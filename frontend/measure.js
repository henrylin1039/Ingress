"use strict";
/* Click-to-measure: aperture photometry of one object on every exposure
 * already in memory, shown as a spectrum per survey pass and as brightness
 * over time compared at matched wavelengths.
 *
 * Uses the globals of app.js (state, cfg, $, crop, passStarts, framesFor,
 * fmtDate, fmtSexa, draw).
 */

const Measure = (() => {
  // Aperture geometry on the sky. SPHEREx pixels are 6.15", so the aperture
  // is about three native pixels across.
  const R_APERTURE = 9.5, R_SKY_IN = 22, R_SKY_OUT = 34; // arcsec
  const MAX_TRACK = 18; // arcsec an object may drift from the clicked point
  const NATIVE_PIXEL = 6.15; // arcsec
  const ERROR_FLOOR = 0.01; // 1% added in quadrature to relative brightness
  // Colour identifies the survey pass in both charts (validated for colour
  // vision deficiency against the panel background); later passes are grey.
  const PASS_COLORS = ["#3987e5", "#d95926", "#199e70", "#c98500", "#d55181", "#0ca30c", "#9085e9", "#e66767"];
  const OTHER_COLOR = "#8a90a0";
  const DIM = "#8e97ab", GRID = "rgba(236,234,228,0.1)", AXIS = "rgba(236,234,228,0.36)", SURFACE = "#05060a";
  const CHART_FONT = '12px "Geist Mono", ui-monospace, Menlo, monospace';

  let pick = null; // {i, j} clicked position in grid pixels
  let result = null;

  const colorOf = (pass) => PASS_COLORS[pass] || OTHER_COLOR;
  const median = (a) => {
    if (!a.length) return NaN;
    const s = Float64Array.from(a).sort();
    return s.length % 2 ? s[s.length >> 1] : (s[s.length / 2 - 1] + s[s.length / 2]) / 2;
  };

  /** Pixel offsets inside a ring r0 <= r < r1 (in grid pixels). */
  function ring(r0, r1) {
    const out = [], R = Math.ceil(r1);
    for (let dy = -R; dy <= R; dy++) for (let dx = -R; dx <= R; dx++) {
      const d = Math.hypot(dx, dy);
      if (d >= r0 && d < r1) out.push([dx, dy]);
    }
    return out;
  }

  /** Follow the object in one pass's stacked image: brightness-weighted
   *  centroid near the clicked point, so moving objects stay in the aperture. */
  function trackCentre(img, n, scale) {
    const rAp = R_APERTURE / scale, offs = ring(0, rAp), sky = ring(R_SKY_IN / scale, R_SKY_OUT / scale);
    let ci = pick.i, cj = pick.j;
    for (let iter = 0; iter < 4; iter++) {
      const bi = Math.round(ci), bj = Math.round(cj);
      const skyVals = [];
      for (const [dx, dy] of sky) {
        const x = bi + dx, y = bj + dy;
        if (x >= 0 && y >= 0 && x < n && y < n && img[y * n + x] === img[y * n + x]) skyVals.push(img[y * n + x]);
      }
      const bkg = median(skyVals);
      const sigma = 1.4826 * median(skyVals.map((v) => Math.abs(v - bkg)));
      let sw = 0, sx = 0, sy = 0, peak = 0;
      for (const [dx, dy] of offs) {
        const x = bi + dx, y = bj + dy;
        if (x < 0 || y < 0 || x >= n || y >= n) continue;
        const w = img[y * n + x] - bkg;
        if (!(w > 0)) continue;
        sw += w; sx += w * x; sy += w * y; peak = Math.max(peak, w);
      }
      if (!(sw > 0) || !(peak > 4 * sigma)) break; // nothing clearly there: stay put
      const ni = sx / sw, nj = sy / sw;
      if (Math.hypot(ni - pick.i, nj - pick.j) * scale > MAX_TRACK) break;
      const moved = Math.hypot(ni - ci, nj - cj);
      ci = ni; cj = nj;
      if (moved < 0.05) break;
    }
    return { i: ci, j: cj };
  }

  function compute() {
    result = null;
    const f = state.field;
    if (!pick || !f || !state.exposures.length) return;
    const n = f.n, scale = f.scale_arcsec, c = (n - 1) / 2;
    const starts = passStarts();
    const passOf = (mjd) => Math.max(0, starts.filter((s) => s <= mjd).length - 1);

    // Where the object is in each pass. While exposures are still arriving
    // the clicked position is used as it is.
    const centres = starts.map(() => ({ i: pick.i, j: pick.j }));
    if (!state.loading) {
      for (const fr of framesFor("all", "pass")) centres[passOf(fr.t0)] = trackCentre(fr.img, n, scale);
    }

    const omega = (scale / 206264.806) ** 2; // steradians per grid pixel
    const toMilliJy = omega * 1e9; // MJy/sr summed over pixels -> mJy
    const correlated = Math.max(1, NATIVE_PIXEL / scale); // resampling shares noise between pixels

    /** Photometry of every exposure in the ring g.in..g.out (arcsec), with
     *  the sky level taken from g.skyIn..g.skyOut. */
    const measureAll = (g) => {
      const ap = ring(g.in / scale, g.out / scale), sky = ring(g.skyIn / scale, g.skyOut / scale);
      const out = [];
      for (const e of state.exposures) {
        const pass = passOf(e.meta.mjd), ctr = centres[pass];
        const bi = Math.round(ctr.i), bj = Math.round(ctr.j), d = e.data;
        const skyVals = [];
        for (const [dx, dy] of sky) {
          const x = bi + dx, y = bj + dy;
          if (x >= 0 && y >= 0 && x < n && y < n && d[y * n + x] === d[y * n + x]) skyVals.push(d[y * n + x]);
        }
        if (skyVals.length < sky.length * 0.5) continue;
        const bkg = median(skyVals);
        const sigma = 1.4826 * median(skyVals.map((v) => Math.abs(v - bkg)));
        let sum = 0, count = 0, peak = -Infinity, atPeak = 0;
        const vals = [];
        for (const [dx, dy] of ap) {
          const x = bi + dx, y = bj + dy;
          if (x < 0 || y < 0 || x >= n || y >= n) continue;
          const v = d[y * n + x];
          if (v !== v) continue;
          sum += v; count++; vals.push(v);
          if (v > peak) { peak = v; atPeak = 1; } else if (v === peak) atPeak++;
        }
        if (count < ap.length * 0.85) continue;
        // In a halo ring the median around the ring is used in place of the
        // sum, so neighbouring stars and diffraction spikes do not count.
        let noise = 1;
        if (g.halo) { sum = median(vals) * count; noise = 1.2533; }
        out.push({
          mjd: e.meta.mjd, det: e.meta.det, pass,
          wave: e.meta.wave + e.meta.dwave_dx * (ctr.i - c) + e.meta.dwave_dy * (ctr.j - c),
          bandwidth: e.meta.bandwidth,
          flux: (sum - count * bkg) * toMilliJy,
          err: noise * sigma * Math.sqrt(count * (1 + count / skyVals.length)) * correlated * toMilliJy,
          // Saturated cores are filled with one constant value by the server.
          saturated: atPeak >= 4,
        });
      }
      return out;
    };

    let geometry = { in: 0, out: R_APERTURE, skyIn: R_SKY_IN, skyOut: R_SKY_OUT, halo: false };
    let points = measureAll(geometry);
    if (points.filter((p) => p.saturated).length > 0.3 * points.length) {
      // A star too bright for its core to be measured. Its halo scales with
      // its brightness, so measure a ring just outside the largest saturated
      // core seen in any exposure instead.
      // Size of each exposure's saturated core, as the radius of a circle
      // of the same area; a high percentile ignores the odd bad exposure.
      const reach = Math.ceil(60 / scale), radii = [];
      for (const e of state.exposures) {
        const ctr = centres[passOf(e.meta.mjd)], bi = Math.round(ctr.i), bj = Math.round(ctr.j), d = e.data;
        let top = -Infinity, area = 0;
        for (let dy = -3; dy <= 3; dy++) for (let dx = -3; dx <= 3; dx++) {
          const x = bi + dx, y = bj + dy;
          if (x >= 0 && y >= 0 && x < n && y < n && d[y * n + x] > top) top = d[y * n + x];
        }
        if (!(top > 0)) continue;
        for (let dy = -reach; dy <= reach; dy++) for (let dx = -reach; dx <= reach; dx++) {
          const x = bi + dx, y = bj + dy;
          if (x >= 0 && y >= 0 && x < n && y < n && d[y * n + x] === top) area++;
        }
        if (area >= 4) radii.push(Math.sqrt(area / Math.PI) * scale);
      }
      radii.sort((p, q) => p - q);
      const rSat = radii.length ? radii[Math.floor(0.9 * (radii.length - 1))] : R_APERTURE;
      geometry = { in: rSat + 9, out: rSat + 24, skyIn: rSat + 45, skyOut: rSat + 60, halo: true };
      points = measureAll(geometry);
    }
    points.sort((a, b) => a.wave - b.wave);

    // Brightness relative to the object's own typical spectrum: each
    // measurement is divided by the median of all measurements (any pass)
    // at nearly the same wavelength.
    const usable = points.filter((p) => !p.saturated);
    for (const p of points) {
      p.rel = p.relErr = NaN;
      if (p.saturated) continue;
      const near = usable.filter((q) => Math.abs(q.wave - p.wave) <= 0.6 * p.bandwidth);
      if (near.length < 2 || new Set(near.map((q) => q.pass)).size < 2) continue;
      const ref = median(near.map((q) => q.flux));
      if (!(ref > 5 * median(near.map((q) => q.err)))) continue; // too faint here
      p.rel = p.flux / ref;
      p.relErr = p.err / ref;
    }

    const combine = (list) => {
      const ok = list.filter((p) => p.rel === p.rel && p.relErr > 0);
      if (ok.length < 2) return null;
      let sw = 0, swx = 0;
      for (const p of ok) { const w = 1 / p.relErr ** 2; sw += w; swx += w * p.rel; }
      const mean = swx / sw;
      const scatter = ok.length >= 3 ? Math.sqrt(ok.reduce((s, p) => s + (p.rel - mean) ** 2, 0) / (ok.length - 1) / ok.length) : 0;
      return { value: mean, err: Math.hypot(Math.max(Math.sqrt(1 / sw), scatter), ERROR_FLOOR), n: ok.length };
    };
    const passes = starts.map((_, k) => {
      const mine = points.filter((p) => p.pass === k);
      return {
        index: k, n: mine.length,
        t0: Math.min(...mine.map((p) => p.mjd)), t1: Math.max(...mine.map((p) => p.mjd)),
        all: combine(mine),
        bands: [1, 2, 3, 4, 5, 6].map((det) => combine(mine.filter((p) => p.det === det))),
        centre: centres[k],
      };
    }).filter((p) => p.n);

    // Largest change between any two passes, in units of its uncertainty.
    let verdict = null;
    const m = passes.filter((p) => p.all && p.all.n >= 5);
    for (let a = 0; a < m.length; a++) for (let b = a + 1; b < m.length; b++) {
      const diff = m[b].all.value - m[a].all.value, err = Math.hypot(m[a].all.err, m[b].all.err);
      if (!verdict || Math.abs(diff / err) > verdict.sigma) verdict = { a: m[a], b: m[b], diff, err, sigma: Math.abs(diff / err) };
    }
    result = { points, passes, verdict, scale, geometry, saturated: points.filter((p) => p.saturated).length };
  }

  // ------------------------------------------------------------- charts ---

  const hits = new Map(); // canvas id -> [{x, y, html}]

  function niceTicks(lo, hi, count) {
    const span = hi - lo || 1, raw = span / count, mag = 10 ** Math.floor(Math.log10(raw));
    const step = [1, 2, 2.5, 5, 10].map((k) => k * mag).find((s) => s >= raw);
    const ticks = [];
    for (let v = Math.ceil(lo / step) * step; v <= hi + step * 1e-9; v += step) ticks.push(Math.abs(v) < step * 1e-9 ? 0 : v);
    return ticks;
  }
  const fmtTick = (v) => (Math.abs(v) >= 100 ? v.toFixed(0) : Math.abs(v) >= 10 ? v.toFixed(1).replace(/\.0$/, "") : String(+v.toPrecision(3)));

  /** Set up a canvas and draw the frame; returns the plotting helpers. */
  function frame(id, { x, y, xTicks, yLabel, xLabel, zero }) {
    const cv = $(id), dpr = window.devicePixelRatio || 1;
    const W = cv.clientWidth || 600, H = cv.clientHeight || 260;
    cv.width = Math.round(W * dpr); cv.height = Math.round(H * dpr);
    const ctx = cv.getContext("2d");
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, W, H);
    const labelled = id === "spec-chart" || id === "lc-chart";
    const box = { l: 66, r: W - 16, t: labelled ? 34 : 12, b: H - 54 };
    const X = (v) => box.l + ((v - x[0]) / (x[1] - x[0])) * (box.r - box.l);
    const Y = (v) => box.b - ((v - y[0]) / (y[1] - y[0])) * (box.b - box.t);
    ctx.font = CHART_FONT;
    ctx.lineWidth = 1;
    ctx.textBaseline = "middle";
    ctx.textAlign = "right";
    for (const v of niceTicks(y[0], y[1], 5)) {
      // Dotted hairlines, with the zero line drawn solid.
      const solid = zero != null && v === zero;
      ctx.strokeStyle = solid ? AXIS : GRID;
      ctx.setLineDash(solid ? [] : [2, 4]);
      ctx.beginPath(); ctx.moveTo(box.l, Math.round(Y(v)) + 0.5); ctx.lineTo(box.r, Math.round(Y(v)) + 0.5); ctx.stroke();
      ctx.setLineDash([]);
      ctx.fillStyle = DIM;
      ctx.fillText(fmtTick(v), box.l - 12, Y(v));
    }
    ctx.textAlign = "center"; ctx.textBaseline = "top";
    ctx.strokeStyle = AXIS;
    for (const t of xTicks) {
      const tx = Math.round(X(t.value)) + 0.5;
      ctx.beginPath(); ctx.moveTo(tx, box.b + 3); ctx.lineTo(tx, box.b + 8); ctx.stroke();
      ctx.fillStyle = DIM; ctx.fillText(t.label, X(t.value), box.b + 14);
    }
    ctx.fillStyle = DIM;
    ctx.fillText(xLabel, (box.l + box.r) / 2, H - 17);
    ctx.save(); ctx.translate(13, (box.t + box.b) / 2); ctx.rotate(-Math.PI / 2); ctx.textBaseline = "middle";
    ctx.fillText(yLabel, 0, 0); ctx.restore();
    hits.set(id, []);
    const mark = (px, py, color, { r = 3, hollow = false, err = 0, diamond = false, html }) => {
      if (err > 0) {
        const a = Math.max(box.t, py - err), b = Math.min(box.b, py + err);
        ctx.strokeStyle = color; ctx.globalAlpha = 0.6; ctx.lineWidth = 1;
        ctx.beginPath(); ctx.moveTo(px, a); ctx.lineTo(px, b);
        ctx.moveTo(px - 2.5, a); ctx.lineTo(px + 2.5, a); ctx.moveTo(px - 2.5, b); ctx.lineTo(px + 2.5, b);
        ctx.stroke();
        ctx.globalAlpha = 1;
      }
      if (!hollow) { // a soft halo, so points read as light on the dark plot
        ctx.globalAlpha = 0.14; ctx.fillStyle = color;
        ctx.beginPath(); ctx.arc(px, py, Math.min(r * 2.1, 9), 0, 2 * Math.PI); ctx.fill();
        ctx.globalAlpha = 1;
      }
      ctx.beginPath();
      if (diamond) { ctx.moveTo(px, py - r); ctx.lineTo(px + r, py); ctx.lineTo(px, py + r); ctx.lineTo(px - r, py); ctx.closePath(); }
      else ctx.arc(px, py, r, 0, 2 * Math.PI);
      ctx.fillStyle = hollow ? SURFACE : color; ctx.fill();
      ctx.strokeStyle = hollow ? color : SURFACE; ctx.lineWidth = hollow ? 1.5 : 1; ctx.stroke();
      hits.get(id).push({ x: px, y: py, html });
    };
    return { ctx, X, Y, box, mark, yScale: (box.b - box.t) / (y[1] - y[0]) };
  }

  function drawSpectrum() {
    const pts = result.points;
    const hiVals = pts.map((p) => p.flux + (p.saturated ? 0 : p.err)), loVals = pts.map((p) => p.flux - (p.saturated ? 0 : p.err));
    let hi = Math.max(...hiVals), lo = Math.min(0, ...loVals);
    const pad = (hi - lo || 1) * 0.06;
    const ticks = [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5].map((v) => ({ value: v, label: String(v) }));
    const ch = frame("spec-chart", { x: [0.7, 5.05], y: [lo - (lo < 0 ? pad : 0), hi + pad], xTicks: ticks, xLabel: "Wavelength (µm)", yLabel: result.geometry.halo ? "Halo brightness (mJy)" : "Flux (mJy)", zero: 0 });
    // The six detector bands: a tinted column each, keyed in the band's colour.
    const { ctx, box } = ch;
    ctx.font = CHART_FONT; ctx.textAlign = "center"; ctx.textBaseline = "top";
    Object.values(cfg.bands).forEach(([lo, hi], k) => {
      const a = Math.max(box.l, ch.X(lo)), b = Math.min(box.r, ch.X(hi)), color = BAND_COLORS[k + 1] || DIM;
      ctx.fillStyle = color;
      ctx.globalAlpha = 0.05; ctx.fillRect(a, box.t, b - a, box.b - box.t);
      ctx.globalAlpha = 0.9; ctx.fillRect(a + 1, box.t - 9, b - a - 2, 2);
      ctx.globalAlpha = 1; ctx.fillStyle = DIM;
      if (b - a > 54) ctx.fillText(`Band ${k + 1}`, (a + b) / 2, box.t - 27);
      else ctx.fillText(String(k + 1), (a + b) / 2, box.t - 27);
    });
    // A faint line through each pass, so the shape of the spectrum reads at a glance.
    for (const pass of new Set(pts.map((p) => p.pass))) {
      const line = pts.filter((p) => p.pass === pass && !p.saturated);
      if (line.length < 2) continue;
      ctx.strokeStyle = colorOf(pass); ctx.globalAlpha = 0.32; ctx.lineWidth = 1; ctx.lineJoin = "round";
      ctx.beginPath();
      line.forEach((p, k) => (k ? ctx.lineTo(ch.X(p.wave), ch.Y(p.flux)) : ctx.moveTo(ch.X(p.wave), ch.Y(p.flux))));
      ctx.stroke();
      ctx.globalAlpha = 1;
    }
    for (const p of pts) {
      ch.mark(ch.X(p.wave), ch.Y(p.flux), colorOf(p.pass), {
        hollow: p.saturated, err: p.saturated ? 0 : p.err * ch.yScale,
        html: `<b>${p.wave.toFixed(3)} µm</b> · band ${p.det}<br>${p.flux.toPrecision(4)} ± ${p.err.toPrecision(2)} mJy${p.saturated ? " (saturated)" : ""}<br>${fmtDate(p.mjd)} · pass ${p.pass + 1}`,
      });
    }
  }

  function drawLightCurve() {
    const passes = result.passes;
    const vals = [];
    for (const p of passes) for (const b of [p.all, ...p.bands]) if (b) vals.push((b.value - 1) * 100 + b.err * 100, (b.value - 1) * 100 - b.err * 100);
    $("lc-block").hidden = !vals.length;
    if (!vals.length) return;
    const lim = Math.max(5, ...vals.map(Math.abs)) * 1.1;
    const ticks = passes.map((p, k) => ({ value: k + 0.5, label: p.t1 - p.t0 < 1 ? fmtDate(p.t0) : `${fmtDate(p.t0).slice(0, 7)}` }));
    const ch = frame("lc-chart", { x: [0, passes.length], y: [-lim, lim], xTicks: ticks, xLabel: "Survey pass (month of observation)", yLabel: "Change from average (%)", zero: 0 });
    // One column per pass, alternately shaded and topped in the pass's colour.
    passes.forEach((p, k) => {
      const a = ch.X(k), b = ch.X(k + 1);
      ch.ctx.fillStyle = "rgb(236,234,228)";
      ch.ctx.globalAlpha = k % 2 ? 0 : 0.035; ch.ctx.fillRect(a, ch.box.t, b - a, ch.box.b - ch.box.t);
      ch.ctx.fillStyle = colorOf(p.index);
      ch.ctx.globalAlpha = 0.9; ch.ctx.fillRect(a + 1, ch.box.t - 9, b - a - 2, 2);
      ch.ctx.globalAlpha = 1; ch.ctx.fillStyle = DIM;
      ch.ctx.font = CHART_FONT; ch.ctx.textAlign = "center"; ch.ctx.textBaseline = "top";
      ch.ctx.fillText(`Pass ${p.index + 1}`, (a + b) / 2, ch.box.t - 27);
    });
    passes.forEach((p, k) => {
      const color = colorOf(p.index), when = `${fmtDate(p.t0)} → ${fmtDate(p.t1)}`;
      p.bands.forEach((b, d) => {
        if (!b) return;
        ch.mark(ch.X(k + 0.14 + d * 0.1), ch.Y((b.value - 1) * 100), color, {
          r: 3.5, err: b.err * 100 * ch.yScale,
          html: `<b>Band ${d + 1}</b> · pass ${p.index + 1}<br>${((b.value - 1) * 100).toFixed(1)} ± ${(b.err * 100).toFixed(1)}%<br>${b.n} exposures · ${when}`,
        });
      });
      if (p.all) {
        ch.mark(ch.X(k + 0.84), ch.Y((p.all.value - 1) * 100), color, {
          r: 6, diamond: true, err: p.all.err * 100 * ch.yScale,
          html: `<b>All bands</b> · pass ${p.index + 1}<br>${((p.all.value - 1) * 100).toFixed(1)} ± ${(p.all.err * 100).toFixed(1)}%<br>${p.all.n} exposures · ${when}`,
        });
      }
    });
  }

  function render() {
    const body = $("measure-body");
    body.hidden = !result || !result.points.length;
    $("measure-hint").hidden = !!pick && !body.hidden;
    if (pick && body.hidden) {
      $("measure-hint").hidden = false;
      $("measure-hint").textContent = state.exposures.length ? "No usable measurements at that position. Click a different object." : "Click any object in the image to measure it.";
      return;
    }
    if (!result || body.hidden) { $("measure-hint").textContent = "Click any object in the image to measure its spectrum and how its brightness has changed."; return; }

    const sky = pixelToSky(pick.i, pick.j);
    const v = result.verdict;
    let line;
    if (state.loading) line = "Still loading exposures; the result will firm up when loading finishes.";
    else if (!v) line = "Too faint, too saturated or seen in only one pass, so brightness cannot be compared between passes.";
    else {
      const pct = `${Math.abs(v.diff * 100).toFixed(1)} ± ${(v.err * 100).toFixed(1)}%`;
      const span = `between ${fmtDate(v.a.t0).slice(0, 7)} and ${fmtDate(v.b.t0).slice(0, 7)}`;
      if (v.sigma >= 5) line = `<b>Variable:</b> ${v.diff > 0 ? "brightened" : "faded"} by ${pct} ${span} (${v.sigma.toFixed(0)}σ).`;
      else if (v.sigma >= 3) line = `<b>Possibly variable:</b> ${v.diff > 0 ? "brighter" : "fainter"} by ${pct} ${span} (${v.sigma.toFixed(1)}σ).`;
      else line = `<b>No significant change:</b> the largest difference is ${pct} ${span} (${v.sigma.toFixed(1)}σ).`;
    }
    let sat = "";
    if (result.geometry.halo) sat = ` <span class="warn">This star is too bright for SPHEREx to measure its centre, so the ring of light around it is measured instead. Changes over time can still be followed this way, less precisely than for an unsaturated star, and the method has been checked on one star only. The flux values are a fraction of the star's true brightness.</span>`;
    else if (result.saturated) sat = ` <span class="warn">${result.saturated} of ${result.points.length} measurements are saturated (hollow points) and are left out of the comparison.</span>`;
    $("measure-method").textContent = result.geometry.halo
      ? `Typical brightness in a ring from ${result.geometry.in.toFixed(0)}″ to ${result.geometry.out.toFixed(0)}″ around the star (solid rings on the image), sky from the dashed rings. Error bars come from the scatter of the surrounding sky.`
      : "Aperture photometry in a 9.5″ circle (solid ring on the image), sky from the dashed rings. Error bars come from the scatter of the surrounding sky. Fluxes are not corrected for light falling outside the aperture: compared with the 2MASS and WISE catalogues they read about 4–13% low.";
    // What the numbers rest on. The band buttons only change the images.
    const band = state.band === "all" ? "" : state.band.startsWith("c:") ? `the ${state.band.slice(2).replace("-", "–")} µm range` : `Band ${state.band.slice(1)}`;
    const scope = `Measured on all loaded wavelengths (0.75–5.0 µm)${band ? `, not only ${band}: the Wavelength buttons change the images, not the measurements` : "; the Wavelength buttons change the images, not the measurements"}. Each ± reflects background noise only, so the true uncertainty is larger for bright stars.`;
    $("measure-summary").innerHTML = `${fmtSexa(sky.ra, sky.dec)} · ${result.points.length} measurements<br>${line}${sat}<span class="scope">${scope}</span>`;

    const legend = $("measure-legend");
    legend.textContent = "";
    for (const p of result.passes) {
      const item = document.createElement("span");
      item.className = "key";
      item.innerHTML = `<i style="background:${colorOf(p.index)}"></i>Pass ${p.index + 1}: ${fmtDate(p.t0)}${p.t1 - p.t0 >= 1 ? " → " + fmtDate(p.t1) : ""}`;
      legend.append(item);
    }
    drawSpectrum();
    drawLightCurve();
  }

  function pixelToSky(i, j) {
    const f = state.field, rad = Math.PI / 180, as = (f.scale_arcsec / 3600) * rad, c = (f.n - 1) / 2;
    const xi = -(i - c) * as, eta = (j - c) * as, d0 = f.dec * rad;
    const den = Math.cos(d0) - eta * Math.sin(d0), dra = Math.atan2(xi, den);
    return { ra: (((f.ra + dra / rad) % 360) + 360) % 360, dec: Math.atan(((Math.sin(d0) + eta * Math.cos(d0)) * Math.cos(dra)) / den) / rad };
  }

  function skyToPixel(ra, dec) {
    const f = state.field, rad = Math.PI / 180, c = (f.n - 1) / 2, as = (f.scale_arcsec / 3600) * rad;
    const d0 = f.dec * rad, d = dec * rad, da = (ra - f.ra) * rad;
    const cosc = Math.sin(d0) * Math.sin(d) + Math.cos(d0) * Math.cos(d) * Math.cos(da);
    const xi = (Math.cos(d) * Math.sin(da)) / cosc;
    const eta = (Math.cos(d0) * Math.sin(d) - Math.sin(d0) * Math.cos(d) * Math.cos(da)) / cosc;
    return { i: c - xi / as, j: c + eta / as };
  }

  function copyLink() {
    const btn = $("measure-link"), done = (msg) => { btn.textContent = msg; setTimeout(() => { btn.textContent = "Copy link to this object"; }, 1800); };
    const fallback = () => { window.prompt("Copy this link:", location.href); };
    if (navigator.clipboard && navigator.clipboard.writeText) navigator.clipboard.writeText(location.href).then(() => done("Link copied"), fallback);
    else fallback();
  }

  /** Show the tooltip of the mark nearest the pointer. `tipId` is the
   *  tooltip element and `hostId` the positioned element that contains it. */
  function hover(id, ev, tipId = "chart-tip", hostId = "measure") {
    const tip = $(tipId), r = $(id).getBoundingClientRect();
    const x = ev.clientX - r.left, y = ev.clientY - r.top;
    let best = null, bd = 18 ** 2;
    for (const h of hits.get(id) || []) { const d = (h.x - x) ** 2 + (h.y - y) ** 2; if (d < bd) { bd = d; best = h; } }
    tip.hidden = !best;
    if (!best) return;
    tip.innerHTML = best.html;
    const host = $(hostId).getBoundingClientRect();
    tip.style.left = `${Math.min(r.left - host.left + best.x + 12, host.width - 190)}px`;
    tip.style.top = `${r.top - host.top + best.y + 12}px`;
  }

  function downloadCsv() {
    if (!result) return;
    const rows = ["mjd,date_utc,pass,band,wavelength_um,flux_mjy,flux_err_mjy,saturated,relative_brightness,relative_err"];
    for (const p of [...result.points].sort((a, b) => a.mjd - b.mjd)) {
      rows.push([p.mjd.toFixed(5), new Date((p.mjd - 40587) * 86400e3).toISOString(), p.pass + 1, p.det, p.wave.toFixed(4), p.flux.toPrecision(6), p.err.toPrecision(4), p.saturated ? 1 : 0,
        p.rel === p.rel ? p.rel.toFixed(4) : "", p.relErr === p.relErr ? p.relErr.toFixed(4) : ""].join(","));
    }
    const sky = pixelToSky(pick.i, pick.j);
    const a = document.createElement("a");
    a.href = URL.createObjectURL(new Blob([rows.join("\n") + "\n"], { type: "text/csv" }));
    a.download = `spherex_${sky.ra.toFixed(4)}_${sky.dec >= 0 ? "+" : ""}${sky.dec.toFixed(4)}.csv`;
    a.click();
    setTimeout(() => URL.revokeObjectURL(a.href), 1000);
  }

  // ------------------------------------------------------------ public ---

  function refresh() { compute(); render(); if (window.Period) Period.update(); if (window.Insight) Insight.update(); }

  return {
    refresh,
    chart: frame, // shared with the period panel
    hover,
    clear() { pick = null; result = null; render(); if (window.Period) Period.update(); if (window.Insight) Insight.update(); },
    /** Sky position of the selected object, or null. */
    picked() { return pick && state.field ? pixelToSky(pick.i, pick.j) : null; },
    /** The current measurements (see compute), or null. */
    result() { return result; },
    /** Select the object at a sky position, e.g. from a shared link. */
    pickSky(ra, dec) { const p = skyToPixel(ra, dec); Measure.pickAt(p.i, p.j); },
    /** Called with the pointer position in grid pixels. */
    pickAt(i, j) {
      const n = state.field && state.field.n;
      if (!n || i < 0 || j < 0 || i > n - 1 || j > n - 1) return;
      pick = { i, j };
      refresh();
      draw();
      syncUrl();
    },
    /** Draw the aperture on the image. `toCanvas` maps grid pixels to canvas. */
    overlay(ctx, toCanvas, unit) {
      if (!pick || !state.field) return;
      const fr = state.frames[state.frame];
      let ctr = pick;
      if (result && fr) {
        const p = result.passes.filter((q) => q.t0 - 1 <= fr.t1).pop() || result.passes[0];
        if (p) ctr = p.centre;
      }
      const [x, y] = toCanvas(ctr.i, ctr.j), k = toCanvas(1, 0)[0] - toCanvas(0, 0)[0], s = state.field.scale_arcsec;
      ctx.save();
      ctx.strokeStyle = "#5fd3c4"; ctx.lineWidth = 1.5 * unit;
      const g = (result && result.geometry) || { in: 0, out: R_APERTURE, skyIn: R_SKY_IN, skyOut: R_SKY_OUT };
      for (const r of [g.in, g.out]) if (r > 0) { ctx.beginPath(); ctx.arc(x, y, (r / s) * k, 0, 2 * Math.PI); ctx.stroke(); }
      ctx.setLineDash([4 * unit, 4 * unit]); ctx.globalAlpha = 0.7;
      for (const r of [g.skyIn, g.skyOut]) { ctx.beginPath(); ctx.arc(x, y, (r / s) * k, 0, 2 * Math.PI); ctx.stroke(); }
      ctx.restore();
    },
    init() {
      for (const id of ["spec-chart", "lc-chart"]) {
        $(id).addEventListener("mousemove", (ev) => hover(id, ev));
        $(id).addEventListener("mouseleave", () => { $("chart-tip").hidden = true; });
      }
      $("measure-csv").onclick = downloadCsv;
      $("measure-clear").onclick = () => { Measure.clear(); draw(); syncUrl(); };
      $("measure-link").onclick = copyLink;
      window.addEventListener("resize", () => { if (result) render(); });
      // Redraw when a chart's own box changes size (not only the window), so a
      // chart is never shown stretched from a drawing made at another width.
      if ("ResizeObserver" in window) {
        let queued = 0;
        const sized = new ResizeObserver(() => {
          if (queued) return;
          queued = requestAnimationFrame(() => {
            queued = 0;
            const stale = ["spec-chart", "lc-chart"].some((id) => {
              const cv = $(id), w = Math.round(cv.clientWidth * (window.devicePixelRatio || 1));
              return cv.clientWidth > 0 && Math.abs(cv.width - w) > 1;
            });
            if (result && stale) render();
          });
        });
        for (const id of ["spec-chart", "lc-chart"]) sized.observe($(id));
      }
      render();
    },
  };
})();
window.Measure = Measure;
