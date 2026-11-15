"use strict";
/* Period search for the selected object, and comparison with a catalogued
 * period where one is known.
 *
 * Works on the measurements made by measure.js: brightness relative to the
 * object's own spectrum, which removes the effect of wavelength. Nightly
 * means are fitted with a sinusoid at each trial period; the best period and
 * its uncertainty come from the chi-square curve. Nothing here is fixed in
 * advance: the result is whatever the loaded exposures give.
 *
 * Uses the globals of app.js and Measure (state, $, fmtDate, Measure).
 */

const Period = (() => {
  const P_MIN = 120, P_MAX = 500, P_STEP = 1; // trial periods, days
  const MIN_NIGHTS = 30; // fewest nights of data for a search
  const MIN_MONTHS = 6; // fewest separate 30-day intervals containing data
  const MIN_PER_NIGHT = 2; // measurements needed for a night to count
  const ERROR_FLOOR = 0.03; // 3% added in quadrature to each nightly mean
  const CLIP = 5; // nights further than this many robust sigma from the median are set aside
  const REPLICATE = 0.1; // the two half-samples must each give the period to within this fraction
  const MIN_CYCLES = 1.5; // the data must span this many cycles of the period found
  const MIN_EXPLAINED = 0.3; // a period is reported only if it explains this much variance
  const MATCH_ARCSEC = 20;
  // Stars with a catalogued period to compare against.
  const REFERENCES = [
    { name: "R Draconis", type: "Mira variable", ra: 248.1675, dec: 66.755, period: 245.6,
      source: "General Catalogue of Variable Stars", url: "http://www.sai.msu.su/gcvs/gcvs/" },
  ];
  const DATA = "#3987e5", MODEL = "#d95926", DIM = "#94a0b8";
  let analysis = null;

  /** Weighted least-squares fit of a + b sin + c cos at period P. */
  function sineFit(t, y, w, P) {
    const A = [[0, 0, 0], [0, 0, 0], [0, 0, 0]], b = [0, 0, 0];
    for (let k = 0; k < t.length; k++) {
      const f = [1, Math.sin((2 * Math.PI * t[k]) / P), Math.cos((2 * Math.PI * t[k]) / P)];
      for (let i = 0; i < 3; i++) { b[i] += w[k] * f[i] * y[k]; for (let j = 0; j < 3; j++) A[i][j] += w[k] * f[i] * f[j]; }
    }
    for (let i = 0; i < 3; i++) for (let j = i + 1; j < 3; j++) {
      const m = A[j][i] / A[i][i];
      for (let k = 0; k < 3; k++) A[j][k] -= m * A[i][k];
      b[j] -= m * b[i];
    }
    const x = [0, 0, 0];
    for (let i = 2; i >= 0; i--) x[i] = (b[i] - A[i].slice(i + 1).reduce((s, a, k) => s + a * x[i + 1 + k], 0)) / A[i][i];
    let chi = 0;
    for (let k = 0; k < t.length; k++) chi += w[k] * (y[k] - model(x, P, t[k])) ** 2;
    return { chi, coef: x, amplitude: Math.hypot(x[1], x[2]) };
  }
  const model = (x, P, t) => x[0] + x[1] * Math.sin((2 * Math.PI * t) / P) + x[2] * Math.cos((2 * Math.PI * t) / P);

  /** Search the measurements of one object for a period. `points` are
   *  Measure's per-exposure results ({mjd, rel, relErr}). */
  function analyse(points) {
    const nights = new Map();
    for (const p of points) {
      if (!(p.rel > 0) || !(p.relErr > 0)) continue;
      const k = Math.floor(p.mjd), n = nights.get(k) || { sw: 0, swx: 0, n: 0 };
      const w = 1 / p.relErr ** 2;
      n.sw += w; n.swx += w * p.rel; n.n++;
      nights.set(k, n);
    }
    const t = [], y = [], w = [], series = [];
    for (const [k, n] of [...nights].sort((a, b) => a[0] - b[0])) {
      if (n.n < MIN_PER_NIGHT) continue;
      const rel = n.swx / n.sw, relErr = Math.hypot(1 / Math.sqrt(n.sw), ERROR_FLOOR);
      const magErr = (relErr / rel) * 1.0857;
      t.push(k); y.push(-2.5 * Math.log10(rel)); w.push(1 / magErr ** 2);
      series.push({ mjd: k, rel, relErr, n: n.n });
    }
    const out = { nights: t.length, series, measurements: points.length };
    if (!t.length) return { ...out, status: "none", why: "This object has no usable brightness measurements to search." };
    out.span = t[t.length - 1] - t[0];
    out.months = new Set(t.map((d) => Math.floor((d - t[0]) / 30))).size;
    if (out.nights < MIN_NIGHTS || out.months < MIN_MONTHS) {
      return { ...out, status: "sparse", why: `A period search needs at least ${MIN_NIGHTS} nights of data spread over at least ${MIN_MONTHS} separate months. This object has ${out.nights} night${out.nights === 1 ? "" : "s"} in ${out.months} month${out.months === 1 ? "" : "s"}. Most of the sky is revisited only every six months, so this is common; a longer date range, or a field near the ecliptic poles, gives more.` };
    }
    // Nights far from the rest (a cosmic ray, a neighbour's scattered light)
    // would steer a least-squares fit, so they are set aside first.
    const sorted = [...y].sort((p, q) => p - q), med = sorted[sorted.length >> 1];
    const mad = [...y].map((v) => Math.abs(v - med)).sort((p, q) => p - q)[y.length >> 1];
    const keep = y.map((v) => Math.abs(v - med) <= CLIP * 1.4826 * mad + 1e-9);
    const pick = (arr, test) => arr.filter((_, k) => keep[k] && test(k));
    const T = pick(t, () => true), Y = pick(y, () => true), W = pick(w, () => true);
    out.clipped = t.length - T.length;
    series.forEach((s, k) => { s.clipped = !keep[k]; });

    const search = (tt, yy, ww) => {
      let best = null;
      const scan = [];
      for (let P = P_MIN; P <= P_MAX; P += P_STEP) {
        const f = sineFit(tt, yy, ww, P);
        scan.push([P, f.chi]);
        if (!best || f.chi < best.chi) best = { P, ...f };
      }
      // Chi-square of a constant (the weighted mean), for comparison.
      let sw = 0, swy = 0, flat = 0;
      for (let k = 0; k < tt.length; k++) { sw += ww[k]; swy += ww[k] * yy[k]; }
      for (let k = 0; k < tt.length; k++) flat += ww[k] * (yy[k] - swy / sw) ** 2;
      return { best, scan, flat, explained: 1 - best.chi / flat };
    };
    const { best, scan, flat, explained } = search(T, Y, W);
    // Replication: the same period must come out of two interleaved halves
    // of the nights analysed separately.
    let n = 0;
    const order = keep.map((ok) => (ok ? n++ : -1));
    const halves = [0, 1].map((h) => search(pick(t, (k) => order[k] % 2 === h), pick(y, (k) => order[k] % 2 === h), pick(w, (k) => order[k] % 2 === h)));
    out.halves = halves.map((h) => ({ period: h.best.P, explained: h.explained }));
    const replicated = halves.every((h) => Math.abs(h.best.P - best.P) <= REPLICATE * best.P && h.explained >= MIN_EXPLAINED);
    // Scale so the best fit has a reduced chi-square of one, then take the
    // periods within delta chi-square of one as the 1-sigma range.
    const scale = best.chi / Math.max(1, T.length - 4);
    const inside = scan.filter(([, chi]) => (chi - best.chi) / scale <= 1).map(([P]) => P);
    Object.assign(out, {
      period: best.P, lo: Math.min(...inside), hi: Math.max(...inside), coef: best.coef,
      amplitude: 2 * best.amplitude, explained,
      scan: scan.map(([P, chi]) => [P, Math.max(0, 1 - chi / flat)]),
    });
    if (best.P <= P_MIN || best.P >= P_MAX) return { ...out, status: "edge", why: `The best fit lies at the edge of the range searched (${P_MIN} to ${P_MAX} days), so no period is reported.` };
    if (out.explained < MIN_EXPLAINED) return { ...out, status: "weak", why: `No convincing periodic variation between ${P_MIN} and ${P_MAX} days: the best sinusoid explains only ${Math.max(0, out.explained * 100).toFixed(0)}% of the variation.` };
    if (out.span < MIN_CYCLES * best.P) return { ...out, status: "short", why: `The best fit is ${best.P} days, but the observations span only ${out.span} days, ${(out.span / best.P).toFixed(1)} cycles of that length. At least ${MIN_CYCLES} cycles must be observed before a period is reported, because a single rise and fall can be fitted by many different periods.` };
    if (!replicated) return { ...out, status: "unconfirmed", why: `A sinusoid of ${best.P} days fits all the nights together, but analysing every other night separately gives ${out.halves[0].period} and ${out.halves[1].period} days. A real period would be found in both halves, so none is reported.` };
    return { ...out, status: "ok" };
  }

  /** Compare a recovered period with a reference value. */
  function compare(a, ref) {
    const diff = a.period - ref.period;
    // The uncertainty on the side facing the reference.
    const sigma = Math.max(P_STEP / 2, diff > 0 ? a.period - a.lo : a.hi - a.period);
    const nSigma = Math.abs(diff) / sigma;
    return {
      diff: Math.abs(diff), percent: (Math.abs(diff) / ref.period) * 100, nSigma,
      level: ref.period >= a.lo && ref.period <= a.hi ? "consistent" : nSigma <= 2 ? "marginal" : "inconsistent",
    };
  }

  function referenceFor(sky) {
    if (!sky) return null;
    return REFERENCES.find((r) => Math.hypot((sky.ra - r.ra) * Math.cos((r.dec * Math.PI) / 180), sky.dec - r.dec) * 3600 <= MATCH_ARCSEC) || null;
  }

  // ------------------------------------------------------------- panel ---

  function drawCharts(a, ref) {
    // Brightness against date, with the fitted sinusoid.
    const t0 = a.series[0].mjd, t1 = a.series[a.series.length - 1].mjd, pad = (t1 - t0) * 0.03 + 2;
    const pct = (rel) => (rel - 1) * 100;
    const curve = [];
    for (let t = t0 - pad; t <= t1 + pad; t += 2) curve.push([t, pct(10 ** (-0.4 * model(a.coef, a.period, t)))]);
    const ys = [...a.series.filter((s) => !s.clipped).flatMap((s) => [pct(s.rel + s.relErr), pct(s.rel - s.relErr)]), ...curve.map((c) => c[1])];
    const lim = Math.max(10, ...ys.map(Math.abs)) * 1.08;
    const ticks = [];
    for (let k = 0; k < 4; k++) { const v = t0 + ((t1 - t0) * (k + 0.5)) / 4; ticks.push({ value: v, label: fmtDate(v).slice(0, 7) }); } // inset so the end labels are not cut off
    let ch = Measure.chart("period-lc", { x: [t0 - pad, t1 + pad], y: [-lim, lim], xTicks: ticks, xLabel: "Date of observation", yLabel: "Change from average (%)", zero: 0 });
    ch.ctx.strokeStyle = MODEL; ch.ctx.lineWidth = 2; ch.ctx.beginPath();
    curve.forEach(([t, v], i) => (i ? ch.ctx.lineTo(ch.X(t), ch.Y(v)) : ch.ctx.moveTo(ch.X(t), ch.Y(v))));
    ch.ctx.stroke();
    for (const s of a.series) {
      if (Math.abs(pct(s.rel)) > lim) continue; // an outlying night beyond the plotted range
      ch.mark(ch.X(s.mjd), ch.Y(pct(s.rel)), DATA, { r: 3, hollow: s.clipped, err: s.clipped ? 0 : s.relErr * 100 * ch.yScale,
        html: `<b>${fmtDate(s.mjd)}</b><br>${pct(s.rel) >= 0 ? "+" : "−"}${Math.abs(pct(s.rel)).toFixed(1)} ± ${(s.relErr * 100).toFixed(1)}%<br>${s.n} exposures that night${s.clipped ? "<br>set aside as an outlier" : ""}` });
    }

    // How well each trial period fits.
    const top = Math.max(0.05, ...a.scan.map((s) => s[1])) * 100 * 1.15;
    ch = Measure.chart("period-scan", { x: [P_MIN, P_MAX], y: [0, top], xTicks: [150, 200, 250, 300, 350, 400, 450, 500].map((v) => ({ value: v, label: String(v) })), xLabel: "Trial period (days)", yLabel: "Variation explained (%)" });
    const c = ch.ctx;
    if (a.status === "ok") { c.fillStyle = "rgba(57,135,229,0.16)"; c.fillRect(ch.X(a.lo), ch.box.t, Math.max(2, ch.X(a.hi) - ch.X(a.lo)), ch.box.b - ch.box.t); }
    c.strokeStyle = DATA; c.lineWidth = 2; c.beginPath();
    a.scan.forEach(([P, v], i) => (i ? c.lineTo(ch.X(P), ch.Y(v * 100)) : c.moveTo(ch.X(P), ch.Y(v * 100))));
    c.stroke();
    c.font = "11px system-ui, sans-serif"; c.textBaseline = "top";
    if (ref) {
      c.strokeStyle = MODEL; c.lineWidth = 1.5; c.setLineDash([5, 4]);
      c.beginPath(); c.moveTo(ch.X(ref.period), ch.box.t); c.lineTo(ch.X(ref.period), ch.box.b); c.stroke();
      c.setLineDash([]);
      c.fillStyle = "#e6e9f0"; c.textAlign = "left";
      c.fillText(`catalogued ${ref.period} d`, ch.X(ref.period) + 6, ch.box.t + 2);
    }
    if (a.status === "ok") {
      c.fillStyle = "#e6e9f0"; c.textAlign = "right";
      c.fillText(`recovered ${a.period} d`, ch.X(a.period) - 6, ch.box.t + (ref ? 16 : 2));
    }
    for (const [P, v] of a.scan) if (P % 5 === 0) ch.mark(ch.X(P), ch.Y(v * 100), "rgba(0,0,0,0)", { r: 0, html: `<b>${P} days</b><br>explains ${(v * 100).toFixed(0)}% of the variation` });
  }

  function update() {
    const box = $("period-body"), charts = $("period-charts");
    const sky = window.Measure && Measure.picked(), res = window.Measure && Measure.result();
    const ref = referenceFor(sky);
    $("period-title").textContent = ref ? "Validation: recovering a known period" : "Period search";
    charts.hidden = true;
    analysis = null;
    if (!state.field || !state.exposures.length) {
      box.innerHTML = `<p class="hint">${state.loading ? "Loading…" : "Load a field, then click a star to search its brightness for a repeating pattern."}</p>`;
      return;
    }
    const total = state.field.total;
    if (state.loading) {
      box.innerHTML = `<p class="status"><span class="spinner" aria-hidden="true"></span> Loading ${ref ? ref.name : "exposures"}: ${state.exposures.length} of ${total}</p>
        <p class="hint">${ref ? `When every exposure has arrived, Ingress will measure ${ref.name} on each one and search for its pulsation period. The catalogued value is ${ref.period} days.` : "The period search runs when every exposure has arrived."}</p>`;
      return;
    }
    if (!sky || !res || !res.points.length) {
      box.innerHTML = `<p class="hint">Click a star in the image to measure it and search its brightness for a repeating pattern. Stars with a catalogued period are compared with it automatically.</p>`;
      return;
    }
    const a = (analysis = analyse(res.points));
    const used = `${a.nights} night${a.nights === 1 ? "" : "s"}${a.span ? ` spanning ${a.span} days` : ""}${a.clipped ? ` (${a.clipped} outlying night${a.clipped > 1 ? "s" : ""} set aside)` : ""}, from ${res.points.length} exposures${state.field.found > state.field.total ? ` (${state.field.total.toLocaleString("en-US")} of ${state.field.found.toLocaleString("en-US")} available were loaded)` : ""}${res.geometry.halo ? ", measured from the star's halo because its centre saturates" : ""}`;
    if (a.status !== "ok") {
      box.innerHTML = `
        <p class="verdict none"><span aria-hidden="true">∅</span> No period reported</p>
        <p>${a.why}</p>
        ${ref ? `<p>${ref.name} has a catalogued period of ${ref.period} days, but these observations cannot test it. Widen the date range to include more of them.</p>` : ""}
        <p class="hint">Data: ${used}.</p>`;
      if (a.scan) { charts.hidden = false; drawCharts(a, ref); }
      return;
    }
    const range = `${a.lo}–${a.hi} days`;
    const unc = `+${a.hi - a.period} ⁄ −${a.period - a.lo} days`;
    if (!ref) {
      box.innerHTML = `
        <p class="verdict info"><span aria-hidden="true">∿</span> Candidate periodic variation</p>
        <div class="figures"><div><small>Recovered period</small><b>${a.period} days</b><i>1σ range ${range}</i></div></div>
        <p>The best-fitting sinusoid explains ${(a.explained * 100).toFixed(0)}% of the variation, with a peak-to-peak amplitude of ${a.amplitude.toFixed(2)} mag, and analysing every other night separately gives ${a.halves[0].period} and ${a.halves[1].period} days. No catalogued period is on file for this object to compare with.</p>
        <p class="caution"><span aria-hidden="true">⚠</span> Treat this as a candidate. Scattered light from a nearby bright star changes as the telescope's orientation turns through the year and can imitate a slow variation. Check the images and the object's surroundings before believing it.</p>
        <p class="hint">Data: ${used}. Periods from ${P_MIN} to ${P_MAX} days were tried.</p>`;
    } else {
      const c = compare(a, ref);
      const verdict = {
        consistent: ["ok", "✓", "Consistent with the catalogued period", `The catalogued ${ref.period} days lies inside the 1σ range of the recovered period (${range}).`],
        marginal: ["warn", "⚠", "Marginally consistent", `The catalogued ${ref.period} days lies ${c.nSigma.toFixed(1)}σ from the recovered period, outside its 1σ range (${range}) but within 2σ.`],
        inconsistent: ["bad", "✕", "Not consistent with the catalogued period", `The catalogued ${ref.period} days lies ${c.nSigma.toFixed(1)}σ from the recovered period (1σ range ${range}).`],
      }[c.level];
      box.innerHTML = `
        <p class="verdict ${verdict[0]}"><span aria-hidden="true">${verdict[1]}</span> ${verdict[2]}</p>
        <div class="figures">
          <div><small>Recovered from SPHEREx</small><b>${a.period} days</b><i>1σ range ${range}</i></div>
          <div><small>Catalogued</small><b>${ref.period} days</b><i><a href="${ref.url}" target="_blank" rel="noopener">${ref.source}</a></i></div>
        </div>
        <table>
          <tr><th>Absolute difference</th><td>${c.diff.toFixed(1)} days</td></tr>
          <tr><th>Percentage difference</th><td>${c.percent.toFixed(2)}%</td></tr>
          <tr><th>Measurement uncertainty (1σ)</th><td>${unc}</td></tr>
          <tr><th>Difference ÷ uncertainty</th><td>${c.nSigma.toFixed(2)}σ</td></tr>
        </table>
        <p>${verdict[3]}</p>
        <p>${ref.name} is a ${ref.type}: a pulsating giant star whose brightness rises and falls on a regular cycle. Ingress computed this just now from the exposures on screen, with no period supplied in advance. Analysing every other night separately gives ${a.halves[0].period} and ${a.halves[1].period} days, so the result does not rest on a few nights. The best sinusoid explains ${(a.explained * 100).toFixed(0)}% of the variation, with a peak-to-peak amplitude of ${a.amplitude.toFixed(2)} mag.</p>
        <p class="hint">Data: ${used}. Trial periods from ${P_MIN} to ${P_MAX} days in ${P_STEP}-day steps. Changing the dates or loading another star recalculates everything. <button type="button" class="link" data-dialog="dialog-method" id="period-method">How this is calculated</button></p>`;
      $("period-method").onclick = () => document.querySelector('.top-actions [data-dialog="dialog-method"]').click();
    }
    charts.hidden = false;
    drawCharts(a, ref);
  }

  function init() {
    for (const id of ["period-lc", "period-scan"]) {
      $(id).addEventListener("mousemove", (ev) => Measure.hover(id, ev, "period-tip", "period"));
      $(id).addEventListener("mouseleave", () => { $("period-tip").hidden = true; });
    }
    window.addEventListener("resize", () => { if (analysis && analysis.scan) update(); });
    update();
  }

  return { analyse, compare, update, init, referenceFor, result: () => analysis, REFERENCES };
})();
window.Period = Period;
