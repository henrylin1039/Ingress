"use strict";
/* Comparing the same spectral region over time.
 *
 * A SPHEREx exposure sees each point of sky at one wavelength, and which
 * wavelengths a field receives differs from one survey pass (epoch) to the
 * next. Averaging "everything in band 3" per epoch therefore mixes a change
 * in wavelength with a change in time. This module builds one frame per
 * epoch from the same set of narrow wavelength slices, weights the slices
 * equally in every epoch, measures how well the epochs line up, and says so
 * when no fair comparison exists.
 *
 * Uses the globals of app.js (state, cfg, $, selected, passStarts, stack,
 * fmtDate).
 */

const Compare = (() => {
  const SLICES_PER_BAND = 6; // each detector band is cut into this many wavelength slices
  const MIN_COVERAGE = 0.5; // an exposure must cover this much of the field to count
  const MIN_SLICES = 0.85; // a pixel needs this fraction of the slices to be shown
  const MAX_NOISE = 3; // exposures noisier than this multiple of typical are left out
  const MIN_SHARED = 0.25; // fewest shared slices, as a fraction of a typical epoch's
  const MIN_STARS = 8; // fewest stars for a trustworthy alignment estimate
  const SHIFT_WARN = 1.5; // arcsec; about a quarter of a SPHEREx pixel
  const WIDTH_WARN = 0.1; // fractional difference in star width between epochs

  const median = (a) => {
    if (!a.length) return NaN;
    const s = [...a].sort((x, y) => x - y);
    return s.length % 2 ? s[s.length >> 1] : (s[s.length / 2 - 1] + s[s.length / 2]) / 2;
  };

  /** Which wavelength slice an exposure falls in. */
  function sliceOf(e) {
    const [lo, hi] = cfg.bands[`D${e.meta.det}`];
    const k = Math.floor(((e.meta.wave - lo) / (hi - lo)) * SLICES_PER_BAND);
    return e.meta.det * 100 + Math.min(SLICES_PER_BAND - 1, Math.max(0, k));
  }

  function build(band) {
    const npx = state.field.n * state.field.n;
    const starts = passStarts();
    const passOf = (mjd) => Math.max(0, starts.filter((s) => s <= mjd).length - 1);
    const epochs = new Map();
    // Slices are weighted equally, so one very noisy exposure (common at the
    // edges of a detector's wavelength range) would spoil a whole frame.
    // Leave out exposures far noisier than is typical for their detector.
    const typicalNoise = {};
    for (const det of [1, 2, 3, 4, 5, 6]) typicalNoise[det] = median(selected(band).filter((e) => e.meta.det === det && e.meta.noise > 0).map((e) => e.meta.noise));
    let noisy = 0;
    for (const e of selected(band)) {
      if (e.meta.coverage < MIN_COVERAGE) continue;
      if (e.meta.noise > MAX_NOISE * typicalNoise[e.meta.det]) { noisy++; continue; }
      const k = passOf(e.meta.mjd);
      if (!epochs.has(k)) epochs.set(k, { index: k, slices: new Map(), exposures: [] });
      const ep = epochs.get(k), s = sliceOf(e);
      if (!ep.slices.has(s)) ep.slices.set(s, []);
      ep.slices.get(s).push(e);
      ep.exposures.push(e);
    }
    const all = [...epochs.values()].sort((a, b) => a.index - b.index);
    for (const ep of all) {
      ep.t0 = Math.min(...ep.exposures.map((e) => e.meta.mjd));
      ep.t1 = Math.max(...ep.exposures.map((e) => e.meta.mjd));
      ep.used = true;
    }
    const union = new Set(all.flatMap((ep) => [...ep.slices.keys()]));

    // Keep the slices present in every epoch. More epochs matter more than
    // more slices, so an epoch is dropped (poorest coverage first) only if
    // keeping it would leave under a quarter of a typical epoch's slices.
    let active = [...all], common = [];
    while (active.length >= 2) {
      common = [...union].filter((s) => active.every((ep) => ep.slices.has(s))).sort((a, b) => a - b);
      const typical = median(active.map((ep) => ep.slices.size));
      if (common.length >= Math.max(1, Math.ceil(typical * MIN_SHARED))) break;
      const worst = active.reduce((w, ep) => (ep.slices.size < w.slices.size || (ep.slices.size === w.slices.size && ep.exposures.length < w.exposures.length) ? ep : w));
      worst.used = false;
      worst.reason = `it covers only ${worst.slices.size} of the ${union.size} wavelength slices seen in this range`;
      active = active.filter((ep) => ep !== worst);
      common = [];
    }

    const report = { band, slicesTotal: union.size, slices: common.length, epochs: all, status: "ok", noisy };
    if (active.length < 2) {
      // Nothing can be compared. Show the best single epoch as a still.
      report.status = "none";
      report.why = all.length < 2
        ? (all.length ? "SPHEREx has observed this field in this wavelength range during only one epoch so far." : "There are no exposures in this wavelength range.")
        : "The epochs that observed this field did not sample overlapping wavelengths in this range.";
      for (const ep of all) ep.used = false;
      const best = all.reduce((b, ep) => (!b || ep.slices.size > b.slices.size ? ep : b), null);
      const frames = [];
      if (best) {
        const img = stack(best.exposures, npx);
        frames.push(frameOf(best, best.exposures, img, false));
        frames[0].mean = img;
        report.still = best.index;
      }
      frames.report = report;
      return frames;
    }

    const frames = [];
    for (const ep of active) {
      // Average within each slice first, then average the slices with equal
      // weight, so every epoch has the same spectral weighting.
      const img = matchedImage(common.map((s) => ep.slices.get(s)), npx);
      const used = common.flatMap((s) => ep.slices.get(s));
      ep.nUsed = used.length;
      ep.wave = common.reduce((t, s) => t + ep.slices.get(s).reduce((a, e) => a + e.meta.wave, 0) / ep.slices.get(s).length, 0) / common.length;
      frames.push(frameOf(ep, used, img, true));
    }
    // Average of the epochs, only where every epoch has data.
    const mean = new Float32Array(npx);
    let shared = 0;
    for (let i = 0; i < npx; i++) {
      let sum = 0;
      for (const f of frames) sum += f.img[i];
      mean[i] = sum / frames.length;
      if (mean[i] === mean[i]) shared++;
    }
    for (const f of frames) f.mean = mean;
    report.sharedPixels = shared / npx;
    const waves = active.map((ep) => ep.wave);
    report.waveSpread = (Math.max(...waves) - Math.min(...waves)) / median(waves);
    const lo = Math.min(...frames.map((f) => f.wmin)), hi = Math.max(...frames.map((f) => f.wmax));
    report.range = [lo, hi];
    diagnose(frames, mean, report);
    frames.report = report;
    return frames;
  }

  /** One epoch's image from its exposures grouped by slice.
   *
   *  Per pixel: values far from the median of all the epoch's exposures
   *  (unflagged hot pixels and cosmic rays) are dropped, each slice is
   *  averaged, and the slices are averaged with equal weight. A pixel is
   *  left empty unless at least MIN_SLICES of the slices observed it. */
  function matchedImage(groups, npx) {
    const data = [], slice = [];
    groups.forEach((g, k) => { for (const e of g) { data.push(e.data); slice.push(k); } });
    const m = data.length, nSlices = groups.length, need = Math.ceil(nSlices * MIN_SLICES);
    const vals = new Float32Array(m), sorted = new Float32Array(m), dev = new Float32Array(m);
    const sum = new Float64Array(nSlices), cnt = new Uint16Array(nSlices);
    const soft = new Float64Array(nSlices), softCnt = new Uint16Array(nSlices);
    const img = new Float32Array(npx);
    for (let i = 0; i < npx; i++) {
      let c = 0;
      for (let e = 0; e < m; e++) { const v = data[e][i]; vals[e] = v; if (v === v) sorted[c++] = v; }
      let med = 0, limit = Infinity, extreme = Infinity;
      if (c >= 5) {
        const s = sorted.subarray(0, c).sort();
        med = s[c >> 1];
        for (let j = 0; j < c; j++) dev[j] = Math.abs(s[j] - med);
        const sigma = 1.4826 * dev.subarray(0, c).sort()[c >> 1];
        limit = 6 * sigma;
        // A star's centre legitimately swings between exposures as it lands
        // on different parts of a pixel. Only a value far beyond the pixel's
        // own level may remove a slice altogether.
        extreme = 8 * (Math.abs(med) + sigma);
      }
      sum.fill(0); cnt.fill(0); soft.fill(0); softCnt.fill(0);
      for (let e = 0; e < m; e++) {
        const v = vals[e];
        if (v !== v) continue;
        const d = Math.abs(v - med);
        if (d <= limit) { sum[slice[e]] += v; cnt[slice[e]]++; }
        else if (d <= extreme) { soft[slice[e]] += v; softCnt[slice[e]]++; }
      }
      let total = 0, have = 0;
      for (let k = 0; k < nSlices; k++) {
        if (cnt[k]) { total += sum[k] / cnt[k]; have++; }
        else if (softCnt[k]) { total += soft[k] / softCnt[k]; have++; }
      }
      img[i] = have >= need ? total / have : NaN;
    }
    return img;
  }

  function frameOf(ep, used, img, matched) {
    let flags = 0;
    for (const e of used) flags |= e.meta.flags;
    const w = used.map((e) => e.meta.wave);
    return {
      img, n: used.length, flags, matched, epoch: ep.index,
      t0: Math.min(...used.map((e) => e.meta.mjd)), t1: Math.max(...used.map((e) => e.meta.mjd)),
      wmin: Math.min(...w), wmax: Math.max(...w), wave: ep.wave || w.reduce((a, b) => a + b, 0) / w.length,
    };
  }

  /** Measure how the epochs differ in alignment and sharpness, using the
   *  positions and widths of ordinary stars in each epoch's frame. */
  function diagnose(frames, mean, report) {
    const n = state.field.n, scale = state.field.scale_arcsec, R = 4;
    const vals = [];
    for (let i = 0; i < mean.length; i += 5) if (mean[i] === mean[i]) vals.push(mean[i]);
    const med = median(vals), sigma = 1.4826 * median(vals.map((v) => Math.abs(v - med)));
    const stars = [];
    for (let j = R + 4; j < n - R - 4; j++) for (let i = R + 4; i < n - R - 4; i++) {
      const v = mean[j * n + i];
      if (!(v > med + 12 * sigma)) continue;
      let peak = true, flat = 0;
      for (let dy = -R; dy <= R && peak; dy++) for (let dx = -R; dx <= R; dx++) {
        const u = mean[(j + dy) * n + i + dx];
        if (u !== u || u > v) { peak = false; break; }
        if (u === v) flat++;
      }
      if (peak && flat < 3) stars.push({ i, j, v }); // a flat top means a saturated core
    }
    stars.sort((a, b) => b.v - a.v);
    const sample = stars.slice(0, 40);
    report.stars = sample.length;
    if (sample.length < MIN_STARS) return;

    const measure = (img, s) => {
      const edge = [];
      for (let d = -R; d <= R; d++) edge.push(img[(s.j - R) * n + s.i + d], img[(s.j + R) * n + s.i + d], img[(s.j + d) * n + s.i - R], img[(s.j + d) * n + s.i + R]);
      const bkg = median(edge.filter((v) => v === v));
      let sw = 0, sx = 0, sy = 0;
      for (let dy = -R; dy <= R; dy++) for (let dx = -R; dx <= R; dx++) {
        const w = img[(s.j + dy) * n + s.i + dx] - bkg;
        if (!(w > 0)) continue;
        sw += w; sx += w * dx; sy += w * dy;
      }
      if (!(sw > 0)) return null;
      const cx = sx / sw, cy = sy / sw;
      let sr = 0;
      for (let dy = -R; dy <= R; dy++) for (let dx = -R; dx <= R; dx++) {
        const w = img[(s.j + dy) * n + s.i + dx] - bkg;
        if (w > 0) sr += w * ((dx - cx) ** 2 + (dy - cy) ** 2);
      }
      return { x: cx, y: cy, width: Math.sqrt(sr / sw) };
    };
    const table = sample.map((s) => frames.map((f) => measure(f.img, s))).filter((row) => row.every(Boolean));
    if (table.length < MIN_STARS) { report.stars = table.length; return; }
    report.stars = table.length;
    frames.forEach((f, k) => {
      const dx = [], dy = [], w = [];
      for (const row of table) {
        const mx = row.reduce((t, m) => t + m.x, 0) / row.length, my = row.reduce((t, m) => t + m.y, 0) / row.length;
        const mw = row.reduce((t, m) => t + m.width, 0) / row.length;
        dx.push(row[k].x - mx); dy.push(row[k].y - my); w.push(row[k].width / mw);
      }
      const ep = report.epochs.find((e) => e.index === f.epoch);
      // East is to the left, so a positive pixel offset in x is a shift west.
      ep.shift = Math.hypot(median(dx), median(dy)) * scale;
      ep.width = median(w);
    });
    const used = report.epochs.filter((e) => e.used);
    report.maxShift = Math.max(...used.map((e) => e.shift));
    report.widthSpread = Math.max(...used.map((e) => e.width)) - Math.min(...used.map((e) => e.width));
  }

  // ------------------------------------------------------------- panel ---

  const rangeText = () => {
    const b = state.band;
    if (b === "all") return "0.75–5.0 µm";
    if (b.startsWith("c:")) return `${b.slice(2).replace("-", "–")} µm`;
    return `${cfg.bands[b][0]}–${cfg.bands[b][1]} µm (band ${b.slice(1)})`;
  };

  const LIMITS = `
    <details><summary>What can still mislead in this view</summary><ul>
      <li><b>Wavelengths are matched in slices, not exactly.</b> Each slice is a sixth of a detector band. Inside a slice the epochs' wavelengths differ a little, so a source with a steep or spiky spectrum can still show a small false change.</li>
      <li><b>Alignment is not corrected.</b> Images are placed using the archive's astrometry. The measured offsets are listed above; an offset of a fraction of a 6.15″ pixel puts a dark and bright pair on every star in the difference view. A real mover shows that pair on one object only.</li>
      <li><b>Star shapes differ between epochs.</b> SPHEREx's sharpness depends on wavelength and on where a star falls on the detector, and the telescope is rotated between epochs. Rings and spikes around bright stars in the difference view are that, not variability.</li>
      <li><b>Sky level is removed.</b> Each exposure has its median brightness subtracted to take out zodiacal light, so a change in overall sky brightness, or in emission that fills the field, is not preserved.</li>
      <li><b>Masked and saturated pixels.</b> Flagged pixels and isolated outliers are left out. A pixel is shown for an epoch only if at least 85% of the matched slices observed it, so at such pixels the wavelength mix can differ slightly between epochs. Saturated star cores are filled with a flat value and carry no brightness information.</li>
      <li><b>Fast movers.</b> An asteroid moves during one epoch and is smeared or lost in an epoch average. Look for those in Explore wavelengths, one frame per day or exposure.</li>
    </ul></details>`;

  const signed = (pct) => { const r = Math.round(pct); return r === 0 ? "0%" : `${r > 0 ? "+" : "−"}${Math.abs(r)}%`; };

  const fmtCount = (n) => Number(n).toLocaleString("en-US");
  /** Headline figures for how much data stands behind what is on screen. */
  const coverage = (items) => `<dl class="coverage">${items.map(([label, value, of]) =>
    `<div><dt>${label}</dt><dd>${value}${of ? ` <small>${of}</small>` : ""}</dd></div>`).join("")}</dl>`;

  function render() {
    const box = $("context");
    const exploring = state.view === "explore";
    $("view-compare").classList.toggle("on", !exploring);
    $("view-explore").classList.toggle("on", exploring);
    $("view-compare").setAttribute("aria-checked", String(!exploring));
    $("view-explore").setAttribute("aria-checked", String(exploring));
    // The frame grouping only applies when exploring. Keep the control in
    // view when comparing, switched off and with a pointer to where it works.
    $("group").disabled = !exploring;
    $("group-note").hidden = exploring;
    $("mode").options[1].textContent = exploring ? "difference from the average (mixes wavelengths)" : "difference from the average of the matched epochs";
    if (!state.field || !state.exposures.length) { box.hidden = true; return; }
    box.hidden = false;
    const diff = state.mode === "diff";

    if (exploring) {
      const unit = { exposure: "exposure", day: "day", pass: "survey pass" }[state.group];
      box.className = "context warn";
      box.innerHTML = `
        <p class="badge"><span aria-hidden="true">⚠</span> Viewing different infrared wavelengths</p>
        ${coverage([
          ["Exposures loaded", fmtCount(state.field.total), state.field.found > state.field.total ? `of ${fmtCount(state.field.found)}` : ""],
          ["Frames in this view", fmtCount(state.frames.length), ""],
        ])}
        <p>Each frame here is one ${unit} within ${rangeText()}, and each was taken at ${state.group === "exposure" ? "a different wavelength" : "a different mix of wavelengths"}, shown on the image.
        An object can brighten or fade from one frame to the next <b>purely because of its colour</b>, with nothing about it changing.
        Use this view to see what the sky looks like across the infrared, or to follow fast movers day by day.</p>
        ${diff ? "<p><b>The difference view here subtracts images taken at different wavelengths.</b> What is left mixes colour with change over time, so it does not isolate variability.</p>" : ""}
        <p><button type="button" id="go-compare">Compare the same wavelengths over time instead</button></p>`;
      $("go-compare").onclick = () => setView("compare");
      return;
    }

    const r = state.frames.report;
    if (state.loading || !r) {
      box.className = "context";
      box.innerHTML = `<p class="badge">Loading…</p><p>Exposures are still arriving. The frames shown are a preview; wavelengths will be matched between epochs once everything has loaded.</p>`;
      return;
    }
    if (r.status === "none") {
      box.className = "context none";
      box.innerHTML = `
        <p class="badge"><span aria-hidden="true">∅</span> No comparison over time is available</p>
        <p>${r.why} ${r.still != null ? "The image shown is a single epoch, not a comparison." : ""}
        Try a wider wavelength range or a longer date range, or browse the exposures individually.</p>
        <p><button type="button" id="go-explore">Explore wavelengths instead</button></p>`;
      $("go-explore").onclick = () => setView("explore");
      return;
    }

    const used = r.epochs.filter((e) => e.used), dropped = r.epochs.filter((e) => !e.used);
    const f = state.field;
    const haveStars = r.maxShift != null;
    const warn = [];
    if (haveStars && r.maxShift > SHIFT_WARN) warn.push(`the epochs are offset by up to ${r.maxShift.toFixed(1)}″`);
    if (haveStars && r.widthSpread > WIDTH_WARN) warn.push(`star images differ in width by ${(r.widthSpread * 100).toFixed(0)}%`);
    if (r.sharedPixels < 0.8) warn.push(`only ${(r.sharedPixels * 100).toFixed(0)}% of the field was observed in every epoch`);
    // One entry per epoch, stacked so it reads in a narrow column.
    const rows = r.epochs.map((e) => e.used
      ? `<li><b>Epoch ${e.index + 1}</b> <span class="when">${fmtDate(e.t0)} → ${fmtDate(e.t1)}</span>
          <dl><div><dt>Exposures</dt><dd>${e.nUsed} of ${e.exposures.length}</dd></div>
          <div><dt>Wavelength</dt><dd>${e.wave.toFixed(3)} µm</dd></div>
          <div><dt>Offset</dt><dd>${e.shift != null ? e.shift.toFixed(1) + "″" : "–"}</dd></div>
          <div><dt>Star width</dt><dd>${e.width != null ? signed((e.width - 1) * 100) : "–"}</dd></div></dl></li>`
      : `<li class="out"><b>Epoch ${e.index + 1}</b> <span class="when">${fmtDate(e.t0)} → ${fmtDate(e.t1)}</span>
          <p>Left out: ${e.reason}.</p></li>`).join("");
    box.className = "context ok";
    box.innerHTML = `
      <p class="badge"><span aria-hidden="true">✓</span> Comparing the same spectral region over time</p>
      ${coverage([
        ["Exposures loaded", fmtCount(f.total), f.found > f.total ? `of ${fmtCount(f.found)}` : ""],
        ["Epochs compared", used.length, `of ${r.epochs.length}`],
        ["Wavelength slices", r.slices, `of ${r.slicesTotal}`],
        ["Field in every epoch", `${(r.sharedPixels * 100).toFixed(0)}%`, ""],
      ])}
      <p>Each frame is one epoch, built from the <b>same ${r.slices} wavelength slice${r.slices > 1 ? "s" : ""}</b> within ${rangeText()}, weighted equally, so a difference between frames is not a difference in wavelength${diff ? "; the difference view subtracts the average of these matched epochs" : ""}.
      Small differences can still come from how well the epochs line up and how sharp their stars are, ${haveStars ? "which are measured below but not corrected" : "which could not be measured in this field"}.
      ${r.slices < r.slicesTotal ? `${r.slicesTotal - r.slices} of ${r.slicesTotal} slices are not used because at least one epoch did not observe them.` : ""}
      ${r.noisy ? `${r.noisy} unusually noisy exposure${r.noisy > 1 ? "s are" : " is"} left out.` : ""}
      ${dropped.length ? `${dropped.length} epoch${dropped.length > 1 ? "s are" : " is"} left out for lack of matching wavelengths.` : ""}</p>
      ${r.epochs.length > 4
        ? `<details class="epochs-more"><summary>Details of all ${r.epochs.length} epochs (${used.length} compared, ${dropped.length} left out)</summary><ol class="epochs">${rows}</ol></details>`
        : `<ol class="epochs">${rows}</ol>`}
      <p class="fine">Mean wavelengths agree to ${r.waveSpread < 0.001 ? "better than 0.1" : "within " + (r.waveSpread * 100).toFixed(1)}%. ${(r.sharedPixels * 100).toFixed(0)}% of the field was observed in every epoch.
      ${haveStars ? `Offset and star width are measured from ${r.stars} stars, relative to the average of the epochs.` : "Too few suitable stars in this field to measure alignment or sharpness."}
      All frames share one brightness scale.</p>
      ${warn.length ? `<p class="caution"><span aria-hidden="true">⚠</span> Take care: ${warn.join("; ")}. Expect residue around stars in the difference view.</p>` : ""}
      ${LIMITS}`;
  }

  return { build, render };
})();
window.Compare = Compare;
