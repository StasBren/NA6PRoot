// ===================================================================================================== SIMPLE MAIN DECK
function bigStat(s, x, y, w, value, label, color) {
  s.addText(value, { x, y, w, h: 0.62, fontSize: 36, bold: true, color: color || INK, margin: 0, isTextBox: true, valign: "top" });
  s.addText(eqRuns(label, { fontSize: 14, color: INK2 }), { x, y: y + 0.66, w, h: 0.62, margin: 0, isTextBox: true, valign: "top" });
}
function flow(s, items, y, h, opt = {}) {
  const n = items.length, gap = opt.gap || 0.42, w = (W - gap * (n - 1)) / n;
  items.forEach(([head, body, hl], i) => {
    const x = X0 + i * (w + gap);
    s.addShape(pres.shapes.RECTANGLE, { x, y, w, h, fill: { color: hl ? TINT : "FFFFFF" }, line: { color: hl === "red" ? RED : hl ? INK : RULE, width: 1.25 } });
    s.addText(head, { x: x + 0.15, y: y + 0.15, w: w - 0.3, h: 0.36, fontSize: 16, bold: true, color: hl === "red" ? RED : INK, margin: 0, isTextBox: true });
    s.addText(eqRuns(body, { fontSize: 13, color: INK2 }), { x: x + 0.15, y: y + 0.58, w: w - 0.3, h: h - 0.7, margin: 0, isTextBox: true, valign: "top" });
    if (i < n - 1) arrow(s, x + w + 0.04, y + h / 2, x + w + gap - 0.04, y + h / 2, INK);
  });
}

// 1 title
pres.addSection({ title: "Main" });
{
  const s = pres.addSlide({ masterName: "DARK", sectionTitle: "Main" });
  s.addShape(pres.shapes.RECTANGLE, { x: X0, y: 2.181, w: 1.111, h: 0.042, fill: { color: RED }, line: { color: RED, width: 0 } });
  s.addText([{ text: "Fast cathode response", options: { breakLine: true } }, { text: "for the NA60+ MWPC" }], { x: X0, y: 2.45, w: 11.9, h: 1.45, fontSize: 46, bold: true, color: DTXT, margin: 0, isTextBox: true, valign: "top" });
  s.addText("Can ten numbers per ion cloud replace ion-by-ion signal calculations?", { x: X0, y: 4.0, w: 11.9, h: 0.45, fontSize: 24, color: DTXT, margin: 0, isTextBox: true });
  s.addText("2100 V · B = 0 · Ar/CO₂ 70:30 · idealized 9-anode cell · 192 Garfield++ muons", { x: X0, y: 4.62, w: 11.9, h: 0.365, fontSize: 18, color: DSUB, margin: 0, isTextBox: true });
  s.addText("Development validation on non-sealed events; the sealed holdout is still closed · November 2026", { x: X0, y: 6.367, w: 11.9, h: 0.272, fontSize: 13, color: MUTED, margin: 0, isTextBox: true });
  N(s, "This talk continues the avalanche study. That deck measured what the ion clouds look like; its open item was the cathode signal, because the archives contain ion births but no strip currents.\n\nThe question today is whether we can compute the U/V strip signal from a few numbers per ion cloud instead of following every ion. Everything is a development test at one setup, 2100 V without magnetic field; the sealed holdout events have not been opened.");
}
// 2 problem
{
  const s = content("THE PROBLEM", "Slow ions make the signal, and following them is expensive", "Main");
  img(s, "signal_anatomy", X0, 1.86, 7.4, 4.5);
  bigStat(s, 8.85, 2.05, 3.6, "≈ 85 s", "per event to follow every ion through the field (detailed reference, S0)");
  bigStat(s, 8.85, 3.65, 3.6, "40–100 µs", "for the ions to reach the cathode: the signal is long and slow");
  N(s, "The plot shows one strip signal from the detailed calculation, which we call S0. When an avalanche happens, the positive ions move from the anode wire to the cathode, and the moving charge induces the signal on the cathode strips. Left: the fast current in the first microsecond, while the ions leave the wire. Right: the charge keeps growing until the ions arrive, after 40 to 100 µs.\n\nS0 follows every sampled ion along its path and evaluates the weighting field of every strip along the way. That takes about 85 seconds per muon on one CPU. For millions of events this is far too slow, which is why we need a fast model.");
}
// 3 idea
{
  const s = content("THE IDEA", "K2A summarises each ion cloud and looks up its signal in a table", "Main");
  flow(s, [["Ion cloud", "all ions of one primary electron on one anode"],
           ["10 numbers", "charge, time, radius, strip position, angular shape (2 harmonics), up/down asymmetry", true],
           ["Response table", "built once from the same ion paths and weighting field as S0", true],
           ["Strip signals", "U and V currents on absolute strips"]], 2.05, 1.95);
  note(s, X0, 4.45, 11.5, 1.6, ["Why these ten numbers: the avalanche study showed that clouds are short in time and small along the wire, that ions are born ≈ 17.5 µm from the anode, and that most clouds sit on one side of it.",
    "Whether they are enough is measured, not assumed: the largest approximation is using one radius per cloud (backup B1)."], 15);
  N(s, "The idea of K2A is simple. Instead of following every ion, we describe each ion cloud, all ions from one primary electron on one wire, by ten numbers: its charge, its time, an effective radius, its position in strip coordinates, two angular harmonics and the up/down asymmetry.\n\nThe response of a cloud is then looked up in a table that is computed once, using exactly the same ion paths and the same weighting field as the detailed calculation. Summing all clouds gives the U and V strip currents.\n\nThe choice of the ten numbers comes from the avalanche study: avalanches are short and small, ions are born very close to the wire, and most clouds are strongly up/down polarised, which is what the histogram of A shows. Whether these numbers are enough is a measured question; the largest single approximation turns out to be using one radius per cloud.");
}
// 4 test
{
  const s = content("THE TEST", "Same avalanches, two calculations, criteria fixed in advance", "Main");
  flow(s, [["Garfield++ avalanche", "archived ion births of one muon (no regeneration)"],
           ["S0: every ion", "detailed reference: drift and weighting field ion by ion"],
           ["K2A: 10 numbers", "fast prediction per cloud, uses nothing from S0", "red"],
           ["Compare", "waveforms, charge, position, timing, threshold decisions"]], 2.05, 1.85);
  bigStat(s, X0, 4.35, 3.6, "192", "independent muons, all stations and channels");
  bigStat(s, 4.85, 4.35, 3.6, "9", "acceptance criteria, frozen before the comparison (27 Nov)");
  bigStat(s, 8.85, 4.35, 3.6, "20", "sealed holdout events, not opened", RED);
  N(s, "To test K2A we use real archived Garfield++ avalanches and compute the strip signal twice: once with the detailed reference S0, once with K2A. Both start from the same ion births, so the only difference is the model.\n\nWe compare what reconstruction would use: the waveforms after a shaping stand-in, the charge per plane, the U/V charge ratio, the reconstructed position from the strip centroids, the timing at threshold, and which strips are above threshold.\n\nNine acceptance criteria were written down and frozen on 27 November, before the comparisons. This talk uses 192 independent muons. A separate set of 20 sealed events is kept untouched for the final test.");
}
// 5 waveforms
{
  const s = content("RESULT · SIGNALS", "The signals agree strip by strip", "Main");
  img(s, "waves_typical", X0, 1.84, 11.556, 4.45, "center");
  note(s, X0, 6.36, 11.5, 0.42, ["A typical event: black S0, red dashed K2A. Difference below 0.3 % of the peak; K2A is slightly late on the rising edge."], 12);
  N(s, "Here is a typical event: the three largest strips in each plane after shaping. Black is the detailed calculation, red dashed is K2A; on this scale you cannot separate them.\n\nThe small panels show the difference relative to the peak. It is below 0.3 % everywhere. There is a small systematic dip on the rising edge: K2A peaks about 0.1 ns later and 0.14 % lower. This comes from describing each cloud with one radius, and it is the same in all events.");
}
// 6 criteria
{
  const s = content("RESULT · CRITERIA", "All 192 events pass; charge is the tightest criterion", "Main");
  img(s, "criteria_margin", X0, 1.86, 7.6, 4.5);
  bigStat(s, 8.85, 1.95, 3.6, "192 / 192", "events pass all frozen criteria");
  bigStat(s, 8.85, 3.3, 3.6, "0.18 %", "typical waveform difference after shaping");
  bigStat(s, 8.85, 4.65, 3.6, "99 %", "of the charge limit used by the worst event: no margin left", RED);
  N(s, "This plot shows, for every criterion, how much of its allowed limit each event uses. Zero means perfect agreement; 100 % means right at the limit.\n\nAll 192 events pass every criterion. Position and timing use only a few percent of their limits, the waveform difference about 20 %.\n\nThe weak point is charge: the charge in the comparison window differs by up to 0.5 %, which is exactly the limit, and the U/V charge ratio uses up to 93 %. They pass, but with no margin, and we do not yet know why. That is the first thing to fix.");
}
// 7 speed
{
  const s = content("RESULT · SPEED", "K2A computes the signal about 50,000 times faster", "Main");
  img(s, "runtime", X0, 1.84, 11.556, 4.25, "center");
  note(s, X0, 6.2, 11.5, 0.55, ["S0 ≈ 85 s per event; K2A ≈ 1.6 ms to evaluate (≈ 12 ms including reading the file). The avalanche itself is not included in either."], 12);
  N(s, "Left: time per event. The detailed calculation takes about 85 seconds; K2A needs 7 ms to read the file, 2.5 ms to build the ten numbers and 1.6 ms to compute the signal.\n\nRight: the speedup for each event. Computing the signal alone is about 54,000 times faster; including file reading and building the numbers it is about 7,600 times faster.\n\nThese numbers are only for the signal calculation. Generating the avalanche is a separate cost, handled later by the fast avalanche model Fast4.");
}
// 8 limits
{
  const s = content("LIMITS", "What is not shown yet", "Main");
  img(s, "selection_bias", X0, 1.84, 6.4, 3.95);
  rows(s, 7.6, 1.95, 4.85, 1.7, [
    ["Small events", "tested up to 26 clouds; 38 % of events are that small"],
    ["First 7 µs", "the table ends there: ≈ 15 % of the final charge"],
    ["One setup", "2100 V, B = 0, ideal cell; other conditions need new tables"],
    ["Charge", "passes, but at the limit; cause unknown"],
    ["Readout", "stand-in shaper, not the real VMM3a"],
  ], { rh: 0.78, fs: 13, lfs: 14 });
  N(s, "Where the result does not yet apply.\n\nThe plot compares the tested events, red, with all events, grey. The test queue runs small events first, so we have only shown it for events with up to 26 clouds, about 38 % of the sample. Two much larger events also pass, but two is not a statistic.\n\nThe response table stops at 8 µs, and in the first 7 µs only about 15 % of the final charge has arrived; the slow tail is not modelled yet. Everything is at one setup, 2100 V without magnetic field. The charge criterion is at its limit. And the shaping is a placeholder for the real electronics.");
}
// 9 summary + next
{
  const s = content("SUMMARY", "Yes, within the tested domain; next, fix charge and widen the domain", "Main");
  s.addText("Shown", { x: X0, y: 1.9, w: 5.5, h: 0.4, fontSize: 18, bold: true, color: INK, margin: 0, isTextBox: true });
  s.addText("Next", { x: 6.95, y: 1.9, w: 5.5, h: 0.4, fontSize: 18, bold: true, color: RED, margin: 0, isTextBox: true });
  const shown = ["Ten numbers per cloud reproduce the detailed signal in all 192 tested events", "Waveforms agree to ≈ 0.2 %, positions to ≈ 1 µm, timing to 0.1 ns",
                 "The signal calculation is ≈ 50,000 times faster", "The largest approximation is known: one radius per cloud"];
  const next = ["Find the cause of the charge difference", "Extend the table: more wires, larger radii, longer times",
                "Test larger events, then the 20 sealed events once", "Connect the fast avalanche model (Fast4) and real electronics"];
  [[shown, X0, INK], [next, 6.95, RED]].forEach(([L, x, c]) => {
    s.addShape(pres.shapes.LINE, { x, y: 2.38, w: 5.5, h: 0, line: { color: c, width: 1.25 } });
    s.addText(L.map((t, i) => ({ text: t, options: { bullet: true, breakLine: i < L.length - 1 } })),
      { x, y: 2.55, w: 5.5, h: 3.6, fontSize: 16, color: INK3, margin: 0, isTextBox: true, valign: "top", paraSpaceAfter: 12 });
  });
  N(s, "To answer the question of the talk: yes, within what we tested. Ten numbers per ion cloud and a precomputed table reproduce the detailed cathode signal in all 192 events: waveforms to about 0.2 %, positions to about a micrometre, timing to a tenth of a nanosecond, at about 50,000 times the speed. We also know which approximation dominates.\n\nNext: first understand the charge difference, because that criterion has no margin. Then extend the table to more wires, larger radii and longer times, test larger events, and run the sealed events once. The last step is connecting the fast avalanche model and the real electronics response.");
}
