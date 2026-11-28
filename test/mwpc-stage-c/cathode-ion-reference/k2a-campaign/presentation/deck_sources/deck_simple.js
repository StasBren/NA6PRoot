// K2A signal-response deck, NA60+ MWPC. Visual system of the avalanche deck (Public Sans, dark title, off-white content).
const pptxgen = require("pptxgenjs");
const fs = require("fs");
const { applyTheme } = require("/root/.claude/skills/synced/40fa6195-ed1c-488b-8003-1ed866607756_7c481322-ac30-49de-9e37-048a12729e85/pptx/scripts/apply_theme.js");

const FIG = "/home/claude/pres/fig/";
const SZ = JSON.parse(fs.readFileSync("/home/claude/pres/an/figsizes.json"));
const THEME = {
  name: "NA60 K2A", headFontFace: "Public Sans", bodyFontFace: "Public Sans",
  colors: { dk1: "15191F", lt1: "FBFBF8", dk2: "3A414A", lt2: "F1F2EF", accent1: "C8102E", accent2: "1F4E9E",
            accent3: "D9822B", accent4: "2E7D38", accent5: "8A9099", accent6: "4A525C", hlink: "1F4E9E", folHlink: "4A525C" },
};
const INK = "15191F", INK2 = "4A525C", INK3 = "3A414A", MUTED = "8A9099", RULE = "C9CED4", RED = "C8102E", LRED = "E57A8A",
  BG = "FBFBF8", DARK = "15191F", DTXT = "F1F2EF", DSUB = "C9CED4", TINT = "F1F2EF";
const X0 = 0.889, W = 11.556, FOOT = "NA60+ MWPC · K2A fast cathode response · 2100 V, B = 0, idealized 9-anode cell";

const pres = new pptxgen();
pres.layout = "LAYOUT_WIDE";
pres.theme = { headFontFace: "Public Sans", bodyFontFace: "Public Sans" };
pres.title = "K2A fast cathode response for the NA60+ MWPC";
pres.author = "NA60+ MWPC simulation";

pres.defineSlideMaster({
  title: "DARK", background: { color: DARK },
  objects: [],
});
pres.defineSlideMaster({
  title: "CONTENT", background: { color: BG },
  objects: [
    { text: { text: FOOT, options: { x: X0, y: 6.856, w: 9.0, h: 0.228, fontFace: "Public Sans", fontSize: 12, color: MUTED, margin: 0 } } },
    { placeholder: { options: { name: "title", type: "title", x: X0, y: 1.172, w: 11.902, h: 0.46, fontFace: "Public Sans", fontSize: 26,
                                bold: true, color: INK, margin: 0, valign: "top", align: "left" }, text: "" } },
  ],
  slideNumber: { x: 11.95, y: 6.856, w: 0.5, h: 0.228, fontFace: "Public Sans", fontSize: 12, color: MUTED, align: "right", margin: 0 },
});

// ---------- helpers
function eqRuns(s, base) {
  // "K_{0}" -> subscript, "x^{2}" -> superscript, **bold** spans
  const runs = []; const re = /(_\{[^}]*\}|\^\{[^}]*\}|\*\*[^*]+\*\*)/g; let last = 0, m;
  while ((m = re.exec(s)) !== null) {
    if (m.index > last) runs.push({ text: s.slice(last, m.index), options: { ...base } });
    const t = m[0];
    if (t.startsWith("_")) runs.push({ text: t.slice(2, -1), options: { ...base, subscript: true } });
    else if (t.startsWith("^")) runs.push({ text: t.slice(2, -1), options: { ...base, superscript: true } });
    else runs.push({ text: t.slice(2, -2), options: { ...base, bold: true } });
    last = m.index + t.length;
  }
  if (last < s.length) runs.push({ text: s.slice(last), options: { ...base } });
  return runs;
}
function para(lines, base) { // array of strings -> runs with breaks
  const out = [];
  lines.forEach((l, i) => {
    const r = eqRuns(l, base);
    if (i < lines.length - 1) r[r.length - 1].options.breakLine = true;
    out.push(...r);
  });
  return out;
}
function kicker(s, text, dark) {
  s.addText(text, { x: X0, y: 0.889, w: 11.902, h: 0.228, fontSize: 12, bold: true, color: dark ? LRED : RED, margin: 0, charSpacing: 1, isTextBox: true, objectName: "kicker" });
}
let APPX = false, BN = 0;
function content(kick, title, sec) {
  if (APPX) sec = "Backup";
  const s = pres.addSlide({ masterName: "CONTENT", sectionTitle: sec });
  if (APPX && !/^A\d/.test(kick)) { BN += 1; kick = `B${BN} · ` + kick.replace(/^\d+ · /, ""); }
  kicker(s, kick);
  s.addText(title, { placeholder: "title" });
  return s;
}
function img(s, name, x, y, maxW, maxH, align) {
  const [w0, h0] = SZ[name]; const k = Math.min(maxW / w0, maxH / h0); const w = w0 * k, h = h0 * k;
  const xx = align === "right" ? x + maxW - w : align === "center" ? x + (maxW - w) / 2 : x;
  s.addImage({ path: FIG + name + ".png", x: xx, y, w, h, altText: name.replace(/_/g, " "), objectName: name });
  return { x: xx, y, w, h };
}
function stat(s, x, y, w, value, label, opts = {}) {
  s.addShape(pres.shapes.LINE, { x, y, w, h: 0, line: { color: opts.rule || INK, width: 1 } });
  s.addText(value, { x, y: y + 0.07, w, h: 0.42, fontSize: opts.vsize || 24, bold: true, color: opts.color || INK, margin: 0, isTextBox: true, valign: "top" });
  s.addText(eqRuns(label, { fontSize: 12, color: INK2 }), { x, y: y + 0.5, w, h: opts.lh || 0.46, margin: 0, isTextBox: true, valign: "top" });
}
function note(s, x, y, w, h, lines, size) {
  s.addShape(pres.shapes.LINE, { x: x + 0.014, y, w: 0, h, line: { color: RED, width: 1.25 } });
  s.addText(para(lines, { fontSize: size || 12, color: INK2 }), { x: x + 0.164, y, w: w - 0.164, h, margin: 0, isTextBox: true, valign: "top", paraSpaceAfter: 2 });
}
function rows(s, x, y, w, labW, items, opt = {}) {
  // items: [label, text] ; horizontal rules between
  const rh = opt.rh || 0.62, fs = opt.fs || 13;
  items.forEach(([lab, txt], i) => {
    const yy = y + i * rh;
    s.addShape(pres.shapes.LINE, { x, y: yy, w, h: 0, line: { color: RULE, width: 0.75 } });
    s.addText(eqRuns(lab, { fontSize: opt.lfs || 14, bold: true, color: opt.lcol || INK }), { x, y: yy + 0.08, w: labW - 0.1, h: rh - 0.12, margin: 0, isTextBox: true, valign: "top" });
    s.addText(eqRuns(txt, { fontSize: fs, color: INK3 }), { x: x + labW, y: yy + 0.08, w: w - labW, h: rh - 0.12, margin: 0, isTextBox: true, valign: "top" });
  });
}
function box(s, x, y, w, h, head, body, opt = {}) {
  s.addShape(pres.shapes.RECTANGLE, { x, y, w, h, fill: { color: opt.fill || "FFFFFF" }, line: { color: opt.line || RULE, width: 1 } });
  s.addText(eqRuns(head, { fontSize: opt.hfs || 12, bold: true, color: opt.hcol || INK }), { x: x + 0.12, y: y + 0.08, w: w - 0.24, h: 0.26, margin: 0, isTextBox: true, valign: "top" });
  if (body) s.addText(para(body, { fontSize: opt.bfs || 11, color: INK2 }), { x: x + 0.12, y: y + 0.36, w: w - 0.24, h: h - 0.42, margin: 0, isTextBox: true, valign: "top" });
}
function arrow(s, x1, y1, x2, y2, col) {
  s.addShape(pres.shapes.LINE, { x: Math.min(x1, x2), y: Math.min(y1, y2), w: Math.abs(x2 - x1), h: Math.abs(y2 - y1),
    flipH: x2 < x1, flipV: y2 < y1, line: { color: col || INK2, width: 1.25, endArrowType: "triangle" } });
}
const N = (s, t) => s.addNotes(t);

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
    "Whether they are enough is measured, not assumed: the largest approximation is using one radius per cloud (backup B2)."], 15);
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

// ===================================================================================================== BACKUP
pres.addSection({ title: "Backup" });
APPX = true;
{
  const s = pres.addSlide({ masterName: "DARK", sectionTitle: "Backup" });
  kicker(s, "BACKUP", true);
  s.addText("Technical detail for questions", { x: X0, y: 1.172, w: 11.6, h: 0.6, fontSize: 26, bold: true, color: DTXT, margin: 0, isTextBox: true });
  const items = ["B1  The ten cloud quantities and their motivation", "B2  K2A construction and the measured approximation budget", "B3  Workflow and what depends on detector conditions",
                 "B4  Validation design and execution gates", "B5  Error distributions and complexity", "B6  Detector observables", "B7  Validity domain in detail",
                 "A1–A11  Criteria table, definitions, response table, reproducibility, thresholds, time collapse, worst event, residual and horizon, runtime, provenance, domain"];
  s.addText(items.map((t, i) => ({ text: t, options: { breakLine: i < items.length - 1 } })), { x: X0, y: 2.2, w: 11.3, h: 4.3, fontSize: 16, color: DSUB, margin: 0, isTextBox: true, paraSpaceAfter: 8, valign: "top" });
  N(s, "Backup slides: the full technical version of the talk.");
}
// =====================================================================================================  6 moments
{
  const s = content("02 · MODEL", "Ten numbers per ion cloud, motivated by the avalanche study", "02 Model");
  img(s, "moments", X0, 1.86, 5.4, 4.2);
  rows(s, 6.55, 1.86, 5.9, 1.65, [
    ["N", "cloud charge (fC); gain is Pólya-like, so N varies strongly"],
    ["t", "T1: the primary's weighted-mean birth time (σ_{t} of one avalanche ≈ 66 ps)"],
    ["r_{eff}", "exp⟨ln r⟩ about the anode; here 17.53 ± 0.17 µm"],
    ["D_{U}, D_{V}", "strip coordinates of the cloud (avalanche ≈ 12 µm along the wire)"],
    ["C_{1}, S_{1}, C_{2}, S_{2}", "⟨cos nθ⟩, ⟨sin nθ⟩, n ≤ 2: angular spread 10–16°"],
    ["A", "⟨sign sin θ⟩: up/down polarisation, RMS 0.88 here (0.87 per avalanche)"],
  ], { rh: 0.57, fs: 12, lfs: 13 });
  note(s, 6.55, 5.38, 5.9, 0.95, ["Cloud = (Heed primary, nearest anode). All moments are charge-weighted over the K256 sample. Stored but not used by the evaluator: t_{rms} (zero by construction after the T1 collapse), raw time RMS, clearance, row count. Derived at evaluation: hemisphere coefficient β, per-primary strip recentring, stereo term."], 10.5);
  N(s, "A cloud is one Heed primary on one anode. For each cloud K2A stores ten numbers: the charge, one time, one effective radius, two strip coordinates, the first two angular harmonics, and the up/down asymmetry A. All are charge-weighted means over the archived births.\n\nThe avalanche study is the motivation for this choice. A single avalanche is short in time, about 66 ps, and small along the wire, about 12 µm, so one time and one strip coordinate per cloud are natural. The ions are born within about 20 µm of the anode: r_eff is 17.5 µm with an RMS of only 0.17 µm across 3,214 clouds. The angular spread of 10 to 16 degrees suggests that low harmonics carry most of the angular information, and the strong up/down polarisation, RMS of A of 0.88 here, is why a separate hemisphere term exists.\n\nThis is motivation, not proof. Whether two harmonics plus the hemisphere term, one radius and one D are enough is decided by the measured comparisons on the next slide and in the accuracy section.\n\nThe moment file also contains diagnostic columns that the evaluator never reads, among them t_rms, which is zero by construction after the time collapse. That column is relevant later for events 0790 and 0877.");
}
// =====================================================================================================  7 construction
{
  const s = content("02 · MODEL", "Exact at the table nodes; each approximation is measured separately", "02 Model");
  const L = [
    "**Green table** (frozen RKF paths and LUT, as in S0):",
    "   G(k, r, θ, D, τ) = φ_{s}(x(t_{0}+τ)) − φ_{s}(x(t_{0}))",
    "   3 anodes × 14 radii × 72 angles; D step 0.17 mm; τ to 8 µs",
    "**Harmonic kernels** (θ_{m}, m = 1…72):",
    "   K_{0} = Σ G_{m}/M,   Kc_{n} = 2Σ cos(nθ_{m}) G_{m}/M,",
    "   Ks_{n} = 2Σ sin(nθ_{m}) G_{m}/M,   K_{A} = Σ s⊥_{m} G_{m}",
    "**Per cloud:**",
    "   Q_{s}(τ) = N [K_{0} + Σ_{n≤2} (C_{n} Kc_{n} + S_{n} Ks_{n}) + β K_{A}]",
    "   at (r_{eff}, D): 4 corners in (ln r, D), linear in τ at the 10 ns bin edges",
  ];
  s.addText(para(L, { fontSize: 12, color: INK3 }), { x: X0, y: 1.88, w: 5.9, h: 2.8, margin: 0, isTextBox: true, valign: "top", paraSpaceAfter: 2 });
  rows(s, X0, 4.72, 5.9, 1.25, [
    ["Exact", "basis sum ≡ sum over the 72 θ nodes; table nodes reproduce S0 (float32 storage)"],
    ["Approximate", "n ≤ 2 + hemisphere; one r_{eff}, one D, one time per cloud; θ, r, D, τ interpolation"],
  ], { rh: 0.62, fs: 11.5, lfs: 12.5 });
  img(s, "closure_budget", 6.95, 1.86, 5.5, 4.1, "right");
  note(s, 6.95, 6.08, 5.5, 0.65, ["Measured on 5 reference muons (S-ladder): one radius per cloud dominates; n ≤ 2 with the hemisphere term costs 0.09 %."], 11);
  N(s, "The construction has a numerical part and an approximation part, and I want to keep them apart.\n\nNumerically, the Green table G is the change of the weighting potential along a frozen RKF ion path started from each table node: three anodes, 14 radii, 72 angles, D in 0.17 mm steps and τ up to 8 µs, with the same field model and the same LUT as S0. At a node the table reproduces S0 up to float32 storage. Projecting the 72 angles onto harmonics is exact linear algebra: the basis sum equals a weighted sum over the angle nodes.\n\nThe approximations are physical: keeping only harmonics up to n = 2 plus the hemisphere term, giving each cloud one radius, one D and one time, and interpolating between nodes. The truncated harmonic representation is not exact for an arbitrary cloud.\n\nThe bar chart measures each approximation on the five reference muons, unshaped ion term, 0 to 7 µs. Table nodes plus interpolation cost about 0.24 %. One D per cloud is negligible. Angular truncation to n ≤ 2 costs 0.14 %, and 0.09 % with the hemisphere term. The single radius r_eff is by far the largest single item at about 0.52 %, and the full K2A comes out at 0.50 % against S0. This is the comparison evidence that n ≤ 2 plus the hemisphere term is sufficient at this level; the avalanche angular spread only motivated it.");
}
// =====================================================================================================  8 workflow
{
  const s = content("02 · MODEL", "One algorithm; voltage, gas, field and geometry enter through tables", "02 Model");
  const yb = 1.95, bw = 2.05, gap = 0.31;
  const chain = [["Heed primaries", ["Geant4 / Heed", "(unchanged)"], "FFFFFF", RULE],
                 ["Avalanche model", ["Fast4: GOAL", "calibrated at 1800 V"], "FFFFFF", RED],
                 ["Cloud moments", ["10 numbers per cloud", "full births or K256"], "FFFFFF", INK],
                 ["K2A evaluator", ["kernel r14t72edge5", "absolute U/V strips"], TINT, INK],
                 ["Readout", ["shaping, thresholds", "provisional stand-in"], "FFFFFF", RULE]];
  chain.forEach(([h, b, f, l], i) => {
    const x = X0 + i * (bw + gap);
    box(s, x, yb, bw, 1.0, h, b, { fill: f, line: l, bfs: 10.5, hcol: i === 1 ? RED : INK });
    if (i < chain.length - 1) arrow(s, x + bw, yb + 0.5, x + bw + gap, yb + 0.5);
  });
  s.addText(eqRuns("Validated today: Garfield++ K256 archive → moments (p3_birth_moments) → compiled K2A. The Fast4 → K2A interface is not yet built or validated.", { fontSize: 12, color: INK2 }),
    { x: X0, y: 3.05, w: 11.5, h: 0.3, margin: 0, isTextBox: true });
  const hdr = { bold: true, color: INK, fontSize: 11.5, fill: { color: TINT } };
  const c = (t, o = {}) => ({ text: t, options: { fontSize: 11, color: INK3, ...o } });
  s.addTable([
    [c("Layer", hdr), c("Content", hdr), c("Depends on", hdr), c("Change of conditions", hdr)],
    [c("Avalanche (Fast4)", { bold: true }), c("LUT500 transport, avalanche library"), c("voltage, gas, B, cell geometry"), c("recalibrate (needed for 2100 V)")],
    [c("Ion transport", { bold: true }), c("RKF Green grid paths"), c("voltage, gas (mobility), geometry"), c("re-run grid paths, re-tabulate")],
    [c("Weighting field", { bold: true }), c("spectral U/V LUT"), c("strip pitch, stereo, gap"), c("regenerate LUT")],
    [c("Response kernel", { bold: true }), c("K(basis, side, wire, r, D, τ)"), c("the two tables above"), c("rebuild; ranges are data")],
    [c("Algorithm", { bold: true }), c("moments → kernel sum"), c("—"), c("unchanged, except: asymmetric gaps, non-contiguous strips, field wires")],
  ], { x: X0, y: 3.55, w: 11.556, colW: [2.1, 2.9, 2.9, 3.656], rowH: 0.36, border: { type: "solid", color: RULE, pt: 0.75 }, fill: { color: "FFFFFF" }, margin: 0.05 });
  note(s, X0, 5.95, 11.5, 0.75, ["Only one configuration is validated: 2100 V, B = 0, Ar/CO₂ 70:30, the idealized nine-anode cell. A magnetic field acts mainly on electron drift and the avalanche position (Fast4's domain); its effect on the ion paths is an estimate, not a validated statement."], 11);
  N(s, "This is the intended chain. Heed primaries go into a fast avalanche model, Fast4, which should produce the cloud moments directly; K2A turns them into U/V strip currents, and the readout model shapes and digitises them.\n\nWhat is validated today is shorter: the moments are built from the archived Garfield++ births, and the compiled K2A evaluates them. Fast4 exists but was calibrated at 1800 V, and the interface to K2A has not been built or validated. So the Fast4 box is a goal, not a result.\n\nThe table is the portability argument. Voltage, gas and geometry enter through tables: the ion-path grid, the weighting-field LUT and the kernel built from them. Changing conditions means regenerating tables, not rewriting the algorithm, with three exceptions where the representation itself breaks: asymmetric gaps, non-contiguous strips and field wires.\n\nNone of the accuracy numbers in this talk transfer automatically to another voltage, gas, field or geometry; each new table set needs its own validation.");
}
// =====================================================================================================  9 validation design
{
  const s = content("03 · ACCURACY", "Paired comparison on identical avalanches; criteria frozen in advance", "03 Accuracy");
  const fl = [["980 non-sealed muons", "977 with preflight moments; 20 sealed events excluded"],
              ["195 processed slots", "5 pilot + 6 extension + 184 from the 800-muon queue"],
              ["192 accuracy set", "VALID execution, independent of the kernel build"]];
  fl.forEach(([h, b], i) => {
    const y = 1.88 + i * 0.95;
    box(s, X0, y, 3.6, 0.78, h, [b], { bfs: 10.5, fill: i === 2 ? TINT : "FFFFFF", line: i === 2 ? INK : RULE });
    if (i < 2) arrow(s, X0 + 1.8, y + 0.78, X0 + 1.8, y + 0.95);
  });
  box(s, X0, 4.8, 3.6, 1.12, "Not in the accuracy set", ["0782: intended out-of-domain, refused", "0790, 0877: execution FAIL on the", "t_{rms} diagnostic (appendix A4)"], { bfs: 10.5, line: RED, hcol: RED });
  s.addText(eqRuns("Sealed 20-event holdout: not opened", { fontSize: 11, bold: true, color: RED }), { x: X0, y: 6.05, w: 3.6, h: 0.3, margin: 0, isTextBox: true });
  const cx = 4.85, cw = 3.6;
  const refs = [["S0_T1", "detailed, primary-mean birth times: decisive reference", INK],
                ["S0_T0", "detailed, true per-ion birth times", INK2],
                ["K2A", "prediction: own production grouping, nothing from S0", RED],
                ["K2Aimp", "diagnostic: S0 grouping imposed", MUTED]];
  s.addText("Signals compared", { x: cx, y: 1.88, w: cw, h: 0.28, fontSize: 13, bold: true, color: INK, margin: 0, isTextBox: true });
  refs.forEach(([h, b, col], i) => {
    const y = 2.25 + i * 0.62;
    s.addShape(pres.shapes.LINE, { x: cx, y, w: cw, h: 0, line: { color: RULE, width: 0.75 } });
    s.addText(h, { x: cx, y: y + 0.07, w: 1.0, h: 0.4, fontSize: 13, bold: true, color: col, margin: 0, isTextBox: true, valign: "top" });
    s.addText(b, { x: cx + 1.05, y: y + 0.07, w: cw - 1.05, h: 0.5, fontSize: 11, color: INK3, margin: 0, isTextBox: true, valign: "top" });
  });
  s.addText(para(["Components: **model** = S0_T1 vs K2A; **total** = S0_T0 vs K2A; **time collapse** = S0_T0 vs S0_T1.",
                  "Ion term + electron step is decisive; ion-only is reported alongside."], { fontSize: 11, color: INK2 }),
    { x: cx, y: 4.85, w: cw, h: 1.1, margin: 0, isTextBox: true, valign: "top", paraSpaceAfter: 4 });
  const gx = 8.8, gw = 3.64;
  s.addText("Gates and definitions", { x: gx, y: 1.88, w: gw, h: 0.28, fontSize: 13, bold: true, color: INK, margin: 0, isTextBox: true });
  rows(s, gx, 2.25, gw, 1.25, [
    ["Identity", "compiled vs Python K2A ≤ 10^{−12}: 194/194"],
    ["Moments", "archive vs export 10^{−9}: 193/195 (2 FAIL)"],
    ["Grouping", "production = S0 grouping, bit-identical: 192/192"],
    ["Time grid", "identical 10 ns bins in every comparison"],
    ["Criteria v1", "frozen 27 Nov 2026; 9 observables, both planes, Tp 100 and 200 ns"],
    ["Windows", "unshaped 0–7 µs; shaped 0–6 µs; CR-RC³ stand-in, not VMM3a"],
  ], { rh: 0.6, fs: 10.5, lfs: 11.5 });
  N(s, "The design isolates the model error. Both sides start from the same archived Garfield++ avalanche, so there is no sampling difference between them.\n\nThe sample: 980 non-sealed muons, of which 195 have been processed, and 192 form the accuracy set. Three are outside it. 0782 was selected on purpose as an out-of-domain event and K2A correctly refuses it. 0790 and 0877 failed an execution-consistency check on a diagnostic moment; I explain that in the appendix, and their frozen outcome stays FAIL. The sealed holdout has not been opened.\n\nFour signals are compared. S0_T1 is the decisive reference. S0_T0 uses true per-ion birth times. K2A is the prediction and uses its own production grouping, receiving nothing from S0. K2Aimp imposes the S0 grouping and is a diagnostic only. Production and S0 grouping turned out identical in all 192 events.\n\nThe gates on the right are consistency checks, not accuracy: identical results from compiled and Python K2A, identical moments from two independent paths, identical time grids. The criteria were frozen on 27 November before these comparisons. Unshaped comparisons cover 0 to 7 µs, shaped ones 0 to 6 µs, with a CR-RC³ stand-in shaper at 100 and 200 ns, which is not the VMM3a response.");
}
// =====================================================================================================  11 distributions
{
  const s = content("03 · ACCURACY", "Shaped errors are 0.18 % per plane, with no trend in the tested range", "03 Accuracy");
  img(s, "dists", X0, 1.84, 6.35, 4.45);
  img(s, "complexity", 7.45, 1.84, 5.0, 4.45, "right");
  note(s, X0, 6.32, 11.5, 0.45, ["Tested range: 7–26 clouds per event, all six stations, three channels. Ion-only shaped L2 (no electron step): mean 0.19 %, max 0.31 %."], 11);
  N(s, "Left, the error distributions over the 192 events for both planes. Shaped L2 at 100 ns: mean 0.178 % for U and 0.179 % for V, RMS across events about 0.023 %; at 200 ns slightly smaller. The unshaped current L2 is larger, about 0.52 %, because the unshaped current keeps the early high-frequency structure where the single-radius closure matters; shaping suppresses it. The window charge difference is centred on zero with an RMS of 0.22 %; that one matters for the charge criteria.\n\nWith the electron step removed, the ion-only shaped L2 is slightly larger, mean 0.19 % and maximum 0.31 %, so the electron step does not hide the ion error.\n\nRight, the dependence on event complexity: shaped L2 against the number of clouds, and the mean per station and channel. There is no visible trend; the rank correlation with the number of clouds is −0.06. But the tested range is narrow, 7 to 26 clouds, because the queue runs small events first. Larger events are a separate statement; the two high-multiplicity stress events in the appendix pass, but they are only two.");
}
// =====================================================================================================  12 observables
{
  const s = content("03 · ACCURACY", "Centroids and timing are far inside the limits; charge is not", "03 Accuracy");
  img(s, "observables", X0, 1.84, 6.75, 4.5);
  const x = 7.95, w = 4.49;
  stat(s, x, 1.89, w, "0.34 µm", "RMS of Δu (max 1.7 µm; limit 50 µm); u amplifies strip-centroid differences 10× relative to w", { lh: 0.55 });
  stat(s, x, 3.0, w, "0.03 µm", "RMS of Δw (max 0.14 µm; limit 10 µm)");
  stat(s, x, 3.98, w, "+0.10 ns", "systematic Δt at the 10 % threshold, RMS 0.01 ns (limit 1 ns)");
  stat(s, x, 4.96, w, "≤ 0.50 %", "window charge per plane, |ΔQ|/Q: limit 0.5 %; 4 events above 80 % of it", { color: RED });
  note(s, X0, 6.36, 11.5, 0.42, ["Δu, Δw are differences of the peak-weighted U/V strip centroid between K2A and S0, not tracking resolution."], 11);
  N(s, "These are the observables that feed reconstruction, as defined in the frozen criteria.\n\nPosition: the peak-weighted strip centroids of U and V are converted to w along the wires and u across them. Because the stereo angle is small, u = 1.7 mm × (V − U)/0.2, so any centroid difference is amplified about ten times more in u than in w. Even so the RMS difference is 0.34 µm in u and 0.03 µm in w, against limits of 50 and 10 µm. These are differences between two simulations, not a position resolution.\n\nTiming: the 10 % threshold crossing is systematically 0.10 ns later in K2A with an RMS of only 0.01 ns, and the lead peak is 0.14 % lower. Both are the single-radius effect seen in the waveforms.\n\nCharge is the weak point. The charge collected per plane in the window differs by up to 0.50 %, essentially at the limit, and four events use more than 80 % of it; the U/V charge ratio reaches 0.93 % against 1 %. All pass, but there is no margin. I do not yet know the cause. Candidates are the D and radial interpolation of the table at late times, which the kernel probe estimates at a few tenths of a percent, and the closure itself; this needs to be resolved before the sealed run.");
}
// =====================================================================================================  14 limits
{
  const s = content("04 · LIMITS", "Validity is shown for small events, 0–7 µs, 2100 V and B = 0", "04 Cost and limits");
  img(s, "selection_bias", X0, 1.84, 6.3, 3.9);
  note(s, X0, 5.85, 6.3, 0.85, ["The queue runs the smallest events of each station × channel first: the 192 events cover ≤ 26 clouds, where 38 % of the non-sealed population lies. Two stress events (141 and 304 clouds) pass but have separate provenance."], 10.5);
  rows(s, 7.45, 1.86, 5.0, 1.32, [
    ["Domain", "929 / 977 events evaluable: 0664, 0599 need anodes outside −1…1; 46 have one cloud above the radial grid (< 10^{−5} of the charge)"],
    ["Horizon", "tables end at 8 µs; 0–7 µs holds ≈ 15 % of the collecting-strip charge; the 43–98 µs tail is not modelled"],
    ["Early residual", "single r_{eff} dominates the unshaped error (0637: 10–50 ns)"],
    ["Charge", "window charge and U/V ratio near their limits; cause not identified"],
    ["Execution", "0790, 0877: t_{rms} diagnostic above its 10^{−5} ns guard (numerical cancellation, likely); FAIL kept"],
    ["References", "electron step approximate in S0 and K2A; CR-RC³ stand-in, not VMM3a"],
  ], { rh: 0.78, fs: 10.5, lfs: 11.5 });
  N(s, "Where the result does not yet apply.\n\nSample: the left plot compares the processed events with all non-sealed events. The queue processes the smallest events of each station and channel first, so the 192 events have at most 26 clouds and at most about 6,700 archived ion samples; 38 % of the population lies in that range and the median event has 32 clouds. The two stress events with 141 and 304 clouds pass, but they were run from a different commit and are reported separately.\n\nDomain: 929 of 977 events are inside the kernel table. Two events, 0664 and 0599, have most of their charge on anodes the table does not cover, because the table holds only anodes −1 to +1. Forty-six events have a single cloud just above the radial grid, carrying less than 10^-5 of the charge; K2A refuses them rather than clamping.\n\nHorizon: the tables stop at 8 µs, and within 7 µs only about 15 % of the collecting-strip charge has been induced; the remaining tail up to about 100 µs is not modelled. Whether that matters depends on the electronics window and the pile-up rate, which are not fixed yet.\n\nEarly residual: the closure anatomy of 0637 shows that the single radius dominates the unshaped error, concentrated in 10 to 50 ns.\n\nCharge criteria: nearly saturated, cause open. Execution: 0790 and 0877 failed a diagnostic check, explained in the appendix. References: the electron term is an instantaneous step on both sides, and the shaper is a stand-in.");
}
const hdrC = { bold: true, color: INK, fontSize: 11, fill: { color: TINT } };
const cell = (t, o = {}) => ({ text: t, options: { fontSize: 10.5, color: INK3, ...o } });
// A1
{
  const s = content("A1 · CRITERIA", "Frozen criteria v1 against the observed values (192 events)", "Appendix");
  const M = JSON.parse(fs.readFileSync("/home/claude/pres/an/margin.json"));
  const pct = (r) => `${(100 * r).toFixed(0)} %`;
  const R = [[cell("Observable (each plane, Tp 100 and 200 ns)", hdrC), cell("Limit", hdrC), cell("Median use", hdrC), cell("Max use", hdrC), cell("Worst event", hdrC)]];
  const rowsA = [["shaped_L2_ceiling", "shaped relative L2", "≤ 1 %"], ["window_charge", "window charge per plane |ΔQ|/Q", "≤ 0.5 %"],
                 ["UV_charge_ratio", "U/V charge ratio", "≤ 1 % rel."], ["threshold_time_f0.1", "time at 10 % threshold", "≤ 1 ns"],
                 ["second_over_lead", "second / lead peak", "≤ 0.005"], ["u_centroid", "u centroid", "≤ 50 µm"], ["w_centroid", "w centroid", "≤ 10 µm"]];
  rowsA.forEach(([k, lab, lim]) => {
    const m = M[k];
    R.push([cell(lab), cell(lim), cell(pct(m.median_ratio)), cell(pct(m.max_ratio), m.max_ratio > 0.8 ? { color: RED, bold: true } : {}), cell(`${m.worst_slot} (${m.worst_key.split("_")[0] === "UV" ? "U/V" : m.worst_key.split("_")[0]}, Tp ${m.worst_Tp})`)]);
  });
  R.push([cell("lead-strip identity"), cell("no change"), cell("—"), cell("0 of 768 decisions"), cell("—")]);
  R.push([cell("hit flips at 5 % and 10 %"), cell("0, marginal exempt"), cell("—"), cell("2 of 1,536 (both exempt)"), cell("0011 U 5 %, 0646 V 10 %, Tp 200")]);
  s.addTable(R, { x: X0, y: 1.9, w: 11.556, colW: [3.9, 1.7, 1.5, 2.3, 2.156], rowH: 0.38, border: { type: "solid", color: RULE, pt: 0.75 }, fill: { color: "FFFFFF" }, margin: 0.05 });
  note(s, X0, 5.95, 11.5, 0.75, ["Use = |K2A − S0_T1| / limit, worst plane and Tp per event. Thresholds fixed from the S0 lead peak (2, 5, 10, 20, 50 %); exemptions are defined from S0 only (S0 peak within ±2 % of the threshold). Coverage requirement for the sealed run: ≥ 17 of 20 events in domain."], 10.5);
  N(s, "The full table of the frozen criteria with what the 192 events actually use. The binding observables are the two charge criteria. All others use at most 30 % of their limits.");
}
// A2
{
  const s = content("A2 · DEFINITIONS", "Definitions and conventions", "Appendix");
  rows(s, X0, 1.88, 11.556, 2.3, [
    ["Cloud", "(Heed primary i, anode k = rint(u/3 mm), half to even); u in the nine-anode patch frame, r and θ = atan2(v, u − 3k) about anode k"],
    ["Strip frame", "U = w − 0.1 u, V = w + 0.1 u; pitch 1.7 mm; D_{U}, D_{V} without the ∓0.1 u term (added per anode by the kernel); absolute strips via per-primary recentring"],
    ["Signals", "ion term = Shockley–Ramo along frozen RKF paths; electron term = instantaneous step at the birth point; 10 ns bins; sign as in the S0 replay"],
    ["Shaping", "unit-peak CR-RC³, τ = Tp/3, Tp = 100 and 200 ns, applied to the 0–6 µs current; output in fC-equivalent peak (stand-in, not VMM3a)"],
    ["L2", "‖y_{K2A} − y_{S0}‖ / ‖y_{S0}‖ over all strips of a plane and the window; unshaped: current, 0–7 µs"],
    ["Centroid", "peak-weighted strip centroid of strips above the lowest threshold; w = 0.85 mm (U + V), u = 8.5 mm (V − U) in strip units"],
    ["T1 / T0", "T1: all ions of a primary at its ion-weighted mean birth time; T0: each ion at its own birth time"],
  ], { rh: 0.64, fs: 11, lfs: 12 });
  N(s, "Reference for the conventions. Note the centroid conversion: in strip units, w = 0.85 mm × (U + V) and u = 8.5 mm × (V − U), which is the factor of ten leverage of u relative to w.");
}
// A3
{
  const s = content("A3 · RESPONSE TABLE", "The response table and its numerical accuracy", "Appendix");
  img(s, "kernel", X0, 1.84, 11.556, 3.4, "center");
  rows(s, X0, 5.32, 11.556, 2.4, [
    ["Axes", "anodes −1, 0, 1 · r 15.0–135.2 µm, 14 log nodes (×1.184) · θ 72 · D −10.2…10.2 mm, step 0.17 · τ 0, 0.01 ns…8 µs, 200 nodes (7.1 %) · 57 MB float32"],
    ["Probe (read-only)", "τ interpolation ≤ 0.11 % of the early bin; radial leave-one-out 1.2 % at 10 ns → 0.02 % at 7 µs (mid-interval ≈ ¼); D leave-one-out ≈ 0.43–0.50 %"],
  ], { rh: 0.68, fs: 10.5, lfs: 11.5 });
  N(s, "Left: the monopole kernel, the induced charge per unit ion charge on the strip above a cloud, for several birth radii. Right: the same kernel relative to its value at 15 µm, as a function of radius, at four times. Early on, the response depends strongly on the birth radius; that is why one effective radius per cloud is the dominant approximation early and becomes harmless later.\n\nThe probe numbers are leave-one-out estimates of the interpolation error, an upper bound; the error at mid-interval is roughly a quarter of it. The D interpolation at late τ, around half a percent leave-one-out, is one candidate for the charge residual.");
}
// A4
{
  const s = content("A4 · REPRODUCIBILITY", "Execution gates, and the two moment-check failures", "Appendix");
  const R = [[cell("Gate", hdrC), cell("Tolerance", hdrC), cell("Result", hdrC)],
             [cell("compiled K2A vs Python K2A"), cell("10⁻¹²"), cell("194 / 194 PASS (max ion L2 7 × 10⁻¹⁵)")],
             [cell("responses from archive vs export moments"), cell("10⁻⁷"), cell("194 / 194 PASS (max 8 × 10⁻¹⁵)")],
             [cell("cloud moments, archive vs export path"), cell("10⁻⁹; t_rms: |t_rms| < 10⁻⁵ ns"), cell("193 / 195 PASS; 0790, 0877 FAIL", { color: RED })],
             [cell("time grid of all compared signals"), cell("identical"), cell("195 / 195 PASS")],
             [cell("production vs S0 grouping"), cell("bit-identical"), cell("192 / 192 (AGREE in all 195 slots)")],
             [cell("repeated K2A evaluations"), cell("bit-identical"), cell("195 / 195")]];
  s.addTable(R, { x: X0, y: 1.9, w: 11.556, colW: [4.3, 3.0, 4.256], rowH: 0.37, border: { type: "solid", color: RULE, pt: 0.75 }, fill: { color: "FFFFFF" }, margin: 0.05 });
  rows(s, X0, 4.45, 11.556, 2.1, [
    ["0790, 0877", "only t_{rms} differs (6.2 × 10⁻⁷ and 8.2 × 10⁻⁷ ns); t ≤ 1.2 × 10⁻¹³ ns, r_{eff} ≤ 8 × 10⁻¹⁶ rel., harmonics ≤ 3 × 10⁻¹⁶; N, D, A identical"],
    ["Why", "t_{rms} is the square root of (mean t² − (mean t)²); after the T1 collapse it is zero by construction, so only rounding cancellation remains, growing with t. In 0790 it reached ≈ 1.07 × 10⁻⁵ ns, above the 10⁻⁵ ns guard. t_{rms} is not an evaluator input; responses from both moment paths agree to 3 × 10⁻¹⁵. Likely numerical cancellation; root cause under investigation"],
    ["Outcome", "frozen execution outcome FAIL kept; not evaluated for accuracy"],
  ], { rh: 0.72, fs: 10.5, lfs: 11.5 });
  N(s, "The gates check that the pipeline computes what it claims, independently of accuracy.\n\nThe two failures are in the comparison of cloud moments built along two independent paths, from the archive and from the frozen export. Only the diagnostic t_rms column differs. After the time collapse every ion of a cloud has the same time, so t_rms is zero by construction; the computed value is the square root of a rounding-level cancellation and grows with the cloud time. In 0790 it reached about 1.07 × 10^-5 ns, just above the absolute guard of 10^-5 ns. 0877 shows the same signature; its guard value is not in the collected files. All physical moments agree to floating-point precision, and the K2A responses computed from the two moment paths agree to 3 × 10^-15. This is most likely numerical cancellation, not an accuracy failure, but the frozen outcome stays FAIL and the events stay out of the accuracy set.");
}
// A5
{
  const s = content("A5 · THRESHOLDS", "Decisions flip only where the reference sits on a threshold", "Appendix");
  const R = [[cell("Event", hdrC), cell("Decision", hdrC), cell("S0 peak / threshold", hdrC), cell("Effect", hdrC), cell("In criteria?", hdrC)],
             [cell("0011"), cell("U hit at 5 %, Tp 200 ns"), cell("0.99984 (0.016 % below)"), cell("one strip crosses in K2A only"), cell("yes: marginal exemption")],
             [cell("0646"), cell("V hit at 10 %, Tp 200 ns"), cell("0.99986 (0.014 % below)"), cell("one strip crosses in K2A only"), cell("yes: marginal exemption")],
             [cell("0960"), cell("centroid membership at 2 %, ion only, Tp 100 ns"), cell("0.99995 (0.005 % below)"), cell("Δu = 230 µm, Δw = 23 µm; with fixed membership Δu = 0.26 µm"), cell("no (ion-only diagnostic)")]];
  s.addTable(R, { x: X0, y: 1.9, w: 11.556, colW: [1.0, 3.0, 2.3, 3.2, 2.056], rowH: 0.5, border: { type: "solid", color: RULE, pt: 0.75 }, fill: { color: "FFFFFF" }, margin: 0.05 });
  note(s, X0, 4.4, 11.5, 1.2, ["Threshold decisions and centroid membership are discontinuous: a 0.01 % amplitude difference can flip them when the reference is within 0.01 % of the threshold. The frozen criteria exempt such cases using S0 alone. With the electron step included (the decisive signal), the largest Δu is 1.7 µm."], 11.5);
  N(s, "Why perfect agreement of every discrete decision is not expected. In all three cases the S0 peak lies within two hundredths of a percent of the threshold, and K2A, being about 0.01 % different, lands on the other side. 0960 is not part of the criteria: in the ion-only signal one strip sits at the 2 % centroid-inclusion threshold, and including it moves the u centroid by 230 µm; with the membership fixed from S0 the difference is 0.26 µm.");
}
// A6
{
  const s = content("A6 · TIME COLLAPSE", "One time per primary costs ≤ 0.1 % after shaping", "Appendix");
  const st = JSON.parse(fs.readFileSync("/home/claude/pres/an/stats.json"));
  const R = [[cell("Comparison (192 events)", hdrC), cell("Shaped L2, Tp 100 ns: median", hdrC), cell("max", hdrC), cell("Unshaped ion L2: median", hdrC), cell("max", hdrC)],
             [cell("S0_T0 vs S0_T1 (time collapse alone)"), cell("0.0011 %"), cell("0.10 %"), cell("0.010 %"), cell("0.18 %")],
             [cell("S0_T1 vs K2A (model, decisive)"), cell("0.17 %"), cell("0.27 %"), cell("0.56 %"), cell("0.92 %")],
             [cell("S0_T0 vs K2A (total)"), cell("0.17 %"), cell("0.27 %"), cell("0.56 %"), cell("0.91 %")]];
  s.addTable(R, { x: X0, y: 1.9, w: 11.556, colW: [3.8, 2.4, 1.4, 2.556, 1.4], rowH: 0.42, border: { type: "solid", color: RULE, pt: 0.75 }, fill: { color: "FFFFFF" }, margin: 0.05 });
  note(s, X0, 4.0, 11.5, 1.0, ["The collapse matters mainly in the unshaped ion+electron signal (median 0.011 %, max 3.5 %), because the electron step moves in time; after shaping it stays ≤ 0.1 %. All 192 events also pass criteria v1 against S0_T0 (information only)."], 11.5);
  N(s, "S0_T0 keeps every ion's own birth time; S0_T1 puts all ions of a primary at one time, like K2A. The difference between them is the cost of the time collapse alone, and it is small after shaping. The unshaped ion-plus-electron difference can reach a few percent in individual events because the instantaneous electron step moves by up to tens of nanoseconds, but the shaper integrates that away.");
}
// A7
{
  const s = content("A7 · WORST EVENT", "The event with the largest shaped L2 (0581) looks like the typical one", "Appendix");
  img(s, "waves_worst", X0, 1.84, 11.556, 4.4, "center");
  note(s, X0, 6.3, 11.5, 0.45, ["U plane L2 0.27 %, V plane 0.20 % (Tp = 100 ns). Same early-edge residual, ≈ 0.3 % of the lead peak."], 11);
  N(s, "For completeness, the event with the largest shaped L2 of the 192. The residual has the same shape as in the typical event, slightly larger. Nothing qualitatively new happens in the tail of the distribution.");
}
// A8
{
  const s = content("A8 · RESIDUAL AND HORIZON", "The early residual is the radial closure; most charge arrives after 7 µs", "Appendix");
  img(s, "early_0637", X0, 1.84, 5.7, 4.45);
  img(s, "horizon", 6.75, 1.84, 5.7, 4.45, "right");
  note(s, X0, 6.32, 11.5, 0.45, ["0637 per-group synthesis: SR vs G (one r_{eff}) ≈ K2A vs S0, both concentrated in 10–50 ns; angular truncation (S3A) small. Right: ion term of three examples, 100 µs replay."], 10.5);
  N(s, "Left: the closure anatomy of 0637, the event with the most visible early residual. Replacing each ion's radius by the cloud's effective radius, SR versus G, gives almost the same error as the full K2A against S0, and both sit in the first 10 to 50 ns. The angular truncation is small. This identifies the source of the early residual; a second radial moment is the natural fix if an observable ever needs it.\n\nRight: cumulative collecting-strip charge for three events over the full 100 µs S0 replay. About 15 % of the final charge is induced by 7 µs. The rest arrives as the ions approach the cathode, which K2A does not model at present.");
}
// A9
{
  const s = content("A9 · RUNTIME", "Runtime definitions (medians over 192 events, seconds)", "Appendix");
  const R = [[cell("Quantity", hdrC), cell("Definition", hdrC), cell("Median", hdrC)],
             [cell("S0_T0 total"), cell("RKF export (archive read + ion drift) + grouping + weighting replay (chosen single pass) + electron step + local sum"), cell("84.6 s")],
             [cell("K2A archive I/O"), cell("read and validate the K256 archive (minimum of 5 repeats)"), cell("7.4 ms")],
             [cell("K2A moments"), cell("moment construction from the births (minimum of 5 repeats)"), cell("2.5 ms")],
             [cell("K2A evaluation"), cell("warm kernel evaluation, all strips, 700 bins (median of 20 repeats, bit-identical)"), cell("1.6 ms")],
             [cell("Kernel open"), cell("verify sha256 and map 57 MB (once per job)"), cell("0.16 s")],
             [cell("Speedup: evaluation"), cell("(S0_T0 total − approximate archive read) / K2A evaluation"), cell("54,180×")],
             [cell("Speedup: + moments"), cell("same numerator / (moments + evaluation)"), cell("21,228×")],
             [cell("Speedup: end to end"), cell("S0_T0 total / (archive I/O + moments + evaluation)"), cell("7,568×")]];
  s.addTable(R, { x: X0, y: 1.9, w: 11.556, colW: [2.6, 7.4, 1.556], rowH: 0.4, border: { type: "solid", color: RULE, pt: 0.75 }, fill: { color: "FFFFFF" }, margin: 0.05 });
  note(s, X0, 5.75, 11.5, 0.95, ["One CPU process per step on the same laptop. Avalanche generation excluded on both sides. A few sub-second S0 steps are recorded with small negative wall times (timing artefact, |t| < 3 s, negligible against 85 s); event 0262 had a 12,000 s RKF wall time (not investigated)."], 10.5);
  N(s, "Exact definitions behind the speedup numbers. The archive read subtracted from S0 is approximate because the frozen exporter reads and drifts in one process. Some small S0 steps have negative recorded times, an artefact of the timing, at most a few seconds against a total of 85 s; it does not change any median.");
}
// A10
{
  const s = content("A10 · PROVENANCE", "Snapshot, new aggregate and stress events are kept separate", "Appendix");
  const R = [[cell("Report", hdrC), cell("Slots", hdrC), cell("Accuracy set", hdrC), cell("Shaped L2 median U / V", hdrC), cell("max |Δu|", hdrC), cell("Speedup eval / e2e", hdrC), cell("Commit", hdrC)],
             [cell("snapshot 2026-11-28 11:51"), cell("155"), cell("153 PASS"), cell("0.173 / 0.176 %"), cell("1.20 µm"), cell("54,500× / 7,610×"), cell("437946a")],
             [cell("new aggregate 14:34 (this talk)", { bold: true }), cell("195"), cell("192 PASS"), cell("0.173 / 0.175 %"), cell("1.49 µm"), cell("54,180× / 7,568×"), cell("437946a")],
             [cell("stress: 0758 (304 clouds)"), cell("1"), cell("PASS"), cell("0.196 / 0.179 %"), cell("—"), cell("S0 1,520 s; K2A eval 27 ms"), cell("b87f710")],
             [cell("stress: 0004 (141 clouds)"), cell("1"), cell("PASS"), cell("0.184 / 0.223 %"), cell("—"), cell("S0 680 s; K2A eval 12 ms"), cell("b87f710")]];
  s.addTable(R, { x: X0, y: 1.9, w: 11.556, colW: [2.75, 0.75, 1.3, 1.9, 1.1, 2.4, 1.356], rowH: 0.42, border: { type: "solid", color: RULE, pt: 0.75 }, fill: { color: "FFFFFF" }, margin: 0.05 });
  note(s, X0, 4.3, 11.5, 1.6, ["All runs: one kernel (kernel_r14t72edge5), one strip-evaluator build, one birth-moment binary, criteria v1 sha256 efb4c833…; the stress events ran from a later harness commit and are not merged with the 192.",
    "Collected upload: 12,742 files, manifest sha256 per file; 0 essential items missing; sealed slots never copied.",
    "Max |Δu| of the new aggregate is at Tp = 100 ns (1.7 µm at Tp = 200 ns)."], 11);
  N(s, "Provenance. The original 153-event snapshot is unchanged and agrees with the new aggregate within the statistical spread. The two high-multiplicity stress events pass but were produced with a later commit of the harness, so they stay separate until equivalence is shown.");
}
// A11
{
  const s = content("A11 · DOMAIN", "48 of 977 events are refused, never clamped", "Appendix");
  const R = [[cell("Cause", hdrC), cell("Events", hdrC), cell("Charge affected", hdrC), cell("Fix", hdrC)],
             [cell("anode outside the table (−1…1)"), cell("2 (0664, 0599)"), cell("78 % and 50 % of the event"), cell("add anodes ±2, ±3 (strict superset)")],
             [cell("r_eff above the radial grid"), cell("46 (one cloud each)"), cell("< 10⁻⁵ of the event; r_eff / r_max 1.07–2.6"), cell("extend radii outward")],
             [cell("evaluated"), cell("929 (95.1 %)"), cell("—"), cell("—")]];
  s.addTable(R, { x: X0, y: 1.9, w: 11.556, colW: [3.3, 2.2, 3.4, 2.656], rowH: 0.42, border: { type: "solid", color: RULE, pt: 0.75 }, fill: { color: "FFFFFF" }, margin: 0.05 });
  note(s, X0, 3.85, 11.5, 1.3, ["0664 (station 5, φ): charge on patch anodes 0–3, 60 % on anode 2 · 0599 (station 1, ω, θ_{u} = 11.9°): 31 % on anode 2, 19 % on anode 3. The kernel's anode range follows from the five development muons used to build it. All 3,024 grid paths end on the physical cathode (|v| = 3 mm within 3 × 10⁻⁶ mm), drift 43–98 µs."], 11);
  N(s, "The domain scan ran K2A alone over all 977 non-sealed events. Every refused event is listed with its cause; nothing is clamped or silently dropped. The two wire-domain events need the strict-superset kernel; the radial outliers are single clouds with negligible charge.");
}

pres.writeFile({ fileName: "/home/claude/pres/deck/K2A_signal_response_simple.pptx" }).then(async (f) => {
  await applyTheme(f, THEME);
  console.log("written", f);
});
