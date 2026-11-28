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
function content(kick, title, sec) {
  const s = pres.addSlide({ masterName: "CONTENT", sectionTitle: sec });
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

// =====================================================================================================  1 title
pres.addSection({ title: "Introduction" });
{
  const s = pres.addSlide({ masterName: "DARK", sectionTitle: "Introduction" });
  s.addShape(pres.shapes.RECTANGLE, { x: X0, y: 2.181, w: 1.111, h: 0.042, fill: { color: RED }, line: { color: RED, width: 0 } });
  s.addText([{ text: "Fast cathode response", options: { breakLine: true } }, { text: "for the NA60+ MWPC" }], { x: X0, y: 2.45, w: 11.9, h: 1.45, fontSize: 46, bold: true, color: DTXT, margin: 0, isTextBox: true, valign: "top" });
  s.addText("K2A: compact ion-cloud moments against ion-by-ion drift", { x: X0, y: 4.0, w: 11.9, h: 0.45, fontSize: 24, color: DTXT, margin: 0, isTextBox: true });
  s.addText("Garfield++ K256 archives · idealized 9-anode cell · Ar/CO₂ 70:30 · 2100 V · B = 0", { x: X0, y: 4.62, w: 11.9, h: 0.365, fontSize: 18, color: DSUB, margin: 0, isTextBox: true });
  s.addText("192 independent development muons against the detailed S0 reference · criteria frozen in advance", { x: X0, y: 5.02, w: 11.9, h: 0.365, fontSize: 18, color: DSUB, margin: 0, isTextBox: true });
  s.addText("Development validation on non-sealed events. The sealed 20-event holdout has not been opened · November 2026", { x: X0, y: 6.367, w: 11.9, h: 0.272, fontSize: 13, color: MUTED, margin: 0, isTextBox: true });
  N(s, "This talk follows the avalanche study. That deck measured what single avalanches and ion clouds look like at 2100 V; its first open item was that the archives contain ion births only, not cathode signals. Here we ask whether the induced U/V strip signal can be computed from a few numbers per ion cloud, instead of drifting every sampled ion through the field and the weighting potential.\n\nEverything shown is a development validation on non-sealed muons at one configuration: 2100 V, B = 0, the idealized nine-anode cell and Ar/CO2 70:30. The sealed 20-event holdout is still closed, and the acceptance criteria were frozen before these comparisons.");
}
// =====================================================================================================  2 questions
{
  const s = pres.addSlide({ masterName: "DARK", sectionTitle: "Introduction" });
  kicker(s, "QUESTIONS FOR THIS DECK", true);
  s.addText("Can ion-by-ion cathode-response calculations be replaced by a compact moment-based kernel, without losing the observables used for reconstruction?",
    { x: X0, y: 1.172, w: 11.6, h: 0.95, fontSize: 24, bold: true, color: DTXT, margin: 0, isTextBox: true, valign: "top" });
  const q = [["01 · SIGNAL", "What has to be computed, and why is the detailed calculation slow?"],
             ["02 · MODEL", "What does K2A keep from each ion cloud, and what does it approximate?"],
             ["03 · ACCURACY", "How close is K2A to the detailed reference, in signals and in detector observables?"],
             ["04 · COST AND LIMITS", "How much faster is it, and where is it not yet valid?"]];
  q.forEach(([k, t], i) => {
    const x = i % 2 ? 6.917 : X0, y = i < 2 ? 2.55 : 3.87;
    s.addShape(pres.shapes.LINE, { x, y, w: 5.528, h: 0, line: { color: RED, width: 1 } });
    s.addText(k, { x, y: y + 0.142, w: 5.583, h: 0.228, fontSize: 12, bold: true, color: LRED, margin: 0, isTextBox: true });
    s.addText(t, { x, y: y + 0.412, w: 5.4, h: 0.653, fontSize: 18, color: DTXT, margin: 0, isTextBox: true, valign: "top" });
  });
  s.addText("The avalanche study established the properties of ion clouds. K2A uses those properties to compute their induced cathode response efficiently; this deck tests whether that is good enough.",
    { x: X0, y: 5.55, w: 11.5, h: 0.6, fontSize: 14, color: DSUB, margin: 0, isTextBox: true, valign: "top" });
  N(s, "Four questions, in the order of the talk. The main question is quantitative: replacing the ion-by-ion calculation is acceptable only if the detector observables that reconstruction uses are preserved, and only within the domain we actually tested.\n\nThe link to the previous deck is direct: the avalanche study told us how big, how extended and how polarised the ion clouds are. K2A summarises each cloud by a few moments that follow from those properties and looks the response up in precomputed tables.");
}
// =====================================================================================================  3 main result
pres.addSection({ title: "Answer" });
{
  const s = content("ANSWER IN ONE SLIDE", "All 192 events pass the frozen criteria; the charge criteria are tightest", "Answer");
  img(s, "criteria_margin", X0, 1.86, 7.25, 4.35);
  const x = 8.45, w = 3.99;
  stat(s, x, 1.89, w, "192 / 192", "independent events pass frozen criteria v1 (K2A vs S0_T1; non-sealed development set)");
  stat(s, x, 2.93, w, "0.18 %", "mean shaped L2 per plane, Tp = 100 ns (largest event 0.27 %; 0.30 % at Tp = 200 ns)");
  stat(s, x, 3.97, w, "99 %", "of the window-charge limit used by the worst event (0740, V plane); U/V charge ratio up to 93 %", { color: RED });
  stat(s, x, 5.01, w, "5.4 × 10⁴", "median speedup of the response evaluation; 7.6 × 10³ including archive I/O and moments", { lh: 0.5 });
  note(s, X0, 6.28, 11.5, 0.45, ["Discrete decisions: 0 lead-strip changes in 768; 2 threshold flips in 1,536, both with the S0 peak within 0.02 % of the threshold (marginal exemption defined in v1). Avalanche generation is not included in any speedup."], 11);
  N(s, "This is the whole answer on one slide; the rest of the talk explains and qualifies it.\n\nThe plot shows, for each frozen v1 criterion, how much of its limit each of the 192 events uses: the worst plane and worst shaping time per event. Shaped waveform L2, timing, second-over-lead and the two centroid criteria use at most 30 % of their limits; the centroids use only a few percent.\n\nThe two charge criteria are different. The window-charge difference per plane reaches 99 % of its 0.5 % limit in event 0740, and the U/V charge ratio reaches 93 % of its 1 % limit. All events pass, but these two criteria have essentially no margin, and the cause is not yet identified. I would not claim more than 'passes' for them.\n\nDiscrete decisions: no lead-strip change in 768 plane-by-shaping decisions. Two hit flips out of 1,536 threshold decisions, both at Tp = 200 ns, and in both the S0 peak lies within 0.02 % of the threshold; the frozen criteria list these as marginal exemptions, so this is not perfect agreement of every decision.\n\nThe speedup is for the cathode-response calculation only, from the same archived avalanche on both sides.");
}
// =====================================================================================================  4 signal
pres.addSection({ title: "01 Signal" });
{
  const s = content("01 · SIGNAL", "Slow ions induce the strip signal; only ≈ 15 % arrives within 7 µs", "01 Signal");
  img(s, "cell_schematic", X0, 1.86, 5.25, 4.25);
  img(s, "signal_anatomy", 6.25, 1.86, 6.2, 4.25, "right");
  note(s, X0, 6.2, 11.5, 0.55, ["Shockley–Ramo: Q_{s}(t) = q [φ_{s}(x(t)) − φ_{s}(x_{0})] for each charge q on its path x(t); φ_{s} is the weighting potential of strip s (sign convention of the S0 replay: positive on the strip the ions approach). Electrons: instantaneous step at birth, an approximation in both S0 and K2A."], 11);
  N(s, "Left: the reference cell. Nine anodes of 30 µm diameter at 3 mm pitch between two cathodes at ±3 mm; the cathodes carry the U and V strips, 1.7 mm pitch, stereo ±0.1. The arrows only illustrate the ion drift; the real paths are computed with RKF in the analytic field.\n\nRight: what one strip actually sees, from the detailed S0 calculation for event 0219. In the first microsecond the ions leave the high-field region around the anode and produce a fast current on the collecting strip and its neighbour. After that the current is small but lasts until the ions reach the cathode after roughly 40 to 100 µs. On the collecting strip only about 16 % of the final induced charge has arrived by 7 µs; the neighbouring strip sees a bipolar signal that integrates to almost zero once the ions land.\n\nThe electron contribution is treated as an instantaneous step at the ion birth position in both S0 and K2A, so neither reference simulates the electron current in time.\n\nAll comparisons in this talk cover 0 to 7 µs unshaped and 0 to 6 µs after shaping. The late ion tail is a separate question; I come back to it in the limitations.");
}
// =====================================================================================================  5 S0
{
  const s = content("01 · SIGNAL", "The S0 reference follows every sampled ion; it costs ≈ 85 s per event", "01 Signal");
  const steps = [["Garfield++ archive", "K256 weighted ion births per Heed primary (frozen)"],
                 ["RKF ion drift", "every sampled ion, analytic nine-wire field (frozen exporter)"],
                 ["Weighting-potential replay", "φ_{s} along each path from the spectral U/V LUT"],
                 ["Electron term", "instantaneous step at the birth point (approximation)"],
                 ["Sum on absolute strips", "per recentred group, 10 ns bins, 100 µs"]];
  steps.forEach(([h, b], i) => {
    const y = 1.88 + i * 0.83;
    box(s, X0, y, 4.75, 0.7, `${i + 1}  ${h}`, [b], { bfs: 10.5 });  // one-line bodies
    if (i < steps.length - 1) arrow(s, X0 + 2.375, y + 0.7, X0 + 2.375, y + 0.83);
  });
  img(s, "s0_cost", 5.95, 1.86, 6.5, 4.1, "right");
  note(s, 5.95, 6.05, 6.5, 0.7, ["Two references from the same births: **S0_T1** puts every ion of a primary at the primary's weighted-mean birth time (decisive reference); **S0_T0** keeps every ion's own birth time."], 11);
  N(s, "S0 is the established detailed workflow. It starts from the frozen Garfield++ archive, which holds 256 weighted ion births per Heed primary. Each sampled ion is drifted with RKF through the analytic nine-wire field; the weighting potential of every strip is then evaluated along each path from the spectral lookup table, and the strip currents are summed on absolute strip numbers.\n\nThe electron term is not simulated in time: it is an instantaneous step at the birth point, the same approximation for all references.\n\nThere are two versions. S0_T1 collapses all ions of one primary to the primary's weighted-mean birth time; this is the decisive reference because K2A uses the same time. S0_T0 keeps every ion's own birth time and tells us what the time collapse itself costs.\n\nThe box plot is measured wall time per event for the 192 events, one CPU process on the laptop: median 85 s, dominated by the weighting-potential replay at about 49 s and the RKF drift at about 30 s. That is the cost K2A has to remove.");
}
// =====================================================================================================  6 moments
pres.addSection({ title: "02 Model" });
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
pres.addSection({ title: "03 Accuracy" });
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
// =====================================================================================================  10 waveforms
{
  const s = content("03 · ACCURACY", "Waveforms agree; K2A peaks are 0.14 % low and 0.1 ns late", "03 Accuracy");
  img(s, "waves_typical", X0, 1.84, 11.556, 4.4, "center");
  note(s, X0, 6.3, 11.5, 0.45, ["Median-L2 event 0219. Residual: a dip of ≈ 0.25 % of the lead peak at ≈ 100 ns on the rising edge, then agreement to ≈ 0.05 %. Worst event (0581) in the appendix."], 11);
  N(s, "This is a typical event, the one with the median shaped L2, 0219. The three largest strips in each plane: S0_T1 in black, K2A dashed red, after the 100 ns stand-in shaper. On this scale the curves are indistinguishable; the plane L2 is 0.17 %.\n\nThe residual panels show the difference relative to the lead-strip peak. The structure is systematic, not noise: K2A is slightly low on the rising edge, about 0.25 % of the peak around 100 ns, and then agrees to a few hundredths of a percent over the long tail.\n\nThe same pattern appears in all 192 events: the lead peak is on average 0.14 % low and the 10 % threshold crossing about 0.1 ns late. That is consistent with the closure study, where the single effective radius produces an early-time residual and a best-fit time shift of about 0.1 ns. It is a known, small, systematic effect of the single-radius approximation.\n\nThe event with the largest L2, 0581, looks the same and is in the appendix.");
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
// =====================================================================================================  13 runtime
pres.addSection({ title: "04 Cost and limits" });
{
  const s = content("04 · COST", "The response evaluation is ≈ 5 × 10⁴ times faster than ion-by-ion drift", "04 Cost and limits");
  img(s, "runtime", X0, 1.84, 11.556, 4.15, "center");
  rows(s, X0, 6.02, 11.556, 2.7, [["Ratios (medians, 192 events)", "evaluation only 54,180× · + moment construction 21,228× · + archive I/O 7,568× · single cold event incl. kernel open 464×"]],
    { rh: 0.62, fs: 11, lfs: 11.5 });
  N(s, "Left: measured wall time per event, one CPU process. S0_T0 takes a median of 85 s. K2A needs 7.4 ms to read and validate the archive, 2.5 ms to build the moments, and 1.6 ms to evaluate the response once the kernel is open. Opening and verifying the 57 MB kernel costs 0.16 s once per job.\n\nRight: the speedup per event for three definitions, against the number of clouds. Evaluation only, which is what a production chain with moments already in memory would pay: median 54,000 times. Including the moment construction: 21,000 times. End to end from the archive, including its I/O: 7,600 times. A single cold event including the kernel opening: 460 times.\n\nThree caveats. The S0 time counts the chosen single pass of its replay; ladder retries are excluded. The archive read subtracted from S0 is an approximation, because the frozen exporter reads and drifts in one process. And avalanche generation is excluded on both sides: none of these numbers says anything about the cost of producing the avalanche, which is Fast4's job. One event, 0262, had an RKF wall time of 12,000 s, probably machine contention; it is why a few points sit very high, and it does not affect the medians.");
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
// =====================================================================================================  15 conclusions
pres.addSection({ title: "Conclusions" });
{
  const s = content("SUMMARY", "Answer: yes, for the observables and the domain tested", "Conclusions");
  rows(s, X0, 1.88, 11.556, 2.65, [
    ["Demonstrated", "192/192 independent events pass frozen criteria v1 against S0_T1: shaped L2 0.18 % (max 0.30 %), centroid differences ≤ 1.7 µm (u) and 0.14 µm (w), Δt = +0.10 ns, no lead-strip change"],
    ["Approximation budget", "one r_{eff} per cloud dominates (≈ 0.5 % unshaped); n ≤ 2 + hemisphere ≈ 0.09 %; one D per cloud negligible; time collapse ≤ 0.1 % shaped"],
    ["Cost", "response evaluation 5.4 × 10⁴ faster than S0 (7.6 × 10³ with archive I/O); avalanche generation not included"],
    ["Provisional", "charge observables at the edge of their limits; development sample of small events; stand-in shaper; electron step approximate"],
    ["Not shown", "sealed holdout · events above 26 clouds (two exceptions) · ion tail beyond 7 µs · other voltages, fields, gases, geometries · Fast4 input"],
  ], { rh: 0.86, fs: 12, lfs: 13 });
  N(s, "To answer the question of the talk: yes, for the observables defined in the frozen criteria and for the domain we tested, the ion-by-ion calculation can be replaced by ten numbers per cloud and a precomputed kernel.\n\nDemonstrated: all 192 independent development events pass against the decisive reference; waveform errors are about 0.2 % after shaping; position differences are micrometres or less; timing shifts by a systematic tenth of a nanosecond; no lead strip changes.\n\nWe also know where the remaining error comes from: the single effective radius is the dominant approximation, the angular truncation with the hemisphere term costs less than a tenth of a percent, and the time collapse is negligible after shaping.\n\nThe cost of the response evaluation drops by more than four orders of magnitude, not counting avalanche generation.\n\nProvisional: the charge observables pass with little margin, the sample is biased to small events, and the electronics and electron term are approximations. Not shown at all: the sealed holdout, large events beyond two examples, the late ion tail, other operating conditions, and Fast4 as the input.");
}
// =====================================================================================================  16 next steps
{
  const s = content("NEXT STEPS", "Next: explain the charge residual, then widen the domain", "Conclusions");
  const steps = [["1", "Charge residual", "Locate the window-charge and U/V-ratio differences (table interpolation at late τ vs closure) without touching criteria v1."],
                 ["2", "Strict-superset kernel", "Add anodes ±2, ±3 and outward radii; old table entries must stay bit-identical; new kernel id."],
                 ["3", "Longer horizon", "All 3,024 grid paths end on the cathode (43–98 µs), so τ can be extended by re-tabulating existing paths."],
                 ["4", "Larger events", "Finish the 800-muon queue and re-check accuracy and cost against the number of clouds."],
                 ["5", "Sealed holdout", "Run the 20 sealed events once, with criteria v1 unchanged."],
                 ["6", "Production chain", "Fast4 → K2A interface at 2100 V; real VMM3a response and readout window."]];
  steps.forEach(([n, h, b], i) => {
    const col = i % 2, r = Math.floor(i / 2);
    const x = X0 + col * 5.9, y = 1.9 + r * 1.32;
    s.addShape(pres.shapes.OVAL, { x, y: y + 0.02, w: 0.42, h: 0.42, fill: { color: i < 1 ? RED : INK }, line: { color: i < 1 ? RED : INK, width: 0 } });
    s.addText(n, { x, y: y + 0.02, w: 0.42, h: 0.42, fontSize: 14, bold: true, color: "FFFFFF", align: "center", valign: "middle", margin: 0, isTextBox: true });
    s.addText(h, { x: x + 0.6, y, w: 5.0, h: 0.32, fontSize: 15, bold: true, color: INK, margin: 0, isTextBox: true });
    s.addText(b, { x: x + 0.6, y: y + 0.36, w: 5.0, h: 0.8, fontSize: 12, color: INK3, margin: 0, isTextBox: true, valign: "top" });
  });
  note(s, X0, 5.95, 11.5, 0.75, ["Open questions: is a second radial moment needed once the electronics window is fixed? Which readout window and pile-up rate decide whether the ion tail matters? Should the t_{rms} guard of the moment check be documented as a numerical tolerance (harness only, not the model)?"], 11);
  N(s, "Six steps, in order.\n\nFirst, before anything else, understand the charge residual. The criteria stay frozen; the point is to know whether it comes from the table interpolation at late times or from the closure, because that decides whether the fix is a finer table or a model change.\n\nSecond, a strict-superset kernel that adds anodes ±2 and ±3 and larger radii, with the existing entries bit-identical, so that the 48 refused events become evaluable without changing any validated result.\n\nThird, the time horizon: every grid path already ends on the cathode, between 43 and 98 µs, so a longer table needs only re-tabulation of existing paths, not new transport.\n\nFourth, finish the 800-muon queue, which brings in the larger events. Fifth, the sealed holdout, once, with the frozen criteria. Sixth, connect Fast4 at 2100 V and replace the stand-in shaper by the real VMM3a response.\n\nOpen questions are listed at the bottom; the most physical one is whether a second radial moment will be needed once we know the real electronics window.");
}

// ===================================================================================================== APPENDIX
pres.addSection({ title: "Appendix" });
{
  const s = pres.addSlide({ masterName: "DARK", sectionTitle: "Appendix" });
  kicker(s, "APPENDIX", true);
  s.addText("Technical definitions, validation criteria and reproducibility", { x: X0, y: 1.172, w: 11.6, h: 0.6, fontSize: 26, bold: true, color: DTXT, margin: 0, isTextBox: true });
  const items = ["A1  Frozen criteria v1 against the observed values", "A2  Definitions and conventions", "A3  Response table and its numerical accuracy",
                 "A4  Reproducibility gates; events 0790 and 0877", "A5  Decisions near a threshold", "A6  Time collapse: S0_T0 vs S0_T1",
                 "A7  Worst-L2 event", "A8  Early-time residual and ion-tail horizon", "A9  Runtime definitions",
                 "A10 Provenance: snapshot, new aggregate and stress events", "A11 Kernel domain"];
  s.addText(items.map((t, i) => ({ text: t, options: { breakLine: i < items.length - 1 } })), { x: X0, y: 2.2, w: 11, h: 4.2, fontSize: 16, color: DSUB, margin: 0, isTextBox: true, paraSpaceAfter: 6, valign: "top" });
  N(s, "Backup material for questions.");
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

pres.writeFile({ fileName: "/home/claude/pres/deck/K2A_signal_response.pptx" }).then(async (f) => {
  await applyTheme(f, THEME);
  console.log("written", f);
});
