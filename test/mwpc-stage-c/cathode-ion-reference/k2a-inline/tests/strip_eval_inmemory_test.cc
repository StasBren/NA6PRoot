// SYNTHETIC in-memory chain test: births -> K2AInline::MomentAccumulator -> K2AInline::Event -> K2AStrip::Evaluator,
// with no file in between, against the same evaluator on the written moments CSV (must be bit-identical: %.17g
// round trip), plus the blocked-input and shift-table paths.
//   strip_eval_inmemory_test KERNEL_DIR KERNEL_ID WORKDIR
// Writes WORKDIR/inmem.moments.csv (for the Python reference) and WORKDIR/inmem/ (in-memory result, replay format).
#include <cstring>
#include <filesystem>
#include <iostream>
#include <random>

#include "../K2AStripEvaluator.hh"

namespace {
int gFail = 0;
void Check(bool ok, const std::string& m) { std::cout << (ok ? "  ok   " : "  FAIL ") << m << "\n"; if (!ok) ++gFail; }
bool Same(const K2AStrip::Result& a, const K2AStrip::Result& b) {
  auto eq = [](const std::vector<double>& x, const std::vector<double>& y) {
    return x.size() == y.size() && (x.empty() || std::memcmp(x.data(), y.data(), x.size() * sizeof(double)) == 0); };
  return a.status == b.status && a.U.strips == b.U.strips && a.V.strips == b.V.strips && eq(a.U.ion, b.U.ion) && eq(a.V.ion, b.V.ion) &&
         eq(a.U.electron, b.U.electron) && eq(a.V.electron, b.V.electron);
}
K2AInline::Event Births(bool insideAnode) {
  std::mt19937_64 g(424242);
  std::uniform_real_distribution<double> U(0, 1);
  const K2AInline::MuonFrame f{0.0123, 1.23, 0.456};
  K2AInline::MomentAccumulator acc(f);
  const double wcen[] = {0.0, 0.52, -0.61, 0.95, -1.1, 0.03};   // primary w offsets (cm): D ~ 10 w + 0.44 mm -> m = 0, +3, -3, +6, -6, 0
  for (int p = 0; p < 6; ++p) {
    acc.BeginPrimary(2 * p + 1);
    const int n = 3000 + 500 * p;
    for (int i = 0; i < n; ++i) {
      const int k = (i % 3) - 1 + (p % 2);                          // wires -1..1 or 0..2: multiwire primaries
      const double r = 0.003 * std::pow(0.1 / 0.003, U(g)), ph = 6.283185307179586 * U(g);
      acc.AddCallbackBirth(0.3 * k + r * std::cos(ph), r * std::sin(ph), wcen[p] + 0.004 * (U(g) - .5), 10 + 30 * p + 3 * U(g), std::uint64_t(i));
    }
    if (insideAnode && p == 2) { acc.AddCallbackBirth(0.3 + 0.0006, 0.0008, wcen[p], 50, std::uint64_t(n)); acc.EndPrimary(std::uint64_t(n + 1)); }
    else acc.EndPrimary(std::uint64_t(n));
  }
  return acc.Finish();
}
}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc != 4) throw std::runtime_error("usage: strip_eval_inmemory_test KERNEL_DIR KERNEL_ID WORKDIR");
    const std::string W = argv[3];
    if (std::filesystem::exists(W)) throw std::runtime_error("WORKDIR exists");
    std::filesystem::create_directories(W);
    const auto k = K2AStrip::KernelTable::Open(argv[1], argv[2]);
    K2AStrip::Evaluator ev(*k, {});
    const auto e = Births(false);
    Check(e.UsableForCathodeSignal(), "synthetic full-birth event usable (" + K2AInline::StatusNames(e.status) + ")");
    std::set<int> wires; for (const auto& c : e.clouds) wires.insert(c.wire);
    Check(wires.size() >= 3 && e.clouds.size() >= 12, "multiwire event: " + std::to_string(e.clouds.size()) + " clouds on " + std::to_string(wires.size()) + " wires");
    K2AInline::WriteOutputs(W + "/inmem", e);
    const auto rm = ev.Evaluate(e);
    const auto rc = ev.Evaluate(K2AStrip::ReadMomentsCsv(W + "/inmem.moments.csv"));
    Check(rm.usable && rm.status == K2AStrip::OK, "in-memory evaluation usable, status " + K2AStrip::StatusNames(rm.status));
    Check(Same(rm, rc), "in-memory Event == moments-CSV evaluation, bit for bit");
    std::set<std::pair<int, int>> ms; for (const auto& s : rm.shifts) ms.insert({s.m_U, s.m_V});
    int pos = 0, neg = 0; for (const auto& s : rm.shifts) { pos += s.m_U > 0 || s.m_V > 0; neg += s.m_U < 0 || s.m_V < 0; }
    Check(pos > 0 && neg > 0 && rm.groups.size() >= 4, "positive and negative displacements: " + std::to_string(rm.groups.size()) + " groups");
    std::map<int, std::pair<int, int>> table; for (const auto& s : rm.shifts) table[s.primary] = {s.m_U, s.m_V};
    const auto rt = ev.Evaluate(e, &table);
    Check(Same(rm, rt) && rt.shifts_differing_from_computed == 0, "imposed shift table equal to the rule -> identical result");
    auto t2 = table; t2.begin()->second.first += 1;
    const auto rx = ev.Evaluate(e, &t2);
    Check(rx.usable && rx.shifts_differing_from_computed == 1, "one imposed shift differing from the rule is reported");
    // provenance: every cloud appears in exactly one group, groups' charge adds up
    double q = 0; int nc = 0; for (const auto& g : rm.groups) { q += g.N_fC; nc += g.clouds; }
    double qc = 0; for (const auto& c : e.clouds) qc += c.N_fC;
    Check(nc == int(e.clouds.size()) && std::abs(q - qc) <= 1e-12 * qc, "group provenance covers every cloud once; charge adds up");
    std::filesystem::create_directories(W + "/inmem_result");
    K2AStrip::WriteResult(W + "/inmem_result", rm, ev.Edges());
    // blocked input: the accumulator's INSIDE_ANODE block must stop the evaluation
    const auto eb = Births(true);
    const auto rb = ev.Evaluate(eb);
    Check(!eb.UsableForCathodeSignal() && !rb.usable && (rb.status & K2AStrip::INPUT_BLOCKED) && rb.U.strips.empty() &&
          rb.blocked_charge_fC == rb.charge_fC, "blocked moment event -> INPUT_BLOCKED, no waveform, all charge reported blocked");
    std::cout << (gFail ? "FAIL" : "PASS") << " in-memory strip evaluation test (SYNTHETIC)\n";
    return gFail ? 1 : 0;
  } catch (const std::exception& x) { std::cerr << "ERROR " << x.what() << "\n"; return 2; }
}
