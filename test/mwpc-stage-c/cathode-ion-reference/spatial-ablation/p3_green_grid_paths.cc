// Green's-function ion paths for the spatial (S-ladder) ablation.
//
// Drifts single Ar+ ions from a polar grid of birth points around chosen anodes with the FROZEN transport code:
// the frozen exporter source p3_export_full_spectral_paths.cc is compiled into this tool unchanged (its main() is
// renamed by the preprocessor), so gas, 9-wire analytic field, 2100 V, cathodes, sensor area, RKF step limit and the
// path validation are byte-identical to the T0/T1 reference exports. Only the birth points differ: a grid instead of
// archived ions. Output is the exporter's format (ions.csv, vertices.f64, complete.txt) plus grid.csv, so the frozen
// replay's read_export()/path_for() validate it unchanged.
//
//   p3_green_grid_paths MOBILITY NEW_DIR --wires -1,0,1 --r-um 16,18,...,120 --ntheta 72
//   p3_green_grid_paths MOBILITY NEW_DIR --points FILE     (validation: start points u_mm,v_mm,t_ns per line, e.g.
//                                                           taken from a frozen export; paths must reproduce it)
//
// Grid: wire k at u = k * 3 mm (patch frame), radius r from the wire centre, azimuth theta measured from +u towards
// +v, theta_m = (m + 0.5) * 360/ntheta (never exactly in the v = 0 plane, where the frozen code requires v != 0).
// w = 0, t = 0, weight 1 per row. Strip columns (wU, wV) are written as 0; the table builder adds the strip phase.
#define main p3_frozen_export_main_unused
#include "p3_export_full_spectral_paths.cc"
#undef main

#include <sstream>

struct Start { int wire; double r_um, theta_deg, u_cm, v_cm, t_ns; };

namespace {
std::vector<double> ParseList(const std::string& s) {
  std::vector<double> v; std::stringstream ss(s); std::string item;
  while (std::getline(ss, item, ',')) v.push_back(std::stod(item));
  return v;
}
}

int main(int argc, char** argv) {
  try {
    Require(argc >= 4, "usage: p3_green_grid_paths MOBILITY NEW_DIR --wires a,b,.. --r-um r1,r2,.. --ntheta N");
    namespace fs = std::filesystem;
    const std::string mobility = argv[1];
    const fs::path out = argv[2];
    std::vector<double> wires{-1, 0, 1}, rum;
    int ntheta = 72;
    std::string pointsFile;
    for (int i = 3; i < argc; ++i) {
      const std::string a = argv[i];
      Require(i + 1 < argc, "missing value for " + a);
      if (a == "--wires") wires = ParseList(argv[++i]);
      else if (a == "--r-um") rum = ParseList(argv[++i]);
      else if (a == "--ntheta") ntheta = std::stoi(argv[++i]);
      else if (a == "--points") pointsFile = argv[++i];
      else Require(false, "unknown option " + a);
    }
    std::vector<Start> starts;
    if (!pointsFile.empty()) {
      std::ifstream in(pointsFile); Require(bool(in), "cannot read " + pointsFile);
      std::string line; std::getline(in, line);   // header
      while (std::getline(in, line)) {
        if (line.empty()) continue;
        const auto x = ParseList(line); Require(x.size() == 3, "points line must be u_mm,v_mm,t_ns");
        const int k = int(std::lround(x[0] / 3.));
        const double du = x[0] - 3. * k;
        starts.push_back({k, 1000. * std::hypot(du, x[1]), std::atan2(x[1], du) * 180 / M_PI, x[0] / 10., x[1] / 10., x[2]});
      }
      Require(!starts.empty(), "no points");
    } else {
      Require(!rum.empty() && std::is_sorted(rum.begin(), rum.end()), "need an increasing --r-um list");
      Require(rum.front() > 1e4 * NA60Plane::wireRadius && rum.back() < 1000., "radii must be in (wire radius, 1 mm)");
      Require(ntheta >= 8 && ntheta <= 720 && ntheta % 4 == 0, "ntheta must be a multiple of 4 in [8,720]");
      for (double k : wires) Require(k == std::round(k) && std::abs(k) <= 3, "wires must be integers in [-3,3]");
      for (double k : wires)
        for (double r : rum)
          for (int m = 0; m < ntheta; ++m) {
            const double th = (m + 0.5) * 2 * M_PI / ntheta;
            starts.push_back({int(k), r, th * 180 / M_PI, k * NA60Plane::pitch + 1e-4 * r * std::cos(th),
                              1e-4 * r * std::sin(th), 0.});
          }
    }
    Require(!fs::exists(out), "output exists");
    Require(fs::create_directories(out), "cannot create output directory");
    PhysicalIonPaths transport(mobility, .005);   // same step as the frozen exporter
    std::ofstream data(out / "vertices.f64", std::ios::binary), meta(out / "ions.csv"), grid(out / "grid.csv");
    Require(bool(data) && bool(meta) && bool(grid), "cannot create output files");
    meta << std::setprecision(17)
         << "row,primary,weight,birth_ns,end_ns,first_vertex,vertices,upper_index,lower_index\n";
    grid << std::setprecision(17) << "row,wire,r_um,theta_deg,u_mm,v_mm,end_u_mm,end_v_mm,end_ns\n";
    long long offset = 0, row = 0;
    for (const auto& st : starts) {
      NA60Plane::Archive::Birth b{};
      b.u = st.u_cm; b.v = st.v_cm; b.w = 0.; b.t = st.t_ns; b.weight = 1.; b.primary = int(row);
      auto path = transport.Run(b);
      path.points.back().v = std::copysign(NA60Plane::gap, path.end.v);   // as in the frozen exporter
      meta << row << "," << row << ",1," << b.t << "," << path.end.t << "," << offset << "," << path.points.size() << ",0,0\n";
      grid << row << "," << st.wire << "," << st.r_um << "," << st.theta_deg << "," << 10 * b.u << "," << 10 * b.v << ","
           << 10 * path.end.u << "," << 10 * path.end.v << "," << path.end.t << "\n";
      for (const auto& p : path.points) {
        Require(std::abs(p.w - b.w) < 1e-8, "unexpected longitudinal motion");
        const double rowv[5] = {p.t, 10 * p.u, 10 * p.v, 0., 0.};
        data.write(reinterpret_cast<const char*>(rowv), sizeof(rowv));
      }
      offset += path.points.size();
      ++row;
    }
    data.close(); meta.close(); grid.close();
    Require(bool(data) && bool(meta) && bool(grid), "output I/O failed");
    Require(fs::file_size(out / "vertices.f64") == static_cast<unsigned long long>(offset) * 40, "binary size mismatch");
    std::ofstream seal(out / "complete.txt");
    seal << std::setprecision(17) << "rows=" << row << "\nexact_ions=" << row << "\nweight_sum=" << row
         << "\nvertices=" << offset << "\nstep_mm=0.005\n"
#ifdef FROZEN_EXPORTER_SHA256
         << "frozen_exporter_source_sha256=" << FROZEN_EXPORTER_SHA256 << "\n"
#endif
         << "mobility=" << mobility << "\n";
    seal.close(); Require(bool(seal), "seal write failed");
    if (pointsFile.empty())
      std::cout << "PASS green grid rows=" << row << " vertices=" << offset << " wires=" << wires.size()
                << " radii=" << rum.size() << " ntheta=" << ntheta << "\n";
    else std::cout << "PASS green points rows=" << row << " vertices=" << offset << "\n";
    return 0;
  } catch (const std::exception& e) { std::cerr << e.what() << "\n"; return 2; }
}
