// SYNTHETIC test driver for BirthMoments.hh (no ROOT): births.csv (primary,u_cm,v_cm,w_cm,t_ns,weight) + frame.csv
// (patch_shift_cm,u_in_mm,w_in_mm) -> moments CSV.  Used by test_birth_moments.sh only.
#include <fstream>
#include <iostream>
#include <sstream>
#include "../BirthMoments.hh"
int main(int argc, char** argv) {
  try {
    BirthMoments::Require(argc == 4, "usage: birth_moments_core_test births.csv frame.csv out.csv");
    std::ifstream b(argv[1]), f(argv[2]);
    std::string line; std::getline(b, line);
    std::vector<BirthMoments::Birth> v;
    while (std::getline(b, line)) {
      if (line.empty()) continue;
      std::stringstream ss(line); std::string x; std::vector<std::string> c;
      while (std::getline(ss, x, ',')) c.push_back(x);
      v.push_back({std::stoi(c[0]), std::stod(c[1]), std::stod(c[2]), std::stod(c[3]), std::stod(c[4]), std::stod(c[5])});
    }
    std::getline(f, line); std::getline(f, line);
    std::stringstream ss(line); std::string x; std::vector<double> fr;
    while (std::getline(ss, x, ',')) fr.push_back(std::stod(x));
    const auto r = BirthMoments::Compute(v, {fr.at(0), fr.at(1), fr.at(2)});
    std::ofstream o(argv[3]); BirthMoments::Write(o, r);
    std::cout << "core: rows=" << r.rows << " clouds=" << r.clouds.size() << "\n";
    return 0;
  } catch (const std::exception& e) { std::cerr << e.what() << "\n"; return 2; }
}
