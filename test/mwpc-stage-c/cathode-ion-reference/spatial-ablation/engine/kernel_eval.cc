// Compiled evaluation of the K-model response (same basis, corners and time interpolation as the validated
// per-cloud loop of s_kernel.py / CathodeResponseEngine._loop). Plain C ABI for ctypes; no Python headers needed.
//   g++ -O3 -std=c++17 -shared -fPIC -o libkernel_eval.so kernel_eval.cc
// Order of operations per (side, cloud, strip): for each of the 4 (r, D) corners, Q[t] += w * sum_b C[b] K[b,..,t]
// (basis sum in index order), then np.interp semantics (slope form, left = 0, right = Q[last]) at edges - t_cloud,
// then the bin difference is added to the strip; the electron step P is added in bin floor(t / dt).
#include <cmath>
#include <cstdint>
#include <vector>

extern "C" int kernel_eval(
    const float* K, const double* KP,            // K (nbasisAll, 2, nwire, nr, nD, ntau), KP (nbasisAll, 2, nwire, nr, nD)
    int nbasisAll, int nwire, int nr, int nD, int ntau,
    const int* bidx, int nb,                     // basis rows used by the model
    const double* tau,                           // ntau
    int ncl, const double* C,                    // (ncl, nb) coefficients x charge
    const int* wk, const int* jr, const double* fr,
    const int* jdU, const double* fdU, const int* jdV, const double* fdV,
    const double* tcl,                           // cloud times
    int nO, const int* offs, int dsub,
    int nE, const double* edges, double dt,
    double* ion, double* ele) {                  // (nE-1, 2, nO), accumulated (caller zeroes); unit: fC per bin
  (void)nbasisAll;
  const long sB = 2L * nwire * nr * nD * ntau, sS = 1L * nwire * nr * nD * ntau, sW = 1L * nr * nD * ntau,
             sR = 1L * nD * ntau, sD = ntau;
  const long pB = 2L * nwire * nr * nD, pS = 1L * nwire * nr * nD, pW = 1L * nr * nD, pR = nD;
  std::vector<double> Q(ntau), qe(nE);
  const int nbins = nE - 1;
  for (int si = 0; si < 2; ++si) {
    const int* jd = si ? jdV : jdU; const double* fd = si ? fdV : fdU;
    for (int c = 0; c < ncl; ++c) {
      const double* Cc = C + long(c) * nb;
      const int b0 = int(std::floor(tcl[c] / dt));
      for (int oi = 0; oi < nO; ++oi) {
        const int d0 = jd[c] - dsub * offs[oi];
        if (d0 < 0 || d0 + 1 >= nD) return 1;
        for (int t = 0; t < ntau; ++t) Q[t] = 0.;
        double P = 0.;
        for (int dr = 0; dr < 2; ++dr) {
          const double wr = dr ? fr[c] : 1 - fr[c];
          for (int dd = 0; dd < 2; ++dd) {
            const double w = wr * (dd ? fd[c] : 1 - fd[c]);
            const long base = si * sS + long(wk[c]) * sW + long(jr[c] + dr) * sR + long(d0 + dd) * sD;
            const long pbase = si * pS + long(wk[c]) * pW + long(jr[c] + dr) * pR + (d0 + dd);
            for (int t = 0; t < ntau; ++t) {
              double s = 0.;
              for (int b = 0; b < nb; ++b) s += Cc[b] * double(K[bidx[b] * sB + base + t]);
              Q[t] += w * s;
            }
            double p = 0.;
            for (int b = 0; b < nb; ++b) p += Cc[b] * KP[bidx[b] * pB + pbase];
            P += w * p;
          }
        }
        int j = 0;
        for (int e = 0; e < nE; ++e) {                 // edges increasing: monotone pointer, np.interp semantics
          const double x = edges[e] - tcl[c];
          if (x < tau[0]) { qe[e] = 0.; continue; }
          if (x > tau[ntau - 1]) { qe[e] = Q[ntau - 1]; continue; }
          while (j < ntau - 2 && tau[j + 1] <= x) ++j;
          if (x == tau[ntau - 1]) { qe[e] = Q[ntau - 1]; continue; }
          const double slope = (Q[j + 1] - Q[j]) / (tau[j + 1] - tau[j]);
          qe[e] = slope * (x - tau[j]) + Q[j];
        }
        for (int e = 0; e < nbins; ++e) ion[(long(e) * 2 + si) * nO + oi] += qe[e + 1] - qe[e];
        if (b0 < nbins && b0 >= 0) ele[(long(b0) * 2 + si) * nO + oi] += P;
      }
    }
  }
  return 0;
}
