# 200 ns TW/TC stereo-u outlier (20 development muons, ion only): initial analysis

Source: the 20-muon `t_ladder_summary.txt` (ion term only, provisional CR-RC^3 stand-in). The original result is kept
as reported: TW (and TC, identical) minus T0 at Tp = 200 ns gives u mean +0.01132 mm, rms 0.04871 mm,
max |du| = 0.2177 mm, against max |du| = 0.0087 mm at 100 ns.

## What the summary alone already shows (arithmetic, not yet confirmed per muon)
* Stereo reconstruction: u = 1.7 mm x (c_V - c_U) / (2 x 0.1) = 8.5 mm per strip of centroid difference, and
  w = 0.85 mm x (c_U + c_V).
* At 200 ns, max |d c_U| = 0.02533 strips -> 8.5 x 0.02533 = 0.215 mm; max |d c_V| = 0.00073 strips -> 0.006 mm.
  Max |dw| = 0.02128 mm ~= 0.85 x 0.0253. So the 0.2177 mm u outlier is almost entirely one U-plane centroid jump.
* U strip-hit flips at 200 ns occur only at f = 0.02 (one flip in one muon: mean 0.05, max 1). The centroid uses
  exactly the strips with peak >= 0.02 x reference lead, so this flip is a strip entering or leaving the centroid.
* At 100 ns the only U flip was at f = 0.1, not at the inclusion threshold, and the U centroid moved by <= 0.0017
  strips: no membership change, no outlier.
* Size check: a strip with peak p ~ 0.02 L at d ~ 3-4 strips from the centroid, in a plane sum S ~ 2.5-3.5 L, moves
  the centroid by p d / S ~ 0.02-0.03 strips, i.e. 0.17-0.25 mm of u. This matches.

**Hypothesis (to confirm):** the 200 ns outlier is a discontinuity of the threshold-dependent centroid definition. A
peripheral U strip sits within a fraction of a percent of the 2% inclusion threshold, and the ~0.2% TW amplitude
change moves it across. It is not a waveform-level effect of TW; the waveform errors stay at the 0.35% level.

## Confirmation (WSL, no new simulation)
`t_shaped_compare.py` from the commit that adds this note writes per-strip peaks for every variant, the membership changes at
the inclusion threshold, and u, w under alternative centroids (membership fixed from the reference, lead +-1,
lead +-2, all positive strips). `t_centroid_diagnostic.py` turns one compare directory into the report.
```bash
cd ~/na6p/wt-validation && git pull
TA=test/mwpc-stage-c/cathode-ion-reference/timing-ablation
D=<the 20-muon T-ladder run dir>                     # contains muon_*/T0,TW,TC/replay/waveforms.csv
python3 $TA/t_shaped_compare.py $D --variants T0,TW,TC --tp 100,200 --no-prompt --out $D/compare_centroid_diag_$(date +%H%M%S)
python3 $TA/t_centroid_diagnostic.py <that compare dir> --tp 200 --variant TW --abs-thr-fC 0.5,1,2,5
python3 $TA/t_centroid_diagnostic.py <that compare dir> --tp 100 --variant TW --abs-thr-fC 0.5,1,2,5
```
This re-runs only the shaping and comparison of existing waveforms. Return the two `centroid_diagnostic_*.md` files
and the compare `t_ladder_summary.txt`. The report gives (1) the muon, (2) every strip's reference and TW peak and
its distance from the 2% threshold, with the membership change marked, (3) u, w with fixed vs threshold-dependent
membership, (4) u, w under lead +-1, +-2, all-positive and absolute (fC-equivalent) thresholds.
The absolute thresholds are in units of the unit-peak stand-in shaper (1 fC delta -> 1); they become real
electronics thresholds only with the VMM gain and transfer function.

## Why it matters for stereo-angle and pitch optimization
The u lever arm is pitch / (2 tan alpha): 8.5 mm per strip now, larger for smaller angles. Any hard strip-membership
rule then converts a sub-percent amplitude change into a 0.1-0.2 mm u jump. Detector-optimization studies must
therefore use the actual cluster algorithm and electronics threshold, or a centroid without hard membership, and
quote u-resolution with the threshold-crossing population included. Agreement in waveform L2 alone does not
guarantee agreement in reconstructed u.
