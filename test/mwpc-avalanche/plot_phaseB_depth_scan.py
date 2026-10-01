#!/usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

RUNS = [
    (-0.2, "phaseB_weighting_2mm_m0p2.csv"),
    (-0.5, "phaseB_weighting_2mm_m0p5.csv"),
    (-1.0, "phaseB_weighting_2mm_m1p0.csv"),
    (-1.5, "phaseB_weighting_2mm_m1p5.csv"),
    (-1.8, "phaseB_weighting_2mm_m1p8.csv"),
]

STRIP_INDICES = np.arange(-8, 9)
STRIP_CENTERS_MM = STRIP_INDICES * 1.7


def strip_column(k: int) -> str:
    if k < 0:
        return f"strip_m{-k}"
    if k > 0:
        return f"strip_p{k}"
    return "strip_0"


def main() -> None:
    plt.figure(figsize=(8, 5))

    summary_rows = []

    for v_mm, filename in RUNS:
        df = pd.read_csv(filename)

        # Use the scan point closest to w = 0.
        row = df.iloc[(df["w_mm"].abs()).argmin()]

        strip_cols = [strip_column(k) for k in STRIP_INDICES]
        phi = row[strip_cols].to_numpy(dtype=float)

        total = phi.sum()
        if total <= 0:
            raise RuntimeError(f"Non-positive weighting-potential sum in {filename}")

        frac = phi / total

        mean_w = np.sum(frac * STRIP_CENTERS_MM)
        rms_w = np.sqrt(np.sum(frac * (STRIP_CENTERS_MM - mean_w) ** 2))

        summary_rows.append(
            {
                "v_mm": v_mm,
                "sum_phi": total,
                "central_fraction": float(frac[STRIP_INDICES == 0][0]),
                "profile_rms_mm": float(rms_w),
            }
        )

        plt.plot(
            STRIP_CENTERS_MM,
            frac,
            marker="o",
            label=f"v = {v_mm:.1f} mm",
        )

    plt.xlabel("Strip center w [mm]")
    plt.ylabel("Normalized weighting-potential share")
    plt.title("2 mm cathode: strip sharing vs depth")
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig("phaseB_2mm_depth_scan.png", dpi=200)
    plt.close()

    summary = pd.DataFrame(summary_rows)
    summary.to_csv("phaseB_2mm_depth_scan_summary.csv", index=False)

    print(summary.to_string(index=False))
    print("\nWrote:")
    print("  phaseB_2mm_depth_scan.png")
    print("  phaseB_2mm_depth_scan_summary.csv")


if __name__ == "__main__":
    main()
