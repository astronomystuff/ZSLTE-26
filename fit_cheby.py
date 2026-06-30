import numpy as np
from numpy.polynomial.chebyshev import Chebyshev
import json

# -----------------------------
# Load raw REBOUND output
# -----------------------------
raw = np.load("zslte_raw_xyz.npz")
jd = raw["jd"]
x = raw["x"]
y = raw["y"]
z = raw["z"]

# -----------------------------
# Slice into segments
# -----------------------------
def slice_segments(jd, x, y, z, window_days=32):
    segments = []
    N = len(jd)
    step = window_days

    for start in range(0, N, step):
        end = start + step
        if end > N:
            break

        seg = {
            "jd0": float(jd[start]),
            "jd1": float(jd[end-1]),
            "t": jd[start:end],
            "x": x[start:end],
            "y": y[start:end],
            "z": z[start:end]
        }
        segments.append(seg)

    return segments

segments_raw = slice_segments(jd, x, y, z, window_days=32)

# -----------------------------
# Fit Chebyshev
# -----------------------------
def fit_cheby_segment(seg, degree=12):
    jd0 = seg["jd0"]
    jd1 = seg["jd1"]

    u = 2 * (seg["t"] - jd0) / (jd1 - jd0) - 1

    cx = Chebyshev.fit(u, seg["x"], degree).coef.tolist()
    cy = Chebyshev.fit(u, seg["y"], degree).coef.tolist()
    cz = Chebyshev.fit(u, seg["z"], degree).coef.tolist()

    return {
        "jd0": jd0,
        "jd1": jd1,
        "cx": cx,
        "cy": cy,
        "cz": cz
    }

segments_fitted = [fit_cheby_segment(seg) for seg in segments_raw]

# -----------------------------
# Save kernel
# -----------------------------
kernel = {
    "version": "ZSLTE-26",
    "segments": segments_fitted
}

with open("ZSLTE-26.json", "w") as f:
    json.dump(kernel, f)

print("Saved zslte-2026.json")
