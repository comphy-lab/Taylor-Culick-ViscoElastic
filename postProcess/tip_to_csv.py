#!/usr/bin/env python3
"""Convert a planar Taylor-Culick tip log into a tip-velocity CSV.

Reads the ``c<CaseNo>-tip.dat`` written by the in-code tip diagnostic of
``simulationCases/TaylorCulickPlanar.c`` and emits

    t,tstar,x_tip,v_tip,v_over_VTC

``x_tip`` is the interface position on the midplane (minimum x over the
reconstructed VOF facets in the midplane band).  ``v_tip`` is d(x_tip)/dt,
obtained from a local first-order least-squares (Savitzky-Golay) slope over a
window of ``--window`` samples, which is far less noisy than a bare two-point
difference at the sub-cell scale of the VOF reconstruction.

Non-dimensionalisation: rho = sigma = h0 = 1 with h0 the FULL sheet thickness
gives V_TC = sqrt(2); the tip log header records the run's actual
V_TC = sqrt(2*sigma/(rho1*h0)), which this script parses rather than
assuming rho1 = h0 = 1.  Savva & Bush (JFM 626, 2009) use the half-thickness,
so a case quoted at Oh is run with mu = sqrt(2) * Oh; their viscous time
is tau_vis = mu/2 and their reduced time is t* = t / tau_vis.

Usage
-----
    python3 postProcess/tip_to_csv.py <case-dir-or-tipfile> [-o out.csv] [--window 21]
"""

from __future__ import annotations

import argparse
import glob
import math
import os
import sys

import numpy as np

FALLBACK_V_TC = math.sqrt(2.0)


def find_tip_file(path: str) -> str:
    """Resolve a directory or explicit file to a tip data file."""
    if os.path.isfile(path):
        return path
    hits = sorted(glob.glob(os.path.join(path, "c*-tip.dat")))
    if not hits:
        raise SystemExit(f"no c*-tip.dat found under {path}")
    if len(hits) > 1:
        raise SystemExit(f"ambiguous tip files under {path}: {hits}")
    return hits[0]


def read_tip(fname: str):
    """Return (t, x_tip, x_tip_vof, tau_vis, v_tc) from a tip log.

    ``v_tc`` is ``None`` if the header predates the ``V_TC`` field, in which
    case the caller must fall back to the rho1 = h0 = 1 default.
    """
    tau_vis = None
    v_tc = None
    with open(fname) as fh:
        for line in fh:
            if line.startswith("#") and "tau_vis" in line and "V_TC" in line:
                parts = line.replace("#", "").split()
                tau_vis = float(parts[parts.index("tau_vis") + 1])
                v_tc = float(parts[parts.index("V_TC") + 1])
                break
    data = np.loadtxt(fname, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[0] < 2:
        raise SystemExit(f"{fname}: fewer than two samples; nothing to do")
    t = data[:, 0]
    x_tip = data[:, 2]
    x_vof = data[:, 3]
    if tau_vis is None:
        tstar = data[:, 1]
        with np.errstate(divide="ignore", invalid="ignore"):
            tau_vis = float(np.nanmedian(np.where(tstar > 0, t / tstar, np.nan)))
    return t, x_tip, x_vof, tau_vis, v_tc


def savgol_slope(t: np.ndarray, x: np.ndarray, window: int) -> np.ndarray:
    """Local linear least-squares slope over a centred window of samples.

    Falls back to a shorter, one-sided window near the ends so the output has
    the same length as the input.
    """
    n = len(t)
    half = max(1, window // 2)
    v = np.empty(n)
    for k in range(n):
        lo = max(0, k - half)
        hi = min(n, k + half + 1)
        if hi - lo < 2:
            lo, hi = max(0, k - 1), min(n, k + 2)
        tt = t[lo:hi]
        xx = x[lo:hi]
        tt = tt - tt.mean()
        denom = float(np.dot(tt, tt))
        v[k] = np.dot(tt, xx - xx.mean()) / denom if denom > 0 else np.nan
    return v


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("path", help="case directory or c<N>-tip.dat file")
    ap.add_argument("-o", "--output", default=None,
                    help="output CSV (default: <case-dir>/tip_velocity.csv)")
    ap.add_argument("--window", type=int, default=21,
                    help="samples in the local slope fit (default: 21)")
    args = ap.parse_args()

    if args.window < 3:
        raise SystemExit("--window must be at least 3")

    tip_file = find_tip_file(args.path)
    case_dir = os.path.dirname(os.path.abspath(tip_file))
    out = args.output or os.path.join(case_dir, "tip_velocity.csv")

    t, x_tip, x_vof, tau_vis, v_tc = read_tip(tip_file)
    if v_tc is None:
        print(f"warning: {tip_file}: no V_TC in header, "
              f"falling back to rho1 = h0 = 1 default ({FALLBACK_V_TC:.6g})",
              file=sys.stderr)
        v_tc = FALLBACK_V_TC

    # Drop duplicated times (a restart can repeat the last recorded sample).
    keep = np.concatenate(([True], np.diff(t) > 0))
    t, x_tip, x_vof = t[keep], x_tip[keep], x_vof[keep]

    v_tip = savgol_slope(t, x_tip, args.window)
    tstar = t / tau_vis if tau_vis and tau_vis > 0 else np.zeros_like(t)

    with open(out, "w") as fh:
        fh.write("t,tstar,x_tip,v_tip,v_over_VTC\n")
        for k in range(len(t)):
            fh.write("%.8e,%.8e,%.8e,%.8e,%.8e\n"
                     % (t[k], tstar[k], x_tip[k], v_tip[k], v_tip[k] / v_tc))

    dev = np.max(np.abs(x_tip - x_vof))
    print(f"{tip_file} -> {out}")
    print(f"  samples {len(t)}  t in [{t[0]:.4g}, {t[-1]:.4g}]  "
          f"tau_vis {tau_vis:.6g}  V_TC {v_tc:.6g}")
    print(f"  x_tip final {x_tip[-1]:.6g}  v_tip final {v_tip[-1]:.6g}  "
          f"v/V_TC final {v_tip[-1] / v_tc:.6g}")
    print(f"  max |x_tip - x_tip_vof| = {dev:.3g} "
          f"(large values flag rim pinch-off or midplane bubble entrainment)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
