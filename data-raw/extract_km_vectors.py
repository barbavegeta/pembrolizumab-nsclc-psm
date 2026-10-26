"""Extract Kaplan-Meier step coordinates from the vector artwork of Figure 2
in Reck et al., J Clin Oncol 2021;39(21):2339-2349 (PMC8280089).

Figure 2 is embedded in the article PDF as vector paths. Each arm's KM
curve is a single stroked path in the arm's colour, so its step corners can
be read exactly from the path coordinates, without pixel tracing. Axis
calibration uses the tick marks drawn on each axis, fitted by least squares.

The PDF itself is not stored in the repository (CC BY-NC-ND 4.0). Download
it from the PMC Open Access dataset:

    https://pmc-oa-opendata.s3.amazonaws.com/PMC8280089.1/PMC8280089.1.pdf
    (md5 4b2d15722a5c6f1228124ce57c1c78f2)

Usage:
    python extract_km_vectors.py PMC8280089.1.pdf OUTPUT_DIR

Requires PyMuPDF (tested with 2.0) and numpy.
"""

import csv
import hashlib
import sys
from pathlib import Path

import numpy as np
import pymupdf

EXPECTED_MD5 = "4b2d15722a5c6f1228124ce57c1c78f2"
PAGE_INDEX = 5  # page 6 of the PDF
PEMBRO_RGB = (0.0, 0.463, 0.659)
CHEMO_RGB = (0.671, 0.086, 0.173)

# Panel layout on page 6: vertical extent of each panel's plotting area (pt),
# and the x-axis range printed on the figure.
PANELS = {
    "os": {"y_range": (60.0, 200.0), "x_max_months": 72},
    "pfs": {"y_range": (270.0, 420.0), "x_max_months": 66},
}


def rgb(d):
    c = d.get("color")
    return tuple(round(v, 3) for v in c) if c else None


def in_panel(d, panel):
    lo, hi = PANELS[panel]["y_range"]
    r = d["rect"]
    return r.y0 >= lo and r.y1 <= hi


def tick_positions(drawings, panel, axis):
    """Positions of the tick marks on the x or y axis of a panel."""
    black = [d for d in drawings if rgb(d) == (0.0, 0.0, 0.0) and in_panel(d, panel)]
    out = []
    for d in black:
        for it in d["items"]:
            if it[0] != "l":
                continue
            p, q = it[1], it[2]
            length = abs(q.x - p.x) + abs(q.y - p.y)
            if length > 5:  # axis lines and the reference line in the legend
                continue
            if axis == "x" and abs(p.x - q.x) < 1e-6:
                out.append(p.x)
            if axis == "y" and abs(p.y - q.y) < 1e-6:
                out.append(p.y)
    return np.unique(np.round(out, 3))


def calibrate(drawings, panel):
    """Least-squares maps from page coordinates to months and to survival."""
    xt = tick_positions(drawings, panel, "x")
    yt = tick_positions(drawings, panel, "y")
    x_max = PANELS[panel]["x_max_months"]
    if len(xt) != x_max:
        raise ValueError(f"{panel}: expected {x_max} x ticks (1 per month), found {len(xt)}")
    if len(yt) != 10:
        raise ValueError(f"{panel}: expected 10 y ticks (10%-100%), found {len(yt)}")
    months = np.arange(1, x_max + 1)
    bx = np.polyfit(xt, months, 1)
    surv = np.arange(1.0, 0.0, -0.1)  # ticks run top (100%) to bottom (10%)
    by = np.polyfit(np.sort(yt), surv, 1)
    resid_x = np.max(np.abs(np.polyval(bx, xt) - months))
    resid_y = np.max(np.abs(np.polyval(by, np.sort(yt)) - surv))
    return bx, by, resid_x, resid_y


def curve_path(drawings, panel, colour):
    """The KM step path: the arm's widest-stroke path in the panel."""
    cands = [d for d in drawings if rgb(d) == colour and in_panel(d, panel)]
    cands.sort(key=lambda d: (d.get("width") or 0, len(d["items"])), reverse=True)
    return cands[0], cands[1:]


def path_vertices(d):
    """Ordered vertices of a path made of connected line segments."""
    pts = []
    for it in d["items"]:
        if it[0] != "l":
            raise ValueError(f"unexpected path item {it[0]}")
        p, q = it[1], it[2]
        if pts and (abs(pts[-1][0] - p.x) > 1e-3 or abs(pts[-1][1] - p.y) > 1e-3):
            raise ValueError("path is not continuous")
        if not pts:
            pts.append((p.x, p.y))
        pts.append((q.x, q.y))
    return pts


def main(pdf_path, out_dir):
    data = Path(pdf_path).read_bytes()
    md5 = hashlib.md5(data).hexdigest()
    if md5 != EXPECTED_MD5:
        raise SystemExit(f"PDF md5 {md5} does not match the expected {EXPECTED_MD5}")
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    drawings = pymupdf.open(pdf_path)[PAGE_INDEX].get_drawings()

    log = []
    for panel in PANELS:
        bx, by, rx, ry = calibrate(drawings, panel)
        log.append(f"{panel}: calibration max residual {rx:.4f} months, {ry:.5f} survival")
        for arm, colour in (("pembro", PEMBRO_RGB), ("chemo", CHEMO_RGB)):
            d, others = curve_path(drawings, panel, colour)
            verts = path_vertices(d)
            rows = [(float(np.polyval(bx, x)), float(np.polyval(by, y))) for x, y in verts]
            with open(out / f"{panel}_{arm}.csv", "w", newline="") as f:
                w = csv.writer(f)
                w.writerow(["time_months", "surv"])
                for t, s in rows:
                    w.writerow([f"{t:.4f}", f"{s:.5f}"])
            n_ticks = sum(len(o["items"]) for o in others)
            log.append(f"{panel}_{arm}: {len(rows)} vertices, "
                       f"last t = {rows[-1][0]:.2f} months, {n_ticks} censoring tick marks")

    with open(out / "axis_limits.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["endpoint", "x_max_months"])
        for panel, spec in PANELS.items():
            w.writerow([panel, spec["x_max_months"]])

    print("\n".join(log))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
