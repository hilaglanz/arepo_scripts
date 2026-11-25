"""Read saved flux CSVs; write a comparison figure and reproducible statistics."""
from pathlib import Path
import csv
import hashlib
import json

import numpy as np
from reportlab.graphics import renderSVG
from reportlab.graphics.charts.lineplots import LinePlot
from reportlab.graphics.shapes import Drawing, Line, Rect, String
from reportlab.lib import colors

SOURCE = Path(r"C:\Users\hila\Dropbox\RMHD\IDORT\AGB\1.5")
OUTPUT = Path(__file__).resolve().parent
FIELDS = ["outward_Msun_per_yr", "inward_Msun_per_yr", "net_Msun_per_yr"]
FILES = ["bowen_flux.csv", "no_bowen_flux_291_295.csv", "no_bowen_flux_344_496.csv"]
data = {}
provenance = {}
for filename in FILES:
    path = SOURCE / filename
    with path.open(newline="", encoding="utf-8-sig") as handle:
        data[filename] = [{key: float(value) for key, value in row.items()}
                          for row in csv.DictReader(handle)]
    rows = data[filename]
    assert all(np.isfinite(list(row.values())).all() for row in rows)
    assert len({(row["snapshot"], row["inner_Rsun"]) for row in rows}) == len(rows)
    assert all(row[FIELDS[0]] >= 0 and row[FIELDS[1]] >= 0 for row in rows)
    assert all(np.isclose(row[FIELDS[0]] - row[FIELDS[1]], row[FIELDS[2]],
                          rtol=1e-13, atol=1e-25) for row in rows)
    provenance[filename] = {"path": str(path), "rows": len(rows),
                            "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def series(filename, inner):
    rows = sorted((row for row in data[filename] if row["inner_Rsun"] == inner),
                  key=lambda row: row["time_yr"])
    time = np.array([row["time_yr"] for row in rows])
    values = np.array([[row[field] for field in FIELDS] for row in rows])
    assert np.all(np.diff(time) > 0)
    return time, values


def interpolate(time, values, requested):
    assert requested[0] >= time[0] and requested[-1] <= time[-1]
    return np.column_stack([np.interp(requested, time, values[:, i]) for i in range(3)])


def average(time, values):
    return (np.diff(time)[:, None] * (values[1:] + values[:-1]) / 2).sum(axis=0) / (time[-1] - time[0])


statistics = []
for inner in (1000, 1500, 2000):
    bt, bv = series("bowen_flux.csv", inner)
    nt, nv = series("no_bowen_flux_344_496.csv", inner)
    start, end = max(bt[0], nt[0]), min(bt[-1], nt[-1])
    for name, begin in (("shared later interval", start), ("last year", end - 1)):
        grid = np.r_[begin, nt[(nt > begin) & (nt < end)], end]
        bmean = average(grid, interpolate(bt, bv, grid))
        nmean = average(grid, interpolate(nt, nv, grid))
        bfull = np.r_[begin, bt[(bt > begin) & (bt < end)], end]
        full_mean = average(bfull, interpolate(bt, bv, bfull))
        statistics.append({"inner_Rsun": inner, "period": name,
                           "time_yr": [float(begin), float(end)],
                           "bowen_mean": dict(zip(FIELDS, bmean.tolist())),
                           "no_bowen_mean": dict(zip(FIELDS, nmean.tolist())),
                           "bowen_full_cadence_mean": dict(zip(FIELDS, full_mean.tolist()))})

(OUTPUT / "statistics.json").write_text(json.dumps({
    "sources": provenance,
    "method": "Linear interpolation of Bowen onto no-Bowen physical times; trapezoidal time weighting. No interpolation across missing snapshots 295-344. Rates are finite-shell transport averages, not an escape criterion.",
    "statistics": statistics,
}, indent=2), encoding="utf-8")

# Standard ReportLab line plots: show only the common, continuously sampled segment.
width, height = 850, 755
drawing = Drawing(width, height)
drawing.add(Rect(0, 0, width, height, fillColor=colors.white, strokeColor=None))
drawing.add(String(55, 723, "Net mass transport with and without Bowen opacity",
                   fontName="Helvetica-Bold", fontSize=18))
drawing.add(String(55, 701, "Shared later interval: 7.14-10.29 yr; positive = outward, negative = inward",
                   fontName="Helvetica", fontSize=11, fillColor=colors.HexColor("#444444")))
blue, orange = colors.HexColor("#0072B2"), colors.HexColor("#D55E00")
for x, color, label in ((55, blue, "Tables + Bowen"), (270, orange, "Tables only")):
    drawing.add(Line(x, 678, x + 25, 678, strokeColor=color, strokeWidth=2.5))
    drawing.add(String(x + 33, 674, label, fontName="Helvetica", fontSize=11))

plot_left, plot_width, plot_height = 92, 693, 144
limits = [(-16, 8, [-15, -10, -5, 0, 5]), (-4, 1, [-4, -3, -2, -1, 0, 1]), (0, 13, [0, 3, 6, 9, 12])]
for index, (inner, outer) in enumerate(((1000, 1500), (1500, 2000), (2000, 3000))):
    bottom = 487 - index * 190
    bt, bv = series("bowen_flux.csv", inner)
    nt, nv = series("no_bowen_flux_344_496.csv", inner)
    begin, end = nt[0], min(nt[-1], bt[-1])
    bgrid = np.r_[begin, bt[(bt > begin) & (bt < end)], end]
    bnet = interpolate(bt, bv, bgrid)[:, 2] * 1e6
    nmask = (nt >= begin) & (nt <= end)
    plot = LinePlot()
    plot.x, plot.y = plot_left, bottom
    plot.width, plot.height = plot_width, plot_height
    plot.data = [list(zip(bgrid.tolist(), bnet.tolist())),
                 list(zip(nt[nmask].tolist(), (nv[nmask, 2] * 1e6).tolist()))]
    plot.lines[0].strokeColor, plot.lines[1].strokeColor = blue, orange
    plot.lines[0].strokeWidth = plot.lines[1].strokeWidth = 2
    plot.xValueAxis.valueMin, plot.xValueAxis.valueMax = 7.1, 10.35
    plot.xValueAxis.valueSteps = [7.2, 7.8, 8.4, 9.0, 9.6, 10.2]
    plot.xValueAxis.labelTextFormat = "%.1f"
    low, high, ticks = limits[index]
    plot.yValueAxis.valueMin, plot.yValueAxis.valueMax = low, high
    plot.yValueAxis.valueSteps = ticks
    plot.yValueAxis.labelTextFormat = "%g"
    for axis in (plot.xValueAxis, plot.yValueAxis):
        axis.labels.fontName = "Helvetica"
        axis.labels.fontSize = 10
        axis.strokeColor = colors.HexColor("#777777")
        axis.gridStrokeColor = colors.HexColor("#E8E8E8")
        axis.visibleGrid = 1
    zero_height = bottom + (-low) / (high - low) * plot_height
    drawing.add(plot)
    drawing.add(Line(plot_left, zero_height, plot_left + plot_width, zero_height,
                     strokeColor=colors.HexColor("#444444"), strokeWidth=0.6,
                     strokeDashArray=[3, 3]))
    drawing.add(String(plot_left, bottom + plot_height + 12,
                       f"{inner}-{outer} solar radii", fontName="Helvetica-Bold", fontSize=12))
    drawing.add(String(25, bottom + plot_height / 2 - 28,
                       "Net flux", fontName="Helvetica", fontSize=10, angle=90))

drawing.add(String(438, 65, "Simulation time [yr]", textAnchor="middle",
                   fontName="Helvetica", fontSize=11))
drawing.add(String(55, 36, "All net-flux axes: 10^-6 solar masses / yr. Shells follow the core.",
                   fontName="Helvetica", fontSize=10, fillColor=colors.HexColor("#444444")))
drawing.add(String(55, 19, "Lines join saved measurements; the missing earlier no-Bowen interval is excluded.",
                   fontName="Helvetica", fontSize=10, fillColor=colors.HexColor("#444444")))
renderSVG.drawToFile(drawing, str(OUTPUT / "net_flux_comparison.svg"))
print(json.dumps(statistics, indent=2))
