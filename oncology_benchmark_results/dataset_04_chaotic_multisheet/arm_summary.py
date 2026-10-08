"""Summarize the 'Arm-CN wide' sheet (rows = arms, columns = patients, values = arm-level total CN)
into one row per patient. Gain/loss is relative to round(median arm CN) of that patient."""
import sys, csv, statistics
sys.path.insert(0, sys.argv[3])
from tableio import read_table
t = read_table(sys.argv[1], sheet="Arm-CN wide")
arms = t.column("chr-arm")
out = csv.writer(open(sys.argv[2], "w", newline=""))
out.writerow(["patient", "arm_CN_median", "arms_gained", "arms_lost", "gained_arms", "lost_arms"])
for p in t.columns[1:]:
    vals = [float(v) for v in t.column(p)]
    base = round(statistics.median(vals))
    g = [a for a, v in zip(arms, vals) if v >= base + 0.5]
    l = [a for a, v in zip(arms, vals) if v <= base - 0.5]
    out.writerow([p, round(statistics.median(vals), 2), len(g), len(l), " ".join(g), " ".join(l)])
