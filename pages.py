"""Highcharts page and CSV helpers shared with StormSense Latest (copied from its products.py)."""
import csv
import json
import os
import re
from datetime import datetime, timezone

import numpy as np

DATA = re.compile(r"var data = (\[.*?\]);\n", re.S)
UTC = re.compile(r"Date\.UTC\((\d+),(\d+),(\d+),(\d+),(\d+),(\d+)(?:,[\d.]+)?\)")
GUIDE = re.compile(r"<script>var ssGuidance = .*?;</script>\n?")


# --- Highcharts pages --------------------------------------------------------
def read_page(path):
    html = open(path, encoding="utf-8").read()
    m = DATA.search(html)
    if not m:  # a truncated or empty Tidewatch page: leave it as it is
        print(f"  no chart data in {os.path.basename(path)}")
        return html, []
    js = UTC.sub(lambda u: str(int(datetime(int(u[1]), int(u[2]) + 1, int(u[3]), int(u[4]), int(u[5]), int(u[6]),
                                            tzinfo=timezone.utc).timestamp() * 1000)), m.group(1))
    return html, json.loads(js)




def write_page(path, html, series, title=None, guidance=None):
    if guidance:  # read by the styler, which writes it under the disclaimer
        html = GUIDE.sub("", html).replace("</head>", "<script>var ssGuidance = %s;</script>\n</head>"
                                           % json.dumps(guidance), 1)
    html = DATA.sub(lambda m: "var data = " + json.dumps(series, separators=(",", ":")) + ";\n", html, count=1)
    if title:
        html = re.sub(r'("title": \{"text": ")[^"]*(")', lambda m: m.group(1) + title + m.group(2), html, count=1)
    open(path, "w", encoding="utf-8").write(html)


def ms(t):
    return int(t.timestamp() * 1000)


def pts(times, vals):
    return [[ms(t), round(float(v), 3)] for t, v in zip(times, vals) if not np.isnan(v)]


# --- CSV ---------------------------------------------------------------------
def read_csv(path):
    with open(path, newline="") as f:
        rows = list(csv.reader(f))
    return rows[0], rows[1:]


def write_csv(path, head, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(head)
        w.writerows(rows)


def fmt(v):
    return "" if v is None or np.isnan(v) else "%.3f" % (round(v, 3) + 0.0)


def iso(t):
    return t.strftime("%Y-%m-%dT%H:%M:%SZ")


def parse_t(s):
    return datetime.strptime(s, "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=timezone.utc)


def offset(f36_dir, st):
    """Tidewatch's MLLW minus NAVD for the station (median over its f36 pair)."""
    try:
        a = dict(r[:4:3] for r in read_csv(os.path.join(f36_dir, f"{st}_2018_f36_MLLW.csv"))[1])
        b = dict(r[:4:3] for r in read_csv(os.path.join(f36_dir, f"{st}_2018_f36_NAVD.csv"))[1])
        d = sorted(float(a[t]) - float(b[t]) for t in a.keys() & b.keys() if a[t] and b[t])
        return d[len(d) // 2] if d else None
    except Exception:
        return None


