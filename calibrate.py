#!/usr/bin/env python3
"""
calibrate.py
============
Anchor-day calibration harness for OUR Belgium hail-size formula.

It does NOT re-implement the formula. It imports `compute()` and the
coefficient constants straight from belgium_hail_map.py, so the numbers it
prints are exactly what the live map would produce. Tune the coefficients in
belgium_hail_map.py, re-run this, read the table.

What it does per anchor day:
  - computes the DAILY-MAX hail field over a convective window (default
    12/15/18 UTC of the 00Z run) -- not a single snapshot. ESWD reports a
    day's *maximum* stone, so we compare day-max model to day-max obs and
    remove the "wrong forecast hour" confound.
  - masks to the Belgium polygon (falls back to the full box if Natural Earth
    isn't available) and reports:
        BE-max   : largest predicted size in Belgium  (compare to obs stone)
        BE-p99   : 99th percentile (is the peak broad or a single hot cell?)
        %>=2cm   : fraction of Belgium at/over severe -- the BASE-RATE check
        %>=1cm   : fraction at/over the plotting floor
  - flags each day against its expectation (big / severe / marginal /
    control / null), so leaks are called out in words.

The truth column (observed_max_cm) is YOURS to fill from ESWD -- a handful of
numbers. Leave it blank and the row still computes; you just won't get an
error figure for it.

Usage:  python calibrate.py            (reads anchor_days.csv)
Output: printed table + a clean results.txt
"""

import csv
import datetime as dt

import numpy as np

import belgium_hail_map as M   # <-- single source of truth for formula + coeffs

# --------------------------------------------------------------------------
# Config
# --------------------------------------------------------------------------
RUN_HOUR    = 0                 # GFS cycle to evaluate (00Z, as on the map)
VALID_HOURS = [12, 15, 18]      # daily-max window, UTC (convective afternoon)
ANCHOR_CSV  = "anchor_days.csv"
RESULTS_TXT = "results.txt"
SEVERE_CM   = 2.0               # ESSL/ESTOFEX severe-hail threshold


# --------------------------------------------------------------------------
# Output helper: collect lines so we can print AND write a clean results.txt
# (keeps Herbie's download chatter out of the saved file)
# --------------------------------------------------------------------------
_LINES = []
def out(s=""):
    print(s)
    _LINES.append(s)


# --------------------------------------------------------------------------
# Belgium polygon mask over the GFS grid (built once; same grid every date)
# --------------------------------------------------------------------------
def belgium_mask(lat, lon):
    try:
        from cartopy.io import shapereader
        from shapely.geometry import Point
        from shapely.prepared import prep
        shp = shapereader.natural_earth(resolution="10m", category="cultural",
                                        name="admin_0_countries")
        geom = None
        for rec in shapereader.Reader(shp).records():
            if rec.attributes.get("NAME") == "Belgium":
                geom = rec.geometry
                break
        if geom is None:
            return None
        pg = prep(geom)
        m = np.zeros((len(lat), len(lon)), dtype=bool)
        for i, la in enumerate(lat):
            for j, lo in enumerate(lon):
                if pg.contains(Point(float(lo), float(la))):
                    m[i, j] = True
        return m if m.any() else None
    except Exception as e:
        out(f"  (Belgium mask unavailable: {e}; using full domain box instead)")
        return None


# --------------------------------------------------------------------------
# Daily-max field for one date
# --------------------------------------------------------------------------
def daily_max(date_str):
    d = dt.datetime.strptime(date_str, "%Y-%m-%d")
    run = dt.datetime(d.year, d.month, d.day, RUN_HOUR)
    stack = lat = lon = None
    for vh in VALID_HOURS:
        fxx = vh - RUN_HOUR
        size, lat, lon = M.compute(run, fxx)
        stack = size if stack is None else np.maximum(stack, size)
    return stack, lat, lon


# --------------------------------------------------------------------------
# Metrics + verdict
# --------------------------------------------------------------------------
def metrics(size, mask):
    vals = size[mask] if mask is not None else size.ravel()
    return {
        "bmax": float(np.nanmax(vals)),
        "b99":  float(np.nanpercentile(vals, 99)),
        "f2":   float(np.mean(vals >= SEVERE_CM)),
        "f1":   float(np.mean(vals >= 1.0)),
    }


def verdict(expect, m):
    """Heuristic flags -- a starting read, not gospel. bmax in cm, f2 a fraction."""
    bmax, f2 = m["bmax"], m["f2"]
    e = (expect or "").strip().lower()
    if e == "null":
        return "GATE LEAK (null day not ~0)" if bmax >= 2.0 else "ok - quiet"
    if e == "marginal":
        if bmax >= 5.0 or f2 > 0.20:
            return "TOO HOT for a marginal day"
        return "ok"
    if e == "control":      # organized/sheared but little big hail expected
        return "SHEAR LEAK (organization carrying it)" if bmax >= 4.0 else "ok"
    if e == "severe":
        return "too weak" if bmax < 2.0 else "ok"
    if e == "big":
        return "too weak" if bmax < 3.0 else "ok"
    return ""


# --------------------------------------------------------------------------
def main():
    # --- read anchors --------------------------------------------------
    rows = []
    with open(ANCHOR_CSV, newline="") as f:
        for r in csv.DictReader(f):
            if not r.get("date"):
                continue
            obs = (r.get("observed_max_cm") or "").strip()
            rows.append({
                "date": r["date"].strip(),
                "expect": (r.get("expect") or "").strip(),
                "obs": float(obs) if obs else None,
                "note": (r.get("note") or "").strip(),
            })

    # --- coefficient banner -------------------------------------------
    out("=" * 78)
    out("COEFFICIENTS IN USE  (from belgium_hail_map.py)")
    out(f"  SIZE_MAX_CM = {M.SIZE_MAX_CM}   LI_SCALE = {M.LI_SCALE}   "
        f"LR_MIN/REF = {M.LR_MIN}/{M.LR_REF}   LR_WEIGHT = {M.LR_WEIGHT}")
    out(f"  SHEAR_SCALE = {M.SHEAR_SCALE}   ORG_FLOOR = {M.ORG_FLOOR}   "
        f"MELT_C = {M.MELT_C}")
    out(f"  Daily-max window (UTC): {VALID_HOURS}   "
        f"Run: {RUN_HOUR:02d}Z   Severe threshold: {SEVERE_CM} cm")
    out("=" * 78)

    # --- header --------------------------------------------------------
    hdr = (f"{'date':<11}{'expect':<10}{'obs':>5} {'BE-max':>7}{'BE-p99':>7}"
           f"{'%>=2':>6}{'%>=1':>6}{'err':>6}  flag")
    out(hdr)
    out("-" * 78)

    mask = None
    mask_built = False
    region = "full box"
    flagged = []

    for row in rows:
        try:
            size, lat, lon = daily_max(row["date"])
        except Exception as e:
            out(f"{row['date']:<11}{row['expect']:<10}{'--':>5}  "
                f"FETCH FAILED: {type(e).__name__}")
            continue

        if not mask_built:
            mask = belgium_mask(lat, lon)
            region = "Belgium polygon" if mask is not None else "full box"
            mask_built = True

        m = metrics(size, mask)
        v = verdict(row["expect"], m)
        obs_s = f"{row['obs']:.1f}" if row["obs"] is not None else "--"
        err_s = (f"{m['bmax'] - row['obs']:+.1f}"
                 if row["obs"] is not None else "--")

        out(f"{row['date']:<11}{row['expect']:<10}{obs_s:>5} "
            f"{m['bmax']:>7.1f}{m['b99']:>7.1f}"
            f"{m['f2']*100:>5.0f}%{m['f1']*100:>5.0f}%{err_s:>6}  {v}")

        if v and not v.startswith("ok"):
            flagged.append((row, v))

    out("-" * 78)
    out(f"Region for stats: {region}.  "
        f"BE-max vs obs = did we reach the stone; err = BE-max - obs.")

    # --- interpretation ------------------------------------------------
    out("")
    out("READ:")
    out("  * Base rate: even on BIG days %>=2 over Belgium should be MODEST.")
    out("    If ordinary/marginal/null days show large %>=2, the gain is too")
    out("    high or a gate is leaking -- pull SIZE_MAX_CM down / LI_SCALE up")
    out("    / MELT_C up before trusting any single-case 'hit'.")
    out("  * err > 0 on big days = over-forecast peak; err < 0 = under.")
    out("  * 'control' day (organized, low hail) going hot = shear is carrying")
    out("    the result -> raise ORG_FLOOR's ceiling effect / re-check coupling.")
    if flagged:
        out("")
        out("FLAGGED:")
        for row, v in flagged:
            out(f"  - {row['date']} ({row['expect']}): {v}"
                + (f"  | {row['note']}" if row["note"] else ""))
    else:
        out("")
        out("No flags -- but fill in observed_max_cm before trusting that.")

    with open(RESULTS_TXT, "w") as f:
        f.write("\n".join(_LINES) + "\n")
    print(f"\nWrote {RESULTS_TXT}")


if __name__ == "__main__":
    main()
