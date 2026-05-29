#!/usr/bin/env python3
"""
belgium_hail_map.py
===================
Custom GFS-derived hail-SIZE map for Belgium (output in cm).

This is OUR formula, not SHIP and not Groenemeijer's. Structure:

    size_surface = MELT(size_aloft, freezing level)             # survival
                 * Capacity(buoyancy)                            # the gate
                 * [ORG_FLOOR + (1-ORG_FLOOR) * Efficiency(shear)]  # trim

Design decisions baked in (from our discussion):
  - Buoyancy is the GATE and the primary axis (capacity). Saturates at the
    top so big CAPE doesn't run away linearly.
  - Organization (deep-layer shear) is a bounded EFFICIENCY multiplier with a
    high floor -- it trims, it never zeroes the result and never dominates.
    Floor is set high (0.55) because we only care in severe-ish environments
    where some shear is already present.
  - Melt survival is size-dependent (big stones survive a deep warm layer,
    small ones melt away) and is the live discriminator within the
    favourable box.
  - Plotted from 1 cm, with the visual breaks placed on the operational
    thresholds: 2 cm (severe, ESSL/ESTOFEX) and 5 cm (very large).

!!! THE COEFFICIENTS BELOW ARE UNCALIBRATED FIRST GUESSES. !!!
They give physically plausible behaviour but have NOT been tuned against
ESWD reports. Treat the absolute numbers with suspicion until calibrated;
the field structure is the trustworthy part for now.

Run target: 00Z GFS, forecast hour 15 (valid 15Z) of the current UTC day.
Override with env vars RUN_HOUR / VALID_HOUR / RUN_DATE (YYYY-MM-DD) if needed.

Data: NOAA GFS 0.25 deg (pgrb2.0p25) via `herbie`.
Output: output/belgium_hail_<date>_<run>z_f<fxx>.png  (+ copy to output/latest.png)
"""

import os
from datetime import datetime, timezone

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from herbie import Herbie

# --------------------------------------------------------------------------
# Domain
# --------------------------------------------------------------------------
LON_MIN, LON_MAX = 1.5, 7.5          # GFS lon is 0..360; Belgium stays positive
LAT_MIN, LAT_MAX = 48.5, 52.5
PLOT_EXTENT = [2.4, 6.5, 49.4, 51.6]  # what the map actually shows

# --------------------------------------------------------------------------
# FORMULA COEFFICIENTS  -- all uncalibrated first guesses, tune these
# --------------------------------------------------------------------------
SIZE_MAX_CM  = 8.0     # max "aloft" size an extreme environment can reach
LI_SCALE     = 4.0     # deg C; growth-zone buoyancy (Best LI) saturates on this
LR_MIN       = 5.5     # deg C/km; below this, mid-level lapse rate adds nothing
LR_REF       = 8.0     # deg C/km; at/above this, full lapse-rate credit
LR_WEIGHT    = 0.5     # capacity floor from buoyancy alone (0..1); rest from lapse rate

SHEAR_SCALE  = 18.0    # m/s; organization efficiency saturates on this scale
ORG_FLOOR    = 0.55    # disorganized storms still realize this fraction

MELT_C       = 0.40    # melting strength; higher = more melt for a given setup

PLOT_FLOOR_CM = 1.0    # blank anything below this on the map

# --------------------------------------------------------------------------
# Herbie field helpers (proven search strings)
# --------------------------------------------------------------------------
def _ds(H, search):
    ds = H.xarray(search)
    if isinstance(ds, list):
        ds = ds[0]
    return ds


def _first(ds):
    return ds[list(ds.data_vars)[0]]


def _best_li(H):
    """
    Best (4-layer) Lifted Index from GFS -- the most-unstable parcel's buoyancy
    referenced to 500 hPa (~ -12 to -18 C over Belgium in convective season),
    i.e. buoyancy weighted toward the hail GROWTH ZONE. This is the GFS-native
    analogue of Groenemeijer's MU_LI_10 and a far better hail-capacity proxy
    than MUCAPE, which is inflated by low-level (sub-growth-zone) buoyancy.
    Negative = unstable. Falls back to surface LFTX if 4LFTX is absent.
    """
    for s in (":4LFTX:180-0 mb above ground:",
              ":4LFTX:",
              ":LFTX:500-1000 mb:",
              ":LFTX:"):
        try:
            return _first(_ds(H, s))
        except Exception:
            continue
    raise RuntimeError("No Lifted Index field found in GFS file.")


# --------------------------------------------------------------------------
# Pull fields + compute the hail-size field
# --------------------------------------------------------------------------
def compute(run, fxx):
    H = Herbie(run.strftime("%Y-%m-%d %H:%M"), model="gfs",
               product="pgrb2.0p25", fxx=fxx)
    if H.grib is None:
        raise RuntimeError(
            f"GFS {run:%Y-%m-%d %H}Z f{fxx:03d} not available yet. "
            f"If you ran this too soon after 00Z, wait until ~05-06Z.")

    li_da = _best_li(H)
    lat = li_da.latitude.values
    lon = li_da.longitude.values
    ila = np.where((lat >= LAT_MIN) & (lat <= LAT_MAX))[0]
    ilo = np.where((lon >= LON_MIN) & (lon <= LON_MAX))[0]

    def C(a):
        return np.asarray(a)[np.ix_(ila, ilo)]

    best_li = C(li_da.values)                                    # deg C, neg = unstable

    t = _first(_ds(H, ":TMP:(500|700) mb:"))
    t500 = C(t.sel(isobaricInhPa=500).values) - 273.15           # deg C
    t700 = C(t.sel(isobaricInhPa=700).values) - 273.15
    gh = _first(_ds(H, ":HGT:(500|700) mb:"))
    z500 = C(gh.sel(isobaricInhPa=500).values)                   # m
    z700 = C(gh.sel(isobaricInhPa=700).values)
    lr75 = (t700 - t500) / (z500 - z700) * 1000.0                # deg C/km

    u10 = C(_first(_ds(H, ":UGRD:10 m above ground:")).values)
    v10 = C(_first(_ds(H, ":VGRD:10 m above ground:")).values)
    u500 = C(_first(_ds(H, ":UGRD:500 mb:")).values)
    v500 = C(_first(_ds(H, ":VGRD:500 mb:")).values)
    shear06 = np.hypot(u500 - u10, v500 - v10)                   # m/s (~0-6 km proxy)

    fz_msl = C(_first(_ds(H, ":HGT:0C isotherm:")).values)       # m MSL
    orog = C(_first(_ds(H, ":HGT:surface:")).values)             # m MSL
    fz_agl = np.maximum(fz_msl - orog, 0.0)                      # m AGL

    # ---- THE FORMULA --------------------------------------------------
    # Capacity: growth-zone buoyancy (Best LI). Gate CLOSES where the column
    # is stable to 500 hPa (LI >= 0); saturates as LI gets strongly negative.
    # Emphasised by mid-level lapse rate.
    buoy_term = 1.0 - np.exp(np.minimum(best_li, 0.0) / LI_SCALE)  # 0..1, 0 if LI>=0
    lr_term = np.clip((lr75 - LR_MIN) / (LR_REF - LR_MIN), 0.0, 1.0)
    capacity = SIZE_MAX_CM * buoy_term * (LR_WEIGHT + (1 - LR_WEIGHT) * lr_term)

    # Efficiency: organization via deep-layer shear, floored + saturating.
    eff = ORG_FLOOR + (1 - ORG_FLOOR) * np.tanh(shear06 / SHEAR_SCALE)

    size_aloft = capacity * eff                                  # cm, pre-melt

    # Melt survival: deeper warm layer melts more; bigger stones survive.
    melt_frac = np.exp(-MELT_C * (fz_agl / 1000.0) / np.maximum(size_aloft, 0.1))
    size_sfc = size_aloft * melt_frac                            # cm at ground

    lat_c = lat[ila]
    lon_c = lon[ilo]
    if lat_c[0] > lat_c[-1]:                 # contourf wants ascending lat
        lat_c = lat_c[::-1]
        size_sfc = size_sfc[::-1, :]
    return size_sfc, lat_c, lon_c


# --------------------------------------------------------------------------
# Plot
# --------------------------------------------------------------------------
def plot(size, lat, lon, run, fxx, outpath):
    levels = [1.0, 1.5, 2.0, 3.0, 5.0, 7.0, 12.0]
    colors = ["#9ecae1", "#41b6c4",        # 1-1.5, 1.5-2  : cool/muted (sub-severe)
              "#fed976", "#fd8d3c",        # 2-3,   3-5    : warm (severe)
              "#e31a1c", "#7a0177"]        # 5-7,   7+     : hot/extreme (very large)
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(levels, cmap.N)

    masked = np.where(size >= PLOT_FLOOR_CM, size, np.nan)

    fig = plt.figure(figsize=(9, 8))
    ax = plt.axes(projection=ccrs.Mercator())   # conformal: fixes the lat squish
    ax.set_extent(PLOT_EXTENT, crs=ccrs.PlateCarree())

    cf = ax.contourf(lon, lat, masked, levels=levels, cmap=cmap, norm=norm,
                     extend="neither", transform=ccrs.PlateCarree())

    # operational threshold contours so the severe read pops
    try:
        cl = ax.contour(lon, lat, size, levels=[2.0, 5.0], colors="black",
                        linewidths=[1.2, 1.6], transform=ccrs.PlateCarree())
        ax.clabel(cl, fmt={2.0: "2 cm", 5.0: "5 cm"}, fontsize=8)
    except Exception:
        pass

    ax.add_feature(cfeature.OCEAN, facecolor="#dfe7ef", zorder=0)
    ax.add_feature(cfeature.BORDERS, linewidth=0.9, edgecolor="black")
    ax.add_feature(cfeature.COASTLINE, linewidth=0.6)
    gl = ax.gridlines(draw_labels=True, linewidth=0.3, color="gray", alpha=0.5)
    gl.top_labels = gl.right_labels = False

    try:
        from cartopy.io import shapereader
        shp = shapereader.natural_earth(resolution="10m", category="cultural",
                                        name="admin_0_countries")
        for rec in shapereader.Reader(shp).records():
            if rec.attributes.get("NAME") == "Belgium":
                ax.add_geometries([rec.geometry], ccrs.PlateCarree(),
                                  facecolor="none", edgecolor="red",
                                  linewidth=1.8, zorder=5)
    except Exception as e:
        print("  (Belgium outline skipped:", e, ")")

    cb = plt.colorbar(cf, ax=ax, shrink=0.85, pad=0.02, ticks=levels)
    cb.set_label("Estimated max hail size at ground [cm]")

    import datetime as _dt
    valid = run + _dt.timedelta(hours=fxx)
    ax.set_title(
        "Belgium - GFS experimental hail size (OUR formula, v0.1 UNCALIBRATED)\n"
        f"Run {run:%Y-%m-%d %H}Z   valid {valid:%Y-%m-%d %H}Z   (f{fxx:03d})",
        fontsize=12)
    ax.text(0.5, -0.13,
            "Buoyancy = Best Lifted Index (growth-zone).  Cool 1-2 cm | warm "
            "2-5 cm (severe) | hot >=5 cm.  Coefficients are first guesses.",
            transform=ax.transAxes, ha="center", fontsize=8, style="italic")

    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    plt.savefig(outpath, dpi=150, bbox_inches="tight")
    print("Saved", outpath)


# --------------------------------------------------------------------------
def main():
    import datetime as _dt
    run_hour = int(os.environ.get("RUN_HOUR", "0"))
    valid_hour = int(os.environ.get("VALID_HOUR", "15"))
    date_str = os.environ.get("RUN_DATE")  # YYYY-MM-DD, optional
    if date_str:
        d = _dt.datetime.strptime(date_str, "%Y-%m-%d").date()
    else:
        d = _dt.datetime.now(timezone.utc).date()

    run = _dt.datetime(d.year, d.month, d.day, run_hour, tzinfo=timezone.utc)
    fxx = valid_hour - run_hour
    if fxx < 0:
        fxx += 24
    run = run.replace(tzinfo=None)  # Herbie wants naive UTC

    print(f"Target: GFS {run:%Y-%m-%d %H}Z  f{fxx:03d}  (valid {valid_hour:02d}Z)")
    size, lat, lon = compute(run, fxx)

    fname = f"belgium_hail_{run:%Y%m%d}_{run.hour:02d}z_f{fxx:03d}.png"
    out = os.path.join("output", fname)
    plot(size, lat, lon, run, fxx, out)

    # also keep a stable filename for easy phone viewing
    import shutil
    shutil.copyfile(out, os.path.join("output", "latest.png"))
    print("Updated output/latest.png")


if __name__ == "__main__":
    main()
