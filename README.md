# Belgium hail-size map (experimental GFS formula)

Generates a hail-**size** map (in cm) over Belgium from the GFS 00Z run,
valid 15Z, using a custom formula. **All compute happens on GitHub's
servers** — your phone only taps a button and views the result.

> **Status: v0.1, coefficients UNCALIBRATED.** The field *structure* is
> physically reasonable; the absolute centimetre values are first guesses
> until tuned against ESWD reports. See "Calibration" below.

## The formula (what it does)

```
size_surface = melt_survival(size_aloft, freezing level)     # survival fraction
             × capacity(buoyancy)                            # the gate, saturating
             × [floor + (1-floor)·efficiency(deep-layer shear)]   # bounded trim
```

- **Capacity** — saturating in MUCAPE, emphasised by 700–500 hPa lapse rate.
  Buoyancy is the gate and the primary axis.
- **Efficiency** — organisation from ~0–6 km bulk shear, with a high floor
  (0.55) so it trims but never zeroes or dominates.
- **Melt survival** — size-dependent: small stones melt out under a deep warm
  layer, big stones survive. The live discriminator in favourable setups.

Plotted from **1 cm**, with visual breaks on the operational thresholds
**2 cm** (severe) and **5 cm** (very large), per ESSL/ESTOFEX definitions.

## One-time setup (doable on a phone)

1. Create a new GitHub repo (e.g. `belgium-hail-map`), public is fine.
2. Add these four files, keeping the paths exactly:
   - `belgium_hail_map.py`
   - `requirements.txt`
   - `README.md`
   - `.github/workflows/hail-map.yml`  ← note the `.github/workflows/` folder.
     On the GitHub web/app, "Create new file" and type the full path with
     slashes; it makes the folders for you.
3. Go to the **Actions** tab and enable workflows if prompted.

## Running it

- **Manual (the button):** Actions → "Belgium hail map" → **Run workflow**.
  Leave the date blank for today, or type a `YYYY-MM-DD` to backfill.
- **Automatic:** it also runs every day at 06:00 UTC (00Z GFS f015 is ready
  well before then).

## Viewing the map on your phone

After a run finishes (1–3 min), open **`output/latest.png`** in the repo —
GitHub renders PNGs inline, so it just shows. Each run also keeps a dated
copy (`output/belgium_hail_YYYYMMDD_00z_f015.png`) and uploads the PNG as a
downloadable **artifact** on the run's summary page as a backup.

## Calibration (the real next step)

The coefficients live at the top of `belgium_hail_map.py`, each labelled.
To tune honestly:

- Score only the **≥2 cm** field against ESWD reports — sub-2 cm hail is
  badly under-reported, so don't trust or tune to the 1–2 cm band.
- Adjust `CAPE_SCALE`, `SIZE_MAX_CM`, `SHEAR_SCALE`, `ORG_FLOOR`, `MELT_C`
  to remove bias, rather than chasing individual cases.

## Known simplifications (hooks for later)

- Capacity uses MUCAPE × lapse rate as a stand-in for true growth-zone
  (≥ −10 °C) buoyancy / LI@−10. Swapping in a real MU-parcel ascent (MetPy)
  is the main physics upgrade.
- Organisation uses bulk shear only; SRH / storm-relative inflow would
  sharpen the mesocyclone proxy.
- Melt uses dry-bulb freezing level; wet-bulb zero + sub-cloud RH is better.
