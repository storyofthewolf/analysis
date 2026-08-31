# Design Note — A Unified ExoCAM Reading Primitive

**Status: design note, not a commitment.** This records an idea for possible future
work. Nothing here is implemented, and none of it is scheduled. The underlying
infrastructure work was explicitly deferred as not aligned with immediate job
duties; the note exists so the idea is not lost.

*Written August 30, 2026. Claims below were checked against the code in
`exocam-trend`, `exocam-analysis`, and `exocam-spectral` as of that date.*

---

## The problem — three independent readers

Three repositories each independently solve the same problem — "read N netCDF
history files into arrays" — and each arrived at a different set of arbitrary
limitations.

### `exocam-trend`

The timeseries package. Its reading path is `trend_core.read_monthly_files`,
driven by `trend.py`.

- **Monthly `h0` cadence only.** `trend.py` states it directly at the top:
  *"Presently it only works with monthly mean outputs (nhtfrq=0)."*
  `read_monthly_files` globs `{root_path}/{case_id}{prefix}*.nc` and filters with
  the regex `\.(\d{4})-(\d{2})\.nc$` — a file that is not stamped `YYYY-MM` is
  skipped outright. The three component prefixes are hardcoded in `trend.py`
  as `.cam.h0.`, `.clm2.h0.`, and `.cice.h.`.
- **2D variables only.** `read_monthly_files` reads `ncid.variables[vname][0, :, :]`
  — first time index, all lats and lons. There is no path for a `(lev, lat, lon)`
  field. `vars.in` documents this as a deliberate, dated deferral:
  *"As of June 1, 2023. only 2D variables are supported. […] I will eventually add
  option for read in full 3D variables when the utility is requested or becomes
  self-evident."*
- **Variable selection through a relative-path `vars.in`.**
  `trend_utils.read_request_var()` (line 32) is literally `with open('vars.in', 'r')`
  — a relative path, no argument, no search path, no fallback. `trend.py` must be
  run from the repository root. The file has three blocks (read / print / plot),
  each three lines (atm / ice / lnd), with print and plot validated as subsets of
  read; the token `energy` is exempt because it is derived rather than read.

  The consequence: **a bespoke variable set means editing a file that is shared,
  version-controlled repo state.** Two analyses wanting different variables cannot
  run concurrently from the same checkout — one silently gets the other's variable
  list. `run_trend_batch.sh` acknowledges this rather than solving it, warning that
  it *"does NOT edit vars.in; verify it before the run."*

**Correction to the original framing:** trend is *not* archive-only. `trend.py`
already supports `--rundir` (`$RUNDIR/<case>/run`) and `--testdir` (a flat
directory), alongside the default archive layout
(`archive/<case>/{atm,ice,lnd}/hist`). Directory generality is the one axis where
trend is ahead; cadence, dimensionality, and variable selection are the fractured
ones.

### `exocam-analysis` (this repo)

Equilibrium-state analysis of small batches. Its reader is `core/reader.py:read_ncfile`,
which reads exactly **one file per call** and returns a raw dict.

- Variables are **hardcoded module-level lists**: `_VARS_2D`, `_VARS_3D`,
  `_VARS_FLUX`, `_VARS_CF_2D`, `_VARS_CF_FLUX`, `_VARS_OPTIONAL_2D`. Every one is
  read from every file. Anything outside those lists reaches the caller only via
  the on-demand `fields_2d` / `fields_3d` cache arguments, which exist to serve
  contour plots, not general analysis.
- There is **no multi-file, no time-axis, and no cadence concept at all** — batch
  mode is a YAML list of separate equilibrium files, each read independently and
  each becoming its own `Diagnostics` object. Time series are simply out of scope.
- Averaging goes through `core/coords.py:area_weighted_avg`, which is a genuinely
  correct implementation but is a *separate* one from trend's.

### `exocam-spectral`

TOA band-resolved spectra. Its reader is `plotspectra_GCM.read_gcm_file`.

- **Single file per call**, `nc.Dataset(fname)`, one `time_idx` (default 0). No
  multi-file path at all; two-run comparison is two independent calls.
- It carries a reassembly step neither other repo has:
  `gather_spectral_bands(ds, prefix, time_idx)` compiles
  `^{prefix}_toa_int(\d+)$`, walks `ds.variables` to discover every match,
  **sorts by the integer band index** (not lexically), and stacks into
  `(nband, nlat, nlon)`. This exists because the GCM writes the spectral dimension
  exploded into one variable per band. Six prefixes are read: `FUL`, `FULC`,
  `FUS`, `FUSC`, `FDS`, `FDSC`.
- Its global mean is a **third** averaging implementation —
  `global_mean_spectrum` takes an arithmetic mean over longitude, then a
  **Gaussian-quadrature-weighted** mean over latitude using the file's own `gw`
  variable, not a cos(lat) or sin-band construction.

Three readers, three weighting implementations, three variable-selection
mechanisms, three different ideas of what "a file" is.

---

## The proposed shape — one primitive

All of this work condenses to a single reading primitive:

> **Read N netCDF history files, given a case, a set of variables, a time
> selection, and a model component. Return the arrays, the coordinates, and the
> area weights.**

Requirements the primitive must meet that no current reader does:

- **Cadence-general.** Not monthly-only. Whatever `nhtfrq` produced, and whatever
  history stream (`h0`, `h1`, …), parsed from the filename rather than assumed.
- **2D and 3D.** `(lat, lon)` and `(lev, lat, lon)` are the same call.
- **Archive or run directory** (trend's existing `--rundir` / `--testdir`
  generality, kept and extended to every component).
- **Explicit variable-set argument.** A Python list passed in by the caller — never
  a file opened by relative path from the working directory.
- **Weights returned with the data** (see below).

Everything else composes on top of it: energy balance, water conservation,
vertical profiles, global and zonal averages, contour plots, spectra. None of
those needs its own file-reading code.

---

## Key design point: weights travel with the data

This is the part that matters most, and it is a correctness argument rather than
a convenience one.

**The correctness layer must be attached to the reader, not left to the caller,
because these failures are silent.** They do not raise. They produce plausible
numbers that are wrong.

Three specific cases:

1. **Area-weighted averaging.** A bare `np.mean` over a CAM lat/lon grid
   overweights the poles — every grid cell counts equally although polar cells
   cover far less area. The result looks like a temperature. It is off by an
   amount you cannot spot by inspection. Both `coords.area_weighted_avg` and
   `trend_core.build_area_weights` / `global_mean_2d` handle this correctly today;
   the risk is any new script that reaches for numpy first.

2. **Gaussian quadrature weights for spectral work.** `exocam-spectral` uses the
   file's own `gw(lat)` variable, not a `cos(lat)` or `sin`-band approximation.
   For a Gaussian grid these are the quadrature weights the model itself
   integrates with; substituting a geometric approximation gives a subtly
   different number with no error and no warning. A unified reader should return
   `gw` when the file has it, and treat it as the preferred latitude weight rather
   than something the caller reconstructs.

3. **`grav` and `mwdry` for non-Earth atmospheres.** This is the worst of the
   three, because **these values are not discoverable from the file.** In this
   repo `core/reader.py` takes them as arguments defaulting to Earth
   (`grav=9.81`, `mwdry=28.966`), and `core/compute.py` turns them into
   `G = raw['grav']` and `R = 8.314462 / (raw['mwdry']/1000.0)`, which then feed
   every mass-weighted and hypsometric diagnostic (`hybrid2height`, column masses,
   all vertical integrals). Run a Titan or a CO2-dominated case and forget the
   flags and you silently get Earth's answer for someone else's atmosphere.
   Nothing in the netCDF corrects you.

If the reader returns weights (and `grav`/`mwdry`, carried as case metadata)
alongside the data, then correct analysis becomes the path of least resistance
and is hard to bypass. A caller would have to work at it to get the wrong answer.

---

## What stays in code vs. what is composed per-question

**The distinguishing test is verifiability, not complexity.**

- Where a wrong answer is **silent** — weighting, spectral band reassembly,
  planetary constants — it belongs in tested library code that is *always* called.
  You cannot catch these by looking at the output.
- Where a wrong answer is **obvious** — which variables, which time window, which
  runs get compared — it can be composed per question. You would see a nonsense
  plot immediately and fix it.

This is why *"a common suite of standard variable averages and plots"* is of only
modest value. The primitives are the reusable part; the composition is the
science, and the composition is genuinely bespoke per campaign. Codifying a fixed
suite of standard plots freezes yesterday's questions into infrastructure while
doing nothing about the failure mode that actually bites.

---

## Migration risk — the reason this is not trivial

**Consolidating three readers is a rewrite that touches every downstream script,
and the failure mode is silent.**

A subtly different weighting convention, or a time window that is off by one
month, or a `-999.0` sentinel handled differently, does not crash. It yields
plausible numbers that quietly disagree with previously published results. The
existing readers already differ in exactly this kind of detail — trend's
`build_area_weights` places cell edges at latitude midpoints with poles at ±90°
and normalizes to sum to 1; `coords.area_weighted_avg` builds a separable
`dlon × d(sin lat)` outer product and masks `-999.0` inline; spectral uses `gw`
directly. These should agree, but "should agree" is exactly the assumption that
produces a silently wrong paper.

**Non-negotiable requirement: reproduce known-good results from each existing
path before retiring it.** Numerical agreement to floating-point precision on a
real case, per repo, per diagnostic. Not a smoke test. This has precedent in this
repo — the April 2026 refactor validated the vectorized `hybrid2height` and
`area_weighted_avg` against the old loop versions exactly this way. The same
discipline, applied across three repos, is most of the cost of the project.

---

## Secondary benefit — run-time analysis and agent-drivability

A cadence-general reader taking an explicit variable-set argument would close most
of the run-time analysis gap as a side effect. Arbitrary variables, arbitrary time
slices, and 3D fields all become the same call with different arguments, instead
of three separate capability requests against three separate readers. The
`vars.in` collision problem disappears entirely, because there is no shared mutable
file to collide over.

It also makes the toolkit far more tractable to drive from an agent: one function
with named arguments is something an agent can call correctly, where "edit this
shared file in the repo root, then run the script from that directory" is
something an agent gets wrong or, worse, gets right while clobbering a concurrent
analysis.

**Scope boundary:** `atm.log` parsing stays outside this. Neither repo does it
today, and it is not netCDF — a log-parsing capability is a separate piece of work
that should not be smuggled into a reader's contract.

---

## Summary

| | Today |
|---|---|
| Readers | 3 independent (`read_monthly_files`, `read_ncfile`, `read_gcm_file`) |
| Weighting implementations | 3 (sin-band normalized, separable outer product, model `gw`) |
| Variable selection | relative-path `vars.in`; hardcoded module lists; hardcoded prefix list |
| Cadence | monthly-only (trend); single-file, no time axis (analysis, spectral) |
| Dimensionality | 2D-only (trend); 2D + 3D (analysis); band × lat × lon (spectral) |

The idea is sound. The work is infrastructure, it is a rewrite with a silent
failure mode, and it was deferred deliberately. Recorded here so it is available
if the priorities change.
