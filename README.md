# Ingress

A web tool for watching any patch of sky change over time in images from
NASA's SPHEREx mission. Built for the 2026 NASA Space Apps SPHEREx challenge.

Enter coordinates and a date range; Ingress loads every SPHEREx exposure of
that field and plays them as an animation, so that objects which move or
change brightness stand out.

[Read the Ingress methodology (PDF)](Ingress-Methodology.pdf)

## Run it

Open **https://ingress-iuxy.onrender.com** in a browser. Nothing to install.

## Using it

- **Coordinates**: decimal degrees (`133.7951 -7.2451`) or sexagesimal
  (`08h55m10.831s -7d14m42.53s`, `08:55:10.8 -07:14:42`). RA must be 0–360°,
  Dec −90° to +90°.
- **Field of view to load**: 1–20 arcmin. After loading you can switch to any
  smaller field instantly.
- **Dates**: from 2025-04-01 to today.
- **Wavelength**: all light, one of the six detector bands, or a custom range
  between 0.75 and 5.0 µm.
- **Two views**, chosen above the image:
  - **Compare over time** (the default): one frame per epoch, each built from
    the same wavelengths, so a change between frames is a change over time.
  - **Explore wavelengths**: one frame per exposure, day or survey pass, each
    at its own wavelength, for seeing the sky across the infrared.
- Either view can show images or the difference from the average.
- Space plays/pauses; the arrow keys step through frames.
- **Measure an object**: click anything in the image to get its spectrum
  (one colour per survey pass) and its brightness over time, both with error
  bars, plus a CSV download.

All exposures of the field are downloaded once. Changing the wavelength
range, field of view, frame grouping or display settings afterwards works on
data already in the browser.

### Linking to a field

```
/?ra=269.4464&dec=4.7686&fov=5&start=2025-04-01&end=2026-11-15&band=D3&zoom=2&group=day
```

`ra` and `dec` (degrees) are required; the rest are optional. `obj` selects
an object to measure.

## How it works

1. **Search.** IRSA's Simple Image Access service lists the Level 2 spectral
   images whose footprint contains the position.
2. **Read.** For each image, only the pixel rows covering the field are read
   from the public `nasa-irsa-spherex` bucket with HTTP range requests (a few
   percent of each 66 MB file). IRSA's cutout service is the fallback.
3. **Resample.** Each exposure is resampled onto one north-up grid centred
   on the target (3.075″ pixels up to 12′ fields, coarser above), flagged
   pixels are dropped and the median background is removed.
4. **Stream.** Exposures are streamed to the browser as they finish, with
   their time, wavelength and flags. The browser stacks them into frames.

A SPHEREx exposure sees each point of sky at one wavelength, which depends
on where the point falls on the detector. Each survey pass covers a point
about 100 times, once per spectral channel, so a frame for one band is the
average of the exposures of that pass whose wavelength falls in the band.

Fields with more than 900 exposures (near the ecliptic poles) are sampled
evenly in time; narrow the date range to see every exposure.

## Comparing over time without mixing in wavelength

A SPHEREx exposure sees each point of sky at one wavelength, and the
wavelengths a field receives differ from one survey pass to the next.
Blinking exposures, or averages of "everything in band 3", therefore mixes
a change of wavelength with a change over time: a steady red star appears
to vary. The two views keep these apart (`frontend/compare.js`).

**Compare over time** builds one frame per epoch from matched wavelengths:

1. Each detector band is cut into six wavelength slices, 36 in all.
2. Only slices observed in every epoch are used. An epoch is dropped (and
   listed as left out) only if keeping it would leave under a quarter of a
   typical epoch's slices.
3. Within an epoch each slice is averaged, then the slices are averaged with
   equal weight, so every epoch has the same spectral weighting. Exposures
   more than three times noisier than is typical for their detector are left
   out, as are isolated outlier pixels.
4. If fewer than two epochs share any wavelengths, the page says that no
   comparison is available and shows a single labelled epoch instead.

A panel under the image states what is being compared: the slices used, and
for each epoch its dates, exposures used, mean wavelength, offset from the
other epochs and relative star width (both measured from field stars), plus
the fraction of the field seen in every epoch. It warns when epochs are
offset by more than 1.5″, when star widths differ by more than 10%, or when
under 80% of the field is shared.

**Explore wavelengths** carries a standing warning that brightness changes
between frames can be colour alone, labels every frame with its wavelength,
and says that its difference view does not isolate variability.

Limits that remain, listed in the page: matching is by slice, not exact
wavelength; alignment and star shape are measured and reported but not
corrected; each exposure's median sky level is removed; fast movers smear
within an epoch.

## Measuring brightness

Clicking an object runs aperture photometry on every exposure already in the
browser (`frontend/measure.js`):

- **Spectrum.** Flux in a 9.5″ aperture against the wavelength of each
  exposure at that position; sky level and error bars from a surrounding ring.
- **Brightness over time.** A SPHEREx exposure has one wavelength, so raw
  flux against time would trace the spectrum, not variability. Each
  measurement is therefore divided by the median of all measurements at
  nearly the same wavelength, and the ratios are averaged per pass and band.
- **Moving objects** are followed: the aperture is re-centred in each pass.
- **Saturated stars** cannot be measured at their centre. For these the
  median brightness of a ring just outside the saturated core is used, which
  tracks changes over time but not the star's total flux.

The summary line reports the largest change between two passes and calls it
variable at 5σ or more.

A measured object can be shared: clicking one adds `&obj=<ra>,<dec>` to the
address, and opening that link selects and measures it again.

**Validation** ([details](validation/README.md)): fluxes of ordinary stars
come out 4–13% below 2MASS and AllWISE, as expected without an aperture
correction, and the period fitted to R Draconis is 245 days against a
catalogued 245.6.

## Test targets

The five targets in the menu (`TARGETS` in `backend/presets.py`) are kept
ready so they open in seconds:

- Their exposure lists are bundled in `backend/presets.json`, so opening one
  never waits on the archive's search service, and works when it is down.
- With `INGRESS_PREWARM=1` (set in the Dockerfile) the server prepares their
  images at 10′ in the background when it starts. This takes a few minutes;
  until it finishes they load at normal speed.
- Every 12 hours the server asks the archive whether they have new exposures.

After a data release, refresh the bundled lists with
`python -m backend.presets` and commit the result.

## Period search and the opening demonstration

Opened with no field in the address, Ingress loads R Draconis, selects it
and runs the period search (`frontend/period.js`) on the result. The panel
beside the image shows the recovered period with its 1σ range next to the
catalogued 245.6 days, their absolute and percentage difference, and whether
the two are consistent. Nothing is fixed in advance: the figures are
computed from the exposures loaded, and change if the dates do.

The same search runs for any object clicked. Stars with a catalogued period
(the `REFERENCES` list in `frontend/period.js`) are compared with it; others
get a period only if the data support one: enough nights spread over at
least six months and 1.5 cycles, outlying nights set aside, 30% of the
variation explained,
and the same period recovered from the even and the odd nights separately.
Most fields have too few nights, and the panel says so.

## Updating for new data releases

The image collections searched are set in one place, `COLLECTIONS` in
`backend/spherex.py`, or with the `SPHEREX_COLLECTIONS` environment variable
(comma separated, most preferred first). New exposures added to the existing
collections appear automatically.

Other settings: `INGRESS_CACHE` (cache directory), `INGRESS_PARALLEL`
(simultaneous bucket reads, default 24), `INGRESS_SEARCH_TTL` (seconds a
search result is reused, default 6 hours), `SPHEREX_BUCKET_URL`.

## Running your own copy

```bash
pip install -r requirements.txt
```

```bash
uvicorn backend.app:app --port 8000
```

Then open http://localhost:8000. Or with Docker:

```bash
docker build -t ingress . && docker run -p 8000:8000 ingress
```

The server needs outbound HTTPS access to `irsa.ipac.caltech.edu` (search)
and `nasa-irsa-spherex.s3.amazonaws.com` (image data).

## Layout

- `backend/spherex.py`: archive search, download, resampling, limits
- `backend/app.py`: API (`/api/stream`, `/api/config`) and static files
- `frontend/`: the page, with no build step and no dependencies

Data: NASA/IPAC Infrared Science Archive, SPHEREx Spectral Image MEFs.
