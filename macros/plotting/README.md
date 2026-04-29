# `getERPlots.C`

`getERPlots.C` reads a GUNDAM output ROOT file from the pre-fit `model` and `data` directories and produces plots and event-rate summaries for `numu` and `nue` selections, split into CC and NC categories.

The macro produces:

- stacked reco/true energy plots by interaction mode;
- stacked reco/true `cos(theta_Z)` plots by interaction mode;
- 2D reco and true energy-vs-`cos(theta_Z)` histograms;
- inclusive and per-mode event-rate summaries;
- separate reco-energy binning for `numu` and `nue` selections.

---

## Input file

The input ROOT file is set inside `getERPlots()`:

```cpp
const char* inputFileName =
    "/gpfs/projects/McGrewGroup/uyevarou/Work/atm/gudi-anoa-sen25a/outputs/"
    "gundamFitter_config_DUNE_Asimov_DryRun.root";
```

To run on a different GUNDAM output file, edit this path.

The macro expects the following ROOT-file structure:

```text
FitterEngine/preFit/model
FitterEngine/preFit/data
```

Each relevant selection directory should contain an `events_TTree`.

---

## Running the macro

Run:

```bash
root -l -b -q getERPlots.C
```

To save the terminal event-rate printout:

```bash
root -l -b -q getERPlots.C > ER_count.txt
```

Use `>` to create a new output file, or `>>` to append to an existing file.

---

## Output

The macro creates one PDF per category:

```text
modeStacks_numuSel_CC.pdf
modeStacks_nueSel_CC.pdf
modeStacks_numuSel_NC.pdf
modeStacks_nueSel_NC.pdf
```

Each PDF contains:

1. reco-energy stack by interaction mode;
2. reco-`cos(theta_Z)` stack by interaction mode;
3. true-energy stack by interaction mode;
4. true-`cos(theta_Z)` stack by interaction mode;
5. 2D model reco histogram;
6. 2D model true histogram;
7. 2D data reco histogram.

The terminal output contains inclusive and per-mode event-rate summaries.

---

## Main options

### Custom reco-energy binning

Custom reco-energy binning is controlled by:

```cpp
static const bool kUseCustomRecoBinning = true;
```

If this option is `true`, the macro uses separate reco-energy binning for `numu` and `nue` selections.

If it is `false`, the macro uses the common uniform binning:

```cpp
static const int    kRecoNBinsUniform = 50;
static const double kRecoMinUniform   = 0.0;
static const double kRecoMaxUniform   = 10.0;
```

### `numu` reco-energy binning

Used for:

```text
numuSel_CC
numuSel_NC
```

Edit:

```cpp
const std::vector<double> kRecoBinEdgesNuMu = {
    0.0, 0.25, 0.5, 0.75, 1.0,
    1.25, 1.5, 1.75, 2.0,
    2.25, 2.5, 2.75, 3.0,
    3.5, 4.0, 5.0, 6.0, 8.0, 10.0
};
```

### `nue` reco-energy binning

Used for:

```text
nueSel_CC
nueSel_NC
```

Edit:

```cpp
const std::vector<double> kRecoBinEdgesNuE = {
    0.0, 0.5, 1.0, 1.5, 2.0,
    2.5, 3.0, 3.5, 4.0,
    5.0, 6.0, 8.0, 10.0
};
```

### Truth-energy binning

Truth-energy binning is common to all selections:

```cpp
static const int    kTrueNBins = 50;
static const double kTrueMin   = 0.0;
static const double kTrueMax   = 10.0;
```

### `cos(theta_Z)` binning

The same `cos(theta_Z)` binning is used for reco and truth histograms:

```cpp
static const int    kCosNBins = 100;
static const double kCosMin   = -1.0;
static const double kCosMax   =  1.0;
```

---

## Fast test mode

To run script for the small number of events:

```cpp
for (Long64_t i = 0; i < nEntries; ++i) {
```

For quick debugging, temporarily replace it with:

```cpp
for (Long64_t i = 0; i < std::min<Long64_t>(nEntries, 100); ++i) {
```