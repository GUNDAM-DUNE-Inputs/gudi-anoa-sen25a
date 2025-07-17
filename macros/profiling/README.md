# Oscillation‑Parameter Profiling with GUNDAM

Build profile scans of oscillation parameters using the GUNDAM framework.

---

## Prerequisites

## 1. Checkout & Compile

Switch to the dedicated GUNDAM profiling branch and initialize submodules:

```bash
git checkout feature/profiling
git pull
git submodule sync
git submodule update --init --recursive
```

Recompile GUNDAM as usual (e.g. `make` ).

---

## 2. Configure Submission Scripts

### 2.1 `submit_DUNE_profiling.sh`

1. Source your GUNDAM environment:
   ```bash
   source /YOURPATH_GUNDAM/setup.sh
   ```
2. Edit the input‑folder variable to point at your root files:
   ```bash
   OA_INPUT_FOLDER="/path/to/your/inputs"
   ```

### 2.2 `submit_slurm_profiling_DUNE.sh`

1. Set the override YAML for your fit configuration (default: MaCh3Comparison300Zx300Esmoothing04Profiling.yaml):
   ```bash
   input_yaml="override/yourProfilingConfig.yaml"
   ```
2. Adjust SLURM directives (e.g. `--time`, `--mem`, `--array`) to match your cluster’s policies.

---

## 3. Prepare the Parameter Injector

Customize the injector template with your prior values (default: NuFit priors):

```bash
./macros/profiling/injector/parameterInjector.yaml
```

---

## 4. Submit the SLURM Array

Determine the total number of jobs:

```
jobs = N_params × (2 × NUM_SIGMA + 1)
```

- For 6 parameters and `NUM_SIGMA=2` (5 points each), you need `6 × 5 = 30` jobs.

Create a folder for output files and submit the array:

```bash
mkdir ../../outputs/profiling/

sbatch submit_slurm_profiling_DUNE.sh
```

---

## 5. Plot the Profiles

After all profiling jobs finish, use a `plotOApar.cpp` to plot oscillation parameters profiles. You need to specify the directory containing your generated `.root` scan files as `rootDir`

```bash
root -l -b plotOApar.cpp
```

This will produce χ²‑vs‑parameter plots for each oscillation parameter in the current working directory.

---