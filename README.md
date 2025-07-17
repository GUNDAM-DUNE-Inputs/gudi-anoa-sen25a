# GUNDAM Inputs and Configuration for DUNE Atmospheric Neutrinos OA

## 1. Clone & Prepare the Repository 📥

```bash
# Clone the main repo
git clone https://github.com/ulyevarou/gudi-anoa-sen25a.git
cd gudi-anoa-sen25a

# Initialize and update submodules
git submodule init
git submodule update
```

## 2. Set Up the GUNDAM Environment ⚙️

To install GUNDAM, follow [GUNDAM Getting Started](https://gundam-organization.github.io/gundam/GettingStarted.html).

Then, in `submit_DUNE.sh`, source your installation. Replace with your actual install path:
This script sets up the environment variables (e.g., `PATH`, `LD_LIBRARY_PATH`) needed for GUNDAM to run correctly.

```bash
source yourpath/gundam/install/setup.sh
```



# NN Home cluster users:
# Recommended GUNDAM environment
```bash
source /home/isgould/work/gundam/gundam-gcc_12-x86_64-linux-gnu/setup.sh
```
# Fallback (if the above doesn't work)
```bash
source /home/isgould/work/gundam/install/setup.sh
```
- `submit_DUNE.sh` (if running manually)
- or `nnhome_slurm.sh` (if using SLURM)

Alternatively, you can add gundam to your path as described in [GUNDAM Getting Started](https://gundam-organization.github.io/gundam/GettingStarted.html):

```bash
export PATH="$INSTALL_DIR/gundam/bin:$PATH"
export LD_LIBRARY_PATH="$INSTALL_DIR/gundam/lib:$LD_LIBRARY_PATH"
```


## 3. Compile gundamOscAnaTools🔧

Run the following command:

```bash
./update-externals.sh
```
Alternatively you can build manually. If needed, configure `./gundamOscAnaTools/resources/TabulateNuOscillator/CMakeLists.txt` and compile by using:

```bash
cd ./gundamOscAnaTools/resources/TabulateNuOscillator/
mkdir build-$(uname -m)
cd build-$(uname -m)
cmake -DCMAKE_INSTALL_PREFIX=${PWD} ..
make install
```
### Troubleshooting: CUDA build errors

If `./update-externals.sh` (or `compile.sh`) fails with an error about `Unsupported gpu architecture 'compute_'` or similar CUDA detection issues, it may be because you’re unintentionally using Conda’s CMake. To force use of the system CMake and CUDA toolchain, deactivate your Conda environment:

```bash
conda deactivate   # run twice if you still see (base) in your prompt
```
Then retry: 

```bash
cd gundamOscAnaTools/resources/TabulateNuOscillator
rm -rf build-*
./compile.sh
# or from repo root:
./update-externals.sh
```

#### 3a. (Optional) Verify your GPU compute capability
Note: you may get a CUDA error, in that case you can check what SM-XX your GPU supports by running:

```bash
nvidia-smi --query-gpu=compute_cap --format=csv,noheader,nounits
```
#### 3b. (Optional) Then you can update compile.sh to match your GPU:

Open: 
gundamOscAnaTools/resources/TabulateNuOscillator/compile.sh
and locate the CMake invocation. Change the architecture flag from "all" to your numeric target, e.g:

-  -DCMAKE_CUDA_ARCHITECTURES="all" \
+  -DCMAKE_CUDA_ARCHITECTURES=52 \



## 4. Set the Input Files Path 📂

In `submit_DUNE.sh`, set the folder containing the input files.


* **On NN home:**

  ```bash
  export OA_INPUT_FOLDER=/storage/shared/DUNE/OA-inputs/atm/gudi-inputs/v3/
  ```

* **On dunegpvm:**

  ```bash
  export OA_INPUT_FOLDER=/pnfs/dune/persistent/users/weishi/OA-inputs/atm_reprocessed_v2/
  ```

## 5. Using NuOscillator Interface (instructions for submodule) 🔧

To add `gundamOscAnaTools` to your path make sure that the following lines appear in `submit_DUNE.sh`:

```bash
export NUOSCILLATOR_ROOT=./gundamOscAnaTools/resources/TabulateNuOscillator/build-x86_64/
source ${NUOSCILLATOR_ROOT}/bin/setup.NuOscillator.sh
```

## 6. **Running GUNDAM  🚀**

Use correctly configured `submit_DUNE.sh` script to perform either a dry run or a full Asimov fit. You can adjust the number of cores according to your cluster using the `-t` option

The `submit_DUNE.sh` script sets up environment variables and runs `gundamFitter` with your desired configuration. You can run it directly or call it from a SLURM job script (e.g., `nnhome_slurm.sh` as you will see below).


### Dry run

```bash
gundamFitter -a -d -c ./config_DUNE.yaml -t 8
```

### Asimov fit (omit `-d`)

```bash
gundamFitter -a -c ./config_DUNE.yaml -t 8
```

Then submit the script:

```bash
./submit_DUNE.sh
```
-------------------
### 📄 The `submit_DUNE.sh` Wrapper Script

The `submit_DUNE.sh` script is a lightweight wrapper that:

- Sets the required environment variables for GUNDAM and NuOscillator
- Chooses the correct input path depending on the cluster
- Calls `gundamFitter` with your configuration and output location

This script allows you to easily manage different run configurations and is also called by `nnhome_slurm.sh` when submitting jobs through SLURM.

Here’s an example of what `submit_DUNE.sh` looks like:

```bash
#!/bin/bash

# Optionally, you can source GUNDAM here if not using a SLURM wrapper (e.g. `nnhome_slurm.sh`):

# at seawulf 
# source /gpfs/home/uyevarouskay/Work/env_gundam_home_main_v1.sh;

# at SBU NN home cluster (you can source GUNDAM if you aren't using nnhome_slurm.sh script):
# source /home/uyevarou/Work/env_dune_lbl.sh

# Set the atmospheric input folder (Current as of June 2025; subject to change): 
export OA_INPUT_FOLDER=/storage/shared/DUNE/OA-inputs/atm/gudi-inputs/v3/

# Set and source the NuOscillator interface (Current as of June 2025; subject to change): 
export NUOSCILLATOR_ROOT=./gundamOscAnaTools/resources/TabulateNuOscillator/build-x86_64/
source ${NUOSCILLATOR_ROOT}/bin/setup.NuOscillator.sh

# Run GUNDAM
gundamFitter -a \
  -c config_DUNE.yaml \
  --out-dir outputs/noOsc-hists-v1 \
  -t 4
```

------------------

### ***sbatch*** 🆚  ***srun***

`sbatch` submits a batch script (e.g., `nnhome_slurm.sh`) to the SLURM scheduler. This is the standard method for running longer GUNDAM jobs or full fits:
```bash
sbatch nnhome_slurm.sh
```
- The job runs in the background.

- Standard output and errors are saved to log files (defined by `#SBATCH --output` and `--error`).

`srun` runs a command **interactively** on an allocated compute node:

- Useful for short test runs, debugging, or launching jobs manually.

- Output appears directly in your terminal.
  
-------------------
## **6.5 Running GUNDAM with SLURM 🧠💻**

If you're working on a SLURM-managed cluster (like NN home), you can submit your GUNDAM jobs using the provided `nnhome_slurm.sh` script.

This script wraps the `submit_DUNE.sh` logic inside a SLURM batch job, and includes GPU usage, memory settings, and email notifications.

-----------------

### ✅ SLURM Submission Script: `nnhome_slurm.sh`

```bash
#!/bin/bash
# Set the name of your SLURM job (this will show up in squeue)
#SBATCH --job-name=ReplaceWithYourName
# Path to write stdout log (feel free to replace path as you wish):
#SBATCH --output=slurm_files/GundamInputDUNE_LBL_%j.out
# Path to write stderr log (feel free to replace path as you wish):
#SBATCH --error=slurm_files/GundamInputDUNE_LBL_%j.err
# Sends an email when the job ends:
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=youremail@yourschool.edu
# Memory required per node:
#SBATCH --mem=50G
# Request 1 GPU; If you want CPU-only, remove this line entirely
#SBATCH --gres=gpu
# Number of CPU cores to use for this job (make sure to adjust this value to match -t in gundamFitter) 
#SBATCH -c 4

# Print system information and timestamp (useful for debugging/logging)
uname -a 
date

# Load GUNDAM environment, if you have your own, source that. Otherwise, you can use this one. 
```bash
source /home/isgould/work/gundam/gundam-gcc_12-x86_64-linux-gnu/setup.sh
```

# Alternate setup (fallback):
```bash
# source /home/isgould/work/gundam/install/setup.sh
```
# Submit the run logic (delegated to submit_DUNE.sh)
```bash 
./submit_DUNE.sh
```
🛠️ uname -a and date log system and timestamp info are helpful for debugging or job tracking.

⚠️ Make sure to replace ReplaceWithYourName and youremail@yourschool.edu with your job name and email.

📚 More on SLURM
For more details on SLURM options and usage, see the [SLURM documentation](https://slurm.schedmd.com/documentation.html).

## 7. Enable/disable/modify Oscillation Parameters ⚙️

**Disable Oscillation Parameters** 

  The Oscillation parameters are provided in the `./inputs/parameters/configParSet.yaml`. To run GUNDAM without oscillation parameters you can use an override file by running the following command

  ```bash
  gundamConfigUnfolder -c config_DUNE.yaml -of overrides/disableOscillationParameters.yaml -o
  ```

  It will create an override configuration file `config_DUNE_With_disableOscillationParameters.json` that can be run by modifying `submit_DUNE.sh` in the following way:

  ```bash
  gundamFitter -a -d -c ./config_DUNE_With_disableOscillationParameters.json -t 8
  ```

  Alternatively you can disable oscillation parameter set in the `./inputs/parameters/configParSet.yaml` by changing the following flag to `false` (not recommended):

  ```yaml
  - name: "Oscillation Parameters"
    isEnabled: false
  ``` 
  And run `submit_DUNE.sh` without modification.

  ```bash
  ./submit_DUNE.sh
  ```

**Change the priors** 


  Different priors can be provided for the oscillation parameters by using corresponding override files. For instance, to set up prior values to PDG recommended use the following override file:

  ```bash
  gundamConfigUnfolder -c config_DUNE.yaml -of overrides/pdg25NormalOscillationParameters.yaml -o
  ```
 
  Then modify `submit_DUNE.sh` to:

  ```bash
  gundamFitter -a -d -c ./config_DUNE_With_pdg25NormalOscillationParameters.json -t 8
  ```

  And run:

  ```bash
  ./submit_DUNE.sh
  ``` 

---

## Using NuOscillator Interface (old instructions for not a submodule case) 🔧

To use the GUNDAM interface for NuOscillator, follow these steps:

**Build the NuOscillator Interface**
Make sure that GUNDAM is set up and in your PATH. You will also need to have the development environment set up. Then navigate to the directory and compile the interface:

```bash
cd gundamOscAnaTools/resources/TabulateNuOscillator/
./compile.sh
```

**On NN home:**

```bash
source /home/uyevarou/Work/env_dune_lbl.sh
cd gundamOscAnaTools/resources/TabulateNuOscillator/
```

**On dunegpvm:**

```bash
source /cvmfs/larsoft.opensciencegrid.org/spack-packages/setup-env.sh
spack load root@6.28.12
spack load cmake@3.27.7%gcc@12.2.0
spack load gcc@12.2.0
# First time compile only to properly get ROOT exports
source /cvmfs/larsoft.opensciencegrid.org/spack-packages/opt/spack/linux-almalinux9-x86_64_v2/gcc-12.2.0/root-6.28.12-sfwfmqorvxttrxgfrfhoq5kplou2pddd/bin/thisroot.sh
```

**Modify the GUNDAM source file and add additional libraries to your PATH** (included in `env_dune_lbl.sh`):

```bash
source /YOURPATH/gundamOscAnaTools/resources/TabulateNuOscillator/build-x86_64/bin/setup-NuOscillator.sh
# or alternatively:
export LD_LIBRARY_PATH="/YOURPATH/gundamOscAnaTools/resources/TabulateNuOscillator/build-x86_64/_deps/oscprob-src/lib/:$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH="/YOURPATH/gundamOscAnaTools/resources/TabulateNuOscillator/build-x86_64/lib/:$LD_LIBRARY_PATH"
```

---

**Modify the Run Line** 📝

In the `submit_DUNE.sh` script, update the run line as follows:

```bash
export NUOSCILLATOR_ROOT=/YOURPATH/gundamOscAnaTools/resources/TabulateNuOscillator/
export NUOSCILLATOR_ROOT_LIB=/YOURPATH/gundamOscAnaTools/resources/TabulateNuOscillator/build-x86_64/
```

---

## Run GUNDAM 🎉

```bash
./submit_DUNE.sh
```

Happy fitting! 😊
