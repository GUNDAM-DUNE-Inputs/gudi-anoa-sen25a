
## Prepare Inputs for GUNDAM

This repository contains a script that prepares files in a GUNDAM-friendly format by copying true and reco variables from the initial `cafTree` to a new `event_tree`. Please note, that the names of the branches are modified. The copied variables and their new names are listed in `fileConvertor_ATM.cpp` (lines 261–302). Additionally, new added branches are listed between lines 252 and 260.

### New Branches Added
The following branches are added in the `event_tree`: 

- `iniNuFlavor`   // nue: 12, numu: 14, nuebar: -12, numubar: -14
- `interNuFlavor` // nue: 12, numu: 14, nutau: 16, nuebar: -12, numubar: -14, nutaubar: -16 
- `detNuFlavor`   //nue selection: 12, numu selection: 16, NC selection: 0
- `sample_idx`    //sample id from 1-36

## How to Build and Run the Script

### Compiling the Script

To compile the `fileConvertor_ATM.cpp` script, use the following command:

```bash
g++ -o fileConvertor_ATM fileConvertor_ATM.cpp -I$(root-config --incdir) $(root-config --libs) -std=c++17
```

This command compiles the C++ code and links it against ROOT libraries. Make sure ROOT is properly installed and configured on your system.

### Running the Script Locally


Once the script is compiled, use the `run_jobs.sh` script to execute it.

Make the script executable and run it:

```bash
chmod +x run_jobs.sh
./run_jobs.sh
```

### Submitting Jobs to SLURM 

1. **Edit the Slurm script**  

   - Update the `#SBATCH` directives (partition, time, memory, CPUs) to match your cluster’s available resources.  
   - Set the `OA_INPUTS` and `OA_OUTPUTS` variables to point at your data directories.

2. **Submit the job**  

```bash
sbatch run_jobs_slurm.sh
```
