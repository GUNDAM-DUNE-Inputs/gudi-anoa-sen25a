# Earth Model Systematics Covariance

The code in this directory is copied from  
`MaCh3_DUNE/Utils/atmospherics/EarthSystematicsCov/` and is used to construct  
Earth model parameter **correlation** and **covariance** matrices.

It has been modified to:

- Produce a ROOT file containing a `TMatrixD` covariance matrix.
- Save the corresponding parameter names and priors as vectors  
  (via the `save_covariance_root` function).

Additionally, the script `merge_oa_and_earth_cov.py` builds a **joint covariance matrix** that also includes oscillation parameter priors from **NuFit 5.2**.  
This combined matrix can be used directly in **GUNDAM** configuration files.

---

## 1. Create the Earth Model correlation/covariance matrix

To create the Earth Model correlation matrix with the default settings and save
the covariance matrix into `EarthCov.root`:

```bash
module load root
pip install --user pandas

# Check that Python can import ROOT and pandas
python -c "import ROOT, pandas; print('OK')"

python example.py
```
Alternatively, you can submit example.py to SLURM:

```bash
sbatch runMtxCov.sh
```
## 2. Merge with oscillation parameter covariance matrix

To merge the Earth Model covariance matrix with the oscillation-parameter
covariance matrix (NuFit 5.2 priors), run:

```bash
python merge_oa_and_earth_cov.py
```

