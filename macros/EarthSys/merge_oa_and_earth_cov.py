#!/usr/bin/env python3
import ROOT
import numpy as np
import matplotlib.pyplot as plt
import builtins

# ============================================================
# USER SETTINGS – EDIT THESE
# ============================================================

# --- OA covariance file and object names ---
OA_FILE          = "oscCovAveragedNufit52.root"   # your OA file
OA_COV_NAME      = "osc_param_cov"                # TMatrixD name
OA_NAMES_NAME    = "osc_param_names"              # TObjArray of TObjString
OA_CENTRAL_NAME  = "osc_param_priors"             # TVectorD with priors (central or sigmas)

# --- Earth covariance file and object names (from EarthSystematics.save_covariance_root) ---
EARTH_FILE         = "EarthCov.root"
EARTH_COV_NAME     = "EarthCovMatrix"
EARTH_NAMES_NAME   = "ParNames"
EARTH_CENTRAL_NAME = "ParCentral"

# --- Output ROOT + PNG ---
OUT_ROOT_FILE     = "Merged_OA_EarthCov.root"
OUT_COV_NAME      = "MergedCovMatrix"
OUT_NAMES_NAME    = "MergedParamNames"
OUT_CENTRAL_NAME  = "MergedParCentral"
OUT_PNG           = "MergedCovMatrix.png"


# ============================================================
# Helpers
# ============================================================

def read_cov_names_priors(filename, cov_name, names_name, central_name):
    """
    Read TMatrixD (cov), TObjArray(TObjString) (names),
    and TVectorD (central values) from ROOT.

    Returns: (cov_np, names_list, central_np)
    """
    f = ROOT.TFile.Open(filename, "READ")
    if not f or f.IsZombie():
        raise RuntimeError(f"Cannot open file: {filename}")

    # Covariance
    cov = f.Get(cov_name)
    if not cov:
        raise RuntimeError(f"Cannot find TMatrixD '{cov_name}' in {filename}")

    n = cov.GetNrows()
    M = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            M[i, j] = cov[i][j]

    # Names
    names_obj = f.Get(names_name)
    if not names_obj:
        raise RuntimeError(f"Cannot find TObjArray '{names_name}' in {filename}")

    names = []
    for i in range(names_obj.GetLast() + 1):
        obj = names_obj.At(i)
        if not obj:
            continue
        names.append(obj.GetName())   # TObjString::GetName() -> content

    # Priors: central values
    v_central = f.Get(central_name)
    if not v_central:
        raise RuntimeError(f"Cannot find TVectorD '{central_name}' in {filename}")

    if v_central.GetNrows() != n:
        raise RuntimeError(
            "Size mismatch between covariance and priors in "
            f"{filename}: cov={n}, central={v_central.GetNrows()}"
        )

    central = np.array([v_central[i] for i in range(n)], dtype=float)

    f.Close()
    return M, names, central


def write_merged_to_root(filename, cov_np, names_list,
                         central_np,
                         cov_name=OUT_COV_NAME,
                         names_name=OUT_NAMES_NAME,
                         central_name=OUT_CENTRAL_NAME):
    """
    Save merged covariance, names, and priors to a ROOT file.
    """
    npar = cov_np.shape[0]
    if len(names_list) != npar:
        raise RuntimeError("names_list length does not match covariance dimension")
    if central_np.shape[0] != npar:
        raise RuntimeError("central vector length does not match covariance dimension")

    f = ROOT.TFile(filename, "RECREATE")

    # 1) Covariance matrix as TMatrixD
    M = ROOT.TMatrixD(npar, npar)
    for i in range(npar):
        for j in range(npar):
            M[i][j] = float(cov_np[i, j])
    M.Write(cov_name)

    # 2) Names as TObjArray of TObjString

    name_array = ROOT.TObjArray(npar)
    #name_array.SetName("ParNames")  # This works in PyROOT when called before AddAt
    for i,name in enumerate(names_list):
        name_array.AddAt(ROOT.TObjString(name), i)
    name_array.Write(names_name, ROOT.TObject.kSingleKey)


    # 3) Central values as TVectorD
    v_central = ROOT.TVectorD(npar)
    for i in range(npar):
        v_central[i] = float(central_np[i])
    v_central.Write(central_name)

    f.Close()
    print(f"Saved merged covariance, names and priors to {filename}")


def plot_covariance_png(cov_np, names_list, png_name):
    """
    Make an annotated PNG of the covariance matrix.
    """
    npar = cov_np.shape[0]
    if npar == 0:
        print("Empty covariance, nothing to plot.")
        return

    vabs = float(np.max(np.abs(cov_np)))
    fig_size = max(6, 0.4 * npar + 4)
    fig, ax = plt.subplots(figsize=(fig_size, fig_size))

    im = ax.imshow(cov_np, vmin=-vabs, vmax=vabs, cmap="PuOr")
    fig.colorbar(im, ax=ax)

    ax.set_xticks(np.arange(npar))
    ax.set_yticks(np.arange(npar))
    ax.set_xticklabels(names_list, rotation=90)
    ax.set_yticklabels(names_list)

    # Dynamically choose font size so text fits cells reasonably
    base_font = builtins.max(6, int(160 / max(1, npar)))

    for i in range(npar):
        for j in range(npar):
            val = cov_np[i, j]
            ax.text(j, i, f"{val:.2e}",
                    ha="center", va="center",
                    fontsize=base_font,
                    color="black" if abs(val) < vabs * 0.6 else "white")

    ax.set_title("Merged Covariance Matrix (OA + Earth)")
    plt.tight_layout()

    if not png_name.lower().endswith(".png"):
        png_name += ".png"
    plt.savefig(png_name, dpi=200)
    print(f"Saved merged covariance plot as {png_name}")
    plt.close(fig)


# ============================================================
# Main
# ============================================================

def main():
    # 1) Read OA inputs
    print(f"Reading OA covariance and priors from {OA_FILE} ...")
    cov_oa, names_oa, central_oa = read_cov_names_priors(
        OA_FILE, OA_COV_NAME, OA_NAMES_NAME, OA_CENTRAL_NAME
    )
    print(f"  OA: {cov_oa.shape[0]} parameters")

    # 2) Read Earth inputs
    print(f"Reading Earth covariance and priors from {EARTH_FILE} ...")
    cov_earth, names_earth, central_earth = read_cov_names_priors(
        EARTH_FILE, EARTH_COV_NAME, EARTH_NAMES_NAME, EARTH_CENTRAL_NAME
    )
    print(f"  Earth: {cov_earth.shape[0]} parameters")

    n_oa = cov_oa.shape[0]
    n_e  = cov_earth.shape[0]
    n_total = n_oa + n_e

    # 3) Build block-diagonal merged covariance
    merged_cov = np.zeros((n_total, n_total), dtype=float)
    merged_cov[:n_oa, :n_oa] = cov_oa
    merged_cov[n_oa:, n_oa:] = cov_earth

    # 4) Concatenate names and priors
    merged_names   = names_oa + names_earth
    merged_central = np.concatenate((central_oa, central_earth))

    # 5) Save to new ROOT file
    write_merged_to_root(
        OUT_ROOT_FILE,
        merged_cov,
        merged_names,
        merged_central,
        cov_name=OUT_COV_NAME,
        names_name=OUT_NAMES_NAME,
        central_name=OUT_CENTRAL_NAME
    )

    # 6) Draw PNG
    plot_covariance_png(merged_cov, merged_names, OUT_PNG)


if __name__ == "__main__":
    main()
