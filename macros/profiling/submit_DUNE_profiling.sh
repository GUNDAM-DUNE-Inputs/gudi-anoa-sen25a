#!/bin/bash

## Set Up the GUNDAM Environment
# at seawulf
source /gpfs/home/uyevarouskay/Work/env_gundam_home_v4.sh; #latest version

# at SBU NN home cluster
#source /home/isgould/work/gundam/install/setup.sh

## Set the Input Files Path
# at seawulf cluster
export OA_INPUT_FOLDER=/gpfs/scratch/uyevarouskay/atm/gundam_files_fin_v13/  #add is bad

# at SBU NN home cluster
#export OA_INPUT_FOLDER=/storage/shared/DUNE/OA-inputs/atm/gudi-inputs/v3/

# dunegpvm
#export OA_INPUT_FOLDER=/pnfs/dune/persistent/users/weishi/OA-inputs/atm_reprocessed_v2/ #might be outdated

## NuOscillator 
export NUOSCILLATOR_ROOT_LIB=../../gundamOscAnaTools/resources/TabulateNuOscillator/build-x86_64/
source ${NUOSCILLATOR_ROOT_LIB}/bin/setup.NuOscillator.sh

#!/usr/bin/env bash
set -euo pipefail
set -x

occurrence="${1:?Need which isFixed occurrence to flip (1..6))}"
newPrior="${2:?Need new priorValue (e.g. 0.403)}"
paramName="${3:?Need parameter basename (e.g. sin23_0)}"
input_yaml="${4:?Need a parameter config (e.g. ./overrides/MaCh3Comparison300Zx300Esmoothing04Profiling.yaml)}"

#paramName="PMNS_SIN_SQUARED_23"
#occurrence=3 #(3rd in the list of the parameters)
#newPrior=0.403
#input_yaml="./overrides/MaCh3Comparison300Zx300Esmoothing04Profiling.yaml"
#input_yaml="./overrides/MaCh3ComparisonFixedH300Zx300Eprofiling.yaml"
input_yaml_toys="./overrides/toysForProfiling.yaml"

safePrior="${newPrior//./_}"

base="${input_yaml##*/}"; base="${base%.*}"
override_dir="./macros/profiling/overrides/"
configs_json_dir="./macros/profiling/configs_json/"

# ensure dirs
cd ../..
mkdir -p "$override_dir" "$configs_json_dir"

override_output_yaml="$override_dir/${paramName}_${safePrior}.yaml"
json_out="$configs_json_dir/config_DUNE_With_${paramName}_${safePrior}.json"

#json_out="configs_json/config_DUNE_With_${base}_fixed.json"
#override_output_yaml="${override_dir}/${base}_fixed.yaml"


# only flip the *first* isFixed line to true
#sed '0,/^[[:space:]]*isFixed:/ s/^\([[:space:]]*\)isFixed:.*$/\1isFixed: true/' \
#    "$input_yaml" \
#  > "$override_output_yaml"

# -- flip the Nth isFixed to true ------------------------------------
awk -v N="$occurrence" -v V="$newPrior" '
  BEGIN       { count=0; changePV=0 }
  # Flip the Nth isFixed line:
  /^[[:space:]]*isFixed:/ && ++count==N {
    match($0, /^[[:space:]]*/)                # grab indent
    printf "%sisFixed: true\n", substr($0, RSTART, RLENGTH)
    changePV=1                                 # flag to change next priorValue
    next                                       # skip default print
  }
  # If flagged, rewrite the very next priorValue:
  /^[[:space:]]*priorValue:/ && changePV {
    match($0, /^[[:space:]]*/)                # grab indent
    printf "%spriorValue: %s\n", substr($0, RSTART, RLENGTH), V
    changePV=0                                 # done with this block
    next
  }
  { print }                                    # everything else unchanged
' "$input_yaml" > "$override_output_yaml"

#produce json
gundamConfigUnfolder \
  -c config_DUNE.yaml \
  -of "$override_output_yaml" \
  -of "$input_yaml_toys" \
  -o "$json_out"

gundamFitter -a -c config_DUNE.yaml -of "$override_output_yaml" -of "$input_yaml_toys" --inject-toy-parameter ./macros/profiling/injector/parameterInjector.yaml --out-dir ./outputs/profiling/ -t 8 --toy 0