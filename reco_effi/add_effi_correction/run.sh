#!/bin/bash
set -e
beam_array=("pPb" "Pbp")
beam_pol=("Up")
#nbr_chunks=60 # for "MC-Sim09k-Multiplicity"
nbr_chunks=30
#MC="MC-Sim09d" #for MC pp
MC="MC-Sim09e"  # for MC of pPb Pbp
#MC="MC-Sim09k-Multiplicity"
for b in "${beam_array[@]}"
do
    for p in "${beam_pol[@]}"
    do
    python /home3/sara.sellam/RpPb_identified_hadrons_project/efficiencies/reco_effi/add_effi_correction/add_effi_correction.py -beam "$b" -pol "$p" -MC "$MC"  -nbr_chunks "$nbr_chunks"
    done
done
    