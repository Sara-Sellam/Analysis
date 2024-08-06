#!/bin/bash
set -e


beam_array=("pPb" "Pbp" "pp")
#MC_array=("MC-Sim09k-Multiplicity" "MC-Sim09k-Multiplicity" "MC-Sim09d")
#beam_array=("pPb" "Pbp")

#MC_array=("MC-Sim09d")
MC_array=("MC-Merged" "MC-Merged" "MC-Sim09d")

#pol_array=("Down")
pol_array=("Merged" "Merged" "Down")

weight_set_array=("1" "2" "7")

probe_array=("pi" "K" "p")
#probe_array=("noPID")

binning_mode="oscar_ana"
#binning_mode=pi_0_ana
#binning_mode="sara_ana"



python /home3/sara.sellam/RpPb_identified_hadrons_project/efficiencies/reco_effi/plot_effi.py  -beam "${beam_array[@]}" -pol "${pol_array[@]}" -particle "${probe_array[@]}" -binning_mode "$binning_mode"  -weights "${weight_set_array[@]}" -MC "${MC_array[@]}" # >> "$output_file" 

