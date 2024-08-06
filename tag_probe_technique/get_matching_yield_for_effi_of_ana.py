
import sys
import glob
from math import *
from numpy import tan, arcsin 
import json
import datetime
import os
from array import array
import multiprocessing
import uproot as up
import statistics

from ROOT import *
from ROOT import RDF
from ROOT import RDataFrame
TH1.SetDefaultSumw2(True)    
TH2.SetDefaultSumw2(True)
gStyle.SetPaintTextFormat("1.3f")
gStyle.SetPalette(1)
sys.path.append("/home3/sara.sellam/RpPb_identified_hadrons_project")

from Binning.Binning_calibration import binning_TagandProbe
from Binning.Binning import *
from efficiencies.pid_effi.pid_selections import _pid_selections
from Ntuple_path.call_calibration_ntuple  import _part,Tuple_Names_for_TagandProbe
from selections.selections import _ETA


from Fit_Model import *
from fit_configuration import fit_opt

from Helpers_function.sys_functions import weighted_std, gen_Toys

from ple_yield.utils import *
ROOT.EnableImplicitMT()


git_in_path="/scratch43/ssellam/results"
git_out_path="/scratch43/ssellam/results"



import argparse

# Create an ArgumentParser object to handle command line arguments
parser = argparse.ArgumentParser()

# Add the command line arguments
parser.add_argument("-b", "--beam", nargs='+', type=str, help="Beam type")
parser.add_argument("-data_type", "--data_type", nargs='+', type=str, help="data_type")
parser.add_argument("-mom", "--mother", nargs='+', type=str, help="The mother of the probe mom")
parser.add_argument("-probe", "--probe", nargs='+', type=str, help="The probe mom")
parser.add_argument("-pid", "--pid", nargs='+', type=str, help="requested pid selection")
parser.add_argument("-binning_mode", "--binning_mode", nargs='+', type=str, help="binning_mode")

# Parse the command line arguments
args = parser.parse_args()

# Extract the values of the arguments
beam = args.beam[0]
mom=args.mother[0]
probe=args.probe[0]
pid=args.pid[0]
binning_mode=args.binning_mode[0]
data_type=args.data_type[0]
effi=probe+"_as_"+pid

tuple_name=_part[beam][mom]['Tuple']
binning_mode_dict={"ana":{"h1_h1":"tight",
                                "h1_h":"croaser"},
                     "Fit_sys":{"h1_h1":"tight",
                                "h1_h":"croaser"},
                        "Bin_sys_v1":{"h1_h1":"sys_tight_v1",
                                      "h1_h":"sys_croaser_v1"},
                        "Bin_sys_v2":{"h1_h1":"sys_tight_v2",
                                      "h1_h":"sys_croaser_v2"},
                        "Bin_sys_v3":{"h1_h1":"sys_tight_v3",
                                      "h1_h":"sys_croaser_v3"},
                        "loose_pid_cuts_sys":{"h1_h1":"tight",
                                            "h1_h":"croaser"},
                        "tight_pid_cuts_sys":{"h1_h1":"tight",
                                            "h1_h":"croaser"}}

if probe == pid:
    binning_name=binning_mode_dict["ana"]["h1_h1"]
    _pt_bins=array('d',binning_TagandProbe[tuple_name][binning_name]["PT"])
    _eta_bins=array('d',binning_TagandProbe[tuple_name][binning_name]["ETA"])
else:
    binning_name=binning_mode_dict["ana"]["h1_h"]
    _pt_bins=array('d',binning_TagandProbe[tuple_name][binning_name]["PT"])
    _eta_bins=array('d',binning_TagandProbe[tuple_name][binning_name]["ETA"])


if "MC" in data_type:
    df_data,cuts=process_TrueID_selected_data(data=data_type, beam=beam, pol="Down", particle=probe, weights="1",test_mode=False)
    cuts = _Kinematic_range["MC"][probe] + "&&" + _ana_selection[data][beam][pol]["TrueID_"+probe]+"&&pi_MC_isPrompt==1&&"+_PID_selections[pid][beam]
else:
    df_data,cuts=process_data(data="Data", beam=beam, pol="Down", particle=probe, test_mode=False)


out_file=TFile(git_out_path+"/efficiencies/pid_effi/tag_probe_technique/rfiles/yield_for_effi_"+effi+"_"+beam+"_"+binning_name+"_"+data_type+"_ana.root","recreate")

cuts=cuts+"&&pi_isMuon==0"

df_data=df_data.Filter(cuts)
if "MC" in data_type:
    df_data=df_data.Define("ETA_lab",_ETA["MC"][probe])
else:
    df_data=df_data.Define("ETA_lab",_ETA["Data"][probe])

print("cuts",cuts)
print("_eta_bins",_eta_bins)
h_data=df_data.Histo2D(RDF.TH2DModel("N_tracks_for_effi","h_data",len(_pt_bins)-1,_pt_bins,len(_eta_bins)-1,_eta_bins),"pi_PT","ETA_lab")
h_data.SetDirectory(0)
out_file.cd()
h_data.Write()
for i in range(len(_eta_bins)-1):
    h_data_1d=h_data.ProjectionX("N_tracks_for_effi_{}_{}".format(_eta_bins[i],_eta_bins[i+1]),i,i+1)
    h_data_1d.SetDirectory(0)
    out_file.cd()
    h_data_1d.Write()


out_file.Close()    
