import sys
import glob
from math import *
from numpy import tan, arcsin 
import json
import datetime
import os
from array import array
import multiprocessing

from ROOT import *
from ROOT import RDF
from ROOT import RDataFrame
TH1.SetDefaultSumw2(True)    
TH2.SetDefaultSumw2(True)
gStyle.SetPaintTextFormat("1.3f")
gStyle.SetPalette(1)
sys.path.append("/home3/sara.sellam/RpPb_identified_hadrons_project")


from Binning.Binning_calibration import binning_TagandProbe
from efficiencies.pid_effi.pid_selections import _pid_selections
from Ntuple_path.call_phi_for_TM  import _path 
from efficiencies.TM_effi.fit import mass_modelling
import argparse

# Create an ArgumentParser object to handle command line arguments
parser = argparse.ArgumentParser()

# Add the command line arguments
parser.add_argument("-beam", "--beam", nargs='+', type=str, help="Beam type")
parser.add_argument("-pol", "--pol", nargs='+', type=str, help="Polarity")

# Parse the command line arguments
args = parser.parse_args()

# Extract the values of the arguments
beam_type = args.beam[0]
pol=args.pol[0]

_pt_bins=
_eta_bins=

in_file_path=_path[beam][pol]
data_list=glob.glob(in_file_path)
in_tree=TChain()

if isinstance(data_list,list):
    for data in [data_list[0]]:
        for tree in in_tree_list:
            print("tree",tree)
            in_tree.Add(data + "/" + tree + "/DecayTree")
            if weights!="0":
                #occupancy_weights.Add(ion_occupancy_weights_before_sWeights[mom][beam]+"/tree")
                print("weights")
else:           
    data=data_list
    for tree in in_tree_list:
        print("tree",tree)
        in_tree.Add(data + "/" + tree + "/DecayTree")
        if weights!="0":
            #occupancy_weights.Add(ion_occupancy_weights_before_sWeights[mom][beam]+"/tree")
            print("weights")
assert in_tree.GetEntries()!=0

df_MM=RDataFrame(in_tree)

cuts=""
        
df_MM=df_MM.Filter(cuts)

if __name__ == '__main__':
    # Define the number of processes to use (adjust as needed)
    num_processes = 4
        
    outFile_name_after=git_out_path+"/efficiencies/TM_effi/rfiles/canvas_"+beam+"_"+effi
    
    pool = multiprocessing.Pool(processes=num_processes)
    for i in range(len(_pt_bins)-1):
        for j in range(len(_eta_bins)-1):
            pt_bin_size=in_range("Kmi_PT",str(_pt_bins[i]),str(_pt_bins[i+1]))
            probe_eta="(-log(tan(asin(Kmi_PT/Kmi_PT)/2)))"
            eta_bin_size=in_range(probe_eta,str(_eta_bins[j]),str(_eta_bins[j+1]))
            bin_size=pt_bin_size+"&&"+eta_bin_size
            df_slice=df_MM.Filter(bin_size)
            df_slice=df_slice.AsNumpy([_part[beam][mom]["mom"]+"_MM"])
            args_list.append((i, j, df_slice, fit_opt, weights, binned, beam,outFile_name))
            # Use the pool to parallelize the calculations
    results_list[mode] = pool.map(calculate_mass_modelling, args_list)
    # Close the pool of processes
    pool.close()
    pool.join()
           
       