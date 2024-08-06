import sys

import numpy as np
import pandas as pd
import sys
from ROOT import *
from ROOT import RDF
from ROOT import RDataFrame
import glob
sys.path.append("/home3/sara.sellam/RpPb_identified_hadrons_project")

from selections.selections import *
from Binning.Binning import *
from selections.selections import _branch_list,_ana_selection,_track_type, _P,_PT, _ETA_cms, _selections,_Kinematic_range, _bad_runs, _ETA
from Ntuple_path.call_ntuple import _path, _corrections,get_data
from Ntuple_path.call_ntuple import  _corrections
ROOT.EnableImplicitMT()
TH1.SetDefaultSumw2(True)    
TH2.SetDefaultSumw2(True)  

git_path="/home3/sara.sellam/RpPb_identified_hadrons_project"
git_out_path="/scratch43/ssellam/results"

def make_ratio_histogram(beam="",particle="",track_type="",weights=None,binning_mode="",ratio="",test_mode=False,merged_mc=True,merged_pol=True):
    _pt_bins=binning[binning_mode][beam]["pt"]
    _eta_bins=binning[binning_mode][beam]["eta"] 
    print("requested binning",binning_mode)
    print("_pt_bins",_pt_bins)
    print("_eta_bins",_eta_bins)
    #frac computations
    if beam =="pp":
        merged_mc=False
        pol="Down"
        data1="MC-Sim09d"
        data2="Gen"+data1
    if weights!="0":
        frac_all=TH2D(ratio+"_"+track_type+"_"+particle,"",len(_pt_bins)-1,_pt_bins,len(_eta_bins)-1,_eta_bins) 
        h_num_all=TH2D("h_num_all","",len(_pt_bins)-1,_pt_bins,len(_eta_bins)-1,_eta_bins) 
        h_denom_all=TH2D("h_denom_all","",len(_pt_bins)-1,_pt_bins,len(_eta_bins)-1,_eta_bins) 
        if merged_mc:
            for data1 in ["MC-Sim09k-Multiplicity","MC-Sim09e"]:
                data2="Gen"+data1
                if merged_pol:
                    for pol in ["Up","Down"]:
                        f=TFile(git_out_path+"/efficiencies/reco_effi/rfiles/"+ratio+"_"+track_type+"_"+particle+"_"+beam+"_"+pol+"_w_"+weights+"_"+data1+"_"+data2+"_"+binning_mode+"_num_denom.root","open")
                        h_num=f.Get("h_num")
                        h_denom=f.Get("h_denom")
                        h_num_all.Add(h_num)
                        h_denom_all.Add(h_denom)
            
            frac_all.Divide(h_num,h_denom,1.,1.,"B")  
            data1="MC-Merged"
            data2="Gen"+data1
            f_out=TFile(git_out_path+"/efficiencies/reco_effi/rfiles/"+ratio+"_"+track_type+"_"+particle+"_"+beam+"_Merged_w_"+weights+"_"+data1+"_"+data2+"_"+binning_mode+".root","recreate")
            h_num.SetDirectory(0)
            h_denom.SetDirectory(0)
            f_out.cd()
            print("Writing")
            h_num.Write()
            h_denom.Write()
            frac_all.Write()
            f_out.Close()
        else:
            f=TFile(git_out_path+"/efficiencies/reco_effi/rfiles/"+ratio+"_"+track_type+"_"+particle+"_"+beam+"_"+pol+"_w_"+weights+"_"+data1+"_"+data2+"_"+binning_mode+"_num_denom.root","open")
            frac_all=TH2D(ratio+"_"+track_type+"_"+particle,"",len(_pt_bins)-1,_pt_bins,len(_eta_bins)-1,_eta_bins) 
            h_num_all=f.Get("h_num")
            h_denom_all=f.Get("h_denom")
            frac_all.Divide(h_num_all,h_denom_all,1.,1.,"B")  
            f_out=TFile(git_out_path+"/efficiencies/reco_effi/rfiles/"+ratio+"_"+track_type+"_"+particle+"_"+beam+"_"+pol+"_w_"+weights+"_"+data1+"_"+data2+"_"+binning_mode+".root","recreate")
            h_num_all.Write()
            h_denom_all.Write()
            frac_all.Write()
            f_out.Close()


       


            
import argparse

# Create an ArgumentParser object to handle command line arguments
parser = argparse.ArgumentParser()

# Add the command line arguments
parser.add_argument("-beam", "--beam", nargs='+', type=str, help="Beam type")
#parser.add_argument("-pol", "--pol", nargs='+', type=str, help="Polarization")
parser.add_argument("-particle", "--particle", nargs='+', type=str, help="particle")
parser.add_argument("-binning_mode", "--binning_mode", nargs='+', type=str, help="binning_mode")
parser.add_argument("-weights", "--weights", nargs='+', type=str, help="Occupancy weight")
#parser.add_argument("-MC", "--MC", nargs='+', type=str, help="Which MC sample")

parser.add_argument("-test_mode", "--test_mode", action="store_true", help="Test Mode")

# Parse the command line arguments
args = parser.parse_args()

# Extract the values of the arguments
beam = args.beam[0]
#pol=args.pol[0]
probe=args.particle[0]
binning_mode=args.binning_mode[0]
weights= args.weights[0]
#MC= args.MC[0]

test_mode = args.test_mode


# Print the values of the arguments
print(f"Beam type: {beam}")
#print(f"Polarization : {pol}")
print(f"Particle : {probe}")
print(f"binning_mode : {binning_mode}")
print(f"Test Mode : {test_mode}")




make_ratio_histogram(beam=beam,particle=probe,track_type="prompt",weights=weights,binning_mode=binning_mode,ratio="reco_effi",test_mode=test_mode)



