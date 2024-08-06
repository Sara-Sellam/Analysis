
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
from uncertainties import ufloat
from uncertainties.umath import * 
from sklearn.linear_model import LinearRegression
import numpy as np


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

from Fit_Model import *
from fit_configuration import fit_opt

from Helpers_function.sys_functions import weighted_std, gen_Toys

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
parser.add_argument("-weights", "--weights", nargs='+', type=str, help="requested weights")

# Parse the command line arguments
args = parser.parse_args()

# Extract the values of the arguments
beam = args.beam[0]
mom=args.mother[0]
probe=args.probe[0]
pid=args.pid[0]
binning_mode=args.binning_mode[0]
weights=args.weights[0]
data_type=args.data_type[0]

effi=probe+"_as_"+pid

def find_closest(arr, target):
    return arr[min(range(len(arr)), key=lambda i: abs(arr[i] - target))]


tuple_name=_part[beam][mom]['Tuple']
binning_mode_dict={"ana":{"h1_h1":"tight_old",
                                "h1_h":"croaser"},
                    "ana2":{"h1_h1":"oscar_ana_tight",
                            "h1_h":"oscar_ana_croaser"},
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



_eta_bins_ana=binning[binning_mode][beam]["eta"]
_pt_bins_ana=binning[binning_mode][beam]["pt"]


h_pid={}


out_file=TFile(git_out_path+"/efficiencies/pid_effi/tag_probe_technique/rfiles/PT_eta/effi_"+effi+"_"+beam+"_"+binning_mode+"_ana_weight_"+weights+"_v1.root","recreate")

for eta in range(len(_eta_bins_ana)-1):
    print("eta",_eta_bins_ana[eta],_eta_bins_ana[eta+1])
    pid_file=TFile(git_in_path+"/efficiencies/pid_effi/tag_probe_technique/rfiles/PT_eta/ratio_"+effi+"_"+beam+"_"+binning_name+"_ana_weight_"+weights+".root","open")
    NTracks_file=TFile(git_in_path+"/efficiencies/pid_effi/tag_probe_technique/rfiles/yield_for_effi_"+effi+"_"+beam+"_"+binning_name+"_"+data_type+"_ana.root","open")


    ratio_pid=pid_file.Get("ratio_ana")
    print("ratio_pid",ratio_pid.GetEntries())
    ratio_nTracks=NTracks_file.Get("N_tracks_for_effi")

    h_pid[_eta_bins_ana[eta],_eta_bins_ana[eta+1]]=TH1D("ratio_ana_"+str(_eta_bins_ana[eta])+"_"+str(_eta_bins_ana[eta+1]),"",len(_pt_bins_ana)-1,_pt_bins_ana)
    
    h_sys=TH1D("h_std_"+str(_eta_bins_ana[eta])+"_"+str(_eta_bins_ana[eta+1]),"",len(_pt_bins_ana)-1,_pt_bins_ana)

    eta_center=(_eta_bins_ana[eta]+_eta_bins_ana[eta+1])/2
    if beam=="pp":
        eta_center=eta_center
    elif beam=="pPb":
        eta_center=eta_center+0.5
    else:
        eta_center=(eta_center+0.5)*-1
    
    effi_array=[]
    effi_errors=[]
    mean_efficiencies = []
    mean_efficiencies_err= []
    
    effi_value_array=[]
    effi_error_array=[]
    nTracks_value_array=[]
    nTracks_error_array=[]

    bin_center=[]
    for pt_bin in range(len(_pt_bins)-1):
        pt_center=(_pt_bins[pt_bin]+_pt_bins[pt_bin+1])/2
        bin=ratio_pid.FindBin(pt_center,eta_center)
        print("pt_center",pt_center,"eta_center",eta_center)
        effi_value=ratio_pid.GetBinContent(bin)
        print("effi_value",effi_value)
        effi_value_err=ratio_pid.GetBinError(bin)
        nTracks_value=ratio_nTracks.GetBinContent(bin)
        nTracks_value_err=ratio_nTracks.GetBinError(bin)

        effi_value_array.append(effi_value)
        effi_error_array.append(effi_value_err)

        nTracks_value_array.append(nTracks_value)
        nTracks_error_array.append(nTracks_value_err)

        bin_center.append(pt_center)
        
    new_effi_value_array=effi_value_array
    new_effi_error_array=effi_error_array


    print("see if it has zero new_effi_value_array",new_effi_value_array)
    if np.any(np.asarray(effi_value_array) != 0) :   
        X = np.array([bin_center[i] for i in range(len(bin_center)) if effi_value_array[i] != 0]).reshape(-1, 1)
        y = np.array([effi_value_array[i] for i in range(len(effi_value_array)) if effi_value_array[i] != 0])
        y_err = np.array([effi_error_array[i] for i in range(len(effi_value_array)) if effi_value_array[i] != 0])

        model = LinearRegression()
        model2 = LinearRegression()
        print("bin_center",bin_center)
        print("effi_value_array",effi_value_array)
        print("X",X)
        print("y",y)
        model.fit(X, y) 
        model2.fit(X, y_err)
        zero_effi_indices = [i for i, e in enumerate(effi_value_array) if e == 0]
        for zero_index in zero_effi_indices:
            predicted_effi = model.predict(np.array([bin_center[zero_index]]).reshape(-1, 1))[0]
            predicted_effi_error = model2.predict(np.array([bin_center[zero_index]]).reshape(-1, 1))[0]

            new_effi_value_array[zero_index] = predicted_effi
            new_effi_error_array[zero_index] = predicted_effi_error
    
    print("effi_value_array",effi_value_array)
    print("new_effi_value_array",new_effi_value_array)



    for pt_bin in range(len(_pt_bins_ana)-1):
        
        lower_bound = _pt_bins_ana[pt_bin]
        upper_bound = _pt_bins_ana[pt_bin+1]


        #indices_within_range = [index for index, value in enumerate(_pt_bins) if lower_bound <= value <= upper_bound]
        pt_within_ana_range = [element for element in _pt_bins if lower_bound <= element <= upper_bound]
        closest_lower = find_closest(_pt_bins, lower_bound)
        closest_upper = find_closest(_pt_bins, upper_bound)
        
        # Find indices of closest values in _pt_bins
        closest_lower_index = _pt_bins.index(closest_lower)
        closest_upper_index = _pt_bins.index(closest_upper)
        pt_within_ana_range = _pt_bins[closest_lower_index:closest_upper_index + 1]


        print("_pt_bins caibration",_pt_bins)
        print("lower_bound",lower_bound,"upper_bound",upper_bound,"pt_within_ana_range",pt_within_ana_range)
        
        effi_array=[]
        effi_errors=[]
        


        if len(pt_within_ana_range)>2:
            #print("*"*10," more than 2 in range")
            mean_eff_err_2=0
            final_effi=ufloat(0,0)
            tot_tracks=ufloat(0,0)

            effi_value_array_within_range=[]
            effi_error_array_within_range=[]
            for indx in range(len(pt_within_ana_range)-1):
                pt_center_within_range=(pt_within_ana_range[indx]+pt_within_ana_range[indx+1])/2
                idx_pt_center_within_range=bin_center.index(pt_center_within_range)

                effi_value=new_effi_value_array[idx_pt_center_within_range]
                effi_value_err=new_effi_error_array[idx_pt_center_within_range]
                nTracks_value=nTracks_value_array[idx_pt_center_within_range]
                nTracks_value_err=nTracks_error_array[idx_pt_center_within_range]
                
                print("final_effi +=",final_effi,"effi_value ",effi_value,"nTracks_value",nTracks_value)
                
                effi_value_array_within_range.append(effi_value)
                effi_error_array_within_range.append(abs(effi_value_err))

                final_effi+=ufloat(effi_value,abs(effi_value_err))*ufloat(nTracks_value,nTracks_value_err)
                tot_tracks+=ufloat(nTracks_value,nTracks_value_err)
            if tot_tracks>0:
                mean_eff = final_effi/tot_tracks
            else:
                mean_eff=ufloat(0,0)
     
            std=weighted_std(effi_value_array_within_range,effi_error_array_within_range)
        
        elif len(pt_within_ana_range)==2:
            std=0
            for indx in range(len(pt_within_ana_range)-1):
                pt_center_within_range=(pt_within_ana_range[indx]+pt_within_ana_range[indx+1])/2
                idx_pt_center_within_range=bin_center.index(pt_center_within_range)
                value=new_effi_value_array[idx_pt_center_within_range]
                value_err=new_effi_error_array[idx_pt_center_within_range]
            
                mean_eff=value
                mean_eff_err=value_err+0.01
                mean_eff=ufloat(mean_eff,abs(mean_eff_err))
                #print("mean_eff",mean_eff)
        else:
            std=0
            print()
            if pt_within_ana_range[0]!= _pt_bins[0]:
                upper_limit=_pt_bins.index(pt_within_ana_range[0])
                lower_limit=upper_limit-1
                pt_center_within_range=(_pt_bins[upper_limit]+_pt_bins[lower_limit])/2
                idx_pt_center_within_range=bin_center.index(pt_center_within_range)
                value=new_effi_value_array[idx_pt_center_within_range]
                value_err=new_effi_error_array[idx_pt_center_within_range]
            else:
                lower_limit=_pt_bins.index(pt_within_ana_range[0])
                upper_limit=lower_limit+1
                pt_center_within_range=(_pt_bins[upper_limit]+_pt_bins[lower_limit])/2
                idx_pt_center_within_range=bin_center.index(pt_center_within_range)
                value=new_effi_value_array[idx_pt_center_within_range]
                value_err=new_effi_error_array[idx_pt_center_within_range]
            mean_eff=value
            mean_eff_err=value_err+0.01
            mean_eff=ufloat(mean_eff,abs(mean_eff_err))
            """
            std=0
            #bin=ratio_pid.FindBin(pt_within_ana_range[0],eta_center)
            bin=ratio_pid.FindBin(lower_bound,eta_center)
            mean_eff=ratio_pid.GetBinContent(bin)
            mean_eff_err=ratio_pid.GetBinError(bin)

            #print("*"*10,"1 in range")
            mean_eff=ufloat(mean_eff,mean_eff_err)
            """
        
        print("mean_eff ",mean_eff)
        mean_efficiencies.append(mean_eff.nominal_value)
        mean_efficiencies_err.append(mean_eff.std_dev)
        
        h_pid[_eta_bins_ana[eta],_eta_bins_ana[eta+1]].SetBinContent(pt_bin+1,mean_eff.nominal_value )
        tot_error=sqrt(mean_eff.std_dev**2+std**2)
        h_pid[_eta_bins_ana[eta],_eta_bins_ana[eta+1]].SetBinError(pt_bin+1,tot_error)

        
        #print("eta",_eta_bins_ana[eta],"pt_bin",_pt_bins[pt_bin],"final ========> effi_value",mean_eff.nominal_value)

        h_sys.SetBinContent(pt_bin+1,sqrt(mean_eff.std_dev**2+std**2))

    h_pid[_eta_bins_ana[eta],_eta_bins_ana[eta+1]].SetDirectory(0)
    h_sys.SetDirectory(0)
    out_file.cd()
    h_sys.Write()
    h_pid[_eta_bins_ana[eta],_eta_bins_ana[eta+1]].Write()
    assert len(mean_efficiencies)==len(_pt_bins_ana)-1
out_file.Close()


out_file=TFile(git_out_path+"/efficiencies/pid_effi/tag_probe_technique/rfiles/PT_eta/effi_"+effi+"_"+beam+"_"+binning_mode+"_ana_weight_"+weights+"_v2.root","recreate")

for eta in range(len(_eta_bins_ana)-1):
    print("eta",_eta_bins_ana[eta],_eta_bins_ana[eta+1])
    pid_file=TFile(git_in_path+"/efficiencies/pid_effi/tag_probe_technique/rfiles/PT_eta/ratio_"+effi+"_"+beam+"_"+binning_name+"_ana_weight_"+weights+".root","open")
    NTracks_file=TFile(git_in_path+"/efficiencies/pid_effi/tag_probe_technique/rfiles/yield_for_effi_"+effi+"_"+beam+"_"+binning_name+"_"+data_type+"_ana.root","open")


    ratio_pid=pid_file.Get("ratio_ana")
    print("ratio_pid",ratio_pid.GetEntries())
    ratio_nTracks=NTracks_file.Get("N_tracks_for_effi")

    h_pid[_eta_bins_ana[eta],_eta_bins_ana[eta+1]]=TH1D("ratio_ana_"+str(_eta_bins_ana[eta])+"_"+str(_eta_bins_ana[eta+1]),"",len(_pt_bins_ana)-1,_pt_bins_ana)
    
    h_sys=TH1D("h_std_"+str(_eta_bins_ana[eta])+"_"+str(_eta_bins_ana[eta+1]),"",len(_pt_bins_ana)-1,_pt_bins_ana)

    eta_center=(_eta_bins_ana[eta]+_eta_bins_ana[eta+1])/2
    if beam=="pp":
        eta_center=eta_center
    elif beam=="pPb":
        eta_center=eta_center+0.5
    else:
        eta_center=(eta_center+0.5)*-1
    
    effi_array=[]
    effi_errors=[]
    mean_efficiencies = []
    mean_efficiencies_err= []
    
    effi_value_array=[]
    effi_error_array=[]
    nTracks_value_array=[]
    nTracks_error_array=[]

    bin_center=[]
    for pt_bin in range(len(_pt_bins)-1):
        pt_center=(_pt_bins[pt_bin]+_pt_bins[pt_bin+1])/2
        bin=ratio_pid.FindBin(pt_center,eta_center)
        print("pt_center",pt_center,"eta_center",eta_center)
        effi_value=ratio_pid.GetBinContent(bin)
        print("effi_value",effi_value)
        effi_value_err=ratio_pid.GetBinError(bin)
        nTracks_value=ratio_nTracks.GetBinContent(bin)
        nTracks_value_err=ratio_nTracks.GetBinError(bin)

        effi_value_array.append(effi_value)
        effi_error_array.append(effi_value_err)

        nTracks_value_array.append(nTracks_value)
        nTracks_error_array.append(nTracks_value_err)

        bin_center.append(pt_center)
        
    new_effi_value_array=effi_value_array
    new_effi_error_array=effi_error_array


    for pt_bin in range(len(_pt_bins_ana)-1):
        
        lower_bound = _pt_bins_ana[pt_bin]
        upper_bound = _pt_bins_ana[pt_bin+1]


        #indices_within_range = [index for index, value in enumerate(_pt_bins) if lower_bound <= value <= upper_bound]
        pt_within_ana_range = [element for element in _pt_bins if lower_bound <= element <= upper_bound]
        closest_lower = find_closest(_pt_bins, lower_bound)
        closest_upper = find_closest(_pt_bins, upper_bound)
        
        # Find indices of closest values in _pt_bins
        closest_lower_index = _pt_bins.index(closest_lower)
        closest_upper_index = _pt_bins.index(closest_upper)
        pt_within_ana_range = _pt_bins[closest_lower_index:closest_upper_index + 1]


        print("_pt_bins caibration",_pt_bins)
        print("lower_bound",lower_bound,"upper_bound",upper_bound,"pt_within_ana_range",pt_within_ana_range)
        
        effi_array=[]
        effi_errors=[]
        


        if len(pt_within_ana_range)>2:
            #print("*"*10," more than 2 in range")
            mean_eff_err_2=0
            final_effi=ufloat(0,0)
            tot_tracks=ufloat(0,0)

            effi_value_array_within_range=[]
            effi_error_array_within_range=[]
            for indx in range(len(pt_within_ana_range)-1):
                pt_center_within_range=(pt_within_ana_range[indx]+pt_within_ana_range[indx+1])/2
                idx_pt_center_within_range=bin_center.index(pt_center_within_range)

                effi_value=new_effi_value_array[idx_pt_center_within_range]
                effi_value_err=new_effi_error_array[idx_pt_center_within_range]
                nTracks_value=nTracks_value_array[idx_pt_center_within_range]
                nTracks_value_err=nTracks_error_array[idx_pt_center_within_range]
                
                print("final_effi +=",final_effi,"effi_value ",effi_value,"nTracks_value",nTracks_value)
                
                effi_value_array_within_range.append(effi_value)
                effi_error_array_within_range.append(abs(effi_value_err))

                final_effi+=ufloat(effi_value,abs(effi_value_err))*ufloat(nTracks_value,nTracks_value_err)
                tot_tracks+=ufloat(nTracks_value,nTracks_value_err)
            if tot_tracks>0:
                mean_eff = final_effi/tot_tracks
            else:
                mean_eff=ufloat(0,0)
     
            std=weighted_std(effi_value_array_within_range,effi_error_array_within_range)
        
        elif len(pt_within_ana_range)==2:
            std=0
            for indx in range(len(pt_within_ana_range)-1):
                pt_center_within_range=(pt_within_ana_range[indx]+pt_within_ana_range[indx+1])/2
                idx_pt_center_within_range=bin_center.index(pt_center_within_range)
                value=new_effi_value_array[idx_pt_center_within_range]
                value_err=new_effi_error_array[idx_pt_center_within_range]
            
                mean_eff=value
                mean_eff_err=value_err+0.01
                mean_eff=ufloat(mean_eff,abs(mean_eff_err))
                #print("mean_eff",mean_eff)
        else:
            std=0
            print()
            if pt_within_ana_range[0]!= _pt_bins[0]:
                upper_limit=_pt_bins.index(pt_within_ana_range[0])
                lower_limit=upper_limit-1
                pt_center_within_range=(_pt_bins[upper_limit]+_pt_bins[lower_limit])/2
                idx_pt_center_within_range=bin_center.index(pt_center_within_range)
                value=new_effi_value_array[idx_pt_center_within_range]
                value_err=new_effi_error_array[idx_pt_center_within_range]
            else:
                lower_limit=_pt_bins.index(pt_within_ana_range[0])
                upper_limit=lower_limit+1
                pt_center_within_range=(_pt_bins[upper_limit]+_pt_bins[lower_limit])/2
                idx_pt_center_within_range=bin_center.index(pt_center_within_range)
                value=new_effi_value_array[idx_pt_center_within_range]
                value_err=new_effi_error_array[idx_pt_center_within_range]
            mean_eff=value
            mean_eff_err=value_err+0.01
            mean_eff=ufloat(mean_eff,abs(mean_eff_err))
            """
            std=0
            #bin=ratio_pid.FindBin(pt_within_ana_range[0],eta_center)
            bin=ratio_pid.FindBin(lower_bound,eta_center)
            mean_eff=ratio_pid.GetBinContent(bin)
            mean_eff_err=ratio_pid.GetBinError(bin)

            #print("*"*10,"1 in range")
            mean_eff=ufloat(mean_eff,mean_eff_err)
            """
        
        print("mean_eff ",mean_eff)
        mean_efficiencies.append(mean_eff.nominal_value)
        mean_efficiencies_err.append(mean_eff.std_dev)
        
        h_pid[_eta_bins_ana[eta],_eta_bins_ana[eta+1]].SetBinContent(pt_bin+1,mean_eff.nominal_value )
        tot_error=sqrt(mean_eff.std_dev**2+std**2)
        h_pid[_eta_bins_ana[eta],_eta_bins_ana[eta+1]].SetBinError(pt_bin+1,tot_error)

        
        #print("eta",_eta_bins_ana[eta],"pt_bin",_pt_bins[pt_bin],"final ========> effi_value",mean_eff.nominal_value)

        h_sys.SetBinContent(pt_bin+1,sqrt(mean_eff.std_dev**2+std**2))

    h_pid[_eta_bins_ana[eta],_eta_bins_ana[eta+1]].SetDirectory(0)
    h_sys.SetDirectory(0)
    out_file.cd()
    h_sys.Write()
    h_pid[_eta_bins_ana[eta],_eta_bins_ana[eta+1]].Write()
    assert len(mean_efficiencies)==len(_pt_bins_ana)-1
out_file.Close()











    
    










    
    