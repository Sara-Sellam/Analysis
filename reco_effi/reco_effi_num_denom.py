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
git_out_path="/scratch43/ssellam/results/"

def make_ratio_histogram(beam="",pol="",data1="",particle="",track_type="",weights=None,data2="",binning_mode="",ratio="",test_mode=False):
    h_num,h_denom,frac={},{},{}
    ######## selections ###########
    num_sel=_Kinematic_range["MC"][particle]+"&&"+_track_type[track_type][particle]
    denom_sel=_Kinematic_range["GenMC"]["gen_"+particle]+"&&"+_track_type[track_type]["gen_"+particle]
    print("numerator selections: ", num_sel)
    print("denominator selections: ", denom_sel)

    ################################


    mc_input_files_list=get_data(_path[data1][beam][pol][0],test_mode)
    ch_num=TChain()
    pseudoIP_ch=TChain()
    track_reco_corr_ch=TChain()
    ch_num_weights=TChain()
    for idx, mc_input_files in enumerate(mc_input_files_list):
        for in_file in mc_input_files:
            ch_num.Add(in_file+"/"+_path[data1][beam][pol][1][particle])
            in_file_name=in_file.split('/')[-1]
            pseudoIP_ch.Add(_corrections["pseudoIP"][data1][beam][pol][0][idx]+"pseudoip_"+in_file_name+"/"+_corrections["pseudoIP"][data1][beam][pol][1])
            track_reco_corr_ch.Add(_corrections["track_reco_corr_sys"][data1][beam][pol][0][idx]+"track_reco_corr_"+in_file_name+"/DecayTree")
            
            if weights!="0":
                ch_num_weights.Add(_corrections["weights"+weights][data1][beam][pol][0][idx]+"w_"+weights+"_"+in_file_name+"/"+_corrections["weights"+weights][data1][beam][pol][1])

    assert ch_num.GetEntries()== pseudoIP_ch.GetEntries()
    assert pseudoIP_ch.GetEntries()== track_reco_corr_ch.GetEntries()
    ch_num.AddFriend(track_reco_corr_ch)
    if weights!="0": 
        if beam=="pp":
            ch_num.AddFriend(ch_num_weights)
            ch_num.AddFriend(track_reco_corr_ch)
            _df_num=RDataFrame(ch_num)
            df_num=_df_num.Filter(num_sel)
        else:
            ch_num.AddFriend(ch_num_weights)
            _df_num=RDataFrame(ch_num)
            df_num=_df_num.Filter(num_sel)
        assert track_reco_corr_ch.GetEntries()==ch_num_weights.GetEntries()
    else:
        _df_num=RDataFrame(ch_num)
        df_num=_df_num.Filter(num_sel)
   

    
    

    if beam=="pp":
        df_num=df_num.Filter("pi_TRACK_GhostProb<0.3")
    #assert  df_num.Count().GetValue() !=0
    P=_P["MC"][particle]
    
    df_num=df_num.Define("pi_ETA_CMS",_ETA_cms[beam]["MC"][particle])

    if weights!="0":
        df_num=df_num.Define("weights1","w").Define("weights2","Trcalib_ratio").Define("tot_weight","weights1*weights2")
    else:
        df_num=df_num.Define("tot_weight","Trcalib_ratio")




    #denominator

    ch_denom=TChain()
    ch_denom_weights=TChain()
    P=_P["GenMC"][particle]
    PT=_PT["GenMC"][particle]
    gen_mc_input_files_list=get_data(_path[data2][beam][pol][0],test_mode)
    for idx, gen_mc_input_files in enumerate(gen_mc_input_files_list):
        for in_file in gen_mc_input_files:
            ch_denom.Add(in_file+"/"+_path[data2][beam][pol][1][particle])
            in_file_name=in_file.split('/')[-1]
            if weights!="0":
                ch_denom_weights.Add(_corrections["weights"+weights][data2][beam][pol][0][idx]+"true_w_"+weights+"_"+in_file_name+"/"+_corrections["weights"+weights][data2][beam][pol][1][particle])
    
    if weights!="0":   
        if beam =="pp":
            ch_denom_weights.AddFriend(ch_denom)
            df_denom=RDataFrame(ch_denom_weights)
            df_denom=df_denom.Define(particle+"_TRUEP",P).Filter(denom_sel)
            df_denom=df_denom.Define(particle+"_TRUE_ETA_CMS",_ETA_cms[beam]["GenMC"][particle]).Define("tot_weight","w").Redefine(particle+"_TRUEPT",PT)
        elif data2=="GenMC-Sim09e":
            ch_denom_weights.AddFriend(ch_denom)
            ch_denom=ch_denom_weights
            
            df_denom=RDataFrame(ch_denom)
            df_denom=df_denom.Define(particle+"_TRUEP",P).Filter(denom_sel)
            df_denom=df_denom.Define(particle+"_TRUE_ETA_CMS",_ETA_cms[beam]["GenMC"][particle]).Define("tot_weight","w")
     
        else:
            ch_denom.AddFriend(ch_denom_weights)
            df_denom=RDataFrame(ch_denom)
            df_denom=df_denom.Define(particle+"_TRUEP",P).Filter(denom_sel)
            df_denom=df_denom.Define(particle+"_TRUE_ETA_CMS",_ETA_cms[beam]["GenMC"][particle]).Define("tot_weight","w")
            
        assert ch_denom.GetEntries()==ch_denom_weights.GetEntries()
    else:
        df_denom=RDataFrame(ch_denom)
        df_denom=df_denom.Define(particle+"_TRUEP",P).Filter(denom_sel)
        df_denom=df_denom.Define(particle+"_TRUE_ETA_CMS",_ETA_cms[beam]["GenMC"][particle])

    _pt_bins=binning[binning_mode][beam]["pt"]
    _eta_bins=binning[binning_mode][beam]["eta"] 
    print("requested binning",binning_mode)
    print("_pt_bins",_pt_bins)
    print("_eta_bins",_eta_bins)
    #frac computations
    if weights!="0": 
        f=TFile(git_out_path+"/efficiencies/reco_effi/rfiles/"+ratio+"_"+track_type+"_"+particle+"_"+beam+"_"+pol+"_w_"+weights+"_"+data1+"_"+data2+"_"+binning_mode+"_num_denom.root","recreate")
        frac_all=TH2D(ratio+"_"+track_type+"_"+particle,"",len(_pt_bins)-1,_pt_bins,len(_eta_bins)-1,_eta_bins)
        h_num=df_num.Histo2D(RDF.TH2DModel("h_num","",len(_pt_bins)-1,_pt_bins,len(_eta_bins)-1,_eta_bins),"pi_PT","pi_ETA_CMS","tot_weight")
        h_denom=df_denom.Histo2D(RDF.TH2DModel("h_denom","",len(_pt_bins)-1,_pt_bins,len(_eta_bins)-1,_eta_bins),particle+"_TRUEPT",particle+"_TRUE_ETA_CMS","tot_weight")  
        h_num.SetDirectory(0)
        h_denom.SetDirectory(0)
        f.cd()
        print("Writing")
        h_num.Write()
        h_denom.Write()
        f.Close()
        """
        print("now dividing ")          
        frac_all.Divide(h_num.GetPtr(),h_denom.GetPtr(),1.,1.,"B")  
        print("I am here ")
        frac={}
        f.cd()
        frac_all.Write()
        """
        """
        for i in range(len(_eta_bins)-1):
            frac[i,i+1]=frac_all.ProjectionX(ratio+"_"+track_type+"_"+particle+"_{}_{}".format(_eta_bins[i],_eta_bins[i+1]),i,i+1)
 
            df_num_eta=df_num.Filter("pi_ETA_CMS>{}&&pi_ETA_CMS<{}".format(_eta_bins[i],_eta_bins[i+1]))
            
            heffi=TH1D("reco_effi_with_syserr_"+track_type+"_"+particle+"_{}_{}".format(_eta_bins[i],_eta_bins[i+1]),"reco_effi_with_syserr_{}_{}".format(_eta_bins[i],_eta_bins[i+1]),len(_pt_bins)-1,_pt_bins)
            hSys=TH1D("reco_effi_sys_"+track_type+"_"+particle+"_{}_{}".format(_eta_bins[i],_eta_bins[i+1]),"reco_effi_sys_{}_{}".format(_eta_bins[i],_eta_bins[i+1]),len(_pt_bins)-1,_pt_bins)
            hSys_100=TH1D("reco_effi_sys_100_"+track_type+"_"+particle+"_{}_{}".format(_eta_bins[i],_eta_bins[i+1]),"reco_effi_sys_100_{}_{}".format(_eta_bins[i],_eta_bins[i+1]),len(_pt_bins)-1,_pt_bins)
            hReleffi=TH1D("reco_effi_Relsys_"+track_type+"_"+particle+"_{}_{}".format(_eta_bins[i],_eta_bins[i+1]),"reco_effi_Relsys_{}_{}".format(_eta_bins[i],_eta_bins[i+1]),len(_pt_bins)-1,_pt_bins)
            hReleffi_100=TH1D("reco_effi_Relsys_100_"+track_type+"_"+particle+"_{}_{}".format(_eta_bins[i],_eta_bins[i+1]),"reco_effi_Relsys_100_{}_{}".format(_eta_bins[i],_eta_bins[i+1]),len(_pt_bins)-1,_pt_bins)

            for j  in range(len(_pt_bins)-1):
                df_num_eta_pt=df_num_eta.Filter("(pi_PT>{}&&pi_PT<{})".format(_pt_bins[j],_pt_bins[j+1]))
                corr_dict=df_num_eta_pt.AsNumpy(["track_corr_total_err"])
                corr_arr=corr_dict["track_corr_total_err"]
                bin=frac[i,i+1].FindBin(_pt_bins[j])
                value=frac[i,i+1].GetBinContent(bin)
                stat_err=frac[i,i+1].GetBinError(bin)
                track_corr_total_err=np.mean(corr_arr)
                total_err=(stat_err**2+track_corr_total_err**2)**0.5
                
                heffi.SetBinContent(bin,value)
                heffi.SetBinError(bin,total_err)

                hSys.SetBinContent(bin,total_err)
                hSys_100.SetBinContent(bin,total_err*100)

                hReleffi.SetBinContent(bin,total_err/value)
                hReleffi_100.SetBinContent(bin,total_err/value*100)
            
            
            f.cd()
            frac[i,i+1].Write()
            heffi.Write()
            hSys.Write()
            hSys_100.Write()
            hReleffi.Write()
            hReleffi_100.Write()
        f.Close()   
        """
    else:
        f=TFile(git_out_path+"/efficiencies/reco_effi/rfiles/"+ratio+"_"+track_type+"_"+particle+"_"+beam+"_"+pol+"_w_"+weights+"_"+data1+"_"+data2+"_"+binning_mode+"_num_denom.root","recreate")
        frac_all=TH2D(ratio+"_"+track_type+"_"+particle,"",len(_pt_bins)-1,_pt_bins,len(_eta_bins)-1,_eta_bins)
        h_num=df_num.Histo2D(RDF.TH2DModel("h_num","",len(_pt_bins)-1,_pt_bins,len(_eta_bins)-1,_eta_bins),"pi_PT","pi_ETA_CMS")
        h_denom=df_denom.Histo2D(RDF.TH2DModel("h_denom","",len(_pt_bins)-1,_pt_bins,len(_eta_bins)-1,_eta_bins),particle+"_TRUEPT",particle+"_TRUE_ETA_CMS")  
        h_num.SetDirectory(0)
        h_denom.SetDirectory(0)
        f.cd()
        print("Writing")
        h_num.Write()
        h_denom.Write()
        f.Close()


        """
           
        f=TFile(git_out_path+"/efficiencies/reco_effi/rfiles/"+ratio+"_"+track_type+"_"+particle+"_"+beam+"_"+pol+"_w_0_"+data1+"_"+data2+"_"+binning_mode+".root","recreate")

        for i in range(len(_eta_bins)-1):
            print(ratio,"_{}_{}".format(_eta_bins[i],_eta_bins[i+1]))
            frac=TH1D(ratio+"_"+track_type+"_"+particle+"_{}_{}".format(_eta_bins[i],_eta_bins[i+1]),ratio+"_{}_{}".format(_eta_bins[i],_eta_bins[i+1]),len(_pt_bins)-1,_pt_bins)
            
            df_num_eta=df_num.Filter("pi_ETA_CMS>{}&&pi_ETA_CMS<{}".format(_eta_bins[i],_eta_bins[i+1]))
            h_num=df_num_eta.Histo1D(RDF.TH1DModel("h_num_{}_{}_".format(_eta_bins[i],_eta_bins[i+1]),"h_num_{}_{}_".format(_eta_bins[i],_eta_bins[i+1]),len(_pt_bins)-1,_pt_bins),"pi_PT","tot_weight")
            #print("h_num entries",h_num.GetEntries())

            df_denom_eta=df_denom.Define("ple_TRUE_ETA_CMS",particle+"_ETA_CMS").Filter("ple_TRUE_ETA_CMS>{}&&ple_TRUE_ETA_CMS<{}".format(_eta_bins[i],_eta_bins[i+1]))
            h_denom=df_denom_eta.Histo1D(RDF.TH1DModel("h_denom_{}_{}_".format(_eta_bins[i],_eta_bins[i+1]),"h_denom_{}_{}_".format(_eta_bins[i],_eta_bins[i+1]),len(_pt_bins)-1,_pt_bins),particle+"_TRUEPT")
            #print("h_denom entries",h_denom.GetEntries())
            
            frac.Divide(h_num.GetPtr(),h_denom.GetPtr(),1.,1.,"B") 
            
            heffi=TH1D("reco_effi_with_syserr_"+track_type+"_"+particle+"_{}_{}".format(_eta_bins[i],_eta_bins[i+1]),"reco_effi_with_syserr_{}_{}".format(_eta_bins[i],_eta_bins[i+1]),len(_pt_bins)-1,_pt_bins)
            
            hSys=TH1D("reco_effi_sys_"+track_type+"_"+particle+"_{}_{}".format(_eta_bins[i],_eta_bins[i+1]),"reco_effi_sys_{}_{}".format(_eta_bins[i],_eta_bins[i+1]),len(_pt_bins)-1,_pt_bins)
            hSys_100=TH1D("reco_effi_sys_100_"+track_type+"_"+particle+"_{}_{}".format(_eta_bins[i],_eta_bins[i+1]),"reco_effi_sys_100_{}_{}".format(_eta_bins[i],_eta_bins[i+1]),len(_pt_bins)-1,_pt_bins)
            
            hReleffi=TH1D("reco_effi_Relsys_"+track_type+"_"+particle+"_{}_{}".format(_eta_bins[i],_eta_bins[i+1]),"reco_effi_Relsys_{}_{}".format(_eta_bins[i],_eta_bins[i+1]),len(_pt_bins)-1,_pt_bins)
            hReleffi_100=TH1D("reco_effi_Relsys_100_"+track_type+"_"+particle+"_{}_{}".format(_eta_bins[i],_eta_bins[i+1]),"reco_effi_Relsys_100_{}_{}".format(_eta_bins[i],_eta_bins[i+1]),len(_pt_bins)-1,_pt_bins)

            for j  in range(len(_pt_bins)-1):
                df_num_eta_pt=df_num_eta.Filter("pi_PT>{}&&pi_PT<{}".format(_pt_bins[i],_pt_bins[i+1]))
                corr_dict=df_num_eta_pt.AsNumpy(["track_corr_total_err"])
                corr_arr=corr_dict["track_corr_total_err"]
                bin=frac.FindBin(_pt_bins[j])
                value=frac.GetBinContent(bin)
                stat_err=frac.GetBinError(bin)
                track_corr_total_err=np.mean(corr_arr)
                total_err=(stat_err**2+track_corr_total_err**2)**0.5

                heffi.SetBinContent(bin,value)
                heffi.SetBinError(bin,total_err)

                hSys.SetBinContent(bin,total_err)
                hSys_100.SetBinContent(bin,total_err*100)

                hReleffi.SetBinContent(bin,total_err/value)
                hReleffi_100.SetBinContent(bin,total_err/value*100)

            f.cd()
            frac.Write()
            heffi.Write()
            hSys.Write()
            hSys_100.Write()
            hReleffi.Write()
            hReleffi_100.Write()
        """
    f.Close() 
        
import argparse

# Create an ArgumentParser object to handle command line arguments
parser = argparse.ArgumentParser()

# Add the command line arguments
parser.add_argument("-beam", "--beam", nargs='+', type=str, help="Beam type")
parser.add_argument("-pol", "--pol", nargs='+', type=str, help="Polarization")
parser.add_argument("-particle", "--particle", nargs='+', type=str, help="particle")
parser.add_argument("-binning_mode", "--binning_mode", nargs='+', type=str, help="binning_mode")
parser.add_argument("-weights", "--weights", nargs='+', type=str, help="Occupancy weight")
parser.add_argument("-MC", "--MC", nargs='+', type=str, help="Which MC sample")

parser.add_argument("-test_mode", "--test_mode", action="store_true", help="Test Mode")

# Parse the command line arguments
args = parser.parse_args()

# Extract the values of the arguments
beam = args.beam[0]
pol=args.pol[0]
probe=args.particle[0]
binning_mode=args.binning_mode[0]
weights= args.weights[0]
MC= args.MC[0]

test_mode = args.test_mode


# Print the values of the arguments
print(f"Beam type: {beam}")
print(f"Polarization : {pol}")
print(f"Particle : {probe}")
print(f"binning_mode : {binning_mode}")
print(f"Test Mode : {test_mode}")




make_ratio_histogram(beam=beam,pol=pol,data1=MC,data2="Gen"+MC,particle=probe,track_type="prompt",weights=weights,binning_mode=binning_mode,ratio="reco_effi",test_mode=test_mode)



