
import numpy as np
import pandas as pd
import sys
sys.path.append("/home3/sara.sellam/RpPb_identified_hadrons_project")
from ROOT import *
from ROOT import RDF
from ROOT import RDataFrame
import glob
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


def make_ratio_histogram(beam="",pol="",data1="",particle="",weights=None,binning_mode="",ratio="",test_mode=False):
    h_num,h_denom,frac={},{},{}

    ################################


    mc_input_files_list=get_data(_path[data1][beam][pol][0],test_mode,num_files=10)
    ch=TChain()
    pseudoIP_ch=TChain()
    track_reco_corr_ch=TChain()
    ch_weights=TChain()
    for idx, mc_input_files in enumerate(mc_input_files_list):
        for in_file in mc_input_files:
            ch.Add(in_file+"/"+_path[data1][beam][pol][1][particle])
            in_file_name=in_file.split('/')[-1]
            pseudoIP_ch.Add(_corrections["pseudoIP"][data1][beam][pol][0][idx]+"pseudoip_"+in_file_name+"/"+_corrections["pseudoIP"][data1][beam][pol][1])
            track_reco_corr_ch.Add(_corrections["track_reco_corr_sys"][data1][beam][pol][0][idx]+"track_reco_corr_"+in_file_name+"/DecayTree")
            
            if weights!="0":
                ch_weights.Add(_corrections["weights"+weights][data1][beam][pol][0][idx]+"w_"+weights+"_"+in_file_name+"/"+_corrections["weights"+weights][data1][beam][pol][1])

    assert ch.GetEntries()== pseudoIP_ch.GetEntries()
    assert pseudoIP_ch.GetEntries()== track_reco_corr_ch.GetEntries()
    pseudoIP_ch.AddFriend(ch)
    ch=pseudoIP_ch
    if weights!="0": 
        if beam=="pp":
            ch.AddFriend(ch_weights)
            df=RDataFrame(ch)
        else:
            ch.AddFriend(ch_weights)
            df=RDataFrame(ch)
            
        assert track_reco_corr_ch.GetEntries()==ch_weights.GetEntries()
    else:
        df=RDataFrame(ch)
   

    
    

   
    df=df.Define("pi_ETA_CMS",_ETA_cms[beam]["MC"][particle])

    if weights!="0":
        df=df.Define("tot_weight","w")
   



    
    _pt_bins=binning[binning_mode][beam]["pt"]
    _eta_bins=binning[binning_mode][beam]["eta"] 
    print("requested binning",binning_mode)
    print("_pt_bins",_pt_bins)
    print("_eta_bins",_eta_bins)

    cuts={"pseudoIP":{"pp":0.368,
                 "pPb":0.348,
                  "Pbp":0.348},
        "GhostP":{"pp":0.078,
                 "pPb":0.103,
                  "Pbp":0.109}}
    
    cuts_num="pi_TRACK_Type==3&&BCType==3&&pi_Reconstructed==1&&pi_ETA>2&&pi_ETA<4.8&&pi_P>7000&&pi_PT>500&&pi_PT<8000&&pi_TRACK_GhostProb<"+str(cuts["GhostP"][beam])+"&&pseudoIP<"+str(cuts["pseudoIP"][beam])#+"&&pi_MC_isPrompt==1"

    cuts_denom="pi_TRUEID!=0&&pi_TRACK_Type==3&&BCType==3&&pi_Reconstructed==1"

    df_num=df.Filter(cuts_num)
    df_denom=df.Filter(cuts_denom)


    #frac computations
    if weights!="0": 
        f=TFile(git_out_path+"/efficiencies/sel_effi/rfiles/"+ratio+"_"+beam+"_"+pol+"_w_"+weights+"_"+data1+"_"+binning_mode+"_num_denom_new2.root","recreate")
        frac_all=TH2D(ratio,"",len(_pt_bins)-1,_pt_bins,len(_eta_bins)-1,_eta_bins)
        h_num=df_num.Histo2D(RDF.TH2DModel("h_num","",len(_pt_bins)-1,_pt_bins,len(_eta_bins)-1,_eta_bins),"pi_PT","pi_ETA_CMS","tot_weight")
        h_denom=df_denom.Histo2D(RDF.TH2DModel("h_denom","",len(_pt_bins)-1,_pt_bins,len(_eta_bins)-1,_eta_bins),"pi_PT","pi_ETA_CMS","tot_weight")
        h_num.SetDirectory(0)
        h_denom.SetDirectory(0)
        print("dividing")
        frac_all.Divide(h_num.GetPtr(),h_denom.GetPtr(),1,1,"B")
        f.cd()
        print("Writing")
        h_num.Write()
        h_denom.Write()
        frac_all.Write()
        f.Close()
       
    else:
        f=TFile(git_out_path+"/efficiencies/sel_effi/rfiles/"+ratio+"_"+beam+"_"+pol+"_w_"+weights+"_"+data1+"_"+binning_mode+"_num_denom.root","recreate")
        frac_all=TH2D(ratio,"",len(_pt_bins)-1,_pt_bins,len(_eta_bins)-1,_eta_bins)
        h_num=df_num.Histo2D(RDF.TH2DModel("h_num","",len(_pt_bins)-1,_pt_bins,len(_eta_bins)-1,_eta_bins),"pi_PT","pi_ETA_CMS","tot_weight")
        h_denom=df_denom.Histo2D(RDF.TH2DModel("h_denom","",len(_pt_bins)-1,_pt_bins,len(_eta_bins)-1,_eta_bins),"pi_PT","pi_ETA_CMS","tot_weight")
        h_num.SetDirectory(0)
        h_denom.SetDirectory(0)
        f.cd()
        print("Writing")
        h_num.Write()
        h_denom.Write()
        frac_all.Write()
        f.Close()


        

"""
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


"""

#make_ratio_histogram(beam=beam,pol=pol,data1=MC,data2="Gen"+MC,particle=probe,track_type="prompt",weights=weights,binning_mode=binning_mode,ratio="sel_effi",test_mode=test_mode)



make_ratio_histogram(beam="pPb",pol="Down",data1="MC-Sim09e",particle="pi",weights="1",binning_mode="sara_ana",ratio="sel_effi",test_mode=False)
