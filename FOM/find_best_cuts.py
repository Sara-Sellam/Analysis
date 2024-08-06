import numpy as np
import sys
sys.path.append("/home3/sara.sellam/RpPb_identified_hadrons_project")
import glob
import numpy as np
import pandas as pd
import time
from ROOT import *
from ROOT import RDF
from ROOT import RDataFrame
import matplotlib.pyplot as plt

import mplhep as hep

from selections.selections import *
from Binning.Binning import *
from selections.selections import _branch_list, _ana_selection_without_PID,_ana_selection,_track_type, _P,_PT, _ETA_cms, _selections,_Kinematic_range
from Ntuple_path.call_ntuple import _path, _corrections,get_data
from ple_yield.utils import *
ROOT.EnableImplicitMT(20)




def significance(s,b):
    return s/sqrt(b)


def significance2(s,b):
    return s/sqrt(s+b)

def significance3(s,b):
    return s/b

def significance_x_purity(s,b):
    s*s/pow((s+b),3/2)
    return s/b


def make_histogram(beam="",pol="",data="MC",particle="",binning_mode="oscar_ana",test_mode="",weights="",_cuts_modes="TrueID_cuts",mode="Significance x Purity"):
    start_time = time.time()
    ch_data=TChain()
    pseudoIP_ch=TChain()
    ch_data_weights=TChain()
    ratio="N_cand"
    df={}
    PIDpK={}
    PIDp={}
    PIDKp={}
    PIDK={}
    if data!="Data":
        mc_name=data.split("-")[0]
        
        if mc_name=="MC":
            if not _cuts_modes:
                df_data,cuts=process_mc_data(data=data,beam=beam, pol=pol, particle=particle, weights=weights,test_mode=test_mode,nbr_files=1)
            if _cuts_modes =="TrueID_pid_cuts":
                df_data,cuts=process_mc_data(data=data,beam=beam, pol=pol, particle=particle, weights=weights,test_mode=test_mode,nbr_files=1)
            if _cuts_modes =="TrueID_cuts":
                df_data,cuts=process_TrueID_selected_data(data=data,beam=beam, pol=pol, particle=particle, weights=weights,test_mode=test_mode,nbr_files=1)
            if  _cuts_modes =="pid_cuts":
                df_data,cuts=process_pid_selected_data(data=data,beam=beam, pol=pol, particle=particle, weights=weights,test_mode=test_mode,nbr_files=1)
            if _cuts_modes =="pid_cuts_all":
                df_data,cuts=process_pid_selected_allcand_data(data=data,beam=beam, pol=pol, particle=particle, weights=weights,test_mode=test_mode,nbr_files=1)
        elif mc_name=="GenMC":
            df_data,cuts=process_gen_mc_data(data=data,beam=beam, pol=pol, particle=particle, weights=weights,test_mode=test_mode,nbr_files=1)
    else:
        df_data,cuts=process_data(data=data, beam=beam, pol=pol, particle=particle, test_mode=test_mode,nbr_files=1)
    df={}
    h={}
    df_data=df_data.Define("pi_PIDKp","pi_PIDK-pi_PIDp")
    df_data=df_data.Define("pi_PIDpK","pi_PIDp-pi_PIDK")
    df_data=df_data.Define("pi_PIDK_negative","pi_PIDK*-1")
    df_data=df_data.Define("pi_PIDp_negative","pi_PIDp*-1")
    _pt_bins=binning[binning_mode][beam]["pt"]




    particle_mapping = {
    "pi": ("K", "p"),
    "K": ("pi", "p"),
    "p": ("K", "pi")}

    particle_pdg_ID = {
    "pi": ("211"),
    "K": ("321"),
    "p": ("2212")}


    """
    PIDh_label= {
    "pi": (r"$\Delta \log \mathcal{L}_{K-\pi}<$", r"$\Delta log \mathcal{L}_{p-\pi}<$"),
    "K":(r"$\Delta \log \mathcal{L}_{K-\pi}>$", r"$\Delta log \mathcal{L}_{K-p}>$"),
    "p": (r"$\Delta \log \mathcal{L}_{p-\pi}>$",r"$\Delta log \mathcal{L}_{p-K}>$")}
    """
    PIDh_label= {
    "pi": (r"$DLL_{K\pi}<$", r"$DLL_{p\pi}<$"),
    "K":(r"$DLL_{K\pi}>$", r"$DLL_{Kp}>$"),
    "p": (r"$DLL_{p\pi}>$",r"$DLL_{pK}>$")}




    PIDh_cuts = {
    "pi": ("pi_PIDK<", "pi_PIDp<"),
    "K":("pi_PIDK>", "pi_PIDKp>"),
    "p": ("pi_PIDp>", "pi_PIDpK>")}


    # Optimization over cut values
    best_fom = 0
    best_cuts = (None, None)
    

    particle2, particle3 = particle_mapping.get(particle, (None, None))
    PIDh_x,PIDh_y=PIDh_label.get(particle, (None, None))
    print("sig",cuts)

    bkg_cuts=_Kinematic_range[mc_name]['pi'] +"&&pi_MC_isPrompt==1 &&(abs(pi_TRUEID)==321 || abs(pi_TRUEID)==2212|| abs(pi_TRUEID)==211)" # &&(("+ _ana_selection[mc_name][beam][pol]["TrueID_"+particle2]+")||("+ _ana_selection[mc_name][beam][pol]["TrueID_"+particle3]+"))"
    
    bkg_cuts=_Kinematic_range[mc_name]['pi'] +"&&pi_MC_isPrompt==1 &&(abs(pi_TRUEID)==321 || abs(pi_TRUEID)==2212|| abs(pi_TRUEID)==211)"



    pid=particle_pdg_ID.get(particle, (None))
    pid_cut1,pid_cut2=PIDh_cuts.get(particle, (None, None))


    cuts=_Kinematic_range[mc_name]['pi'] +"&&pi_MC_isPrompt==1 && abs(pi_TRUEID)=="+pid

    df["sig"]=df_data.Filter(cuts)  
    df["bkg"]= df_data.Filter(bkg_cuts) 


    
    cut1_range = np.linspace(-20, 10, 20)
    cut2_range = np.linspace(-20, 10, 20)
    f=TFile("rfiles/"+mode+"_"+particle+"_FoM_"+beam+"_"+data+".root","recreate")
    FOM=np.zeros((len(cut1_range), len(cut2_range)))
    max_fom = 0
    best_cut1 = None
    best_cut2 = None

    for i,cut1 in enumerate(cut1_range):
        for j,cut2 in enumerate(cut2_range):
            """
            h["bkg"]=df["bkg"].Filter(PIDh_x+"<"+str(cut1)+"&&"+PIDh_y+"<"+str(cut2)).Histo1D(RDF.TH1DModel("h_bkg","",100,400,5000),"pi_PT","weights")
            h["sig"]=df["sig"].Filter(PIDh_x+"<"+str(cut1)+"&&"+PIDh_y+"<"+str(cut2)).Histo1D(RDF.TH1DModel("h_sig","",100,400,5000),"pi_PT","weights")
            """
            h["bkg"]=df["bkg"].Filter(pid_cut1+str(cut1)+"&&"+pid_cut2+str(cut2)).Histo1D(RDF.TH1DModel("h_bkg","",100,400,5000),"pi_PT","weights")
            h["sig"]=df["sig"].Filter(pid_cut1+str(cut1)+"&&"+pid_cut2+str(cut2)).Histo1D(RDF.TH1DModel("h_sig","",100,400,5000),"pi_PT","weights")

            h["bkg"].SetDirectory(0)
            h["sig"].SetDirectory(0)
            f.cd()
            h["bkg"].Write("bkg_"+pid_cut1+str(cut1)+"&&"+pid_cut2+str(cut2))
            h["sig"].Write("sig_"+pid_cut1+str(cut1)+"&&"+pid_cut2+str(cut2))
            sum_s_b=h["bkg"].Clone()
            sum_s_b.Add(h["sig"].GetPtr())
            purity=h["sig"].Clone()
            #purity.Divide(h["sig"].GetPtr(),sum_s_b)
            purity.Divide(h["sig"].GetPtr(),h["bkg"].GetPtr())
            purity.SetDirectory(0)
            f.cd()
            purity.Write("purity_"+pid_cut1+str(cut1)+"&&"+pid_cut2+str(cut2))
            # Calculate the integral of the signal histogram
            total_signal = h["sig"].Integral()
            total_bkg = h["bkg"].Integral()


            # Calculate the integral of the combined signal and background histogram
            if mode=="Purity":
                # Calculate overall purity

                fom = total_signal / total_bkg if total_bkg > 0 else 0


                FOM[i, j] = fom*100
                if fom > max_fom:
                    max_fom = fom
                    best_cut1 = cut1
                    best_cut2 = cut2

            if mode=="Significance x Purity":
                fom=significance_x_purity(total_signal,total_bkg)
                FOM[i, j] = fom

                if fom > max_fom:
                    max_fom = fom
                    best_cut1 = cut1
                    best_cut2 = cut2


    plt.style.use(hep.style.LHCb2)
    plt.rcParams["text.usetex"] =True
    plt.figure(figsize=(10, 8))
    plt.contourf(cut1_range, cut2_range, FOM.T, levels=50, cmap='viridis')
    plt.colorbar()
    plt.xlabel(PIDh_x,fontsize=20)
    plt.ylabel(PIDh_y,fontsize=20)
    plt.title(mode+' of PID Cuts[\%]')

    print("best cuts",pid_cut1+str(best_cut1),pid_cut2+str(best_cut2))
    #plt.text(best_cut1, best_cut2, f'[{best_cut1:.2f}, {best_cut2:.2f}]', 
    #    ha='center', va='center', color='white', fontsize=8, bbox=dict(facecolor='black', alpha=0.5))

    plt.savefig("Plots/"+mode+"_"+particle+"_FoM_"+beam+"_"+data+".pdf")
 
    f.Close()

for b in ["pp"]:
    for ple in ["pi","K","p"]:
        make_histogram(beam=b,pol="Down",data="MC-Sim09d",particle=ple,binning_mode="oscar_ana",test_mode=True,weights="1",_cuts_modes="TrueID_cuts",mode="Purity")


