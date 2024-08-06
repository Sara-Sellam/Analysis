
# run with lb-conda pidcalib bash

import sys
import subprocess 
sys.path.append("/home3/sara.sellam/RpPb_identified_hadrons_project")
from Binning.Binning import *

from selections.pid_selections import _pid_param,_params
import glob


#subprocess.call(["lb-conda pidcalib bash"], shell=True)


ple={"pi":"Pi","K":"K","p":"P"}
selections=["DLLK<"+_pid_param["pp"]["pi"][0]+"&DLLp<"+_pid_param["pp"]["pi"][1],"DLLK>"+_pid_param["pp"]["K"][0]+"&(DLLK-DLLp)>"+_pid_param["pp"]["K"][1],"DLLp>"+_pid_param["pp"]["p"][0]+"&(DLLp-DLLK)>"+_pid_param["pp"]["p"][1]]
extra_selections="TRACK_GHOSTPROB<"+_params["GhostP"]["pp"]+"&P>2000&PT>400&PT<8000"
selection_for_ple_to_ghost="TRACK_GHOSTPROB>"+_params["GhostP"]["pp"]
#selection_for_ple_to_ghost="TRACK_GHOSTPROB > 0.2 & P>2000&PT>400&PT<8000" 

if 0: 
    for k in ple.keys():
        for cut in selections:
            index_ple=selections.index(cut)
            index_cut=list(ple.keys()).index(k)
            if 0:#index_ple==index_cut:
                subprocess.call(["pidcalib2.make_eff_hists --sample Turbo15 --magnet down --particle "+ple[k]+" --pid-cut '"+cut+"' --cut '"+extra_selections+"' --bin-var PT --bin-var ETA --binning-file binning_tight.json --output-dir pidcalib_output"],shell=True)
            else:
                subprocess.call(["pidcalib2.make_eff_hists --sample Turbo15 --magnet down --particle "+ple[k]+" --pid-cut '"+cut+"' --cut '"+extra_selections+"' --bin-var PT --bin-var ETA --binning-file binning_croaser.json --output-dir pidcalib_output"],shell=True)

                
    for i in range(len(ple.keys())): 
        k=list(ple.keys())[i]
        cut="HasRich==1"  
        subprocess.call(["pidcalib2.make_eff_hists --sample Turbo15 --magnet down --particle "+ple[k]+" --pid-cut '"+selection_for_ple_to_ghost+"' --bin-var PT --bin-var ETA --binning-file binning.json --output-dir pidcalib_for_ghost_misID_output"],shell=True)
    pkl_list=glob.glob("pidcalib_output/*.pkl")
    pkl_list2=glob.glob("pidcalib_for_ghost_misID_output/*.pkl")
    for f in pkl_list :
        subprocess.call(["pidcalib2.pklhisto2root '"+f+"'"],shell=True)
    for f2 in pkl_list2:     
        subprocess.call(["pidcalib2.pklhisto2root '"+f2+"'"],shell=True)




from ROOT import *
gROOT.ProcessLine(".x ~/lhcbStyle.C")
gStyle.SetPaintTextFormat("1.3f")
gStyle.SetOptTitle(1)
gStyle.SetPadRightMargin(0.13)
import uproot as up
name={"pi":"#pi","K":"K","p":"#font[12]{p}"}
binning_mode="sara_ana"
beam="pp"
_pt_bins=binning[binning_mode]["pp"]["pt"]
_eta_bins=binning[binning_mode]["pp"]["eta"] 

for k in ple.keys():
    for cut in selections:
        f=up.open("pidcalib_output/effhists-Turbo15-down-"+ple[k]+"-"+cut+"-PT.ETA.root")    
        i=selections.index(cut)
        h_name=f.keys()[0]
        key=list(ple.keys())
        in_file=TFile("pidcalib_output/effhists-Turbo15-down-"+ple[k]+"-"+cut+"-PT.ETA.root","open")
        h1=in_file.Get(h_name)
        c1=TCanvas("c1","c1",500,400)
        pad1 = TPad("pad","pad", 0.01, 0.01, 0.95, 0.95)
        pad1.Draw()
        pad1.cd()
        gStyle.SetPalette(1)
        h1.SetStats(0)
        h1.SetYTitle("#eta")
        h1.SetTitle("#varepsilon_{ "+name[k]+"#rightarrow "+name[key[int(i)]]+"}")
        h1.SetXTitle("p_{T}[MeV/c]")
        h1.Draw("colz,text,e")
        c1.SaveAs("/home3/sara.sellam/RpPb_identified_hadrons_project/efficiencies/pid_effi/pp/Plots/"+k+"_to_"+key[int(i)]+"_pt_eta.pdf")
        f_out=TFile("/home3/sara.sellam/RpPb_identified_hadrons_project/efficiencies/pid_effi/rfiles/pid_"+k+"_as_"+key[int(i)]+"_"+beam+"_"+binning_mode+".root","recreate")
        for e in range(len(_eta_bins)-1):
            h_pid=TH1D(k+"_as_"+key[int(i)]+"_{}_{}".format(_eta_bins[e],_eta_bins[e+1]),k+"_as_"+key[int(i)]+"_{}_{}".format(_eta_bins[e],_eta_bins[e+1]),len(_pt_bins)-1,_pt_bins)
            for j in range(len(_pt_bins)-1):
                _pt_centers=(_pt_bins[j]+_pt_bins[j+1])/2
                _eta_centers=(_eta_bins[e]+_eta_bins[e+1])/2
                bin=h1.FindBin(_pt_centers,_eta_centers)
                content=h1.GetBinContent(bin)
                print(k+"_as_"+key[int(i)],"_pt_centers",_pt_centers,"_eta_centers",_eta_centers,"content",content)
                if content:
                    h_pid.SetBinContent(j+1,h1.GetBinContent(bin))
                    h_pid.SetBinError(j+1,h1.GetBinError(bin))
                else:
                    h_pid.SetBinContent(j+1,0)
                    h_pid.SetBinError(j+1,0)
            f_out.cd()
            h_pid.Write()
        f_out.Close()

"""
for i in range(3):
    cut=selection_for_ple_to_ghost
    p=list(ple.keys())[i]
    f=up.open("pidcalib_for_ghost_misID_output/effhists-Turbo15-down-"+ple[p]+"-"+cut+"-PT.ETA.root")    
    h_name=f.keys()[0]
    in_file=TFile("pidcalib_for_ghost_misID_output/effhists-Turbo15-down-"+ple[p]+"-"+cut+"-PT.ETA.root","open")
    h1=in_file.Get(h_name)
    c1=TCanvas("c1","c1",500,400)
    pad1 = TPad("pad","pad", 0.01, 0.01, 0.95, 0.95)
    pad1.Draw()
    pad1.cd()
    gStyle.SetPalette(1)
    h1.SetStats(0)
    h1.SetYTitle("#eta")
    h1.SetTitle("#varepsilon_{ "+name[p]+"#rightarrow g}")
    h1.SetXTitle("p_{T}[MeV/c]")
    h1.Draw("colz,text,e")
    c1.SaveAs("/home3/sara.sellam/RpPb_identified_hadrons_project/efficiencies/pid_effi/pp/Plots/"+p+"_to_g_pt_eta.pdf")
    f_out=TFile("/home3/sara.sellam/RpPb_identified_hadrons_project/efficiencies/pid_effi/rfiles/pid_"+p+"_as_g_"+beam+"_"+binning_mode+".root","recreate")
    for e in range(len(_eta_bins)-1):
        h_pid=TH1D(p+"_as_g_{}_{}".format(_eta_bins[e],_eta_bins[e+1]),p+"_as_g_{}_{}".format(_eta_bins[e],_eta_bins[e+1]),len(_pt_bins)-1,_pt_bins)
        for j in range(len(_pt_bins)-1):
            _pt_centers=(_pt_bins[j]+_pt_bins[j+1])/2
            _eta_centers=(_eta_bins[e]+_eta_bins[e+1])/2
            bin=h1.FindBin(_pt_centers,_eta_centers)
            content=h1.GetBinContent(bin)
            print(content)
            if content:
                h_pid.SetBinContent(j+1,h1.GetBinContent(bin))
                h_pid.SetBinError(j+1,h1.GetBinError(bin))
            else:
                h_pid.SetBinContent(j+1,0)
                h_pid.SetBinError(j+1,0)
        f_out.cd()
        h_pid.Write()
    f_out.Close()
"""