from ROOT import *
from ROOT import RDF
from ROOT import RDataFrame
TH1.SetDefaultSumw2(True)    
TH2.SetDefaultSumw2(True)
gStyle.SetPaintTextFormat("1.3f")
gStyle.SetPalette(1)
import sys 
import json 
from math import *
import uproot as up
sys.path.append("/home3/sara.sellam/RpPb_identified_hadrons_project")
from Binning.Binning import *
from Ntuple_path.call_calibration_ntuple  import _part,Tuple_Names_for_TagandProbe
from Helpers_function.import_dependencies import *

git_out_path="/scratch43/ssellam/results/Bkg_fractions"

git_path="/scratch43/ssellam/results"


binning_mode="oscar_ana"
beam="pPb"
probe="p"
pid="p"
mom="Lambda"
weights="5"

if probe==pid:
    mode="tight"
else:
    mode="croaser"

_pt_bins=binning[binning_mode][beam]["pt"]
_eta_bins=binning[binning_mode]["pp"]["eta"] 

effi=probe+"_as_"+pid
binning_name="oscar_ana"
m="ana2"


f=TFile(git_path+"/efficiencies/pid_effi/tag_probe_technique/rfiles/P_eta/ratio/"+beam+"/ratio_"+effi+"_"+beam+"_"+binning_name+"_"+mode+"_"+m+"_weight_"+weights+".root","open")
h_pid=f.Get("ratio")
f2=TFile(git_path+"/efficiencies/pid_effi/tag_probe_technique/rfiles/PT_eta/ratio/"+beam+"/ratio_"+effi+"_"+beam+"_"+binning_name+"_"+mode+"_"+m+"_weight_"+weights+".root","open")
_h_pid2=f2.Get("ratio")

h_pid2=TH2D("pid2","",len(_pt_bins)-1,_pt_bins,len(_eta_bins)-1,_eta_bins)

for eta in range(len(_eta_bins)-1):
    for pt in range(len(_pt_bins)-1):
        pt_center=(_pt_bins[pt]+_pt_bins[pt+1])/2
        eta_center=(_eta_bins[eta]+_eta_bins[eta+1])/2
        bin=h_pid2.FindBin(pt_center,eta_center)
        h_pid2.SetBinContent(bin,_h_pid2.GetBinContent(bin))
        h_pid2.SetBinError(bin,_h_pid2.GetBinError(bin))









file=open("/home3/sara.sellam/RpPb_identified_hadrons_project/Ntuple_path/calibration_samples.json")

data_list = json.load(file)
data_list=data_list[beam]["Data"]
tuple_name=_part[beam][mom]['Tuple']
in_tree_list=Tuple_Names_for_TagandProbe[beam][tuple_name]
in_tree=TChain()      
for data in data_list:
    for tree in in_tree_list:
        in_tree.Add(data + "/" + tree + "/DecayTree")
    assert in_tree.GetEntries()!=0
max=in_tree.GetEntries()
df=RDataFrame(in_tree)
df=df.Range(max)
dic=df.AsNumpy(columns=["probe_PT","probe_P","probe_ETA"])

pid_effi=0
trk=0


import numpy as np
pid_effi_array=np.zeros(max)
for i in range(max):#in_tree.GetEntires()):
    in_tree.GetEntry(i)
    bin=h_pid.FindBin(in_tree.probe_P,in_tree.probe_ETA)
    pid_value=h_pid.GetBinContent(bin)
    pid_effi_array[i]=pid_value

new_dict =dic
new_dict["pid_effi"]=pid_effi_array
new_df=RDF.FromNumpy(new_dict)
h_num=new_df.Histo2D(RDF.TH2DModel("h_num","",len(_pt_bins)-1,_pt_bins,len(_eta_bins)-1,_eta_bins),"probe_PT","probe_ETA","pid_effi")
h_denom=new_df.Histo2D(RDF.TH2DModel("h_denom","",len(_pt_bins)-1,_pt_bins,len(_eta_bins)-1,_eta_bins),"probe_PT","probe_ETA")
f=TFile("/home3/sara.sellam/RpPb_identified_hadrons_project/efficiencies/pid_effi/tag_probe_technique/effi_"+probe+"_as_"+pid+"_"+beam+".root","recreate")
h=TH2D("effi_with_p","",len(_pt_bins)-1,_pt_bins,len(_eta_bins)-1,_eta_bins)
h.Divide(h_num.GetPtr(),h_denom.GetPtr())
h_dev=TH2D("dev","",len(_pt_bins)-1,_pt_bins,len(_eta_bins)-1,_eta_bins)

for eta in range(len(_eta_bins)-1):
    for pt in range(len(_pt_bins)-1):
        pt_center=(_pt_bins[pt]+_pt_bins[pt+1])/2
        eta_center=(_eta_bins[eta]+_eta_bins[eta+1])/2

        bin=h.FindBin(pt_center,eta_center)
        pid_with_p=h.GetBinContent(bin)
        pid_with_p_err=h.GetBinError(bin)
        pid_with_pt=h_pid2.GetBinContent(bin)
        pid_with_pt_err=h_pid2.GetBinError(bin)
        if pid_with_pt>0 and pid_with_p>0:
            std=abs(pid_with_pt-pid_with_p)/((pid_with_pt+pid_with_p)/2)
            print("std",std)
            std_err=sqrt(pid_with_p_err**2+pid_with_pt_err**2)
        else:
            std=0
            std_err=0
        h_dev.SetBinContent(bin,std)
        h_dev.SetBinError(bin,std_err)


h.SetDirectory(0)
f.cd()
#h_num.Write()
#h_denom.Write()
h.Write()
h_pid2.Write("effi_with_pt")
h_dev.Write("dev")
#h_dev2.Write("dev2")
for i in range(len(_eta_bins)-1):
    h=h_dev.ProjectionX("dev_{}_{}".format(_eta_bins[i],_eta_bins[i+1]),i+1,i+1)
    f.cd()
    h.SetDirectory(0)
    h.Write()    
f.Close()



ylabel={"pi":r"$P^{sim}(\pi)[\%]$","K":r"$P^{sim}(K)[\%]$","p":r"$P^{sim}(p)[\%]$","noPID":r"$P^{sim}[\%]$"}
xlabel=r"$p_{\mathrm{T}}$[\mathrm{MeV}/$c$]"
beam_name={"pp":r"$pp$","pPb":r"$p\mathrm{Pb}$","Pbp":r"$\mathrm{Pb}p$"}
ple_name={"pi":r"$\pi$","K":r"$K$","p":r"$p$","noPID":"noPID"}
particles_color={"pi":"red","K":"blue","p":"green"}
particles_marke={"pi":"o","K":"^","p":"s"}
removed_bins={"pp":{"pi":{"2.0_2.5":[0,0,0,0,0,0,1,1,1],
                         "2.5_3.0":[0,0,0,0,1,1,1,1,1],
                         "3.0_3.5":[0,0,1,1,1,1,1,1,1],
                         "3.5_4.0":[0,1,1,1,1,1,1,1,1],
                         "4.0_4.3":[1,1,1,1,1,1,1,1,1],
                         "4.3_4.8":[1,1,1,1,1,1,1,0,0]},

                    "K":{"2.0_2.5":[0,0,0,0,0,0,1,1,1],
                         "2.5_3.0":[0,0,0,0,1,1,1,1,1],
                         "3.0_3.5":[0,0,1,1,1,1,1,1,1],
                         "3.5_4.0":[0,1,1,1,1,1,1,1,1],
                         "4.0_4.3":[1,1,1,1,1,1,1,1,1],
                         "4.3_4.8":[1,1,1,1,1,1,1,0,0]},

                    "p":{"2.0_2.5":[0,0,0,0,0,0,1,1,1],
                         "2.5_3.0":[0,0,0,0,1,1,1,1,1],
                         "3.0_3.5":[0,0,1,1,1,1,1,1,1],
                         "3.5_4.0":[0,1,1,1,1,1,1,1,1],
                         "4.0_4.3":[1,1,1,1,1,1,1,1,1],
                         "4.3_4.8":[1,1,1,1,1,1,1,0,0]}},

            "pPb":{"pi":{"1.6_2.0":[0,0,0,0,0,0,1,1,1],
                         "2.0_2.5":[0,0,0,0,1,1,1,1,1],
                         "2.5_3.0":[0,0,1,1,1,1,1,1,1],
                         "3.0_3.5":[0,1,1,1,1,1,1,1,1],
                         "3.5_4.0":[1,1,1,1,1,1,1,1,1],
                         "4.0_4.3":[1,1,1,1,1,1,1,0,0]},
                    
                    "K":{"1.6_2.0":[0,0,0,0,0,0,1,1,1],
                         "2.0_2.5":[0,0,0,0,1,1,1,1,1],
                         "2.5_3.0":[0,0,1,1,1,1,1,1,1],
                         "3.0_3.5":[0,1,1,1,1,1,1,1,1],
                         "3.5_4.0":[1,1,1,1,1,1,1,1,1],
                         "4.0_4.3":[1,1,1,1,1,1,1,0,0]},

                    "p":{"1.6_2.0":[0,0,0,0,0,0,1,1,1],
                         "2.0_2.5":[0,0,0,0,1,1,1,1,1],
                         "2.5_3.0":[0,0,1,1,1,1,1,1,1],
                         "3.0_3.5":[0,1,1,1,1,1,1,1,1],
                         "3.5_4.0":[1,1,1,1,1,1,1,1,1],
                         "4.0_4.3":[1,1,1,1,1,1,1,0,0]}},
                         
                "Pbp":{"pi":{"-3.0_-2.5":[0,0,0,0,0,0,1,1,1],
                            "-3.5_-3.0":[0,0,0,0,1,1,1,1,1],
                            "-4.0_-3.5":[0,0,1,1,1,1,1,1,1],
                            "-4.5_-4.0":[0,1,1,1,1,1,1,1,1],
                            "-4.8_-4.5":[1,1,1,1,1,1,1,1,1],
                            "-5.2_-4.8":[1,1,1,1,1,1,1,0,0]},

                        "K":{"-3.0_-2.5":[0,0,0,0,0,0,1,1,1],
                            "-3.5_-3.0":[0,0,0,0,1,1,1,1,1],
                            "-4.0_-3.5":[0,0,1,1,1,1,1,1,1],
                            "-4.5_-4.0":[0,1,1,1,1,1,1,1,1],
                            "-4.8_-4.5":[1,1,1,1,1,1,1,1,1],
                            "-5.2_-4.8":[1,1,1,1,1,1,1,0,0]},
                            
                        "p":{"-3.0_-2.5":[0,0,0,0,0,0,1,1,1],
                            "-3.5_-3.0":[0,0,0,0,1,1,1,1,1],
                            "-4.0_-3.5":[0,0,1,1,1,1,1,1,1],
                            "-4.5_-4.0":[0,1,1,1,1,1,1,1,1],
                            "-4.8_-4.5":[1,1,1,1,1,1,1,1,1],
                            "-5.2_-4.8":[1,1,1,1,1,1,1,0,0]}}}
                        

position={"Pbp":{"-3.0_-2.5":[0,0],
                     "-3.5_-3.0":[0,1],
                     "-4.0_-3.5":[0,2],
                     "-4.5_-4.0":[1,0],
                     "-4.8_-4.5":[1,1],
                     "-5.2_-4.8":[1,2]}}


_eta_bins=binning[binning_mode]["pp"]["eta"]
fig, axarr = plot_2_3_subplot(ylabel="Percentage difference",xlabel=r"$p_{\mathrm{T}}$[MeV/c]")

#for probe in ["pi","K","p"]:
if 1:
    out_file_name="/home3/sara.sellam/RpPb_identified_hadrons_project/efficiencies/pid_effi/tag_probe_technique/effi_"+probe+"_as_"+pid+"_"+beam+".root"
    in_file=up.open(out_file_name)
    s=0
    for i in  range(axarr.ndim):
        for j in range(3):
            if beam == "Pbp":
                i, j = position["Pbp"]["{}_{}".format(_eta_bins[s], _eta_bins[s+1])]
                ax = axarr[i][j]
            else:
                ax = axarr[i][j]
            removed_array=np.array(removed_bins["pp"][probe]["{}_{}".format(_eta_bins[s], _eta_bins[s+1])])
            center=in_file["dev_{}_{}".format(_eta_bins[s],_eta_bins[s+1])].to_boost().axes.centers[0]
            width=in_file["dev_{}_{}".format(_eta_bins[s],_eta_bins[s+1])].to_boost().axes.widths[0]/2
            y=in_file["dev_{}_{}".format(_eta_bins[s],_eta_bins[s+1])].to_numpy()[0]*100*removed_array
            yerr=in_file["dev_{}_{}".format(_eta_bins[s],_eta_bins[s+1])].errors()*100*removed_array
                       
            ax.errorbar(x=center,y=y,xerr=width,yerr=yerr,marker=particles_marke[probe], color=particles_color[probe], markersize=15, linestyle='none',label=probe+"_as_"+pid)
            ax.set_xscale("log")
            ax.set_ylim(0,70)

            ax.set_xlim(400,4619)
            ax.grid(True)
            
            ax.tick_params(axis='x', labelsize=50)  
            ax.tick_params(axis='y', labelsize=50) 
            s=s+1

ax= axarr[0][0]
leg = ax.legend()
plt.savefig("sys_p_pt_"+probe+"_as_"+pid+"_"+beam+"_weight_"+weights+".pdf")
         



