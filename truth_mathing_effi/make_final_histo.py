from ROOT  import *
from ROOT import RDataFrame
from ROOT import RDF
import ROOT
import numpy as np

import glob
import sys
sys.path.append("/home3/sara.sellam/RpPb_identified_hadrons_project")
from Binning.Binning import *
from selections.selections import _branch_list, _PID_selections,_ana_selection,_track_type, _P,_PT, _ETA_cms,_ETA , _selections,_Kinematic_range


#ROOT.EnableImplicitMT()
gStyle.SetPaintTextFormat("1.2f")
gStyle.SetPalette(1)
gROOT.ProcessLine(".x /home3/sara.sellam/RpPb_identified_hadrons_project/lhcbstyle.C++")



ROOT.gInterpreter.Declare("""
float getBinContentFromTH2D(float x, float y) 
{
    #include "TH2D.h"
    #include "TFile.h"
    #include <iostream>
    using namespace std;

    
    TFile* file = TFile::Open("/home3/sara.sellam/RpPb_identified_hadrons_project/efficiencies/truth_mathing_effi/effi_TM_pp.root");
    TH2D* hist = dynamic_cast<TH2D*>(file->Get("effi_TM"));
    int bin = hist->FindBin(x,y);
    Double_t value= hist->GetBinContent(bin);
    file ->Close();
    //cout<<value;
    return value;
}                   
""")

ROOT.gInterpreter.Declare("""
float getBinErrorFromTH2D(float x, float y) 
{
    #include "TH2D.h"
    TFile* file = new TFile("/home3/sara.sellam/RpPb_identified_hadrons_project/efficiencies/truth_mathing_effi/effi_TM_pp.root", "READ");
    TH2D* hist = dynamic_cast<TH2D*>(file->Get("effi_TM"));
    int bin = hist->FindBin(x,y);
    Double_t value = hist->GetBinError(bin);
    file ->Close();
    return value;
}                   
""")

ch={}
dic={}
df={}
h={}
binning_mode="oscar_ana"


beam="pp"

_pt_bins=binning[binning_mode][beam]["pt"]
_eta_bins=binning[binning_mode][beam]["eta"] 



_eta_bins2_dic={"oscar_ana":{"pp":{"eta":array('d',[2,2.5,3,3.5,4,4.3,4.8])},
                            "pPb":{"eta":array('d',[2.065, 2.465, 2.965, 3.465, 3.965, 4.465,4.765])},
                            "Pbp":{"eta":array('d',[2.035, 2.535, 3.035, 3.535, 3.835, 4.035, 4.335, 4.735])},
                            }}
        
    
    
_eta_bins2=_eta_bins2_dic[binning_mode][beam]["eta"] 


list_data={}

path={"oscar":{"pp":"/scratch41/oscar/pp/RealData/RealData_v20_pp_MD_5TeV_Job*_PseudoIPBranch_withDOCA.root",
               "pPb":"/scratch41/oscar/pA/RealData/RealData_v20_pA_MD_5TeV_Job*_PseudoIPBranch_withDOCA.root",
               "Pbp":"/scratch41/oscar/pA/RealData/RealData_v20_Ap_MD_5TeV_Job*_PseudoIPBranch_withDOCA.root"},
       "sara":{"pp":"/scratch38/HITuple/Data/reduced_Data_2013/pp_Down/*_HITuple_pp_magdown.root",
               "pPb":"/scratch38/HITuple/Data/reduced_Data_2013/pPb_Down/*_HITuple_pA_magdown.root",
                "Pbp":"/scratch38/HITuple/Data/reduced_Data_2013/Pbp_Down/*_HITuple_Ap_magdown.root"}}


cuts={"pseudoIP":{"pp":0.368,
                 "pPb":0.348,
                  "Pbp":0.348},
        "GhostP":{"pp":0.078,
                 "pPb":0.103,
                  "Pbp":0.109}}
list_data["oscar"]=glob.glob(path["oscar"][beam])

ch["oscar"]=TChain("")

for l in list_data["oscar"]:#[list_data["oscar"][0],list_data["oscar"][1],list_data["oscar"][2],list_data["oscar"][3],list_data["oscar"][4],list_data["oscar"][5],list_data["oscar"][6]]:
#for l in [list_data["oscar"][0],list_data["oscar"][1],list_data["oscar"][2],list_data["oscar"][3],list_data["oscar"][4],list_data["oscar"][5],list_data["oscar"][6]]:
    ch["oscar"].Add(l+"/pAanalysis/pAanalysis")
df["oscar"]=RDataFrame(ch["oscar"])




ch["sara"]=TChain("")
ch_pseudoIP=TChain("")
list_data["sara"]=glob.glob(path["sara"][beam])
#for l in list_data["sara"]:#[list_data["sara"][0],list_data["sara"][1],list_data["sara"][2],list_data["sara"][3],list_data["sara"][4]]:#,list_data["sara"][5],list_data["sara"][6]]:
for l in [list_data["sara"][0]]:#,list_data["sara"][1],list_data["sara"][2],list_data["sara"][3],list_data["sara"][4]]:#,list_data["sara"][5],list_data["sara"][6]]:
    ch["sara"].Add(l+"/DecayTree")
    file_name = l.split('/')[-1]
    ch_pseudoIP.Add("/scratch38/HITuple/Data/pseudoIP/"+beam+"_Down/pseudoip_"+file_name+"/tree")

assert ch["sara"].GetEntries()==ch_pseudoIP.GetEntries()
ch_pseudoIP.AddFriend(ch["sara"])
df["sara"]=RDataFrame(ch_pseudoIP)


f=TFile("/scratch43/ssellam/results/efficiencies/TM_effi/effi_TM_candidates_"+beam+".root","recreate")



if beam =="pp": 
    sara_cuts="pi_TRACK_Type==3&&pi_ETA>2&&pi_ETA<4.8&&pi_P>2000&&pi_PT>400&&pi_PT<8000&&pi_Hlt1NoBiasLeadingCrossingDecision_Dec==1&&BCType==3"
else:
    sara_cuts="pi_TRACK_Type==3&&pi_ETA>2&&pi_ETA<4.8&&pi_P>2000&&pi_PT>400&&pi_PT<8000&&pi_Hlt1MBMicroBiasVeloDecision_Dec==1&&BCType==3&&nPVs==1"

sara_cuts=sara_cuts+"&&pi_TRACK_GhostProb<"+str(cuts["GhostP"][beam])+"&&pseudoIP<"+str(cuts["pseudoIP"][beam])


df["sara"]=df["sara"].Define("pi_ETA","-log(tan(asin(pi_PT/pi_P)/2))")
if beam =="pPb":
    df["sara"]=df["sara"].Define("pi_ETA_cms","pi_ETA-0.5")
elif beam =="Pbp":
    df["sara"]=df["sara"].Define("pi_ETA_cms","(pi_ETA*-1)-0.5")


df["sara"]=df["sara"].Define("TM_effi_values","getBinContentFromTH2D(pi_PT,pi_ETA)")
df["sara"]=df["sara"].Define("TM_effi_errors","getBinErrorFromTH2D(pi_PT,pi_ETA)")

#print("I am here")
h["sara-weighted"]=df["sara"].Filter(sara_cuts).Histo2D(RDF.TH2DModel("h_weighted","",len(_pt_bins)-1,_pt_bins,len(_eta_bins2)-1,_eta_bins2),"pi_PT","pi_ETA","TM_effi_values")
#print("I am here 2")

h["sara"]=df["sara"].Filter(sara_cuts).Histo2D(RDF.TH2DModel("h","",len(_pt_bins)-1,_pt_bins,len(_eta_bins2)-1,_eta_bins2),"pi_PT","pi_ETA")
#print("I am here 3")
"""

for eta in range(len(_eta_bins2)-1):
    for pt in range(len(_pt_bins)-1):
        bin=h["sara"].FindBin(_pt_bins[pt],_eta_bins2[eta])
        print(h["sara"].GetBinContent(bin))
"""



h["sara-weighted"].SetDirectory(0)
f.cd()
h["sara-weighted"].Write()

"""
print("I am here 4")
h["sara-weighted"].SetDirectory(0)
f.cd()
h["sara-weighted"].Write()
print("I am here 5")
"""
f.Close()
"""


h["sara"]=df["sara"].Filter(sara_cuts).Histo2D(RDF.TH2DModel("h_sara_bins2","",len(_pt_bins)-1,_pt_bins,len(_eta_bins2)-1,_eta_bins2),"pi_PT","pi_ETA")
if beam !="pp":
    h["sara_cms"]=df["sara"].Filter(sara_cuts).Histo2D(RDF.TH2DModel("h_sara","",len(_pt_bins)-1,_pt_bins,len(_eta_bins)-1,_eta_bins),"pi_PT","pi_ETA_cms")


h["sara"].SetDirectory(0)
f.cd()
h["sara"].Write("sara_cand")

if beam !="pp":
    h["sara_cms"].SetDirectory(0)
    f.cd()
    h["sara_cms"].Write("sara_cand_cms")


scale = 1/(h["sara"].Integral())

assert h["sara"].GetEntries()!=0
assert h["sara"].Integral()!=0
h["sara"].Scale(scale)
h["sara"].SetDirectory(0)
f.cd()
h["sara"].Write("sara_cand_scaled")

for pid in ["pi","K","p"]:
    sara_cuts_pid=sara_cuts+"&&"+_PID_selections[pid][beam]
    h["sara_"+pid]=df["sara"].Filter(sara_cuts_pid).Histo2D(RDF.TH2DModel("h_sara_"+pid,"",len(_pt_bins)-1,_pt_bins,len(_eta_bins2)-1,_eta_bins2),"pi_PT","pi_ETA")   
    if beam !="pp":
        h["sara_"+pid+"_cms"]=df["sara"].Filter(sara_cuts_pid).Histo2D(RDF.TH2DModel("h_sara_"+pid+"_cms","",len(_pt_bins)-1,_pt_bins,len(_eta_bins)-1,_eta_bins),"pi_PT","pi_ETA_cms")   
        h["sara_"+pid+"_cms"].SetDirectory(0)
        f.cd()
        h["sara_"+pid+"_cms"].Write("sara_cand_"+pid+"_cms")
    h["sara_"+pid].SetDirectory(0)
    f.cd()
    h["sara_"+pid].Write("sara_cand_"+pid)


"""


