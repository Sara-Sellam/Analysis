import ROOT
from ROOT import *
import numpy as np 
import os,sys
import glob
from ROOT import RDataFrame
Minimizer = ROOT.RooFit.Minimizer
ROOT.EnableImplicitMT()

gStyle.SetPaintTextFormat("1.3f")
gStyle.SetPalette(1)
gROOT.ProcessLine('.x /home3/sara.sellam/RpPb_identified_hadrons_project/efficiencies/sel_effi/fit_mass/Bkg.cxx++')

###########################################
# in this Ntuples:
# Kpl is the tag 
# Kmi is the prob 


############################################


beam="pPb"

ch=TChain()
if beam =="pPb":
    ch.Add("/scratch43/ssellam/calibration/MC/Phi_KK/"+beam+"/*"+beam+"*.root/PhiMCTuple/DecayTree")
else:
    list_file=glob.glob("/scratch43/ssellam/fMC/phi_for_TM/"+beam+"_Up/*"+beam+"*.root")
    for i in range(len(list_file)):
        ch.Add(list_file[i]+"/PhiTuple/DecayTree")

assert ch.GetEntries()!=0
if beam =="pPb":
    df=RDataFrame(ch)
else:
    df=RDataFrame(ch)
    df=df.Define("Kmi_ETA","-log(tan(asin(Kmi_PT/Kmi_P)/2))")

df=df.Filter("abs(phi_TRUEID)==333&&abs(Kpl_TRUEID)==321&&abs(Kmi_TRUEID)==321&&abs(Kpl_MC_MOTHER_ID)==333&&abs(Kmi_MC_MOTHER_ID)==333")
probe_cuts="Kmi_TRACK_Type==3&&Kmi_P>7000"
tag_cuts="Kpl_ProbNNk>0.5"
df=df.Filter(tag_cuts)





pt_bins=[500,961,1249,4619]
eta_bins=[2,3.5,4,4.8] 
sig_dict={}
mode="ana" #sys #ana

cuts={"pseudoIP":{"pp":0.368,
                 "pPb":0.348,
                  "Pbp":0.348},
        "GhostP":{"pp":0.078,
                 "pPb":0.103,
                  "Pbp":0.109}}


#sel_cuts="Kpl_Reconstructed==1&&Kpl_TRACK_Type==3&&Kpl_P>7000&&BCType==3&&nPVs==1&&Kpl_TRACK_GhostProb<"+str(cuts["GhostP"][beam])#+"&&Kpl_IPCHI2_ORIVX<"+str(cuts["pseudoIP"][beam])
#selection="sara_cuts"
#sel_cuts="Kpl_TRACK_GhostProb<"+str(cuts["GhostP"][beam])+"&&Kpl_IPCHI2_ORIVX<"+str(cuts["pseudoIP"][beam])

selection="oscar_cuts"
#sel_cuts="Kpl_Reconstructed==1&&Kpl_TRACK_Type==3&&Kpl_P>2000&&BCType==3&&nPVs==1&&Kpl_TRACK_GhostProb<"+str(cuts["GhostP"][beam])#+"&&Kpl_IPCHI2_ORIVX<"+str(cuts["pseudoIP"][beam])

sel_cuts=probe_cuts+"&&Kmi_TRACK_GhostProb<"+str(cuts["GhostP"][beam])#+"&&Kpl_IPCHI2_ORIVX<"+str(cuts["pseudoIP"][beam])

for i in range(len(eta_bins)-1):
    sig_dict[str(eta_bins[i])+"_"+str(eta_bins[i+1])]={}
    for j in range(len(pt_bins)-1):
        sig_dict[str(eta_bins[i])+"_"+str(eta_bins[i+1])][str(pt_bins[j])+"_"+str(pt_bins[j+1])]={}

for i in range(len(eta_bins)-1):
    for j in range(len(pt_bins)-1):
        sig_dict[str(eta_bins[i])+"_"+str(eta_bins[i+1])][str(pt_bins[j])+"_"+str(pt_bins[j+1])]["before"]={}


        df_sel=df.Filter("Kmi_PT>"+str(pt_bins[j])+"&&Kmi_PT<"+str(pt_bins[j+1])+"&&Kmi_ETA>"+str(eta_bins[i])+"&&Kmi_ETA<"+str(eta_bins[i+1]))
    
        
        mass = ROOT.RooRealVar( "phi_MM","#phi(1020)[MeV/c^{2}]",1000,1040) 
        
        _mass_branch= df_sel.AsNumpy(["phi_MM"])
        arr_mass_branch= _mass_branch["phi_MM"]
        bins = [np.linspace( 1000,1040, 100)]
        counts, _edges = np.histogramdd([arr_mass_branch], bins=bins)
        d = ROOT.RooDataHist.from_numpy(counts, [mass], bins=bins)


        mean  = ROOT.RooRealVar ( 'mean', '', 1020,  1016 ,1023)
        sigma = ROOT.RooRealVar ( 'sigma', '', 2,0.1,10)
        if mode =="sys":
           
            sig1   = ROOT.RooGaussian('sig1', 'sig1', mass, mean, sigma)
        
            sigma2 = ROOT.RooRealVar("sigma2","", 1,0.1,6)
            sig2 = ROOT.RooGaussian("sig2", "", mass, mean, sigma2)
            w12 = ROOT.RooRealVar("w12", "", 0.001, 0., 1.)
            sig = ROOT.RooAddPdf("sig", "", ROOT.RooArgList(sig1, sig2), ROOT.RooArgList(w12) )
            coeff = ROOT.RooRealVar("coeff", "coeff", 0, -1, 1)
            bkg = ROOT.RooChebychev("bkg", "bkg", mass, ROOT.RooArgList( coeff) ) 

        elif mode =="ana":

            resol = ROOT.RooRealVar('resol','resol',1,0.1, 6)
            sig   = ROOT.RooVoigtian('sig', 'sig', mass, mean, sigma, resol) 

            coeff = ROOT.RooRealVar("coeff", "coeff", 0, -1, 1)
            bkg = ROOT.RooChebychev("bkg", "bkg", mass, ROOT.RooArgList( coeff) )
            
            a1 = ROOT.RooRealVar('a1','a1',-1.,1.)
            a2 = ROOT.RooRealVar('a2','a2',0.,1.)
            bkg = Bkg('Bkg','Bkg',mass,a1,a2)


        numtot=d.sumEntries()

        N_sig=ROOT.RooRealVar("N_sig" ,"",0.8*numtot,0,numtot*1.3)
        N_bkg=ROOT.RooRealVar("N_bkg","",0.2*numtot,0,numtot*1.3)

        SIG=ROOT.RooExtendPdf("SIG","",sig,N_sig)
        BKG=ROOT.RooExtendPdf("BKG","",bkg,N_bkg)
        model=ROOT.RooAddPdf("model","sig+bkg",ROOT.RooArgSet(SIG,BKG))
        
    

        r = model.fitTo(d,ROOT.RooFit.Save(True),Minimizer("Minuit2"),ROOT.RooFit.PrintLevel(-1))
        frame = mass.frame (ROOT.RooFit.Name('frame'),ROOT.RooFit.Title(' '))
        d.plotOn ( frame , ROOT.RooFit.Name('data'))#,ROOT.RooFit.Binning(60))
        model.plotOn ( frame , ROOT.RooFit.Components ('sig'), ROOT.RooFit.LineColor( ROOT.kRed ) )
        model.plotOn ( frame , ROOT.RooFit.Components ('bkg'), ROOT.RooFit.LineColor( ROOT.kOrange ) )
        model.plotOn ( frame, ROOT.RooFit.Name('PDF_mass') )
        chi2 = frame.chiSquare()
        r.Print()
        can=ROOT.TCanvas("can","can")
        can.cd()
        frame.Draw()
        can.Update()
        ptext = TPaveText(0.65,0.65,0.85,0.85,"NDC")
        ptext.SetFillStyle(4000)
        ptext.SetBorderSize(0)
        ptext.AddText("{} < #eta <{}".format(eta_bins[i],eta_bins[i+1]))
        ptext.AddText("{} < p_{{T}} < {}".format(pt_bins[j],pt_bins[j+1]))
        ptext.AddText("N_{{sig}} ={:1.3f} #pm {:1.3f} ".format(N_sig.getValV(),N_sig.getError()))
        ptext.AddText("#chi^{{2}}/ndf={:1.3f} ".format(chi2))

        ptext.SetFillStyle(0)
        ptext.Draw()        
        can.SaveAs("/home3/sara.sellam/RpPb_identified_hadrons_project/efficiencies/sel_effi/fit_mass/Plots/can_mc_"+beam+"_before_"+str(eta_bins[i])+"_"+str(eta_bins[i+1])+"_"+str(pt_bins[j])+"_"+str(pt_bins[j+1])+"_"+mode+"_"+selection+".pdf")
        sig_dict[str(eta_bins[i])+"_"+str(eta_bins[i+1])][str(pt_bins[j])+"_"+str(pt_bins[j+1])]["before"]["value"]=N_sig.getValV()
        sig_dict[str(eta_bins[i])+"_"+str(eta_bins[i+1])][str(pt_bins[j])+"_"+str(pt_bins[j+1])]["before"]["error"]=N_sig.getError()
        
for i in range(len(eta_bins)-1):
    for j in range(len(pt_bins)-1):
        sig_dict[str(eta_bins[i])+"_"+str(eta_bins[i+1])][str(pt_bins[j])+"_"+str(pt_bins[j+1])]["after"]={}


        df_sel=df.Filter("Kmi_PT>"+str(pt_bins[j])+"&&Kmi_PT<"+str(pt_bins[j+1])+"&&Kmi_ETA>"+str(eta_bins[i])+"&&Kmi_ETA<"+str(eta_bins[i+1]))
        df_sel=df_sel.Filter(sel_cuts)
        
        mass = ROOT.RooRealVar( "phi_MM","#phi(1020)[MeV/c^{2}]",1000,1040) 
        
        _mass_branch= df_sel.AsNumpy(["phi_MM"])
        arr_mass_branch= _mass_branch["phi_MM"]
        bins = [np.linspace( 1000,1040, 100)]
        counts, _edges = np.histogramdd([arr_mass_branch], bins=bins)
        d = ROOT.RooDataHist.from_numpy(counts, [mass], bins=bins)


        mean  = ROOT.RooRealVar ( 'mean', '', 1020,  1016 ,1023)
        sigma = ROOT.RooRealVar ( 'sigma', '', 2,0.1,10)
        if mode =="sys":
           
            sig1   = ROOT.RooGaussian('sig1', 'sig1', mass, mean, sigma)
        
            sigma2 = ROOT.RooRealVar("sigma2","", 1,0.1,6)
            sig2 = ROOT.RooGaussian("sig2", "", mass, mean, sigma2)
            w12 = ROOT.RooRealVar("w12", "", 0.001, 0., 1.)
            sig = ROOT.RooAddPdf("sig", "", ROOT.RooArgList(sig1, sig2), ROOT.RooArgList(w12) )
            coeff = ROOT.RooRealVar("coeff", "coeff", 0, -1, 1)
            bkg = ROOT.RooChebychev("bkg", "bkg", mass, ROOT.RooArgList( coeff) ) 

        elif mode =="ana":

            resol = ROOT.RooRealVar('resol','resol',1,0.1, 6)
            sig   = ROOT.RooVoigtian('sig', 'sig', mass, mean, sigma, resol) 

            coeff = ROOT.RooRealVar("coeff", "coeff", 0, -1, 1)
            bkg = ROOT.RooChebychev("bkg", "bkg", mass, ROOT.RooArgList( coeff) )
            
            a1 = ROOT.RooRealVar('a1','a1',-1.,1.)
            a2 = ROOT.RooRealVar('a2','a2',0.,1.)
            bkg = Bkg('Bkg','Bkg',mass,a1,a2)


        numtot=d.sumEntries()

        N_sig=ROOT.RooRealVar("N_sig" ,"",0.8*numtot,0,numtot*1.3)
        N_bkg=ROOT.RooRealVar("N_bkg","",0.2*numtot,0,numtot*1.3)

        SIG=ROOT.RooExtendPdf("SIG","",sig,N_sig)
        BKG=ROOT.RooExtendPdf("BKG","",bkg,N_bkg)
        model=ROOT.RooAddPdf("model","sig+bkg",ROOT.RooArgSet(SIG,BKG))

        r = model.fitTo(d,ROOT.RooFit.Save(True),Minimizer("Minuit2"),ROOT.RooFit.PrintLevel(-1))
      
        frame = mass.frame (ROOT.RooFit.Name('frame'),ROOT.RooFit.Title(' '))
        d.plotOn ( frame , ROOT.RooFit.Name('data'))#,ROOT.RooFit.Binning(60))
        model.plotOn ( frame , ROOT.RooFit.Components ('sig'), ROOT.RooFit.LineColor( ROOT.kRed ) )
        model.plotOn ( frame , ROOT.RooFit.Components ('bkg'), ROOT.RooFit.LineColor( ROOT.kOrange ) )
        model.plotOn ( frame, ROOT.RooFit.Name('PDF_mass') )
        chi2 = frame.chiSquare()

        can=ROOT.TCanvas("can","can")
        can.cd()
        frame.Draw()
        can.Update()
        ptext = TPaveText(0.65,0.65,0.85,0.85,"NDC")
        ptext.SetFillStyle(4000)
        ptext.SetBorderSize(0)
        ptext.AddText("{} < #eta <{}".format(eta_bins[i],eta_bins[i+1]))
        ptext.AddText("{} < p_{{T}} < {}".format(pt_bins[j],pt_bins[j+1]))
        ptext.AddText("N_{{sig}} ={:1.3f} #pm {:1.3f} ".format(N_sig.getValV(),N_sig.getError()))
        ptext.AddText("#chi^{{2}}/ndf={:1.3f} ".format(chi2))

        ptext.SetFillStyle(0)
        ptext.Draw()
        can.SaveAs("/home3/sara.sellam/RpPb_identified_hadrons_project/efficiencies/sel_effi/fit_mass/Plots/can_mc_"+beam+"_after_"+str(eta_bins[i])+"_"+str(eta_bins[i+1])+"_"+str(pt_bins[j])+"_"+str(pt_bins[j+1])+"_"+mode+"_"+selection+".pdf")

        sig_dict[str(eta_bins[i])+"_"+str(eta_bins[i+1])][str(pt_bins[j])+"_"+str(pt_bins[j+1])]["after"]["value"]=N_sig.getValV()
        sig_dict[str(eta_bins[i])+"_"+str(eta_bins[i+1])][str(pt_bins[j])+"_"+str(pt_bins[j+1])]["after"]["error"]=N_sig.getError()


from uncertainties import ufloat
import json

for i in range(len(eta_bins)-1):
    for j in range(len(pt_bins)-1):
        sig_dict[str(eta_bins[i])+"_"+str(eta_bins[i+1])][str(pt_bins[j])+"_"+str(pt_bins[j+1])]["sel_effi_mc"]={}
        value_num=ufloat(sig_dict[str(eta_bins[i])+"_"+str(eta_bins[i+1])][str(pt_bins[j])+"_"+str(pt_bins[j+1])]["after"]["value"],sig_dict[str(eta_bins[i])+"_"+str(eta_bins[i+1])][str(pt_bins[j])+"_"+str(pt_bins[j+1])]["after"]["error"])
        value_denom=ufloat(sig_dict[str(eta_bins[i])+"_"+str(eta_bins[i+1])][str(pt_bins[j])+"_"+str(pt_bins[j+1])]["before"]["value"],sig_dict[str(eta_bins[i])+"_"+str(eta_bins[i+1])][str(pt_bins[j])+"_"+str(pt_bins[j+1])]["before"]["error"])
        ratio=value_num/value_denom
        sig_dict[str(eta_bins[i])+"_"+str(eta_bins[i+1])][str(pt_bins[j])+"_"+str(pt_bins[j+1])]["sel_effi_mc"]["value"]=ratio.nominal_value
        sig_dict[str(eta_bins[i])+"_"+str(eta_bins[i+1])][str(pt_bins[j])+"_"+str(pt_bins[j+1])]["sel_effi_mc"]["error"]=ratio.std_dev

# Define the output file path
output_file_path = "/home3/sara.sellam/RpPb_identified_hadrons_project/efficiencies/sel_effi/sel_effi_mc_"+beam+"_"+mode+"_"+selection+".json"

# Write the dictionary to the JSON file
with open(output_file_path, 'w') as json_file:
    json.dump(sig_dict, json_file, indent=4)