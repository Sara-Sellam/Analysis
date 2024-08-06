import ROOT
from ROOT import *
import numpy as np 
from ROOT import RDataFrame
Minimizer = ROOT.RooFit.Minimizer
gStyle.SetPaintTextFormat("1.3f")
gStyle.SetPalette(1)

beam="pPb"

ch=TChain()
ch.Add("/scratch43/ssellam/calibration/MC/Ks_pipi/"+beam+"/*"+beam+"*.root/KsTuple/DecayTree")
df=RDataFrame(ch)

df=df.Filter("abs(Ks_TRUEID)==310&&abs(pi_pl_TRUEID)==211&&abs(pi_mi_TRUEID)==211&&abs(pi_pl_MC_MOTHER_ID)==310&&abs(pi_mi_MC_MOTHER_ID)==310")

pt_bins=[500,1249, 4619]
eta_bins=[2,3,4,4.8]
sig_dict={}
mode="sys" #sys #ana

cuts={"pseudoIP":{"pp":0.368,
                 "pPb":0.348,
                  "Pbp":0.348},
        "GhostP":{"pp":0.078,
                 "pPb":0.103,
                  "Pbp":0.109}}

for i in range(len(eta_bins)-1):
    sig_dict[str(eta_bins[i])+"_"+str(eta_bins[i+1])]={}
    for j in range(len(pt_bins)-1):
        sig_dict[str(eta_bins[i])+"_"+str(eta_bins[i+1])][str(pt_bins[j])+"_"+str(pt_bins[j+1])]={}

for i in range(len(eta_bins)-1):
    #sig_dict[str(eta_bins[i])+"_"+str(eta_bins[i+1])]={}
    for j in range(len(pt_bins)-1):
        sig_dict[str(eta_bins[i])+"_"+str(eta_bins[i+1])][str(pt_bins[j])+"_"+str(pt_bins[j+1])]["before"]={}


        df_sel=df.Filter("pi_pl_PT>"+str(pt_bins[j])+"&&pi_pl_PT<"+str(pt_bins[j+1])+"&&pi_pl_ETA>"+str(eta_bins[i])+"&&pi_pl_ETA<"+str(eta_bins[i+1]))
    
        
        mass = ROOT.RooRealVar( "Ks_MM","",475,520) 
        
        _mass_branch= df_sel.AsNumpy(["Ks_MM"])
        arr_mass_branch= _mass_branch["Ks_MM"]
        bins = [np.linspace( 475,520, 100)]
        counts, _edges = np.histogramdd([arr_mass_branch], bins=bins)
        d = ROOT.RooDataHist.from_numpy(counts, [mass], bins=bins)


        mean  = ROOT.RooRealVar ( 'mean', '', 497.6,  496 ,500)
        sigma = ROOT.RooRealVar ( 'sigma', '', 2,0.1,8)
        if mode =="ana":
           
            sig1   = ROOT.RooGaussian('sig1', 'sig1', mass, mean, sigma)
        
            sigma2 = ROOT.RooRealVar("sigma2","", 1,0.1,5)
            sig2 = ROOT.RooGaussian("sig2", "", mass, mean, sigma2)
            w12 = ROOT.RooRealVar("w12", "", 0.001, 0., 1.)
            sig = ROOT.RooAddPdf("sig", "", ROOT.RooArgList(sig1, sig2), ROOT.RooArgList(w12) )


            coeff = ROOT.RooRealVar("coeff", "coeff", 0, -1, 1)
            bkg = ROOT.RooChebychev("bkg", "bkg", mass, ROOT.RooArgList( coeff) ) 
        
        if mode =="sys":

            resol = ROOT.RooRealVar('resol','resol',1,0.1, 3)
            sig   = ROOT.RooVoigtian('sig', 'sig', mass, mean, sigma, resol) 

            coeff = ROOT.RooRealVar("coeff", "coeff", 0, -1, 1)
            bkg = ROOT.RooChebychev("bkg", "bkg", mass, ROOT.RooArgList( coeff) )

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
        can.SaveAs("Plots/can_mc_before_"+str(eta_bins[i])+"_"+str(eta_bins[i+1])+"_"+str(pt_bins[j])+"_"+str(pt_bins[j+1])+"_"+mode+".pdf")
        sig_dict[str(eta_bins[i])+"_"+str(eta_bins[i+1])][str(pt_bins[j])+"_"+str(pt_bins[j+1])]["before"]["value"]=N_sig.getValV()
        sig_dict[str(eta_bins[i])+"_"+str(eta_bins[i+1])][str(pt_bins[j])+"_"+str(pt_bins[j+1])]["before"]["error"]=N_sig.getError()

for i in range(len(eta_bins)-1):
    for j in range(len(pt_bins)-1):
        sig_dict[str(eta_bins[i])+"_"+str(eta_bins[i+1])][str(pt_bins[j])+"_"+str(pt_bins[j+1])]["after"]={}


        df_sel=df.Filter("pi_pl_PT>"+str(pt_bins[j])+"&&pi_pl_PT<"+str(pt_bins[j+1])+"&&pi_pl_ETA>"+str(eta_bins[i])+"&&pi_pl_ETA<"+str(eta_bins[i+1]))
        sel_cuts="pi_pl_Reconstructed==1&&pi_pl_TRACK_Type==3&&pi_pl_P>7000&&BCType==3&&pi_pl_TRACK_GhostProb<"+str(cuts["GhostP"][beam])+"&&pi_pl_IPCHI2_ORIVX<"+str(cuts["pseudoIP"][beam])
        df_sel=df_sel.Filter(sel_cuts)
        
        mass = ROOT.RooRealVar( "Ks_MM","",475,520) 
        
        _mass_branch= df_sel.AsNumpy(["Ks_MM"])
        arr_mass_branch= _mass_branch["Ks_MM"]
        bins = [np.linspace( 475,520, 100)]
        counts, _edges = np.histogramdd([arr_mass_branch], bins=bins)
        d = ROOT.RooDataHist.from_numpy(counts, [mass], bins=bins)


        mean  = ROOT.RooRealVar ( 'mean', '', 497.6,  496 ,500)
        sigma = ROOT.RooRealVar ( 'sigma', '', 2,0.1,8)
        if mode =="ana":
           
            sig1   = ROOT.RooGaussian('sig1', 'sig1', mass, mean, sigma)
        
            sigma2 = ROOT.RooRealVar("sigma2","", 1,0.1,5)
            sig2 = ROOT.RooGaussian("sig2", "", mass, mean, sigma2)
            w12 = ROOT.RooRealVar("w12", "", 0.001, 0., 1.)
            sig = ROOT.RooAddPdf("sig", "", ROOT.RooArgList(sig1, sig2), ROOT.RooArgList(w12) )


            coeff = ROOT.RooRealVar("coeff", "coeff", 0, -1, 1)
            bkg = ROOT.RooChebychev("bkg", "bkg", mass, ROOT.RooArgList( coeff) ) 
        
        if mode =="sys":

            resol = ROOT.RooRealVar('resol','resol',1,0.1, 3)
            sig   = ROOT.RooVoigtian('sig', 'sig', mass, mean, sigma, resol) 

            coeff = ROOT.RooRealVar("coeff", "coeff", 0, -1, 1)
            bkg = ROOT.RooChebychev("bkg", "bkg", mass, ROOT.RooArgList( coeff) ) 

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
        can.SaveAs("Plots/can_mc_after_"+str(eta_bins[i])+"_"+str(eta_bins[i+1])+"_"+str(pt_bins[j])+"_"+str(pt_bins[j+1])+"_"+mode+".pdf")
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
output_file_path = "sel_effi_mc_"+beam+"_"+mode+".json"

# Write the dictionary to the JSON file
with open(output_file_path, 'w') as json_file:
    json.dump(sig_dict, json_file, indent=4)