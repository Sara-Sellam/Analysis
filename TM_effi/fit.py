
import ROOT
from ROOT import RDF 
#Script containing te fit model

import ROOT
import numpy as np
import uproot as up
import pandas as pd
#import uproot3 as uproot
from Ntuple_path.call_calibration_ntuple  import _part
#import root_numpy as rnp
import sys
sys.path.append("/home3/sara.sellam/RpPb_identified_hadrons_project")
sys.path.append("/home3/sara.sellam/RpPb_identified_hadrons_project/efficiencies/pid_effi/p-ion")
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)

Minimizer = ROOT.RooFit.Minimizer
nottokill = list() #python and RooFit don't like themselves in terms of ownership
def prepareAndDrawPads (mypad1,mypad2):
    mypad1.SetBorderMode(0)
    mypad1.SetBottomMargin( 0.1 )
    mypad2.SetBottomMargin( 0.30  )
    mypad2.SetTopMargin( 0.1 )
    mypad2.SetBorderMode(0)
    mypad1.Draw()
    mypad2.Draw()
def MakeResidual(mypdf,mydata,frame,t):
   PlusTwosigma = ROOT.TLine(t.getMin() , 2 ,t.getMax() , 2 )
   zero_line    = ROOT.TLine(t.getMin() , 0 ,t.getMax() , 0 )
   MinusTwosigma= ROOT.TLine(t.getMin() ,-2 ,t.getMax() ,-2 )
   zero_line.SetLineColor(ROOT.kBlack)
   PlusTwosigma.SetLineColor(ROOT.kRed)
   MinusTwosigma.SetLineColor(ROOT.kRed)
   hpull= frame.pullHist( mydata ,mypdf, False)
   print("frame.pullHist( mydata ,mypdf, False)",hpull)
   hpull.SetMarkerSize(0.8)
   frame_pull  = t.frame()
   frame_pull.addObject(PlusTwosigma)
   frame_pull.addObject(MinusTwosigma)
   frame_pull.addObject(zero_line)
   frame_pull.addObject(hpull,"P" )
   frame_pull.SetAxisRange(-11 , 11,  "Y")
   frame_pull.GetXaxis().SetTitle("")
   frame_pull.GetYaxis().SetTitle("")
   frame_pull.GetXaxis().SetLabelSize(0.1)
   frame_pull.GetYaxis().SetLabelSize(0.1)
   frame_pull.GetYaxis().SetNdivisions(504)
   frame_pull.SetTitle("")
   return frame_pull, hpull


def mass_modelling( df, params, normOnly = True, binned=True,weights=False,fit_results_name="",beam="",mom="",outFile_name="",i="",j=""):


    mass = ROOT.RooRealVar( _part[beam][mom]["mom"] + "_MM", _part[beam][mom]["mom"] + ' candidates mass [MeV/c^{2}]',
                            params["massMin"] - params["shift"], params["massMax"] - params["shift"] ) 
 

    weight=ROOT.RooRealVar("w","",0,1)
    if isinstance(df,dict):
        df=RDF.FromNumpy(df)
    if weights:
        data_dict=df.AsNumpy([_part[beam][mom]["mom"]+"_MM","w"])
        d = ROOT.RooDataSet.from_numpy(data_dict,[mass,weight])
    else:
        data_dict=df.AsNumpy([_part[beam][mom]["mom"]+"_MM"])
        d = ROOT.RooDataSet.from_numpy(data_dict,[mass])
    
    if binned:
        if weights:
            _mass_branch= df.AsNumpy([_part[beam][mom]["mom"]+"_MM","w"])
            arr_mass_branch= _mass_branch[_part[beam][mom]["mom"] + "_MM"]
            bins = [np.linspace( params["massMin"] - params["shift"],params["massMax"] - params["shift"], 100)]
            counts, _ = np.histogramdd([arr_mass_branch], bins=bins)
            d = ROOT.RooDataHist.from_numpy(counts, [mass,weight], bins=bins)
        else:
            _mass_branch= df.AsNumpy([_part[beam][mom]["mom"]+"_MM"])
            arr_mass_branch= _mass_branch[_part[beam][mom]["mom"] + "_MM"]
            bins = [np.linspace( params["massMin"] - params["shift"],params["massMax"] - params["shift"], 100)]
            counts, _ = np.histogramdd([arr_mass_branch], bins=bins)
            d = ROOT.RooDataHist.from_numpy(counts, [mass], bins=bins)


    pars = list()
    params={"mass":1019.5,
            "massMin":1010,
            "massMax":1028,
            
            "sigma":2,
            "sigmaMin":0,
            "sigmaMax":10,
            
            "sigma2":1.5,
            "sigma2Min":0,
            "sigma2Max":6}

    #Signal: Gaussian/Gaussian composition/Voigtian
    mean  = ROOT.RooRealVar ( 'mean', 'mean1', params['mass'], params["massMin"]  ,params["massMax"] )
    sigma = ROOT.RooRealVar ( 'sigma', 'sigma1', params['sigma'], params['sigmaMin'], params['sigmaMax'])

    if params["sig"] == "VG":
        resol = ROOT.RooRealVar('resol','resol',params['sigma2'], params['sigma2Min'], params['sigma2Max'])
        sig   = ROOT.RooVoigtian('sig', 'sig', mass, mean, sigma, resol) 

        pars += [resol ]        
        nottokill.extend( [resol])
    """
    elif params["sig"] == "1G":
        sig   = ROOT.RooGaussian('sig', 'sig', mass, mean, sigma)

    elif params["sig"] == "2G":
        sig1   = ROOT.RooGaussian('sig1', 'sig1', mass, mean, sigma)
        sigma2 = ROOT.RooRealVar("sigma2","", params["sigma2"], params['sigma2Min'], params['sigma2Max'])
        sig2 = ROOT.RooGaussian("sig2", "", mass, mean, sigma2)
        w12 = ROOT.RooRealVar("w12", "", params["w12"], 0., 1.)
        sig = ROOT.RooAddPdf("sig", "", ROOT.RooArgList(sig1, sig2), ROOT.RooArgList(w12) )

        pars += [sigma2, w12 ]        
        nottokill.extend( [sig1, sig2])

    elif params["sig"] == "3G":
        sig1   = ROOT.RooGaussian('sig1', 'sig1', mass, mean, sigma)
        sigma2 = ROOT.RooRealVar("sigma2","", params["sigma2"], params['sigma2Min'], params['sigma2Max'])
        sig2 = ROOT.RooGaussian("sig2", "", mass, mean, sigma2)
        sigma3 = ROOT.RooRealVar("sigma3","", params["sigma3"], params['sigma3Min'], params['sigma3Max'])
        sig3 = ROOT.RooGaussian("sig3", "", mass, mean, sigma3)
        w12 = ROOT.RooRealVar("w12", "", params["w12"],  0., 1.)
        w23 = ROOT.RooRealVar("w23", "", params["w23"], 0., 1.)
        sig = ROOT.RooAddPdf("sig", "", ROOT.RooArgList(sig1, sig2, sig3), ROOT.RooArgList(w12, w23) )
        pars += [sigma2, sigma3, w12, w23]        
        nottokill.extend( [sig1, sig2, sig3])
    """


    """
    if params["bkg"] == "pol":
        coeff = ROOT.RooRealVar("coeff", "coeff",   0, -1, 1)
        coeff2 = ROOT.RooRealVar("coeff2", "coeff2",  0, -1, 1)
        bkg = ROOT.RooGenericPdf("bkg", "bkg", "coeff2 + coeff * {}".format(mass.GetName()), 
                                 ROOT.RooArgList(mass, coeff, coeff2) ) 
        pars += [coeff, coeff2]

    elif params["bkg"] == "parab":
        coeff = ROOT.RooRealVar("coeff", "coeff", params["coeff"], 50., 5000)
        coeff2 = ROOT.RooRealVar("coeff2", "coeff2", params["coeff2"], 50., 5000.) 
        bkg = ROOT.RooGenericPdf("bkg", "bkg", "-({} - coeff)*({} - coeff2)".format(mass.GetName(), mass.GetName()), 
                                 ROOT.RooArgList( mass, coeff, coeff2))
        pars += [coeff, coeff2]
    """
    if params["bkg"] == "Cheby1":
        coeff = ROOT.RooRealVar("coeff", "coeff", 0.15 , -1, 1)
        bkg = ROOT.RooChebychev("bkg", "bkg", mass, ROOT.RooArgList( coeff) ) 
        pars += [coeff]
    """
    elif params["bkg"] == "Cheby2":
        coeff = ROOT.RooRealVar("coeff", "coeff", params["coeff"], -1, 1)
        coeff2 = ROOT.RooRealVar("coeff2", "coeff2", params["coeff2"], -1, 1) 
        coeff3 = ROOT.RooRealVar("coeff3", "coeff3", params["coeff3"], -1, 1) 
        bkg = ROOT.RooChebychev("bkg", "bkg", mass, ROOT.RooArgList( coeff,coeff2 ) ) 
        pars += [coeff, coeff2]
    """
    #Signal + Background Addition
    f_sig = ROOT.RooRealVar( 'f_sig', 'f_sig', params["fsig"], 0., 1. )
    numtot=d.sumEntries()
    
    N_sig=ROOT.RooRealVar("N_sig" ,"",0.8*numtot,0,numtot*1.3)
    N_bkg=ROOT.RooRealVar("N_bkg","",0.2*numtot,0,numtot*1.3)

    SIG=ROOT.RooExtendPdf("SIG","",sig,N_sig)
    BKG=ROOT.RooExtendPdf("BKG","",bkg,N_bkg)
    model=ROOT.RooAddPdf("model","sig+bkg",ROOT.RooArgSet(SIG,BKG))

    #model=ROOT.RooAddPdf('model', 'model',ROOT.RooArgList( sig, bkg),ROOT.RooArgList(N_sig,N_bkg))
    r = model.fitTo(d,ROOT.RooFit.Save(True),Minimizer("Minuit2"),ROOT.RooFit.PrintLevel(-1))
    frame = mass.frame (ROOT.RooFit.Name('frame'),ROOT.RooFit.Title(' '))
    d.plotOn ( frame , ROOT.RooFit.Name('data'))#,ROOT.RooFit.Binning(60))
    model.plotOn ( frame , ROOT.RooFit.Components ('sig'), ROOT.RooFit.LineColor( ROOT.kRed ) )
    model.plotOn ( frame , ROOT.RooFit.Components ('bkg'), ROOT.RooFit.LineColor( ROOT.kOrange ) )
    model.plotOn ( frame, ROOT.RooFit.Name('PDF_mass') )
    
    f=ROOT.TFile(outFile_name+"_"+str(i)+"_"+str(j)+".root","recreate")
    d.Write("data")
    frame.Write("frame")
    model.Write("PDF_mass")
    f.Close()
    r.Print()

    # Calculate the chi-square value
    chi2 = frame.chiSquare()

    # Calculate the degrees of freedom
    num_data_points = d.numEntries()
    num_fitted_parameters = model.getParameters(d).getSize()
    NDF = num_data_points - num_fitted_parameters
    status=r.status()
    print("Chi-square value:", chi2)
    print("Degrees of freedom (NDF):", NDF)
    print("Status",r.status())
    out={
        "N_sig"      : N_sig.getValV(), 
        "N_sig_error":N_sig.getError(), 
        "N_bkg"      : N_bkg.getValV(),
        "chi2"     :chi2,
        "status":status,
        }

    return out    

