#Script containing te fit model

from ROOT import RDF
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

_part_latex={"Ks":"m(\pi^{+}\pi^{-})",
             "Lambda":"m(p \pi)",
             "phi":"m(K^{+} K^{-})"}

def mass_modelling( df, params, normOnly = True, binned=True,weights=False,fit_results_name="",beam="",mom="",outFile_name="",i="",j="",sig_bkg=["sig","bkg"]):


    mass = ROOT.RooRealVar( _part[beam][mom]["mom_mass"], _part_latex[mom]+ '[MeV/c^{2}]',
                            params["massMin"], params["massMax"]  ) 
                            
    print("I am here 1")
    if weights in ["0_pi","0_K","0_p","1_pi","1_K","1_p","2_pi","2_K","2_p","3_pi","3_K","3_p"]:
        print("I am here 2")
        weight=ROOT.RooRealVar("w","",0,1)
    else:
        weight=ROOT.RooRealVar("probe_before_sw_w"+weights,"",0,1)
    if isinstance(df,dict):
        df=RDF.FromNumpy(df)
    
    
    if binned:
        if weights!="0":
            if weights in ["0_pi","0_K","0_p","1_pi","1_K","1_p","2_pi","2_K","2_p","3_pi","3_K","3_p"]:
                print("I am here 3")
                _mass_branch= df.AsNumpy([_part[beam][mom]["mom_mass"],"w"])
                arr_weights=_mass_branch["w"]

            else:
                _mass_branch= df.AsNumpy([_part[beam][mom]["mom_mass"],"probe_before_sw_w"+weights])
                arr_weights=_mass_branch["probe_before_sw_w"+weights]

            arr_mass_branch= _mass_branch[_part[beam][mom]["mom_mass"]]
            bins = [np.linspace( params["massMin"] - params["shift"],params["massMax"] - params["shift"], 100)]
            counts, _ = np.histogramdd([arr_mass_branch], bins=bins ,weights=arr_weights)
            weights_squared_sum = np.sqrt(np.sum(arr_weights**2)) * np.ones_like(counts)
            weights_squared_sum = np.sqrt(np.histogramdd([arr_mass_branch], bins=bins ,weights=arr_weights**2)[0])

            d = ROOT.RooDataHist.from_numpy(counts, [mass], weights_squared_sum=weights_squared_sum,bins=bins)
            #d = ROOT.RooDataHist.from_numpy(counts, [mass,weight], bins=bins)
        else:
            _mass_branch= df.AsNumpy([_part[beam][mom]["mom_mass"]])
            arr_mass_branch= _mass_branch[_part[beam][mom]["mom_mass"]]
            bins = [np.linspace( params["massMin"] - params["shift"],params["massMax"] - params["shift"], 100)]
            counts, _edges = np.histogramdd([arr_mass_branch], bins=bins)
            d = ROOT.RooDataHist.from_numpy(counts, [mass], bins=bins)
            bin_widths = np.diff(_edges[0])
            #errors_per_unit = np.sqrt(counts) / bin_widths
    else:
        counts=0
        if weights!="0":
            if weights in ["0_pi","0_K","0_p","1_pi","1_K","1_p","2_pi","2_K","2_p"]:
                data_dict=df.AsNumpy([_part[beam][mom]["mom_mass"],"w"])

            else:
                data_dict=df.AsNumpy([_part[beam][mom]["mom_mass"],"probe_before_sw_w"+weights])
            
            d = ROOT.RooDataSet.from_numpy(data_dict,[mass,weight])
        else:
            data_dict=df.AsNumpy([_part[beam][mom]["mom_mass"]])
            d = ROOT.RooDataSet.from_numpy(data_dict,[mass])

    #d = ROOT.RooDataSet( "ds", 'Mass', ROOT.RooArgSet(mass),ROOT.RooFit.Import(tree))


    
    pars = list()

    #Signal: Gaussian/Gaussian composition/Voigtian
    mean  = ROOT.RooRealVar ( 'mean', 'mean1', params['mass'], params["massMin"]  ,params["massMax"] )
    #mean  = ROOT.RooRealVar ( 'mean', 'mean1', params['mass'])#, params["massMin"]  ,params["massMax"] )
    sigma = ROOT.RooRealVar ( 'sigma', 'sigma1', params['sigma'], params['sigmaMin'], params['sigmaMax'])
    
    sig_mode=sig_bkg[0]
    bkg_mode=sig_bkg[1]
    if params[sig_mode] == "VG":

        resol = ROOT.RooRealVar('resol','resol',params['sigma2'], params['sigma2Min'], params['sigma2Max'])
        if mom =="phi":
            sigma = ROOT.RooRealVar ('sigma', 'sigma1',4.26)
        sig   = ROOT.RooVoigtian('sig', 'sig', mass, mean, sigma, resol) 

        pars += [resol ]        
        nottokill.extend( [resol])

    elif params[sig_mode] == "1G":
        sig   = ROOT.RooGaussian('sig', 'sig', mass, mean, sigma)

    elif params[sig_mode] == "2G":
        sig1   = ROOT.RooGaussian('sig1', 'sig1', mass, mean, sigma)
        sigma2 = ROOT.RooRealVar("sigma2","", params["sigma2"], params['sigma2Min'], params['sigma2Max'])
        sig2 = ROOT.RooGaussian("sig2", "", mass, mean, sigma2)
        w12 = ROOT.RooRealVar("w12", "", params["w12"], 0., 1.)
        sig = ROOT.RooAddPdf("sig", "", ROOT.RooArgList(sig1, sig2), ROOT.RooArgList(w12) )
        if mom=="Ks":
            sigma.setVal(2.45)
            sigma.setConstant()
            sigma2.setVal(5)
            sigma2.setConstant()
            
        pars += [sigma2, w12 ]        
        nottokill.extend( [sig1, sig2])

    elif params[sig_mode] == "3G":
        if mom =="Ks":
            sigma = ROOT.RooRealVar ( 'sigma', 'sigma1',3.6)
        sig1   = ROOT.RooGaussian('sig1', 'sig1', mass, mean, sigma)
        sigma2 = ROOT.RooRealVar("sigma2","", params["sigma2"], params['sigma2Min'], params['sigma2Max'])
        if mom =="Ks":
            sigma2 = ROOT.RooRealVar ( 'sigma2', '',6.6)
        sig2 = ROOT.RooGaussian("sig2", "", mass, mean, sigma2)
        sigma3 = ROOT.RooRealVar("sigma3","", params["sigma3"], params['sigma3Min'], params['sigma3Max'])
        if mom =="Ks":
            sigma3 = ROOT.RooRealVar ( 'sigma3', '',1.97)
        sig3 = ROOT.RooGaussian("sig3", "", mass, mean, sigma3)
        w12 = ROOT.RooRealVar("w12", "", params["w12"],  0., 1.)
        w23 = ROOT.RooRealVar("w23", "", params["w23"], 0., 1.)
        sig = ROOT.RooAddPdf("sig", "", ROOT.RooArgList(sig1, sig2, sig3), ROOT.RooArgList(w12, w23) )
        pars += [sigma2, sigma3, w12, w23]        
        nottokill.extend( [sig1, sig2, sig3])
    elif params[sig_mode] == "CB_1G" :
        # Define the parameters for the two Crystal Ball functions
        #mean1 = ROOT.RooRealVar("mean1", "mean1", 1.115, 1.10, 1.13)
        sigma1 = ROOT.RooRealVar("sigma1", "sigma1", params["sigma"], params['sigmaMin'], params['sigmaMax'])

        #mean2 = ROOT.RooRealVar("mean2", "mean2", 1.117, 1.10, 1.13)
        sigma2 = ROOT.RooRealVar("sigma2", "sigma2", params["sigma2"], params['sigma2Min'], params['sigma2Max'])
        alpha2 = ROOT.RooRealVar("alpha2", "alpha2",  -1.0, -5.0, -0.1)
        n2 = ROOT.RooRealVar("n2", "n2",1.0, 0.1, 10.0)
        # Create the two Crystal Ball functions
        cb1 = ROOT.RooGaussian("sig1", "sig1", mass, mean, sigma1)
        cb2 = ROOT.RooCBShape("cb2", "cb2", mass, mean, sigma2, alpha2, n2)
        mean = {}
        sigma = {}
        gaus = {}
        w12 = ROOT.RooRealVar("w12", "", params["w12"],  0., 1.)
        sig = ROOT.RooAddPdf("sig", "sig", ROOT.RooArgList(cb1, cb2), ROOT.RooArgList(w12))
        nottokill.extend( [cb1, cb2])
        pars += [sigma1,sigma2,alpha2,n2]    

    elif params[sig_mode] == "2CB" :
        # Define the parameters for the two Crystal Ball functions
        #mean1 = ROOT.RooRealVar("mean1", "mean1", 1.115, 1.10, 1.13)
        sigma1 = ROOT.RooRealVar("sigma1", "sigma1", params["sigma"], params['sigmaMin'], params['sigmaMax'])
        alpha1 = ROOT.RooRealVar("alpha1", "alpha1", params["alpha1"],params['alpha1Min'], params['alpha1Max'])
        n1 = ROOT.RooRealVar("n1", "n1", 1.0, 0.1, 10.0)

        #mean2 = ROOT.RooRealVar("mean2", "mean2", 1.117, 1.10, 1.13)
        sigma2 = ROOT.RooRealVar("sigma2", "sigma2", params["sigma2"], params['sigma2Min'], params['sigma2Max'])
        alpha2 = ROOT.RooRealVar("alpha2", "alpha2", params["alpha2"],params['alpha2Min'], params['alpha2Max'])
        n2 = ROOT.RooRealVar("n2", "n2", 1., 1.0, 10)
        # Create the two Crystal Ball functions
        cb1 = ROOT.RooCBShape("cb1", "cb1", mass, mean, sigma1, alpha1, n1)
        cb2 = ROOT.RooCBShape("cb2", "cb2", mass, mean, sigma2, alpha2, n2)
        mean = {}
        sigma = {}
        gaus = {}
        frac1 = ROOT.RooRealVar("frac1", "frac1", 0.5, 0.0, 1.0)
        #sig = ROOT.RooAddPdf("sig", "sig", ROOT.RooArgList(cb1, cb2), ROOT.RooArgList(frac1))

       
        for iGaus in [str(i) for i in range(0, 3)]:
            m = ROOT.RooRealVar("mean" + iGaus, "mean" + iGaus, 1116.3, 1115, 1117)
            s = ROOT.RooRealVar("sigma" + iGaus, "sigma" + iGaus, 1.5, 0.5, 30)
            gaus[iGaus] =ROOT. RooGaussian("gaus" + iGaus, "gaus" + iGaus, mass, m, s)
            mean[iGaus] = m
            sigma[iGaus] = s
        cf1 = ROOT.RooRealVar("cf1", "cf1", 0.2, 0, 1.)
        cf2 = ROOT.RooRealVar("cf2", "cf2", 1., 0, 1.)
        signalCB = ROOT.RooAddPdf("signalL0", "signalGaus", cb1, cb2, cf1)
        sig = ROOT.RooAddPdf("sig", "", signalCB, gaus['0'], cf2)

        #signal =ROOT.RooAddPdf("signalL0", "signal", cb1, gaus['0'], cf1)
        #sig =ROOT.RooAddPdf("sig", "signal", cb1, gaus['1'], cf1)
        #sig = ROOT.RooAddPdf("sig", "", signal,gaus['1'], cf2)
        #cf2.setConstant()
        nottokill.extend( [cb1, cb2])
        pars += [sigma1,sigma2,alpha1,alpha2,n2,n1]
        
    elif params[sig_mode] == "1CB" :
        
        mean1 = ROOT.RooRealVar("mean1", "mean1", 1.115, 1.10, 1.13)
        sigma1 = ROOT.RooRealVar("sigma1", "sigma1", params["sigma"], params['sigmaMin'], params['sigmaMax'])
        alpha1 = ROOT.RooRealVar("alpha1", "alpha1", 1.0, 0.1, 5.0)
        n1 = ROOT.RooRealVar("n1", "n1", 1.0, 0.1, 20.0)

        # Define the fraction of events in each Crystal Ball function
        frac1 = ROOT.RooRealVar("frac1", "frac1", 0.5, 0.0, 1.0)
        sig = ROOT.RooCBShape("sig", "sig", mass, mean, sigma1, alpha1, n1)

        # Create the sum of the two Crystal Ball functions
        #sig = ROOT.RooAddPdf("sig", "sig", ROOT.RooArgList(cb1))
        #nottokill.extend( [cb1]
        
        pars += [sigma1,alpha1,n1]
    elif params[sig_mode] == "2CB_test" :
        
        mean1 = ROOT.RooRealVar("mean1", "mean1", 1.115, 1.10, 1.13)
        sigma1 = ROOT.RooRealVar("sigma1", "sigma1", params["sigma"], params['sigmaMin'], params['sigmaMax'])
        alpha1 = ROOT.RooRealVar("alpha1", "alpha1", 1.0, 0.1, 5.0)
        n1 = ROOT.RooRealVar("n1", "n1", 1.0, 0.1, 15.0)

        mean2 = ROOT.RooRealVar("mean2", "mean2", 1.117, 1.10, 1.13)
        sigma2 = ROOT.RooRealVar("sigma2", "sigma2",params["sigma"], params['sigmaMin'], params['sigmaMax'])
        alpha2 = ROOT.RooRealVar("alpha2", "alpha2", -1.0, -5.0, -0.1)
        n2 = ROOT.RooRealVar("n2", "n2", 1.0, 0.1, 5)

        # Create the two Crystal Ball functions
        cb1 = ROOT.RooCBShape("cb1", "cb1", mass, mean, sigma1, alpha1, n1)
        cb2 = ROOT.RooCBShape("cb2", "cb2", mass, mean, sigma2, alpha2, n2)

        # Define the fraction of events in each Crystal Ball function
        frac1 = ROOT.RooRealVar("frac1", "frac1", 0.5, 0.0, 1.0)

        # Create the sum of the two Crystal Ball functions
        sig = ROOT.RooAddPdf("sig", "sig", ROOT.RooArgList(cb1, cb2), ROOT.RooArgList(frac1))
        nottokill.extend( [cb1, cb2])
        pars += [sigma1,sigma2,alpha1,alpha2,n2,n1]
    elif params[sig_mode] == "RBW" :
        mean = ROOT.RooRealVar("mean", "mean", 1019.46, 1015, 1025)    
        width = ROOT.RooRealVar("width", "width", 4.26, 0.1, 10)

        # Define the RooGenericPdf object pointing to the BreitWigner function
        sig1 = ROOT.RooGenericPdf("sig1", "sig1", "BreitWigner(mass, m0, gamma)", {"mass": mass, "m0": mean, "gamma": width})
        # Define the RooGenericPdf object pointing to myPdfFunction
        #sig1 = RooGenericPdf("sig1", "My sig1", "BreitWigner({m,s}, {m0,gamma})", {"m": mass, "s": s, "m0": m0, "gamma": gamma})
        #sig1   = BreitWigner('sig1', 'sig1', mass, mean, sigma)


    #Background: Argus/polynomial function/Chebychev polynomial function
    if params[bkg_mode] == "Argus": 
        m0 = ROOT.RooRealVar("m0", "m0", params["Argus_end"], 60, 1100)
        c = ROOT.RooRealVar("c", "c", params["Argus_c"], -1, 1)
        p = ROOT.RooRealVar("p", "p", params["Argus_p"], -100., 40)
        bkg = ROOT.RooArgusBG ('bkg', 'bkg', mass, m0, c, p)    
        pars += [m0, c, p]

    elif params[bkg_mode] == "pol":
        coeff = ROOT.RooRealVar("coeff", "coeff",   0, -1, 1)
        coeff2 = ROOT.RooRealVar("coeff2", "coeff2",  0, -1, 1)
        bkg = ROOT.RooGenericPdf("bkg", "bkg", "coeff2 + coeff * {}".format(mass.GetName()), 
                                 ROOT.RooArgList(mass, coeff, coeff2) ) 
        pars += [coeff, coeff2]

    elif params[bkg_mode] == "parab":
        coeff = ROOT.RooRealVar("coeff", "coeff", params["coeff"], 50., 5000)
        coeff2 = ROOT.RooRealVar("coeff2", "coeff2", params["coeff2"], 50., 5000.) 
        bkg = ROOT.RooGenericPdf("bkg", "bkg", "-({} - coeff)*({} - coeff2)".format(mass.GetName(), mass.GetName()), 
                                 ROOT.RooArgList( mass, coeff, coeff2))
        pars += [coeff, coeff2]

    elif params[bkg_mode] == "Cheby1":
        coeff = ROOT.RooRealVar("coeff", "coeff", params["coeff"], -1, 1)
        coeff2 = ROOT.RooRealVar("coeff2", "coeff2", params["coeff2"], -1, 1) 
        coeff3 = ROOT.RooRealVar("coeff3", "coeff3", params["coeff3"], -1, 1) 
        bkg = ROOT.RooChebychev("bkg", "bkg", mass, ROOT.RooArgList( coeff) ) 
        pars += [coeff]
    
    elif params[bkg_mode] == "Cheby2":
        coeff = ROOT.RooRealVar("coeff", "coeff", params["coeff"], -1, 1)
        coeff2 = ROOT.RooRealVar("coeff2", "coeff2", params["coeff2"], -1, 1) 
        coeff3 = ROOT.RooRealVar("coeff3", "coeff3", params["coeff3"], -1, 1) 
        bkg = ROOT.RooChebychev("bkg", "bkg", mass, ROOT.RooArgList( coeff,coeff2 ) ) 
        pars += [coeff, coeff2]
    #Signal + Background Addition
    f_sig = ROOT.RooRealVar( 'f_sig', 'f_sig', params["fsig"], 0., 1. )
    numtot=d.sumEntries()
    #model = ROOT.RooAddPdf ( 'model', 'model', ROOT.RooArgList( sig, bkg), ROOT.RooArgList( f_sig )) 

    
    #N_sig = ROOT.RooFormulaVar("N_sig", "f_sig * {}".format(d.sumEntries()), ROOT.RooArgList(f_sig))
    #N_bkg = ROOT.RooFormulaVar("N_bkg", "(1. - f_sig ) * {}".format(d.sumEntries()), ROOT.RooArgList(f_sig))
    N_sig=ROOT.RooRealVar("N_sig" ,"",0.8*numtot,0,numtot*1.3)
    N_bkg=ROOT.RooRealVar("N_bkg","",0.2*numtot,0,numtot*1.3)

    SIG=ROOT.RooExtendPdf("SIG","",sig,N_sig)
    BKG=ROOT.RooExtendPdf("BKG","",bkg,N_bkg)
    model=ROOT.RooAddPdf("model","sig+bkg",ROOT.RooArgSet(SIG,BKG))

    if params[sig_mode] == "2G":   
        if mom=="Ks":
            if beam=="pPb":
                mean.setVal(497.82)
                sigma.setVal(2.54)
                sigma2.setVal(5.2)
                mean.setConstant()
                sigma.setConstant()
                sigma2.setConstant()
            elif beam=="Pbp":
                mean.setVal(497.82)
                sigma.setVal(2.54)
                sigma2.setVal(5.2)
                mean.setConstant()
                sigma.setConstant()
                sigma2.setConstant()
        elif mom =="Lambda":
            
            if beam=="ppPb":
                mean.setVal(1115.75)
                sigma.setVal(2.03)
                sigma2.setVal(0.82)
                mean.setConstant()
                sigma.setConstant()
                sigma2.setConstant()
            
                
            elif beam=="Pbp":
                mean.setVal(1115.76)
                sigma.setVal(2.19)
                sigma2.setVal(0.88)
                mean.setConstant()
                sigma.setConstant()
                sigma2.setConstant()




    #model=ROOT.RooAddPdf('model', 'model',ROOT.RooArgList( sig, bkg),ROOT.RooArgList(N_sig,N_bkg))
    
    r = model.fitTo(d,ROOT.RooFit.Save(True),Minimizer("Minuit2"),ROOT.RooFit.PrintLevel(-1))
    #,ROOT.RooFit.Parallelize(20))#,ROOT.RooFit.BatchMode("cpu"))#, ROOT.RooFit.Save( True ), ROOT.RooFit.Extended(False))
    #r.Print()
    frame = mass.frame (ROOT.RooFit.Name('frame'),ROOT.RooFit.Title(' '))
    d.plotOn ( frame , ROOT.RooFit.Name('data'))#,ROOT.RooFit.Binning(60))
    model.plotOn ( frame , ROOT.RooFit.Components ('sig'), ROOT.RooFit.LineColor( ROOT.kRed ) )
    model.plotOn ( frame , ROOT.RooFit.Components ('bkg'), ROOT.RooFit.LineColor( ROOT.kOrange ) )
    model.plotOn ( frame, ROOT.RooFit.Name('PDF_mass') )
    can=ROOT.TCanvas("can","can")
    can.cd()
    frame.Draw()
    f=ROOT.TFile(outFile_name+"_"+str(i)+"_"+str(j)+".root","recreate")
    d.Write("data")
    frame.Write("frame")
    model.Write("PDF_mass")
    can.Write("can")
    f.Close()


    #f=ROOT.TFile("/scratch43/ssellam/calibration/p-ion/sWeights/Plots/"+part+"_MM_fit.root","recreate")
    #f.Close()
    """
    
    pad_M= ROOT.TPad("pad_M","pad_M",0 ,0.25 ,1 ,1) 
    pad_pull=ROOT.TPad("pad_M_pull","pad_M_pull",0 ,0,1 ,0.25)
    prepareAndDrawPads(pad_M ,pad_pull )
    frame_Mass_pull,hpull = MakeResidual("PDF_mass","data",frame,mass)
    frame_Mass_pull.SetTitle("")
    pad_M.cd()
    frame.Draw()
    pad_pull.cd()
    frame_Mass_pull.Draw()
    can.Close()
    """ 
    #can.Print("/scratch43/ssellam/calibration/p-ion/sWeights/Plots/"+part+"_MM_fit.pdf")

    #for ipar in pars: ipar.setConstant() 
    #nottokill.extend ( pars ) 
    #nottokill.extend([sig, bkg])
    #r = model.fitTo( d, ROOT.RooFit.Save( True ), ROOT.RooFit.Extended(False) ) 
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
    print("N sig:", N_sig.getValV())
    print("Status",r.status())
    #ROOT.gPad.SetLogy()
    #ROOT.gPad.SetLogy(False)
    """
    data_sig_values, data_bkg_values = list(), list()
       
    #print("N_sig",N_sig.getVal(),"N_bkg",N_bkg.getVal(),"numtot:",d.sumEntries())
    tot_sig = 0
    tot_bkg = 0
    
    dict_MM=df.AsNumpy(columns=[_part[beam][mom]["mom_mass"]])
    df_MM=pd.DataFrame(dict_MM)
    
    for ientry in df_MM[ _part[beam][mom]["mom"]+ "_MM" ]:
        M = mass
        M.setVal( ientry )
        #print(M)
        isig = sig.getValV( ROOT.RooArgSet( M )) 
        ibkg = bkg.getValV( ROOT.RooArgSet( M ) )
        #print("isig",isig)
        data_sig_values.append( isig ) 
        data_bkg_values.append( ibkg ) 
        tot_sig += sig.getValV( ROOT.RooArgSet( M ))
        tot_bkg += bkg.getValV( ROOT.RooArgSet( M ) )
    
    """
    out={
        #"mass"       : mass, 
        #"Signal"        : sig,  
        #"Background"        : bkg,  
        #"model":    model,
        #"data_sig_values":data_sig_values,
        #"data_bkg_values":data_bkg_values,
        #"tot_sig":tot_sig,
        #"tot_bkg":tot_bkg,
    
        "N_sig"      : N_sig.getValV(), 
        "N_sig_error":N_sig.getError(), 
        "N_bkg"      : N_bkg.getValV(),
        "N_bkg_error":N_bkg.getError(),
        "chi2"     :chi2,
        
        "status":status,
        #"uncert": uncert,
        "counts":counts
        }
    


    # Calculate the Poisson error for each bin
    errors = np.sqrt(counts)

    # Calculate the relative error for each bin
    relative_errors = errors / counts

    # Define a threshold for a "large" uncertainty (e.g., 0.2 for 20%)
    threshold = 350

    # Check if any bin has a relative error larger than the threshold
    has_large_uncertainty = np.any(relative_errors > threshold)
    if 0:#N_sig.getValV()<threshold:
        out={
        #"mass"       : mass, 
        #"Signal"        : sig,  
        #"Background"        : bkg,  
        #"model":    model,
        #"data_sig_values":data_sig_values,
        #"data_bkg_values":data_bkg_values,
        #"tot_sig":tot_sig,
        #"tot_bkg":tot_bkg,
    
        "N_sig"      : 0, 
        "N_sig_error":0, 
        "N_bkg"      : 0,
        "N_bkg_error":0,
        "chi2"     :chi2,
        
        "status":status,
        #"uncert": uncert,
        "counts":counts
        }


    return out    

