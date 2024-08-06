from ROOT import * 
from array import array

_pt_binning=array('d',[200,1249,8000])

_eta_binning=array('d',[2,3.2,4.8])
beam="Pbp"
h=TH2D("effi_TM","",len(_pt_binning)-1,_pt_binning,len(_eta_binning)-1,_eta_binning)
for eta in range(len(_eta_binning)-1):
    for pt in range(len(_pt_binning)-1):
        pt_center=(_pt_binning[pt]+_pt_binning[pt+1])/2
        eta_center=(_eta_binning[eta]+_eta_binning[eta+1])/2
        bin=h.FindBin(pt_center,eta_center)
        if beam =="pp":
            if (pt_center >=200 and pt_center <= 1249) and (eta_center>=2 and eta_center<= 3.2):
                h.SetBinContent(bin,0.9944)
                h.SetBinError(bin,0.0015)
            elif (pt_center >=200 and pt_center <= 1249) and (eta_center>=3.2 and eta_center<= 4.8):
                h.SetBinContent(bin,0.9913)
                h.SetBinError(bin,0.0012)
            elif (pt_center >=1249 and pt_center <= 8000) and (eta_center>=2 and eta_center<= 3.2): 
                h.SetBinContent(bin,0.9903)
                h.SetBinError(bin,0.0017)
            elif (pt_center >=1249 and pt_center <= 8000) and (eta_center>=3.2 and eta_center<= 4.8):  
                h.SetBinContent(bin,0.9870)
                h.SetBinError(bin,0.0021)
        elif beam =="pPb":
            if (pt_center >=200 and pt_center <= 1249) and (eta_center>=2 and eta_center<= 3.2):
                h.SetBinContent(bin,0.99029)
                h.SetBinError(bin,0.00056)
            elif (pt_center >=200 and pt_center <= 1249) and (eta_center>=3.2 and eta_center<= 4.8):
                h.SetBinContent(bin,0.98873)
                h.SetBinError(bin,0.00059)
            elif (pt_center >=1249 and pt_center <= 8000) and (eta_center>=2 and eta_center<= 3.2): 
                h.SetBinContent(bin,0.98708)
                h.SetBinError(bin,0.00014)
            elif (pt_center >=1249 and pt_center <= 8000) and (eta_center>=3.2 and eta_center<= 4.8):  
                h.SetBinContent(bin,0.985514)
                h.SetBinError(bin,0.000082)
        elif beam =="Pbp":
            if (pt_center >=200 and pt_center <= 1249) and (eta_center>=2 and eta_center<= 3.2):
                h.SetBinContent(bin,0.98942)
                h.SetBinError(bin,0.00070)
            elif (pt_center >=200 and pt_center <= 1249) and (eta_center>=3.2 and eta_center<= 4.8):
                h.SetBinContent(bin,0.98785)
                h.SetBinError(bin,0.00062)
            elif (pt_center >=1249 and pt_center <= 8000) and (eta_center>=2 and eta_center<= 3.2): 
                h.SetBinContent(bin,0.98358)
                h.SetBinError(bin,0.00082)
            elif (pt_center >=1249 and pt_center <= 8000) and (eta_center>=3.2 and eta_center<= 4.8):  
                h.SetBinContent(bin,0.98267)
                h.SetBinError(bin,0.00012)

f_out=TFile("/home3/sara.sellam/RpPb_identified_hadrons_project/efficiencies/truth_mathing_effi/effi_TM_"+beam+".root","recreate")
h.SetDirectory(0)
f_out.cd()
h.Write()
f_out.Close()
      
            