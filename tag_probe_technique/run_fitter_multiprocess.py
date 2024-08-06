 
import sys
import glob
from math import *
from numpy import tan, arcsin 
import json
import datetime
import os
from array import array
import multiprocessing


from ROOT import *
from ROOT import RDF
from ROOT import RDataFrame
TH1.SetDefaultSumw2(True)    
TH2.SetDefaultSumw2(True)
gStyle.SetPaintTextFormat("1.3f")
gStyle.SetPalette(1)
sys.path.append("/home3/sara.sellam/RpPb_identified_hadrons_project")

from Binning.Binning_calibration import binning_TagandProbe
from Binning.Binning import *
from efficiencies.pid_effi.pid_selections import _pid_selections
from Ntuple_path.call_calibration_ntuple  import _part,Tuple_Names_for_TagandProbe
from selections.selections import _params,_pid_param
from Ntuple_path.call_ntuple import _corrections



from Fit_Model import *
from fit_configuration import fit_opt

from Helpers_function.sys_functions import weighted_std, gen_Toys
         

git_out_path="/scratch43/ssellam/results"
git_path="/home3/sara.sellam/RpPb_identified_hadrons_project"
ROOT.EnableImplicitMT()
in_range = lambda iq, imin, imax : "{} > {} && {} < {}".format(iq,imin, iq, imax)

def process_files_from_json(file_path, beam ="",data_type="",num_files_to_use=1):
    data_list = get_dictionary_from_json(file_path)

    
    if data_list is None:
        return

    if beam not in data_list:
        print(f"Error: '{beam}' not found in the JSON data for beam values.")
        return

    files_list = data_list[beam][data_type]
    if len(files_list)>2:
        desired_file = files_list[0:num_files_to_use]
    else:
        desired_file=files_list
    return desired_file

#useful functions 
def get_dictionary_from_json(file_path):
    try:
        with open(file_path) as file:
            data = json.load(file)
        return data
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
        return None
    except json.JSONDecodeError:
        print(f"Error: Invalid JSON format in file '{file_path}'.")
        return None


def calculate_mass_modelling(args):
    i, j, df_slice, fit_opt, weights, binned, beam, mom,outFile_name,sig_bkg = args
    result= mass_modelling(df_slice, fit_opt.loc[mom], weights=weights, binned=binned, beam=beam, mom=mom,outFile_name=outFile_name,i=i,j=j,sig_bkg=sig_bkg)
    
    return i,j,result



_part_latex={"Ks":"m(\pi^{+}\pi^{-})",
             "Lambda":"m(p \pi)",
             "phi":"m(K^{+} K^{-})"}

def mass(beam="",mom="",probe="",pid_cut="",nbr_weight="1",cuts="",xvar="",binning_mode="oscar_ana",binned=False,weights=False,x_y_used=""):
    print(" start calculations for : "+beam+"  probe : "+probe+" pid_cut : "+pid_cut)
    particle_name={"pi":"#pi","K":"K","p":"p"}

    tuple_name=_part[beam][mom]['Tuple']

    in_tree_list=Tuple_Names_for_TagandProbe[beam][tuple_name]
    
    data_list=process_files_from_json(git_path+"/Ntuple_path/calibration_samples.json",data_type="Data",beam=beam,num_files_to_use=20 )
   
    _pt_bins_ana=binning[binning_mode][beam]["pt"]
    _eta_bins_ana=binning[binning_mode][beam]["eta"]

    in_tree=TChain()
   

    if len(data_list)>1:
        reduced_data=[data_list[0],data_list[1],data_list[2],data_list[3],data_list[4],data_list[5],data_list[6],data_list[7],data_list[8],data_list[9],data_list[10],data_list[11],data_list[12],data_list[13],data_list[14],data_list[15]]
        for data in reduced_data:
            for tree in in_tree_list:
                print("tree",tree)
                in_tree.Add(data + "/" + tree + "/DecayTree")
    else:           
        for data in data_list:
            for tree in in_tree_list:
                print("tree",tree)
                in_tree.Add(data + "/" + tree + "/DecayTree")
               
    assert in_tree.GetEntries()!=0

    generate_array = lambda *nums: [f"{num}{suffix}" for num in nums for suffix in ["_pi", "_K", "_p"]]
    pid_weight_array=generate_array(0,1, 2,3,4,5,6,7)


    if probe!="p":

        df_raw=RDataFrame(in_tree)
        df_raw=df_raw.Filter("probe_P>7000&&probe_PT>500")
        dict_raw=df_raw.AsNumpy(["probe_P","probe_PT","probe_ETA","probe_TRACK_GhostProb","probe_hasRich","probe_PIDK","probe_PIDp",_part[beam][mom]["mom_mass"]])
        pd_raw=pd.DataFrame({"probe_P":dict_raw["probe_P"],
                         "probe_PT":dict_raw["probe_PT"],
                         "probe_ETA":dict_raw["probe_ETA"],
                         "probe_TRACK_GhostProb":dict_raw["probe_TRACK_GhostProb"],
                         "probe_hasRich":dict_raw["probe_hasRich"],
                         "probe_PIDK":dict_raw["probe_PIDK"],
                         "probe_PIDp":dict_raw["probe_PIDp"],
                         _part[beam][mom]["mom_mass"]:dict_raw[_part[beam][mom]["mom_mass"]]})

        if weights in pid_weight_array:
            data_list=glob.glob(data_list[0])[0]
            file_name=data_list.split("/")[-1]
            
            if weights in ["0_pi","0_K","0_p"]:
                pd_raw["w"]=np.ones(len(pd_raw))
            
            else:
                ch_weights=TChain()
                ch_weights.Add(_corrections[f"weights_for_pid_"+probe]["Calibration"][beam][0][0]+"w_"+weights+"_"+file_name+"/tree")
                ch_weights.Print()
                df_w_raw=RDataFrame(ch_weights)
                dict_w_raw=df_w_raw.AsNumpy(["w"])

                pd_raw["w"]=dict_w_raw["w"]

            dict_all=pd_raw.to_dict("list")
            
        else: 
            print("I am here")



    elif probe =="p":
        Kinematic_cuts="p_P>7000&&p_PT>500"
        ch_target1=TChain()
        ch_target1.Add("/scratch43/ssellam/calibration/Data/original_data/Lambda_ppi/La_pA_forPID.root/LambdaTuple_Line1/DecayTree")
        print("ch_target1",ch_target1.GetEntries())
        df_target1=RDataFrame(ch_target1)
        df_target1=df_target1.Filter(Kinematic_cuts+"&&p_P<40000")

        df_target1=df_target1.Define("p_ETA","-log(tan(asin(p_PT/p_P)/2))")
        df_target1_arr=df_target1.AsNumpy(["p_P","p_PT","p_ETA","p_TRACK_GhostProb","p_hasRich","p_PIDK","p_PIDp",_part[beam][mom]["mom_mass"]])
        pd_target1=pd.DataFrame(df_target1_arr)
        target1_arr_w=np.ones(len(pd_target1))
        target1_arr_w=target1_arr_w*(1/0.18)

        ch_target2=TChain()
        ch_target2.Add("/scratch43/ssellam/calibration/Data/original_data/Lambda_ppi/La_pA_forPID.root/LambdaTuple_Line2/DecayTree")
        print("ch_target2",ch_target2.GetEntries())
        df_target2=RDataFrame(ch_target2)
        df_target2=df_target2.Filter(Kinematic_cuts+"&&p_P>40000")
        df_target2=df_target2.Define("p_ETA","-log(tan(asin(p_PT/p_P)/2))")
        df_target_arr2=df_target2.AsNumpy(["p_P","p_PT","p_ETA","p_TRACK_GhostProb","p_hasRich","p_PIDK","p_PIDp",_part[beam][mom]["mom_mass"]])
        pd_target2=pd.DataFrame(df_target_arr2)
        target2_arr_w=np.ones(len(pd_target2))
        target = pd.concat([pd_target1, pd_target2])
        target_arr_w = np.concatenate((target1_arr_w, target2_arr_w))
        target["scale_w"]=target_arr_w
        print("pd_target1",len(pd_target1),"ch_target1",ch_target1.GetEntries(Kinematic_cuts+"&&p_P<40000")," all ",ch_target1.GetEntries())
        print("pd_target2",len(pd_target2),"ch_target2",ch_target2.GetEntries(Kinematic_cuts+"&&p_P>40000")," all ",ch_target2.GetEntries())
        print("calib",len(target)) 

        if weights in pid_weight_array:
            if weights in ["0_p"]:
                target["w"]=target["scale_w"]
              
            else:
                data_list=glob.glob("/scratch43/ssellam/calibration/Data/calibration2013/pPb_sample.pidcalib.root")[0]
                file_name=data_list.split("/")[-1]
                ch_weights=TChain()
                ch_weights.Add(_corrections[f"weights_for_pid_"+probe]["Calibration"][beam][0][0]+"w_"+weights+"_"+file_name+"/tree")
                print(_corrections[f"weights_for_pid_"+probe]["Calibration"][beam][0][0]+"w_"+weights+"_"+file_name)
                df_w_raw=RDataFrame(ch_weights)
                dict_w_raw=df_w_raw.AsNumpy(["w"])
                print("weights",len(dict_w_raw["w"]))
                print("scale_w",len(target["scale_w"]))


                target["w"]=dict_w_raw["w"]*target["scale_w"]
            target["probe_P"]=target["p_P"]
            target["probe_PT"]=target["p_PT"]
            target["probe_ETA"]=target["p_ETA"]
            target["probe_TRACK_GhostProb"]=target["p_TRACK_GhostProb"]
            target["probe_hasRich"]=target["p_hasRich"]
            target["probe_PIDK"]=target["p_PIDK"]
            target["probe_PIDp"]=target["p_PIDp"]
        


        dict_all=target.to_dict("list")
            





    df_denom=RDF.FromNumpy({"probe_P":np.array(dict_all["probe_P"]),
                            "probe_PT":np.array(dict_all["probe_PT"]),
                            "probe_ETA":np.array(dict_all["probe_ETA"]),
                            "w":np.array(dict_all["w"]),
                            "probe_TRACK_GhostProb":np.array(dict_all["probe_TRACK_GhostProb"]),
                            "probe_hasRich":np.array(dict_all["probe_hasRich"]).astype(int),
                            "probe_PIDK":np.array(dict_all["probe_PIDK"]),
                            "probe_PIDp":np.array(dict_all["probe_PIDp"]),
                            _part[beam][mom]["mom_mass"]:np.array(dict_all[_part[beam][mom]["mom_mass"]])})




    pid_selection_mode="ana"
    
    if len(cuts) > 0:
        df_denom=df_denom.Filter(cuts)
   

    
    results_before={}
    results_after={}
    ns_before,ns_after={},{}
    ns_before_err,ns_after_err={},{}
    
    ptext={}
    effi=probe+"_as_"+pid_cut
    MT=True
    if __name__ == '__main__':
        # Define the number of processes to use (adjust as needed)
        num_processes = 4
        

        
        # Create a pool of processes

        # Create a list of arguments for the function
        args_list = []
        result={}
        results_list={}
        results={}
        h,h_sys,chi2,status={},{},{},{}
        data_frame={}
        ns,ns_err,chi2_value,fit_status_value={},{},{},{}
        ratio,can,pad1,can_effi={},{},{},{}
        can_chi2,can_fitstatus={},{}

        for key in ["Before","After"]:
            ns[key],ns_err[key],chi2_value[key],fit_status_value[key]={},{},{},{}
            results_list[key]={}
            data_frame[key]={}
            can[key]={}
            can[key,"raw"]={}


        results["Before"],results["After"]={},{}
        
        
        fit_mode_dict={"ana":["sig","bkg"],
                       "ana2":["sig","bkg"],
                       "ana3":["sig","bkg"],
                 "Fit_sys":["sig_sys","bkg_sys"],
                 "Bin_sys_v1":["sig","bkg"],
                 "Bin_sys_v2":["sig","bkg"],
                 "Bin_sys_v3":["sig","bkg"],
                 "loose_pid_cuts_sys":["sig","bkg"],
                 "tight_pid_cuts_sys":["sig","bkg"],
                 "tight_old":["sig","bkg"]}
        
        binning_mode_dict={"ana":{"h1_h1":"tight",
                                "h1_h":"croaser"},
                          "ana2":{"h1_h1":"oscar_ana_tight",
                                    "h1_h":"oscar_ana_croaser"},
                           "ana3":{"h1_h1":"sara_ana_tight",
                                    "h1_h":"sara_ana_croaser"},
                     "Fit_sys":{"h1_h1":"tight",
                                "h1_h":"croaser"},
                        "Bin_sys_v1":{"h1_h1":"sys_tight_v1",
                                      "h1_h":"sys_croaser_v1"},
                        "Bin_sys_v2":{"h1_h1":"sys_tight_v2",
                                      "h1_h":"sys_croaser_v2"},
                        "Bin_sys_v3":{"h1_h1":"sys_tight_v3",
                                      "h1_h":"sys_croaser_v3"},
                        "loose_pid_cuts_sys":{"h1_h1":"tight",
                                            "h1_h":"croaser"},
                        "tight_pid_cuts_sys":{"h1_h1":"tight",
                                            "h1_h":"croaser"},
                        "tight_old":{"h1_h1":"tight_old",
                                    "h1_h":"croaser"},
                        }
        
        pid_selection_mode_dict={"ana":"ana",
                                 "ana2":"ana",
                                 "ana3":"ana",
                                "Fit_sys":"ana",
                                "Bin_sys_v1":"ana",
                                "Bin_sys_v2":"ana",
                                "Bin_sys_v3":"ana",
                                "loose_pid_cuts_sys":"loose",
                                "tight_pid_cuts_sys":"tight",
                                "tight_old":"ana",
                                }

        
        fit_mode=["ana2"]#,"Fit_sys","Bin_sys_v1","Bin_sys_v2"]
        #fit_mode=["ana","Fit_sys","Bin_sys_v1","Bin_sys_v2","Bin_sys_v3","loose_pid_cuts_sys","tight_pid_cuts_sys"]

        for m in fit_mode:
            pid_selection_mode=pid_selection_mode_dict[m]
            df_num=df_denom.Filter(_pid_selections[pid_selection_mode][pid_cut][beam])
       
            if probe == pid_cut:
                binning_name=binning_mode_dict[m]["h1_h1"]
                _pt_bins=array('d',binning_TagandProbe[tuple_name][binning_name][xvar])
                _eta_bins=array('d',binning_TagandProbe[tuple_name][binning_name]["ETA"])
            else:
                binning_name=binning_mode_dict[m]["h1_h"]
                _pt_bins=array('d',binning_TagandProbe[tuple_name][binning_name][xvar])
                _eta_bins=array('d',binning_TagandProbe[tuple_name][binning_name]["ETA"])
            
            out_file=TFile(git_out_path+"/efficiencies/pid_effi/tag_probe_technique/rfiles/"+x_y_used+"/ratio/"+beam+"/ratio_"+effi+"_"+beam+"_"+binning_name+"_"+m+"_weight_"+weights+".root","recreate")
            for mode in ["Before","After"]:
                pool = multiprocessing.Pool(processes=num_processes)
                for key in ["Before","After"]:
                    results[key][m]={}
                    ns[key][m],ns_err[key][m],chi2_value[key][m],fit_status_value[key][m]={},{},{},{}
                
                data_frame["Before"][m]=df_denom
                data_frame["After"][m]=df_num
                sig_bkg=fit_mode_dict[m]
                outFile_name=git_out_path+"/efficiencies/pid_effi/tag_probe_technique/rfiles/"+x_y_used+"/"+mode+"/"+mode+"_canvas_"+beam+"_"+effi+"_"+m+"_weight_"+weights
                df_in_bins={}
                for i in range(len(_pt_bins)-1):
                    for j in range(len(_eta_bins)-1):

                        pt_bin_size=in_range("probe_"+xvar,str(_pt_bins[i]),str(_pt_bins[i+1]))
                        probe_eta="(-log(tan(asin(probe_PT/probe_P)/2)))"
                        eta_bin_size=in_range(probe_eta,str(_eta_bins[j]),str(_eta_bins[j+1]))
                        bin_size=pt_bin_size+"&&"+eta_bin_size
                        df_slice=data_frame[mode][m].Filter(bin_size)
                        df_in_bins[i,j]=df_slice

                        if weights!="0":
                            if weights in pid_weight_array:
                                df_slice=df_slice.AsNumpy([_part[beam][mom]["mom_mass"],"w"])
                            else:
                                df_slice=df_slice.AsNumpy([_part[beam][mom]["mom_mass"],"probe_before_sw_w"+weights])
                        else:
                            df_slice=df_slice.AsNumpy([_part[beam][mom]["mom_mass"]])
                        
                        args_list.append((i, j, df_slice, fit_opt, weights, binned, beam, mom,outFile_name,sig_bkg))
                
            
                # Use the pool to parallelize the calculations
                results_list[mode][m] = pool.map(calculate_mass_modelling, args_list)
                
                # Close the pool of processes
                pool.close()
                pool.join()
                can[mode][m]=TCanvas("can","",0,0,700,1000)
                
                can[mode][m].Divide(len(_eta_bins)-1, len(_pt_bins)-1)
                fit_var={}
                frame,data,PDF_mass,c_file={},{},{},{}
                for i in range(len(_pt_bins)-1):
                    for j in range(len(_eta_bins)-1):
                        c_file[i,j]=TFile(outFile_name+"_"+str(i)+"_"+str(j)+".root","open")
                        frame[i,j] = c_file[i,j].Get("frame")
                        data[i,j] = c_file[i,j].Get("data")
                        PDF_mass[i,j] = c_file[i,j].Get("PDF_mass")
                        PDF_mass[i,j].plotOn(frame[i,j], ROOT.RooFit.Name('PDF_mass'), ROOT.RooFit.MarkerColor(ROOT.kRed))
                        data[i,j].plotOn(frame[i,j], ROOT.RooFit.Name('data'))
                        can[mode][m].cd(1+(len(_eta_bins)-1)*i+j)
                        frame[i,j].Draw()
                        can[mode][m].Update()
                        ptext[1+(len(_eta_bins)-1)*i+j] = TPaveText(0.65,0.65,0.85,0.85,"NDC")
                        ptext[1+(len(_eta_bins)-1)*i+j].SetFillStyle(4000)
                        ptext[1+(len(_eta_bins)-1)*i+j].SetBorderSize(0)
                        ptext[1+(len(_eta_bins)-1)*i+j].AddText("{} < #eta <{}".format(_eta_bins[j],_eta_bins[j+1]))
                        if x_y_used=="PT_eta":
                            ptext[1+(len(_eta_bins)-1)*i+j].AddText("{} < p_{{T}} < {}".format(_pt_bins[i],_pt_bins[i+1]))
                        else:
                            ptext[1+(len(_eta_bins)-1)*i+j].AddText("{} < p < {}".format(_pt_bins[i],_pt_bins[i+1]))
                        
                        ptext[1+(len(_eta_bins)-1)*i+j].SetFillStyle(0)
                        ptext[1+(len(_eta_bins)-1)*i+j].Draw()
                        os.remove(outFile_name+"_"+str(i)+"_"+str(j)+".root")

                        can_inv=ROOT.TCanvas("can_"+str(i)+"_"+str(j),"", 800, 600)

                        box = ROOT.TPaveText(0.6, 0.65, 0.85, 0.8, 'NDC')

                        box.SetFillColor(0)
                        box.SetBorderSize(1)
                        box.AddText("{} < #eta< {}".format(_eta_bins[j],_eta_bins[j+1]))
                        if x_y_used=="PT_eta":
                            box.AddText("${} < p_{{T}} < {}$".format(_pt_bins[i], _pt_bins[i+1]))
                        else:
                            box.AddText("{} < p < {}".format(_pt_bins[i],_pt_bins[i+1]))
                   


                        can_inv.cd()
                        frame[i,j].Draw()
                        box.Draw()
                        can_inv.SaveAs(git_out_path+"/efficiencies/pid_effi/tag_probe_technique/Plots/"+x_y_used+"/"+mode+"/"+beam+"/individual_canvas/"+mode+"_canvas_"+beam+"_"+effi+"_"+m+"_weight_"+weights+"_"+"_eta_"+str(_eta_bins[j])+"_"+str(_eta_bins[j+1])+"_pt_"+str(_pt_bins[i])+"_"+str(_pt_bins[i+1])+".pdf")

                    """
                    for i in range(len(_pt_bins)-1):
                        for j in range(len(_eta_bins)-1):
                            if int(weights) !=0:
                                h_p=df_in_bins[i,j].Histo1D(("h_p_"+str(i)+"_"+str(j),"",100,7000,10000),"probe_P","w")
                                h_pt=df_in_bins[i,j].Histo1D(("h_pt_"+str(i)+"_"+str(j),"",20,_pt_bins[i],_pt_bins[i+1]),"probe_PT","w")
                                h_p.SetStats(0)
                                h_pt.SetStats(0)

                                can_inv=ROOT.TCanvas("can_"+str(i)+"_"+str(j),"", 800, 600)

                                box = ROOT.TPaveText(0.6, 0.65, 0.85, 0.8, 'NDC')

                                box.SetFillColor(0)
                                box.SetBorderSize(1)
                                box.AddText("{} < #eta< {}".format(_eta_bins[j],_eta_bins[j+1]))
                                if x_y_used=="PT_eta":
                                    box.AddText("${} < p_{{T}} < {}$".format(_pt_bins[i], _pt_bins[i+1]))
                                else:
                                    box.AddText("{} < p < {}".format(_pt_bins[i],_pt_bins[i+1]))
                                can_inv.cd()
                                h_p.Draw()
                                box.Draw()
                                
                                
                                os.makedirs(git_out_path+"/efficiencies/pid_effi/tag_probe_technique/Plots/"+x_y_used+"/"+mode+"/"+beam+"/individual_canvas/P", exist_ok=True)
                                can_inv.SaveAs(git_out_path+"/efficiencies/pid_effi/tag_probe_technique/Plots/"+x_y_used+"/"+mode+"/"+beam+"/individual_canvas/P/"+mode+"_canvas_"+beam+"_"+effi+"_"+m+"_weight_"+weights+"_"+"_eta_"+str(_eta_bins[j])+"_"+str(_eta_bins[j+1])+"_pt_"+str(_pt_bins[i])+"_"+str(_pt_bins[i+1])+".pdf")
                                
                                can_inv=ROOT.TCanvas("can_"+str(i)+"_"+str(j),"", 800, 600)

                                box = ROOT.TPaveText(0.6, 0.65, 0.85, 0.8, 'NDC')

                                box.SetFillColor(0)
                                box.SetBorderSize(1)
                                box.AddText("{} < #eta< {}".format(_eta_bins[j],_eta_bins[j+1]))
                                if x_y_used=="PT_eta":
                                    box.AddText("${} < p_{{T}} < {}$".format(_pt_bins[i], _pt_bins[i+1]))
                                else:
                                    box.AddText("{} < p < {}".format(_pt_bins[i],_pt_bins[i+1]))
                                can_inv.cd()
                                h_pt.Draw()
                                box.Draw()
                                os.makedirs(git_out_path+"/efficiencies/pid_effi/tag_probe_technique/Plots/"+x_y_used+"/"+mode+"/"+beam+"/individual_canvas/PT", exist_ok=True)
                                can_inv.SaveAs(git_out_path+"/efficiencies/pid_effi/tag_probe_technique/Plots/"+x_y_used+"/"+mode+"/"+beam+"/individual_canvas/PT/"+mode+"_canvas_"+beam+"_"+effi+"_"+m+"_weight_"+weights+"_"+"_eta_"+str(_eta_bins[j])+"_"+str(_eta_bins[j+1])+"_pt_"+str(_pt_bins[i])+"_"+str(_pt_bins[i+1])+".pdf")
                    """

                    
                    


                
                can[mode][m].SaveAs(git_out_path+"/efficiencies/pid_effi/tag_probe_technique/Plots/"+x_y_used+"/"+mode+"/"+beam+"/"+mode+"_canvas_"+beam+"_"+effi+"_"+m+"_weight_"+weights+".pdf")

                for i in range(len(_pt_bins)-1):
                    for j in range(len(_eta_bins)-1):      
                        h[mode]=TH2D("h_"+mode+"_"+m+"_eta_"+str(_eta_bins[j])+"_"+str(_eta_bins[j+1])+"_pt_"+str(_pt_bins[i])+"_"+str(_pt_bins[i+1])+"".format(),"",len(_pt_bins)-1,_pt_bins,len(_eta_bins)-1,_eta_bins)
                        h[mode].SetDirectory(0)
                        h[mode,"raw"]=TH2D("h_raw_"+mode+"_"+m+"_eta_"+str(_eta_bins[j])+"_"+str(_eta_bins[j+1])+"_pt_"+str(_pt_bins[i])+"_"+str(_pt_bins[i+1])+"".format(),"",len(_pt_bins)-1,_pt_bins,len(_eta_bins)-1,_eta_bins)
                        h[mode,"raw"].SetDirectory(0)
                        chi2[mode]=TH2D("h_chi2_"+mode+"_"+m+"_eta_"+str(_eta_bins[j])+"_"+str(_eta_bins[j+1])+"_pt_"+str(_pt_bins[i])+"_"+str(_pt_bins[i+1])+"".format(),"",len(_pt_bins)-1,_pt_bins,len(_eta_bins)-1,_eta_bins)
                        chi2[mode].SetDirectory(0)
                        status[mode]=TH2D("h_status_"+mode+"_"+m+"_eta_"+str(_eta_bins[j])+"_"+str(_eta_bins[j+1])+"_pt_"+str(_pt_bins[i])+"_"+str(_pt_bins[i+1])+"".format(),"",len(_pt_bins)-1,_pt_bins,len(_eta_bins)-1,_eta_bins)
                        status[mode].SetDirectory(0)
                        c_file[i,j].Close()
                output_dict={}
                for i, j, result_dict in results_list[mode][m]:
                    

                    results[mode][m][i, j] = result_dict
                    ns[mode][m][i, j] = result_dict["N_sig"]
                    ns_err[mode][m][i, j] = result_dict["N_sig_error"]
                    chi2_value[mode][m][i, j]=result_dict["chi2"]
                    fit_status_value[mode][m][i, j]=result_dict["status"]
                    key_str = f"{i}_{j}" 
                    output_dict[key_str]={} 
                    output_dict[key_str]["ns_err"]=str(ns_err[mode][m][i, j])
                    output_dict[key_str]["ns"]=str(ns[mode][m][i, j])
                    output_dict[key_str]["nbkg_err"]=result_dict["N_bkg_error"]
                    output_dict[key_str]["nbkg"]=result_dict["N_bkg"]
                    output_dict[key_str]["counts"]=np.sum(result_dict["counts"])
                    # some tests
                    
                    combined_uncertainty = np.sqrt(ns_err[mode][m][i, j]**2 + result_dict["N_bkg_error"]**2)
                    output_dict[key_str]["combined_uncertainty"]=combined_uncertainty

                    mean,mean_err,sigma,sigma_err,canvas=gen_Toys(ns[mode][m][i, j],ns_err[mode][m][i, j],beam,_eta_bins,_pt_bins,pol="Up+Down",can_name="test",xlabel="test")
                    output_dict[key_str]["toy_sigma"]=sigma
                    tot_unc=np.sqrt(sigma**2+ns_err[mode][m][i, j]**2)
                    output_dict[key_str]["toy_sigma+ns_err"]=np.sqrt(sigma**2+ns_err[mode][m][i, j]**2)



                    
                    if ns[mode][m][i, j]<100:
                        h[mode].SetBinContent(i + 1, j + 1, 0)
                        h[mode].SetBinError(i + 1, j + 1,0)
                    else:
                        h[mode].SetBinContent(i + 1, j + 1, ns[mode][m][i, j])
                        h[mode].SetBinError(i + 1, j + 1, ns_err[mode][m][i, j])

                    h[mode,"raw"].SetBinContent(i + 1, j + 1, ns[mode][m][i, j])
                    h[mode,"raw"].SetBinError(i + 1, j + 1, ns_err[mode][m][i, j])
                    
                    """
                    h[mode].SetBinContent(i + 1, j + 1, ns[mode][m][i, j])

                    )
                    """
                    #h[mode].SetBinError(i + 1, j + 1, combined_uncertainty+(0.01*ns[mode][m][i, j]))
                    #h[mode].SetBinError(i + 1, j + 1,tot_unc)
                    if isnan(chi2_value[mode][m][i, j]):
                        chi2_value[mode][m][i, j] = -1
                    chi2[mode].SetBinContent(i + 1, j + 1, chi2_value[mode][m][i, j])
                    status[mode].SetBinContent(i + 1, j + 1, fit_status_value[mode][m][i, j] + 1)
                out_file.cd()
                h[mode].Write()
                chi2[mode].Write()
                status[mode].Write()    
                """
                with open("output_stat_uncert_new.json", "w") as json_file:
                    json.dump(output_dict, json_file, indent=4)     
                """   
            ratio[m]=TH2D("ratio","",len(_pt_bins)-1,_pt_bins,len(_eta_bins)-1,_eta_bins)
            ratio[m,"raw"]=TH2D("ratio_raw","",len(_pt_bins)-1,_pt_bins,len(_eta_bins)-1,_eta_bins)

            #ratio[m]=TEfficiency(h["After"],h["Before"])
            ratio[m].Divide(h["After"],h["Before"],1.,1.,"B")
            ratio[m,"raw"].Divide(h["After","raw"],h["Before","raw"],1.,1.,"B")

            #ratio[m].SetStatisticOption(ROOT.TEfficiency.kFWilson)

            frac={}
            
            for eta_bin in range(len(_eta_bins)-1):
                frac[eta_bin,eta_bin+1]=ratio[m].ProjectionX("ratio_"+m+"_{}_{}".format(_eta_bins[eta_bin],_eta_bins[eta_bin+1]),eta_bin,eta_bin+1)
                out_file.cd()
                frac[eta_bin,eta_bin+1].Write()
            
            can_effi[m]={}
            
            can_effi[m]["ratio"] = TCanvas("can_effi","",500,400)
            pad1[m] = TPad("pad","pad", 0.01, 0.01, 0.95, 0.95)
            pad1[m].Draw()
            pad1[m].cd()  
            ratio[m].SetXTitle("p_{T} [MeV/c]")
            ratio[m].SetYTitle("#eta")
            ratio[m].SetTitle(particle_name[probe]+" #rightarrow "+particle_name[pid_cut])
            ratio[m].Draw("colz,text,e")
            ratio[m].SetStats(0)
            can_effi[m]["ratio"].SaveAs(git_out_path+"/efficiencies/pid_effi/tag_probe_technique/Plots/"+x_y_used+"/ratio/"+beam+"/ratio_"+effi+"_"+beam+"_"+binning_name+"_"+m+"_"+weights+".pdf")
            
            can_effi[m,"raw"]={}
            can_effi[m,"raw"]["ratio"] = TCanvas("can_effi","",500,400)
            pad1[m,"raw"] = TPad("pad","pad", 0.01, 0.01, 0.95, 0.95)
            pad1[m,"raw"].Draw()
            pad1[m,"raw"].cd()  
            ratio[m,"raw"].SetXTitle("p_{T} [MeV/c]")
            ratio[m,"raw"].SetYTitle("#eta")
            ratio[m,"raw"].SetTitle(particle_name[probe]+" #rightarrow "+particle_name[pid_cut])
            ratio[m,"raw"].Draw("colz,text,e")
            ratio[m,"raw"].SetStats(0)
            can_effi[m,"raw"]["ratio"].SaveAs(git_out_path+"/efficiencies/pid_effi/tag_probe_technique/Plots/"+x_y_used+"/ratio/"+beam+"/ratio_raw_"+effi+"_"+beam+"_"+binning_name+"_"+m+"_"+weights+".pdf")
            
            
            
            for h_key in ["Before","After"]:
                can_effi[m][h_key] = TCanvas("can_effi","",500,400)
                pad1[m] = TPad("pad","pad", 0.01, 0.01, 0.95, 0.95)
                pad1[m].Draw()
                pad1[m].cd()  
                if x_y_used=="PT_eta":
                    h[h_key].SetXTitle("p_{T} [MeV/c]")
                else:
                    h[h_key].SetXTitle("p [MeV/c]")
                h[h_key].SetYTitle("#eta")
                h[h_key].Draw("colz,text,e")
                h[h_key].SetStats(0)
                can_effi[m][h_key].SaveAs(git_out_path+"/efficiencies/pid_effi/tag_probe_technique/Plots/"+x_y_used+"/ratio/"+beam+"/"+h_key+"_"+effi+"_"+beam+"_"+binning_name+"_"+m+"_"+weights+".pdf")
                


            out_file.cd()
            ratio[m].Write()
            ratio[m,"raw"].Write()



            can_chi2[m] = TCanvas("can_chi2","",500,400)
            pad1[m] = TPad("pad","pad", 0.01, 0.01, 0.95, 0.95)
            pad1[m].Draw()
            pad1[m].cd()  
            chi2[mode].SetXTitle("p_{T} [MeV/c]")
            chi2[mode].SetYTitle("#eta")
            chi2[mode].Draw("colz,text,e")
            chi2[mode].SetStats(0)
            chi2[mode].SetMaximum(10)
            can_chi2[m].SaveAs(git_out_path+"/efficiencies/pid_effi/tag_probe_technique/Plots/"+x_y_used+"/ratio/"+beam+"/chi2_"+effi+"_"+beam+"_"+binning_name+"_"+m+"_"+weights+".pdf")
            out_file.cd()
            chi2[mode].Write()

            can_fitstatus[m] = TCanvas("can_fitstatus","",500,400)
            pad1[m] = TPad("pad","pad", 0.01, 0.01, 0.95, 0.95)
            pad1[m].Draw()
            pad1[m].cd()  
            status[mode].SetXTitle("p_{T} [MeV/c]")
            status[mode].SetYTitle("#eta")
            status[mode].Draw("colz,text,e")
            status[mode].SetStats(0)
            can_fitstatus[m].SaveAs(git_out_path+"/efficiencies/pid_effi/tag_probe_technique/Plots/"+x_y_used+"/ratio/"+beam+"/fit_status_"+effi+"_"+beam+"_"+binning_name+"_"+m+"_"+weights+".pdf")
            out_file.cd()
            status[mode].Write()
            out_file.Close()
            """
            pid_file=TFile(git_out_path+"/efficiencies/pid_effi/tag_probe_technique/rfiles/"+x_y_used+"/ratio_"+effi+"_"+beam+"_"+binning_mode+"_"+m+"_"+weights+".root","recreate")
            for eta in range(len(_eta_bins_ana)-1):
                h_pid=TH1D("ratio_ana_"+str(_eta_bins_ana[eta])+"_"+str(_eta_bins_ana[eta+1]),"",len(_pt_bins_ana)-1,_pt_bins_ana)
                eta_center=(_eta_bins_ana[eta]+_eta_bins_ana[eta+1])/2
                if beam=="pp":
                   eta_center=eta_center
                elif beam=="pPb":
                    eta_center=eta_center+0.5
                else:
                    eta_center=(eta_center+0.5)*-1
                for pt in range(len(_pt_bins_ana)-1):
                    pt_center=(_pt_bins_ana[pt]+_pt_bins_ana[pt+1])/2
                    print("pt_center",pt_center)
                    bin=ratio[m].FindBin(pt_center,eta_center)
                    print("bin",bin)
                    value=ratio[m].GetBinContent(bin)
                    value_err=ratio[m].GetBinError(bin)
                    print("value",value)
                    print("value",value_err)
                    bin2=h_pid.FindBin(pt_center)
                    print("bin2",bin2)
                    h_pid.SetBinContent(bin2,value)
                    h_pid.SetBinError(bin2,value_err)
                pid_file.cd()
                h_pid.Write()
            pid_file.Close()
            """

        
        





    
import argparse

# Create an ArgumentParser object to handle command line arguments
parser = argparse.ArgumentParser()

# Add the command line arguments
parser.add_argument("-b", "--beam", nargs='+', type=str, help="Beam type")
parser.add_argument("-mom", "--mother", nargs='+', type=str, help="The mother of the probe mom")
parser.add_argument("-probe", "--probe", nargs='+', type=str, help="The probe mom")
parser.add_argument("-pid", "--pid", nargs='+', type=str, help="requested pid selection")
parser.add_argument("-weights", "--weights", nargs='+', type=str, help="requested weights")

# Parse the command line arguments
args = parser.parse_args()

# Extract the values of the arguments
beam_type = args.beam[0]
mom=args.mother[0]
probe=args.probe[0]
pid=args.pid[0]
weights=args.weights[0]

# Print the values of the arguments
print(f"Beam type: {beam_type}")
print(f"Mother of the probe : {mom}")
print(f"The probe : {probe}")
print(f"requested PID selection : {pid}")
print(f"requested weights : {weights}")
print("efficiency name "+probe+" as "+pid)
Kinematic_cuts="probe_P>7000&&probe_PT>500&&probe_TRACK_GhostProb<"+_params["GhostP"][beam_type]
#mass(beam=beam_type,mom=mom,probe=probe,cuts=Kinematic_cuts,pid_cut=pid,xvar="PT",binned=True,weights=weights,x_y_used="PT_eta")
mass(beam=beam_type,mom=mom,probe=probe,cuts=Kinematic_cuts,pid_cut=pid,xvar="P",binned=True,weights=weights,x_y_used="P_eta")
