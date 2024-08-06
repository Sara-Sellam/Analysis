

from ROOT import RDF
from ROOT  import *
from ROOT import RDataFrame


import time
import concurrent.futures
from multiprocessing import Process,Queue
from subprocess import call
from array import array
from math import *
import numpy as np
from numpy import log,tan,sqrt
import glob
import os.path


import sys
sys.path.append("/home3/sara.sellam/RpPb_identified_hadrons_project")
from Ntuple_path.call_ntuple import _path,_corrections_v2,get_data
from Binning.Binning import *

opts = RDF.RSnapshotOptions()
opts.fMode = "UPDATE"                                                                                                        
def chunks(in_list, nbr_chunks):
    chunk_size = int(len(in_list) / nbr_chunks)
    remainder = len(in_list) % nbr_chunks
    chunks_set = [in_list[i:i + chunk_size] for i in range(0, len(in_list) - remainder, chunk_size)]
    # Adjust the last chunk to include the remainder elements
    chunks_set[-1] = in_list[-(chunk_size + remainder):] 
    return  chunks_set
def check_root_file(file_path):
    if os.path.exists(file_path):
        root_file = TFile(file_path)
        if root_file.IsZombie():
            print("Error: The file is a zombie root file.")
            return False
        root_file.Close()
        return True
    else:
        print("Error: File does not exist.")
        return False

def add_corr(beam,pol,ratio,in_file,path_out):
    output_file_name=(in_file.split(beam+"_"+pol+"/"))[1]
    if check_root_file(path_out+"track_reco_corr_"+output_file_name):
        print("file exists ")
    else:
        MC_in_tree=TChain("MC_in_tree","MC_in_tree")
        MC_in_tree.Add(in_file+"/HITuple/DecayTree")
        entries=MC_in_tree.GetEntries()
        momentum=np.zeros(entries)

        Trcalib_ratio=np.zeros(entries)
        Trcalib_ratio_err=np.zeros(entries)
        Trcalib_Sys=np.zeros(entries)
        had_Sys=np.zeros(entries)
        total_err=np.zeros(entries)
    
        Sys = {"pPb": 0.004,"Pbp": 0.004, "pp":0.008 }
        hadSys = 0.014

        for i in range( MC_in_tree.GetEntries()):
            MC_in_tree.GetEntry(i)
            eta =(-log(tan(asin(MC_in_tree.pi_PT/MC_in_tree.pi_P))))
            if beam=="pp":
                bin =ratio.FindBin(MC_in_tree.pi_P,eta)   
            else:       
                bin =ratio.FindBin(MC_in_tree.pi_P/1000,eta)

            momentum[i]=MC_in_tree.pi_P
            Trcalib_ratio[i]= ratio.GetBinContent(bin) 
            Trcalib_ratio_err[i] =ratio.GetBinError(bin)
            if Trcalib_ratio[i]==0:
                Trcalib_ratio[i]=1
                Trcalib_ratio_err[i]=0.05
            Trcalib_Sys[i]=Sys[beam]
            had_Sys[i]=hadSys
            total_err[i]=( Trcalib_ratio_err[i]**2+(Sys[beam]*Trcalib_ratio[i] )**2+( hadSys*Trcalib_ratio[i] )**2 )**0.5 
        arr_Trcalib_ratio=np.asarray(Trcalib_ratio)
        arr_Trcalib_ratio_err=np.asarray(Trcalib_ratio_err)
        arr_Trcalib_Sys=np.asarray(Trcalib_Sys)
        arr_had_Sys=np.asarray(had_Sys)
        arr_total_err=np.asarray(total_err)
        dict_Tr_corr_sys={"Trcalib_ratio":arr_Trcalib_ratio,
                      "Trcalib_ratio_err":arr_Trcalib_ratio_err,
                      "Trcalib_Sys":arr_Trcalib_Sys,
                      "had_Sys":arr_had_Sys,
                      "track_corr_total_err":arr_total_err
                      }
        df=RDF.FromNumpy(dict_Tr_corr_sys)
        df.Snapshot("DecayTree",path_out+"track_reco_corr_"+output_file_name)



track_corr_file_dict={"pPb":"/scratch38/HITuple/MC/track_reco_corr/Tr_reco_corr_table/ratio2012S20.root",
                "Pbp":"/scratch38/HITuple/MC/track_reco_corr/Tr_reco_corr_table/ratio2012S20.root",
                "pp":"/scratch38/HITuple/MC/track_reco_corr/Tr_reco_corr_table/Ratio_2015_5TeV_GP03_Sim09b_Long.root"}
def add_reco_effi_correction(beam="",pol="",nbr_chunks=40,MC=""):
    for in_files in _path[MC][beam][pol][0]:
        print("in_files",in_files)
        in_list=glob.glob(in_files)
    chunks_lists=list(chunks(in_list,nbr_chunks))

    file=TFile(track_corr_file_dict[beam],"open")
    ratio=file.Get("Ratio")
    path_out=_corrections_v2["track_reco_corr_sys"][MC][beam][pol][0][0]                                                                                                                                    
    queue=Queue()
    for i in range(nbr_chunks):
        print("########################################")
        print("###Begin with mini list number",i," ####")
        print("########################################")
        processes = []
        mini_list=chunks_lists[i]
        for i in range(len(mini_list)):
            process=Process(target=add_corr,args=[beam,pol,ratio,mini_list[i],path_out])
            processes.append(process)

        for process in processes:
            process.start()
            print("process start")

        time.sleep(30)
        for proc in processes:
            proc.join(timeout=0)
            while proc.is_alive():
                print("Job is not finished!")
                time.sleep(30)

        for proc in processes:
            proc.terminate()    
        




import sys
sys.path.append("/home3/sara.sellam/RpPb_identified_hadrons_project")

import argparse

# Create an ArgumentParser object to handle command line arguments
parser = argparse.ArgumentParser()

# Add the command line arguments
parser.add_argument("-beam", "--beam", nargs='+', type=str, help="Beam type")
parser.add_argument("-pol", "--pol", nargs='+', type=str, help="Polarization")
parser.add_argument("-nbr_chunks", "--nbr_chunks", nargs='+', type=int, help="Number of chunks")
parser.add_argument("-MC", "--MC", nargs='+', type=str, help="Which MC sample")

# Parse the command line arguments
args = parser.parse_args()

# Extract the values of the arguments
beam = args.beam[0]
pol=args.pol[0]
nbr_chunks=args.nbr_chunks[0]
MC= args.MC[0]
# Print the values of the arguments
print(f"Beam type: {beam}")
print(f"Polarization : {pol}")
print(f"MC : {MC}")

add_reco_effi_correction(beam=beam,pol=pol,nbr_chunks=nbr_chunks,MC=MC)
