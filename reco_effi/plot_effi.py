import sys
sys.path.append("/home3/sara.sellam/RpPb_identified_hadrons_project")
import uproot
import matplotlib.pyplot as plt
import mplhep as hep
from Binning.Binning import *
from Helpers_function.import_dependencies import *
from Helpers_function.labels import * 
git_out_path="/scratch43/ssellam/results"

import argparse

# Create an ArgumentParser object to handle command line arguments
parser = argparse.ArgumentParser()

# Add the command line arguments
parser.add_argument("-beam", "--beam", nargs='+', type=str, help="Beam type")
parser.add_argument("-pol", "--pol", nargs='+', type=str, help="Polarization")
parser.add_argument("-particle", "--particle", nargs='+', type=str, help="particle")
parser.add_argument("-binning_mode", "--binning_mode", nargs='+', type=str, help="binning_mode")
parser.add_argument("-weights", "--weights", nargs='+', type=str, help="Occupancy weight")
parser.add_argument("-MC", "--MC", nargs='+', type=str, help="Which MC sample")

# Parse the command line arguments
args = parser.parse_args()

# Extract the values of the arguments
beam = args.beam
pol=args.pol
particles=args.particle
binning_mode=args.binning_mode[0]
weights= args.weights

MC= args.MC

# Print the values of the arguments
print(f"Beam type: {beam}")
print(f"Polarization : {pol}")
print(f"Particles : {particles}")
print(f"binning_mode : {binning_mode}")
   





def plot_efficiency(beam="", pol="", data1="", particle=[], track_type="", weights="", data2="", binning_mode="", ratio="", xlabel="$p_{T}$[MeV/c]", ylabel=r"$\varepsilon_{reco}$", legend=""):
    
    fig, axarr = plot_2_3_subplot(ylabel=ylabel, xlabel=pt_label)

    colors, markers = listof_colors_markers()
    position={"Pbp":{"-3.0_-2.5":[0,0],
                     "-3.5_-3.0":[0,1],
                     "-4.0_-3.5":[0,2],
                     "-4.5_-4.0":[1,0],
                     "-4.8_-4.5":[1,1],
                     "-5.2_-4.8":[1,2]}}
    if isinstance(particle, list):
        _pt_bins = binning[binning_mode][beam]["pt"]
        _eta_bins = binning[binning_mode][beam]["eta"] 
        for particle in particles:
            
            index = particles.index(particle)
            print(ratio+"_"+track_type+"_"+particle+"_"+beam)
            in_file = uproot.open(git_out_path+"/efficiencies/reco_effi/rfiles/"+ratio+"_"+track_type+"_"+particle+"_"+beam+"_"+pol+"_w_"+weights+"_"+data1+"_"+data2+"_"+binning_mode+".root")
            print(git_out_path+"/efficiencies/reco_effi/rfiles/"+ratio+"_"+track_type+"_"+particle+"_"+beam+"_"+pol+"_w_"+weights+"_"+data1+"_"+data2+"_"+binning_mode+".root")
            s = 0
            for i in range(axarr.ndim):
                for j in range(3):
                    removed_bins_array=np.array(removed_bins[beam][particle]["{}_{}".format(_eta_bins[s], _eta_bins[s+1])])
                    if beam == "Pbp":
                        i, j = position["Pbp"]["{}_{}".format(_eta_bins[s], _eta_bins[s+1])]
                        ax = axarr[i][j]
                    else:
                        ax = axarr[i][j]
                    #h_name = ratio + "_" + track_type + "_" + particle + "_{}_{}".format(_eta_bins[s], _eta_bins[s+1])
                    h_name="reco_effi_prompt_"+particle
                    h_errors = in_file[h_name].to_boost().variances()[:,s]*removed_bins_array
                    h_values = in_file[h_name].to_boost().counts()[:,s]*removed_bins_array
  
                    pt_bins= np.asarray(_pt_bins)
                    h_centers = (pt_bins[:-1] + pt_bins[1:]) / 2
                    h_edges = pt_bins[1:] - pt_bins[:-1]
                    #ax.text(0.05, 0.95, beam_name[beam] + r": $ {} <  \eta_{{cms}} < {}$".format(_eta_bins[s], _eta_bins[s+1]), transform=ax.transAxes, fontsize=30, verticalalignment='top')
                    
                    ax.text(0.40, 0.3, beam_name[beam] + r": $ {} <  \eta_{{cms}} < {}$".format(_eta_bins[s], _eta_bins[s+1]), transform=ax.transAxes, fontsize=35, verticalalignment='top')

                    ax.errorbar(x=h_centers, y=h_values, xerr=h_edges*0.5, yerr=h_errors, marker=particles_marker[particle], color=particles_color[particle], markersize=16, linestyle='none', label=particles_name[particle])
                    ax.set_xscale("log")
                    ax.set_ylim(0.2,0.95)
                    
                    ax.set_xlim(400, 4619)
                    ax.grid(True)
                    ax.axhline(y=1.0, color='Black', linestyle='--')
                    ax.tick_params(axis='x', labelsize=50)  
                    ax.tick_params(axis='y', labelsize=50) 
                    
                    s += 1
            if beam == "Pbp":
                i, j = position["Pbp"]["{}_{}".format(_eta_bins[0], _eta_bins[1])]
                ax = axarr[i][j]
            else:
                ax= axarr[0][0]
            #ax.text(0.6, 0.98,"Sim : "+data1, transform=ax.transAxes, fontsize=20, verticalalignment='top')
            #ax.text(0.6, 0.935,"Pol : "+pol, transform=ax.transAxes, fontsize=20, verticalalignment='top')
            #ax.text(0.6, 0.89,"Weight set : "+weights, transform=ax.transAxes, fontsize=20, verticalalignment='top')        
            leg = ax.legend(loc='upper left',fontsize=35)
            if weights !="0":
                title_text = "Sim: {}, Pol: {}, Weight set: {}".format(data1, pol, weights)
            else:
                title_text = "Sim: {}, Pol: {}, Unweighted".format(data1, pol)

            fig.suptitle(title_text, fontsize=50)

        plt.savefig(git_out_path+"/efficiencies/reco_effi/Plots/"+beam+"/reco_"+beam+"_"+pol+"_w_"+weights+"_"+data1+".pdf")
    """
    elif isinstance(beam, list):
        for b in beam:
            _pt_bins = binning[binning_mode][b]["pt"]
            _eta_bins = binning[binning_mode][b]["eta"] 
            index = beam.index(b)
            _eta_bins = binning[binning_mode][b]["eta"] 
            in_file = uproot.open(git_out_path+"/efficiencies/reco_effi/rfiles/"+ratio+"_"+track_type+"_"+particle+"_"+b+"_"+pol+"_w_"+weights+"_"+data1+"_"+data2+"_"+binning_mode+".root")
            s = 0
            for i in range(axarr.ndim):
                for j in range(3):  # Use shape[i] to determine the range dynamically
                    if b == "Pbp":
                        i, j = position["Pbp"]["{}_{}".format(_eta_bins[s], _eta_bins[s+1])]
                        ax = axarr[i][j]
                    else:
                        ax = axarr[i][j]
                    h_name = ratio + "_" + track_type + "_" + particle + "_{}_{}".format(_eta_bins[s], _eta_bins[s+1])
                    h_errors = in_file[h_name].errors()
                    h_values = in_file[h_name].to_numpy()[0]
                    h_centers = in_file[h_name].to_boost().axes.centers[0]
                    x_bin_size = in_file[h_name].to_boost().axes.widths[0] / 2
                    ax.errorbar(x=h_centers, y=h_values, xerr=x_bin_size, yerr=h_errors, marker=markers[index], color=colors[index], markersize=15, linestyle='none', label=beam_name[b] + ": $ {}<  $\eta_{cms}$ < {}$".format(_eta_bins[s], _eta_bins[s+1]))
                    ax.set_xscale("log")
                    ax.set_ylim(0.2,0.95)
                    ax.set_xlim(400, 4619)
                    ax.grid(True)
                    ax.axhline(y=1.0, color='Black', linestyle='--')
                    leg = ax.legend()
                    s += 1
            if beam == "Pbp":
                i, j = position["Pbp"]["{}_{}".format(_eta_bins[0], _eta_bins[1])]
                ax = axarr[i][j]
            else:
                ax= axarr[0][0]
            ax.text(0.6, 0.98,"Sim : "+data1, transform=ax.transAxes, fontsize=20, verticalalignment='top')
            ax.text(0.6, 0.935,"Pol : "+pol, transform=ax.transAxes, fontsize=20, verticalalignment='top')
            if weights !="0":
                ax.text(0.6, 0.89,"Weight set : "+weights, transform=ax.transAxes, fontsize=20, verticalalignment='top')        
            else:
                ax.text(0.6, 0.89,"Unweighted",transform=ax.transAxes, fontsize=20, verticalalignment='top')        

            leg = ax.legend()

        plt.savefig(git_out_path+"/efficiencies/reco_effi/Plots/reco_effi_"+particle+"_"+weights+"_all_beams_w_"+weights+"_"+data1+".pdf")
    elif isinstance(weights, list) and isinstance(beam, str) and isinstance(particle, str) :
        for w in weights:
            index = weights.index(w)
            _eta_bins = binning[binning_mode][beam]["eta"] 
            in_file = uproot.open(git_out_path+"/efficiencies/reco_effi/rfiles/"+ratio+"_"+track_type+"_"+particle+"_"+beam+"_"+pol+"_w_"+w+"_"+data1+"_"+data2+"_"+binning_mode+".root")
            s = 0
            for i in range(axarr.ndim):
                for j in range(3):  
                    if beam == "Pbp":
                        i, j = position["Pbp"]["{}_{}".format(_eta_bins[s], _eta_bins[s+1])]
                        ax = axarr[i][j]
                    else:
                        ax = axarr[i][j]
                    h_name = ratio + "_" + track_type + "_" + particle + "_{}_{}".format(_eta_bins[s], _eta_bins[s+1])
                    h_errors = in_file[h_name].errors()
                    h_values = in_file[h_name].to_numpy()[0]
                    h_centers = in_file[h_name].to_boost().axes.centers[0]
                    
                    
                    x_bin_size = in_file[h_name].to_boost().axes.widths[0] / 2
                    ax.text(0.40, 0.3, beam_name[beam] + r": $ {} <  \eta_{{cms}} < {}$".format(_eta_bins[s], _eta_bins[s+1]), transform=ax.transAxes, fontsize=40, verticalalignment='top')
                    if w !="0":
                        ax.errorbar(x=h_centers, y=h_values, xerr=x_bin_size, yerr=h_errors, marker=markers[index], color=colors[index], markersize=15, linestyle='none', label="Weight set :"+w)
                    else:
                        ax.errorbar(x=h_centers, y=h_values, xerr=x_bin_size, yerr=h_errors, marker=markers[index], color=colors[index], markersize=15, linestyle='none', label="Unweighted")

                    ax.set_xscale("log")
                    ax.set_ylim(0.2,0.95)
                    ax.set_xlim(400, 4619)
                    ax.grid(True)
                    ax.axhline(y=1.0, color='Black', linestyle='--')
                    ax.tick_params(axis='x', labelsize=50)  
                    ax.tick_params(axis='y', labelsize=50) 
                    #leg = ax.legend()
                    s += 1
            if beam == "Pbp":
                i, j = position["Pbp"]["{}_{}".format(_eta_bins[0], _eta_bins[1])]
                ax = axarr[i][j]
            else:
                ax= axarr[0][0]
            #ax.text(0.6, 0.96,"Sim : "+data1, transform=ax.transAxes, fontsize=20, verticalalignment='top')
            #ax.text(0.6, 0.92,"Pol : "+pol, transform=ax.transAxes, fontsize=20, verticalalignment='top')                    
            leg = ax.legend(loc='upper left',fontsize=22)
            title_text = "Sim: {}, Pol: {}".format(data1, pol)
            fig.suptitle(title_text, fontsize=50)
       
        plt.savefig(git_out_path+"/efficiencies/reco_effi/Plots/reco_"+particle+"_all_weights_beam_"+beam+"_"+data1+".pdf")
       
    """

beam_list=beam
ple_list=particles
pol_list=pol




print("MC",MC)
for pol in  pol_list:
    for i,beam in  enumerate(beam_list):  
        for w in weights:
            plot_efficiency(beam=beam,pol=pol,data1=MC[i],particle=particles,track_type="prompt",weights=w,data2="Gen"+MC[i],binning_mode="oscar_ana",ratio="reco_effi",ylabel=ylabes_reco_effi_beam[beam])

"""
for pol in  pol_list:
    for ple in particles:  
        for w in weights:
            plot_efficiency(beam=beam,pol=pol,data1=MC,particle=ple,track_type="prompt",weights=w,data2="Gen"+MC,binning_mode="oscar_ana",ratio="reco_effi",ylabel=ylabes_reco_effi_ple[ple])

for pol in  pol_list:
    for i,beam in  enumerate(beam_list):  
        for particles in ple_list:
            plot_efficiency(beam=beam[i],pol=pol,data1=MC[i],particle=particles,track_type="prompt",weights=weights,data2="Gen"+MC,binning_mode="oscar_ana",ratio="reco_effi",ylabel=ylabes_reco_effi_ple[particles])

"""