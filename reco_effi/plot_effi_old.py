import sys
sys.path.append("/home3/sara.sellam/RpPb_identified_hadrons_project")
import uproot
import matplotlib.pyplot as plt
import mplhep as hep
from Binning.Binning import *
from Helpers_function.import_dependencies import *

git_out_path="/scratch43/ssellam/results/"

def plot_beam_effi(beam="",pol="",data1="",particle=[],track_type="",weights="",data2="",binning_mode="",ratio="",xlabel="$p_{T}$[MeV/$c^{2}$]",ylabel=r"$\varepsilon_{reco}$",legend=""):
    _pt_bins=binning[binning_mode][beam]["pt"]
    _eta_bins=binning[binning_mode][beam]["eta"] 
    f, axarr = plot_2_3_subplot(ylabel=ylabel,xlabel="$p_{T}$[MeV/$c^{2}$]")
  
    colors,markers=listof_colors_markers()
    #particles=particle
    beam_name={"pp":r"$pp$","pPb":r"$pPb$","Pbp":r"$Pbp$"}
    particles_name={"pi":r"$\pi$","K":"K","p":r"$p$"}
    position={"Pbp":{"-3.0_-2.5":[0,0],
                     "-3.5_-3.0":[0,1],
                     "-4.0_-3.5":[0,2],
                     "-4.5_-4.0":[1,0],
                     "-4.8_-4.5":[1,1],
                     "-5.2_-4.8":[1,2]}}
    
    #print("particles",particles)
    for ple in particles:
        index = particles.index(ple)
        particle=ple

        print("particle",particle)
        in_file=uproot.open(git_out_path+"/efficiencies/reco_effi/rfiles/"+ratio+"_"+track_type+"_"+particle+"_"+beam+"_"+pol+"_w_"+weights+"_"+data1+"_"+data2+"_"+binning_mode+".root")
        s=0
        for i in  range(axarr.ndim):
            for j in range(3):
                if beam =="Pbp":
                    i,j= position["Pbp"]["{}_{}".format(_eta_bins[s],_eta_bins[s+1])]
                    ax=axarr[i][j]
                else :
                    ax=axarr[i][j]
                h_name=ratio+"_"+track_type+"_"+particle+"_{}_{}".format(_eta_bins[s],_eta_bins[s+1])
                #print(in_file.keys())
                h_errors=in_file[h_name].errors()
                h_values=in_file[h_name].to_numpy()[0]
                h_centers=in_file[h_name].to_boost().axes.centers[0]
                x_bin_size=in_file[h_name].to_boost().axes.widths[0]/2
                ax.text(0.05, 0.95,beam_name[beam]+": $ {}< \eta < {}$".format(_eta_bins[s],_eta_bins[s+1]), transform=ax.transAxes,fontsize=40,verticalalignment='top')
                ax.text(0.05, 0.85,data1, transform=ax.transAxes,fontsize=40,verticalalignment='top')
                print("h_values",h_values)
                ax.errorbar(x=h_centers,y=h_values,xerr=x_bin_size,yerr=h_errors,marker=markers[index],color=colors[index], markersize=15,linestyle='none',label=particles_name[particle])
                print("label",particles_name[particle])
                ax.set_xscale("log")
                #ax.set_yscale("log")
                ax.set_ylim(0.0001,1.1999)
                ax.set_xlim(400,4619)
                ax.grid(True)
                ax.axhline(y=1.0, color='Black', linestyle='--')
                leg = ax.legend()
                s=s+1
    plt.savefig(git_out_path+"/efficiencies/reco_effi/Plots/reco_"+beam+"_"+pol+"_w_"+weights+".pdf")
    
def plot_ple_effi(beam="",pol="",data1="",particle=[],track_type="",weights=[],data2="",binning_mode="",ratio="",xlabel="$p_{T}$[MeV/$c^{2}$]",ylabel=r"$\varepsilon_{reco}$",legend=""): 
    f, axarr = plot_2_3_subplot(ylabel=ylabel,xlabel="$p_{T}$[MeV/$c^{2}$]")
    colors,markers=listof_colors_markers()

    position={"Pbp":{"-3.0_-2.5":[0,0],
                     "-3.5_-3.0":[0,1],
                     "-4.0_-3.5":[0,2],
                     "-4.5_-4.0":[1,0],
                     "-4.8_-4.5":[1,1],
                     "-5.2_-4.8":[1,2]}}
    beam_name={"pp":r"$pp$","pPb":r"$pPb$","Pbp":r"$Pbp$"}

    if 1:
        for b in beam:
            index = beam.index(b)
            _eta_bins=binning[binning_mode][b]["eta"] 
            in_file=uproot.open(git_out_path+"/efficiencies/reco_effi/rfiles/"+ratio+"_"+track_type+"_"+particle+"_"+b+"_"+pol+"_w_"+weights[b]+"_"+data1+"_"+data2+"_"+binning_mode+".root")
            s=0
            for i in  range(axarr.ndim):
                for j in range(3):
                    if b =="Pbp":
                        i,j= position["Pbp"]["{}_{}".format(_eta_bins[s],_eta_bins[s+1])]
                        ax=axarr[i][j]
                    else:
                        ax=axarr[i][j]
                    h_name=ratio+"_"+track_type+"_"+particle+"_{}_{}".format(_eta_bins[s],_eta_bins[s+1])
                    h_errors=in_file[h_name].errors()
                    h_values=in_file[h_name].to_numpy()[0]
                    h_centers=in_file[h_name].to_boost().axes.centers[0]
                    x_bin_size=in_file[h_name].to_boost().axes.widths[0]/2
                    ax.text(0.05, 0.85,data1, transform=ax.transAxes,fontsize=40,verticalalignment='top')
                    ax.text(0.05, 0.95,beam_name[beam]+": $ {}< \eta < {}$".format(_eta_bins[s],_eta_bins[s+1]), transform=ax.transAxes,fontsize=40,verticalalignment='top')
                    ax.errorbar(x=h_centers,y=h_values,xerr=x_bin_size,yerr=h_errors,marker=markers[index],color=colors[index], markersize=15,linestyle='none',label=beam_name[b]+": $ {}< \eta < {}$".format(_eta_bins[s],_eta_bins[s+1]))
                    ax.set_xscale("log")
                    ax.set_ylim(0.0001,1.1999)
                    ax.set_xlim(400,4619)
                    ax.grid(True)
                    ax.axhline(y=1.0, color='Black', linestyle='--')
                    leg = ax.legend()
                    s=s+1
        plt.savefig(git_out_path+"/efficiencies/reco_effi/Plots/reco_"+particle+"_all_beams.pdf")



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
weights= args.weights[0]

MC= args.MC[0]

# Print the values of the arguments
print(f"Beam type: {beam}")
print(f"Polarization : {pol}")
print(f"Particles : {particles}")
print(f"binning_mode : {binning_mode}")
   


ylabes_beam_effi={"pp":r"$\varepsilon_{reco}$",
                "pPb":r"$\varepsilon_{reco}$",
                "Pbp":r"$\varepsilon_{reco}$"}

beam_name={"pp":r"$pp$","pPb":r"$pPb$","Pbp":r"$Pbp$"}


particles_name={"pi":r"$\pi$","K":r"$K$","p":r"$p$"}
ylabes_ple_effi={"pi":r"$\varepsilon_{reco}^{\pi}$",
                "K":r"$\varepsilon_{reco}^{K}$",
                "p":r"$\varepsilon_{reco}^{p}$"}

for b in beam:
    for p in pol: 
            plot_beam_effi(beam=b,pol=p,data1=MC,particle=particles,track_type="prompt",weights=weights,data2="Gen"+MC,binning_mode="oscar_ana",ratio="reco_effi",ylabel=ylabes_beam_effi[b])
               

for ple in particles:
    for p in pol: 
        plot_ple_effi(beam=["pPb","Pbp","pp"],pol=p,data1=MC,particle=ple,track_type="prompt",weights=weights,data2="Gen"+MC,binning_mode="oscar_ana",ratio="reco_effi",ylabel=ylabes_ple_effi[ple])
               
