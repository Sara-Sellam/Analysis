import sys
sys.path.append("/home3/sara.sellam/RpPb_identified_hadrons_project")
import uproot as up
import matplotlib.pyplot as plt
import mplhep as hep
from Binning.Binning import *
from Helpers_function.import_dependencies import *
from Helpers_function.labels import * 
plt.style.use(hep.style.LHCb2)

beam="pPb"

in_file=up.open("/scratch43/ssellam/results/efficiencies/sel_effi/rfiles/sel_effi_"+beam+"_Down_w_1_MC-Sim09e_sara_ana_num_denom.root")

fig, axarr = plot_2_3_subplot(ylabel="effi", xlabel=pt_label)
colors, markers = listof_colors_markers()
position={"Pbp":{"-3.0_-2.5":[0,0],
                     "-3.5_-3.0":[0,1],
                     "-4.0_-3.5":[0,2],
                     "-4.5_-4.0":[1,0],
                     "-4.8_-4.5":[1,1],
                     "-5.2_-4.8":[1,2]}}
binning_mode="sara_ana"
_pt_bins = binning[binning_mode][beam]["pt"]
_eta_bins = binning[binning_mode][beam]["eta"] 

s = 0
for i in range(axarr.ndim):
    for j in range(3):
        if beam == "Pbp":
            i, j = position["Pbp"]["{}_{}".format(_eta_bins[s], _eta_bins[s+1])]
            ax = axarr[i][j]
        else:
            ax = axarr[i][j]
        h_errors = in_file["sel_effi"].to_boost().variances()[:,s]
        h_values = in_file["sel_effi"].to_boost().counts()[:,s]
        pt_bins= np.asarray(_pt_bins)
        h_centers = (pt_bins[:-1] + pt_bins[1:]) / 2
        h_edges = pt_bins[1:] - pt_bins[:-1]
        ax.text(0.40, 0.3, beam_name[beam] + r": $ {} <  \eta_{{cms}} < {}$".format(_eta_bins[s], _eta_bins[s+1]), transform=ax.transAxes, fontsize=35, verticalalignment='top')
        ax.errorbar(x=h_centers, y=h_values, xerr=h_edges*0.5, yerr=h_errors, marker=".", color="blue", markersize=16, linestyle='none', label="all")
        ax.set_xscale("log")
        ax.set_ylim(0.2,1.2)
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
leg = ax.legend(loc='upper left',fontsize=35)

plt.savefig("/scratch43/ssellam/results/efficiencies/sel_effi/Plots/sel_effi_"+beam+"_Down_w_1_MC-Sim09e_sara_ana_num_denom.pdf")
