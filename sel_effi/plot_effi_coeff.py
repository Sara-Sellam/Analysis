
import sys
sys.path.append("/home3/sara.sellam/RpPb_identified_hadrons_project")
import json
import json
import matplotlib.pyplot as plt
import numpy as np
import mplhep as hep
from Helpers_function.labels import * 

plt.style.use(hep.style.LHCb2)

beam = "pPb"

# Load JSON data from file
with open("/home3/sara.sellam/RpPb_identified_hadrons_project/efficiencies/sel_effi/sel_effi_coeffi_" + beam + "_ana_oscar_cuts.json", 'r') as f:
    data = json.load(f)

etas_str = list(data.keys())
_eta = [_eta_str.split('_') for _eta_str in etas_str]
etas = [float(item) for sublist in _eta for item in sublist]
etas = list(dict.fromkeys(etas))

bins_str = list(data[etas_str[0]].keys())
_bins = [bin_str.split('_') for bin_str in bins_str]
bins = [int(item) for sublist in _bins for item in sublist]
bins = list(dict.fromkeys(bins))

# Prepare data for plotting
x_bins = [bin for bin in bins]
indices = np.arange(len(bins))

fig, ax = plt.subplots(figsize=(12, 8))

for i, eta in enumerate(etas_str):
    values = [data[eta][bin_str]["sel_effi_coeffi"]["value"] for bin_str in bins_str]
    errors = [data[eta][bin_str]["sel_effi_coeffi"]["error"] for bin_str in bins_str]
    if eta=="2_3":
        values = np.array(values)*np.array([0,1])
        errors = np.array(errors)*np.array([0,1])
    bin_centers = [(bins[j] + bins[j+1]) / 2 for j in range(len(bins) - 1)]
    bin_widths = [(bins[j+1] - bins[j]) / 2 for j in range(len(bins) - 1)]

    color = plt.cm.plasma(i / len(etas_str))  # Assign color based on eta index
    ax.errorbar(bin_centers, values, xerr=bin_widths, yerr=errors, fmt='o', label=r'$ {} <  \eta < {}$'.format(etas[i], etas[i+1]), color=color)

ax.axhline(y=1, color='black', linestyle='--')
ax.set_xlabel(pt_label, fontsize=50)
ax.set_ylabel(r"$C_{\mathrm{sel}}$", fontsize=50)
ax.set_ylim(0.88, 1.02)
ax.set_xlim(bins[0], bins[-1])
ax.grid(True)
ax.tick_params(axis='x', labelsize=30)
ax.tick_params(axis='y', labelsize=30)
ax.legend()
# Add legend outside the plot
#handles, labels = ax.get_legend_handles_labels()
#ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5), fontsize=12)

plt.tight_layout()
plt.savefig("/home3/sara.sellam/RpPb_identified_hadrons_project/efficiencies/sel_effi/sel_effi_coeffi_"+beam+"_ana_sara_cuts.pdf")

""""
import sys
sys.path.append("/home3/sara.sellam/RpPb_identified_hadrons_project")
import json
import matplotlib.pyplot as plt
import numpy as np
from Helpers_function.labels import * 
import mplhep as hep
plt.style.use(hep.style.LHCb2)


beam="pPb"
# Load JSON data from file
with open("/home3/sara.sellam/RpPb_identified_hadrons_project/efficiencies/sel_effi/sel_effi_coeffi_"+beam+"_ana_oscar_cuts.json", 'r') as f:
    data = json.load(f)



etas_str = list(data.keys())
_eta=[_eta_str.split('_')for _eta_str in etas_str]
etas = [float(item) for sublist in _eta for item in sublist]
etas=list(dict.fromkeys(etas))


bins_str = list(data[etas_str[0]].keys())  
_bins = [bin_str.split('_')for bin_str in bins_str]
bins = [int(item) for sublist in _bins for item in sublist]
bins=list(dict.fromkeys(bins))
print("bins",bins)


# Prepare data for plotting
x_bins = [bin for bin in bins]
indices = np.arange(len(bins))

fig, axs = plt.subplots(1, len(etas_str), figsize=(20, 8), sharey=True)
fig.supylabel(r"$C_{sel}$",x=0,y=0.55,fontsize=50)
for i, eta in enumerate(etas_str):
    ax = axs[i]
    values = [data[eta][bin_str]["sel_effi_coeffi"]["value"] for bin_str in bins_str]
    errors = [data[eta][bin_str]["sel_effi_coeffi"]["error"] for bin_str in bins_str]
    
    bin_centers = [(bins[j] + bins[j+1]) / 2 for j in range(len(bins) - 1)]
    bin_widths = [(bins[j+1] - bins[j])/2 for j in range(len(bins)-1)]

    print("bin_centers",bin_centers)
    
    print("errors",errors)
    print("values",values)
    ax.text(0.10, 0.9, beam_name[beam] + r": $ {} <  \eta < {}$".format(etas[i], etas[i+1]), transform=ax.transAxes, fontsize=35, verticalalignment='top')
    ax.errorbar(bin_centers, values,xerr=bin_widths, yerr=errors, fmt='o', label='Coefficient Value')

    ax.axhline(y=1, color='black', linestyle='--') 
    ax.set_xlabel(pt_label,fontsize=50)
    #ax.set_xscale("log")
    ax.set_ylim(0.8,1.1)
    ax.set_xlim(bins[0], bins[-1])
    ax.grid(True)
    ax.tick_params(axis='x', labelsize=50)  
    ax.tick_params(axis='y', labelsize=50) 
   
    #ax.legend()


plt.tight_layout()
plt.savefig("/home3/sara.sellam/RpPb_identified_hadrons_project/efficiencies/sel_effi/sel_effi_coeffi_pPb_ana_sara_cuts.pdf")
"""