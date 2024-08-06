import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import mplhep as hep

# Define a new range for the Delta log L values centered around zero for the difference
delta_log_L_difference = np.linspace(-3, 3, 1000)

# Assume a single distribution for the difference in log-likelihood
# This could represent, for example, the distribution of differences where positive values favor particle X
# and negative values favor pions
difference_dist = norm.pdf(delta_log_L_difference, loc=0, scale=1)

# Define a cut point for demonstration
cut_point_for_difference = 0.5  # Arbitrary choice, positive favoring particle X


#plt.style.use(hep.style.LHCb2)
# Plotting the difference in log-likelihood distribution
plt.figure(figsize=(10, 6))
plt.plot(delta_log_L_difference, difference_dist, label=r'$\Delta \log \mathcal{L}_{X-\pi}$ Distribution', color='purple')
plt.axvline(cut_point_for_difference, color='orange', linestyle='--', label='Cut Point for Particle X')

plt.fill_between(delta_log_L_difference, difference_dist, where=(delta_log_L_difference < cut_point_for_difference), color='purple', alpha=0.2, label='Region favoring Pions')
plt.fill_between(delta_log_L_difference, difference_dist, where=(delta_log_L_difference > cut_point_for_difference), color='orange', alpha=0.2, label='Region favoring Particle X')

plt.xlabel(r'$\Delta \log \mathcal{L}_{X-\pi}$')
plt.ylabel('Probability Density')
plt.title(r'Distribution of $\Delta \log \mathcal{L}_{X-\pi}$ with Cut Point')
plt.legend()
plt.savefig("delta_likelyhood.pdf")

