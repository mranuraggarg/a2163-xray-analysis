from sherpa.astro import ui
import matplotlib.pyplot as plt

ui.set_stat("chi2") 

# Load radial profile data: columns 3 = RMID, SUR_BRI, SUR_BRI_ERR
ui.load_data(1, "radial_profile.fits", 3, ["RMID", "SUR_BRI", "SUR_BRI_ERR"])

# Define a double beta model as the source
ui.set_source("beta1d.core + beta1d.halo")

# Set initial parameter values for core component
core.xpos = 0
core.r0 = 50.0
core.beta = 0.7
core.ampl = 1.0

# Set initial parameter values for halo component
halo.xpos = 0
halo.r0 = 200.0
halo.beta = 1.0
halo.ampl = 0.5

# Freeze the center (assumed same for both components)
core.xpos.freeze()
halo.xpos.freeze()

# Set bounds to improve fitting stability
core.r0.min = 1e-1
core.r0.max = 1e4
core.beta.min = 0.5
core.beta.max = 10.0
core.ampl.min = 1e-6
core.ampl.max = 1e3

halo.r0.min = 1e-1
halo.r0.max = 1e4
halo.beta.min = 0.5
halo.beta.max = 10.0
halo.ampl.min = 1e-6
halo.ampl.max = 1e3

# Fit the model to the data
ui.fit()
ui.covar()

# # Run MCMC sampling to explore parameter space
# stats, accept, params = ui.get_draws(niter=10000)

# # Extract beta samples (beta 1 for src.beta), discarding burn-in
# import numpy as np
# burn = 1000
# beta_samples = params[1, burn:]  # param beta 1 corresponds to src.beta

# # Print posterior statistics
# print("Posterior mean beta:", np.mean(beta_samples))
# print("Posterior std dev beta:", np.std(beta_samples))

# # Plot the fit with logarithmic axes
# ui.plot_fit(xlog=True, ylog=True)
ui.plot_fit()
plt.xlabel("R_MID (pixel)")
plt.ylabel("Surface Brightness (photons/cm²/pixel²/s)")
plt.title("Beta Model Fit to Radial Profile")
plt.grid(True)
plt.tight_layout()
plt.show()
