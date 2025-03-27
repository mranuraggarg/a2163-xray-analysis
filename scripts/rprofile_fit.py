import numpy as np
import matplotlib.pyplot as plt
import emcee
from pycrates import read_file
from scipy.optimize import minimize
from scipy.stats import chi2 as chi2_dist

# Load the radial profile data
tab = read_file("radial_profile.fits")
xx = tab.get_column("rmid").values
yy = tab.get_column("sur_bri").values
ye = tab.get_column("sur_bri_err").values

# Filter NaNs from the input arrays
mask = ~np.isnan(xx) & ~np.isnan(yy) & ~np.isnan(ye)
r = xx[mask]
y = yy[mask]
yerr = ye[mask]

# Define the double beta model
def double_beta_model(r, r01, beta1, A1, r02, beta2, A2):
    sb1 = A1 * (1 + (r / r01)**2)**(-3 * beta1 + 0.5)
    sb2 = A2 * (1 + (r / r02)**2)**(-3 * beta2 + 0.5)
    return sb1 + sb2

# Log-likelihood
def lnlike(theta, r, y, yerr):
    r01, beta1, A1, r02, beta2, A2 = theta
    model = double_beta_model(r, r01, beta1, A1, r02, beta2, A2)
    inv_sigma2 = 1.0 / (yerr ** 2)
    return -0.5 * np.sum((y - model) ** 2 * inv_sigma2)

# Log-prior with physical constraints
def lnprior(theta):
    r01, beta1, A1, r02, beta2, A2 = theta
    if (10 < r01 < 500 and 0.4 < beta1 < 1.0 and 0 < A1 < 100 and
        r01 < r02 < 1500 and 0.6 < beta2 < 1.5 and 0 < A2 < 100):
        return 0.0
    return -np.inf

# Log-probability
def lnprob(theta, r, y, yerr):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, r, y, yerr)

# Initial guess
initial = np.array([50, 0.6, 1.0, 200, 1.0, 0.5])
result = minimize(lambda *args: -lnprob(*args), initial, args=(r, y, yerr))
best = result.x

# MCMC setup
ndim, nwalkers = 6, 60
pos = best + 1e-4 * np.random.randn(nwalkers, ndim)

sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(r, y, yerr))
sampler.run_mcmc(pos, 5000, progress=True)

samples = sampler.get_chain(discard=1000, thin=10, flat=True)
medians = np.percentile(samples, 50, axis=0)
errors = np.std(samples, axis=0)

param_names = ["r01", "beta1", "A1", "r02", "beta2", "A2"]
for i, name in enumerate(param_names):
    print(f"{name} = {medians[i]:.2f} ± {errors[i]:.2f}")

# Compute fit curves
rfit = np.linspace(min(r), max(r), 500)

# Extract individual components from medians
r01, beta1, A1, r02, beta2, A2 = medians
y_core = A1 * (1 + (rfit / r01) ** 2) ** (-3 * beta1 + 0.5)
y_halo = A2 * (1 + (rfit / r02) ** 2) ** (-3 * beta2 + 0.5)

# Plot best fit
yfit = double_beta_model(rfit, *medians)

# Compute predicted surface brightness from best-fit model
y_model = double_beta_model(r, *medians)

# Compute chi-square and reduced chi-square
chi2 = np.sum(((y - y_model) / yerr) ** 2)
dof = len(r) - len(medians)
reduced_chi2 = chi2 / dof
p_value = chi2_dist.sf(chi2, dof)

print(f"\nChi² = {chi2:.2f}")
print(f"DOF = {dof}")
print(f"Reduced Chi² = {reduced_chi2:.3f}")
print(f"P-value = {p_value:.3e}")

# === Posterior-based chi2 and p-value sampling ===
chi2_chain = []
pval_chain = []

for sample in samples:
    model_pred = double_beta_model(r, *sample)
    chi2_i = np.sum(((y - model_pred) / yerr) ** 2)
    dof_i = len(r) - len(sample)
    red_chi2_i = chi2_i / dof_i
    pval_i = chi2_dist.sf(chi2_i, dof_i)
    chi2_chain.append(red_chi2_i)
    pval_chain.append(pval_i)

chi2_chain = np.array(chi2_chain)
pval_chain = np.array(pval_chain)

print(f"\nPosterior Median Reduced Chi² = {np.median(chi2_chain):.3f} ± {np.std(chi2_chain):.3f}")
print(f"Posterior Median P-value = {np.median(pval_chain):.3f} ± {np.std(pval_chain):.3f}")

# === Add histogram plots to current figure ===
fig, axs = plt.subplots(3, 2, figsize=(12, 14))

# Top-left: Data
axs[0, 0].errorbar(r, y, yerr=yerr, fmt="o")
axs[0, 0].set_title("Surface Brightness Data")
axs[0, 0].set_xscale("log")
axs[0, 0].set_yscale("log")
axs[0, 0].set_xlabel("R_MID (pixel)")
axs[0, 0].set_ylabel("SUR_BRI")
axs[0, 0].grid(True)

# Top-right: Fit + components
axs[0, 1].errorbar(r, y, yerr=yerr, fmt="o", label="Data", alpha=0.5)
axs[0, 1].plot(rfit, yfit, label="Total Fit", color="black")
axs[0, 1].plot(rfit, y_core, label="Core", linestyle="--", color="blue")
axs[0, 1].plot(rfit, y_halo, label="Halo", linestyle="--", color="green")
axs[0, 1].set_title("Model Fit with Components")
axs[0, 1].set_xscale("log")
axs[0, 1].set_yscale("log")
axs[0, 1].set_xlabel("R_MID (pixel)")
axs[0, 1].set_ylabel("SUR_BRI")
axs[0, 1].legend()
axs[0, 1].grid(True)

# Middle-left: Residuals
residuals = (y - double_beta_model(r, *medians)) / yerr
axs[1, 0].axhline(0, color='gray', linestyle='--')
axs[1, 0].errorbar(r, residuals, yerr=np.ones_like(yerr), fmt='o', color="black")
axs[1, 0].set_title("Fit Residuals")
axs[1, 0].set_xlabel("R_MID (pixel)")
axs[1, 0].set_ylabel("Residuals (σ)")
axs[1, 0].grid(True)

# Middle-right: Posterior for beta1
axs[1, 1].hist(samples[:, 1], bins=30, density=True, alpha=0.7, color="skyblue")
axs[1, 1].set_title("Posterior: $\\beta_1$")
axs[1, 1].set_xlabel("$\\beta_1$")
axs[1, 1].set_ylabel("Density")

# Bottom-left: Reduced Chi² Distribution
axs[2, 0].hist(chi2_chain, bins=30, density=True, alpha=0.7, color="tomato")
axs[2, 0].set_title("Posterior: Reduced Chi²")
axs[2, 0].set_xlabel("Reduced Chi²")
axs[2, 0].set_ylabel("Density")

# Bottom-right: P-value Distribution
axs[2, 1].hist(pval_chain, bins=30, density=True, alpha=0.7, color="seagreen")
axs[2, 1].set_title("Posterior: P-value")
axs[2, 1].set_xlabel("P-value")
axs[2, 1].set_ylabel("Density")

plt.tight_layout()
plt.savefig("graphs/surface_brightness.png", dpi=600)
plt.show()