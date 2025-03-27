import numpy as np
from scipy.special import gamma
from scipy.integrate import quad
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme(style="ticks", palette="deep")

# Reading the beta fitting parameters from the file fit_results
# Read and parse the fit result file safely
fit_params = {}
with open("results/fit_result") as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 5 and parts[1] == '=':
            try:
                key = parts[0]
                value = float(parts[2])
                err = float(parts[4])
                fit_params[key] = np.array([value, err])
            except ValueError:
                continue  # skip lines with invalid float conversion

## Adding spectral fit parameters to fit_params directory
fit_params['norm1'] = np.array([0.0017175, 5.88e-5])
fit_params['Abundanc'] = np.array([0.355905, 0.125528])
fit_params['kT'] = np.array([12.4685, 1.97064])

# Cosmological and geometric constants
pix_size_arcsec = 0.492            # arcsec per pixel
arcsec_to_kpc = 3.387              # conversion at z=0.203
angular_scale = pix_size_arcsec * arcsec_to_kpc

# Convert core radius from pixels to kpc
R0_kpc = fit_params['r01'] * angular_scale  # R0 in kpc
Physical_parameters = {
        'R_0': R0_kpc
}

print(f"Core radius R₀ = {R0_kpc[0]:.2f} ± {R0_kpc[1]:.2f} kpc")

# Calculate central electron density n_e0
# Using XSPEC normalization: 
# norm = 10^-14 / (4π DA^2 (1+z)^2) × ∫ n_e n_H dV
# For a β-model and assuming n_e ≈ 1.2 n_H, we solve for n_e0:

# Input values
kpc_to_cm = 3.0857e21  # 1 kpc = 3.0857e21 cm
z = 0.203
Da_kpc = 683000  # Angular diameter distance in kpc for z=0.203
Da_cm = Da_kpc * kpc_to_cm  # Convert distances
rc_cm = R0_kpc * kpc_to_cm  # core radius in cm

# Load parameters with uncertainties
norm_val = fit_params['norm1']        # [value, uncertainty]
beta = fit_params['beta1']            # [value, uncertainty]

# Compute numerator of n_e0 expression
numerator = norm_val[0] * 4 * np.pi * (Da_cm * (1 + z))**2 / 1e-14

# Compute denominator: includes beta model volume integral
gamma_num = gamma(3 * beta[0] - 1.5)
gamma_den = gamma(3 * beta[0])
denominator = 1.2 * np.pi**1.5 * rc_cm[0]**3 * gamma_num / gamma_den

# Estimate fractional errors
d_norm = norm_val[1] / norm_val[0]
d_rc = 3 * (rc_cm[1] / rc_cm[0])  # due to rc^3 term

# Numerical derivative approximation for gamma ratio w.r.t. beta
delta = 1e-5
gamma_up_num = gamma(3 * beta[0] - 1.5 + delta)
gamma_up_den = gamma(3 * beta[0] + delta)
gamma_ratio = gamma_num / gamma_den
d_gamma_ratio = (
    3 * gamma_num * gamma_den**-2 *
    (gamma_den * gamma_up_num - gamma_num * gamma_up_den)
) / gamma_ratio
d_beta = abs(d_gamma_ratio) * beta[1] / beta[0]

# Combine denominator errors
d_den = np.sqrt(d_rc**2 + d_beta**2)

# Total fractional uncertainty in n_e0
d_total = np.sqrt(d_norm**2 + d_den**2)

# Final n_e0 value and uncertainty
ne0_val = np.sqrt(numerator / denominator)
ne0_err = 0.5 * ne0_val * d_total
ne0 = np.array([ne0_val, ne0_err])

Physical_parameters['n_e0'] = ne0

print(f"Central electron density nₑ₀ ≈ {ne0[0]:.2e} ± {ne0[1]:.2e} cm⁻³")

# Compute the gas mass profile
# Define constants
mp = 1.6726e-24  # proton mass in grams
mu_e = 1.14      # mean molecular weight per electron
r_max_kpc = 1000  # maximum radius to integrate to in kpc
r0_cm = rc_cm[0]
beta_val = beta[0]

# Define electron density profile function
def ne_r(r):
    return ne0_val * (1 + (r / r0_cm) ** 2) ** (-1.5 * beta_val)

# Define integrand for gas mass profile
def integrand(r):
    return 4 * np.pi * r**2 * ne_r(r)

# Integrate from 0 to r_max
r_max_cm = r_max_kpc * kpc_to_cm
gas_mass, _ = quad(integrand, 0, r_max_cm)

# Convert to solar masses
M_gas = mu_e * mp * gas_mass / 1.989e33  # in solar masses

Physical_parameters['M_gas(<1000 kpc)'] = M_gas

from numpy.random import normal

# Monte Carlo sampling for gas mass error estimation
N_samples = 1000
ne0_samples = normal(loc=ne0_val, scale=ne0_err, size=N_samples)

gas_mass_samples = []
for ne0_s in ne0_samples:
    def ne_r_sampled(r):
        return ne0_s * (1 + (r / r0_cm) ** 2) ** (-1.5 * beta_val)
    def integrand_sampled(r):
        return 4 * np.pi * r**2 * ne_r_sampled(r)
    mass_s, _ = quad(integrand_sampled, 0, r_max_cm)
    M_s = mu_e * mp * mass_s / 1.989e33
    gas_mass_samples.append(M_s)

M_gas_err = np.std(gas_mass_samples)
Physical_parameters['M_gas_error'] = M_gas_err
print(f"Gas mass within 1000 kpc: {M_gas:.2e} ± {M_gas_err:.2e} M_sun")

# Define radii and pressure profile
radii_kpc = np.linspace(0.1, r_max_cm / kpc_to_cm, 30)
radii_cm = radii_kpc * kpc_to_cm
pe_profile = []
pe_profile_err = []
for r_cm in radii_cm:
    ne_samples = normal(loc=ne0_val, scale=ne0_err, size=1000)
    ne_at_r = ne_samples * (1 + (r_cm / r0_cm) ** 2) ** (-1.5 * beta_val)
    pe_at_r = ne_at_r * fit_params["kT"][0]  # Use kT from fit_params
    pe_profile.append(np.mean(pe_at_r))
    pe_profile_err.append(np.std(pe_at_r))
pe_profile = np.array(pe_profile)
pe_profile_err = np.array(pe_profile_err)

# Plot and save the electron pressure profile using improved Seaborn aesthetics
plt.figure(figsize=(10, 7))

# Use Seaborn's default color palette
color_palette = sns.color_palette("colorblind")

# Plot electron pressure with error bars
plt.errorbar(
    radii_kpc, pe_profile, yerr=pe_profile_err,
    fmt='o', capsize=4, label="Electron Pressure (data)",
    color=color_palette[0], alpha=0.8
)

# Fit curve (e.g., theoretical model) over a smooth range
r_fit = np.linspace(0.1, r_max_kpc, 300)
pe_fit = ne0_val * (1 + (r_fit * kpc_to_cm / r0_cm)**2)**(-1.5 * beta_val) * fit_params["kT"][0]

plt.plot(
    r_fit, pe_fit,
    label="Model Fit",
    color=color_palette[2],
    linewidth=2.5
)

plt.xlabel("Radius (kpc)", fontsize=14)
plt.ylabel(r"Pressure $P_e$ (keV cm$^{-3}$)", fontsize=14)
plt.title("Electron Pressure Profile", fontsize=16)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.grid(True, linestyle="--", alpha=0.5)
plt.legend(frameon=False, fontsize=12)
sns.despine()
plt.tight_layout()
plt.savefig("graphs/electron_pressure_profile.png", dpi=600)
plt.close()

# -----------------------------
# Hydrostatic Mass Profile
# -----------------------------
from scipy.constants import G, m_p, k as k_B

# Constants
mu = 0.61  # mean molecular weight (fully ionized plasma)
T_keV = fit_params['kT'][0]
T_erg = T_keV * 1.60218e-9  # convert keV to erg
G_cgs = G * 1e3  # convert m^3/kg/s^2 to cm^3/g/s^2

# Compute hydrostatic mass using isothermal beta-model
def hydrostatic_mass(r_cm):
    return (3 * beta_val * T_erg * r_cm**3) / (G_cgs * mu * m_p * (r_cm**2 + r0_cm**2))

# Radii array (in cm)
radii_cm_mass = np.linspace(10 * kpc_to_cm, r_max_kpc * kpc_to_cm, 200)
radii_kpc_mass = radii_cm_mass / kpc_to_cm

# Calculate mass profile and convert to solar mass
mass_profile = hydrostatic_mass(radii_cm_mass) / 1.989e33

# Store mass profile for further use
Physical_parameters['M_hydrostatic_profile'] = {
    'r_kpc': radii_kpc_mass,
    'M_solar': mass_profile
}

# Plot hydrostatic mass profile using Seaborn
plt.figure(figsize=(10, 6))
sns.lineplot(x=radii_kpc_mass, y=mass_profile, color="crimson", linewidth=2.5)
plt.xlabel("Radius (kpc)", fontsize=14)
plt.ylabel(r"Hydrostatic Mass $M(<r)$ ($M_\odot$)", fontsize=14)
plt.title("Hydrostatic Mass Profile (Isothermal β-Model)", fontsize=16)
plt.grid(True, linestyle="--", alpha=0.6)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
sns.despine()
plt.tight_layout()
plt.savefig("graphs/hydrostatic_mass_profile.png", dpi=600)
plt.close()

import json

# Convert numpy arrays to regular lists for JSON serialization
def convert_for_json(d):
    def convert_value(val):
        if isinstance(val, np.ndarray):
            return val.tolist()
        elif isinstance(val, dict):
            return convert_for_json(val)
        return val
    return {k: convert_value(v) for k, v in d.items()}

with open("../SZ_effect/fit_params.json", "w") as f:
    json.dump(convert_for_json(fit_params), f, indent=4)

with open("../SZ_effect/physical_parameters.json", "w") as f:
    json.dump(convert_for_json(Physical_parameters), f, indent=4)