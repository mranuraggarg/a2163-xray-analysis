from ciao_contrib.runtool import dmextract
from pycrates import read_file
import matplotlib.pylab as plt

# Reset any previous parameter settings
dmextract.punlearn()

# Set parameters as per your shell commands
dmextract.infile = "acis_1653_evt2.fits[bin sky=@annuli.reg]"
dmextract.outfile = "radial_profile.fits"
dmextract.bkg = "acis_1653_evt2.fits[bin sky=@annuli_bkg.reg]"
# dmextract.bkg = "blanksky_reproj.fits[bin sky=@annuli_bkg.reg]"
dmextract.opt = "generic"
dmextract.clobber = "yes"
dmextract.verbose = 1  # Must be integer, not string!

# Run the tool
dmextract()

# Read radial profile
tab = read_file("radial_profile.fits")
xx = tab.get_column("rmid").values
yy = tab.get_column("sur_bri").values
ye = tab.get_column("sur_bri_err").values

# Setup plotting and plot radial profile
plt.figure(figsize=(8, 6))  # Width = 8 inches, Height = 6 inches

plt.errorbar(xx, yy, yerr=ye, marker="o")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("R_MID (pixel)")
plt.ylabel("SUR_BRI (photons/cm**2/pixel**2/s)")
plt.title("G21.5-0.9 [Chip S3, T=120, Offsets=-1,0,1]")

plt.tight_layout()  # Optional: adjusts spacing to fit labels
plt.show()