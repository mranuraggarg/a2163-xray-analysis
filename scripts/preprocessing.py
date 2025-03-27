from ciao_contrib.runtool import chandra_repro, dmcopy

# Reset parameters
chandra_repro.punlearn()

# Set inputs
chandra_repro.indir = "../../1653"
chandra_repro.outdir = "data"
chandra_repro.clobber = "yes"
chandra_repro.cleanup = "no"
chandra_repro.verbose = 1

# Run chandra_repro
# chandra_repro()

# Resent paramerters
dmcopy.punlearn()

# Set inputs
dmcopy.infile = "data/acisf01653_repro_evt2.fits[energy=300:8000]"
dmcopy.outfile = "acis_1653_evt2.fits"
dmcopy.clobber = "yes"
dmcopy.verbose = 1

# Run dmcopy
dmcopy()