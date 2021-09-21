# make MCMC PT profile figure 

import numpy as np
import os, sys
import argparse, ConfigParser
import scipy.interpolate as si
from scipy.ndimage.filters import gaussian_filter1d as gaussf
import matplotlib
import matplotlib.pyplot as plt
plt.ion()

sys.path.append("../BART/code")
import makeatm as mat
import reader as rd
import makeatm as mat
import PT as pt
import wine as w
import readtransit as rt
import bestFit as bf



############################ INPUTS ########################################

# 4 mol 4 specs
direct = './4species_4opac_uniform/'
cfile =  './BART_4species-4opac-uniform.cfg'
filename = './4species_4opac_uniform/bestFit.atm'

######################### END OF INPUTS ####################################

def parray(string):
  """
  Convert a string containin a list of white-space-separated (and/or
  newline-separated) values into a numpy array
  """
  if string == 'None':
    return None
  try:    # If they can be converted into doubles, do it:
    return np.asarray(string.split(), np.double)
  except: # Else, return a string array:
    return string.split()


def read_cfg(cfile):
  # Parse the config file from the command line:
  cparser = argparse.ArgumentParser(description=__doc__, add_help=False,
                         formatter_class=argparse.RawDescriptionHelpFormatter)
  # Add config file option:
  cparser.add_argument("-c", "--config_file",
                       help="Configuration file", metavar="FILE")

  # Parser for the MCMC arguments:
  parser = argparse.ArgumentParser(parents=[cparser])

  # Directories and files options:
  group = parser.add_argument_group("Directories and files")
  group.add_argument("--tep_name", dest="tep_name",
           help="Transiting exoplanet file name.",
           type=str, action="store", default=None)
  group.add_argument("--logfile", dest="logfile",
           help="MCMC log file [default: %(default)s]",
           type=str, action="store", default="MCMC.log")

  # Atmospheric model options:
  group.add_argument("--atmfile", dest="atmfile",
           help="Atmospheric model file [default: %(default)s]",
           type=str, action="store", default="")
  group.add_argument("--uniform", dest="uniform",
           help="If not None, set uniform abundances with the specified "
                "values for each species in out_spec [default: %(default)s]",
           type=parray, action="store", default=None)

  # MCMC options:
  group = parser.add_argument_group("MCMC")
  group.add_argument("--params",  dest="params",
           help="Model-fitting parameters [default: %(default)s]",
           type=parray, action="store", default=None)
  group.add_argument("--molfit",  dest="molfit", 
           help="Molecules fit [default: %(default)s]",
           type=parray, action="store", default=None)
  group.add_argument("--stepsize", dest="stepsize",
           help="Parameters stepsize",
           type=parray, action="store", default=None)
  group.add_argument("--burnin", dest="burnin",
           help="Number of burn-in iterations per chain",
           type=parray, action="store", default=None)
  group.add_argument("--data", dest="data",
           help="Transit or eclipse depths",
           type=parray, action="store", default=None)
  group.add_argument("--uncert", dest="uncert",
           help="Uncertanties on transit or eclipse depths",
           type=parray, action="store", default=None)

  # Output-Converter Options:
  group = parser.add_argument_group("Output Converter Options")
  group.add_argument("--filters",                 action="store",
           help="Waveband filter name [default: %(default)s]",
           dest="filters",   type=parray, default=None)
  group.add_argument("--kurucz_file",           action="store",
           help="Stellar Kurucz file [default: %(default)s]",
           dest="kurucz",   type=str,       default=None)
  group.add_argument("--solution",                    action="store",
           help="Solution geometry [default: %(default)s]",
           dest="solution", type=str,       default="None",
           choices=('transit', 'eclipse'))

  # Transit options:
  group.add_argument("--outspec", dest="outspec",
           help="Output with flux values [default: %(default)s]",
           type=str, action="store", default=None)
  group.add_argument("--outmod", dest="outmod",
           help="Output with modulation values [default: %(default)s]",
           type=str, action="store", default=None)

  # Remaining_argv contains all other command-line-arguments:
  cargs, remaining_argv = cparser.parse_known_args()
  # Get only the arguments defined above:
  known, unknown = parser.parse_known_args(remaining_argv)

  # Read values from configuration file:
  config = ConfigParser.SafeConfigParser()
  config.optionxform = str  # This one enable Uppercase in arguments
  config.read([cfile])
  defaults = dict(config.items("MCMC"))

  # Set the defaults from the configuration file:
  parser.set_defaults(**defaults)
  # Set values from command line:
  args, unknown = parser.parse_known_args(remaining_argv)

  # Unpack the variables from args:
  variables = dir(args)
  for var in dir(known):
    if not var.startswith("_"):
      exec("{:s} = args.{:s}".format(var, var))

  return stepsize, molfit, params, burnin, atmfile, logfile, filters, kurucz, tep_name, solution, outspec, data, uncert

# read config file
stepsize, molfit, params, burnin, atmfile, logfile, filters, kurucz, tep_name, solution, outspec, data, uncert= read_cfg(cfile)


def bestFit_abun(init_abun, species_names, filename):
    mol_fin, pres_fin, temp_fin, abun_fin = mat.readatm(filename)

    # Open the atmospheric file and read
    f = open(filename, 'r')
    lines = np.asarray(f.readlines())
    f.close()

    # Get molecules names
    imol = np.where(lines == "#SPECIES\n")[0][0] + 1
    molecules = lines[imol].split()
    nmol = len(molecules)
    for m in np.arange(nmol):
        molecules[m] = molecules[m].partition('_')[0]

    # Take user input for species and split species strings into separate strings 
    #      convert the list to tuple
    species_names = tuple(species_names.split(','))
    nspec = len(species_names)

    # Populate column numbers for requested species and 
    #          update list of species if order is not appropriate
    columns = []
    spec    = []
    for i in np.arange(nmol):
        for j in np.arange(nspec):
            if molecules[i] == species_names[j]:
                columns.append(i+3)   # defines p, T +2 or rad, p, T +3
                spec.append(species_names[j])

    # Convert spec to tuple
    spec = tuple(spec)

    # Concatenate spec with pressure for data and columns
    data    = tuple(np.concatenate((['p'], spec)))
    usecols = tuple(np.concatenate(([1], columns))) # defines p as 0 columns, or p as 1 columns

    # Load all data for all interested species
    data = np.loadtxt(filename, dtype=float, comments='#', delimiter=None,    \
                    converters=None, skiprows=13, usecols=usecols, unpack=True)

    return data


def multi_hist(init_abun, params, species, species_names, filename):
    # get MCMC data:
    MCMCdata = direct + "/output.npy"
    data = np.load(MCMCdata)
    nchains, npars, niter = np.shape(data)

    # this is data from first chain zero
    burnin = 100000
    data_stack = data[0,:,burnin:]
    bad_chain = 0
    bad_chain2 = 0
    # now he stack all other chains from 1 
    for c in np.arange(1, nchains):
        if c != bad_chain and c != bad_chain2:
            print c
        else:
            continue

    nrows =1
    ncolumns = len(params)
    plt.figure(100, figsize=(10, 2))
    plt.clf()

    data = bestFit_abun(init_abun, species_names, filename)
    best_abun = np.zeros(7)
    best_abun[0] = data[4][0]
    best_abun[1] = data[2][0]
    best_abun[2] = data[1][0]
    best_abun[3] = data[3][0]

    if len(init_abun)>4:
        best_abun[4] = data[5][0]
        best_abun[5] = data[6][0]
        best_abun[6] = data[7][0]


    nabun_params = len(init_abun)

    color = 'b', 'orangered', 'm', 'cadetblue', 'firebrick', 'darkgreen', 'palevioletred' 
    maxylim = 0  # Max Y limit
    for i in np.arange(npars-3):
        ax = plt.subplot(nrows, ncolumns, i+1)
        a  = plt.xticks(size=14, rotation=90)
        num = params[i]
        plt.xlabel(species[i], size=15)
        a = plt.hist(np.log10(10**data_stack[num]*init_abun[i]), 20, color= color[i], normed=False, alpha=0.5)
        # add best fit values as dashed vertical lines
        plt.axvline(np.log10(best_abun[i]), color= color[i], linestyle='dashed', linewidth=2)
        # sets uniform height
        maxylim = np.amax((maxylim, ax.get_ylim()[1]))
        ax.set_ylim(0, maxylim)
        if i==0:
            a = plt.yticks(size=14)
        else:
            a = plt.yticks(size=14, visible=False)  

    plt.subplots_adjust(bottom=0.50)

    plt.savefig('4mol-4specs-Histograms.png') 
    plt.savefig('4mol-4specs-Histograms.ps') 


# 4 species and their order
species_names = 'CO,CO2,CH4,H2O'
init_abun = np.array([5e-4, 1e-7, 3e-4, 1e-6])
params = [3, 4, 5, 6]
species = [r"${\rm H_{2}O}$", r"${\rm CO_{2}}$", r"${\rm CO}$", r"${\rm CH_{4}}$"]
multi_hist(init_abun, params, species, species_names, 
filename)






