# when running in ipython session MUST BE in the BART/code/ directory


import numpy as np
import os, sys
import argparse, ConfigParser
import scipy.interpolate as si
from scipy.ndimage.filters import gaussian_filter1d as gaussf
import matplotlib
import matplotlib.pyplot as plt
plt.ion()

sys.path.append("../BART/code")
import reader as rd
import makeatm as mat
import PT as pt
import wine as w
import readtransit as rt
import bestFit as bf

filters = ['../WASP43b_filt/spitzer_irac1_sa.dat',
           '../WASP43b_filt/spitzer_irac2_sa.dat',
           '../WASP43b_filt/spitzer_irac1_sa.dat',
           '../WASP43b_filt/spitzer_irac2_sa.dat',
           '../WASP43b_filt/Wang-Hband.dat',
           '../WASP43b_filt/Wang-Kband.dat',
           '../WASP43b_filt/VLT_1190.dat',
           '../WASP43b_filt/VLT_2090.dat',
           '../WASP43b_filt/GROND_K_JB.dat',
           '../WASP43b_filt/GROND_i_JB.dat',
           '../WASP43b_filt/Zhou_Ks.dat',
           '../WASP43b_filt/grism141-filter0.dat',
           '../WASP43b_filt/grism141-filter1.dat',
           '../WASP43b_filt/grism141-filter2.dat',
           '../WASP43b_filt/grism141-filter3.dat',
           '../WASP43b_filt/grism141-filter4.dat',
           '../WASP43b_filt/grism141-filter5.dat',
           '../WASP43b_filt/grism141-filter6.dat',
           '../WASP43b_filt/grism141-filter7.dat',
           '../WASP43b_filt/grism141-filter8.dat',
           '../WASP43b_filt/grism141-filter9.dat',
           '../WASP43b_filt/grism141-filter10.dat',
           '../WASP43b_filt/grism141-filter11.dat',
           '../WASP43b_filt/grism141-filter12.dat',
           '../WASP43b_filt/grism141-filter13.dat',
           '../WASP43b_filt/grism141-filter14.dat']

tep_name = '../tepfile/WASP-43b.tep'

kurucz = '../kurucz/WASP43b-fp00ak2odfnew.pck'

data = np.array([ 0.00347 ,  0.00382 ,  0.0033  ,  0.003827,  0.00103 ,  0.00194 ,
        0.00079 ,  0.00156 ,  0.00197 ,  0.00037 ,  0.00181 ,  0.000367,
        0.000431,  0.000414,  0.000482,  0.00046 ,  0.000473,  0.000353,
        0.000313,  0.00032 ,  0.000394,  0.000439,  0.000458,  0.000595,
        0.000614,  0.000732])

uncert = np.array([  1.30000000e-04,   1.50000000e-04,   8.90000000e-05,
         8.40000000e-05,   1.70000000e-04,   2.90000000e-04,
         3.20000000e-04,   1.40000000e-04,   4.20000000e-04,
         2.20000000e-04,   2.70000000e-04,   4.50000000e-05,
         3.90000000e-05,   3.80000000e-05,   3.60000000e-05,
         3.70000000e-05,   3.30000000e-05,   3.40000000e-05,
         3.00000000e-05,   3.60000000e-05,   3.60000000e-05,
         3.30000000e-05,   3.50000000e-05,   3.60000000e-05,
         3.70000000e-05,   4.20000000e-05])

burnin = 100000
molfit = ['H2O', 'CO2', 'CO', 'CH4', 'NH3', 'HCN', 'C2H2']

plt.figure(1, figsize=(4, 5))
plt.ylabel("Pressure  (bar)", fontsize=14)
plt.xlabel("Temperature  (K)", fontsize=14)

plt.rcParams["mathtext.default"] = 'rm'
matplotlib.rcParams.update({'mathtext.default':'rm'})
matplotlib.rcParams.update({'axes.labelsize': 16,
                                'xtick.labelsize': 14,
                                'ytick.labelsize': 14,})

# 4 mol uniform
atmfile ='./7species_7opac_uniform/bestFit.atm'
params= np.array([-0.6 , -0.4 ,  0.0  ,  0.0  ,  1.09,  -1.07,   -1.168,   1.784,  -1.749,   -3.0,     3.0,    -4.5])
logfile='./7species_7opac_uniform/MCMC.log'
stepsize= np.array([ 0.01,  0.01,  0.0  ,  0.0,    0.01,     0.01,    0.01,    0.01,    0.01,     0.01,    0.01,    0.01])
MCMCdata= './7species_7opac_uniform/output.npy'


# read atmfile
molecules, pressure, temp, abundances = mat.readatm(atmfile)
# get surface gravity
grav, Rp = mat.get_g(tep_name)
# convert gravity to cm/s^2
grav = grav*1e2
# get star data
R_star, T_star, sma, gstar = bf.get_starData(tep_name)
# get best parameters
bestP, uncer = bf.read_MCMC_out(logfile)
# get all params
allParams = bf.get_params(bestP, stepsize, params)
# get PTparams and abundances factors
nparams = len(allParams)
nmol = len(molfit)
nPTparams = nparams - nmol
PTparams  = allParams[:nPTparams]
kappa,  gamma1, gamma2, alpha, beta = PTparams
# HARDCODED !
T_int = 100 # K
T_int_type = 'const'
# call PT line profile to calculate temperature
best_T = pt.PT_line(pressure, kappa,  gamma1, gamma2, alpha, beta, R_star, T_star, T_int, sma, grav, T_int_type)

# get MCMC data:
data = np.load(MCMCdata)
nchains, npars, niter = np.shape(data)

# this is data fro first chain zero
data_stack = data[0,:,burnin:]
# pick a chain that you do not like
bad_chain = 0
bad_chain2 = 0
# now he stack all other chains from 1, you can start from 5 here 
for c in np.arange(1, nchains):
    if c != bad_chain and c!= bad_chain2:
        data_stack = np.hstack((data_stack, data[c, :, burnin:]))
        print c
    else:
        continue

# create array of PT profiles
PTprofiles = np.zeros((np.shape(data_stack)[1], len(pressure)))
# current PT parameters for each chain, iteration
curr_PTparams = PTparams
# fill-in PT profiles array
for k in np.arange(0, np.shape(data_stack)[1]):
    j = 0
    for i in np.arange(len(PTparams)):
        if stepsize[i] != 0.0:
            curr_PTparams[i] = data_stack[j,k]
            j +=1
        else:
            pass
    kappa,  gamma1, gamma2, alpha, beta = curr_PTparams
    PTprofiles[k] = pt.PT_line(pressure, kappa,  gamma1, gamma2, alpha, beta, R_star, T_star,
                                                          T_int, sma, grav, T_int_type)

# get percentiles (for 1,2-sigma boundaries):
low1 = np.percentile(PTprofiles, 16.0, axis=0)
hi1  = np.percentile(PTprofiles, 84.0, axis=0)
median = np.median(PTprofiles, axis=0)

############################### plot figure

plt.fill_betweenx(pressure, low1, hi1, facecolor="orchid",edgecolor="palevioletred")
plt.semilogy(median, pressure, "-", lw=2, label='Median',color="deepskyblue")
plt.legend(loc="best", fontsize=12)
plt.ylim(pressure[0], pressure[-1])

plt.xlim(1250,2000)
plt.xticks([1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000], ["1300", "", "1500", "", "1700", "", "1900", ""])   

plt.subplots_adjust(left=0.20)

plt.savefig('7mol-PT.png')
plt.savefig('7mol-PT.ps')




