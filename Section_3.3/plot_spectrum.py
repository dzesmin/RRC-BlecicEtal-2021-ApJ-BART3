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

# get star data
R_star, T_star, sma, gstar = bf.get_starData(tep_name)
# get surface gravity
grav, Rp = mat.get_g(tep_name)
# convert Rp to m
Rp = Rp * 1000
# ratio planet to star
rprs = Rp/R_star
# read kurucz file
starfl, starwn, tmodel, gmodel = w.readkurucz(kurucz, T_star, gstar)

###################### plot figure #############################

plt.rcParams["mathtext.default"] = 'rm'
matplotlib.rcParams.update({'mathtext.default':'rm'})
matplotlib.rcParams.update({'axes.labelsize': 16,
                                'xtick.labelsize': 14,
                                'ytick.labelsize': 14,})

fig= plt.figure(figsize=(8.5, 5))
plt.ylabel(r"$F_p/F_s$ (10$^{-3}$)", fontsize=14)
plt.xlabel(r"${\rm Wavelength\ \ (um)}$", fontsize=14)

######################################### for each spectrum separately

outflux = './4species_4opac_uniform/4species_4opac_uniform-flux.dat'

# read best-fit spectrum output file, take wn and spectra values
specwn, bestspectrum = rt.readspectrum(outflux, wn=True)
# convert wn to wl
specwl = 1e4/specwn
# number of filters
nfilters = len(filters)
# read and resample the filters:
nifilter  = [] # Normalized interpolated filter
istarfl   = [] # interpolated stellar flux
wnindices = [] # wavenumber indices used in interpolation
meanwn    = [] # Filter mean wavenumber
for i in np.arange(nfilters):
    # read filter:
    filtwaven, filttransm = w.readfilter(filters[i])
    meanwn.append(np.sum(filtwaven*filttransm)/sum(filttransm))
    # resample filter and stellar spectrum:
    nifilt, strfl, wnind = w.resample(specwn, filtwaven, filttransm,
                                            starwn,    starfl)
    nifilter.append(nifilt)
    istarfl.append(strfl)
    wnindices.append(wnind)

# convert mean wn to mean wl
meanwl = 1e4/np.asarray(meanwn)

# band-integrate the flux-ratio or modulation:
bandflux = np.zeros(nfilters, dtype='d')
bandmod  = np.zeros(nfilters, dtype='d')
for i in np.arange(nfilters):
    fluxrat = (bestspectrum[wnindices[i]]/istarfl[i]) * rprs*rprs
    bandflux[i] = w.bandintegrate(fluxrat, specwn, nifilter[i], wnindices[i])
    bandmod[i]  = w.bandintegrate(bestspectrum[wnindices[i]],
                                            specwn, nifilter[i], wnindices[i])

# stellar spectrum on specwn:
sinterp = si.interp1d(starwn, starfl)
sflux = sinterp(specwn)
frat = bestspectrum/sflux * rprs * rprs

# plot eclipse spectrum
gfrat = gaussf(frat, 1)
plt.semilogx(specwl, gfrat*1e3, "slateblue", lw=1.5, label="Best-fit", alpha=0.6)
plt.errorbar(meanwl, data*1e3, uncert*1e3, fmt="ro", label="Data", alpha=0.8)
plt.plot(meanwl, bandflux*1e3, "ko", label="Model", alpha=1.0)
leg = plt.legend(loc="upper left")
leg.get_frame().set_alpha(0.5)

nfilters = len(filters)
# plot filter bandpasses
for i in np.arange(nfilters-15):
    (head, tail) = os.path.split(filters[i])
    lbl = tail[:-4]
    # read filter:
    wn, respons = w.readfilter(filters[i])
    respons = respons/3 -0.4
    wl = 10000.0/wn
    if lbl == 'spitzer_irac1_sa' or lbl== 'spitzer_irac2_sa':
        respons = respons*2 +0.4
        plt.plot(wl, respons, color='grey', linewidth =1, alpha=0.3)
    elif lbl == 'Wang-Hband' or lbl == 'Wang-Kband':
        plt.plot(wl, respons, 'grey', linewidth =1, alpha=0.3)
    elif lbl == 'VLT_1190' or lbl == 'VLT_2090':
        plt.plot(wl, respons, color='grey', linewidth =2, alpha=0.3)
    elif lbl == 'GROND_K_JB' or lbl == 'GROND_i_JB':
        plt.plot(wl, respons, 'grey', linewidth =1, alpha=0.3)
    elif lbl == 'Zhou_Ks':
        plt.plot(wl, respons, 'grey', linewidth =1, alpha=0.3)

plt.gca().xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
plt.xticks([0.6, 0.8, 1.0, 2.0, 3.0, 4.0, 5.0],["0.6", "0.8", "1.0", "2.0", "3.0", "4.0", "5.0"], fontsize=12) 
  
plt.ylim(-0.4, 7)
plt.xlim(0.61, 5.5)

###################### INSET FIGURE ####################

b = plt.axes([.36, .52, .33, .23])
# plot eclipse spectrum
gfrat = gaussf(frat, 1)
plt.semilogx(specwl, gfrat*1e3, "slateblue", lw=1.5, label="Best-fit", alpha=0.6)
plt.errorbar(meanwl, data*1e3, uncert*1e3, fmt="ro", label="Data", alpha=0.8)
plt.plot(meanwl, bandflux*1e3, "ko", label="Model", alpha=1.0)
xticks = [1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7]
xlabels=["1.1", "1.2", "1.3", "1.4", "1.5", "1.6", 1.7]
plt.xticks(xticks, xlabels, fontsize=8 )  
yticks = [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
ylabels=["0.0", "", "0.4", "", "0.8", "", 1.2]
plt.yticks(yticks, ylabels, fontsize=8) 
plt.xlim(1.1, 1.78)
plt.ylim(-0.1, 1.3)

for i in np.arange(nfilters):
    (head, tail) = os.path.split(filters[i])
    lbl = tail[:-4]
    # read filter:
    wn, respons = w.readfilter(filters[i])
    respons = respons/3.6 -0.3
    wl = 10000.0/wn
    if lbl.startswith("grism"):
        plt.plot(wl, respons, color='grey', linewidth =1, alpha=0.3)
    elif lbl == 'Wang-Hband' or lbl == 'Wang-Kband':
        respons = respons*2+0.1
        plt.plot(wl, respons, 'grey', linewidth =1, alpha=0.3)
    elif lbl == 'VLT_1190' or lbl == 'VLT_2090':
        respons = respons *2 +0.12
        plt.plot(wl, respons*2+0.2, color='grey', linewidth =1, alpha=0.3)       

    ###################### INSET FIGURE ####################

plt.subplots_adjust(bottom=0.14)

######################################### for each speactrum separately

plt.savefig('4mol-Spectrum.png')
plt.savefig('4mol-Spectrum.ps')





