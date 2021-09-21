# this code plots the Fp/F* for any output

import numpy as np
import os, sys
import scipy.interpolate as si
from scipy.ndimage.filters import gaussian_filter1d as gaussf
import matplotlib
import matplotlib.pyplot as plt
plt.ion()

sys.path.append("../BART/code")
import makeatm as mat
import wine as w
import readtransit as rt
import bestFit as bf


############################ INPUTS ########################################


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


# 4 mol 4 specs 
direct ='.'
atmfile ='./4species_4opac_uniform/bestFit.atm'


######################### END OF INPUTS ####################################


def plot_bestFit_Spectrum(low, high, xtic, xlab, spec, clr, specs, atmfile, filters, kurucz, tepfile, outflux, data, uncert, direct):
    '''
    Plot Transit spectrum
    '''
    # get star data
    R_star, T_star, sma, gstar = bf.get_starData(tepfile)

    # get surface gravity
    grav, Rp = mat.get_g(tepfile)

    # convert Rp to m
    Rp = Rp * 1000

    # ratio planet to star
    rprs = Rp/R_star
  
    # read kurucz file
    starfl, starwn, tmodel, gmodel = w.readkurucz(kurucz, T_star, gstar)

    # read best-fit spectrum output file, take wn and spectra values
    (head, tail) = os.path.split(outflux)
    specwn, bestspectrum = rt.readspectrum(direct + '/' + tail, wn=True)

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
        bandflux[i] = w.bandintegrate(fluxrat, specwn, nifilter[i],
                                                                 wnindices[i])
        bandmod[i]  = w.bandintegrate(bestspectrum[wnindices[i]],
                                            specwn, nifilter[i], wnindices[i])

    # stellar spectrum on specwn:
    sinterp = si.interp1d(starwn, starfl)
    sflux = sinterp(specwn)
    frat = bestspectrum/sflux * rprs * rprs

    ###################### plot figure #############################

    plt.rcParams["mathtext.default"] = 'rm'
    matplotlib.rcParams.update({'mathtext.default':'rm'})
    matplotlib.rcParams.update({'axes.labelsize': 16,
                                'xtick.labelsize': 20,
                                'ytick.labelsize': 20,})

    plt.figure(2, (8.5, 5))
    plt.clf()

    # smooth the spectrum a fit
    gfrat = gaussf(frat, 2)

    # plot eclipse spectrum
    plt.semilogx(specwl, gfrat*1e3, clr, lw=1.5, label="Spectrum", linewidth=4)
    plt.errorbar(meanwl, data*1e3, uncert*1e3, fmt="ko", label="Data", alpha=0.7)
    plt.ylabel(r"$F_p/F_s$ (10$^{-3}$)", fontsize=26)

    leg = plt.legend(loc="upper left", fontsize=16)
    leg.get_frame().set_alpha(0.5)
    ax = plt.subplot(111)
    ax.set_xscale('log')
    plt.xlabel(r"${\rm Wavelength\ \ (um)}$", fontsize=26)  
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    plt.gca().xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xticks([0.7, 0.8, 0.9, 1.0, 2.0, 3.0, 4.0, 5.0])   
    ax.set_xticklabels(["0.7", "", "", "1.0", "2.0", "3.0", "4.0", "5.0"])
    plt.xlim(min(specwl),max(specwl))

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
            plt.plot(wl, respons, color='grey', linewidth =1, alpha=0.5)
        elif lbl == 'Wang-Hband' or lbl == 'Wang-Kband':
            plt.plot(wl, respons, 'grey', linewidth =1, alpha=0.5)
        elif lbl == 'VLT_1190' or lbl == 'VLT_2090':
            plt.plot(wl, respons, color='grey', linewidth =2, alpha=0.5)
        elif lbl == 'GROND_K_JB' or lbl == 'GROND_i_JB':
            plt.plot(wl, respons, 'grey', linewidth =1, alpha=0.5)
        elif lbl == 'Zhou_Ks':
            plt.plot(wl, respons, 'grey', linewidth =1, alpha=0.5)

    plt.ylim(-0.4, 7)
    plt.text(1.9, 3, specs, color = clr, fontsize = 26)

    plt.subplots_adjust(bottom=0.20)

    plt.savefig(spec +"_BestFit-transSpec.png")
    plt.savefig(spec +"_BestFit-transSpec.ps")


# H2O
low, high = 1e-5, 1e-1
xtic = [1e-4, 1e-3, 1e-2] 
xlab = ["10$^{-4}$", "10$^{-3}$", "10$^{-2}$"]
outflux = 'H2O-ExoMol-1e100-flux.dat'
plot_bestFit_Spectrum(low, high, xtic, xlab, 'H2O', 'b', '$H_2O$', atmfile, filters, kurucz, tep_name, outflux, data, uncert, direct)

# CO2
low, high = 1e-13, 1e-9
xtic = [1e-12, 1e-11, 1e-10] 
xlab = ["10$^{-12}$", "10$^{-11}$", "10$^{-10}$"]
outflux = 'CO2-ExoMol-1e100-flux.dat'
plot_bestFit_Spectrum(low, high, xtic, xlab, 'CO2', 'orangered', '$CO_2$', atmfile, filters, kurucz, tep_name, outflux, data, uncert, direct)

# CO
low, high = 1e-3, 1e1
xtic = [1e-2, 1e-1, 1e0] 
xlab = ["10$^{-2}$", "10$^{-1}$", "10$^{0}$"]
outflux = 'HITRAN2019-CO-1e100-flux.dat'
plot_bestFit_Spectrum(low, high, xtic, xlab, 'CO', 'm', '$CO$', atmfile, filters, kurucz, tep_name, outflux, data, uncert, direct)

# CH4
low, high = 1e-10, 1e-6
xtic = [1e-9, 1e-8, 1e-7] 
xlab = ["10$^{-9}$", "10$^{-8}$", "10$^{-7}$"]
outflux = 'HITEMP2020_CH4-1e100-flux.dat'
plot_bestFit_Spectrum(low, high, xtic, xlab, 'CH4', 'cadetblue', '$CH_4$', atmfile, filters, kurucz, tep_name, outflux, data, uncert, direct)

