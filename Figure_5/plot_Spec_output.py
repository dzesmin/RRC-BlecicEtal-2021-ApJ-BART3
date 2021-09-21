# this code plots the Fp/F* for any output

import numpy as np
import os
import scipy.interpolate as si
from scipy.ndimage.filters import gaussian_filter1d as gaussf
import matplotlib
import matplotlib.pyplot as plt
plt.ion()

import sys
sys.path.append("../BART/code")
import makeatm as mat
import wine as w
import readtransit as rt
import bestFit as bf


############################ INPUTS ########################################

direct = '../Figure_5'

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

kurucz = '../kurucz/WASP43b-fp00ak2odfnew.pck'

tepfile = '../tepfile/WASP-43b.tep'


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

spec = np.array(['H2O',\
                 'CO2',\
                 'CO',\
                 'CH4',\
                 'NH3',\
                 'HCN',\
                 'C2H2',\
                 'C2H4',\
                 'H2S',\
                 'TiO',\
                 'VO'])



outflux = np.array(['H2O-ExoMol-1e100-flux.dat',\
                    'CO2-ExoMol-1e100-flux.dat',\
                    'HITRAN2019-CO-1e100-flux.dat',\
                    'HITEMP2020_CH4-1e100-flux.dat',\
                    'NH3-ExoMol-1e100-flux.dat',\
                    'HCN-ExoMol-1e100-flux.dat',\
                    'C2H2-ExoMol-1e100-flux.dat',\
                    'C2H4-ExoMol-1e100-flux.dat',\
                    'H2S-ExoMol-1e100-flux.dat',\
                    'TiO-ExoMol-1e100-flux.dat',\
                    'VO-ExoMol-1e100-flux.dat'])


clr = np.array(['b',\
                'orangered',\
                'm',\
                'cadetblue',\
                'firebrick',\
                'darkgreen',\
                'palevioletred',\
                'mediumpurple',\
                'mediumturquoise',\
                'grey',\
                'darkorange'])


atmfile = '../atmfiles/WASP43b_11mol_bestPT.tea'


specs = np.array(['$H_2O$',\
                  '$CO_2$',\
                  '$CO$',\
                  '$CH_4$',\
                  '$NH_3$',\
                  '$HCN$',\
                  '$C_2H_2$',\
                  '$C_2H_4$',\
                  '$H_2S$',\
                  '$TiO$',\
                  '$VO$'])





######################### END OF INPUTS ####################################


def plot_bestFit_Spectrum(spec, clr, specs, atmfile, filters, kurucz, tepfile, outflux, data, uncert, direct):
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
    #matplotlib.rcParams.update({'fontsize': 10,})
    matplotlib.rcParams.update({'axes.labelsize': 16,
                                #'text.fontsize':   10,
                                'legend.fontsize': 14,
                                'xtick.labelsize': 20,
                                'ytick.labelsize': 20,})

    plt.figure(2, (8.5, 5))
    plt.clf()
    #plt.xlim(0.60, 5.5)

    #plt.xlim(min(specwl),max(specwl))

    # plot eclipse spectrum
    #gfrat = gaussf(frat, 0)
    plt.semilogx(specwl, frat*1e3, clr, lw=1.5, label="Spectrum", linewidth=4)
#cornflowerblue, lightskyblue
    #plt.errorbar(meanwl, data*1e3, uncert*1e3, fmt="ko", label="Data", alpha=0.7)
    plt.errorbar(meanwl, data*1e3, uncert*1e3, fmt=".", color = 'k', zorder=100, capsize=2, capthick=1, label="Data", alpha=0.7)
    #plt.plot(meanwl, bandflux*1e3, "ko", label="model")
    plt.ylabel(r"$F_p/F_s$ (10$^{-3}$)", fontsize=24)

    leg = plt.legend(loc="upper left")
    #leg.draw_frame(False)
    leg.get_frame().set_alpha(0.5)

    ax = plt.subplot(111)
    ax.set_xscale('log')

    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_xticks([0.7, 0.8, 0.9, 1.0, 2.0, 3.0, 4.0, 5.0])   
    ax.set_xticklabels(["0.7", "", "", "1.0", "2.0", "3.0", "4.0", "5.0"])

    plt.xlabel(r"${\rm Wavelength\ \ (um)}$", fontsize=24)  

    nfilters = len(filters)
    # plot filter bandpasses
    for i in np.arange(nfilters-15):
        (head, tail) = os.path.split(filters[i])
        lbl = tail[:-4]
        # read filter:
        wn, respons = w.readfilter(filters[i])
        respons = respons/3 -0.4
        wl = 10000.0/wn
        #plt.plot(wl, respons, color='crimson', linewidth =1)
        if lbl == 'spitzer_irac1_sa' or lbl== 'spitzer_irac2_sa':
            respons = respons*2 +0.4
            plt.plot(wl, respons, color='grey', linewidth =1, alpha=0.5)
            #plt.plot(wl, respons*2, color='orangered', linewidth =1)
        elif lbl == 'Wang-Hband' or lbl == 'Wang-Kband':
            plt.plot(wl, respons, 'grey', linewidth =1, alpha=0.5)
        elif lbl == 'VLT_1190' or lbl == 'VLT_2090':
            plt.plot(wl, respons, color='grey', linewidth =2, alpha=0.5)
            #plt.plot(wl, respons, color='firebrick', linewidth =2)
        elif lbl == 'GROND_K_JB' or lbl == 'GROND_i_JB':
            plt.plot(wl, respons, 'grey', linewidth =1, alpha=0.5)
        elif lbl == 'Zhou_Ks':
            plt.plot(wl, respons, 'grey', linewidth =1, alpha=0.5)

    plt.ylim(-0.4, 7)

    plt.text(1.9, 3, specs, color = clr, fontsize = 26)



    ###################### INSET PT and ABUN FIGURE ####################
    b = plt.axes([.21, .45, .14, .24])

    # read atmfile
    molecules, pres, temp, abundances = mat.readatm(atmfile)

    plt.semilogy(temp, pres, color='r', linewidth=3)
    plt.xlim(1000, 2200)
    plt.ylim(max(pres), min(pres))
    b.minorticks_off()
    yticks = [1e2, 1e1, 1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5]
    ylabels=["10$^{2}$", "", "10$^{0}$", "", "10$^{-2}$", "", "10$^{-4}$", ""]
    plt.yticks(yticks, ylabels, fontsize=8)
    xticks = [1000, 1200, 1400, 1600, 1800, 2000, 2200]
    xlabels=["", "1200", "", "", "1800", ""]
    plt.xticks(xticks, xlabels, fontsize=12)
    plt.xlabel('T (K)', fontsize =12)
    plt.ylabel('P (bar)', fontsize=12)


    # ############################## SECOND INSET ABUN
    c = plt.axes([.35, .45, .14, .24])

    # Sets the second argument given as the species names
    species = spec

    # Open the atmospheric file and read
    f = open(atmfile, 'r')
    lines = np.asarray(f.readlines())
    f.close()

    # Get molecules names
    imol = np.where(lines == "#SPECIES\n")[0][0] + 1
    molecules = lines[imol].split()
    nmol = len(molecules)
    for m in np.arange(nmol):
        molecules[m] = molecules[m].partition('_')[0]

    nspec = 1

    # Populate column numbers for requested species and 
    #          update list of species if order is not appropriate
    columns = []
    spec    = []
    for i in np.arange(nmol):
        if molecules[i] == species:
            columns.append(i+3)   # defines p, T +2 or rad, p, T +3
            spec.append(species)

    # Convert spec to tuple
    spec = tuple(spec)

    # Concatenate spec with pressure for data and columns
    data    = tuple(np.concatenate((['p'], spec)))
    usecols = tuple(np.concatenate(([1], columns))) # defines p as 0 columns, or p as 1 columns

    # Load all data for all interested species
    data = np.loadtxt(atmfile, dtype=float, comments='#', delimiter=None,    \
                converters=None, skiprows=13, usecols=usecols, unpack=True)

    plt.loglog(data[1], data[0], '-', color=clr, \
                                    linewidth=3)

    plt.ylim(max(pres), min(pres))
    c.minorticks_off()
    yticks = [1e2, 1e1, 1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5]
    ylabels=[]
    plt.yticks(yticks, ylabels)
    plt.xlim(1e-12, 1e-2) 
    xticks = [1e-11, 1e-9, 1e-7,1e-5, 1e-3] 
    xlabels = ["10$^{-11}$", "", "10$^{-5}$", "", "10$^{-3}$"]
    plt.xticks(xticks, xlabels, fontsize=12)
    plt.xlabel('Mix. fraction', fontsize=12)

    plt.subplots_adjust(bottom=0.16)

    spec = spec[0]
    print(spec)    

    plt.savefig(spec +"_transSpec_new.png")
    plt.savefig(spec +"_transSpec_new.ps")


for i in np.arange(len(spec)):
    plot_bestFit_Spectrum(spec[i], clr[i], specs[i], atmfile, filters, kurucz, tepfile, outflux[i], data, uncert, direct)
