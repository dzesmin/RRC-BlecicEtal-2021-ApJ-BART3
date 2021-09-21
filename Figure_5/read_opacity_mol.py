import numpy as np
import string
import matplotlib.pylab as plt
import matplotlib
import os
import sys

sys.path.append("../BART/code")
import wine as w
plt.ion()

opacity_file = np.array(['opacity_file_H2O-ExoMol_0.61-5.5um_300_3000K_1e100.dat',\
                         'opacity_file_CO2-ExoMol_0.61-5.5um_300_3000K_1e100.dat',\
                         'opacity_file_HITRAN2019_CO_061-55_300_3000K_1e100.dat',\
                         'opacity_file_HITEMP2020_CH4_0.61-5.5um_300_3000K_1e100.dat',\
                         'opacity_file_NH3-ExoMol_0.61-5.5um_300_3000K_1e100.dat',\
                         'opacity_file_HCN-ExoMol_0.61-5.5um_300_3000K_1e100.dat',\
                         'opacity_file_C2H2-ExoMol_0.61-5.5um_300_3000K_1e100.dat',\
                         'opacity_file_C2H4-ExoMol_0.61-5.5um_300_3000K_1e100.dat',\
                         'opacity_file_H2S-ExoMol_0.61-5.5um_300_3000K_1e100.dat',\
                         'opacity_file_TiO-ExoMol_0.61-5.5um_300_3000K_1e100.dat',\
                         'opacity_file_VO-ExoMol_0.61-5.5um_300_3000K_1e100.dat'])


clr = np.array(['b', 'orangered', 'm', 'cadetblue', 'firebrick', 'darkgreen', 'palevioletred', 'mediumpurple', 'mediumturquoise', 'grey', 'darkorange'])

loc = np.array([[0.8, 1e3], [1.45, 1e3], [1.45, 1e3], [1.45, 1e3], [1.3, 1e3], [1.5, 1e4], [1.67, 1e3], [2.2, 1e3], [0.9, 5e-1], [0.9, 1e-2], [1.8, 3e1]]) 

specs = np.array(['$H_2O$', '$CO_2$', '$CO$', '$CH_4$', '$NH_3$', '$HCN$', '$C_2H_2$', '$C_2H_4$', '$H_2S$', '$TiO$', '$VO$'])

spec_names = np.array(['H2O', 'CO2', 'CO', 'CH4', 'NH3', 'HCN', 'C2H2', 'C2H4', 'H2S', 'TiO', 'VO'])


direct = '../Figure_4'

def plot_opacity(direct, file, clr, spec, spec_names, loc):
   

    # Read the binary opacity file
    f = open(direct + '/' + file, "r")

    # read file dimension sizes, Nmol, Ntemp, Nlayer, Nwave
    a = np.fromfile(f, dtype=np.uint64, count=4)
    Nmol, Ntemp, Nlayer, Nwave  = a[0], a[1], a[2], a[3]

    # read the molecular IDs, molID
    molID = np.fromfile(f, dtype=np.uint32, count=a[0])

    # read temperatures, temp
    temp = np.fromfile(f, dtype=np.float64, count=a[1])

    # read pressures in each layer, pres
    pres = np.fromfile(f, dtype=np.float64, count=a[2])
    # convert barie to bars
    pres = pres/1e6
    # revert to the decreasing order (from the top to the bottom of the atmosphere)
    pres = pres[::-1]

    # read wavenumber sample, wns
    wns = np.fromfile(f, dtype=np.float64, count=a[3])

    # read opacities
    opacities = np.fromfile(f, dtype=np.float64, count=a[0]*a[1]*a[2]*a[3])

    # correctly reshape opacities Nlayer, Ntemp, Nmol, Nwave, 
    # see transit/src/opacity.c
    opac_matrix = opacities.reshape(Nlayer, Ntemp, Nmol, Nwave)

    # opacities per molecule, H2O, CH4, CO, CO2
    mol1_opac = opac_matrix[:,:,0,:]

    # find index of 1bar and 1500K pre and temp
    inx_1bar = np.where(pres>1.0)[0][0]
    inx_1500K = np.where(temp==1500.0)[0][0]

    ##################### PLOTTING

    # plot opacities devided by abundances
    plt.figure(2, (8.5, 5))

    plt.clf()

    plt.loglog(10000/wns, mol1_opac[inx_1bar, inx_1500K, :], color= clr, linewidth=3)
    plt.text(loc[0], loc[1], spec, color= clr, fontsize = 26)


    plt.rcParams["mathtext.default"] = 'rm'
    matplotlib.rcParams.update({'mathtext.default':'rm'})
    #matplotlib.rcParams.update({'font size': 10,})
    matplotlib.rcParams.update({'axes.labelsize': 16,
                                #'text.fontsize':   10,
                                'legend.fontsize': 14,
                                'xtick.labelsize': 20,
                                'ytick.labelsize': 20,})


    #plt.legend(loc='best', prop={'size':10})
    #plt.xlim(0.61, 5.5)
    labels =["0.7", "", "", "1.0", "2.0", "3.0", "4.0", "5.0"]
    plt.xticks([0.7, 0.8, 0.9, 1.0, 2.0, 3.0, 4.0, 5.0], labels, fontsize=20)
    plt.ylim(1e-12, 1e7)
    labels =["10$^{-11}$", "10$^{-7}$", "10$^{-3}$", "10$^{1}$", "10$^{5}$"]
    plt.yticks([1e-11, 1e-7, 1e-3, 1e1, 1e5], labels, fontsize=20)
  

    plt.xlabel(r"${\rm Wavelength\ \ (um)}$", fontsize=24)
    plt.ylabel(r"${\rm Opacity/density\ \ (cm^{2}/g)}$", fontsize=24)



    ########################### FILTERS

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

    plt.subplots_adjust(left=0.16, bottom=0.16)

    # change the name of every plot to save
 
    plt.savefig(spec_names + "-opacity.png")
    plt.savefig(spec_names + "-opacity.ps")


for i in np.arange(len(specs)):
    plot_opacity(direct, opacity_file[i], clr[i], specs[i], spec_names[i], loc[i])






