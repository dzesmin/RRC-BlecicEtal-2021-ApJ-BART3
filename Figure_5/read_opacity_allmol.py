import numpy as np
import string
import matplotlib.pylab as plt
import matplotlib
import os
import sys

sys.path.append("../BART/code")
plt.ion()

opacity_file = np.array(['../opacity_files/output/opacity_file_ExoMol_H2O-CO2_NH3_HCN_C2H2_C2H4_H2S_TiO_VO_HITEMP_CO-CH4_11mol_061-55_300_3000K_1e100.dat'])


def plot_opacity(opacity_file):
   

    # Read the binary opacity file
    f = open(opacity_file, "r")

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
    #pres = pres/1e6
    # revert to the decreasing order (from the top to the bottom of the atmosphere)
    pres = pres[::-1]

    # read wavenumber sample, wns
    wns = np.fromfile(f, dtype=np.float64, count=a[3])

    # read opacities
    opacities = np.fromfile(f, dtype=np.float64, count=a[0]*a[1]*a[2]*a[3])

    # correctly reshape opacities Nlayer, Ntemp, Nmol, Nwave, 
    # see transit/src/opacity.c
    opac_matrix = opacities.reshape(Nlayer, Ntemp, Nmol, Nwave)

    # molID
    # array([101, 103, 104, 102, 106, 113, 111, 112, 115, 107, 108], dtype=uint32)
    # opacities per molecule
    # 101    H2O       18.01528     3.2       01         Water           
    # 103    CO        28.0101      2.8       01         Carbon Monoxide
    # 104    CO2       44.0095      2.8       01         Carbon Dioxide
    # 102    CH4       16.0425      4.0       01         Methane  
    # 106    NH3       17.03052     3.6       01         Ammonia
    # 113    HCN       27.02534     5.0       05         Hydrogen Cyanide
    # 111    C2H2      26.0373      5.26      05         Acetylene
    # 112    C2H4      28.0532      5.69      05         Ethylene
    # 115    H2S       34.0809      3.6       02         Hydrogen Sulfide
    # 107    TiO       63.8664      3.45      03         Titanium Monoxide
    # 108    VO        66.94090     3.32      03         Vanadium Monoxide

    molH2O_opac   = opac_matrix[:,:,0,:]
    molCO_opac    = opac_matrix[:,:,1,:]
    molCO2_opac   = opac_matrix[:,:,2,:]
    molCH4_opac   = opac_matrix[:,:,3,:]
    molNH3_opac   = opac_matrix[:,:,4,:]
    molHCN_opac   = opac_matrix[:,:,5,:]
    molC2H2_opac  = opac_matrix[:,:,6,:]
    molC2H4_opac  = opac_matrix[:,:,7,:]
    molH2S_opac   = opac_matrix[:,:,8,:]
    molTiO_opac   = opac_matrix[:,:,9,:]
    molVO_opac    = opac_matrix[:,:,10,:]

    mol_opac = np.array([molH2O_opac, molCO2_opac, molCO_opac, molCH4_opac, molNH3_opac, molHCN_opac, molC2H2_opac, molC2H4_opac, molH2S_opac, molTiO_opac, molVO_opac])

    clr = np.array(['b', 'orangered', 'm', 'cadetblue', 'firebrick', 'darkgreen', 'palevioletred', 'mediumpurple', 'mediumturquoise', 'grey', 'darkorange'])

    loc = np.array([[0.8, 1e3], [1.45, 1e3], [1.45, 1e3], [1.45, 1e3], [1.3, 1e3], [1.5, 1e4], [1.67, 1e3], [2.2, 1e3], [0.9, 5e-1], [0.9, 1e-2], [1.8, 3e1]]) 

    specs = np.array(['$H_2O$', '$CO_2$', '$CO$', '$CH_4$', '$NH_3$', '$HCN$', '$C_2H_2$', '$C_2H_4$', '$H_2S$', '$TiO$', '$VO$'])

    spec_names = np.array(['H2O', 'CO2', 'CO', 'CH4', 'NH3', 'HCN', 'C2H2', 'C2H4', 'H2S', 'TiO', 'VO'])

    # find index of 1bar and 1500K pre and temp
    inx_1bar = np.where(pres>1.0)[0][0]
    inx_1500K = np.where(temp==1500.0)[0][0]

    ##################### PLOTTING

    # plot opacities devided by abundances

    for i in np.arange(len(clr)):
        plt.figure(i, (8.5, 5))

        plt.clf()

        plt.loglog(10000/wns, mol_opac[i][inx_1bar, inx_1500K, :], color= clr[i], linewidth=3)
        plt.text(loc[i,0], loc[i,1], specs[i], color= clr[i], fontsize = 26)


        plt.rcParams["mathtext.default"] = 'rm'
        matplotlib.rcParams.update({'mathtext.default':'rm'})
        matplotlib.rcParams.update({'axes.labelsize': 16,
                                #'text.fontsize':   10,
                                'legend.fontsize': 14,
                                'xtick.labelsize': 20,
                                'ytick.labelsize': 20,})


        #plt.legend(loc='best', prop={'size':10})
        plt.xlim(0.61, 5.5)
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

        # removes double tick labels from the plot
        plt.gca().xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

        plt.subplots_adjust(left=0.16, bottom=0.16)

        # change the name of every plot to save
 
        plt.savefig(spec_names[i] + "-opacity_new.png")
        plt.savefig(spec_names[i] + "-opacity_new.ps")


plot_opacity(opacity_file)






