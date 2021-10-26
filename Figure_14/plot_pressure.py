# reads tau_toomuch

import numpy as np
import string
import matplotlib.pyplot as plt
plt.ion()


def read_file(fileName):
     # opens the file to read
     f = open(fileName, 'r')

     # allocates array of pressures
     wl = []
     tau_max = []
     rad = []
     rad_ind = []


     # finds the line of interest in the header
     for line in f:
          if string.find(line, '   (microns) ') > -1:

               # reads the file from the row below the header
               for line in f:
                    data = line.split()

                    # appends the data in every row from character 11th to 20th
                    wl      = np.append(wl     , data[0])
                    tau_max = np.append(tau_max, data[1])
                    rad     = np.append(rad    , data[2])
                    rad_ind = np.append(rad_ind, data[3])

               # converts strings to floats
               wl      = wl.astype(float)
               tau_max = tau_max.astype(float)
               rad     = rad.astype(float)
               rad_ind = rad_ind.astype(float)

     return wl, tau_max, rad, rad_ind 


def plot_pressure(fileName, plotName):
    wl, tau_max, rad, rad_ind = read_file(fileName)
    pressure = np.array([  1.00000000e+02,   8.49750000e+01,   7.22080000e+01,
         6.13590000e+01,   5.21400000e+01,   4.43060000e+01,
         3.76490000e+01,   3.19930000e+01,   2.71860000e+01,
         2.31010000e+01,   1.96300000e+01,   1.66810000e+01,
         1.41750000e+01,   1.20450000e+01,   1.02350000e+01,
         8.69750000e+00,   7.39070000e+00,   6.28030000e+00,
         5.33670000e+00,   4.53490000e+00,   3.85350000e+00,
         3.27450000e+00,   2.78260000e+00,   2.36450000e+00,
         2.00920000e+00,   1.70740000e+00,   1.45080000e+00,
         1.23280000e+00,   1.04760000e+00,   8.90220000e-01,
         7.56460000e-01,   6.42810000e-01,   5.46230000e-01,
         4.64160000e-01,   3.94420000e-01,   3.35160000e-01,
         2.84800000e-01,   2.42010000e-01,   2.05650000e-01,
         1.74750000e-01,   1.48500000e-01,   1.26190000e-01,
         1.07230000e-01,   9.11160000e-02,   7.74260000e-02,
         6.57930000e-02,   5.59080000e-02,   4.75080000e-02,
         4.03700000e-02,   3.43050000e-02,   2.91510000e-02,
         2.47710000e-02,   2.10490000e-02,   1.78860000e-02,
         1.51990000e-02,   1.29150000e-02,   1.09750000e-02,
         9.32600000e-03,   7.92480000e-03,   6.73420000e-03,
         5.72240000e-03,   4.86260000e-03,   4.13200000e-03,
         3.51120000e-03,   2.98360000e-03,   2.53540000e-03,
         2.15440000e-03,   1.83070000e-03,   1.55570000e-03,
         1.32190000e-03,   1.12330000e-03,   9.54550000e-04,
         8.11130000e-04,   6.89260000e-04,   5.85700000e-04,
         4.97700000e-04,   4.22920000e-04,   3.59380000e-04,
         3.05390000e-04,   2.59500000e-04,   2.20510000e-04,
         1.87380000e-04,   1.59230000e-04,   1.35300000e-04,
         1.14980000e-04,   9.77010000e-05,   8.30220000e-05,
         7.05480000e-05,   5.99480000e-05,   5.09410000e-05,
         4.32880000e-05,   3.67840000e-05,   3.12570000e-05,
         2.65610000e-05,   2.25700000e-05,   1.91790000e-05,
         1.62980000e-05,   1.38490000e-05,   1.17680000e-05,
         1.00000000e-05])

    pres = pressure[::-1]
    press = np.zeros(len(wl), dtype=float)
    for i in np.arange(len(wl)):
        press[i] = pres[int(rad_ind[i])]

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 5))

    ax.set_xlabel('Wavelength (um)', fontsize=22)
    ax.set_ylabel('Pressure (bar)', fontsize=22) 
    ax.loglog(wl, press, color='darkseagreen')
    ax.set_xlim(0.57, 6) 
    ax.set_ylim(3e2, 1e-5)
    plt.yticks(fontsize=14)  
    ax.set_xticks([0.7, 0.8, 0.9, 1.0, 2.0, 3.0, 4.0, 5.0])  
    ax.set_xticklabels(["0.7", "", "", "1.0", "2.0", "3.0", "4.0", "5.0"], fontsize=14)  
    plt.subplots_adjust(left=None, bottom=0.2, right=None, top=None)

    plt.savefig(plotName + '-pressure.png')
    plt.savefig(plotName + '-pressure.ps')


# 4mol 4 opac
fileName = '../Section_3.3/4species_4opac_uniform/4species_4opac_uniform-toom.dat'
plotName = '4mol_4opac'
plot_pressure(fileName, plotName)


