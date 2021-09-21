
####### IMPORTANT NOTE: 
# the code works for Python 3.6 


##################################
# Source files used in this paper
##################################

'''
1. BART source files are provided in BART/ folder

2. Pressure file is provided in press_file/ folder

3. Kurucz file is provided in kurucz/ folder

4. CIA files are provided in CIA/ folder

5. Tep file carrying all the system parameters 
   is provided in tepfile/ folder

6. All filter files are provided in the WASP43b_filt/ folder
'''

##################################
# ExoMol data and Partition functions
##################################

'''
1. Download the HITRAN/Exomol/repack data (to be uploaded) in the corresponding molecule directory
2. For Exomol data the partition functions are provided in the same directory
'''

# CO
cd ../EXOMOL/HITRAN-HITEMP/CO
wget https://hitran.org/hitemp/data/bzip2format/05_HITEMP2019.par.bz2

# CH4
cd ../EXOMOL/HITRAN-HITEMP/CH4
wget https://hitran.org/hitemp/data/bzip2format/06_HITEMP2020.par.bz2

# CO2
cd ../EXOMOL/CO2
wget https://zenodo.org/record/3768504/files/pcubillos/CO2_exomol_ucl4000_0.5-500.0um_100-3500K_threshold_0.01_lbl.dat
../EXOMOL/CO2/PF_exomol_CO2.dat

# H2O
cd ../EXOMOL/H2O
wget https://zenodo.org/record/3768504/files/pcubillos/H2O_exomol_pokazatel_0.24-500.0um_100-3500K_threshold_0.01_lbl.dat
../EXOMOL/H2O/PF_exomol_H2O.dat

# C2H2
cd ../EXOMOL/C2H2
wget https://zenodo.org/record/3768504/files/pcubillos/C2H2_exomol_acety_1.0-500.0um_100-3500K_threshold_0.03_lbl.dat
../EXOMOL/C2H2/PF_exomol_C2H2.dat

# C2H4
cd ../EXOMOL/C2H4
wget https://zenodo.org/record/3768504/files/pcubillos/C2H4_exomol_mayty_1.4-500um_100-3000K_threshold_0.03_lbl.dat
../EXOMOL/C2H4/PF_hitran_C2H4.dat

# HCN
cd ../EXOMOL/HCN
wget https://zenodo.org/record/3768504/files/pcubillos/HCN_exomol_harris-larner_0.56-500um_100-3500K_threshold_0.01_lbl.dat
../EXOMOL/HCN/PF_exomol_HCN.dat

# NH3
cd ../EXOMOL/NH3
wget https://zenodo.org/record/3768504/files/pcubillos/NH3_exomol_coyute-byte_0.5-500.0um_100-3000K_threshold_0.03_lbl.dat
../EXOMOL/NH3/PF_hitran_NH3.dat

# TiO
cd ../EXOMOL/TiO
wget https://zenodo.org/record/3768504/files/pcubillos/TiO_exomol_toto_0.33-500um_100-3500K_threshold_0.01_lbl.dat
../EXOMOL/TiO/PF_exomol_TiO.dat

# VO
cd ../EXOMOL/VO
wget https://zenodo.org/record/3768504/files/pcubillos/VO_exomol_vomyt_0.29-500um_100-3500K_threshold_0.01_lbl.dat
../EXOMOL/VO/PF_exomol_VO.dat

# H2S
cd ../EXOMOL/H2S
wget https://zenodo.org/record/3768504/files/pcubillos/H2S_exomol_ayt2_0.28-500.0um_100-3000K_threshold_0.01_lbl.dat
../EXOMOL/H2S/PF_exomol_H2S.dat


#######################
# MAKE TLI FILES
#######################

cd tli_files

# CO
../BART/modules/transit/pylineread/src/pylineread.py -c pyline_HITEMP2019_CO_0.61-5.5um.cfg

# CH4
../BART/modules/transit/pylineread/src/pylineread.py -c pyline_HITEMP2020_CH4_0.61-5.5um.cfg

# CO2
../BART/modules/transit/pylineread/src/pylineread.py -c pyline_CO2-repack_0.61-5.5um.cfg

# H2O
../BART/modules/transit/pylineread/src/pylineread.py -c pyline_H2O-repack_0.61-5.5um.cfg

# C2H2
../BART/modules/transit/pylineread/src/pylineread.py -c pyline_C2H2-repack_0.61-5.5um.cfg

# C2H4
../BART/modules/transit/pylineread/src/pylineread.py -c pyline_C2H4-repack_0.61-5.5um.cfg

# HCN
../BART/modules/transit/pylineread/src/pylineread.py -c pyline_HCN-repack_0.61-5.5um.cfg

# NH3
../BART/modules/transit/pylineread/src/pylineread.py -c pyline_NH3-repack_0.61-5.5um.cfg

# TiO
../BART/modules/transit/pylineread/src/pylineread.py -c pyline_TiO-repack_0.61-5.5um.cfg

# VO
../BART/modules/transit/pylineread/src/pylineread.py -c pyline_VO-repack_0.61-5.5um.cfg

# H2S
../BART/modules/transit/pylineread/src/pylineread.py -c pyline_H2S-repack_0.61-5.5um.cfg

# 4 molecules for Section 3.3
../BART/modules/transit/pylineread/src/pylineread.py -c pyline_4mol_0.61-5.5um.cfg

# 7 molecules for Section 3.4
../BART/modules/transit/pylineread/src/pylineread.py -c pyline_7mol_0.61-5.5um.cfg


#####################
# MAKE OPACITY FILES
#####################

cd opacity_files

# make opacity file using transit_4species.cfg file
../BART/modules/transit/transit/transit -c transit_4species.cfg

# make opacity file using transit_7species.cfg file
../BART/modules/transit/transit/transit -c transit_7species.cfg

####################
# FIGURE 5
####################

cd Figure_5
# use provided transit_[moleculeName].cfg files 
# run them as:

# for all the tli files and all the molecules do

../BART/modules/transit/transit/transit -c transit_[moleculeName].cfg

# run code to make plots of spectra
# plot_Spec_output.py
ipython
copy paste the code

# this run will also make the corresponding opacity files
# and place them in the same directory

# run code to make plots of opacities
# read_opacity_mol.py
ipython
copy paste the code


#######################
# SECTION 3.3 FIGURES
#######################

cd Section_3.3
# execute BART.cfg file as
../BART/BART.py -c BART_4species-4opac-uniform.cfg 

# execute all transit_[moleculeName].cfg
../BART/modules/transit/transit/transit -c transit_[moleculeName].cfg

# then copy and paste the code from each of the following modules
# to produce the plots from Section 3.3
plot_spectrum.py
plot_Hist.py
plot_PTprofiles.py
plot_Spec_output-indivSpecs.py

# Figure 9 is given in /4species_4opac_uniform/ folder
# it is named output_pairwise.png


#######################
# SECTION 3.4 FIGURES
#######################

cd Section_3.4
# execute BART.cfg file as
../BART/BART.py -c BART_7species_7opac_uniform.cfg 

# execute all transit_[moleculeName].cfg
../BART/modules/transit/transit/transit -c transit_[moleculeName].cfg

# then copy and paste the code from each of the following modules
# to produce the plots from Section 3.4
plot_spectrum.py
plot_Hist.py
plot_PTprofiles.py
plot_Spec_output-indivSpecs.py

# Figure 13 is given in /7species_7opac_uniform/ folder
# it is named output_pairwise.png


####################
# FIGURE 14
####################

cd Figure_14
# left panel is placed in Section 3.3/4species_4opac_uniform/
# it is named MCMC_PTprofiles_cf.png

# right panel
plot_toomuch.py





