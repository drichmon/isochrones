import numpy as np
import scipy.signal
import re
import os
from collections import OrderedDict
from astropy.io import fits
import scipy.signal
from scipy import interpolate
from scipy.special import erfc
import bisect


inpath = '/home/richmond/Isochrones/MIST/MIST_v1.0_vvcrit0.4_UBVRIplus/' #path to isochrone files
outpath = '/home/richmond/Isochrones/FITS/' #desired output for FITS files


def _cdf_lowmass_Kroupa(mass):
    a_1 = 0.3; A = 1.785;
    return A * np.power(mass, 1-a_1) / (1 - a_1)

def _cdf_midmass_Kroupa(mass):
    m_0 = 0.08; a_2 = 1.3; A = 1.785; B = 0.4352; k_1 = 0.08;
    return B + A * k_1 / (1-a_2) * (np.power(mass, 1-a_2) - np.power(m_0, 1-a_2))
    
def _cdf_highmass_Kroupa(mass):
    m_1 = 0.5; a_3 = 2.3; A = 1.785; C = 0.86469; k_2 = 0.04;
    return C + A * k_2 / (1-a_3) * (np.power(mass, 1-a_3) - np.power(m_1, 1-a_3))

def cdf_mass(mass):
    if mass < 0.08:
        return _cdf_lowmass_Kroupa(mass)
    elif mass < 0.5:
        return _cdf_midmass_Kroupa(mass)
    else:
        return _cdf_highmass_Kroupa(mass)

def find_Metallicity(MIST_file):
    i = 0
    with MIST_file as data:
        for line in data:
            if '[Fe/H]' in line:
                x = data.next()
                metallicity = x[25:31]
                break
    MIST_file.close()
    return float(metallicity.strip())

def Weights(allAges, listOfAges):
    weightslist = []
    binborders = []
    weights = []
    for i, item in enumerate(listOfAges):
        try:
            binborders.append(np.average([listOfAges[i], listOfAges[i+1]]))
        except (IndexError):
            pass
    binborders.insert(0, listOfAges[0]-(np.average([listOfAges[0], listOfAges[1]])-listOfAges[0]))
    binborders.append(listOfAges[len(listOfAges)-1]-(np.average([listOfAges[len(listOfAges)-1], listOfAges[len(listOfAges)-2]])-listOfAges[len(listOfAges)-1]))
    binborders = np.power(10, binborders)
    for i, item in enumerate(binborders):
        try:
            weightslist.append(binborders[i+1]-binborders[i])
        except (IndexError):
            pass
    for item in allAges:
        weights.append(weightslist[listOfAges.index(item)])
    return weights
        

def massWeights(initial_mass, allAges): #mass grids in PARSEC isochrones differ between individual isochrones. There is not a standard grid, which is slightly annoying.
    single_Isochrone_Mass_List = []
    Isochrone_Masses = {}
    weights = []
    print len(initial_mass), len(allAges)
    for i, mass in enumerate(initial_mass):
        try:
            if allAges[i] == allAges[i+1]:
                single_Isochrone_Mass_List.append(mass)
            else:
                single_Isochrone_Mass_List.append(mass)
                Isochrone_Masses[i] = single_Isochrone_Mass_List
                single_Isochrone_Mass_List = []
        except (IndexError):
            single_Isochrone_Mass_List.append(mass)
            Isochrone_Masses[i] = single_Isochrone_Mass_List
            pass
    for isochrone in sorted(Isochrone_Masses):
        masses = Isochrone_Masses[isochrone]
        binborders = []
        for i, item in enumerate(masses):
            try:
                binborders.append(np.average([masses[i], masses[i+1]]))
            except (IndexError):
                pass
        binborders.insert(0, masses[0] - (np.average([masses[0], masses[1]])-masses[0])) #sets the lower bin border
        binborders.append(masses[len(masses)-1]-(np.average([masses[len(masses)-1], masses[len(masses)-2]])-masses[len(masses)-1]))#sets the upper bin border 
        for single_mass in masses:
            right_index = bisect.bisect_left(binborders, single_mass)  #indicies of bin borders of which mass falls between 
            left_index = right_index -1 
            left_mass = binborders[left_index]
            right_mass = binborders[right_index]
            left_CDF = cdf_mass(left_mass)
            right_CDF = cdf_mass(right_mass)
            weights.append(right_CDF - left_CDF)
    print len(weights)
    return weights

def fitsTable(logAge, metallicitylist, initial_mass, Teff, log_g, _2Mass_J, _2Mass_H, _2Mass_Ks, evo_stage, outpath, listOfAges, age_weights, metallicityWeights, mass_weights):
    logAge = np.array(logAge)
    initial_mass = np.array(initial_mass)
    Teff = np.array(Teff)
    log_g = np.array(log_g)
    _2Mass_J = np.array(_2Mass_J)
    _2Mass_H = np.array(_2Mass_H)
    _2Mass_Ks = np.array(_2Mass_Ks)
    evo_stage = np.array(evo_stage)
    metallicitylist = np.array(metallicitylist)
    final_weights = np.array(age_weights) * np.array(mass_weights)
    col1 = fits.Column(name='age', format='E', array=logAge)
    col2 = fits.Column(name='feh', format='E', array=metallicitylist)
    col3 = fits.Column(name='T', format='E', array=Teff)
    col4 = fits.Column(name='logg', format='E', array=log_g)
    col5 =  fits.Column(name='mass', format='E', array=initial_mass)
    col6 =  fits.Column(name='J', format='E', array=_2Mass_J)
    col7 =  fits.Column(name='H', format='E', array=_2Mass_H)
    col8 =  fits.Column(name='K', format='E', array=_2Mass_Ks)
    col9 =  fits.Column(name='stage', format='K', array=evo_stage)
    col10 =  fits.Column(name='age weights', format='E', array=age_weights)
    col11 =  fits.Column(name='feh weights', format='E', array=metallicityWeights)
    col12 =  fits.Column(name='mass weights', format='E', array=mass_weights)
    col13 =  fits.Column(name='weight', format='E', array=final_weights)
    cols = fits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13])
    tbhdu = fits.BinTableHDU.from_columns(cols)
    ageGridhdu = fits.BinTableHDU.from_columns([fits.Column(name='Age Grid', format='E', array=listOfAges)])
    agelisthdr = fits.Header()
    agelisthdr['Ages'] = 'Age grid'
    agelisthdu = fits.PrimaryHDU(header=agelisthdr)
    tbhdulist = fits.HDUList([agelisthdu, tbhdu, ageGridhdu])
    tbhdulist.writeto(outpath + 'MIST.fits', clobber=True) #fits outfile


logAge = []
initial_mass = []
log_Teff = []
log_g = []
_2Mass_J = []
_2Mass_H = []
_2Mass_Ks = []
evo_stage = []
metallicitylist = []
Teff = []
listOfAges = []


for file in sorted(os.listdir(inpath)):
    if 'cmd' in file: # makes sure it only iterates over isochrone files
        columns = np.loadtxt(inpath + file)
        #print columns.shape[1]
        if columns.shape[1] == 24:
            logAge.extend(columns[:, 1])
            initial_mass.extend(columns[:,2])
            log_Teff.extend(columns[:,3])
            log_g.extend(columns[:,4])
            _2Mass_J.extend(columns[:,12])
            _2Mass_H.extend(columns[:,13])
            _2Mass_Ks.extend(columns[:,14])
            evo_stage.extend(columns[:, 23])
            MIST_file = open(inpath +file, 'r')
            metallicity = find_Metallicity(MIST_file)
            for item in columns:
                metallicitylist.append(metallicity)
            #metallicityarray.extend(metallicitylist[file])
            print file + '\t' + str(metallicity)
        else:
            logAge.extend(columns[:, 1])
            initial_mass.extend(columns[:,2])
            log_Teff.extend(columns[:,4])
            log_g.extend(columns[:,5])
            _2Mass_J.extend(columns[:,14])
            _2Mass_H.extend(columns[:,15])
            _2Mass_Ks.extend(columns[:,16])
            evo_stage.extend(columns[:, 25])
            MIST_file = open(inpath +file, 'r')
            metallicity = find_Metallicity(MIST_file)
            for item in columns:
                metallicitylist.append(metallicity)
            print file + '\t' + str(metallicity)

mass_weights = massWeights(initial_mass, logAge)
listOfAges = np.unique(logAge).tolist()
age_weights = Weights(logAge, listOfAges)
metallicities = np.unique(metallicitylist).tolist()
metallicityWeights = Weights(metallicitylist, metallicities)

for n, i in enumerate(evo_stage): #converts MIST evolutionary stages into stages UniDAM can comprehend
    if i < 3:
        evo_stage[n] = 1 #pre CHeB
    elif i == 3:
        evo_stage[n] = 2 #CHeB
    else:
        evo_stage[n] = 3 #post CHeB

Teff.extend(np.power(10, log_Teff)) # log_teff_to_teff(log_Teff)
fitsTable(logAge, metallicitylist, initial_mass, Teff, log_g, _2Mass_J, _2Mass_H, _2Mass_Ks, evo_stage, outpath, listOfAges, age_weights, metallicityWeights, mass_weights)
print 'FITS file created for MIST'

