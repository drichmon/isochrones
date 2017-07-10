import sys
import os
import numpy as np
from scipy import interpolate
from scipy.special import erfc
from astropy.io import fits
import bisect
import math

inpath = '/home/richmond/Isochrones/PARSEC/' #path to isochrone files
outpath = '/home/richmond/Isochrones/FITS/' #desired output for FITS files
'''#if you ever decide to use Chabrier IMF
# Chabrier CDF taken from Jumper & Fisher (2013) w/ params from Chabrier (2005)
def _cdf_lowmass_Chabrier(mass):
    A_s = 0.724;  sig = 0.55; m_c = 0.2;
    return A_s * sig * np.sqrt(np.pi/2.) * erfc((np.log10(m_c) - np.log10(mass)) / (sig * np.sqrt(2.0)))
def _cdf_highmass_Chabrier(mass):
    B_s = 0.323; x = 1.35; m_o = 1.; C = 0.896;
    return C + (B_s / np.log(10.)) * (np.power(m_o,-x) / x) * (1. - np.power((mass / m_o),-x))
'''
# Using Kroupa (2001) w/ CDF calculated by hand
# k values used to connect each part of the Kroupa IMF
def _cdf_lowmass_Kroupa(mass):
    a_1 = 0.3; A = 1.785;
    return A * np.power(mass, 1-a_1) / (1 - a_1)

def _cdf_midmass_Kroupa(mass):
    m_0 = 0.08; a_2 = 1.3; A = 1.785; B = 0.4352; k_1 = 0.08;
    return B + A * k_1 / (1-a_2) * (np.power(mass, 1-a_2) - np.power(m_0, 1-a_2))
    
def _cdf_highmass_Kroupa(mass):
    m_1 = 0.5; a_3 = 2.3; A = 1.785; C = 0.86469; k_2 = 0.04;
    return C + A * k_2 / (1-a_3) * (np.power(mass, 1-a_3) - np.power(m_1, 1-a_3))
'''
# Takes a 2 element array representing the mass range for the provided IMF 
# distribution and produces n_samples number of stars in the given mass range
def generate_masses(mass_range,distrub='Chabrier',n_samples=10000):
    rang = np.arange(mass_range[0],mass_range[1],0.001)
    if distrub == 'Chabrier':
        m_0 = 1.;
        cdf = _cdf_lowmass_Chabrier(rang[rang <= m_0])
        cdf = np.append(cdf, _cdf_highmass_Chabrier(rang[rang > m_0]))
        inv_cdf = interpolate.interp1d(cdf, rang)
    elif distrub == 'Kroupa':
        m_0 = 0.08; m_1 = 0.5;
        cdf = _cdf_lowmass_Kroupa(rang[rang <= m_0])
        cdf = np.append(cdf, _cdf_midmass_Kroupa(rang[(rang > m_0) & (rang <= m_1)]))
        cdf = np.append(cdf, _cdf_highmass_Kroupa(rang[rang > m_1]))
        inv_cdf = interpolate.interp1d(cdf, rang)
    else:
        raise NameError("The " + distrub + " IMF is not provided in this method")
        return None 
    r = np.random.uniform(np.min(cdf),np.max(cdf),n_samples) 
    return inv_cdf(r)
'''
def cdf_mass(mass):
    if mass < 0.08:
        return _cdf_lowmass_Kroupa(mass)
    elif mass < 0.5:
        return _cdf_midmass_Kroupa(mass)
    else:
        return _cdf_highmass_Kroupa(mass)


def find_Metallicity(PARSEC_file):
    i = 0
    with PARSEC_file as data:
        for line in data:
            if '[M/H]' in line:
                metallicity = float(line[45:51])
    PARSEC_file.close()
    return metallicity
            
def log_teff_to_teff(log_Teff):
    i = 0
    for item in log_Teff:
        item = 10 ** item
        np.put(log_Teff, i, item)
        i +=1
    return log_Teff
        
def getAge(PARSEC_file, x):
    Ages = []
    i = 0
    term = "Age ="
    with PARSEC_file as infile:
        for line in infile:
            if term in line:
                #print line[72:78] + '\t' + line[80:83]
                Age = line[72:83]#float(line[72:78]) * 10 ** int(line[80:83]) #PREVIOUSLY +1 TO ALL INDICIES IN THIS LINE
                Ages.append(float(Age))
    return Ages
    PARSEC_file.close()

def get_skiprows(inputAge, PARSEC_file):
    i = 1
    term = "Age = \t" + str(inputAge)
    for line in PARSEC_file:
        if term not in line:
            i += 1
        elif term in line:
            i += 1
            return i
    PARSEC_file.close()

def CompiledAges(listOfAges, initial_mass, columns):
    finalAges = []
    i = 0
    j = 1
    k = 0
    for item in columns:
        np.set_printoptions(precision=14)
        try:
            if columns[j, 1] != columns[i, 1]:
                k += 1
                nextAge = listOfAges[k]
                j += 1
                #print nextAge
                finalAges.append(nextAge)
                #metallicityList.append(metallicity)
                i += 1
            elif columns[j, 1] == columns[i, 1]:
                nextAge = listOfAges[k]
                j += 1
                #print nextAge
                finalAges.append(nextAge)
                #metallicityList.append(metallicity)
                i += 1
            else:
                finalAges.append(nextAge)
                break
        except (IndexError):
            pass    
    return finalAges

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

def fitsTable(allAges, allmetallicities, initial_mass, Teff, log_g, _2Mass_J, _2Mass_H, _2Mass_Ks, evo_stage, outpath, listOfAges, age_weights, metallicityWeights, mass_weights):
    metallicityList = np.array(allmetallicities)  
    x = np.array(allAges) 
    initial_mass = np.array(initial_mass)
    Teff = np.array(Teff)
    log_g = np.array(log_g)
    _2Mass_J = np.array(_2Mass_J)
    _2Mass_H = np.array(_2Mass_H)
    _2Mass_Ks = np.array(_2Mass_Ks)
    metallicityArray = np.array(metallicityList)
    evo_stage = np.array(evo_stage)
    final_weights = np.array(age_weights) * np.array(mass_weights)
    col1 = fits.Column(name='age', format='E', array=x)
    col2 = fits.Column(name='feh', format='E', array=metallicityList)
    col3 = fits.Column(name='T', format='E', array=Teff)
    col4 = fits.Column(name='logg', format='E', array=log_g)
    col5 =  fits.Column(name='mass', format='E', array=initial_mass)
    col6 =  fits.Column(name='J', format='E', array=_2Mass_J)
    col7 =  fits.Column(name='H', format='E', array=_2Mass_H)
    col8 =  fits.Column(name='K', format='E', array=_2Mass_Ks)
    col9 = fits.Column(name='stage', format='K', array=evo_stage)
    col10 =  fits.Column(name='age weights', format='E', array=age_weights)
    col11 =  fits.Column(name='feh weights', format='E', array=metallicityWeights)
    col12 =  fits.Column(name='mass weights', format='E', array=mass_weights)
    col13 =  fits.Column(name='weight', format='E', array=final_weights)
    cols = fits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13])
    tbhdu = fits.BinTableHDU.from_columns(cols)
    ageGridhdu = fits.BinTableHDU.from_columns([fits.Column(name='Age Grid', format='E', array=listOfAges)])
    agelisthdr = fits.Header()
    agelisthdr['Ages'] = 'Age grid' #don't really know what I'm doing here/if this is necessary
    agelisthdu = fits.PrimaryHDU(header=agelisthdr)
    tbhdulist = fits.HDUList([agelisthdu, tbhdu, ageGridhdu])
    tbhdulist.writeto(outpath + 'PARSEC.fits', clobber=True) #fits outfile

initial_mass = [] 
log_Teff = []
log_g = []
_2Mass_J = []
_2Mass_H = []
_2Mass_Ks = []
Teff = []
allmetallicities = []
allAges = []
evo_stage = []

for file in sorted (os.listdir(inpath)): 
    if '.dat' in file: #makes sure we only iterate over the isochrone files
        x = 5
        columns = np.loadtxt(inpath + file)
        initial_mass.extend(columns[:, 2])
        log_Teff.extend(columns[:, 5])
        log_g.extend(columns[:, 6])
        _2Mass_J.extend(columns[:, 13])
        _2Mass_H.extend(columns[:, 14])
        _2Mass_Ks.extend(columns[:, 15])
        allAges.extend(columns[:, 1])
        evo_stage.extend(columns[:, 17])
        PARSEC_file = open(inpath +file, 'r')
        metallicity = find_Metallicity(PARSEC_file)
        for item in columns:
            allmetallicities.append(metallicity)
        print file + '\t' + str(metallicity)

mass_weights = massWeights(initial_mass, allAges)
listOfAges = np.unique(allAges).tolist()
age_weights = Weights(allAges, listOfAges)
metallicities = np.unique(allmetallicities).tolist()
metallicityWeights = Weights(allmetallicities, metallicities) #What is the distribution of metallicities? May have to change this if most stars do not have metallicities of 0.5...


for n, i in enumerate(evo_stage): #converts PARSEC evolutionary stages into stages UniDAM can comprehend
    if i <= 3:
        evo_stage[n] = 1 #pre CHeB
    elif 4 <= i <= 6:
        evo_stage[n] = 2 #CHeB
    else:
        evo_stage[n] = 3 #post CHeB

Teff.extend(np.power(10, log_Teff))

fitsTable(allAges, allmetallicities, initial_mass, Teff, log_g, _2Mass_J, _2Mass_H, _2Mass_Ks, evo_stage, outpath, listOfAges, age_weights, metallicityWeights, mass_weights)

print 'FITS file created for PARSEC'
