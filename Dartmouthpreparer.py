import sys
import os
import numpy as np
import argparse
from astropy.io import fits
import scipy.signal
from scipy import interpolate
from scipy.special import erfc
import bisect


inpath = '/home/richmond/Isochrones/Dartmouth/interpolated/' #path to isochrone files
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

def find_Metallicity(Dartmouth_file):
    i = 0
    with Dartmouth_file as data:
        for line in data:
            if '[Fe/H]' in line:
                x = data.next()
                metallicity = x[40:46]
    Dartmouth_file.close()
    return metallicity.strip()
            
def log_teff_to_teff(log_Teff):
    i = 0
    for item in log_Teff:
        item = 10 ** item
        np.put(log_Teff, i, item)
        i +=1
    return log_Teff


def getAge(Dartmouth_file, x):
    Ages = []
    i = 0
    y = x - 5 # didn't really feel like figuring out the exact difference, 5 is good enough
    with Dartmouth_file as infile:
        for line in infile:
            if i < y:
                i += 1
            elif i >= y:
                if "AGE=" in line:
                    Ages.append(float(line[5:12]))
    return Ages
    Dartmouth_file.close()

def get_skiprows(inputAge, Dartmouth_file):
    i = 1
    term = "#AGE= " + inputAge
    for line in Dartmouth_file:
        if term not in line:
            i += 1
        elif term in line:
            i += 1
            return i
    Dartmouth_file.close()

def CompiledAges(listOfAges, initial_mass, columns):
    finalAges = []
    i = 0
    j = 1
    k = 0
    for item in columns:
        np.set_printoptions(precision=14)
        try:
            if columns[j, 0] - columns[i, 0] <0:
                nextAge = listOfAges[k]
                k += 1
                j += 1
                finalAges.append(nextAge)
                i += 1
            elif columns[j, 0] - columns[i, 0] >0 :
                nextAge = listOfAges[k]
                j += 1
                finalAges.append(nextAge)
                i += 1
            else:
                finalAges.append(nextAge)
                break
        except (IndexError):
            finalAges.append(nextAge)
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
        for i, item in enumerate(sorted(masses)):
            try:
                binborders.append(np.average([masses[i], masses[i+1]]))
            except (IndexError):
                pass
                
        binborders.insert(0, masses[0] - (np.average([masses[0], masses[1]])-masses[0])) #sets the lower bin border
        binborders.append(sorted(masses)[len(masses)-1]-(np.average([sorted(masses)[len(masses)-1], sorted(masses)[len(masses)-2]])-sorted(masses)[len(masses)-1]))#sets the upper bin border 
        for single_mass in masses:
            try:
                right_index = bisect.bisect_left(binborders, single_mass)  #indicies of bin borders of which mass falls between 
                left_index = right_index -1 
                left_mass = binborders[left_index]
                right_mass = binborders[right_index]
                left_CDF = cdf_mass(left_mass)
                right_CDF = cdf_mass(right_mass)
                weights.append(right_CDF - left_CDF)
            except (IndexError): #if something goes wrong it will report the troublesome value
                print single_mass
                print len(binborders)
                print left_index, right_index
                print binborders[left_index], binborders[right_index]
                break
    print len(weights)
    return weights

def fitsTable(allAges, metallicitylist, initial_mass, Teff, log_g, _2Mass_J, _2Mass_H, _2Mass_Ks, listOfAges, outpath, age_weights, mass_weights):
    metallicityList = np.array(metallicitylist)   
    x = np.array(allAges) 
    initial_mass = np.array(initial_mass)
    Teff = np.array(Teff)
    log_g = np.array(log_g)
    _2Mass_J = np.array(_2Mass_J)
    _2Mass_H = np.array(_2Mass_H)
    _2Mass_Ks = np.array(_2Mass_Ks)
    evo_stage = np.array([])
    ones = np.full(len(allAges), 1)
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
    #col11 =  fits.Column(name='feh weights', format='E', array=metallicityWeights)#giving troubles and not really needed 
    col12 =  fits.Column(name='mass weights', format='E', array=mass_weights)
    col13 =  fits.Column(name='weight', format='E', array=final_weights)
    cols = fits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col12, col13])
    tbhdu = fits.BinTableHDU.from_columns(cols)
    ageGridhdu = fits.BinTableHDU.from_columns([fits.Column(name='Age Grid', format='E', array=listOfAges)])
    agelisthdr = fits.Header()
    agelisthdr['Ages'] = 'Age grid'
    agelisthdu = fits.PrimaryHDU(header=agelisthdr)
    tbhdulist = fits.HDUList([agelisthdu, tbhdu, ageGridhdu])
    tbhdulist.writeto(outpath + 'Dartmouth.fits', clobber=True) #fits outfile

EEP = []
initial_mass = [] 
log_Teff = []
log_g = []
_2Mass_J = []
_2Mass_H = []
_2Mass_Ks = []
Teff = []
metallicitylist = []
allAges = []

#inputAge = str(1.000)
for file in sorted(os.listdir(inpath)):
    file1 = file
    x = 5
    if 'dat' in file1:
        columns = np.loadtxt(inpath + file1)
        EEP.extend(columns[:, 0])
        initial_mass.extend(columns[:, 1])
        log_Teff.extend(columns[:, 2])
        log_g.extend(columns[:, 3])
        _2Mass_J.extend(columns[:, 10])
        _2Mass_H.extend(columns[:, 11])
        _2Mass_Ks.extend(columns[:, 12])
        Dartmouth_file = open(inpath +file1, 'r')
        metallicity = find_Metallicity(Dartmouth_file)
        Dartmouth_file = open(inpath +file1, 'r')
        for item in columns:
            metallicitylist.append(metallicity)
        #print len(metallicitylist)
        Dartmouth_file = open(inpath +file1, 'r')
        listOfAges = getAge(Dartmouth_file, x)
        print file1 + '\t' + metallicity
        Ages = CompiledAges(listOfAges, initial_mass, columns)
        allAges.extend(Ages)


allAges = np.array(allAges).dot(1000000000)
allAges = np.log10(allAges)
mass_weights = massWeights(initial_mass, allAges)
listOfAges = np.unique(allAges).tolist()
age_weights = Weights(allAges, listOfAges)
metallicities = np.unique(metallicitylist).tolist()
#metallicityWeights = Weights(metallicitylist, metallicities)
Teff.extend(np.power(10, log_Teff))
#print len(Teff)
fitsTable(allAges, metallicitylist, initial_mass, Teff, log_g, _2Mass_J, _2Mass_H, _2Mass_Ks, listOfAges, outpath, age_weights, mass_weights)

print 'FITS file created for Dartmouth'
