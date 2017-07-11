import sys
import os
import numpy as np
import argparse
import re
from scipy import interpolate
from scipy.special import erfc
from astropy.io import fits
import bisect
import math

''' Be sure to "fix" the yale files before running this script. 
DO NOT FIX MORE THAN ONCE OR ELSE SCRIPT WILL PRODUCE INCORRECT AGES
'''


inpath = '/home/richmond/V2/interpolated/' #path to isochrone files
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

def find_Metallicity(Yale_file):
    i = 0
    with Yale_file as data:
        for line in data:
            if '[Fe/H]' in line:
                metallicity = float(line[52:61])
    Yale_file.close()
    return metallicity

def getAge(Yale_file):
    Ages = []
    i = 0
    term = "age(Gyr)="
    with Yale_file as infile:
        for line in infile:
            if term in line:
                #print line.index('=')
                m= re.search('\d', line) #find the first integer in the line
                x=  m.start()
                Age = line[x:17] #make that the starting index for age
                try:
                    if float(Age) > 0.000:
                        Ages.append(float(Age))
                    else:
                        Ages.append(float(line[10:16]))
                except (ValueError):
                    pass
    return Ages
    Yale_file.close()

def CompiledAges(listOfAges, initial_mass, columns):
    finalAges = []
    i = 0
    j = 1
    k = 0
    for item in columns:
        np.set_printoptions(precision=14)
        try:
            if columns[j, 0] - columns[i, 0] <-0.2: #sometimes the mass values do not increase within an isochrone as assumed, we are only looking for when the differences are large to know when to jump to the next age
                nextAge = listOfAges[k]
                k += 1
                j += 1
                finalAges.append(float(nextAge))
                i += 1
            elif columns[j, 0] - columns[i, 0] >-0.2:
                nextAge = listOfAges[k]
                j += 1
                finalAges.append(float(nextAge))
                i += 1
            else:
                finalAges.append(float(nextAge))
                break
        except (IndexError):
            finalAges.append(float(nextAge))
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
    j = 0
    #print len(initial_mass), len(allAges)
    i = 0 
    for mass in initial_mass:
        if i < 139: #banking on every isochrone having 140 masses
            single_Isochrone_Mass_List.append(mass)
            i += 1
        else:
            single_Isochrone_Mass_List.append(mass)
            Isochrone_Masses[j] = single_Isochrone_Mass_List
            single_Isochrone_Mass_List = []
            i = 0 
            j += 1
    for isochrone in sorted(Isochrone_Masses):
        masses = Isochrone_Masses[isochrone]
        binborders = []
        for i, item in enumerate(sorted(masses)):
            try:
                if masses[i] != masses[i+1]: #some masses repeat in Yale... super annoying
                    binborders.append(np.average([masses[i], masses[i+1]]))
            except (IndexError):
                pass
        binborders.insert(0, masses[0] - (np.average([masses[0], masses[1]])-masses[0])) #sets the lower bin border
        binborders.append(sorted(masses)[len(masses)-1]-(np.average([sorted(masses)[len(masses)-1], sorted(masses)[len(masses)-2]])-sorted(masses)[len(masses)-1]))#sets the upper bin border.needs to be sorted!
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
                #print len(binborders)
                print left_index, right_index
                print binborders[left_index], binborders[right_index]
                break
    #print len(weights)
    return weights

def fitsTable(allAges, metallicitylist, initial_mass, Teff, log_g, _2Mass_J, _2Mass_H, _2Mass_Ks, outpath, listOfAges, age_weights, mass_weights):
    metallicitylist = np.array(metallicitylist)
    x = np.array(allAges) 
    initial_mass = np.array(initial_mass)
    Teff = np.array(Teff)
    log_g = np.array(log_g)
    _2Mass_J = np.array(_2Mass_J)
    _2Mass_H = np.array(_2Mass_H)
    _2Mass_Ks = np.array(_2Mass_Ks)
    evo_stage = np.array([])
    final_weights = np.array(age_weights) * np.array(mass_weights)
    col1 = fits.Column(name='age', format='E', array=x)
    col2 = fits.Column(name='feh', format='E', array=metallicitylist)
    col3 = fits.Column(name='T', format='E', array=Teff)
    col4 = fits.Column(name='logg', format='E', array=log_g)
    col5 =  fits.Column(name='mass', format='E', array=initial_mass)
    col6 =  fits.Column(name='J', format='E', array=_2Mass_J)
    col7 =  fits.Column(name='H', format='E', array=_2Mass_H)
    col8 =  fits.Column(name='K', format='E', array=_2Mass_Ks)
    col9 =  fits.Column(name='stage', format='K', array=evo_stage) #'K' --> 64-bit integer
    col10 =  fits.Column(name='age weights', format='E', array=age_weights)
    #col11 =  fits.Column(name='feh weights', format='E', array=metallicityWeights) #wasn't working for some reason but parameter is not used 
    col12 =  fits.Column(name='mass weights', format='E', array=mass_weights)
    col13 =  fits.Column(name='weight', format='E', array=final_weights)
    cols = fits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col12, col13])
    tbhdu = fits.BinTableHDU.from_columns(cols)
    ageGridhdu = fits.BinTableHDU.from_columns([fits.Column(name='Age Grid', format='E', array=listOfAges)])
    agelisthdr = fits.Header()
    agelisthdr['Ages'] = 'Age grid'
    agelisthdu = fits.PrimaryHDU(header=agelisthdr)
    tbhdulist = fits.HDUList([agelisthdu, tbhdu, ageGridhdu])
    tbhdulist.writeto(outpath + 'Yale.fits', clobber=True) #fits outfile
    
 
#inputAge = str(0.150)# raw_input('')
initial_mass = [] 
log_Teff = []
log_g = []
_2Mass_J = []
_2Mass_H = []
_2Mass_Ks = []
Teff = []
metallicitylist = []
allAges = []

for file in sorted(os.listdir(inpath)):
    file1 = file
    columns = np.loadtxt(inpath + file1)
    initial_mass.extend(columns[:, 0])
    log_Teff.extend(columns[:, 1])
    log_g.extend(columns[:, 3])
    _2Mass_J.extend(columns[:, 4] - columns[:, 9])
    _2Mass_H.extend(columns[:, 4] - columns[:, 10])
    _2Mass_Ks.extend(columns[:, 4] - columns[:, 11])
    Yale_file = open(inpath +file1, 'r')
    metallicity = find_Metallicity(Yale_file)
    for item in columns:
        metallicitylist.append(metallicity)
    Yale_file = open(inpath +file1, 'r')
    listOfAges = getAge(Yale_file) #this only really needs to be performed once
    Yale_file = open(inpath +file1, 'r')
    Ages = CompiledAges(listOfAges, initial_mass, columns)
    allAges.extend(Ages)
    print file1 + '\t' + str(metallicity)

allAges = np.array(allAges).dot(1000000000)
allAges = np.log10(allAges)
listOfAges = np.unique(allAges).tolist()
#print listOfAges
mass_weights = massWeights(initial_mass, allAges)
age_weights = Weights(allAges, listOfAges)
#print len(age_weights), len(mass_weights)
metallicities = np.unique(metallicitylist).tolist()
#metallicityWeights = Weights(metallicitylist, metallicities) 
Teff.extend(np.power(10, log_Teff))
fitsTable(allAges, metallicitylist, initial_mass, Teff, log_g, _2Mass_J, _2Mass_H, _2Mass_Ks, outpath, listOfAges, age_weights, mass_weights)
print 'FITS file created for Yale'

'''
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
'''
