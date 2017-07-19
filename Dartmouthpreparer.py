import sys
import os
import numpy as np
import argparse
from astropy.io import fits
import scipy.signal
from scipy import interpolate
from scipy.special import erfc
import bisect
import preparer

inpath = '/home/richmond/Isochrones/Dartmouth/interpolated/test/' #path to isochrone files
outpath = '/home/richmond/Isochrones/FITS/' #desired output for FITS files

def find_Metallicity(Dartmouth_file):
    with Dartmouth_file as data:
        for line in data:
            if '[Fe/H]' in line:
                x = data.next()
                metallicity = x[40:46]
    Dartmouth_file.close()
    return metallicity.strip()

def getAge(Dartmouth_file, x):
    Ages = []
    i = 0
    y = x - 5 
    with Dartmouth_file as infile:
        for line in infile:
            if i < y:
                i += 1
            elif i >= y:
                if "AGE=" in line:
                    Ages.append(float(line[5:12]))
    return Ages
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
    binborders.append(listOfAges[-1]-((listOfAges[-1] + listOfAges[-2])/2.-listOfAges[-1]))
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
        binborders.append(sorted(masses)[-1] - ((sorted(masses)[-1] + sorted(masses)[-2])/2. -sorted(masses)[-1]))#sets the upper bin border 
        for single_mass in masses:
            try:
                right_index = bisect.bisect_left(binborders, single_mass)  #indicies of bin borders of which mass falls between 
                left_index = right_index -1 
                left_mass = binborders[left_index]
                right_mass = binborders[right_index]
                left_CDF = preparer.cdf_mass(left_mass)
                right_CDF = preparer.cdf_mass(right_mass)
                weights.append(right_CDF - left_CDF)
            except (IndexError): #if something goes wrong it will report the troublesome value
                pass
    return weights

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
        Dartmouth_file = open(inpath +file1, 'r')
        listOfAges = getAge(Dartmouth_file, x)
        Ages = CompiledAges(listOfAges, initial_mass, columns)
        allAges.extend(Ages)
        print file1 + '\t' + metallicity

allAges = np.array(allAges).dot(1000000000)
allAges = np.log10(allAges)
mass_weights = massWeights(initial_mass, allAges)
listOfAges = np.unique(allAges).tolist()
age_weights = Weights(allAges, listOfAges)
metallicities = np.unique(metallicitylist).tolist()
evo_stage = np.array([])
Teff.extend(np.power(10, log_Teff))
preparer.fitsTable(age_weights, allAges, metallicitylist, Teff, log_g, initial_mass, _2Mass_J, _2Mass_H, _2Mass_Ks, evo_stage, mass_weights, listOfAges, outpath, 'Dartmouth')

