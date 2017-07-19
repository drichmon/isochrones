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
import preparer

''' Be sure to "fix" the yale files before running this script. 
DO NOT FIX MORE THAN ONCE OR ELSE SCRIPT WILL PRODUCE INCORRECT AGES
'''


inpath = '/home/richmond/V2/test/' #path to isochrone files
outpath = '/home/richmond/Isochrones/FITS/' #desired output for FITS files

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
    binborders.append(listOfAges[-1]-(np.average([listOfAges[-1], listOfAges[-2]])-listOfAges[-1]))
    binborders = np.power(10, binborders)
    weightslist = binborders[1:] - binborders[:-1]
    for item in allAges:
        weights.append(weightslist[listOfAges.index(item)])
    return weights

def massWeights(initial_mass, allAges): #mass grids in PARSEC isochrones differ between individual isochrones. There is not a standard grid, which is slightly annoying.
    single_Isochrone_Mass_List = []
    Isochrone_Masses = {}
    weights = []
    j = 0
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
        binborders.append(sorted(masses)[-1]-(np.average([sorted(masses)[-1], sorted(masses)[-2]])-sorted(masses)[-1]))#sets the upper bin border.needs to be sorted!
        for single_mass in masses:
            try:
                right_index = bisect.bisect_left(binborders, single_mass)  #indicies of bin borders of which mass falls between 
                left_index = right_index -1
                left_mass = binborders[left_index]
                right_mass = binborders[right_index]
                left_CDF = preparer.cdf_mass(left_mass)
                right_CDF = preparer.cdf_mass(right_mass)
                weights.append(right_CDF - left_CDF)
            except (IndexError): 
                print single_mass
                print left_index, right_index
                print binborders[left_index], binborders[right_index]
                break
    return weights

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
    listOfAges = getAge(Yale_file) 
    Yale_file = open(inpath +file1, 'r')
    Ages = CompiledAges(listOfAges, initial_mass, columns)
    allAges.extend(Ages)
    print file1 + '\t' + str(metallicity)

allAges = np.array(allAges).dot(1000000000)
allAges = np.log10(allAges)
listOfAges = np.unique(allAges).tolist()
mass_weights = massWeights(initial_mass, allAges)
age_weights = Weights(allAges, listOfAges)
metallicities = np.unique(metallicitylist).tolist()
Teff.extend(np.power(10, log_Teff))
evo_stage = []
preparer.fitsTable(age_weights, allAges, metallicitylist, Teff, log_g, initial_mass, _2Mass_J, _2Mass_H, _2Mass_Ks, evo_stage, mass_weights, listOfAges, outpath, 'Yale')
