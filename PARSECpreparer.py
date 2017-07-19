import sys
import os
import numpy as np
from scipy import interpolate
from scipy.special import erfc
from astropy.io import fits
import bisect
import math
import preparer


inpath = '/home/richmond/Isochrones/PARSEC/test/' #path to isochrone files
outpath = '/home/richmond/Isochrones/FITS/' #desired output for FITS files

def find_Metallicity(PARSEC_file):
    i = 0
    with PARSEC_file as data:
        for line in data:
            if '[M/H]' in line:
                metallicity = float(line[45:51])
    PARSEC_file.close()
    return metallicity

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
                finalAges.append(nextAge)
                i += 1
            elif columns[j, 1] == columns[i, 1]:
                nextAge = listOfAges[k]
                j += 1
                finalAges.append(nextAge)
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
    binborders.append(listOfAges[-1]-((listOfAges[-1] + listOfAges[-2])/2.-listOfAges[-1]))
    binborders = np.power(10, binborders)
    weightslist = binborders[1:] - binborders[:-1] # put this line in the other codes
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
        for i, item in enumerate(masses):
            try:
                binborders.append(np.average([masses[i], masses[i+1]]))
            except (IndexError):
                pass
        binborders.insert(0, masses[0] - (np.average([masses[0], masses[1]])-masses[0])) #sets the lower bin border
        binborders.append(masses[-1]-((masses[-1] + masses[-2])/2.-masses[-1]))#sets the upper bin border 
        for single_mass in masses:
            right_index = bisect.bisect_left(binborders, single_mass)  #indicies of bin borders of which mass falls between 
            left_index = right_index -1 
            left_mass = binborders[left_index]
            right_mass = binborders[right_index]
            left_CDF = preparer.cdf_mass(left_mass)
            right_CDF = preparer.cdf_mass(right_mass)
            weights.append(right_CDF - left_CDF)
    return weights

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

for n, i in enumerate(evo_stage): #converts PARSEC evolutionary stages into stages UniDAM can comprehend
    if i <= 3:
        evo_stage[n] = 1 #pre CHeB
    elif 4 <= i <= 6:
        evo_stage[n] = 2 #CHeB
    else:
        evo_stage[n] = 3 #post CHeB

Teff.extend(np.power(10, log_Teff))
preparer.fitsTable(age_weights, allAges, allmetallicities, Teff, log_g, initial_mass, _2Mass_J, _2Mass_H, _2Mass_Ks, evo_stage, mass_weights, listOfAges, outpath, 'PARSEC')

