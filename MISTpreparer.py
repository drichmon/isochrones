import numpy as np
import scipy.signal
import re
import os
from collections import OrderedDict
import scipy.signal
from scipy import interpolate
from scipy.special import erfc
import bisect
import preparer
import time
start_time = time.time()


inpath = '/home/richmond/Isochrones/MIST/MIST_v1.0_vvcrit0.4_UBVRIplus/interpolated/' #path to isochrone files
outpath = '/home/richmond/Isochrones/FITS/' #desired output for FITS files

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
    binborders.append(listOfAges[-1] - ((listOfAges[-1] + listOfAges[-2])/2.-listOfAges[-1]))
    binborders = np.power(10, binborders)
    for i, item in enumerate(binborders):
        try:
            weightslist.append(binborders[i+1]-binborders[i])
        except (IndexError):
            pass
    for item in allAges:
        weights.append(weightslist[listOfAges.index(item)])
    return weights
        

def massWeights(initial_mass, allAges): 
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
        binborders.append(masses[-1]-((masses[-1]+masses[-2])/2-masses[-1]))#sets the upper bin border 
        for single_mass in masses:
            right_index = bisect.bisect_left(binborders, single_mass)  #indicies of bin borders of which mass falls between 
            left_index = right_index -1 
            left_mass = binborders[left_index]
            right_mass = binborders[right_index]
            left_CDF = preparer.cdf_mass(left_mass)
            right_CDF = preparer.cdf_mass(right_mass)
            weights.append(right_CDF - left_CDF)
    return weights

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
age_weights = Weights(logAge, listOfAges)
listOfAges = np.unique(logAge).tolist()

for n, i in enumerate(evo_stage): #converts MIST evolutionary stages into stages UniDAM can comprehend
    if i < 3:
        evo_stage[n] = 1 #pre CHeB
    elif i == 3:
        evo_stage[n] = 2 #CHeB
    else:
        evo_stage[n] = 3 #post CHeB

Teff.extend(np.power(10, log_Teff)) # log_teff_to_teff(log_Teff)
preparer.fitsTable(age_weights, logAge, metallicitylist, Teff, log_g, initial_mass, _2Mass_J, _2Mass_H, _2Mass_Ks, evo_stage, mass_weights, listOfAges, outpath, 'MIST')
print 'run time = %s' % (time.time() - start_time)
