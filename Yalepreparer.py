import sys
import os
import numpy as np
import argparse
import re
from astropy.io import fits

''' Be sure to "fix" the yale files before running this script. 
DO NOT FIX MORE THAN ONCE OR ELSE SCRIPT WILL PRODUCE INCORRECT AGES
'''


inpath = '/home/richmond/V2/interpolated/' #path to isochrone files
outpath = '/home/richmond/Isochrones/FITS/' #desired output for FITS files


def find_Metallicity(Yale_file):
    i = 0
    with Yale_file as data:
        for line in data:
            if '[Fe/H]' in line:
                metallicity = float(line[52:61])
    Yale_file.close()
    return metallicity

def getAge(Yale_file, x):
    Ages = []
    i = 0
    y = x - 2 # didn't really feel like figuring out the exact difference
    term = "age(Gyr)="
    with Yale_file as infile:
        for line in infile:
            if i < y:
                i += 1
            elif i >= y:
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
'''
def get_skiprows(inputAge, Yale_file): #currently only works for age values < 10
    i = 1
    term = "age(Gyr)= " + str(inputAge)
    #term2 = "age(Gyr)=" + str(inputAge)
    for line in Yale_file:
        if term not in line:
            i += 1
        elif term in line:
            i += 1
            return i
    Yale_file.close()
'''
def CompiledAges(listOfAges, initial_mass, columns):
    finalAges = []
    i = 0
    j = 1
    k = 0
    for item in columns:
        np.set_printoptions(precision=14)
        try:
            if columns[j, 0] - columns[i, 0] <-0.5: #sometimes the mass values do not increase within an isochrone as assumed, we are only looking for when the differences are large to know when to jump to the next age
                nextAge = listOfAges[k]
                k += 1
                j += 1
                finalAges.append(float(nextAge))
                i += 1
            elif columns[j, 0] - columns[i, 0] >-0.5:
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

def fitsTable(allAges, metallicitylist, initial_mass, Teff, log_g, _2Mass_J, _2Mass_H, _2Mass_Ks, outpath, listOfAges):
    metallicitylist = np.array(metallicitylist)
    x = np.array(allAges) 
    initial_mass = np.array(initial_mass)
    Teff = np.array(Teff)
    log_g = np.array(log_g)
    _2Mass_J = np.array(_2Mass_J)
    _2Mass_H = np.array(_2Mass_H)
    _2Mass_Ks = np.array(_2Mass_Ks)
    evo_stage = np.array([])
    ones = np.full(len(allAges), 1)
    col1 = fits.Column(name='age', format='E', array=x)
    col2 = fits.Column(name='feh', format='E', array=metallicitylist)
    col3 = fits.Column(name='T', format='E', array=Teff)
    col4 = fits.Column(name='logg', format='E', array=log_g)
    col5 =  fits.Column(name='mass', format='E', array=initial_mass)
    col6 =  fits.Column(name='J', format='E', array=_2Mass_J)
    col7 =  fits.Column(name='H', format='E', array=_2Mass_H)
    col8 =  fits.Column(name='K', format='E', array=_2Mass_Ks)
    col9 =  fits.Column(name='stage', format='K', array=evo_stage) #'K' --> 64-bit integer
    col10 = fits.Column(name='weight', format='E', array=ones)
    cols = fits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8, col9, col10])
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
    x = 5
    columns = np.loadtxt(inpath + file1)
    initial_mass.extend(columns[:, 0])
    log_Teff.extend(columns[:, 1])
    log_g.extend(columns[:, 3])
    _2Mass_J.extend(columns[:, 9])
    _2Mass_H.extend(columns[:, 10])
    _2Mass_Ks.extend(columns[:, 11])
    Yale_file = open(inpath +file1, 'r')
    metallicity = find_Metallicity(Yale_file)
    for item in columns:
        metallicitylist.append(metallicity)
    Yale_file = open(inpath +file1, 'r')
    listOfAges = getAge(Yale_file, x)
    Yale_file = open(inpath +file1, 'r')
    Ages = CompiledAges(listOfAges, initial_mass, columns)
    allAges.extend(Ages)
    print file1 + '\t' + str(metallicity)
    i += 1

#print listOfAges

allAges = np.array(allAges).dot(1000000000)
allAges = np.log10(allAges)
listOfAges = np.array(listOfAges).dot(1000000000)
listOfAges = np.log10(listOfAges)
Teff.extend(np.power(10, log_Teff))
fitsTable(allAges, metallicitylist, initial_mass, Teff, log_g, _2Mass_J, _2Mass_H, _2Mass_Ks, outpath, listOfAges)
print 'FITS file created for Yale'

