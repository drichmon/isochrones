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
            
def log_teff_to_teff(log_Teff):
    i = 0
    for item in log_Teff:
        item = 10 ** item
        np.put(log_Teff, i, item)
        i +=1
    return log_Teff


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
                            Ages.append(Age)
                        else:
                            Ages.append(line[10:16])
                    except (ValueError):
                        pass
    return Ages
    Yale_file.close()

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

def fitsTable(listOfAges, metallicity, initial_mass, Teff, log_g, _2Mass_J, _2Mass_H, _2Mass_Ks, Yale_file, outpath):
    i = 0
    j = 1
    k = 0
    metallicityList = []
    Ages = listOfAges
    finalAges = []
    with Yale_file as infile:
        for line in infile:
            try:
                if len(line.strip()) == 0:
                    #print ''
                    k += 1
                elif len(line.strip()) != 0:
                    if '#' in line:
                        x = '##age(Gyr)=' #do something pointless
                        #print x + str(nextAge)
                    elif '#' not in line:
                        nextAge = Ages[k]
                        finalAges.append(float(nextAge))
                        metallicityList.append(metallicity)
                        i += 1
                else:
                    break
            except (IndexError):
                pass
        Yale_file.close()
    x = np.array(finalAges)
    metallicityArray = np.array(metallicityList)
    col1 = fits.Column(name='logAge', format='E', array=x)
    col2 = fits.Column(name='Metallicity', format='E', array=metallicityList)
    col3 = fits.Column(name='Teff', format='E', array=Teff)
    col4 = fits.Column(name='Log G', format='E', array=log_g)
    col5 =  fits.Column(name='Mass', format='E', array=initial_mass)
    col6 =  fits.Column(name='J', format='E', array=_2Mass_J)
    col7 =  fits.Column(name='H', format='E', array=_2Mass_H)
    col8 =  fits.Column(name='Ks', format='E', array=_2Mass_Ks)
    cols = fits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8])
    tbhdu = fits.BinTableHDU.from_columns(cols)
    tbhdu.writeto(outpath + 'Yale_%s.fits' % str(metallicity), clobber=True) #fits outfile
    
 
inputAge = str(0.150)# raw_input('')

for file in sorted(os.listdir(inpath)):
    file1 = file
    Yale_file = open(inpath +file1, 'r')
    num_skip_rows = get_skiprows(inputAge, Yale_file)
    x = num_skip_rows -1
    initial_mass = np.loadtxt(inpath +file1, skiprows = x, usecols = (0,))
    log_Teff = np.loadtxt(inpath +file1, skiprows = x, usecols = (1,))
    log_g = np.loadtxt(inpath +file1, skiprows = x, usecols = (3,))
    _2Mass_J = np.loadtxt(inpath +file1, skiprows = x, usecols = (9,))
    _2Mass_H = np.loadtxt(inpath +file1, skiprows = x, usecols = (10,))
    _2Mass_Ks = np.loadtxt(inpath +file1, skiprows = x, usecols = (11,))
    Teff = log_teff_to_teff(log_Teff)
    Yale_file = open(inpath +file1, 'r')
    metallicity = find_Metallicity(Yale_file)
    Yale_file = open(inpath +file1, 'r')
    listOfAges = getAge(Yale_file, x)
    Yale_file = open(inpath +file1, 'r')
    print file1 + '\t' + str(metallicity)
    fitsTable(listOfAges, metallicity, initial_mass, Teff, log_g, _2Mass_J, _2Mass_H, _2Mass_Ks, Yale_file, outpath)
