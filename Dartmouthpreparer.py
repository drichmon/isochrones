import sys
import os
import numpy as np
import argparse
from astropy.io import fits


inpath = '/home/richmond/Isochrones/Dartmouth/interpolated/' #path to isochrone files
outpath = '/home/richmond/Isochrones/FITS/' #desired output for FITS files

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
                    Ages.append(line[5:12])
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

def fitsTable(listOfAges, metallicity, initial_mass, Teff, log_g, _2Mass_J, _2Mass_H, _2Mass_Ks, outpath):
    i = 0
    j = 1
    k = 0
    metallicityList = []
    Ages = listOfAges
    finalAges = []
    for item in initial_mass:
        np.set_printoptions(precision=14)
        try:
            if EEP[j] < EEP[i]:
                nextAge = Ages[k]
                k += 1
                j += 1
                finalAges.append(nextAge)
                metallicityList.append(metallicity)
                i += 1
            elif EEP[j] > EEP[i]:
                nextAge = Ages[k]
                j += 1
                finalAges.append(nextAge)
                metallicityList.append(metallicity)
                i += 1
            else:
                break
        except (IndexError):
            pass    
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
    tbhdu.writeto(outpath + 'Dartmouth_%s.fits' % metallicity, clobber=True) #fits outfile


inputAge = str(1.000) #raw_input('')
for file in sorted(os.listdir(inpath)):
    file1 = file
    Dartmouth_file = open(inpath +file1, 'r')
    num_skip_rows = get_skiprows(inputAge, Dartmouth_file)
    x = num_skip_rows - 1
    EEP = np.loadtxt(inpath +file1, skiprows = x, usecols = (0,))
    initial_mass = np.loadtxt(inpath +file1, skiprows = x, usecols = (1,))
    log_Teff = np.loadtxt(inpath +file1, skiprows = x, usecols = (2,))
    log_g = np.loadtxt(inpath +file1, skiprows = x, usecols = (3,))
    _2Mass_J = np.loadtxt(inpath +file1, skiprows = x, usecols = (10,))
    _2Mass_H = np.loadtxt(inpath +file1, skiprows = x, usecols = (11,))
    _2Mass_Ks = np.loadtxt(inpath +file1, skiprows = x, usecols = (12,))
    Teff = log_teff_to_teff(log_Teff)
    Dartmouth_file = open(inpath +file1, 'r')
    metallicity = find_Metallicity(Dartmouth_file)
    Dartmouth_file = open(inpath +file1, 'r')
    listOfAges = getAge(Dartmouth_file, x)
    Dartmouth_file = open(inpath +file1, 'r')
    print file1 + '\t' + metallicity
    fitsTable(listOfAges, metallicity, initial_mass, Teff, log_g, _2Mass_J, _2Mass_H, _2Mass_Ks, outpath)
