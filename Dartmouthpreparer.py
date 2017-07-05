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


def fitsTable(allAges, metallicitylist, initial_mass, Teff, log_g, _2Mass_J, _2Mass_H, _2Mass_Ks, listOfAges, outpath):
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
    col1 = fits.Column(name='age', format='E', array=x)
    col2 = fits.Column(name='feh', format='E', array=metallicityList)
    col3 = fits.Column(name='T', format='E', array=Teff)
    col4 = fits.Column(name='logg', format='E', array=log_g)
    col5 =  fits.Column(name='mass', format='E', array=initial_mass)
    col6 =  fits.Column(name='J', format='E', array=_2Mass_J)
    col7 =  fits.Column(name='H', format='E', array=_2Mass_H)
    col8 =  fits.Column(name='K', format='E', array=_2Mass_Ks)
    col9 = fits.Column(name='stage', format='K', array=evo_stage)
    col10 = fits.Column(name='weight', format='E', array=ones)
    cols = fits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8, col9, col10])
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
listOfAges = np.array(listOfAges).dot(1000000000)
listOfAges = np.log10(listOfAges)
Teff.extend(np.power(10, log_Teff))
#print len(Teff)
fitsTable(allAges, metallicitylist, initial_mass, Teff, log_g, _2Mass_J, _2Mass_H, _2Mass_Ks, listOfAges, outpath)

print 'FITS file created for Dartmouth'
