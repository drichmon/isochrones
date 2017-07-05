import sys
import os
import numpy as np
from astropy.io import fits

inpath = '/home/richmond/Isochrones/PARSEC/' #path to isochrone files
outpath = '/home/richmond/Isochrones/FITS/' #desired output for FITS files

def find_Metallicity(PARSEC_file):
    i = 0
    with PARSEC_file as data:
        for line in data:
            if '[M/H]' in line:
                metallicity = line[45:51]
    PARSEC_file.close()
    return metallicity
            
def log_teff_to_teff(log_Teff):
    i = 0
    for item in log_Teff:
        item = 10 ** item
        np.put(log_Teff, i, item)
        i +=1
    return log_Teff
        
def getAge(PARSEC_file, x):
    Ages = []
    i = 0
    term = "Age ="
    with PARSEC_file as infile:
        for line in infile:
            if term in line:
                #print line[72:78] + '\t' + line[80:83]
                Age = line[72:83]#float(line[72:78]) * 10 ** int(line[80:83]) #PREVIOUSLY +1 TO ALL INDICIES IN THIS LINE
                Ages.append(float(Age))
    return Ages
    PARSEC_file.close()

def get_skiprows(inputAge, PARSEC_file):
    i = 1
    term = "Age = \t" + str(inputAge)
    for line in PARSEC_file:
        if term not in line:
            i += 1
        elif term in line:
            i += 1
            return i
    PARSEC_file.close()

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
                #print nextAge
                finalAges.append(nextAge)
                #metallicityList.append(metallicity)
                i += 1
            elif columns[j, 1] == columns[i, 1]:
                nextAge = listOfAges[k]
                j += 1
                #print nextAge
                finalAges.append(nextAge)
                #metallicityList.append(metallicity)
                i += 1
            else:
                finalAges.append(nextAge)
                break
        except (IndexError):
            pass    
    return finalAges

def fitsTable(allAges, metallicitylist, initial_mass, Teff, log_g, _2Mass_J, _2Mass_H, _2Mass_Ks, evo_stage, outpath, listOfAges):
    metallicityList = np.array(metallicitylist)  
    x = np.array(allAges) 
    initial_mass = np.array(initial_mass)
    Teff = np.array(Teff)
    log_g = np.array(log_g)
    _2Mass_J = np.array(_2Mass_J)
    _2Mass_H = np.array(_2Mass_H)
    _2Mass_Ks = np.array(_2Mass_Ks)
    metallicityArray = np.array(metallicityList)
    evo_stage = np.array(evo_stage)
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
    col10 =  fits.Column(name='weight', format='E', array=ones)
    cols = fits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8, col9, col10])
    tbhdu = fits.BinTableHDU.from_columns(cols)
    ageGridhdu = fits.BinTableHDU.from_columns([fits.Column(name='Age Grid', format='E', array=listOfAges)])
    agelisthdr = fits.Header()
    agelisthdr['Ages'] = 'Age grid' #don't really know what I'm doing here/if this is necessary
    agelisthdu = fits.PrimaryHDU(header=agelisthdr)
    tbhdulist = fits.HDUList([agelisthdu, tbhdu, ageGridhdu])
    tbhdulist.writeto(outpath + 'PARSEC.fits', clobber=True) #fits outfile

initial_mass = [] 
log_Teff = []
log_g = []
_2Mass_J = []
_2Mass_H = []
_2Mass_Ks = []
Teff = []
metallicitylist = []
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
            metallicitylist.append(metallicity)
        print file + '\t' + metallicity

listOfAges = getAge(open(inpath + 'bv_iso_0.000100.dat', 'r'), x) #not necessary to loop this      

Teff.extend(np.power(10, log_Teff))
fitsTable(allAges, metallicitylist, initial_mass, Teff, log_g, _2Mass_J, _2Mass_H, _2Mass_Ks, evo_stage, outpath, listOfAges)

print 'FITS file created for PARSEC'
