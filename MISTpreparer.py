import numpy as np
import scipy.signal
import re
import os
from collections import OrderedDict
from astropy.io import fits
import scipy.signal

inpath = '/home/richmond/Isochrones/MIST/MIST_v1.0_vvcrit0.4_UBVRIplus/' #path to isochrone files
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
'''
def findMin(log_g, endIndex):
    i = 0
    minList = []
    for index in endIndex:
        try:
            temp_log_g = log_g[endIndex[i]:endIndex[i+1]]
            #print temp_log_g
            numpy_temp_log_g = np.array(temp_log_g)
            minimum = np.argmin(numpy_temp_log_g)
            #print minimum
            minList.append(minimum)
            i += 1
        except (IndexError):
            pass
    #print minList
    return minList
'''
'''
def cut(minIndex, columns, endIndex, Ages):
    i = 0
    FinalcutColumns = []
    for key in sorted(Ages):
        x = []
        try:
            for line in Ages[key]:
                x.append(line)
                del x[minIndex[i]:len(Ages[key])]
                i += 1
        except (IndexError):
            pass
        FinalcutColumns.extend(x)
    FinalcutColumns = np.asarray(FinalcutColumns)
    np.savetxt('/home/richmond/extrastuff/MIST_%stemp.dat' % metallicity, FinalcutColumns, fmt='%.4e')

    j = 0
    i = 0
    FinalcutColumns = []
    for item in endIndex:
        try:
            cutColumns = []
            cutColumns.append(columns[endIndex[i]:endIndex[i+1]])
            cutColumns = np.array(cutColumns)
            np.delete(cutColumns, range(0, minIndex[i]))
            #print cutColumns
            cutColumns.tolist()
            x =  cutColumns[0][:]
            for item in x:
                FinalcutColumns.append(item)
            i += 1
        except (IndexError):
            pass
    FinalcutColumns = np.asarray(FinalcutColumns)
    #print FinalcutColumns
    np.savetxt('/home/richmond/Isochrones/temp/MIST_%stemp.dat' % metallicity, FinalcutColumns, fmt='%.4e')
'''

def fitsTable(logAge, metallicitylist, initial_mass, Teff, log_g, _2Mass_J, _2Mass_H, _2Mass_Ks, outpath):
    logAge = np.array(logAge)
    initial_mass = np.array(initial_mass)
    Teff = np.array(Teff)
    log_g = np.array(log_g)
    _2Mass_J = np.array(_2Mass_J)
    _2Mass_H = np.array(_2Mass_H)
    _2Mass_Ks = np.array(_2Mass_Ks)
    metallicitylist = np.array(metallicitylist)
    col1 = fits.Column(name='logAge', format='E', array=logAge)
    col2 = fits.Column(name='Metallicity', format='E', array=metallicitylist)
    col3 = fits.Column(name='Teff', format='E', array=Teff)
    col4 = fits.Column(name='Log G', format='E', array=log_g)
    col5 =  fits.Column(name='Mass', format='E', array=initial_mass)
    col6 =  fits.Column(name='J', format='E', array=_2Mass_J)
    col7 =  fits.Column(name='H', format='E', array=_2Mass_H)
    col8 =  fits.Column(name='Ks', format='E', array=_2Mass_Ks)
    cols = fits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8])
    tbhdu = fits.BinTableHDU.from_columns(cols)
    tbhdu.writeto(outpath + 'MIST.fits', clobber=True) #fits outfile


logAge = []
initial_mass = []
log_Teff = []
log_g = []
_2Mass_J = []
_2Mass_H = []
_2Mass_Ks = []
metallicitylist = []
Teff = []

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
            MIST_file = open(inpath +file, 'r')
            metallicity = find_Metallicity(MIST_file)
            for item in columns:
                metallicitylist.append(metallicity)
            print file + '\t' + str(metallicity)


Teff.extend(np.power(10, log_Teff)) # log_teff_to_teff(log_Teff)
fitsTable(logAge, metallicitylist, initial_mass, Teff, log_g, _2Mass_J, _2Mass_H, _2Mass_Ks, outpath)
print 'FITS file created for MIST'

