'''
def endIndexx(x, MIST_file):#group all the same ages by their range of indicies in the file, then find the minimum. After the minimum, skip lines until the next age starts
    i = 0
    endIndex = []
    with MIST_file as infile:
        for line in infile:
            if i >= x:
                if len(line.strip()) == 0:
                    endIndex.append(i)
                    i +=1
                else:
                    i += 1
            else:
                i += 1
    endIndex = endIndex[::2]
    endIndex.append(i)
    endIndex.append(x)
    return sorted(endIndex)
'''
'''
def group_Ages(Age_dict, MIST_file, endIndex, x):
    i = 0
    j = 0 
    k = 1
    lines = []
    dictionary = {}
    MIST_file = open(inpath +file1, 'r')
    with MIST_file as infile:
        for line in infile:
            try:
                if i >= (x-5):
                    if i == endIndex[k]:
                        k += 1 
                        i += 1
                        lines.append(line)
                        dictionary[Age_dict[j]] = lines[4:len(lines)-1]
                        j += 1
                        lines = []
                    elif i in range(endIndex[j], endIndex[k]):
                        lines.append(line)
                        i += 1
                    else:
                        lines.append(line)
                        i += 1
                else:
                    i += 1
            except (IndexError):
                pass
    MIST_file.close()
    return dictionary
'''
'''
def findMin(groupAges, log_g): #old method
    minIndex = []
    for item in groupAges:
        a = groupAges[item] 
        i = 0
        log_g = []
        for i in range(len(a)):
            log_g.append(a[i][141:165])
            i += 1
        x = []
        for item in log_g:
            if '+' in item:
                i = float(item[:item.index('E')]) * 10 ** int(item[item.index('+'):])
                x.append(i)
            elif '-' in item:
                i = float(item[:item.index('E')]) * 10 ** int(item[item.index('E')+1:])
                x.append(i)
        numpyarray = np.array(x)
        print len(numpyarray)
        minimum  = np.argmin(numpyarray) #scipy.signal.argrelmin(numpyarray, order=1000)[0][0]
        #print minimum
        minIndex.append(minimum)
    #print minIndex
    return minIndex
'''
'''
def cut(minIndex, groupAges): #old method
    dummylist = []
    f = open('/home/richmond/extrastuff/MIST_%stemp.dat' % metallicity, 'w')
    #i = 0
    print minIndex
    for i, item in enumerate(sorted(groupAges)):
        a = groupAges[item]
        dummylist.extend(a[:minIndex[i]+1])
        #for line in a:
        #    if j <= minIndex[i]:                
        #        dummylist.append(line)
        #        j += 1
        #i += 1 #very important that this stays in this indentation
    print 'dummy', len(dummylist)
    for line in dummylist:
        print >>f, line, 
    f.close()
'''
'''
def get_skiprows(inputAge, MIST_file):
    i = 0 
    with MIST_file as infile:
        for line in infile:
            if inputAge not in line:
                i += 1
            if inputAge in line:
                print i
                return i
    MIST_file.close()
'''
import numpy as np
import scipy.signal
import re
import os
from collections import OrderedDict
from astropy.io import fits
import scipy.signal

inpath = '/home/richmond/Isochrones/MIST/MIST_v1.0_vvcrit0.4_UBVRIplus/interpolated/' #path to isochrone files
outpath = '/home/richmond/extrastuff/' #desired output for FITS files

def find_Metallicity(MIST_file):
    i = 0
    with MIST_file as data:
        for line in data:
            if '[Fe/H]' in line:
                x = data.next()
                metallicity = x[25:31]
                break
    MIST_file.close()
    return metallicity.strip()
            
def log_teff_to_teff(log_Teff):
    i = 0
    for item in log_Teff:
        item = 10 ** item
        np.put(log_Teff, i, item)
        i +=1
    return log_Teff

def Age_dict(columns):
    age_dict = {}
    separation = np.where(columns[1:, 0] - columns[:-1, 0] < 0)
    old = 0
    for ind in separation[0]:
        age_dict[columns[ind, 1]] = columns[old:ind+1]
        old = ind + 1
    return age_dict


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

def cut(minIndex, columns, endIndex, Ages):
    '''i = 0
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
'''
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
    np.savetxt('/home/richmond/extrastuff/MIST_%stemp.dat' % metallicity, FinalcutColumns, fmt='%.4e')


def fitsTable(logAge, metallicitylist, initial_mass, Teff, log_g, _2Mass_J, _2Mass_H, _2Mass_Ks, outpath):
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



metallicitylist = []

for file1 in sorted(os.listdir(inpath)): #this loop can be commented out to save time if temp files are already created
    if '.cmd' in file1:
        if 'plus' in file1:
            columns = np.loadtxt(inpath + file1)
            log_g = columns[:, 4]
            logAge = columns[:, 1]
            Ages = Age_dict(columns) # list of the ages
            #print Ages
            MIST_file = open(inpath +file1, 'r')
            metallicity = find_Metallicity(MIST_file)
            metallicitylist.append(metallicity)
            endIndex = np.where(columns[1:, 0] - columns[:-1, 0] < 0)[0] #endIndexx(x, MIST_file)
            #print endIndex
            minIndex = findMin(log_g, endIndex)
            #print minIndex
            cut(minIndex, columns, endIndex, Ages)
            print 'temp file created for: ' + file1 + "\t" + metallicity
        else:
            columns = np.loadtxt(inpath + file1)
            log_g = columns[:, 5]
            logAge = columns[:, 1]
            Ages = Age_dict(columns)
            MIST_file = open(inpath +file1, 'r')
            metallicity = find_Metallicity(MIST_file)
            metallicitylist.append(metallicity)
            endIndex = np.where(columns[1:, 0] - columns[:-1, 0] < 0)[0] #endIndexx(x, MIST_file)
            minIndex = findMin(log_g, endIndex)
            #print minIndex
            cut(minIndex, columns, endIndex, Ages)
            print 'temp file created for: ' + file1 + "\t" + metallicity

newinpath = '/home/richmond/extrastuff/' #location of temp files

logAge = []
initial_mass = []
log_Teff = []
log_g = []
_2Mass_J = []
_2Mass_H = []
_2Mass_Ks = []
metallicityarray = []
met_index = 0
for file in sorted(os.listdir(newinpath)):
    if '.dat' in file: # makes sure it only iterates over isochrone files
        columns = np.loadtxt(newinpath + file)
        if columns.shape[1] == 24:
            for i in columns:
                logAge.append(i[1])
                initial_mass.append(i[2])
                log_Teff.append(i[3])
                log_g.append(i[4])
                _2Mass_J.append(i[12])
                _2Mass_H.append(i[13])
                _2Mass_Ks.append(i[14])
                metallicityarray.append(metallicitylist[met_index])
        else:    
            for i in columns:
                logAge.append(i[1])
                initial_mass.append(i[2])
                log_Teff.append(i[4])
                log_g.append(i[5])
                _2Mass_J.append(i[14])
                _2Mass_H.append(i[15])
                _2Mass_Ks.append(i[16])
                metallicityarray.append(metallicitylist[met_index])
        met_index += 1

logAge = np.array(logAge)
initial_mass = np.array(initial_mass)
log_Teff = np.array(log_Teff)
log_g = np.array(log_g)
_2Mass_J = np.array(_2Mass_J)
_2Mass_H = np.array(_2Mass_H)
_2Mass_Ks = np.array(_2Mass_Ks)
metallicitylist = np.array(metallicityarray)
Teff = np.power(10, log_Teff) # log_teff_to_teff(log_Teff)
fitsTable(logAge, metallicitylist, initial_mass, Teff, log_g, _2Mass_J, _2Mass_H, _2Mass_Ks, outpath)
print 'FITS file created for MIST'


