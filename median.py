from astropy.io import fits
import numpy as np
import pylab as plt
from scipy import stats

from unidam.utils import mathematics

"""This code is meant to find the median of differences between calculations made by UniDAM for distance, age, and mass. Plots are also created to investigate if the amount of age discrepancy is related to a star's location on the HR diagram. The fits files opened by this script are the results from running UniDAM on a set of data."""

def values(param):
    x = file1_cols[param][np.in1d(file1_id, in_both)]
    y = file2_cols[param][np.in1d(file2_id, in_both)]
    return x, y

def plot(var1, var2):
    hr = plt.scatter(var1, var2, c=age_mean_differences, s=9, picker=True, edgecolors='none')
    cbar = plt.colorbar()
    cbar.set_label('Difference in Age Estimates between {0} and {1}'.format(code1, code2)) 
    plt.title('{0} observations'.format(dataset), color='k')
    ax = plt.gca()
    ax.set_ylim(ax.get_ylim()[::-1])
    ax.set_xlim(ax.get_xlim()[::-1])
    if var2 is logg:
        #h = stats.binned_statistic_2d(var1, var2, age_mean_differences, statistic='median', bins=20)
        #x, y = np.meshgrid(var1, var2)
        #plt.pcolormesh(x,y,h[0])
        plt.xlabel('Temperature')
        plt.ylabel('logg')
        plt.savefig(outpath + '{0}_{1}_{2}_{3}v{4}.png'.format(dataset, code1, 
                                                               code2, 'temp', 
                                                               'logg'))
    else:
        #h = stats.binned_statistic_2d(var1, var2, age_mean_differences, statistic='median', bins=20)
        #x, y = np.meshgrid(var1, var2)
        #plt.pcolormesh(x,y,h[0])
        plt.xlabel('Temperature')
        plt.ylabel('Fe/H')
        plt.savefig(outpath + '{0}_{1}_{2}_{3}v{4}.png'.format(dataset, code1, 
                                                               code2, 'temp', 
                                                               'feh'))
    plt.clf()

dataset = 'GAIA_ESO'
code1 = 'PARSEC'
code2 = 'Yale'

if dataset == 'APOKASC':
    file1 = fits.open('/home/richmond/Isochrones/results/APOKASC_results/%sresults.fits' % code1)[1]
    file2 = fits.open('/home/richmond/Isochrones/results/APOKASC_results/%sresults.fits' % code2)[1]
    file3 = fits.open('/home/richmond/Isochrones/results/APOKASC_results/APOKASC_2014.fits')[1]
    outpath = '/home/richmond/Isochrones/plots/%s/' % dataset
elif dataset == 'GAIA_ESO':
    file1 = fits.open('/home/richmond/Isochrones/results/GAIA_results/GAIA_%sresults.fits' % code1)[1]
    file2 = fits.open('/home/richmond/Isochrones/results/GAIA_results/GAIA_%sresults.fits' % code2)[1]
    file3 = fits.open('/home/richmond/Isochrones/results/GAIA_results/GAIA_ESO.fits')[1]
    outpath = '/home/richmond/Isochrones/plots/%s/' % dataset

file1_cols = file1.data
file2_cols = file2.data
file3_cols = file3.data
file1_cols = file1_cols[file1_cols['uspdf_priority']==0]
file2_cols = file2_cols[file2_cols['uspdf_priority']==0]
file1_id = file1_cols['id']
file2_id = file2_cols['id']
file3_id = file3_cols['id']

if len(file1_id) < len(file2_id): #find the common id's where fits were found
    in_both = np.intersect1d(file1_id, file2_id)
else:
    in_both = np.intersect1d(file2_id, file1_id)

file1_age, file2_age = values('age_mean')
file1_mass, file2_mass = values('mass_mean')
file1_distance, file2_distance = values('distance_modulus_mean')
age_mean_differences = np.subtract(file1_age, file2_age)
mass_mean_differences = np.subtract(file1_mass, file2_mass)
distance_mean_differences = np.subtract(file1_distance, file2_distance)
age_median, age_mad = mathematics.median_mad(age_mean_differences)
mass_median, mass_mad = mathematics.median_mad(mass_mean_differences)
distance_median, distance_mad = mathematics.median_mad(distance_mean_differences)
print age_median, age_mad
print mass_median, mass_mad
print distance_median, distance_mad
temperature = file3_cols['T'][np.in1d(file3_id, in_both)]
logg = file3_cols['logg'][np.in1d(file3_id, in_both)]
feh = file3_cols['feh'][np.in1d(file3_id, in_both)]
plot(temperature, logg)
plot(temperature, feh)
