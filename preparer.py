import numpy as np
from astropy.io import fits

def _cdf_lowmass_Kroupa(mass):
    a_1 = 0.3; A = 1.785;
    return A * np.power(mass, 1-a_1) / (1 - a_1)

def _cdf_midmass_Kroupa(mass):
    m_0 = 0.08; a_2 = 1.3; A = 1.785; B = 0.4352; k_1 = 0.08;
    return B + A * k_1 / (1-a_2) * (np.power(mass, 1-a_2) - np.power(m_0, 1-a_2))
    
def _cdf_highmass_Kroupa(mass):
    m_1 = 0.5; a_3 = 2.3; A = 1.785; C = 0.86469; k_2 = 0.04;
    return C + A * k_2 / (1-a_3) * (np.power(mass, 1-a_3) - np.power(m_1, 1-a_3))

def cdf_mass(mass):
    if mass < 0.08:
        return _cdf_lowmass_Kroupa(mass)
    elif mass < 0.5:
        return _cdf_midmass_Kroupa(mass)
    else:
        return _cdf_highmass_Kroupa(mass)


def fitsTable(age_weights, logAge, metallicitylist, Teff, log_g, initial_mass, _2Mass_J, _2Mass_H, _2Mass_Ks, evo_stage, mass_weights, listOfAges, outpath, filetype):
    final_weights = np.array(age_weights) * np.array(mass_weights)
    col1 = fits.Column(name='age', format='E', array=logAge)
    col2 = fits.Column(name='feh', format='E', array=metallicitylist)
    col3 = fits.Column(name='T', format='E', array=Teff)
    col4 = fits.Column(name='logg', format='E', array=log_g)
    col5 =  fits.Column(name='mass', format='E', array=initial_mass)
    col6 =  fits.Column(name='J', format='E', array=_2Mass_J)
    col7 =  fits.Column(name='H', format='E', array=_2Mass_H)
    col8 =  fits.Column(name='K', format='E', array=_2Mass_Ks)
    col9 =  fits.Column(name='stage', format='K', array=evo_stage)
    col10 =  fits.Column(name='age weights', format='E', array=age_weights)
    col11 =  fits.Column(name='mass weights', format='E', array=mass_weights)
    col12 =  fits.Column(name='weight', format='E', array=final_weights)
    cols = fits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12])
    tbhdu = fits.BinTableHDU.from_columns(cols)
    ageGridhdu = fits.BinTableHDU.from_columns([fits.Column(name='Age Grid', format='E', array=listOfAges)])
    agelisthdr = fits.Header()
    agelisthdr['Ages'] = 'Age grid'
    agelisthdu = fits.PrimaryHDU(header=agelisthdr)
    tbhdulist = fits.HDUList([agelisthdu, tbhdu, ageGridhdu])
    if filetype == 'MIST':
        tbhdulist.writeto(outpath + 'MISTtest.fits', clobber=True) #fits outfile
    elif filetype == 'Dartmouth':
        tbhdulist.writeto(outpath + 'Dartmouthtest.fits', clobber=True) #fits outfile
    elif filetype == 'Yale':
        tbhdulist.writeto(outpath + 'Yaletest.fits', clobber=True) #fits outfile
    elif filetype == 'PARSEC':
        tbhdulist.writeto(outpath + 'PARSECtest.fits', clobber=True) #fits outfile
    print 'FITS file created for %s' % filetype

