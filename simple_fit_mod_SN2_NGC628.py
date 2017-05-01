import numpy as np
import pylab as pl

import pyfits
import orb.fit
import orb.utils.spectrum
import orb.utils.io
from orb.core import Lines

# PARAMETERS
order = 6 # In cube header
step = 1644.511699# In cube header
axis_corr = 1.03758024972 # In cube header
shift_guess = 700 # mean velocity guesss in km/s
velocity_range = 200
lines_nm = Lines().get_line_nm(['Hbeta','[OIII]4959','[OIII]5007']) # lines to fit
print lines_nm
cov_pos = [1,1,1] # groups of velocity, [1,2,1] means that the first and third lines share the same velocity
cov_fwhm = [1,2,2]
#signal_range = [483,512] # range of the signal to fit (in nm)
fwhm_guess = [0.9,0.94,0.94] # FWHM guess of the lines in nm
	
# computed params
nm_laser = 543.5
nm_laser_obs = nm_laser*axis_corr

# Read HII Regions files
mask = pyfits.getdata('NGC628_HIIPhot3_All_GROW_SN2.fits')
dist = pyfits.getdata('NGC628_distance_GROW_SN2.fits') 
sig = np.loadtxt('NGC628_sigma.txt')
cube = pyfits.getdata('NGC628_SN2.merged.nm.1.0_conv3_ciel_pop_tout.fits')
hdr = pyfits.open('NGC628_SN2.merged.nm.1.0_conv3_ciel_pop_tout.fits')
h = hdr[0].header

dump_file = np.zeros([2064,2048])
spectrum_nm = np.zeros([134,4285])
sp = np.zeros([134])
spectrum_fit_cube = np.zeros([134,4285])
spectrum_fit = np.zeros([134])
vel_hb = np.zeros([4285])
vel_hb_err = np.zeros([4285])
vel_oiii = np.zeros([4285])
vel_oiii_err = np.zeros([4285])
amp_hb = np.zeros([4285])
amp_oiii1 = np.zeros([4285])
amp_oiii2 = np.zeros([4285])
amp_hei = np.zeros([4285])
amp_sii1 = np.zeros([4285])
amp_sii2 = np.zeros([4285])
fwhm_hb = np.zeros([4285])
fwhm_oiii1 = np.zeros([4285])
fwhm_oiii2 = np.zeros([4285])
fwhm_hei = np.zeros([4285])
fwhm_sii1 = np.zeros([4285])
fwhm_sii2 = np.zeros([4285])
heig_hb = np.zeros([4285])
heig_oiii1 = np.zeros([4285])
heig_oiii2 = np.zeros([4285])
heig_hei = np.zeros([4285])
heig_sii1 = np.zeros([4285])
heig_sii2 = np.zeros([4285])
snr_hb = np.zeros([4285])
snr_oiii1 = np.zeros([4285])
snr_oiii2 = np.zeros([4285])
snr_hei = np.zeros([4285])
snr_sii1 = np.zeros([4285])
snr_sii2 = np.zeros([4285])
amp_hb_err = np.zeros([4285])
amp_oiii1_err = np.zeros([4285])
amp_oiii2_err = np.zeros([4285])
amp_hei_err = np.zeros([4285])
amp_sii1_err = np.zeros([4285])
amp_sii2_err = np.zeros([4285])
fwhm_hb_err = np.zeros([4285])
fwhm_oiii1_err = np.zeros([4285])
fwhm_oiii2_err = np.zeros([4285])
fwhm_hei_err = np.zeros([4285])
fwhm_sii1_err = np.zeros([4285])
fwhm_sii2_err = np.zeros([4285])
heig_hb_err = np.zeros([4285])
heig_oiii1_err = np.zeros([4285])
heig_oiii2_err = np.zeros([4285])
heig_hei_err = np.zeros([4285])
heig_sii1_err = np.zeros([4285])
heig_sii2_err = np.zeros([4285])

#indicex = []
#indicey = []

#for i in range(1,4286):
#    indicex.append(np.where((mask == i) & (dist < sig[i-1]))[0])
#    indicey.append(np.where((mask == i) & (dist < sig[i-1]))[1])

#with open("indicex_regions.txt", 'w') as f:
#    for i in indicex:
#        f.write(str(i) + '\n')
#
#with open("indicey_regions.txt", 'w') as f:
#    for i in indicey:
#        f.write(str(i) + '\n')
#
#with open("indicex_regions.txt", 'r') as f:
#    indicex = [line.rstrip('\n') for line in f]
#with open("indicey_regions.txt", 'r') as f:
#    indicey = [line.rstrip('\n') for line in f]

#for j in range(0,134):
#    dump_file[0:2064,0:2048] = cube[j,0:2064,0:2048] 	
#    print j

#    for i in range(0,4285):
#        spectrum_nm[j,i] = np.sum(dump_file[indicex[i],indicey[i]])
   
#pyfits.writeto('NGC628_SN2.merged.nm.1.0_conv3_ciel_pop_tout_regions.fits', spectrum_nm, h)    
spectrum_nm = pyfits.getdata('NGC628_SN2.merged.nm.1.0_conv3_ciel_pop_tout_regions.fits')

for i in range(0,4285):
    sp[0:134] = spectrum_nm[0:134,i]
    # fit loop
    
    fit_results = orb.fit.fit_lines_in_spectrum(
        sp, lines_nm, step, order, nm_laser,
        nm_laser_obs=nm_laser_obs, wavenumber=False,
        fwhm_guess=fwhm_guess, cont_guess=None,
        shift_guess=shift_guess, sigma_guess=0.,
        fix_fwhm=False, cov_fwhm=cov_fwhm, cov_pos=cov_pos,
        fix_pos=False, cov_sigma=False,
        fit_tol=1e-10, poly_order=0,
        fmodel='sinc', signal_range=None,
        filter_file_path=None, fix_filter=False,
        apodization=1., velocity_range=velocity_range,
        compute_mcmc_error=False, no_error=False)
        
    

    print i
    if fit_results != []: # if the fit is ok 
	spectrum_fit = fit_results['fitted-vector']
    
        for j in range(0,134):		
            spectrum_fit_cube[j,i] = float(spectrum_fit[j])

    	try :
            vel_hb[i] = fit_results['velocity'][0]
            vel_oiii[i] = fit_results['velocity'][2]
            amp_hb[i] = fit_results['lines-params'][0,1]
            amp_oiii1[i] = fit_results['lines-params'][1,1]
            amp_oiii2[i] = fit_results['lines-params'][2,1]
            fwhm_hb[i] = fit_results['lines-params'][0,3]
            fwhm_oiii1[i] = fit_results['lines-params'][1,3]
            fwhm_oiii2[i] = fit_results['lines-params'][2,3]
            heig_hb[i] = fit_results['lines-params-err'][0,0]
            heig_oiii1[i] = fit_results['lines-params-err'][1,0]
            heig_oiii2[i] = fit_results['lines-params-err'][2,0]

        except :
	    pass
            
    		
	try :

            vel_hb_err[i] = fit_results['velocity-err'][0]
            vel_oiii_err[i] = fit_results['velocity-err'][2]
            amp_hb_err[i] = fit_results['lines-params-err'][0,1]
            amp_oiii1_err[i] = fit_results['lines-params-err'][1,1]
            amp_oiii2_err[i] = fit_results['lines-params-err'][2,1]
            fwhm_hb_err[i] = fit_results['lines-params-err'][0,3]
            fwhm_oiii1_err[i] = fit_results['lines-params-err'][1,3]
            fwhm_oiii2_err[i] = fit_results['lines-params-err'][2,3]
            heig_hb_err[i] = fit_results['lines-params-err'][0,0]
            heig_oiii1_err[i] = fit_results['lines-params-err'][1,0]
            heig_oiii2_err[i] = fit_results['lines-params-err'][2,0]
       
   	except :
	    pass
        
# Write results
    
pyfits.writeto('NGC628_SN2.regions.Hb4861_velocity.fits', vel_hb, h)    
pyfits.writeto('NGC628_SN2.regions.OIII5007_velocity.fits', vel_oiii, h)     
pyfits.writeto('NGC628_SN2.regions.Hb4861_velocity_error.fits', vel_hb_err, h)    
pyfits.writeto('NGC628_SN2.regions.OIII5007_velocity_error.fits', vel_oiii_err, h)
pyfits.writeto('NGC628_SN2.regions.fit.spectrum.fits', spectrum_fit_cube, h)
    
pyfits.writeto('NGC628_SN2.regions.Hb4861_amplitude.fits', amp_hb, h)    
pyfits.writeto('NGC628_SN2.regions.OIII4959_amplitude.fits', amp_oiii1, h)    
pyfits.writeto('NGC628_SN2.regions.OIII5007_amplitude.fits', amp_oiii2, h)    
pyfits.writeto('NGC628_SN2.regions.Hb4861_fwhm.fits', fwhm_hb, h)    
pyfits.writeto('NGC628_SN2.regions.OIII4959_fwhm.fits', fwhm_oiii1, h)    
pyfits.writeto('NGC628_SN2.regions.OIII5007_fwhm.fits', fwhm_oiii2, h)    
pyfits.writeto('NGC628_SN2.regions.Hb4861_heigt.fits', heig_hb, h)    
pyfits.writeto('NGC628_SN2.regions.OIII4959_heigt.fits', heig_oiii1, h)    
pyfits.writeto('NGC628_SN2.regions.OIII5007_heigt.fits', heig_oiii2, h)     
    
pyfits.writeto('NGC628_SN2.regions.Hb4861_amplitude_error.fits', amp_hb_err, h)    
pyfits.writeto('NGC628_SN2.regions.OIII4959_amplitude_error.fits', amp_oiii1_err, h)    
pyfits.writeto('NGC628_SN2.regions.OIII5007_amplitude_error.fits', amp_oiii2_err, h)    
pyfits.writeto('NGC628_SN2.regions.Hb4861_fwhm_error.fits', fwhm_hb_err, h)    
pyfits.writeto('NGC628_SN2.regions.OIII4959_fwhm_error.fits', fwhm_oiii1_err, h)    
pyfits.writeto('NGC628_SN2.regions.OIII5007_fwhm_error.fits', fwhm_oiii2_err, h)      
pyfits.writeto('NGC628_SN2.regions.Hb4861_heigt_error.fits', heig_hb_err, h)    
pyfits.writeto('NGC628_SN2.regions.OIII4959_heigt_error.fits', heig_oiii1_err, h)    
pyfits.writeto('NGC628_SN2.regions.OIII5007_heigt_error.fits', heig_oiii2_err, h)         