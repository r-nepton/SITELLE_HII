import numpy as np
import pylab as pl

import pyfits
import orb.fit
import orb.utils.spectrum
import orb.utils.io
from orb.core import Lines

# PARAMETERS
order = 8 # In cube header
step = 1612.269411# In cube header
axis_corr = 1.03664546289 # In cube header
shift_guess = 1000 # mean velocity guesss in km/s
velocity_range = 200
lines_nm = [372.7] #Lines().get_line_nm(['[OII]3726']) # lines to fit
print lines_nm
cov_pos = [1] # groups of velocity, [1,2,1] means that the first and third lines share the same velocity
cov_fwhm = [1]
#signal_range = [365,385] # range of the signal to fit (in nm)
fwhm_guess = [0.67] # FWHM guess of the lines in nm
	
# computed params
nm_laser = 543.5
nm_laser_obs = nm_laser*axis_corr

# Read HII Regions files
mask = pyfits.getdata('NGC628_HIIPhot3_All_GROW_SN1.fits')
dist = pyfits.getdata('NGC628_distance_GROW_SN1.fits') 
sig = np.loadtxt('NGC628_sigma.txt')
cube = pyfits.getdata('NGC628_SN1.merged.nm.1.0_conv3_ciel_pop_tout2.fits')
hdr = pyfits.open('NGC628_SN1.merged.nm.1.0_conv3_ciel_pop_tout2.fits')
h = hdr[0].header

dump_file = np.zeros([2064,2048])
spectrum_nm = np.zeros([105,4285])
sp = np.zeros([105])
spectrum_fit_cube = np.zeros([105,4285])
spectrum_fit = np.zeros([105])
vel_oii = np.zeros([4285])
vel_oii_err = np.zeros([4285])
vel_oiii = np.zeros([4285])
vel_oiii_err = np.zeros([4285])
amp_oii = np.zeros([4285])
amp_oiii1 = np.zeros([4285])
amp_oiii2 = np.zeros([4285])
amp_hei = np.zeros([4285])
amp_sii1 = np.zeros([4285])
amp_sii2 = np.zeros([4285])
fwhm_oii = np.zeros([4285])
fwhm_oiii1 = np.zeros([4285])
fwhm_oiii2 = np.zeros([4285])
fwhm_hei = np.zeros([4285])
fwhm_sii1 = np.zeros([4285])
fwhm_sii2 = np.zeros([4285])
heig_oii = np.zeros([4285])
heig_oiii1 = np.zeros([4285])
heig_oiii2 = np.zeros([4285])
heig_hei = np.zeros([4285])
heig_sii1 = np.zeros([4285])
heig_sii2 = np.zeros([4285])
snr_oii = np.zeros([4285])
snr_oiii1 = np.zeros([4285])
snr_oiii2 = np.zeros([4285])
snr_hei = np.zeros([4285])
snr_sii1 = np.zeros([4285])
snr_sii2 = np.zeros([4285])
amp_oii_err = np.zeros([4285])
amp_oiii1_err = np.zeros([4285])
amp_oiii2_err = np.zeros([4285])
amp_hei_err = np.zeros([4285])
amp_sii1_err = np.zeros([4285])
amp_sii2_err = np.zeros([4285])
fwhm_oii_err = np.zeros([4285])
fwhm_oiii1_err = np.zeros([4285])
fwhm_oiii2_err = np.zeros([4285])
fwhm_hei_err = np.zeros([4285])
fwhm_sii1_err = np.zeros([4285])
fwhm_sii2_err = np.zeros([4285])
heig_oii_err = np.zeros([4285])
heig_oiii1_err = np.zeros([4285])
heig_oiii2_err = np.zeros([4285])
heig_hei_err = np.zeros([4285])
heig_sii1_err = np.zeros([4285])
heig_sii2_err = np.zeros([4285])

#indicex = []
#indicey = []
#
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

#for j in range(0,105):
#    dump_file[0:2064,0:2048] = cube[j,0:2064,0:2048] 	
#    print j
#
#    for i in range(0,4285):
#        spectrum_nm[j,i] = np.sum(dump_file[indicex[i],indicey[i]])
  
#pyfits.writeto('NGC628_SN1.merged.nm.1.0_conv3_ciel_pop_tout_regions.fits', spectrum_nm, h)    
spectrum_nm = pyfits.getdata('NGC628_SN1.merged.nm.1.0_conv3_ciel_pop_tout_regions.fits')

for i in range(0,4285):
    sp[0:105] = spectrum_nm[0:105,i]
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
    
        for j in range(0,105):		
            spectrum_fit_cube[j,i] = float(spectrum_fit[j])

    	try :
            vel_oii[i] = fit_results['velocity'][0]
            amp_oii[i] = fit_results['lines-params'][0,1]
            fwhm_oii[i] = fit_results['lines-params'][0,3]
            heig_oii[i] = fit_results['lines-params-err'][0,0]

        except :
	    pass
            
	try :

            vel_oii_err[i] = fit_results['velocity-err'][0]
            amp_oii_err[i] = fit_results['lines-params-err'][0,1]
            fwhm_oii_err[i] = fit_results['lines-params-err'][0,3]
            heig_oii_err[i] = fit_results['lines-params-err'][0,0]
     
        except :
	    pass
     
# Write results
    
pyfits.writeto('NGC628_SN1.regions.OII3727_velocity.fits', vel_oii, h)    
pyfits.writeto('NGC628_SN1.regions.OII3727_velocity_error.fits', vel_oii_err, h)    
pyfits.writeto('NGC628_SN1.regions.fit.spectrum.fits', spectrum_fit_cube, h)
    
pyfits.writeto('NGC628_SN1.regions.OII3727_amplitude.fits', amp_oii, h)      
pyfits.writeto('NGC628_SN1.regions.OII3727_fwhm.fits', fwhm_oii, h)      
pyfits.writeto('NGC628_SN1.regions.OII3727_heigt.fits', heig_oii, h)     
    
pyfits.writeto('NGC628_SN1.regions.OII3727_amplitude_error.fits', amp_oii_err, h)    
pyfits.writeto('NGC628_SN1.regions.OII3727_fwhm_error.fits', fwhm_oii_err, h)    
pyfits.writeto('NGC628_SN1.regions.OII3727_heigt_error.fits', heig_oii_err, h)    
