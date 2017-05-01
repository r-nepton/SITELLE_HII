import numpy as np
import pylab as pl

import pyfits
import orb.fit
import orb.utils.spectrum
import orb.utils.io
from orb.core import Lines

# PARAMETERS
order = 8 # In cube header
step = 2880.879435# In cube header
axis_corr = 1.03664546289 # In cube header
shift_guess = 700 # mean velocity guesss in km/s
velocity_range = 200
lines_nm = Lines().get_line_nm(['[NII]6548','Halpha','[NII]6583', 'HeI6678', '[SII]6716', '[SII]6731']) # lines to fit
print lines_nm
cov_pos = [1,1,1,1,1,1] # groups of velocity, [1,2,1] means that the first and third lines share the same velocity
cov_fwhm = [1,1,2,2,3,3]
signal_range = [650,676] # range of the signal to fit (in nm)
fwhm_guess = [0.36,0.36,0.38,0.38,0.40,0.40] # FWHM guess of the lines in nm
	
# computed params
nm_laser = 543.5
#hene = pyfits.getdata('LASER_None.cam1.calibration_laser_map.fits')
nm_laser_obs = nm_laser*axis_corr

# Read HII Regions files
mask = pyfits.getdata('NGC628_HIIPhot3_All_GROW_SN3.fits')
dist = pyfits.getdata('NGC628_distance_grow.fits') 
sig = np.loadtxt('NGC628_sigma.txt')
# sig = sig//14.1443 # to get the pixel values instead of the pc scale
cube = pyfits.getdata('NGC628_SN3.merged.nm.1.0_conv3_ciel_pop_tout.fits')
hdr = pyfits.open('NGC628_SN3.merged.nm.1.0_conv3_ciel_pop_tout.fits')
h = hdr[0].header

dump_file = np.zeros([2064,2048])
spectrum_nm = np.zeros([323,4285])
sp = np.zeros([323])
spectrum_fit_cube = np.zeros([323,4285])
spectrum_fit = np.zeros([323])
vel_ha = np.zeros([4285])
vel_ha_err = np.zeros([4285])
vel_nii = np.zeros([4285])
vel_nii_err = np.zeros([4285])
vel_sii = np.zeros([4285])
vel_sii_err = np.zeros([4285])
amp_ha = np.zeros([4285])
amp_nii1 = np.zeros([4285])
amp_nii2 = np.zeros([4285])
amp_hei = np.zeros([4285])
amp_sii1 = np.zeros([4285])
amp_sii2 = np.zeros([4285])
fwhm_ha = np.zeros([4285])
fwhm_nii1 = np.zeros([4285])
fwhm_nii2 = np.zeros([4285])
fwhm_hei = np.zeros([4285])
fwhm_sii1 = np.zeros([4285])
fwhm_sii2 = np.zeros([4285])
heig_ha = np.zeros([4285])
heig_nii1 = np.zeros([4285])
heig_nii2 = np.zeros([4285])
heig_hei = np.zeros([4285])
heig_sii1 = np.zeros([4285])
heig_sii2 = np.zeros([4285])
snr_ha = np.zeros([4285])
snr_nii1 = np.zeros([4285])
snr_nii2 = np.zeros([4285])
snr_hei = np.zeros([4285])
snr_sii1 = np.zeros([4285])
snr_sii2 = np.zeros([4285])
amp_ha_err = np.zeros([4285])
amp_nii1_err = np.zeros([4285])
amp_nii2_err = np.zeros([4285])
amp_hei_err = np.zeros([4285])
amp_sii1_err = np.zeros([4285])
amp_sii2_err = np.zeros([4285])
fwhm_ha_err = np.zeros([4285])
fwhm_nii1_err = np.zeros([4285])
fwhm_nii2_err = np.zeros([4285])
fwhm_hei_err = np.zeros([4285])
fwhm_sii1_err = np.zeros([4285])
fwhm_sii2_err = np.zeros([4285])
heig_ha_err = np.zeros([4285])
heig_nii1_err = np.zeros([4285])
heig_nii2_err = np.zeros([4285])
heig_hei_err = np.zeros([4285])
heig_sii1_err = np.zeros([4285])
heig_sii2_err = np.zeros([4285])

#indicex = []
#indicey = []
#
#nm_laser_obs = np.zeros([4285])
#for i in range(1,4286):
#    indicex.append(np.where((mask == i) & (dist < sig[i-1]))[0])
#    indicey.append(np.where((mask == i) & (dist < sig[i-1]))[1])
#    nm_laser_obs[i-1] = np.mean(hene[indicex[i-1],indicey[i-1]])
#
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

#for j in range(0,323):
#    dump_file[0:2064,0:2048] = cube[j,0:2064,0:2048] 	
#    
#    for i in range(0,4285):
#	print i
#        spectrum_nm[j,i] = np.sum(dump_file[indicex[i],indicey[i]])
    
#pyfits.writeto('NGC628_SN3.merged.nm.1.0_conv3_ciel_pop_tout_regions.fits', spectrum_nm, h)    
spectrum_nm = pyfits.getdata('NGC628_SN3.merged.nm.1.0_conv3_ciel_pop_tout_regions.fits')

for i in range(0,4285):
    sp[0:323] = spectrum_nm[0:323,i]
    # fit loop
    
    fit_results = orb.fit.fit_lines_in_spectrum(
        sp, lines_nm, step, order, nm_laser,
        nm_laser_obs=nm_laser_obs, wavenumber=False,
        fwhm_guess=fwhm_guess, cont_guess=None,
        shift_guess=shift_guess, sigma_guess=0.,
        fix_fwhm=False, cov_fwhm=cov_fwhm, cov_pos=cov_pos,
        fix_pos=False, cov_sigma=False,
        fit_tol=1e-10, poly_order=0,
        fmodel='sinc', signal_range=[650,676],
        filter_file_path=None, fix_filter=False,
        apodization=1., velocity_range=velocity_range,
        compute_mcmc_error=False, no_error=False)
        
    spectrum_fit = fit_results['fitted-vector']
    
    for j in range(0,323):		
        spectrum_fit_cube[j,i] = float(spectrum_fit[j])


    print i
    if fit_results != []: # if the fit is ok 
    	try :
            vel_ha[i] = fit_results['velocity'][1]
            vel_nii[i] = fit_results['velocity'][2]
            vel_sii[i] = fit_results['velocity'][3]
            amp_ha[i] = fit_results['lines-params'][1,1]
            amp_nii1[i] = fit_results['lines-params'][0,1]
            amp_nii2[i] = fit_results['lines-params'][2,1]
            amp_hei[i] = fit_results['lines-params'][3,1]
            amp_sii1[i] = fit_results['lines-params'][4,1]
            amp_sii2[i] = fit_results['lines-params'][5,1]
            fwhm_ha[i] = fit_results['lines-params'][1,3]
            fwhm_nii1[i] = fit_results['lines-params'][0,3]
            fwhm_nii2[i] = fit_results['lines-params'][2,3]
            fwhm_hei[i] = fit_results['lines-params'][3,3]
            fwhm_sii1[i] = fit_results['lines-params'][4,3]
            fwhm_sii2[i] = fit_results['lines-params'][5,3]
            heig_ha[i] = fit_results['lines-params-err'][1,0]
            heig_nii1[i] = fit_results['lines-params-err'][0,0]
            heig_nii2[i] = fit_results['lines-params-err'][2,0]
            heig_hei[i] = fit_results['lines-params-err'][3,0]
            heig_sii1[i] = fit_results['lines-params-err'][4,0]
            heig_sii2[i] = fit_results['lines-params-err'][5,0]
            snr_ha[i] = fit_results['snr'][1,0]
            snr_nii1[i] = fit_results['snr'][0,0]
            snr_nii2[i] = fit_results['snr'][2,0]
            snr_hei[i] = fit_results['snr'][3,0]
            snr_sii1[i] = fit_results['snr'][4,0]
            snr_sii2[i] = fit_results['snr'][5,0]

        except :
	    pass
            
    		
	try :

            vel_ha_err[i] = fit_results['velocity-err'][1]
            vel_nii_err[i] = fit_results['velocity-err'][2]
            vel_sii_err[i] = fit_results['velocity-err'][3]
            amp_ha_err[i] = fit_results['lines-params-err'][1,1]
            amp_nii1_err[i] = fit_results['lines-params-err'][0,1]
            amp_nii2_err[i] = fit_results['lines-params-err'][2,1]
            amp_hei_err[i] = fit_results['lines-params-err'][3,1]
            amp_sii1_err[i] = fit_results['lines-params-err'][4,1]
            amp_sii2_err[i] = fit_results['lines-params-err'][5,1]
            fwhm_ha_err[i] = fit_results['lines-params-err'][1,3]
            fwhm_nii1_err[i] = fit_results['lines-params-err'][0,3]
            fwhm_nii2_err[i] = fit_results['lines-params-err'][2,3]
            fwhm_hei_err[i] = fit_results['lines-params-err'][3,3]
            fwhm_sii1_err[i] = fit_results['lines-params-err'][4,3]
            fwhm_sii2_err[i] = fit_results['lines-params-err'][5,3]
            heig_ha_err[i] = fit_results['lines-params-err'][1,0]
            heig_nii1_err[i] = fit_results['lines-params-err'][0,0]
            heig_nii2_err[i] = fit_results['lines-params-err'][2,0]
            heig_hei_err[i] = fit_results['lines-params-err'][3,0]
            heig_sii1_err[i] = fit_results['lines-params-err'][4,0]
            heig_sii2_err[i] = fit_results['lines-params-err'][5,0]
        
   	except :
	    pass
        
# Write results
    
pyfits.writeto('NGC628_SN3.regions.Ha6563_velocity.fits', vel_ha, h)    
pyfits.writeto('NGC628_SN3.regions.NII6583_velocity.fits', vel_nii, h)    
pyfits.writeto('NGC628_SN3.regions.SII6716_velocity.fits', vel_sii, h)    
pyfits.writeto('NGC628_SN3.regions.Ha6563_velocity_error.fits', vel_ha_err, h)    
pyfits.writeto('NGC628_SN3.regions.NII6583_velocity_error.fits', vel_nii_err, h)
pyfits.writeto('NGC628_SN3.regions.SII6716_velocity_error.fits', vel_sii_err, h)
pyfits.writeto('NGC628_SN3.regions.fit.spectrum.fits', spectrum_fit_cube, h)
    
pyfits.writeto('NGC628_SN3.regions.Ha6563_amplitude.fits', amp_ha, h)    
pyfits.writeto('NGC628_SN3.regions.NII6548_amplitude.fits', amp_nii1, h)    
pyfits.writeto('NGC628_SN3.regions.NII6583_amplitude.fits', amp_nii2, h)    
pyfits.writeto('NGC628_SN3.regions.HeI6678_amplitude.fits', amp_hei, h)    
pyfits.writeto('NGC628_SN3.regions.SII6716_amplitude.fits', amp_sii1, h)    
pyfits.writeto('NGC628_SN3.regions.SII6731_amplitude.fits', amp_sii2, h)    
pyfits.writeto('NGC628_SN3.regions.Ha6563_fwhm.fits', fwhm_ha, h)    
pyfits.writeto('NGC628_SN3.regions.NII6548_fwhm.fits', fwhm_nii1, h)    
pyfits.writeto('NGC628_SN3.regions.NII6583_fwhm.fits', fwhm_nii2, h)    
pyfits.writeto('NGC628_SN3.regions.HeI6678_fwhm.fits', fwhm_hei, h)    
pyfits.writeto('NGC628_SN3.regions.SII6716_fwhm.fits', fwhm_sii1, h)    
pyfits.writeto('NGC628_SN3.regions.SII6731_fwhm.fits', fwhm_sii2, h)  
pyfits.writeto('NGC628_SN3.regions.Ha6563_heigt.fits', heig_ha, h)    
pyfits.writeto('NGC628_SN3.regions.NII6548_heigt.fits', heig_nii1, h)    
pyfits.writeto('NGC628_SN3.regions.NII6583_heigt.fits', heig_nii2, h)    
pyfits.writeto('NGC628_SN3.regions.HeI6678_heigt.fits', heig_hei, h)    
pyfits.writeto('NGC628_SN3.regions.SII6716_heigt.fits', heig_sii1, h)    
pyfits.writeto('NGC628_SN3.regions.SII6731_heigt.fits', heig_sii2, h)  
pyfits.writeto('NGC628_SN3.regions.Ha6563_snr.fits', snr_ha, h)    
pyfits.writeto('NGC628_SN3.regions.NII6548_snr.fits', snr_nii1, h)    
pyfits.writeto('NGC628_SN3.regions.NII6583_snr.fits', snr_nii2, h)    
pyfits.writeto('NGC628_SN3.regions.HeI6678_snr.fits', snr_hei, h)    
pyfits.writeto('NGC628_SN3.regions.SII6716_snr.fits', snr_sii1, h)    
pyfits.writeto('NGC628_SN3.regions.SII6731_snr.fits', snr_sii2, h)  
    
pyfits.writeto('NGC628_SN3.regions.Ha6563_amplitude_error.fits', amp_ha_err, h)    
pyfits.writeto('NGC628_SN3.regions.NII6548_amplitude_error.fits', amp_nii1_err, h)    
pyfits.writeto('NGC628_SN3.regions.NII6583_amplitude_error.fits', amp_nii2_err, h)    
pyfits.writeto('NGC628_SN3.regions.HeI6678_amplitude_error.fits', amp_hei_err, h)    
pyfits.writeto('NGC628_SN3.regions.SII6716_amplitude_error.fits', amp_sii1_err, h)    
pyfits.writeto('NGC628_SN3.regions.SII6731_amplitude_error.fits', amp_sii2_err, h)    
pyfits.writeto('NGC628_SN3.regions.Ha6563_fwhm_error.fits', fwhm_ha_err, h)    
pyfits.writeto('NGC628_SN3.regions.NII6548_fwhm_error.fits', fwhm_nii1_err, h)    
pyfits.writeto('NGC628_SN3.regions.NII6583_fwhm_error.fits', fwhm_nii2_err, h)    
pyfits.writeto('NGC628_SN3.regions.HeI6678_fwhm_error.fits', fwhm_hei_err, h)    
pyfits.writeto('NGC628_SN3.regions.SII6716_fwhm_error.fits', fwhm_sii1_err, h)    
pyfits.writeto('NGC628_SN3.regions.SII6731_fwhm_error.fits', fwhm_sii2_err, h)  
pyfits.writeto('NGC628_SN3.regions.Ha6563_heigt_error.fits', heig_ha_err, h)    
pyfits.writeto('NGC628_SN3.regions.NII6548_heigt_error.fits', heig_nii1_err, h)    
pyfits.writeto('NGC628_SN3.regions.NII6583_heigt_error.fits', heig_nii2_err, h)    
pyfits.writeto('NGC628_SN3.regions.HeI6678_heigt_error.fits', heig_hei_err, h)    
pyfits.writeto('NGC628_SN3.regions.SII6716_heigt_error.fits', heig_sii1_err, h)    
pyfits.writeto('NGC628_SN3.regions.SII6731_heigt_error.fits', heig_sii2_err, h)      
    
    
    