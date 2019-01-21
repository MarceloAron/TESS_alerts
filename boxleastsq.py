"""
Script for doing BLS on TESS alerts.

@author: emilk, marcelo
"""
# =============================================================================
# Basic packages 
# =============================================================================
<<<<<<< HEAD
<<<<<<< HEAD
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
=======
>>>>>>> marcelo
=======
>>>>>>> 5d66aeec13a16aa6d7550a1edc1d585e8a3f22c6
from astropy import units as u
import numpy as np

# =============================================================================
# Box least-squares of light curves 
# =============================================================================
<<<<<<< HEAD
<<<<<<< HEAD
from scipy.signal import medfilt
from bls import BLS

def boxleastsq(fits,pipe,mindur=0.5,maxdur=10.0):
=======
from bls import BLS

def boxleastsq(BJD,flux,mindur=0.5,maxdur=10.0):
>>>>>>> marcelo
=======
from bls import BLS

def boxleastsq(BJD,flux,errflux,mindur=0.5,maxdur=10.0):
>>>>>>> 5d66aeec13a16aa6d7550a1edc1d585e8a3f22c6
	'''
	Box least squares module that uses bls.py from:
		https://github.com/dfm/bls.py

	Input:
<<<<<<< HEAD
<<<<<<< HEAD
		fits    : string - the .fits file containing the light curve.
		pipe    : string - the pipeline used for data reduction.
=======
		BJD     : array  - barycentric julian dates.
		flux    : array  - the normalized flux.
>>>>>>> marcelo
=======
		BJD     : array  - barycentric julian dates.
		errflux : array  - the error of the flux.
		flux    : array  - the normalized flux.
>>>>>>> 5d66aeec13a16aa6d7550a1edc1d585e8a3f22c6
		mindur  : float  - the minimum duration of the transist.
		maxdur  : float  - the maximum duration of the transist.

	Output:
<<<<<<< HEAD
<<<<<<< HEAD
		BJD     : array  - barycentric julian dates.
		flux    : array  - the normalized flux.
		errflux : array  - the error of the flux.
=======
>>>>>>> 5d66aeec13a16aa6d7550a1edc1d585e8a3f22c6
		BLSdict : dict   - dictionary containing period, mid transit time, transit duration, transit
						   depth, and the error on the depth.
		results : object - the results from the BLS.
	'''
<<<<<<< HEAD
	## Extract TIC ID
	TIC = str(int(fits.split('-')[-3]))
	
	## Open .fits file
	file = pyfits.open(fits)
	hdr = file[0].header
	data = file[1].data

	## Extract data
	flux = data['SAP_FLUX']
	errflux = data['SAP_FLUX_ERR']
	relerr = errflux/flux
	BJD = data['TIME']

	## Remove NaNs and zeros
	idx = np.isfinite(flux) & np.isfinite(BJD) & (flux != 0)
	flux = flux[idx]
	BJD = BJD[idx]
	relerr = relerr[idx]
	errflux = errflux[idx]

	## Remove jitter from the spacecraft
	jitter = (BJD < 1346) | (1350 < BJD)
	flux = flux[jitter]
	BJD = BJD[jitter]
	relerr = relerr[jitter]
	errflux = errflux[jitter]
	
	if pipe == 'spoc':
		## Window length ~5% of the number of points
		win = round((len(BJD)*0.05)/1000)*1000 + 1 
		filt = medfilt(flux,win)
		flux = flux/filt
		## Sigma clip
		win = 101
		filt = medfilt(flux,101)
		flux_clip = flux/filt - 1.0
		MAD = np.median(np.abs(flux_clip))
		idx_out = (flux_clip > -5*MAD) & (flux_clip < 5*MAD)
		flux = flux[idx_out]
		BJD = BJD[idx_out]
		relerr = relerr[idx_out]
		errflux = errflux[idx_out]
=======
		BLSdict : dict   - dictionary containing period, mid transit time, transit duration, transit
						   depth, and the error on the depth.
		results : object - the results from the BLS.
	'''
>>>>>>> marcelo
=======
>>>>>>> 5d66aeec13a16aa6d7550a1edc1d585e8a3f22c6

	## BOX LEASTSQUARE
	durations = np.linspace(mindur, maxdur, 10)*u.hour
	model = BLS(BJD*u.day, flux)
	results = model.autopower(durations, minimum_n_transit = 2,
		                    frequency_factor=5.0)

	period = results.period[np.argmax(results.power)].value
	t0 = results.transit_time[np.argmax(results.power)].value
	dur = results.duration[np.argmax(results.power)].value
	dep = results.depth[np.argmax(results.power)]
	errdep = results.depth_err[np.argmax(results.power)]

<<<<<<< HEAD
<<<<<<< HEAD
=======
	dep_even = model.compute_stats(period,dur,t0)['depth_even'][0]
	dep_odd = model.compute_stats(period,dur,t0)['depth_odd'][0]

>>>>>>> 5d66aeec13a16aa6d7550a1edc1d585e8a3f22c6
	BLSdict = {'period' : period, 'midtransit_time' : t0,
			   'duration' : dur, 'depth' : dep, 'errdepth' : errdep,
			   'depth_even' : dep_even, 'depth_odd' : dep_odd}

<<<<<<< HEAD
	return BJD, flux, errflux, BLSdict
=======
	dep_even = model.compute_stats(period,dur,t0)['depth_even'][0]
	dep_odd = model.compute_stats(period,dur,t0)['depth_odd'][0]

	BLSdict = {'period' : period, 'midtransit_time' : t0,
			   'duration' : dur, 'depth' : dep, 'errdepth' : errdep,
			   'depth_even' : dep_even, 'depth_odd' : dep_odd}

	return BLSdict, results
>>>>>>> marcelo
=======
	return BLSdict, results
>>>>>>> 5d66aeec13a16aa6d7550a1edc1d585e8a3f22c6
