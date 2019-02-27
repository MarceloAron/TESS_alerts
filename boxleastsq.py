"""
Script for doing BLS on TESS alerts.

@author: emilk, marcelo
"""
# =============================================================================
# Basic packages 
# =============================================================================
from astropy import units as u
import numpy as np

# =============================================================================
# Box least-squares of light curves 
# =============================================================================
from bls import BLS

def boxleastsq(BJD,flux,mindur=0.5,maxdur=10.0):
	'''
	Box least squares module that uses bls.py from:
		https://github.com/dfm/bls.py

	Input:
		BJD     : array  - barycentric julian dates.
		flux    : array  - the normalized flux.
		mindur  : float  - the minimum duration of the transist.
		maxdur  : float  - the maximum duration of the transist.

	Output:
		BLSdict : dict   - dictionary containing period, mid transit time, transit duration, transit
						   depth, and the error on the depth.
		results : object - the results from the BLS.
	'''

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

	dep_even = model.compute_stats(period,dur,t0)['depth_even'][0]
	dep_odd = model.compute_stats(period,dur,t0)['depth_odd'][0]

	BLSdict = {'period' : period, 'midtransit_time' : t0,
			   'duration' : dur, 'depth' : dep, 'errdepth' : errdep,
			   'depth_even' : dep_even, 'depth_odd' : dep_odd}

	return BLSdict, results
