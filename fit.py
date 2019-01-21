"""
Script for fitting light curves.

@author: emilk, marcelo
"""
# =============================================================================
# Basic packages 
# =============================================================================
import numpy as np 
import lmfit
import astropy.io.fits as pyfits
from scipy.signal import medfilt

# =============================================================================
# Light curve generator
# =============================================================================
import orbit

def lc_model_fit(parameters, bjd, flux, errflux, LD):
	'''
	Module that returns the chi-squared statistics.
	Input:
		bjd        : float array - dates.
		flux       : float array - flux.
		errflux    : float array - the uncertainty for the flux.
		parameters : object - orbital parameters.
		LD         : string - limb darkening law, either "uni" or "quad".

	Output:
		res_lc     : float array - the chi-squared statistics.
	'''
	### Unpack vallues from parameters
	P = parameters['P'].value
	T0 = parameters['T0'].value
	a = parameters['a'].value
	Rp = parameters['Rp'].value
	inc = parameters['inc'].value
	ecc = parameters['ecc'].value
	omega = parameters['omega'].value

	### Get model
	d = orbit.true_dist(bjd, T0, P, a, inc, ecc, omega)
	lc = orbit.LCuni(d, Rp)

	### Check if LD is included
	if LD == 'quad':
		c1 = parameters['c1'].value
		c2 = parameters['c2'].value
		### Get LD'ed model
		dummy = np.zeros(len(lc))
		model = orbit.LCquad(d, c1, c2, Rp, dummy, lc)
	else:
		model = lc
	
	res_lc = (flux - model)**2/errflux**2
	return res_lc

def fit_LC(bjd, flux, errflux, P=1., T0=0., a=8., Rp=0.1, inc=90., ecc=0., omega=90., LD=None, c1=0.3, c2=0.3):
	'''
	Module that utilizes the Levenberg-Marquardt least-squares routine to
	optimize the light curve.
	Input:
		bjd                           : float array - dates.
		flux                          : float array - flux.
		errflux                       : float array - the uncertainty for the flux.
		P, T0, a, Rp, inc, ecc, omega : floats - orbital parameters.
		LD                            : string - limb darkening law, either "uni" or "quad".
		c1, c2                        : floats - limb darkening coefficients.

	Output:
		model_fit                     : object - the fitted parameters.
	'''
	### Parameters class from lmfit
	parameters = lmfit.Parameters()

	### Orbital parameters
	parameters.add('P', value=P, min=0., vary=True)
<<<<<<< HEAD
	parameters.add('T0', value=T0, min=T0-0.1, max=T0+0.1, vary=True)
	parameters.add('a', value=a, min=0., vary=True)
	parameters.add('Rp', value=Rp, min=0., max=1., vary=True)
	parameters.add('inc', value=inc, min=80., max=90., vary=True)
	parameters.add('ecc', value=ecc, min=0., max=1., vary=False)
	parameters.add('omega', value=omega, min=0., max=360., vary=False)
=======
	parameters.add('T0', value=T0, vary=True)
	parameters.add('a', value=a, min=0., vary=True)
	parameters.add('Rp', value=Rp, min=0., max=1., vary=True)
	parameters.add('inc', value=inc, min=0., max=90., vary=True)
	parameters.add('ecc', value=ecc, min=0., max=1., vary=True)
	parameters.add('omega', value=omega, min=0., max=360., vary=True)
>>>>>>> 5d66aeec13a16aa6d7550a1edc1d585e8a3f22c6

	### Limb darkening coefficients
	if LD not in ['uni','quad']:
		print('Limb darkening law must be either "uni" or "quad".\n')
		print('Setting LD="uni".\n')
		LD = 'uni'
	elif LD == 'quad':
		parameters.add('c1', value=c1, min=c1-0.1, max=c1+0.1, vary=True)
		parameters.add('c2', value=c2, min=c2-0.1, max=c2+0.1, vary=True)

	model_fit = lmfit.minimize(lc_model_fit, parameters, args=(bjd, flux, errflux, LD), ftol=1e-15, maxfev=3000)
	return model_fit
