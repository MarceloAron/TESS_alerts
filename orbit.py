import numpy as np
import glob

# =============================================================================
<<<<<<< HEAD
<<<<<<< HEAD
# Parameters 
# =============================================================================
class OrbitalParams(object):
	'''
	The orbital parameters:
		Rp    : float       - planet-to-star ratio (in units of stellar radius).
		ecc   : float       - eccentricity of orbit.
		per   : float       - period of orbit (in days).
		omega : float       - argument of periastron (in degrees).
		T0    : float       - time of inferior conjunction (in days).
		inc   : float       - inclination (in degrees).
		a     : float       - semi-major axis (in units of stellar radius).
		LD    : str         - Limb-darkening law see StellarParams, default 'uni'.
		cs    : list        - Limb-darkening coefficients.

	Set the parameters by calling
		orbparams = OrbitalParams()
		orbparams.ecc = 0.0
	'''
	def __init__(self):
		self.Rp = 0.09
		self.ecc = 0.05
		self.per = 4.8
		self.omega = 272.
		self.T0 = 0.0
		self.a = 12.
		self.inc = 88.
		self.LD = 'uni'
		self.cs = None


=======
# Parameter class
# =============================================================================
>>>>>>> marcelo
=======
# Parameter class
# =============================================================================
>>>>>>> 5d66aeec13a16aa6d7550a1edc1d585e8a3f22c6

class StellarParams(object):
	'''
	The stellar parameters:
		Teff  : float       - effective temperature (in K).
		logg  : float       - surface gravity (in cm/s2).
		MeH   : float       - metallicity [Fe/H] (in dex).
		xi    : string      - micro-turbulende (in km/s);
								'0.0','2.0','4.0','6.0','8.0'
								Default is '2.0'.
		LD    : string      - limb darkening law used; 
								'uni'  - uniform, no LD.
							    'quad' - quadratic, default.
								'nl'   - non-linear.
								'small'- LD for a small planet. 

	Set the parameters by calling
		stelparams = StellarParams()
		stelparams.Teff = 5000.0

	Default is a sun-like star in terms of Teff, logg, and [Fe/H]. 
	'''
	def __init__(self):
		self.Teff = 5750
		self.logg = 4.5
		self.MeH = 0.0
		self.xi = '2.0'
		self.LD = 'quad'


# =============================================================================
# Limb darkening 
# =============================================================================
def magol(LD):
	'''
	Module that uses f2py to compile the IDL limb darkening routines
	by Mandel & Agol described in ADS:2002ApJ...580L.171M.
	
	Input:
		LD : string - the desired limb darkening routine, see StellarParams.
	'''
	from subprocess import call
	LDroutine = glob.glob('./f_src/occult*{}.f'.format(LD))[0]
	name = LDroutine.split('/f_src/occult')[1][:-2]
	prog = glob.glob('./{}*.so'.format(name))
	if not prog or name.endswith('nl.f'):
		call(['f2py', '-c', '-m', '{}'.format(name), '{}'.format(LDroutine), '--quiet'])
	return

def get_LDcoeff(stelpars):
	'''
	Module that collects limb darkening coefficients from Vizier.
<<<<<<< HEAD
<<<<<<< HEAD
	Calculated by A. Claret using ATLAS atmospheres for TESS.
=======
	Calculated by A. Claret using ATLAS atmospheres for TESS
>>>>>>> marcelo
=======
	Calculated by A. Claret using ATLAS atmospheres for TESS
>>>>>>> 5d66aeec13a16aa6d7550a1edc1d585e8a3f22c6
	in ADS:2017A&A...600A..30C.

	The limb darkening law is decided by the one specified in stelpars.LD.

	Input:
		stelpars : object - stellar parameters from class StellarParams
	
	Output:
		coeffs   : list   - list of LD coefficients in ascending order.
	'''
	Teff, logg, MeH = stelpars.Teff, stelpars.logg, stelpars.MeH
	xi, LD = stelpars.xi, stelpars.LD

	from astroquery.vizier import Vizier

	LDs = {'lin' : 'table24', 'quad' : 'table25', 'sqrt' : 'table26',
	       'log' : 'table27', 'nl' : 'table28', 'small' : 'table28'}
	LDs['vals'] = {'Teff' : [3500,250,50000], 'logg' : [0.0,0.5,5.0], 
	               'MeH' : [-5.0,0.5,1.0]}
<<<<<<< HEAD
<<<<<<< HEAD

=======
	
>>>>>>> marcelo
=======
	
>>>>>>> 5d66aeec13a16aa6d7550a1edc1d585e8a3f22c6
	for par in LDs['vals']:
		vals = np.arange(LDs['vals'][par][0],LDs['vals'][par][-1]+0.01,LDs['vals'][par][1])
		if par == 'Teff':
			minimum = np.argmin(np.sqrt((vals-stelpars.Teff)**2))
			Teff = vals[minimum]
<<<<<<< HEAD
<<<<<<< HEAD
=======
=======
>>>>>>> 5d66aeec13a16aa6d7550a1edc1d585e8a3f22c6
		if par == 'logg':
			minimum = np.argmin(np.sqrt((vals-stelpars.logg)**2))
			logg = vals[minimum]
		if par == 'MeH':
			minimum = np.argmin(np.sqrt((vals-stelpars.MeH)**2))
			MeH = vals[minimum]
<<<<<<< HEAD
>>>>>>> marcelo
=======
>>>>>>> 5d66aeec13a16aa6d7550a1edc1d585e8a3f22c6

	catalog = Vizier.query_constraints(catalog='J/A+A/600/A30/{}'.format(LDs[LD]),
		                               Teff='{}'.format(Teff), 
		                               logg='{}'.format(logg),
		                               Z='{}'.format(MeH), xi=xi)

	try:
		cols = catalog[0][:][0].colnames
	except IndexError:
		print('\nWARNING! No LD coefficients found for star with')
		print('Teff = {} K, logg = {} cm/s2, [Fe/H] = {}\n'.format(Teff,logg,MeH))
		stelpars = StellarParams()
		print('Using Teff = {} K, logg = {} cm/s^2, [Fe/H] = {}\n'
			.format(stelpars.Teff,stelpars.logg,stelpars.MeH))
		catalog = Vizier.query_constraints(catalog='J/A+A/600/A30/{}'.format(LDs[LD]),
			                               Teff='{}'.format(stelpars.Teff), 
		    	                           logg='{}'.format(stelpars.logg),
		        	                       Z='{}'.format(stelpars.MeH), 
		        	                       xi='{}'.format(stelpars.xi))
		cols = catalog[0][:][0].colnames
	
	coeffs = []
	for name in cols:
		if name.endswith('LSM'):
			coeff = catalog[0][:][0][name]
			coeffs.append(coeff)

	return coeffs

# =============================================================================
# Keplerian motion 
# =============================================================================
def solve_keplers_eq(mean_anomaly, ecc, tolerance=1.e-5):
	'''
	Module that solves Kepler's equation:
	M = E - sin(E) ,
	where M is the mean anomaly and E the eccentric anomaly.

	This is done following the Newton-Raphson method as described by
	Carl D. Murray and Alexandre C. M. Correia in arXiv:1009.1738v2

	Input:
		mean_anomaly    : float array - the mean anomaly.
	
	Output:
		new_ecc_anomaly : float array - the new eccentric anomaly.
	'''
	## Circular orbit
	if ecc == 0: return mean_anomaly 

	new_ecc_anomaly = mean_anomaly
	converged = False

	for ii in range(100):
		old_ecc_anomaly = new_ecc_anomaly

		new_ecc_anomaly = old_ecc_anomaly - (old_ecc_anomaly - ecc*np.sin(old_ecc_anomaly) - mean_anomaly)/(1.0 - ecc*np.cos(old_ecc_anomaly))

		if np.max(np.abs(new_ecc_anomaly - old_ecc_anomaly)/old_ecc_anomaly) < tolerance:
			converged = True
			break

	if not converged:
		print('Calculation of the eccentric anomaly did not converge!')

	return new_ecc_anomaly

def true_dist(time, T0, per, a, inc, ecc, omega):
	'''
	Module that returns the projected sky distance calculated 
	from the true anomaly.
	Again following:
	Carl D. Murray and Alexandre C. M. Correia in arXiv:1009.1738v2.
	
	Input:
		time   : float array - times of observations.
		orbpars: floats      - see orbital parameters from class OrbitalParams.

	Output:
		d      : float array - the projected distance.
	'''
	n = 2.0*np.pi/per
	mean_anomaly = n*(time-T0)

	ecc_anomaly = solve_keplers_eq(mean_anomaly,ecc)

	cos_E = np.cos(ecc_anomaly)
	sin_E = np.sin(ecc_anomaly)

	## Cosine and sine of the true anomaly
	cos_f = (cos_E - ecc)/(1.0 - ecc*cos_E)
	sin_f = (np.sqrt(1 - ecc**2)*sin_E)/(1.0 - ecc*cos_E)
	
	## Convert angle from degrees to radians
	omega = np.pi/180.*omega
	inc = np.pi/180.*inc
	nn = len(cos_f)
	d = np.zeros(nn)
	
	## NOTE: Expressions like sin(omega + f) are expanded to stay clear of arctan
	for ii in range(nn):
		## Huge value for separation to make sure not to model planet passing behind star
<<<<<<< HEAD
<<<<<<< HEAD
		if np.sin(inc)*(np.sin(omega)*cos_f[ii] + np.cos(omega)*sin_f[ii]) > 0:
=======
		if np.sin(inc)*(np.sin(omega)*cos_f[ii] + np.cos(omega)*sin_f[ii]) < 0:
>>>>>>> marcelo
=======
		if np.sin(inc)*(np.sin(omega)*cos_f[ii] + np.cos(omega)*sin_f[ii]) < 0:
>>>>>>> 5d66aeec13a16aa6d7550a1edc1d585e8a3f22c6
			d[ii] = 1000.
		else:
			nom = a*(1.0 - ecc**2)
			nom *= np.sqrt(1.0 - (np.sin(omega)*cos_f[ii] + np.cos(omega)*sin_f[ii])**2*np.sin(inc)**2)
			den = 1.0 + ecc*cos_f[ii]
			d[ii] = nom/den
	return d

# =============================================================================
# Light curve 
# =============================================================================
def LCuni(d, p):
	'''
	Module that returns a uniform light curve modeled as
	described in ADS:2002ApJ...580L.171M by Mandel and Agol.
	Input:
		d        : float array - the projected distance.
		p        : float       - the planet-to-star radius ratio.

	Output:
		flux_uni : float array - the flux from a uniform source.
<<<<<<< HEAD
<<<<<<< HEAD
		flux     : float array - dummy array containing ones. 
=======
>>>>>>> marcelo
=======
>>>>>>> 5d66aeec13a16aa6d7550a1edc1d585e8a3f22c6
	'''
	nz = len(d)
	flux = np.ones(nz)
	lam = np.zeros(nz)
	for i, z in enumerate(d):
		if 1 + p < z:
			lam[i] = 0
		elif abs(1-p) < z and z <= 1 + p:
			k0 = np.arccos((p**2 + z**2 - 1)/(2*p*z))
			k1 = np.arccos((1 - p**2 + z**2)/(2*z))
			lam[i] = (p**2*k0 + k1 - np.sqrt((4*z**2 - (1 + z**2 - p**2)**2)/4))/np.pi
		elif z <= 1 - p:
			lam[i] = p**2
		elif z <= p - 1:
			lam[i] = 1
	
	flux_uni = flux - lam
<<<<<<< HEAD
<<<<<<< HEAD
	return flux_uni, flux
=======
	return flux_uni
>>>>>>> 5d66aeec13a16aa6d7550a1edc1d585e8a3f22c6


<<<<<<< HEAD
class TransitModel(object):
	def __init__(self, orbpars, time, sample_factor=1, exp_time=0.):
		'''
		Input:
			orbpars  : object      - orbital parameters from class OrbitalParams.
			time     : float array - times of observations.
		'''
		self.time = time

		# Super sampling
		self.sample_factor = sample_factor
		self.exp_time = exp_time
		if self.sample_factor > 1:
			ts = np.linspace(-self.exp_time/2., self.exp_time/2., self.sample_factor)
			self.t_sampled = (ts + self.time.reshape(self.time.size, 1)).flatten()
		else: self.t_sampled = self.time
		
		## Orbital parameters
		self.T0 = orbpars.T0
		self.ecc = orbpars.ecc
		self.per = orbpars.per
		self.omega = orbpars.omega
		self.a = orbpars.a
		self.inc = orbpars.inc

		## Get the projected seperation
		self.d = true_dist(self.t_sampled, self.T0, self.per, self.a, self.inc, self.ecc, self.omega)

	def get_LC(self, orbpars):
		'''
		Module that returns the light curve given of the orbit modeled
		as described in ADS:2002ApJ...580L.171M by Mandel and Agol.
		
		Input:
			orbpars : object      - orbital parameters from class OrbitalParams.

		Output:
			lc      : float array - light curve.

		'''
		if orbpars.per != self.per or orbpars.T0 != self.T0 or orbpars.ecc != self.ecc or orbpars.omega != self.omega or orbpars.a != self.a or orbpars.inc != self.inc:
			self.d = true_dist(self.t_sampled, self.T0, self.per, self.a, self.inc, self.ecc, self.omega)
		
		self.T0 = orbpars.T0
		self.ecc = orbpars.ecc
		self.per = orbpars.per
		self.omega = orbpars.omega
		self.a = orbpars.a
		self.inc = orbpars.inc
		self.p = orbpars.Rp
		self.LD = orbpars.LD
		self.coeffs = orbpars.cs		

		flux_uni, flux = LCuni(self.d, self.p)

		if self.LD == 'uni':
			lc = flux_uni
		elif self.LD == 'quad':
			magol(self.LD)
			import quad
			nz = len(flux)
			quad.occultquad(self.d,self.coeffs[0],self.coeffs[1],self.p,flux,flux_uni,nz)
			lc = flux

		if self.sample_factor == 1: return lc
		else: return np.mean(lc.reshape(-1, self.sample_factor), axis=1)
=======
	return flux_uni


=======
>>>>>>> 5d66aeec13a16aa6d7550a1edc1d585e8a3f22c6
def LCquad(d, c1, c2, p, flux, flux_uni):
	'''
	Module that returns the light curve given of the orbit modeled
	as described in ADS:2002ApJ...580L.171M by Mandel and Agol.
	
	Input:
		d        : float array  - projected sky distance.
		c1       : float        - linear ld coefficient.
		c2       : float        - quadratic ld coefficient.
		flux     : float array  - array to be stored.
		flux_uni : float array  - light curve from uniform source.

	Output:
		lc       : float array - light curve.

	'''
	magol('quad')
	import quad
	nz = len(flux)
	quad.occultquad(d,c1,c2,p,flux,flux_uni,nz)
	lc = flux
	return lc
<<<<<<< HEAD
>>>>>>> marcelo
=======
>>>>>>> 5d66aeec13a16aa6d7550a1edc1d585e8a3f22c6
