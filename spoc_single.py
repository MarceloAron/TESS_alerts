"""
Script for fitting light curves and producing a pdf file with the result.

@author: emilk, marcelo
"""
# =============================================================================
# Basic packages 
# =============================================================================
import numpy as np 
import lmfit
import glob
import astropy.io.fits as pyfits
from scipy.signal import medfilt
import medianfilter

# =============================================================================
# Plotting
# =============================================================================
import matplotlib.pyplot as plt
import seaborn as sns
cols = sns.color_palette()
from matplotlib.backends.backend_pdf import PdfPages

# =============================================================================
# Scripts 
# =============================================================================
import orbit
import fit
import boxleastsq

# -----------------------------------------------------------------------------

pdf = PdfPages('290131778_2.pdf')

### Open .fits file
TIC = '290131778'
fname = glob.glob('fits/*{}*s_llc.fits'.format(TIC))[0]
file = pyfits.open(fname)
hdr = file[0].header
data = file[1].data

print('\n')
### Extract information
### Physical parameters
try: 
	Teff = float(hdr['TEFF'])
	Teff_L = int(Teff)
except TypeError:
	print('No effective temperature available for TIC {}.'.format(TIC))
	print('Using 5750 K to obtain LD coefficients.')
	Teff = 5750
	Teff_L = '\mathrm{N/A}'
try: 
	logg = float(hdr['logg'])
	logg_L = '{:0.2f}'.format(logg)
except TypeError:
	print('No surface gravity available for TIC {}.'.format(TIC))
	print('Using log g = 4.5 dex to obtain LD coefficients.')
	logg = 4.5
	logg_L = '\mathrm{N/A}'
try: 
	MeH = float(hdr['MH'])
	MeH_L ='{:0.2f}'.format(MeH)
except TypeError:
	print('No metallicity available for TIC {}.'.format(TIC))
	print('Using [M/H] = 0.0 dex to obtain LD coefficients.')
	MeH = 0.0	
	MeH_L = '\mathrm{N/A}'
try: 
	rstar = float(hdr['RADIUS'])
	rstar_L ='{:0.2f}'.format(rstar)
except TypeError:
	rstar_L = '\mathrm{N/A}'
print('\n')

### Observational parameters
Tmag = float(hdr['TESSMAG'])
RA = float(hdr['RA_OBJ'])
DEC = float(hdr['DEC_OBJ'])
SEC = int(hdr['SECTOR'])
CAM = hdr['CAMERA']

### Extract data
bjd = data['TIME']
flux = data['PDCSAP_FLUX']
errflux = data['PDCSAP_FLUX_ERR']

### Remove NaNs and zeros
idx = np.isfinite(flux) & np.isfinite(bjd) & (flux != 0)
flux = flux[idx]
bjd = bjd[idx]
errflux = errflux[idx]
relerr = errflux/flux

## Separate data into chunks
start = np.amin(bjd)
end = np.amax(bjd)
down = 13.94
chunks = int(round((end-start)/down))

## Raw flux and filter for pĺots
raw_flux = np.ndarray(shape=len(flux))
raw_flux[:] = flux
time = np.ndarray(shape=len(bjd))
time[:] = bjd
filts = np.ndarray(shape=len(flux))

## Looping over all chunks and applying median filter
## with window size of 5% the number of points
for cc in range(chunks):
	idxs = np.where((bjd >= (start + cc*down)) & (bjd <= (start + (cc+1)*down)))[0]
	win = round(len(bjd[idxs])*0.05)
	if win%2 == 0: win += 1
	filt = medianfilter.my_median_filter(flux[idxs],win)
	filts[idxs] = filt
	flux[idxs] = flux[idxs]/filt
errflux = relerr*flux

## Second median filter, now with smaller window
win = round(len(bjd)*0.001)
if win%2 == 0: win += 1
filt = medianfilter.my_median_filter(flux,win)

## Sigma clip
flux_clip = flux/filt - 1.0
MAD = np.median(np.abs(flux_clip))
idx_out = (flux_clip > -5*MAD) & (flux_clip < 5*MAD)
flux = flux[idx_out]
bjd = bjd[idx_out]
relerr = relerr[idx_out]
errflux = errflux[idx_out]

## Applying the BLS method
BLS,results = boxleastsq.boxleastsq(bjd,flux)

## Flag the candidate based on the periodogram
flag = 'Good'
rr = results.period.value > 20
rr_power = results.power[rr]
if np.amax(results.power)/np.amax(rr_power) < 1.1: flag = 'Bad'

### Good estimates of the orbital parameters are 
### obtained from ISBN 978-0-521-76559 by M. Perryman
fold = BLS['midtransit_time'] # Folded midtransit time
per = BLS['period'] # Period in days
tT = BLS['duration'] # Duration in days
depth = BLS['depth'] # Depth
dep_even = BLS['depth_even']
dep_odd = BLS['depth_odd']

dep_rat = dep_even/dep_odd
binary = 'Likely'
if 0.75 < dep_rat < 1.25:
	binary = 'Unlikely'

mid = np.argmin(np.sqrt((bjd%per - fold)**2))
T0_bls = float(bjd[mid]) # actual midtransit time

a = per/(np.pi*tT) # eq. (6.2), with i ~ 90 deg, Rp << Rstar, tT << P
Rp = np.sqrt(depth) # eq. (6.13)

### Set stellar parameters to get LD coefficients
stelpars = orbit.StellarParams()
stelpars.Teff = Teff
stelpars.logg = logg
stelpars.MeH = MeH
limb = stelpars.LD
c1, c2 = orbit.get_LDcoeff(stelpars)

### Get best fitting model
lc_fit = fit.fit_LC(bjd, flux, errflux, P=per, T0=T0_bls, a=a, Rp=Rp, inc=89.9, ecc=0., omega=90., LD=limb, c1=c1, c2=c2)
print(lmfit.fit_report(lc_fit, show_correl=False))
per = lc_fit.params['P'].value
T0 = lc_fit.params['T0'].value
a = lc_fit.params['a'].value
Rp = lc_fit.params['Rp'].value
inc = lc_fit.params['inc'].value
ecc = lc_fit.params['ecc'].value
omega = lc_fit.params['omega'].value


### Generate light curve for best fitting model
sep = orbit.true_dist(bjd, T0, per, a, inc, ecc, omega)
dummy = np.zeros(len(flux))
lc_uni = orbit.LCuni(sep, Rp)
if limb == 'quad':
	c1, c2 = lc_fit.params['c1'].value, lc_fit.params['c2'].value
	lc = orbit.LCquad(sep, c1, c2, Rp, dummy, lc_uni)
else:
	lc = lc_uni
idx = np.argmin(lc)
phase = (bjd - T0_bls + 0.5*per)%per - 0.5*per
shift = np.argsort(phase)

### Fit statistics
chi2 = lc_fit.redchi

### Plotting
plt.rc('text',usetex=True)
plt.rc('text.latex',preamble=r'\usepackage{color}')
fig = plt.figure()
fig.suptitle(r'$\rm TIC \ {}$'.format(TIC))

ax = fig.add_subplot(1,1,1)
ax.set_title(r'$\rm Sector \ {:d}, \ Camera \ {}$'.format(SEC,CAM))
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['right'].set_color('none')
ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')

axt = fig.add_subplot(3,2,2)
axt.spines['top'].set_color('none')
axt.spines['bottom'].set_color('none')
axt.spines['left'].set_color('none')
axt.spines['right'].set_color('none')
axt.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')	

first = [r'$\rm {\bf Stellar}$',
r'$\rm Tmag={:.2f}$'.format(Tmag),
r'$\rm RA={:.2f}^\circ$'.format(RA),
r'$\rm DEC={:.2f}^\circ$'.format(DEC),
r'$\ $',
r'$T_{} = {} \ \rm K$'.format('\mathrm{eff}',Teff_L),
r'$\log g = {} \ \rm dex$'.format(logg_L), 
r'$\rm [M/H] = {} \ \rm dex$'.format(MeH_L),
r'$\rm R_\star = {} \ R_\odot$'.format(rstar_L)]

font = 6
for ff, val in enumerate(first): axt.text(-6,ff,val,fontsize=font)

second = [r'$\rm {\bf BLS}$',
r'$P = {:.4f} \ \rm days$'.format(BLS['period']), 
r'$T_0 = {:.4f} \ \rm BTJD$'.format(T0_bls),
r'$t_{}={:0.4f} \ \rm days$'.format('\mathrm{T}',tT),
r'$\delta={:0.6f}$'.format(depth),
r'$\delta_{}={:0.6f}$'.format('\mathrm{odd}',dep_odd),
r'$\delta_{}={:0.6f}$'.format('\mathrm{even}',dep_even),
r'$\rm Flag: \ {}$'.format(flag),
r'$\rm Binary: \ {}$'.format(binary)]

for ss, val in enumerate(second): 
	col = 'k'
	if val.endswith('Likely$') or val.endswith('Bad$'): col = cols[3]
	axt.text(-2,ss,val,fontsize=font,color=col)

third = [r'$\rm {\bf LMfit}$',
r'$\chi^2_\nu = {:.2f}$'.format(chi2),
r'$P = {:.4f} \ \rm days$'.format(per), 
r'$T_0 = {:.4f} \ \rm BTJD$'.format(T0),
r'$R_{} = {:.3f} \ R_\star$'.format('\mathrm{p}',Rp),
r'$a = {:.3f} \ R_\star$'.format(a),
r'$i = {:.2f}^\circ$'.format(inc)]

for tt, val in enumerate(third): axt.text(2.5,tt,val,fontsize=font)

axt.set_xlim(-5,5)
axt.set_ylim(len(first)-2,-1)

### Plot light curve
ax1 = fig.add_subplot(3,2,1)
ax1.plot(phase*24,flux,'.',markersize=1.0,color=cols[7])
ax1.plot(phase[shift]*24,lc[shift],'-',color=cols[0],linewidth=1)
dur_h = (tT*24)/2
ax1.set_xlim(-1*(dur_h+2),dur_h+2)
ax1.set_xlabel(r'$\rm Hours \ from \ Midtransit$')
ax1.set_ylabel(r'$\rm Relative \ Brightness$')

### Plot periodogram from BLS
ax2 = fig.add_subplot(3,1,2)
ax2.axvline(per,color=cols[0],lw=1)
ax2.plot(results.period,results.power,'k',lw=0.5)
ax2.set_xlabel(r'$\rm Period \ (days)$')
ax2.set_ylabel(r'$\rm \log \ likelihood$')

### Plot raw flux
ax3 = fig.add_subplot(3,1,3)
ax3.plot(time,raw_flux,'.',markersize=1.0,color=cols[7])
ax3.plot(time,filts,'.',ms=1.,color='black')
for n in range(-30,30):
	ax3.axvline(T0+n*per, ls='dashed', color=cols[0], lw=0.5)
ax3.set_xlim(time.min()-1,time.max()+1)
ax3.set_xlabel(r'$\rm BTJD$')
ax3.set_ylabel(r'$\rm Raw \ Flux$')

plt.subplots_adjust(hspace=0.5)

pdf.savefig(fig,orientation='portrait')
plt.close()
pdf.close()