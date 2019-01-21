# TESS Alerts

Script to download TESS data and get information from the system by fitting a model light curve. 

## Prerequisites

This routine needs Python 3.6, SciPy, Astropy, NumPy and also the BLS module and the LM fit. SciPy, Astropy and NumPy can be downloaded with

```
conda install astropy numpy scipy
```

whereas the BLS module and the LM fit can be downloaded using

```
pip install bls
```

```
pip install lmfit
```

## Installation

To get the routine, just clone the repository:

```
git clone git@gitlab.com:tessalerts/tess-alerts.git
```

## Usage

First, we need to download the files from the TESS Alert Data website, https://tev.mit.edu/user/login/. Use 

```
ipython download.py 
```

and login typing in your username and your password. download.py downloads the .fits files and a .csv file containing the targets and its parameters. Next, we are ready to run main.py, simply typing

```
ipython main.py 
```

## main.py

Now we will explain step by step what this script does. 

First, we loop through all the candidates downloaded and get their file names and pipeline. We get the flux, the time and the flux error, along with stellar parameters such as effective temperature, *log g* and metallicity (whenever present). The stellar parameters are used to obtain the limb darkening coefficients.

We divide the data into chunks, separating it from before and after the downlink break (this should also separate data from different sectors, although the link to "download all" in the TESS Alert data website seems to be outdated in the time this was written, so we only worked with data from sector 1). We fit a median filter to the data in order to de-trend it and normalize it. Then, a second median filter is applied (with a smaller window) to mimic the real data and bring it all to one. We shift the data divided by this second median filter to zero and apply a cut to take away the outliers. The rule for the cut is that we take out any points outside 5MAD, where MAD is the median of the absolute flux data. 

After removing the outliers, we apply the BLS module to extract parameters from the light curve. We get the period, the duration and the depth. More about the BLS module can be found at https://github.com/dfm/bls.py. 

With the period, duration and depth, we are able to estimate some other parameters for the system, such as the star-planet radius ratio and the ratio between the semi-major axis and the star radius. With a first estimate of these parameters, and assuming initial values for the eccentricity, the inclination and the argument of periapse, along with the limb darkening coefficients obtained before, we generate a light curve using our own light curve generator. 

The parameters for the generated light curve are optimized to give the best chi-squared, and the optimization is done by applying a Levenberg-Marquardt algorithm taken from https://lmfit.github.io/lmfit-py/. The eccentricity and the argument of periapse were maintained fixed. This can be changed in the *fit.py* file, along with the minimum and maximum values the parameters can take. 

When the best fit is obtained, we plot the results, with a plot showing the phase folded light curve with the best model, the parameters obtained, a likelihood plot of the period, obtained from the BLS method, and the raw data plotted with the first median filter applied and dashed lines to show the position of the transits.

## qlp.py

This program does the exact same thing as main.py, but for the data processed by the Quick Look Pipeline (QLP). The reason for this division is that the .fits files from SPOC and from QLP are rather different, so to avoid having a number of *if* statements in our program, we decided to separate it into two different scripts. 

## spoc_single.py and qlp_single.py

These two scripts do the same thing as *main.py* and *qlp.py*, except that they don't loop over all the files, but rather you can manually choose a single candidate to work with and save some time. 