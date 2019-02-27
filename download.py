'''
Script for downloading TESS alert targets:
	Will download the .csv file from the alerts 
	and the .fits files	to your current directory

@author: emilk, marcelo
'''
# =============================================================================
# Basic packages 
# =============================================================================
import glob
from subprocess import call
import os

# =============================================================================
# Downloading TESS alerts
# =============================================================================
import requests
from lxml import html

def get_csv(session_requests,csv_url='https://tev.mit.edu/toi/alerts/all/all/csv/',name='./alerts.csv'):
	'''
	Module to create a .csv file of all TESS alert targets.

	Input:
		csv_url : string - link to the .csv file.
		name    : string - name for the output .csv file. 
	'''
	result = session_requests.get(csv_url, headers=dict(referer=csv_url),timeout=5.)
	if result.ok == True:
		with open(name,'wb') as file:
			file.write(result.content)
	else:
		print('something went wrong - could not download .csv file\n')
	return

def get_fits(session_requests,fits_url='https://tev.mit.edu/toi/alerts/all/all/pkg/fits/'):
	'''
	Module to create the .fits files for all the TESS alert targets.

	Input:
		fits_url: string - link to the .fits files.
	'''
	result = session_requests.get(fits_url, headers=dict(referer=fits_url),timeout=5.)
	tree = html.fromstring(result.text)
	dl_link = list(set(tree.xpath("//a[@class='btn btn-primary']/@href")))[0]
	dl_url = 'https://tev.mit.edu' + dl_link
	result = session_requests.get(dl_url, headers=dict(referer=dl_url),timeout=5.)
	if result.ok == True:
		with open('./tess.tar.bz2','wb') as file:
			file.write(result.content)
	else:
		print('something went wrong - could not download .fits files\n')
	return

def login(login_url='https://tev.mit.edu/user/login/',csv=False,csv_name='./alerts.csv',fits=False,fits_dir='./fits'):
	'''
	Module that logs into the TESS alert site in order to download either the .csv file containing the targets
	or the actual light curves of the targets.

	Input:
		login_url : string - the url to the TESS webpage that requires login.
		csv       : bool   - if True .csv file will be downloaded.
		csv_name  : string - the name of the .csv file.		
		fits       : bool   - if True .fits files will be downloaded.
		fits_dir  : string - name of the directory with the .fits files.
	'''
	## User login, password is hidden
	import getpass
	print('Login to {}\n'.format(login_url))
	username = input('Please type in your username: ')
	password = getpass.getpass('Please type in your password: ')
	login = {'username' : username, 'password' : password, 
			 'csrfmiddlewaretoken': None}

	session_requests = requests.session()
	result = session_requests.get(login_url,timeout=5.)
	
	## This is something that is needed for some reason 
	tree = html.fromstring(result.text)
	authenticity_token = list(set(tree.xpath("//input[@name='csrfmiddlewaretoken']/@value")))[0]
	login['csrfmiddlewaretoken'] = authenticity_token
	
	logged_in = False
	while logged_in == False:
		result = session_requests.post(login_url,data=login,headers=dict(referer=login_url))
		if 'Please enter a correct username and password.' in result.text:
			print('Please enter a correct username and password.')
			login['username'] = input('Username: ')
			login['password'] = getpass.getpass('Password: ')
		else:
			print('Login successful')
			logged_in = True

	if csv:
		get_csv(session_requests,name=csv_name)
	if fits:
		call(['mkdir',fits_dir])
		cwd = os.getcwd()
		os.chdir(fits_dir)
		get_fits(session_requests)
		call(['tar','-xvjf','./tess.tar.bz2'])
		os.remove('./tess.tar.bz2')
		gzfiles = glob.glob('*.fits.gz')
		for gzfile in gzfiles:
			call(['gzip','-d',gzfile])
		os.chdir(cwd)

	return

def download(login_url='https://tev.mit.edu/user/login/',csv_name='./alerts.csv', fits_dir='./fits'):
	'''
		Module that runs the previous ones and downloads the .fits files.

		Input
	'''
	fits_files = glob.glob('{}/*.fits'.format(fits_dir))
	if not fits_files:
		login(fits=True,fits_dir=fits_dir)

download()