'''
	Script for applying a median filter to a one dimensional array

@author: emilk and marcelo
'''

import numpy as np

def my_median_filter(arr,win):
	'''
	Applies a median filter to the one dimensional array given.

	Input:
		arr 	: array       - array containing the data points.
		win 	: odd integer - the size of the window that the median will be calculated.

	Output:
		myfilt  : array 	  - array with the median calculated from the data points. Has
								the same size as "arr".
	'''
	mat = np.ndarray(shape=3*len(arr))
	mat[:] = 0
	mat[0:len(arr)] = arr[::-1] # Putting the array backwards before
	mat[len(arr):2*len(arr)] = arr
	mat[2*len(arr):3*len(arr)] = arr[::-1] #Putting the array backwards after
	myfilt = np.ndarray(shape=len(arr))
	n = int((win-1)/2)
	for i in range(0,len(arr)):
		myfilt[i] = np.median(mat[len(arr)+i-n:len(arr)+i+n])
	return myfilt