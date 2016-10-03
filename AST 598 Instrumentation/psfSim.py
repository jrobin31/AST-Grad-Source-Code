"""
AST 598 Instrumentation
Homework #4
Generated PSFs from synthetic apertures
"""
import numpy as np
from numpy.fft import fftshift, fftn
from matplotlib import pyplot
from matplotlib import cm
from scipy.spatial.distance import pdist

arcsecPerPx = 0.01
wavelength = 0.5e-6 # Wavelength in m
nPixels = 2048
mPerPx = wavelength/nPixels/arcsecPerPx*206265
dishDiamPx = 0.3/mPerPx # 0.3 m converted to px

imgCenter = nPixels/2 - 0.5
coordsX,coordsY = np.meshgrid(range(nPixels),range(nPixels))

## Add dishes at given x,y centers
def addDishes(arr, centersX, centersY):
	for (centerX,centerY) in zip(centersX, centersY):
		arr[(coordsX-centerX)**2 + (coordsY-centerY)**2 < (0.5*dishDiamPx)**2] = 1
	return arr

## Calculate the longest baseline in the dish system
def calcLongestBaseline(centersX, centersY):
	if(len(centersX) < 2 or len(centersY) < 2):
		return dishDiamPx
	else:
		return dishDiamPx + np.sqrt(np.max(pdist(np.vstack((centersX,centersY)).T, 
											metric='sqeuclidean')))

## Calculate the Rayleigh Criterion in pixels, given diameter in pixels
def calcRayleighPx(diam):
	return 1.22 * wavelength / (diam * mPerPx) * 206265 / arcsecPerPx

## Calculated the FFT power image, normed to 1
def fftPower(arr):
	pow = fftshift(fftn(arr))
	pow = np.abs(pow)**2
	pow *= 1.0/np.max(pow)
	return pow
	
## Plot a sqrt-stretched image with a colorbar
def plotImage(img, dispRect=None):
	if(dispRect == None):
		dispRect = [0, img.shape[1], 0, img.shape[1]]
	psp = pyplot.imshow(np.sqrt(img), cmap=cm.gray_r)
	pyplot.axis(dispRect)
	pyplot.colorbar(psp)
	pyplot.show()
	
## Plot horizontal and vertical slice through the given image
def plotSlices(img, lines=[]):
	sz = img.shape[0]
	pyplot.plot(range(sz), img[1024,:], drawstyle='steps-mid', 
			label='Horizontal', color='RoyalBlue')
	pyplot.plot(range(sz), img[:,1024], drawstyle='steps-mid', 
			label='Vertical', color='SeaGreen')
	## Plot vertical lines (for Rayleigh Criterion) and text label
	pyplot.vlines(lines, 0.0, 1.0, color='coral')
	for lineX in lines:
		lineXArcsec = (lineX-imgCenter)*arcsecPerPx
		pyplot.text(lineX+2, 0.95, r'$\theta_R=%.2f^{\prime\prime}$' % lineXArcsec)
	## Width of axes is twice the width of 1-dish Rayleigh angle
	axExtent = 2*np.max(lines) - sz/2 - 0.5 
	pyplot.axis((sz - axExtent, axExtent, 0.0, 1.0))
	pyplot.legend(loc='upper left')
	pyplot.show()

## We will make a dictionary that looks like:
## name => (dish center x positions, dish center y positions)

centers = dict()
## Single dish aperture
centers['single'] = (imgCenter*np.ones((1,1)), 
					imgCenter*np.ones((1,1)))

## 2 dishes separated by 1m along x
centers['double'] = (np.linspace(-0.5,0.5,2)/mPerPx + imgCenter, 
					imgCenter*np.ones((2)))

## 7 dishes separated by 1m each along x
centers['7line'] = (np.linspace(-3.0,3.0,7)/mPerPx + imgCenter,
					imgCenter*np.ones((7)))

## 13 dishes in a '+' shape
radii = np.linspace(0.0,3.0,4)/mPerPx
angles = np.linspace(0.0,3./2.*np.pi,4)
centers['13plus'] = (np.outer(radii, np.cos(angles)).flatten() + imgCenter,
					np.outer(radii, np.sin(angles)).flatten() + imgCenter)

## 13 dishes in a 'Y' shape
radii = np.linspace(0.0,4.0,5)/mPerPx
angles = np.linspace(0.0,4./3.*np.pi,3)
centers['13y'] = (np.outer(radii, np.cos(angles)).flatten() + imgCenter,
					np.outer(radii, np.sin(angles)).flatten() + imgCenter)

## 13 dishes in a gaussian distribution
centers['13gaussian'] = (np.random.normal(0.0, 1.0, size=(13))/mPerPx + imgCenter,
						np.random.normal(0.0, 1.0, size=(13))/mPerPx + imgCenter)

diamOneDish = calcLongestBaseline(centers['single'][0], centers['single'][1])

## Synthesize aperture, do fft, and make plots for each case above
for (centersX,centersY) in centers.values():
	aper = np.zeros((nPixels,nPixels), dtype='int')
	
	aper = addDishes(aper, centersX, centersY)
	
	diams = np.unique([diamOneDish, calcLongestBaseline(centersX, centersY)])
	
	plotImage(aper)
	
	psf = fftPower(aper)
	
	plotImage(psf, dispRect=[896,1153,896,1153])
	plotSlices(psf, lines=calcRayleighPx(diams) + imgCenter)
