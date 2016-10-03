import numpy as np
from numpy.linalg import norm
from matplotlib import pyplot

nfwDispersion = 200.0 ## km/s
nfwRadius = 5.0 ## kpc
cellSize = 0.5 ## kpc

cgs_G = 6.674e-8 ## G in cgs (cm^3 g^-1 s^-2). Doesn't matter, factors out.

nfwRho = nfwDispersion**2 / (4 * np.pi * cgs_G * nfwRadius**2) 

def nfwPotential(position):
	'''
	Calculates Navarro-Frenk-White gravitational potential at a given position 
	in cell coordinates.
	Note: This is slow since the looping is not numpyized. But it's a homework
	problem so that's not particularly important
	@param position: 3-tuple cell index to calculate potential
	'''
	position = np.array(position) ## Turn tuple into vector
	total = 0.0
	for cell in np.ndindex((61,61,61)):
		cell = np.array(cell)
		cellDist = physicalDist(position,cell)
		## Only accumulate for cells that aren't the test cell
		if(cellDist > 0.0):
			total += nfwDensity(cell) / cellDist
	return -cgs_G * cellSize**3 * total

def nfwDensity(position):
	'''
	Calculates the Navarro-Frenk-White density at a given position in cell 
	coordinates
	@param position: 3-tuple cell index to calculate density
	'''
	radMag = physicalDist(position, (30,30,30))
	if(radMag > 0.0):
		return  nfwRho * nfwRadius / (radMag * (1 + radMag / nfwRadius)**2)
	else:
		return 0.0
	
def analyticNFWPotential(position):
	'''
	Calculates the Navarro-Frenk-White potential analytically
	@param position: 3-tuple cell index to calculate potential 
	'''
	radMag = physicalDist(np.array(position), (30,30,30))
	return -nfwDispersion**2 * np.log(1+radMag/nfwRadius) * nfwRadius / radMag

def physicalDist(cell1, cell2):
	'''
	Finds the physical (kpc) distance between cell1 and cell2
	@param cell1: numpy array of first cell index
	@param cell2: numpy array of second cell index
	'''
	return cellSize * norm(cell1 - cell2, 2) ## explicity use 2-norm

radii = np.arange(30,61)
potentials = np.zeros(radii.shape, dtype=float)
analPotentials = np.zeros(radii.shape, dtype=float)
for i,rad in enumerate(radii):
	print 'calculating radius %i...' % rad,
	potentials[i] = nfwPotential((rad,30,30))
	analPotentials[i] = analyticNFWPotential((rad,30,30))
	print '%2f (numeric) %2f (analytic)' % (potentials[i], analPotentials[i])

## Plot everything up
#fig = pyplot.figure()
#ax = pyplot.subplot(111)
#ax.plot((radii-30)*cellSize,potentials, label='Numeric $\Phi$')
#ax.plot((radii-30)*cellSize,analPotentials, label='Analytic $\Phi$')
#ax.plot((radii-30)*cellSize,potentials-8650, label='Numeric $\Phi-8650$')
#ax.legend()
#pyplot.show()