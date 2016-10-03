import numpy as np
from numpy.linalg import norm
from scipy.spatial.distance import cdist
from matplotlib import pyplot

## Horizontal and vertical cell sizes
hCellSize = 0.5
vCellSize = 0.167

astro_G = 4.302e-6 ## G in kpc MSun^-1 (km/s)^2

## Density parameters
r_0Disk = 3.0
z_0Disk = 0.3
rho_0Disk = 6e10 / (4 * np.pi * r_0Disk**2 * z_0Disk)
r_0Bulge = 1.4
rho_0Bulge = 1e10 / (2 * np.pi * r_0Bulge**3)

## 61*61*181 by 3 array containing physical x,y,z positions of each cell 
allPositions = np.array([(x,y,z) for x in range(61) for y in range(61) 
						for z in range(181)], dtype=float)
allPositions = (allPositions - (30,30,90))*(hCellSize,hCellSize,vCellSize)

def calcPotential(position):
	'''
	Calculates gravitational potential at a given position in cell coordinates.
	@param position: 3-tuple cell index to calculate potential
	'''
	physPos = (np.array(position) - (30,30,90))*(hCellSize,hCellSize,vCellSize)
	physPos.shape = (1,3)
	distances = cdist(physPos,allPositions,'euclidean')
	mask = (distances.flat > 0.0)
	return -astro_G * hCellSize**2 * vCellSize * np.sum(density.flat[mask]/distances.flat[mask])

def vCirc(radius, potentials):
	'''
	Calculates circular velocity at radius, given list of potentials
	@param radius: radius where we calculate Vcirc, in cell coordinates (should be n+0.5)
	@param potentials: list of potentials at each radius. must have at least n+1 entries
	'''
	radius = int(radius - 0.5)
	return np.sqrt((potentials[radius + 1] - potentials[radius]) * radius)

# I don't buy that this actually takes much memory. Couple MB or so.
# So we'll try this a more Pythonic way this time
print 'Creating density array...'
density = np.zeros((61,61,181))

## This takes a while since I'm looping over numpy arrays
for i in np.ndindex(density.shape):
	rCyl = norm(i[0:2] - np.array((30,30)), 2) * hCellSize
	zCyl = (i[2] - 90) * vCellSize
	rSph = norm((i - np.array((30,30,90))) * (hCellSize,hCellSize,vCellSize), 2)
	density[i] += rho_0Disk * np.exp(-rCyl / r_0Disk) * np.exp(-np.abs(zCyl) / z_0Disk)
	density[i] += rho_0Bulge / (rSph / r_0Bulge * (rSph / r_0Bulge + 1)**3)

# Replace singular point at the center with one of its neighbors 
density[30,30,90] = density[30,30,91]

totalMass = hCellSize**2 * vCellSize * sum(density) # 6.87e10

radii = np.arange(30,61)
potentials = np.zeros(radii.shape, dtype=float)
for i,rad in enumerate(radii):
	potentials[i] = calcPotential((rad,30,90))

pyplot.plot(radii,potentials)
pyplot.xlabel('Radius')
pyplot.ylabel('Potential')
pyplot.show()	

vRadii = np.linspace(0,29,30) + 0.5
vCircs = np.zeros(vRadii.shape, dtype=float) 
for i,rad in enumerate(vRadii):
	vCircs[i] = vCirc(rad, potentials)

pyplot.plot(vRadii,vCircs)
pyplot.xlabel('Radius')
pyplot.ylabel('Circular Velocity')
pyplot.show()

heights = np.arange(90,181)
zPotentials = np.zeros(heights.shape, dtype=float)
for i,height in enumerate(heights):
	zPotentials[i] =  calcPotential((46,30,height))

pyplot.plot(heights,zPotentials)
pyplot.xlabel('Height Above Disk')
pyplot.ylabel('Potential')
pyplot.show()