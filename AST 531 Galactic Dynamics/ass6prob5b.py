import numpy as np
from scipy.integrate import simps
from matplotlib import pyplot as plt
from progressbar import ProgressBar

def phaseDensity(energy):
	return np.exp(-energy)/((2*np.pi)**1.5)

def radPotential(radius):
	return (radius**2)/8.

def radVelocity(momentum, energy, radius, potential):
	return np.sqrt(2*(energy - potential) - momentum**2/radius**2)

def maxMomentum(energy, radius):
	return np.sqrt(2*radius**2*(energy - radPotential(radius)))

def diffDensity(momentum, energy, radius, potential):
	radVel = radVelocity(momentum, energy, radius, potential)
	return 4*np.pi*momentum*phaseDensity(energy)/(radVel*radius**2)

radii = np.linspace(0.1,10,100)
potentials = radPotential(radii)

densities = np.zeros_like(radii)

## Add a progress bar in the terminal so we know how far along the calculation is
pBar = ProgressBar(maxval=densities.shape[0]).start()

for step,(radius,potential) in enumerate(zip(radii,potentials)):
	energies = potential + np.linspace(0.01, 49.99, 100)
	maxMomenta = maxMomentum(energies, radius)
	
	LGrid = np.outer(maxMomenta, np.linspace(0.01, 0.99, 10))
	EGrid = np.outer(energies, np.ones(10))
		
	diffs = diffDensity(LGrid, EGrid, radius, potential)
	
	densities[step] = simps(simps(diffs, LGrid, axis=1), energies)
	
	pBar.update(step)
pBar.finish()

plt.loglog(radii, densities)
plt.xlabel('Radius')
plt.ylabel('Density')
plt.savefig('ass6prob5b.pdf')
plt.show()