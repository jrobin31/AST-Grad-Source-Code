import numpy as np
from scipy.integrate import simps
from matplotlib import pyplot as plt
from progressbar import ProgressBar

## The phase density of star particles (i.e. fraction of stars at the given energy) 
def phaseDensity(energy):
	return np.exp(-energy)/((2*np.pi)**1.5)

## The radial potential at the given radius
def radPotential(radius):
	return (radius**2)/8.

## The square of radial velocity for the given parameters
def radVelocity(radius, momentum, energy):
	sqrRadVel = 2*(energy - radPotential(radius)) - momentum**2/radius**2
	return np.sqrt(np.clip(sqrRadVel, 0, np.inf))

## The upper limit of integration for momentums
def maxMomentum(energy, radius):
	return np.sqrt(2*radius**2*(energy - radPotential(radius)))

## The upper limit of integration for radial action
def maxRadAction(momentum, energy):
	return 1.1*np.sqrt(4*energy + 2*np.sqrt(4*energy**2 - momentum**2))

def energyForRadAction(radAction, momentum):
	return 0.5*(radAction + momentum)

## Radial action for given parameters
def radAction(momentum, energy):
	momenta = np.outer(momentum, np.ones(100))
	energies = np.outer(energy, np.ones(100))
	rads = np.outer(maxRadAction(momentum, energy), np.linspace(0.01, 0.99, 100))
	rads.shape = momenta.shape = energies.shape = (momentum.shape) +(100,)
	return 2/np.pi*simps(radVelocity(rads, momenta, energies), rads)

## Differential density
def diffDensity(momentum, energy, radius, potential):
	phD = phaseDensity(energyForRadAction(radAction(momentum, energy), momentum))
	radVel = radVelocity(radius, momentum, energy) 
	return 4*np.pi*momentum*phD/(radius**2*radVel)

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
plt.savefig('ass6prob5c.pdf')
plt.show()