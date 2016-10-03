import numpy as np
## simps is Simpson's Rule integration
from scipy.integrate import simps
from matplotlib import pyplot as plt
from progressbar import ProgressBar

bhMass = 5.0

## The phase density of star particles (i.e. fraction of stars at the given energy) 
def phaseDensity(energy):
	return np.exp(-energy)/((2*np.pi)**1.5)

## The radial potential at the given radius
def radPotential(radius):
	return (radius**2)/8. - bhMass/radius

## The radial velocity for the given parameters. Returns sqrRadVel < 0
def radVelocity(radius, momentum, energy):
	sqrRadVel = 2*(energy - radPotential(radius)) - momentum**2/radius**2
	return np.sqrt(np.clip(sqrRadVel, 0, np.inf))

## The upper limit of integration for momenta
def maxMomentum(energy, radius):
	return np.sqrt(2*(energy - radPotential(radius))*radius**2)

## The upper limit of integration for radial action
def calcApoapsis(momentum, energy):
	# Note: last term, 40*bhMass, causes significant overestimation of actual 
	# apoapsis. Exact solution to L^2 = 2*r^2*(E-phi) is ridiculous for this 
	# potential. Instead I calculate this as a first guess, then look for max 
	# radius where v_r(E,L,r) is nonzero to get a more accurate estimation
	return 1.1*np.sqrt(4*energy + 2*np.sqrt(4*energy**2 - momentum**2 + 40*bhMass))

## Radial action for given momentum and energy
def radAction(momentum, energy):
	EGrid = np.outer(energy, np.ones(100))
	LGrid = np.outer(momentum, np.ones(100))
	## Calculate first guess at apoapses
	radGrid = np.outer(calcApoapsis(momentum, energy), np.linspace(0.01, 0.99, 100))
	radGrid.shape = LGrid.shape = EGrid.shape = momentum.shape +(100,)
	radVels = radVelocity(radGrid, LGrid, EGrid)
	## Pick new apoapses based on max value where velocity is nonzero
	maxRadii = np.amax(radGrid*np.sign(radVels), axis=2)
	radGrid = np.outer(maxRadii, np.linspace(0.01, 1.1, 100))
	radGrid.shape = EGrid.shape
	return 2/np.pi*simps(radVelocity(radGrid, LGrid, EGrid), radGrid)

## Back out the energy for a given radial action
def energyForRadAction(radAction, momentum):
	return 0.5*(radAction + momentum)

## Differential density
def diffDensity(momentum, energy, radius):
	phDens = phaseDensity(energyForRadAction(radAction(momentum, energy), momentum))
	radVel = radVelocity(radius, momentum, energy) 
	return 4*np.pi*momentum*phDens/(radVel*radius**2)

## Ensure exception is raised for e.g. bad values in sqrt or divide by 0
np.seterr(all='raise')

radii = np.linspace(0.1,10,100)
potentials = radPotential(radii)
densities = np.zeros_like(radii)

## Add a progress bar in the terminal so we know how far along the calculation is
pBar = ProgressBar(maxval=densities.shape[0]).start()

for step,(radius,potential) in enumerate(zip(radii,potentials)):
	energies = potential + np.linspace(0.01, 49.99, 100)
	maxMomenta = maxMomentum(energies, radius)
	
	EGrid = np.outer(energies, np.ones(10))
	LGrid = np.outer(maxMomenta, np.linspace(0.01, 0.99, 10))
		
	diffs = diffDensity(LGrid, EGrid, radius)
	
	densities[step] = simps(simps(diffs, LGrid), energies)
	
	pBar.update(step)
pBar.finish()

## Save out data file so we can plot all on same axes later
np.savetxt('ass6prob5d_%.1f.txt' % bhMass, np.vstack((radii,densities)).T)

plt.loglog(radii, densities)
plt.xlabel('Radius')
plt.ylabel('Density')
plt.savefig('ass6prob5d.pdf')
plt.show()