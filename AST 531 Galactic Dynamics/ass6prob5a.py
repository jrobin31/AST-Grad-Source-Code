import numpy as np
from scipy.integrate import simps
from matplotlib import pyplot as plt
from progressbar import ProgressBar

def phaseDensity(energy):
	return np.exp(-energy)/((2*np.pi)**1.5)

def radPotential(radius):
	return (radius**2)/8.

def diffDensity(energy, potential):
	return 4*np.pi*(2*(energy - potential))**0.5*phaseDensity(energy)

radii = np.linspace(0.1,10,100)
potentials = radPotential(radii)

densities = np.zeros_like(radii)

## Add a progress bar in the terminal so we know how far along the calculation is
pBar = ProgressBar(maxval=densities.shape[0]).start()

for step,potential in enumerate(potentials):
	energies = np.linspace(potential, potential+50, 100)
	densities[step] = simps(diffDensity(energies, potential), energies)
	pBar.update(step)
pBar.finish()

plt.loglog(radii, densities)
plt.xlabel('Radius')
plt.ylabel('Density')
plt.savefig('ass6prob5a.pdf')
plt.show()