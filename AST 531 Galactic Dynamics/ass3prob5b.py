import numpy as np
from matplotlib import pyplot
from progressbar import ProgressBar

stepSize = 1e-3
steps = 20000

phiStart = 0.0
rStart = 1 - 1e-6

## Angular momenta of orbits we will calculate
angMomentaDivisors = np.array([4., 2., 1.5, 1.2, 1.1, 1.05]) 
angMomenta = 1/angMomentaDivisors

## Equation for dPhi/dt
def phiDot(angMomentum, radius):
	return angMomentum / radius**2

## Equation for dr/dt
def rDot(angMomentum, radius):
	return -np.sqrt(2*(1/radius - 1 - angMomentum**2/2 * (1/radius**2 - 1)))

phis = np.zeros((steps,len(angMomenta)))
radii = np.zeros((steps,len(angMomenta)))

phis[0,:] = phiStart
radii[0,:] = rStart

## Add a progress bar in the terminal so we know how far along the calculation is
pBar = ProgressBar(maxval=len(phis.flat)).start()

for amStep, angMomentum in enumerate(angMomenta):
	for step in xrange(1,steps):
		radii[step,amStep] = radii[step-1,amStep] + stepSize*rDot(angMomentum, radii[step-1,amStep])
		phis[step,amStep] = phis[step-1,amStep] + stepSize*phiDot(angMomentum, radii[step-1,amStep])
		pBar.update(step*amStep)
pBar.finish()

## Plot r vs. phi
fig = pyplot.figure()
ax = pyplot.subplot(111)
for amStep,divisor in enumerate(angMomentaDivisors):
	ax.plot(phis[:,amStep], radii[:,amStep], label='$\~L = 1/%0.3g$' % divisor)
ax.axvline(np.pi)
ax.set_xlabel('$\phi$')
ax.set_ylabel('$\~r$')
ax.set_title('Keplerian Orbital Parameters')
ax.legend(loc='lower left')
pyplot.show()

## Plot y vs. x
ax = pyplot.subplot(111, aspect='equal')
for amStep,divisor in enumerate(angMomentaDivisors):
	x = radii[:,amStep]*np.cos(phis[:,amStep])
	y = radii[:,amStep]*np.sin(phis[:,amStep])
	ax.plot(x, y, label='$\~L = 1/%0.3g$' % divisor)
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.set_title('Keplerian Orbits')
ax.legend(loc='upper right')
pyplot.show()