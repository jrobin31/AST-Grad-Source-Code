import numpy as np
from matplotlib import pyplot
from progressbar import ProgressBar

stepSize = 1e-3
steps = 20000

phi0 = 0.0
vel0 = 0.0
rad0 = 1.0

## Angular momenta of orbits we will calculate
angMomentaDivisors = np.array([4., 2., 1.5, 1.2, 1.1, 1.05]) 
angMomenta = 1./angMomentaDivisors

## Equation for dPhi/dt in Kelperian potential
def phiDot(angMomentum, radius):
	return angMomentum / radius**2

## Equation for dr/dt in Kelperian potential
def velDot(angMomentum, radius):
	return angMomentum**2/radius**3 - 1./radius**2 

## Using leapfrog integration, find velocity at t=1/2
velHalfs = vel0 + stepSize/2.*velDot(angMomenta, rad0)

## Rows are timesteps, columns are values of L
phis = np.zeros((steps,len(angMomenta)))
radii = np.zeros((steps,len(angMomenta)))
vels = np.zeros((steps,len(angMomenta)))

## Set initial conditions for each column (L value)
phis[0,:] = phi0
radii[0,:] = rad0
vels[0,:] = velHalfs

## Add a progress bar in the terminal so we know how far along the calculation is
pBar = ProgressBar(maxval=len(phis.flat)).start()

for amStep, angMomentum in enumerate(angMomenta):
	for step in xrange(1,steps):
		phis[step,amStep] = phis[step-1,amStep] + stepSize*phiDot(angMomentum, radii[step-1,amStep])
		radii[step,amStep] = radii[step-1,amStep] + stepSize*vels[step-1,amStep]
		vels[step,amStep] = vels[step-1,amStep] + stepSize*velDot(angMomentum, radii[step,amStep])
		
		## Early out once we have completed 2*pi in angle
		if(phis[step,amStep] > 2*np.pi):
			phis[step+1:,amStep] = np.NAN
			break
		pBar.update(step*amStep)
		
pBar.finish()

## Plot r vs. phi
fig = pyplot.figure()
ax = pyplot.subplot(111)
for amStep,divisor in enumerate(angMomentaDivisors):
	ax.plot(phis[:,amStep], radii[:,amStep], label='$\~L = \~L_c/%0.3g$' % divisor)
ax.set_xlabel('$\phi$')
ax.set_ylabel('$\~r$')
ax.set_title('Orbital Parameters, apocenter $\~r=1$')
ax.legend(loc='lower right')
pyplot.show()

## Plot y vs. x
ax = pyplot.subplot(111, aspect='equal')
for amStep,divisor in enumerate(angMomentaDivisors):
	x = radii[:,amStep]*np.cos(phis[:,amStep])
	y = radii[:,amStep]*np.sin(phis[:,amStep])
	ax.plot(x, y, label='$\~L = \~L_c/%0.3g$' % divisor)
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.set_title('Orbits, apocenter $\~r=1$')
ax.legend(loc='upper right')
pyplot.show()