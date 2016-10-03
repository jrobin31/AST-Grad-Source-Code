import numpy as np
from matplotlib import pyplot
from progressbar import ProgressBar

stepSize = 1e-3
steps = int(10/stepSize)

phi0 = 0.0
velR0 = -0.2
velZ0 = 0.2
rad0 = 0.5
height0 = 0.7
q = 1.0

## Angular momenta of orbits we will calculate
Lz = 0.2

## Equation for dPhi/dt in Axisymmetric potential
def phiDot(angMomentum, radius):
	return angMomentum / radius**2

## Equation for dv_r/dt in Axisymmetric potential
def velRDot(angMomentum, radius, height):
	return angMomentum**2/radius**3 - radius/(radius**2 + height**2)

## Equation for dv_z/dt in Axisymmetric potential
def velZDot(q, radius, height):
	return -height/q/(radius**2 + height**2)

## Using leapfrog integration, find velocity at t=1/2
velRHalf = velR0 + stepSize/2.*velRDot(Lz, rad0, height0)
velZHalf = velZ0 + stepSize/2.*velZDot(q, rad0, height0)

## Rows are timesteps, columns are values of L
phis = np.zeros((steps,))
radii = np.zeros((steps,))
heights = np.zeros((steps,))
velrs = np.zeros((steps,))
velzs = np.zeros((steps,))

## Set initial conditions for each column (L value)
phis[0] = phi0
radii[0] = rad0
heights[0] = height0
velzs[0] = velZHalf
velrs[0] = velRHalf

## Add a progress bar in the terminal so we know how far along the calculation is
pBar = ProgressBar(maxval=len(phis.flat)).start()

for step in xrange(1,steps):
	phis[step] = phis[step-1] + stepSize*phiDot(Lz, radii[step-1])
	radii[step] = radii[step-1] + stepSize*velrs[step-1]
	heights[step] = heights[step-1] + stepSize*velzs[step-1]/q
	velrs[step] = velrs[step-1] + stepSize*velRDot(Lz, radii[step], heights[step])
	velzs[step] = velzs[step-1] + stepSize*velZDot(q, radii[step], heights[step])
	
	pBar.update(step)
		
pBar.finish()

## find zero crossings
zcUpMask = (np.diff(np.sign(heights), 1) > 0)

plotTitle = r'''$\~v_{r,\tau=0} = %0.2g$, $\~v_{z,\tau=0} = %0.2g$, 
$\~r_{\tau=0} = %0.2g$, $\~z_{\tau=0} = %0.2g$''' % (velR0, velZ0, rad0, height0)

## Plot r vs. phi
fig = pyplot.figure()
ax = pyplot.subplot(111)
ax.plot(radii, heights, label='$\~L_z = %0.2g$, $q = %0.2g$' % (Lz, q))
#pyplot.xlim(0.0,0.5)
ax.set_xlabel('$\~r$')
ax.set_ylabel('$\~z$')
ax.legend(loc='lower right')
ax.set_title(plotTitle)
pyplot.show()

## Plot y vs. x
ax = pyplot.subplot(111, aspect='equal', adjustable='datalim')
x = radii*np.cos(phis)
y = radii*np.sin(phis)
ax.plot(x, y, label='$\~L_z = %0.2g$, $q = %0.2g$' % (Lz, q))
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.legend(loc='upper right')
ax.set_title(plotTitle)
pyplot.show()

## Plot upward zero crossings, r vs. v_r
ax = pyplot.subplot(111)
ax.plot(radii[zcUpMask], velrs[zcUpMask], 'o', label='$\~L_z = %0.2g$, $q = %0.2g$' % (Lz, q))
ax.set_xlabel('$\~r$')
ax.set_ylabel('$\~v_r$')
ax.legend(loc='lower right')
ax.set_title(plotTitle)
pyplot.show()