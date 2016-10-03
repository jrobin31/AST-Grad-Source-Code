import numpy as np
from matplotlib import pyplot
from matplotlib import cm

stepSize = 0.005
plotTimes = np.linspace(0.0, 8.0, 9)
endTime = plotTimes[-1]
steps = int(endTime / stepSize) + 1

## Calculates the radius array R_{i,j} for given x,y position arrays 
def radius(posX, posY):
	return np.sqrt(posX**2 + posY**2)

## Calculates Omega(R) array for given radius array
def omega(radius):
	return np.sqrt(10.0 / radius)

## Calculates X acceleration array for given x,y position arrays  
def accelX(posX, posY):
	return -omega(radius(posX, posY))**2 * posX

## Calculates Y acceleration array for given x,y position arrays	
def accelY(posX, posY):
	return -omega(radius(posX, posY))**2 * posY

## Color stars from a certain region differently
colors = np.zeros((61,61), dtype='float')
colors[25:37,40:51] = 0.9

## Make a plot of all stars within radius 25.0. Note: velX/Y not currently used
def plotStars(posX, posY, velX, velY, time):
	radMask = (radius(posX, posY) < 25.0)
	
	ax = pyplot.subplot(111, aspect='equal', adjustable='datalim')
	ax.scatter(posX[radMask], posY[radMask], c=colors[radMask], s=5, 
			cmap=cm.spectral, vmax=1.0, edgecolors='none')
	ax.set_xlabel('x position')
	ax.set_ylabel('y position')
	ax.set_title(r'$\tau=%.1f$' % time)
	pyplot.show()

## Construct initial positions and velocities
posX,posY = np.meshgrid(np.linspace(-30, 30, 61), np.linspace(-30, 30, 61))
posX += np.random.uniform(-0.25, 0.25, size=(61,61))
posY += np.random.uniform(-0.25, 0.25, size=(61,61))
velX = -omega(radius(posX, posY)) * posY
velY = omega(radius(posX, posY)) * posX

## Advance velocities a half step
velX += stepSize/2. * accelX(posX, posY)
velY += stepSize/2. * accelY(posX, posY)

for step in xrange(steps):
	## First see if we need to plot this step
	if(step*stepSize in plotTimes):
		plotStars(posX, posY, velX, velY, step * stepSize)
	
	## Now advance positions and velocities
	posX += stepSize * velX
	posY += stepSize * velY
	velX += stepSize * accelX(posX, posY)
	velY += stepSize * accelY(posX, posY)
