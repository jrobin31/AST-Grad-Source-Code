"""
Investigates the effect of Malmquist Bias on apparent mag and distance
"""
import numpy as np
import math

def dmod(distance=10):
	return 5*math.log10(d/10)
def idmod(dm):
	return 10*10**(dm/5)
def dist(x,y,z):
	return math.sqrt(x**2 + y**2 + z**2)
def randdeltam():
	return -0.87 * math.log10((np.random.randint(1,7) + 0.5)/6.0)
def calcmetals(deltam):
	return 10**(deltam/-0.87)*0.018 ## Z of the sun

print "Generating 1000 random stars in 260x260x260pc cube..."

## Generate random abs magnitudes. 4.8 is M_V of sun
absmag = np.array(0.3 * np.random.randn(1000) + 4.8)

## Generate random positions
pos = zip(130 * np.random.sample(1000), 130 * np.random.sample(1000), 130 * np.random.sample(1000))

## Calculate distances
r = np.array([dist(*tup) for tup in pos])

## Generate region masks, print region counts
rega = (70 < r) & (r < 90)
regb = (90 < r) & (r < 110)
regc = (110 < r) & (r < 130)
regall = rega | regb | regc
print("Region A:\t%(a)i\nRegion B:\t%(b)i\nRegion C:\t%(c)i\nAll Regions:\t%(all)i" % {'a': len(r[rega]),'b': len(r[regb]),'c': len(r[regc]), 'all': len(r[regall])})

## Calculate apparent magnitudes
appmag = np.array([dmod(d) + absmag[i] for i,d in enumerate(r)])

## Sample all stars with apparent mag < 10 that fall in region A, B, or C
sample = (appmag < 10) & regall

## Print mean abs mags
print("Mean abs mag (all regions):\t\t%(m)f mag" % {'m': np.mean(absmag[regall])})
print("Mean abs mag (sample):\t\t\t%(m)f mag" % {'m': np.mean(absmag[sample])})

## Calculate apparent distance, assuming all sample stars have M_V = 4.8
appdist = np.array([idmod(m - 4.8) for m in appmag])

## Print distances
print("Mean sample dist:\t\t\t%(m)f pc" % {'m': np.mean(r[sample])})
print("Mean apparent dist:\t\t\t%(d)f pc" % {'d': np.mean(appdist[sample])})

## Generate metal-poor magnitudes
mpabsmag = np.array([mag + randdeltam() for mag in absmag])
mpappmag = np.array([dmod(d) + mpabsmag[i] for i,d in enumerate(r)])

## Sample all stars with apparent mag < 10 that fall in region A, B, or C
mpsample = (mpappmag < 10) & regall

## Calculate metallicities
metals = np.array([calcmetals(mpabsmag[i] - absmag[i]) for i in range(len(absmag))])

## Print metallicities
print("Metal fraction (all regions):\t\t%(m)f" % {'m': np.mean(metals[regall])})
print("Metal fraction (sample):\t\t%(m)f" % {'m': np.mean(metals[mpsample])})
print("Metal fraction (sample, regions B&C):\t%(m)f" % {'m': np.mean(metals[mpsample & (regb | regc)])})
