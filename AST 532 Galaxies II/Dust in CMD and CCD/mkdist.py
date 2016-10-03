from pylab import *
import readcol

Av = 1.185

## Minimize the cumulative distance between the two data sets
def minimizeDist():
	bestdist = 100000
	bestDmod = 0
	for i in range(0,500):
		dmod = i/500.0*12.0
		dists = getDists(dmod)
		dist = sum(dists)
		if(dist < bestdist):
			bestDmod = dmod
			bestdist = dist
			movable.set_ydata(Vmag - dmod)
			textAnnot.set_text("Dist Mod: %0.2f\nDist: %00.2f pc" % (dmod, 10**((dmod+5)/5)))
			draw()
	
	print(bestDmod)
	return

## Calculate a list of minimum distances between two data sets
def getDists(dmod):
	dists = []
	for x1,y1 in zip(colBV, Vmag):
		minDist = 100000
		for x2,y2 in zip(colBV0[mask], Vmag0[mask]):
			dist = sqrt((y1 - y2 - dmod)**2 + (x1-x2)**2)
			if(dist < minDist):
				minDist = dist
		dists.append(minDist)
	
	return dists

stars = readcol.readcol("cluster.txt", asStruct=True)
stars0 = readcol.readcol("intrinsic_props.txt", asStruct=True, skipline=2)

eBV = Av/3.1

colBV = (stars.Bn - stars.Vn)
Vmag = stars.Vn

colBV0 = stars0.BV
Vmag0 = stars0.M_v

mask = (colBV0 != -99) & (Vmag0 != -99)

#subplot(121)

movable, = plot(colBV, Vmag, marker="o", linestyle="none", color=(1.0,0.5,0.5))
plot(colBV0[mask], Vmag0[mask], marker="o", linestyle="none", color=(0.5,0.5,1.0))

title("Color-Mag Diagram")
xlabel(r"$B-V$", fontsize="large")
ylabel(r"V Magnitude", fontsize="large")

## Flip y axis
ylim(ylim()[::-1])

textAnnot = text(1,-3,"Dist Mod: %0.2f\nDist: %00.2f pc" % (0, 0))

savefig("ColMagUncorrected.pdf")

## Apply extinction to values we read in
colBV = colBV - eBV
Vmag = Vmag - Av
movable.set_xdata(colBV)
movable.set_ydata(Vmag)
title("Extinction Corrected ($A_V=%0.2f, R_V=3.1$)" % Av)

savefig("ColMagExtCorrected.pdf")

title("Distance Corrected ($A_V=%0.2f, R_V=3.1$)" % Av)
minimizeDist()

savefig("ColMagDistCorrected.pdf")

show()