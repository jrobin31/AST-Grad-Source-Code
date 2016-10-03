from pylab import *
import readcol

## Minimize the cumulative distance between the two data sets
def minimizeDist():
	bestdist = 100000
	bestAv = 0
	for i in range(0,200):
		Av = i/200.0*3.0
		dists = getDists(Av)
		dist = sum(dists)
		if(dist < bestdist):
			bestAv = Av
			bestdist = dist
			movable.set_xdata(colBV - eBV*Av)
			movable.set_ydata(colUB - eUB*Av)
 			plothist(dists, Av)
			draw()
	
	print(bestAv)
	return

## Calculate a list of minimum distances between two data sets
def getDists(Av):
	dists = []
	for x1,y1 in zip(colBV, colUB):
		minDist = 10000
		for x2,y2 in zip(colBV0[mask],colUB0[mask]):
			dist = sqrt((x1 - eBV*Av - x2)**2 + (y1 - eUB*Av - y2)**2)
			if(dist < minDist):
				minDist = dist
		dists.append(minDist)
	
	return dists

## Plot a histogram of the distance data
def plothist(data, Av):
	cla()
	hist(data, 100,range=(0,1))
	title("Distance to nearest MS star")
	xlabel("Distance")
	ylabel("Frequency")
	text(0.6, 9, '$A_V=%0.2f$' % Av)
	axis([0,1,0,10])
	return	

## Read in data files
stars = readcol.readcol("cluster.txt", asStruct=True)
stars0 = readcol.readcol("intrinsic_props.txt", asStruct=True, skipline=2)
redlaw = readcol.readcol("extlaw.txt", asStruct=True, skipline=1)

## Renormalize reddeding law, calculate color excesses
rv31 = redlaw.alam_aj_31 / redlaw.alam_aj_31[where(redlaw.wl == 0.55)][0]
eBV = rv31[where(redlaw.wl == 0.44)][0] - rv31[where(redlaw.wl == 0.55)][0]
eUB = rv31[where(redlaw.wl == 0.365)][0] - rv31[where(redlaw.wl == 0.44)][0]

## U-Band values were wrong in the file, should be 0.1 mag fainter
colUB = (stars.Un + 0.1) - stars.Bn
colBV = stars.Bn - stars.Vn

colUB0 = stars0.UB
colBV0 = stars0.BV

mask = (colBV0 != -99) & (colUB0 != -99)

ax1 = subplot(121)

movable, = plot(colBV, colUB, marker="o", linestyle="none", color=(1.0,0.5,0.5))
plot(colBV0[mask], colUB0[mask], marker="o", linestyle="none", color=(0.5,0.5,1.0))

title("UV-Visible Colors")
xlabel(r"$B-V$", fontsize="large")
ylabel(r"$U-B$", fontsize="large")

## Flip y axis
ylim(ylim()[::-1])

## Reddening vector
text(0.75, -1, '$R_V=3.1$', ha="right", va="bottom")
annotate('', xytext=(0.75,-1), xy=(0.75+eBV,-1+eUB), 
		arrowprops=dict(arrowstyle="-|>", color='red'))

ax2 = subplot(122)
plothist(getDists(0), 0)

savefig("ColColExtUncorrected.pdf")

ax1.set_title("Extinction Corrected")
minimizeDist()

savefig("ColColExtCorrected.pdf")

show()