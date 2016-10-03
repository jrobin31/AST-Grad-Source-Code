from pylab import *
import readcol

labelOS = [(-0.1,-0.1),(0.05,0.06),(-0.15,-0.05),(0.05,0.1),(-0.15,-0.05),(0.0,0.12),
		(-0.15,-0.05),(0.05,0.1),(-0.15,-0.05),(0.05,0.1),(-0.15,-0.05),(0.05,0.1),
		(-0.15,0.0),(0.05,0.1),(0.05,0.1),(0.05,0.1),(0.05,0.1),(0.05,0.1),
		(0.05,0.1),(0.05,0.0)]

## Read in intrinsic colors, and reddening law
stars = readcol.readcol("intrinsic_props.txt", asStruct=True, skipline=2)
redlaw = readcol.readcol("extlaw.txt", asStruct=True, skipline=1)

## set which color to use for X and Y axes
colX = stars.VR - stars.VI
colY = stars.VR
labelX = r"$R-I$"
labelY = r"$V-R$"
wlShort = 0.55
wlMid = 0.7
wlLong = 0.9

## Renormalize to A_V = 1.0
rv31 = redlaw.alam_aj_31 / redlaw.alam_aj_31[where(redlaw.wl == 0.55)]
rv50 = redlaw.alam_aj_50 / redlaw.alam_aj_50[where(redlaw.wl == 0.55)]

## Calculate color excesses
eColX31 = rv31[where(redlaw.wl == wlMid)][0] - rv31[where(redlaw.wl == wlLong)][0]
eColY31 = rv31[where(redlaw.wl == wlShort)][0] - rv31[where(redlaw.wl == wlMid)][0]
eColX50 = rv50[where(redlaw.wl == wlMid)][0] - rv50[where(redlaw.wl == wlLong)][0]
eColY50 = rv50[where(redlaw.wl == wlShort)][0] - rv50[where(redlaw.wl == wlMid)][0]

red31 = [(-1.5, 0.0), (-1.5 + eColX31, 0.0 + eColY31)]
red50 = [(-1.5, 0.0), (-1.5 + eColX50, 0.0 + eColY50)]

mask = (colX != -99) & (colY != -99)

line = plot(colX[mask], colY[mask], marker="o", linestyle="none", color=(0.5,0.5,1.0))
for x,y,sptype,offset in zip(colX[mask],colY[mask],stars.SpType[mask], labelOS):
	text(x+offset[0], y+offset[1], sptype)

#plot(*red31[0], color='green', marker='o')
#plot(*red50[1], color='green', marker='o')

annotate('', xy=red31[0], xytext=red31[1], 
		arrowprops=dict(arrowstyle="<|-", color='red'))
annotate('', xy=red50[0], xytext=red50[1], 
		arrowprops=dict(arrowstyle="<|-", color='red'))
text(red31[1][0]+0.05, red31[1][1]+0.05, '$R_V = 3.1$', va="bottom")
text(red50[1][0], red50[1][1]-0.05, '$R_V = 5.0$', va="bottom")

title("Main Sequence Visible-NIR Colors")
xlabel(labelX, fontsize="large")
ylabel(labelY, fontsize="large")

ylim(ylim()[::-1])

show()
