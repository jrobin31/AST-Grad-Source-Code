from pylab import *
s50myr = genfromtxt('./50Myr.spectrum1', names=True,dtype=None)
s1Gyr3600 = genfromtxt('./1Gyr3600.spectrum1', names=True,dtype=None)
s1Gyr9000 = genfromtxt('./1Gyr9000.spectrum1', names=True,dtype=None)
line50 = plot(s50myr['WAVELENGTH'], s50myr['LOGTOTAL'], color='blue', linestyle='solid',
	label='$\mathsf{50\/Myr\/1.00\\times10^6 M_{\odot}}$')
line1k36 = plot(s1Gyr3600['WAVELENGTH'], s1Gyr3600['LOGTOTAL'], color='red', 
	linestyle='solid', label='$\mathsf{1\/Gyr\/24.0\\times10^6 M_{\odot}}$')
line1k90 = plot(s1Gyr9000['WAVELENGTH'], s1Gyr9000['LOGTOTAL'], color='green', 
	linestyle='solid',label='$\mathsf{1\/Gyr\/5.35\\times10^6 M_{\odot}}$')
axis([3000,10000,36,37.5])
plot([3600,3600],[35,38],color='black',linestyle='dashed')
plot([9000,9000],[35,38],color='black',linestyle='dashed')
title('Starburst 99 Single-Burst Models')
legend(fancybox=True, shadow=True)
xlabel('Wavelength ($\AA$)')
ylabel('log(Flux) (erg/s/$\AA$)')
show()
