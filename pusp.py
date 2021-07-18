import matplotlib.pyplot as plt
import numpy as np
import numpy.matlib
import math

lmin = 2.0          # source min wavelength
lmax = 5.0          # source max wavelength
fmin = 1/lmax       # source min frequency
fmax = 1/lmin       # source max frequency
fcen = 0.5*(fmin+fmax)

thetas = range(0,35,5)
kx = [fcen*math.sin(t) for t in [math.radians(float(t)) for t in thetas]]
Refl = np.empty((50,thetas.size))
Abs = np.empty((50,thetas.size))
theta_out = np.empty((50,thetas.size))

for k in range(thetas.size):
    f0 = np.genfromtxt("flux0_a4.3_theta{}.dat".format(thetas[k]), delimiter=",")
    f = np.genfromtxt("flux_a4.3_r1.72_theta{}.dat".format(thetas[k]), delimiter=",")
    Refl[:,k] = -f[:,1]/f0[:,1]
    theta_out[:,k] = np.asarray([math.degrees(math.asin(kx[k]/f0[j,0])) for j in range(50)])
    Abs[:,k] = np.asarray([(1-Refl[j,k])*math.cos(math.radians(theta_out[j,k])) for j in range(50)])

Abs[Abs<0] = 0
wvl = 1/f0[:,0]
wvls = np.matlib.repmat(np.reshape(wvl,(wvl.size,1)),1,thetas.size)

plt.figure(dpi=100)
plt.pcolormesh(theta_out, wvls, Abs, cmap='hot_r', shading='gouraud', vmin=0, vmax=Abs.max())
plt.axis([theta_out.min(), theta_out.max(), wvl[-1], wvl[0]])
plt.xlabel("emission angle (degrees)")
plt.xticks([t for t in range(0,60,10)])
plt.ylabel("wavelength (um)")
plt.yticks([t for t in np.arange(2.0,5.5,0.5)])
plt.title(r"emissivity: a=4.3 um, r=1.72 um")
cbar = plt.colorbar()
cbar.set_ticks([t for t in np.arange(0,1.0,0.2)])
cbar.set_ticklabels(["{:.1f}".format(t) for t in np.arange(0,1.0,0.2)])
plt.show()
