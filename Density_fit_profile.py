import numpy as np
from matplotlib import pyplot as plt
import numpy as np  # load numpy
import h5py    # hdf5 format 
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from scipy.optimize import curve_fit

# We are interested in the last Snapshot because it is where it is virialized and we can calculate the density profile well
# Snapshot reading
filename='snapshot_027.hdf5'

f = h5py.File(filename, "r")
group = f['PartType1']
data = group["Coordinates"][()]

x=data[:,0]
y=data[:,1]
z=data[:,2]

# Centering all coordinates
x=x-np.median(x)
y=y-np.median(y)
z=z-np.median(z)

#fig = plt.figure()
#ax = plt.axes(projection='3d')
#ax.scatter(x, y, z, color="g", s=0.1)
#plt.show()

# We calculate the distance of each particle with respect to the center
dist=[]
for i in range(len(x)):
	r2=np.sqrt(x[i]**2+y[i]**2+z[i]**2)
	dist.append(r2)
dist=np.array(dist)

# We generate shell spherical from the center to the last external particle and calculate the number of particles contained in that shell volume
shellNumber = 100                               # Number of shell we want
shellWidth = np.mean(dist) / shellNumber        # Size of each shell [Kpc]
shellDensity = np.array([])                     # Density [M_sol/kpc**3]
shellRadius = np.array([])                      # Radius [Kpc]
Mtot= 1e11                                      # Total mass [M_sol]
Mpart=Mtot/len(dist)                            # Particle mass

for i in range(shellNumber):
	numberOfParticles = 0
	for j in range(len(dist)):
		if ((dist[j] > (i * shellWidth)) and (dist[j] < ((i+1) * shellWidth))):
			numberOfParticles = numberOfParticles + 1
	shellRadius = np.append(shellRadius, shellWidth*i+(shellWidth)/2)
	massOfParticles = numberOfParticles * Mpart
	shellDensity = np.append(shellDensity, (massOfParticles)/((4./3)*np.pi*(((i+1)*shellWidth)**3 - (i*shellWidth)**3)))


# Fitting NFW
def NFW(r,rho_0,R_s):
	rho=rho_0/(r/R_s*(1+r/R_s)**2)
	return(rho)


sigma1=1/(shellRadius) # We weigh with more error the values furthest from the center
param1,cov1=curve_fit(NFW, shellRadius, shellDensity, p0=[5,30],sigma=sigma1, maxfev=50000)
fitNFW=NFW(shellRadius,param1[0],param1[1])
print(param1[0],param1[1])
print(np.sqrt(np.diag(cov1))) # errors


# Fitting general double-power law model
def GDPL(r,rho_s,R_s,alpha,beta,gamma):
	rho=rho_s/( (r/R_s)**gamma*((1+r/R_s)**alpha)**((beta-gamma)/alpha) )
	return(rho)



param,cov=curve_fit(GDPL, shellRadius, shellDensity, p0=[100000,30,1,3,1],sigma=sigma1,maxfev=50000,bounds=[0,[1e10,100,4,5,4]])

print(param[0],param[1],param[2],param[3],param[4])
print(np.sqrt(np.diag(cov))) #errors of each parameter


# Plotting the profiles
textstr = '\n'.join((
	r'   NFW profile',
	r'--------------------',
	r'$\rho=%.2e$' % (param1[0], ),
	r' $[\mathrm{M}_{\odot}/\mathrm{kpc}^3]$',
	r'$r_s=%.2f$ [kpc]' % (param1[1], ),
))
te = '\n'.join((
	r'     GDPL profile',
	r'-------------------------',
	r'$\rho=%.2e$' % (param[0], ),
	r' $[\mathrm{M}_{\odot}/\mathrm{kpc}^3]$',
	r'$r_s=%.2f$ [kpc]' % (param[1], ),
	r'$\alpha = %.2f $'%(param[2], ),
	r'$\beta = %.2f $'%(param[3], ),
	r'$\gamma = %.2f $'%(param[4], ),
))
fig=plt.figure()
props = dict(boxstyle='round', facecolor='white', alpha=0.5)
plt.text(0.15, 0.18, textstr, fontsize=10,transform=plt.gcf().transFigure, bbox = props)
plt.text(0.34, 0.18, te, fontsize=10,transform=plt.gcf().transFigure, bbox = props)


fitGDPL=GDPL(shellRadius,param[0],param[1],param[2],param[3],param[4])
plt.gca().set_xscale('log')
plt.gca().set_yscale('log')
plt.plot(shellRadius,shellDensity,'r.',label='Datos')
plt.plot(shellRadius,fitNFW,label='Ajuste NFW')
plt.plot(shellRadius,fitGDPL,label='Ajuste GDPL')
plt.title('Density Profile')
plt.xlabel('log (r) [Kpc]')
plt.ylabel(r'log($\rho$) [$\mathrm{M}_{\odot}$/Kpc³]')


plt.grid()
plt.legend()
#plt.savefig('Perfil5.png')
plt.show()

