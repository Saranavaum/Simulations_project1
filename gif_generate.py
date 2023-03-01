# Packages
import h5py    
import numpy as np    
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import statistics
from scipy.optimize import curve_fit
from scipy.stats import gaussian_kde

# We read the HDF5 knowing that our simulations reach Snapshop27
for i in range(28):
	if i<10:
		filename="snapshot_00"+str(i)+".hdf5"
	if i>=10:
		filename="snapshot_0"+str(i)+".hdf5"

	f = h5py.File(filename, "r")
	group = f['PartType1']
	data = group["Coordinates"][()]

	x=data[:,0]
	y=data[:,1]
	z=data[:,2]
  
  # We make sure to center our data on 0.0
	x=x-np.median(x)
	y=y-np.median(y)
	z=z-np.median(z)
  
  # With this we get the density
	xyz = np.vstack([x, y, z])
	gauss = gaussian_kde(xyz)(xyz)
	
  # We take out the speeds from the HDF5
	data2 = group["Velocities"][()]
	vx=data2[:,0]
	vy=data2[:,1]
	vz=data2[:,2]
  
  # We get the distances from each particle with respect to the center
	dist=[]
	for k in range(len(x)):
		r2=np.sqrt(x[k]**2+y[k]**2+z[k]**2)
		dist.append(r2)
	dist=np.array(dist)


  # Plotting and saving the figures
	fig = plt.figure()
	ax = plt.axes(projection='3d')
	m=ax.scatter(x, y, z, c=vy, s=1)
  #m=ax. scatter (x, y, z, c=gauss , s =1) #If we are interested in representing the density
	ax.set_xlabel('X (kpc)')
	ax.set_ylabel('Y (kpc)')
	ax.set_zlabel('Z (kpc)')
	ax.set_xlim(-100,100)
	ax.set_ylim(-100,100)
	ax.set_zlim(-100,100)
	cbar=plt.colorbar(mappable=m,format='%.0e')
	cbar.ax.set_title('Vy (km/s)')
  # cbar .ax. set_title (’ Densidad ( Msolar /Kpc ) ’) #If we are interested in representing the density
	ax.set_title(filename[:-5] + '  Time='+str(i*0.5)+'  Gyrs')
	plt.show(block=False)

	plt.savefig("./Gif5/p310kGif" + str(i))
	plt.pause(1)
	plt.close()
