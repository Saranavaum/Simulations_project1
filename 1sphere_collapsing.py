# load libraries
import numpy as np  # load numpy
import h5py    # hdf5 format 
import random

#---- Initial condition parameters ----#

Mass = 1e11
G=4.3*10**(-6)                               #kpc*(km/s)²/M_sun
H_o=70*10**(-3)                              #(km/s²)/kpc
God_m=0.26                                   # Matter density of the universe
God_l=1-God_m                                # Dark energy
Z_init=15                                    # Initial redshift
Z_col=4                                      # Colapse redshift 
A_init=1./(1.+Z_init)                        # Initial expansion factor
A_col=1./(1.+Z_col)                          # Colapse expansion factor
Rho_crit=3.0*(H_o**2)/(8.0*np.pi*G)          # Average density of the Universe according to the Friedmann equations
H_init=H_o*((God_m/A_init**3)-God_l)**(1/2)  # Hubble constant to initial redshift
Rho_crit_init=3.0*(H_init**2)/(8.0*np.pi*G)  # Initial average density of the Universe according to the Friedmann equations
Lambda_init=1.686*A_init/A_col               # Overdensity needed to collapse
Rho_col=Rho_crit_init*(Lambda_init+1) 
R_col=(3.*Mass/(4.*np.pi*Rho_col))**(1/3)    # Collapse radio

numberOfPart = 10000
# numberOfPart = 1000  # According to the case in which we are
rMax = R_col
boxSize = rMax * 2                           # Box size containing the initial sphere of random particles
mpart=Mass/numberOfPart
x = np.zeros(numberOfPart)
y = np.zeros(numberOfPart)
z = np.zeros(numberOfPart)


# Then we generate random particles that are contained within a Rcol radio sphere
count = 0
i = 0
while count < numberOfPart:
    xTemp = random.random() * boxSize
    yTemp = random.random() * boxSize
    zTemp = random.random() * boxSize

    if (np.sqrt((xTemp - rMax)**2 + (yTemp - rMax)**2 + (zTemp - rMax)**2) < R_col):
        x[i] = xTemp
        y[i] = yTemp
        z[i] = zTemp
        count += 1
        i += 1

#fig = plt.figure()
#ax = plt.axes(projection='3d')
#ax.scatter(x,y,z, color = "g", s=0.8)
#plt.show()


#---- Initial velocities ----#

vx = np.zeros(numberOfPart)
vy = np.zeros(numberOfPart)
vz = np.zeros(numberOfPart)

#-- Generation of HDF5 with the initial conditions established --#

hf = h5py.File("myInitialConditions.hdf5", 'w')
mass=np.zeros(3)
mass[1]=mpart
number_part=np.zeros(3)
number_part[1]=numberOfPart

header = hf.create_group('Header')
header.attrs['NumPart_ThisFile']    = number_part
header.attrs['MassTable']           = mass
header.attrs['Time']                = 0
header.attrs['Redshift']            = 15
header.attrs['NumPart_Total']       = number_part
header.attrs['NumFilesPerSnapshot'] = 1
header.attrs['BoxSize']             = 1.0
header.attrs['Omega0']              = 1.0
header.attrs['OmegaLambdda']        = 0.
header.attrs['HubbleParam']         = 0.7
header.attrs['Flag_Entropy_ICs']    = 0
header.attrs['NumPart_Total_HighWord'] = np.zeros(3)

g2 = hf.create_group('Parameters')

g3 = hf.create_group('PartType1')
g3.create_dataset("Coordinates", data = np.vstack([x, y, z]).T)
g3.create_dataset("Velocities", data = np.vstack([vx, vy, vz]).T)
ids_d=np.arange(numberOfPart)
g3.create_dataset("ParticleIDs",(numberOfPart, ),dtype='f',data=ids_d)

hf.close()
