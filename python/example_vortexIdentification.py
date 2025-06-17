simport numpy as np
import readData
import vortexTools
import matplotlib.pyplot as plt



inputFile = '../../../array.mean0D_UAvg.vtk'








# Read in the structured VTK data file.
dataSetName,dims,origin,spacing,x,y,z,nFields,fieldName,fieldDim,field = readData.structuredVTK(inputFile)
ni = dims[0]
nj = dims[1]
nk = dims[2]
u = field[0]

#xx,yy = np.mgrid[x,y]
yy,zz = np.meshgrid(y,z)
vcross = np.sqrt(np.square(u[1,100,:,:]) + np.square(u[2,100,:,:]))

vx,vy,vz = vortexTools.findVortexCenterline(x,y,z,u)


plt.figure(1)
plt.pcolor(yy,zz,np.transpose(vcross))
plt.plot(vy,vz,'m.')