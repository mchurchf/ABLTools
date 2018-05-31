# exampleMMCSetUpSOWFA.py
#
#
# Matt Churchfield
# National Renewable Energy Laboratory
# 25 August 2017
#





from inflowPrepMMC import inflowPrepMMC
import numpy as np
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt




# Latitude and longitude of physics site.
latPhysicsSite = 45.638004
lonPhysicsSite = -120.642973
rotation = 0.0
rootDir = '../../fine/boundaryData/'
boundaries = ['north','south','east','west','lower','upper']






# Declare instance of class for inflow boundaries.
physicsSite = inflowPrepMMC();
westBoundary = inflowPrepMMC();
  

[xPhysicsSite,yPhysicsSite,UTMzonePhysicsSite] = physicsSite.LatLonToUTM_elem(latPhysicsSite,lonPhysicsSite)
[lat,lon] = physicsSite.UTMtoLatLon(xPhysicsSite,yPhysicsSite,UTMzonePhysicsSite)
print xPhysicsSite
print yPhysicsSite
print lat
print lon
    
                            

# Create a rotation matrix about z.
theta = rotation * np.pi / 180.0
R = np.zeros((3,3))
R[2,2] = 1.0
R[0,0] = np.cos(theta)
R[0,1] = -np.sin(theta)
R[1,0] = np.sin(theta)
R[1,1] = np.cos(theta)




for i in range(len(boundaries)):
    boundaryDir = rootDir + '/' + boundaries[i] + '/'
    boundary = inflowPrepMMC()
    
    boundary.readBoundaryDataPoints(boundaryDir + 'points')
    
    boundary.rotateXYZ(R)
    
    boundary.writeBoundaryDataPoints(boundaryDir + 'points.rotated')
    
    boundary.writeDataWRFtoFOAM(boundaryDir + boundaries[i] +'_bc.dat',xPhysicsSite,yPhysicsSite,UTMzonePhysicsSite)

    fig = plt.figure(10*(i+1))
    ax = fig.gca(projection='3d')
    ax.plot(boundary.xyz[:,0],boundary.xyz[:,1],boundary.xyz[:,2],'k.')
    plt.xlabel('x')
    plt.ylabel('y')