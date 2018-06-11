# viewBoundaryData.py
#
#
# Matt Churchfield
# National Renewable Energy Laboratory
# 7 June 2018
#



from inflowPrepMMC import inflowPrepMMC
import numpy as np
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt




rootDir = '../boundaryData.old'
boundaries = ['west','lower','east','upper']
fieldName = 'U'
pointsFileName = 'points.30'
time = 900







#for i in range(len(boundaries)):
for i in range(1):
    boundaryDir = rootDir + '/' + boundaries[i] + '/'
    
    boundary = inflowPrepMMC()
    
    boundary.readBoundaryDataPoints(boundaryDir + pointsFileName)
    
    f = boundary.readBoundaryDataField(boundaryDir + '/' + str(time) + '/' + fieldName)

    plt.figure(2)
    #plt.scatter(boundary.xyz[:,1], boundary.xyz[:,2], boundary.xyz[:,2], cmap='nipy_spectral')
    plt.scatter(boundary.xyz[:,1], boundary.xyz[:,2], s=0.75, c=f[:,0], cmap='nipy_spectral')
    plt.colorbar()
    plt.axis('equal')
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.show()