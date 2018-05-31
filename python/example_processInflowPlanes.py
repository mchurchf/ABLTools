# driverExample.py
#
# Matt Churchfield
# National Renewable Energy Laboratory
# 8 May 2017
#
# This is an example python driver script that uses the modules contained in ABLTools 
# to postprocess a set of structured VTK data, each sample being contained in a separate
# file in a separe directory.



# Import the necessary modules/classes.
import getDataLayout
import readData
from statistics import timeStatistics
from statistics import spatialCorrelations

import matplotlib
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

import numpy as np





# Specify the case to analyze.
case = '../../run.neutral.8/array.1'

# Mean start time.
meanTimeStart = 10400.0

# Variance start time.
varianceTimeStart = 11400.0

# Spatial correlation start time.
spatialCorrelationTimeStart = 12505.0

# For plotting.
nLevels = 51;
matplotlib.rcParams['font.serif'] = "Courier New"
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.rm'] = 'serif'
matplotlib.rcParams['text.usetex'] = 'false'
matplotlib.rcParams['font.size'] = 7
matplotlib.rcParams['figure.figsize'] = [10.0,6.0]
default_cmap = plt.cm.viridis





















# Get the data directory layout.
nTimes, outputTimes = getDataLayout.getOutputTimes(case)


# Read in the first sample of velocity so that we can set up the basic time statistics.
fileName = case + '/' + outputTimes[0] + '/' + 'array.1_U.vtk'
dataSetName, dims, origin, spacing, x, y, z, nFields, fieldName, fieldDim, field = readData.structuredVTK(fileName)


# Reshape the velocity data into a plane (i.e. remove the x-axis because data
# here is in y and z.
#u = np.array(field[0][:,0,:,:,].reshape((fieldDim[0],dims[2],dims[1])))
u = np.array(field[0])


# Create the time statistics class object.
uTimeStatistics = timeStatistics(u)
#uTimeStatistics = timeStatistics(field[0])


# Create the spatial correlations statistics class object.
uSpatialCorrelations = spatialCorrelations(u,x,y,z,1500.0,1000.0,99.72)


# Loop over time and accumuate statistics, but disregard anything before the time at which
# you want to start taking the mean.
for m in range(nTimes):
    if (float(outputTimes[m]) >= meanTimeStart):
        # Tell where we are.
        print('Processing ' + outputTimes[m] + '...')
    
    
        # Get the time step size.
        if (m > 0):
             dt = float(outputTimes[m]) - float(outputTimes[m-1])
        else:
            dt = float(outputTimes[m+1]) - float(outputTimes[m])
        
        
        # Read in the data.    
        fileName = case + '/' + outputTimes[m] + '/' + 'array.1_U.vtk'
        dataSetName, dims, origin, spacing, x, y, z, nFields, fieldName, fieldDim, field = readData.structuredVTK(fileName)
    
    
        # Reshape the data into a plane (i.e. remove the x-axis because data here
        # is in y and z.
        #u = np.array(field[0][:,0,:,:,].reshape((fieldDim[0],dims[2],dims[1])))
        u = np.array(field[0])
    
    
        # Acculate the mean.
        uTimeStatistics.accumulateMean(u,dt)
    


        # Accumulate the variance.
        if (float(outputTimes[m]) >= varianceTimeStart):
            uTimeStatistics.accumulateVariance(u,uTimeStatistics.meanField,dt)
            
        
        # Accumulate the spatial correlations.
        if (float(outputTimes[m]) >= spatialCorrelationTimeStart):
            uSpatialCorrelations.accumulateCorrelation(u,uTimeStatistics.meanField,dt)
            

print 'Computed statistics with:'
print '   ' + str(uTimeStatistics.meanAccumulatedTime) + 's of mean accumulation'
print '   ' + str(uTimeStatistics.varianceAccumulatedTime) + 's of variance accumulation'
print '   ' + str(uSpatialCorrelations.correlationAccumulatedTime) + 's of correlation accumulation'



# Plot an instantaneous slice of velocity
#u = u[0,:,:];
u = u[0,0,:,:].reshape((dims[2],dims[1]))
plt.figure(1)
plt.contourf(y,z,u,nLevels,cmap=default_cmap)
plt.gca().set_aspect('equal', adjustable='box')
cbar = plt.colorbar(shrink=0.6, extend='both')
cbar.ax.set_ylabel('$u_x$' + ' (m/s)',fontweight='bold')
plt.xlabel('y (m)',fontweight='bold')
plt.ylabel('z (m)',fontweight='bold')
plt.tight_layout()


# Plot the mean velocity.
uAvg = uTimeStatistics.meanField[0,:].reshape((dims[2],dims[1]));
plt.figure(2)
plt.contourf(y,z,uAvg,nLevels,cmap=default_cmap)
plt.gca().set_aspect('equal', adjustable='box')
cbar = plt.colorbar(shrink=0.6, extend='both')
cbar.ax.set_ylabel('$\overline{u}_x$' + ' (m/s)',fontweight='bold')
plt.xlabel('y (m)',fontweight='bold')
plt.ylabel('z (m)',fontweight='bold')
plt.tight_layout()



# Plot the variance field.  The components map like this;
# 0 = xx, 1 = xy, 2 = xz
# 3 = yx, 4 = yy, 5 = yz
# 6 = zx, 7 = zy, 8 = zz
cmpt = 2
var = uTimeStatistics.varianceField[cmpt,:].reshape((dims[2],dims[1]));
plt.figure(3)
plt.contourf(y,z,var,nLevels,cmap=default_cmap)
plt.gca().set_aspect('equal', adjustable='box')
cbar = plt.colorbar(shrink=0.6, extend='both')
cbar.ax.set_ylabel('$\overline{u\'_x u\'_y}$' + ' (m/s)',fontweight='bold')
plt.xlabel('y (m)',fontweight='bold')
plt.ylabel('z (m)',fontweight='bold')
plt.tight_layout()


cmpt = 0
var = uTimeStatistics.varianceField[cmpt,:].reshape((dims[2],dims[1]));
plt.figure(4)
plt.contourf(y,z,var,nLevels,cmap=default_cmap)
plt.gca().set_aspect('equal', adjustable='box')
cbar = plt.colorbar(shrink=0.6, extend='both')
cbar.ax.set_ylabel('$\overline{u\'_x u\'_x}$' + ' (m/s)',fontweight='bold')
plt.xlabel('y (m)',fontweight='bold')
plt.ylabel('z (m)',fontweight='bold')
plt.tight_layout()



# Plot the spatial correlation field.  The components map like this;
# 0 = xx, 1 = xy, 2 = xz
# 3 = yx, 4 = yy, 5 = yz
# 6 = zx, 7 = zy, 8 = zz
cmpt = 2
var = uSpatialCorrelations.correlationField[cmpt,:].reshape((dims[2],dims[1]))
var = var / uSpatialCorrelations.autoCorrelation[cmpt]
plt.figure(5)
plt.contourf(y,z,var,nLevels,cmap=default_cmap)
plt.gca().set_aspect('equal', adjustable='box')
cbar = plt.colorbar(shrink=0.6, extend='both')
cbar.ax.set_ylabel('$C_{xz}$' + ' (m/s)',fontweight='bold')
plt.xlabel('y (m)',fontweight='bold')
plt.ylabel('z (m)',fontweight='bold')
plt.tight_layout()


cmpt = 0
var = uSpatialCorrelations.correlationField[cmpt,:].reshape((dims[2],dims[1]))
var = var / uSpatialCorrelations.autoCorrelation[cmpt]
plt.figure(6)
plt.contourf(y,z,var,nLevels,cmap=default_cmap)
plt.gca().set_aspect('equal', adjustable='box')
cbar = plt.colorbar(shrink=0.6, extend='both')
cbar.ax.set_ylabel('$C_{xx}$' + ' (m/s)',fontweight='bold')
plt.xlabel('y (m)',fontweight='bold')
plt.ylabel('z (m)',fontweight='bold')
plt.tight_layout()