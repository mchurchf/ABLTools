# virtualRemoteSensing.py
#
# Matt Churchfield
# National Renewable Energy Laboratory
# 20 May 2017
#
# This is a class to deal with virtual remote sensing performed in CFD.




import numpy as np




class scanningLidarSOWFA:
    # initialize the class.
    def __init__(self,name):
        self.timeList = []
        self.beamVectorList = []
        self.beamDistribution = []
        self.maxBeamDist = 1.0
        self.timeBetweenScans = 0.0
        self.beamOrigin = [0.0,0.0,0.0]
        self.rotationAxis = [0.0,0.0,1.0]
        self.elevationAxis = [0.0,-1.0,0.0]
        self.beamAnglesVsTime = []
        self.Uname = "U"
        self.name = name
        
        
        

    # set other scan parameters.
    def setScanParameters(self,timeBetweenScans,beamOrigin,rotationAxis,elevationAxis,beamAnglesVsTime):
        self.timeBetweenScans = timeBetweenScans
        self.beamOrigin = beamOrigin
        self.rotationAxis = rotationAxis
        self.elevationAxis = elevationAxis
        self.beamAnglesVsTime = beamAnglesVsTime
        



    # set the scan pattern.
    def setPattern(self,azimuthMin,azimuthMax,dAzimuth,elevationMin,elevationMax,dElevation,dt):
        nAzimuth = int((azimuthMax - azimuthMin)/dAzimuth) + 1
        nElev = int((elevationMax - elevationMin)/dElevation) + 1
        
        self.timeList = np.zeros((nAzimuth*nElev,1))
        self.beamVectorList = np.zeros((nAzimuth*nElev,3))
        
        ii = 0
        t = 0.0
        
        for i in range(nElev):
            for j in range(nAzimuth):
                elevation = elevationMin+(i*dElevation)
                azimuth = azimuthMin+(j*dAzimuth)
                
                v = self.rotateVector([1.0,0.0,0.0],[0.0,1.0,0.0],-elevation)
                axis = self.rotateVector([0.0,0.0,1.0],[0.0,1.0,0.0],-elevation)
                v = self.rotateVector(v,axis,azimuth)
                
                self.beamVectorList[ii,0] = v[0]
                self.beamVectorList[ii,1] = v[1]
                self.beamVectorList[ii,2] = v[2]
                
                self.timeList[ii] = t
                t = t + dt
                
                ii = ii + 1
                
                
                
                
    # set the beam distribution.
    def setBeam(self,ds,maxDist):
        
        nPoint = int(maxDist/ds)
        
        self.beamDistribution = np.zeros((nPoint,1))
        self.maxBeamDist = maxDist
        
        for i in range(nPoint):
            self.beamDistribution[i] = (0.5*ds + i*ds)/maxDist
                                 
                                 
                                 
    
    # write the SOWFA input file for the sampling.
    def writeInputFile(self,fileName):
        f = open(fileName,'w')
        
        f.write('   ' + self.name + '\n')
        f.write('   {\n')
        
        f.write('      type                  scanningLidar;\n')
        
        f.write('   }\n')
        
        
        
        f.close()
                                 
                                 
    
    
    # rotate vector given angle and rotation axis.
    def rotateVector(self,v,axis,angle):
        RM = np.zeros((3,3))
        RM[0,0] = np.square(axis[0]) + (1.0 - np.square(axis[0])) * np.cos(np.deg2rad(angle));
        RM[0,1] = axis[0] * axis[1] * (1.0 - np.cos(np.deg2rad(angle))) - axis[2] * np.sin(np.deg2rad(angle));
        RM[0,2] = axis[0] * axis[2] * (1.0 - np.cos(np.deg2rad(angle))) + axis[1] * np.sin(np.deg2rad(angle));
        RM[1,0] = axis[0] * axis[1] * (1.0 - np.cos(np.deg2rad(angle))) + axis[2] * np.sin(np.deg2rad(angle));
        RM[1,1] = np.square(axis[1]) + (1.0 - np.square(axis[1])) * np.cos(np.deg2rad(angle));
        RM[1,2] = axis[1] * axis[2] * (1.0 - np.cos(angle)) - axis[0] * np.sin(angle);
        RM[2,0] = axis[0] * axis[2] * (1.0 - np.cos(np.deg2rad(angle))) - axis[1] * np.sin(np.deg2rad(angle));
        RM[2,1] = axis[1] * axis[2] * (1.0 - np.cos(np.deg2rad(angle))) + axis[0] * np.sin(np.deg2rad(angle));
        RM[2,2] = np.square(axis[2]) + (1.0 - np.square(axis[2])) * np.cos(np.deg2rad(angle));
        
        v = np.matmul(RM,v)
        return v