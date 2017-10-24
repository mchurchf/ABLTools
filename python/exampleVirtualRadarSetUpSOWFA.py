# exampleVirtualRadarSetUpSOWFA.py
#
#
# Matt Churchfield
# National Renewable Energy Laboratory
# 20 May 2017
#
# radarName                  name of this radar
# inputFileName              file name for SOWFA sampling input
# azimuthMin                 miniumum azimuth angle of radar (deg)
# azimuthMax                 maximum azimuth angle of radar (deg)
# elevationMin               minimum elevation angle of radar (deg)
# elevationMax               maximum elevation angle of radar (deg)
# deltaAzimuth               azimuth angle increment between each sample beam (deg)
# deltaElevation             elevation angle increment between each sample plane (deg)
# timeBetweenEachBeam        time taken to sample a single beam (s)
# timeBetweenFullScan        time between end of one scan and start of next scan (s)
# beamResolution             beam resolution along beam path (m)
# beamDistanceMin            minimum scanning distance (m)
# beamDistanceMax            maximum scanning distance (m)
# rotationAxis               axis of azimuth rotation (unit vector) before global rotation.
# elevationAxis              axis of elevation rotation (unit vector) before global rotation.
# beamOrigin                 origin of beam (m, m, m)
# tiltOption                 1: change of elevation tilts entire azimuth rotation axis (scan surfaces are always planes)
#                            2: change of elevation does not tilt azimuth rotation axis (scan surfaces are conical for non-zero elevations)
# beamAngleVsTime            [time0, azimuth0, elevation0]
#                             ...
#                            [timeN, azimuthN, elevationN]
#                            the time varying global rotation of the entire radar (in addition to scanning rotation)





from virtualRemoteSensing import scanningLidarSOWFA
import numpy as np




# Radar 1:
radarName = "radar.1"
inputFileName = "radar.1"
azimuthMin = -35.0
azimuthMax =  35.0
elevationMin = 0.3
elevationMax = 2.9
deltaAzimuth = 0.25
deltaElevation = 0.18571428571428572
timeBetweenEachBeam = 0.01652262328418912
timeBetweenFullScan = 3.0
beamResolution = 5.0
beamDistanceMin = 1000.0
beamDistanceMax = 3000.0
rotationAxis = [0.0,0.0,1.0]
elevationAxis = [0.0,-1.0,0.0]
tiltOption = 2;
beamOrigin = [3750.0,1500.0,0.0]
beamAngleVsTime = [[0.0,180.0,0.0],[90000.0,180.0,0.0]]


radar1 = scanningLidarSOWFA(radarName)
radar1.setPattern(azimuthMin,azimuthMax,deltaAzimuth,elevationMin,elevationMax,deltaElevation,timeBetweenEachBeam,tiltOption)
radar1.setBeam(beamResolution,beamDistanceMin,beamDistanceMax)
radar1.setScanParameters(timeBetweenFullScan,beamOrigin,rotationAxis,elevationAxis,beamAngleVsTime)
radar1.writeInputFile(inputFileName)




# Radar 2:
radarName = "radar.2"
inputFileName = "radar.2"
azimuthMin = -35.0
azimuthMax =  35.0
elevationMin = 0.3
elevationMax = 2.9
deltaAzimuth = 0.25
deltaElevation = 0.18571428571428572
timeBetweenEachBeam = 0.01652262328418912
timeBetweenFullScan = 3.0
beamResolution = 5.0
beamDistanceMin = 1000.0
beamDistanceMax = 3000.0
rotationAxis = [0.0,0.0,1.0]
elevationAxis = [0.0,-1.0,0.0]
tiltOption = 2;
beamOrigin = [3750.0,1500.0,0.0]
beamAngleVsTime = [[0.0,0.0,0.0],[90000.0,0.0,0.0]]


radar2 = scanningLidarSOWFA(radarName)
radar2.setPattern(azimuthMin,azimuthMax,deltaAzimuth,elevationMin,elevationMax,deltaElevation,timeBetweenEachBeam,tiltOption)
radar2.setBeam(beamResolution,beamDistanceMin,beamDistanceMax)
radar2.setScanParameters(timeBetweenFullScan,beamOrigin,rotationAxis,elevationAxis,beamAngleVsTime)
radar2.writeInputFile(inputFileName)




# Lidar 1:
lidarName = "StuttgartLidar"
inputFileName = "StuttgartLidar"
azimuthMin = -15.0
azimuthMax =  15.0
elevationMin = -10.0
elevationMax =  10.0
deltaAzimuth = 5.0
deltaElevation = 3.33333333333333
timeBetweenEachBeam = 1.0
timeBetweenFullScan = 1.0
beamResolution = 1.25
beamDistanceMin = 0.0
beamDistanceMax = 5.0*77.0
rotationAxis = [0.0,0.0,1.0]
elevationAxis = [0.0,-1.0,0.0]
tiltOption = 2;
beamOrigin = [1500.0,2500.0,80.0]
beamAngleVsTime = [[0.0,25.0,0.0],[90000.0,25.0,0.0]]


lidar1 = scanningLidarSOWFA(lidarName)
lidar1.setPattern(azimuthMin,azimuthMax,deltaAzimuth,elevationMin,elevationMax,deltaElevation,timeBetweenEachBeam,tiltOption)
lidar1.setBeam(beamResolution,beamDistanceMin,beamDistanceMax)
lidar1.setScanParameters(timeBetweenFullScan,beamOrigin,rotationAxis,elevationAxis,beamAngleVsTime)
lidar1.writeInputFile(inputFileName)