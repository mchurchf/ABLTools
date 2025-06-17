cimport numpy as np
import readData
from statistics import profileFit
import matplotlib.pyplot as plt



avgStart = 20000.0
avgEnd = 22000.0
refHeight = 99.27
powerLawStartHeight = refHeight-77.0
powerLawEndHeight = refHeight+77.0
averagingDirectory = ['../averaging']
lineColors = ['k',\
              'b',\
              'r']
caseName = ['8.0 m/s',\
            '10.4 m/s',\
            '15.0 m/s']










z = []              # height above the flat surface (m)
t = []              # time (s)

u = []              # instantaneous planar-averaged x-component velocity (m/s)
v = []              # instantaneous planar-averaged y-component velocity (m/s)
w = []              # instantaneous planar-averaged y-component velocity (m/s)

uu = []             # instantaneous planar-averaged xx velocity variance correlation (m^2/s^2)
uv = []             # instantaneous planar-averaged xy velocity variance correlation (m^2/s^2)
uw = []             # instantaneous planar-averaged xz velocity variance correlation (m^2/s^2)
vv = []             # instantaneous planar-averaged yy velocity variance correlation (m^2/s^2)
vw = []             # instantaneous planar-averaged yz velocity variance correlation (m^2/s^2)
ww = []             # instantaneous planar-averaged zz velocity variance correlation (m^2/s^2)

uMean = []          # time/planar-averaged x-component velocity (m/s)
vMean = []          # time/planar-averaged x-component velocity (m/s)
wMean = []          # time/planar-averaged x-component velocity (m/s)
windSpeedMean = []  # time/planar-averaged horizontal magnitude of wind speed (m/s)
windDirMean = []    # time/planar-averaged direction of horizontal wind component (degrees)

uuMean = []         # time/planar-averaged xx velocity variance correlation (m^2/s^2)
uvMean = []         # time/planar-averaged xy velocity variance correlation (m^2/s^2)
uwMean = []         # time/planar-averaged xz velocity variance correlation (m^2/s^2)
vvMean = []         # time/planar-averaged yy velocity variance correlation (m^2/s^2)
vwMean = []         # time/planar-averaged yz velocity variance correlation (m^2/s^2)
wwMean = []         # time/planar-averaged zz velocity variance correlation (m^2/s^2)

uuRMean = []        # time/planar-averaged xx velocity variance correlation in local wind-aligned coordinates (m^2/s^2)
uvRMean = []        # time/planar-averaged xy velocity variance correlation in local wind-aligned coordinates (m^2/s^2)
uwRMean = []        # time/planar-averaged xz velocity variance correlation in local wind-aligned coordinates (m^2/s^2)
vvRMean = []        # time/planar-averaged yy velocity variance correlation in local wind-aligned coordinates (m^2/s^2)
vwRMean = []        # time/planar-averaged yz velocity variance correlation in local wind-aligned coordinates (m^2/s^2)
wwRMean = []        # time/planar-averaged zz velocity variance correlation in local wind-aligned coordinates (m^2/s^2)

TIMean = []         # time/planar-averaged turbulence intensity based on Cartesian coordinates (%)
TIRMean = []        # time/planar-averaged turbulence intensity in local wind-aligned coordinates (%)

windSpeedRefHeight = []
windDirRefHeight = []
TIRefHeight = []
TIRRefHeight = []
zPowerLaw = []
windSpeedPowerLaw = []
coeffs = []
fit = []



for m in range(len(averagingDirectory)):
    print 'Case ' + caseName[m] + ':'
    
    # Read in the velocity data.
    z_,t_,dt_,u_ = readData.planarAverages(averagingDirectory[m],'U_mean')
    z_,t_,dt_,v_ = readData.planarAverages(averagingDirectory[m],'V_mean')
    z_,t_,dt_,w_ = readData.planarAverages(averagingDirectory[m],'W_mean')
    z_,t_,dt_,uu_ = readData.planarAverages(averagingDirectory[m],'uu_mean')
    z_,t_,dt_,uv_ = readData.planarAverages(averagingDirectory[m],'uv_mean')
    z_,t_,dt_,uw_ = readData.planarAverages(averagingDirectory[m],'uw_mean')
    z_,t_,dt_,vv_ = readData.planarAverages(averagingDirectory[m],'vv_mean')
    z_,t_,dt_,vw_ = readData.planarAverages(averagingDirectory[m],'vw_mean')
    z_,t_,dt_,ww_ = readData.planarAverages(averagingDirectory[m],'ww_mean')
    
    
    z.append(z_)
    t.append(t_)
    u.append(u_)
    v.append(v_)
    w.append(w_)
    uu.append(uu_)
    uv.append(uv_)
    uw.append(uw_)
    vv.append(vv_)
    vw.append(vw_)
    ww.append(ww_)
    


    # Average the data in time.
    indexAvgStart = np.argmax(np.abs(t[m] >= avgStart))
    indexAvgEnd = np.argmax(np.abs(t[m] >= avgEnd))
    
    print '   Time avgeraging from index ' + str(indexAvgStart) + ' to ' + str((indexAvgEnd))
    
    uMean.append(np.mean(u[m][indexAvgStart:indexAvgEnd,:],axis=0))
    vMean.append(np.mean(v[m][indexAvgStart:indexAvgEnd,:],axis=0))
    wMean.append(np.mean(w[m][indexAvgStart:indexAvgEnd,:],axis=0))
    
    windSpeedMean.append(np.sqrt(np.square(uMean[m]) + np.square(vMean[m])))
    windDirMean.append(np.arctan2(vMean[m],uMean[m])*(180.0/np.pi))
    
    uuMean.append(np.mean(uu[m][indexAvgStart:indexAvgEnd,:],axis=0))
    uvMean.append(np.mean(uv[m][indexAvgStart:indexAvgEnd,:],axis=0))
    uwMean.append(np.mean(uw[m][indexAvgStart:indexAvgEnd,:],axis=0))
    vvMean.append(np.mean(vv[m][indexAvgStart:indexAvgEnd,:],axis=0))
    vwMean.append(np.mean(vw[m][indexAvgStart:indexAvgEnd,:],axis=0))
    wwMean.append(np.mean(ww[m][indexAvgStart:indexAvgEnd,:],axis=0))
    
    uuRMean.append(1.0*uuMean[m])
    uvRMean.append(1.0*uvMean[m])
    uwRMean.append(1.0*uwMean[m])
    vvRMean.append(1.0*vvMean[m])
    vwRMean.append(1.0*vwMean[m])
    wwRMean.append(1.0*wwMean[m])
    
    TIMean.append(1.0*uuMean[m])
    TIRMean.append(1.0*uuMean[m])
    
    
    
    # Get the correlation tensor in the wind-aligned coordinate system.
    for i in range(len(z[m])):
        T = np.zeros((3,3))
        T[0,0] =  np.cos(windDirMean[m][i]*(np.pi/180.0))
        T[0,1] =  np.sin(windDirMean[m][i]*(np.pi/180.0))
        T[0,2] =  0.0
        T[1,0] = -np.sin(windDirMean[m][i]*(np.pi/180.0))
        T[1,1] =  np.cos(windDirMean[m][i]*(np.pi/180.0))
        T[1,2] =  0.0
        T[2,0] =  0.0
        T[2,1] =  0.0
        T[2,2] =  1.0
        
        R = np.zeros((3,3))
        R[0,0] = uuMean[m][i]
        R[0,1] = uvMean[m][i]
        R[0,2] = uwMean[m][i]
        R[1,0] = uvMean[m][i]
        R[1,1] = vvMean[m][i]
        R[1,2] = vwMean[m][i]
        R[2,0] = uwMean[m][i]
        R[2,1] = vwMean[m][i]
        R[2,2] = wwMean[m][i]
        
        Rp = np.matmul(np.matmul(T,R),np.transpose(T))
        
        uuRMean[m][i] = Rp[0,0]
        uvRMean[m][i] = Rp[0,1]
        uwRMean[m][i] = Rp[0,2]
        vvRMean[m][i] = Rp[1,1]
        vwRMean[m][i] = Rp[1,2]
        wwRMean[m][i] = Rp[2,2]
        
        TIMean[m][i] = 100.0*np.sqrt(uuMean[m][i])/uMean[m][i]
        TIRMean[m][i] = 100.0*np.sqrt(uuRMean[m][i])/windSpeedMean[m][i]
    
    
    # Get the values at the reference height.
    index = np.argmax(np.abs(z[m] >= refHeight))
    
    wim1 = (z[m][index] - refHeight)/(z[m][index] - z[m][index-1])
    wi = (refHeight - z[m][index-1])/(z[m][index] - z[m][index-1])
    
    windSpeedRefHeight.append(wi*windSpeedMean[m][index] + wim1*windSpeedMean[m][index-1])
    windDirRefHeight.append(wi*windDirMean[m][index] + wim1*windDirMean[m][index-1])
    TIRefHeight.append(wi*TIMean[m][index] + wim1*TIMean[m][index-1])
    TIRRefHeight.append(wi*TIRMean[m][index] + wim1*TIRMean[m][index-1])
    
    
    print '   Reference-height mean wind speed ' + str(windSpeedRefHeight[m]) + ' m/s'
    print '   Reference-height mean wind direction ' + str(windDirRefHeight[m]) + ' deg'
    print '   Reference-height mean turbulence intensity ' + str(TIRRefHeight[m]) + '%'
    
    
    
    # Fit a power law to this profile.
    fit.append(profileFit(z[m],windSpeedMean[m]))
    fit[m].computePowerLawFit(refHeight,windSpeedRefHeight[m],powerLawStartHeight,powerLawEndHeight)
    coeffs.append(fit[m].exponentValue)
    windSpeedPowerLaw.append(fit[m].uFit)
    zPowerLaw.append(fit[m].zFit)
    print '   Power law fit exponent ' + str(coeffs[m])



    # Plot the time history of the velocity profiles at different heights.
    plt.figure(m+1)
    plt.plot(t[m],ux[m][:,0::20])
    plt.xlabel('time (s)')
    plt.ylabel('u (m/s)')
    
    
    
    
# Plot the mean profiles of wind speed.
for m in range(len(averagingDirectory)):
    plt.figure(10)
    plt.plot(windSpeedMean[m][:],z[m],str(lineColors[m] + '-'))
    plt.xlabel('wind speed (m/s)')
    plt.ylabel('z (m)')
    if (m == len(averagingDirectory)-1):
        plt.legend(caseName)

# Plot the power law fits of u.
for m in range(len(averagingDirectory)):
    plt.figure(10)
    plt.plot(windSpeedPowerLaw[m],zPowerLaw[m],str(lineColors[m] + '.'))
    
    
# Plot the mean profiles of wind directions.
for m in range(len(averagingDirectory)):
    plt.figure(20)
    plt.plot(windDirMean[m][:],z[m],str(lineColors[m] + '-'))
    plt.xlabel('wind direction (deg)')
    plt.ylabel('z (m)')
    if (m == len(averagingDirectory)-1):
        plt.legend(caseName)
        
        
    
    
# Plot the mean profiles of wind directions.
for m in range(len(averagingDirectory)):
    plt.figure(30)
    plt.plot(TIRMean[m][:],z[m],str(lineColors[m] + '-'))
    plt.xlabel('turbulence intensity (%)')
    plt.ylabel('z (m)')
    if (m == len(averagingDirectory)-1):
        plt.legend(caseName)
