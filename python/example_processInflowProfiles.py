import numpy as np
import readData
from statistics import profileFit
import matplotlib.pyplot as plt



avgStart = 16000.0
avgEnd = 18000.0
hubHeight = 99.27
powerLawStartHeight = 5.0
powerLawEndHeight = 500.0
averagingDirectory = ['../postProcessing/averaging']
lineColors = ['k',\
              'b',\
              'r']
caseName = ['8.0 m/s',\
            '10.4 m/s',\
            '15.0 m/s']
varName = 'U'





z = []
t = []
ux = []
uy = []
uMean = []
uMeanHubHeight = []
zPowerLaw = []
uPowerLaw = []
coeffs = []
fit = []

for m in range(len(averagingDirectory)):
    print 'Case ' + caseName[m] + ':'
    # Read the data.
    z_,t_,dt_,ux_ = readData.planarAverages(averagingDirectory[m],'U_mean')
    z_,t_,dt_,uy_ = readData.planarAverages(averagingDirectory[m],'V_mean')
    z.append(z_)
    t.append(t_)
    ux.append(ux_)
    uy.append(uy_)
    

    # Average the data in time.
    indexAvgStart = np.argmax(np.abs(t[m] >= avgStart))
    indexAvgEnd = np.argmax(np.abs(t[m] >= avgEnd))
    print '   Time avgeraging from index ' + str(indexAvgStart) + ' to ' + str((indexAvgEnd))
    uMean.append(np.sqrt(np.square(np.mean(ux[m][indexAvgStart:indexAvgEnd,:],axis=0)) + 
                         np.square(np.mean(uy[m][indexAvgStart:indexAvgEnd,:],axis=0))))
    
    
    # Get the value at hub height.
    index = np.argmax(np.abs(z[m] >= hubHeight))
    uMeanHubHeight.append(uMean[m][index])
    print '   Hub-height mean wind speed ' + str(uMeanHubHeight[m])
    
    
    # Fit a power law to this profile.
    fit.append(profileFit(z[m],uMean[m]))
    fit[m].computePowerLawFit(hubHeight,uMeanHubHeight[m],powerLawStartHeight,powerLawEndHeight)
    coeffs.append(fit[m].exponentValue)
    uPowerLaw.append(fit[m].uFit)
    zPowerLaw.append(fit[m].zFit)
    #coeffs.append(np.polyfit(np.log(z[m][indexPowerLawStart[m]:indexPowerLawEnd[m]]), np.log(uMean[m][indexPowerLawStart[m]:indexPowerLawEnd[m]]),1))
    #uPowerLaw.append(uMeanHubHeight[m] * np.power((z[m][indexPowerLawStart[m]:indexPowerLawEnd[m]]/hubHeight),coeffs[m][0]))
    print '   Power law fit exponent ' + str(coeffs[m])


    # Plot the time history of the velocity profiles at different heights.
    plt.figure(m+1)
    plt.plot(t[m],ux[m][:,0::20])
    plt.xlabel('time (s)')
    plt.ylabel('u (m/s)')
    
    
    
    
# Plot the mean profiles of u.
for m in range(len(averagingDirectory)):
    plt.figure(10)
    plt.plot(uMean[m][:],z[m],str(lineColors[m] + '-'))
    plt.xlabel('<u> (m/s)')
    plt.ylabel('z (m)')
    if (m == len(averagingDirectory)-1):
        plt.legend(caseName)


# Plot the power law fits of u.
for m in range(len(averagingDirectory)):
    plt.figure(10)
    plt.plot(uPowerLaw[m],zPowerLaw[m],str(lineColors[m] + '--'))