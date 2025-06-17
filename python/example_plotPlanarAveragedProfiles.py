# rotateDataPlanes.py
#
# Matt Churchfield
# National Renewable Energy Laboratory
# 24 October 2017
#
# This reads in plane averaged ABL simulation data and makes plots.



case = ['./zero',
        './elevated']

tAvgStart = 0.0
tAvgEnd = 1000.0

variables = ['U_mean','V_mean','T_mean']
varNames = ['U (m/s)','V (m/s)','\theta (K)']

lineType = ['k-','r-','b-','g-','c-','m-','y-','k--','b--','r--','g--','c--','m--','y--']
legendText = ['1 outer, 3 corrections', '2 outer, 1 correction']
lineWeight = [0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5]

modulePath = './ABLTools/python'




# Import the necessary modules/classes.
import readData
import writeData
from statistics import profileFit
import numpy as np
import matplotlib.pyplot as plt




# Load in data set 0 and average in time.
[z0,t0,dt0,U0] = readData.planarAverages(case[0]+'/averaging','U_mean')
[z0,t0,dt0,V0] = readData.planarAverages(case[0]+'/averaging','V_mean')
[z0,t0,dt0,T0] = readData.planarAverages(case[0]+'/averaging','T_mean')
[z0,t0,dt0,uu0] = readData.planarAverages(case[0]+'/averaging','uu_mean')
[z0,t0,dt0,vv0] = readData.planarAverages(case[0]+'/averaging','vv_mean')
[z0,t0,dt0,ww0] = readData.planarAverages(case[0]+'/averaging','ww_mean')

iStart = np.argmax(np.abs(t0 >= tAvgStart))
iEnd = np.argmax(np.abs(t0 >= tAvgEnd))
print('Average Start Time = ' + str(t0[iStart]))
print('Average End Time = ' + str(t0[iEnd]))

U0Avg = np.average(U0[iStart:iEnd,:],axis=0,weights=dt0[iStart:iEnd])
V0Avg = np.average(V0[iStart:iEnd,:],axis=0,weights=dt0[iStart:iEnd])
T0Avg = np.average(T0[iStart:iEnd,:],axis=0,weights=dt0[iStart:iEnd])
U0MagAvg = np.sqrt(U0Avg**2 + V0Avg**2)
uu0Avg = np.average(uu0[iStart:iEnd,:],axis=0,weights=dt0[iStart:iEnd])
vv0Avg = np.average(vv0[iStart:iEnd,:],axis=0,weights=dt0[iStart:iEnd])
ww0Avg = np.average(ww0[iStart:iEnd,:],axis=0,weights=dt0[iStart:iEnd])


plt.figure(1)
plt.plot(U0Avg,z0,'b-',lineWidth=lineWeight[0])
plt.plot(V0Avg,z0,'r-',lineWidth=lineWeight[0])
plt.plot(U0MagAvg,z0,'k-',lineWidth=lineWeight[0])
plt.legend(['U','V','Umag'])
plt.xlabel('wind speed (m/s)')
plt.ylabel('height (m)')

plt.figure(2)
plt.plot(T0Avg,z0,'k-',lineWidth=lineWeight[0])
plt.xlabel('potential temperature (K)')
plt.ylabel('height (m)')

plt.figure(3)
plt.plot(uu0Avg,z0,'b-',lineWidth=lineWeight[0])
plt.plot(vv0Avg,z0,'r-',lineWidth=lineWeight[0])
plt.plot(ww0Avg,z0,'k-',lineWidth=lineWeight[0])
plt.legend(['uu','vv','ww'])
plt.xlabel('velocity variance (m/s)^2')
plt.ylabel('height (m)')




# Load in data set 1 and average in time.
[z1,t1,dt1,U1] = readData.planarAverages(case[1]+'/averaging','U_mean')
[z1,t1,dt1,V1] = readData.planarAverages(case[1]+'/averaging','V_mean')
[z1,t1,dt1,T1] = readData.planarAverages(case[1]+'/averaging','T_mean')
[z1,t1,dt1,uu1] = readData.planarAverages(case[1]+'/averaging','uu_mean')
[z1,t1,dt1,vv1] = readData.planarAverages(case[1]+'/averaging','vv_mean')
[z1,t1,dt1,ww1] = readData.planarAverages(case[1]+'/averaging','ww_mean')

iStart = np.argmax(np.abs(t1 >= tAvgStart))
iEnd = np.argmax(np.abs(t1 >= tAvgEnd))
print('Average Start Time = ' + str(t1[iStart]))
print('Average End Time = ' + str(t1[iEnd]))

U1Avg = np.average(U1[iStart:iEnd,:],axis=0,weights=dt1[iStart:iEnd])
V1Avg = np.average(V1[iStart:iEnd,:],axis=0,weights=dt1[iStart:iEnd])
T1Avg = np.average(T1[iStart:iEnd,:],axis=0,weights=dt1[iStart:iEnd])
U1MagAvg = np.sqrt(U1Avg**2 + V1Avg**2)
uu1Avg = np.average(uu1[iStart:iEnd,:],axis=0,weights=dt1[iStart:iEnd])
vv1Avg = np.average(vv1[iStart:iEnd,:],axis=0,weights=dt1[iStart:iEnd])
ww1Avg = np.average(ww1[iStart:iEnd,:],axis=0,weights=dt1[iStart:iEnd])


plt.figure(1)
plt.plot(U1Avg,z1,'b--',lineWidth=lineWeight[1])
plt.plot(V1Avg,z1,'r--',lineWidth=lineWeight[1])
plt.plot(U1MagAvg,z1,'k--',lineWidth=lineWeight[1])
plt.legend(['U','V','Umag'])
plt.xlabel('wind speed (m/s)')
plt.ylabel('height (m)')

plt.figure(2)
plt.plot(T1Avg,z1,'k--',lineWidth=lineWeight[1])
plt.xlabel('potential temperature (K)')
plt.ylabel('height (m)')

plt.figure(3)
plt.plot(uu1Avg,z1,'b--',lineWidth=lineWeight[1])
plt.plot(vv1Avg,z1,'r--',lineWidth=lineWeight[1])
plt.plot(ww1Avg,z1,'k--',lineWidth=lineWeight[1])
plt.legend(['uu','vv','ww'])
plt.xlabel('velocity variance (m/s)^2')
plt.ylabel('height (m)')
