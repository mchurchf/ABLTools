# statistics.py
#
# Matt Churchfield
# National Renewable Energy Laboratory
# 8 May 2017
#
# This is a class that contains methods to compute basic time statistics, 
# spatial correlations, and spectra on a sampled flow field volume time
# history.



import numpy as np



class timeStatistics:
    # initialize the class.
    def __init__(self,field):
        dimsMean = np.asarray(field.shape)
        self.componentsMean = dimsMean[0]
        self.componentsVariance = self.componentsMean**2
        dimsVariance = 1*dimsMean
        dimsVariance[0] = self.componentsVariance
        
        self.meanField = np.zeros(dimsMean)
        self.varianceField = np.zeros(dimsVariance)
        self.meanAccumulatedTime = 0.0
        self.varianceAccumulatedTime = 0.0
        
        
        
    # accumulate the mean.
    def accumulateMean(self,fieldUpdate,dt):
        # take the current estimate of the mean, multiply by the previous 
        # total averaging time, add on the new contribution times its dt weight,
        # add the dt weight to the total averaging time, and divide by the new
        # total averaging time.
        self.meanField = self.meanField * self.meanAccumulatedTime
        self.meanField = self.meanField + (dt * fieldUpdate)
        self.meanAccumulatedTime = self.meanAccumulatedTime + dt
        self.meanField = self.meanField / self.meanAccumulatedTime
        
        
        
    # accumulate the variance tensor using the current estimate of the mean.
    def accumulateVariance(self,fieldUpdate,dt):
        # take the current estimate of the variances and multiply by the
        # previous total averaging time
        self.varianceField = self.varianceField * self.varianceAccumulatedTime
        
        # subtract out the current mean field from the new contribution field.
        fluctuationField = fieldUpdate - self.meanField
        
        # Add dt to the total averaging time.
        self.varianceAccumulatedTime = self.varianceAccumulatedTime + dt
        
        # Compute the new variance field contributions. Add the new variance 
        # contribution onto the current variance weighted by its dt, and divide
        # by the new total averaging time.
        for i in range(self.componentsMean):
            for j in range(self.componentsMean):
                v = np.multiply(fluctuationField[i,:],fluctuationField[j,:])
                index = (i*self.componentsMean) + j
                self.varianceField[index,:] = self.varianceField[index,:] + (dt * v)
                
        self.varianceField = self.varianceField / self.varianceAccumulatedTime
        
        
        
    # accumulate the variance tensor using a given mean.
    def accumulateVariance(self,fieldUpdate,fieldMean,dt):
        # take the current estimate of the variances and multiply by the
        # previous total averaging time
        self.varianceField = self.varianceField * self.varianceAccumulatedTime
        
        # subtract out the mean field from the new contribution field.
        fluctuationField = fieldUpdate - fieldMean
        
        # Add dt to the total averaging time.
        self.varianceAccumulatedTime = self.varianceAccumulatedTime + dt
        
        # Compute the new variance field contributions. Add the new variance 
        # contribution onto the current variance weighted by its dt, and divide
        # by the new total averaging time.
        for i in range(self.componentsMean):
            for j in range(self.componentsMean):
                v = np.multiply(fluctuationField[i,:],fluctuationField[j,:])
                index = (i*self.componentsMean) + j
                self.varianceField[index,:] = self.varianceField[index,:] + (dt * v)
                
        self.varianceField = self.varianceField / self.varianceAccumulatedTime
        
        
        
# class spatialCorrelations:
# I envision a class here to do spatial correlations.  Maybe you give it the 
# location which is the center point of the correlation.  Like, you give it
# the center of your data set, and then it loops through all points and does
# computes the correlation tensor with the center point.  In the end you get
# a tensor field (if you give it a vector) of correlation

            
# class spectra:
# I envision a class here to do compute the spectra.  You give it a time
# history of a quantity, and it returns the spectrum.  Maybe you are able to
# give it multiple time histories in a general array, and it does the spectra
# of each point, and averages them together.  Like, maybe you have the time
# history of velocity along the whole span of a sampling plane.  It would be 
# nice, too, if when you give it a vector time history, it returns a tensor of
# spectra.