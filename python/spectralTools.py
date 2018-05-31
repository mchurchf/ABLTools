#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May  9 15:40:15 2017

@author: pdoubraw
"""
import numpy as np
#==============================================================================
# 
#==============================================================================
def smoothSpectrum(freqs,psd,dl=0.1):    
    """
    Smooth a noisy spectrum.
    """
    nStart  = 0
    nEnd    = 1
    l       = 0
    freqs_smooth = []
    psd_smooth = []
    while nEnd<len(psd):
        freqs_smooth.append(np.mean(freqs[nStart:nEnd+1]))
        psd_smooth.append(np.mean(psd[nStart:nEnd+1]))
        nStart = nEnd ; l += dl ; m = int(np.round(2**l)) ; nEnd += m
    return np.array(freqs_smooth),np.array(psd_smooth)     