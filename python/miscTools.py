#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May  9 15:39:03 2017

@author: pdoubraw
"""
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
#==============================================================================
# 
#==============================================================================
def parse(string):
    """
    Convert strings to float.
    """
    try:
        return float(string)
    except:
        return np.nan
#==============================================================================
# 
#==============================================================================
def matlab_datenum_to_python_datetime(datenum):
    """
    Parameters
    ----------
    datenum : int
    """
    return datetime.fromordinal(int(datenum)) +\
             timedelta(days=datenum%1) - timedelta(days = 366)
#==============================================================================
#              
#==============================================================================
def make_multi_index(tuple_1=None,tuple_2=None,tuple_3=None):
    """
    Returns a multi_index to be used in a data frame. Super useful!
    
    Parameters
    ----------
    
    tuple_1 : tuple or list,
        values for the leftmost index
    tuple_2 : tuple or list,
        values for the next index to the right
    tuple_3 : tuple or list,
        values for the next index to the right
    
    Returns
    -------    
    multi_index : a pandas index instance
    
    Examples
    --------
    multi_index = make_multi_index(std_roughness,std_heights,['A','k'])    
    """    
    if (tuple_1 is None) or (tuple_2 is None): 
        return
        
    t1 = []
    t2 = []
    t3 = []
    for x1 in tuple_1:
        for x2 in tuple_2:
            if tuple_3 is not None:
                for x3 in tuple_3:
                    t1.append(x1)
                    t2.append(x2)
                    t3.append(x3)
            else:
                t1.append(x1)
                t2.append(x2)
    if tuple_3 is not None:
        arrays = [t1,t2,t3]
    else:
        arrays = [t1,t2]        
    tuples = list(zip(*arrays))
    
    return pd.MultiIndex.from_tuples(tuples) 