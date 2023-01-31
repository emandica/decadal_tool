#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 22 17:05:41 2023

@author: ema
"""
import numpy as np

import xarray as xr

import constants as c

def block_bootstrap(dataset, k):
    """
    dataset (member,startdate)
    k: resample size
    """
    ancillary_a = []
    for i in range(k):
#        seleziona a caso con sostituzione n start dates
        r_time = np.random.choice(dataset.time.size, size=round(len(dataset.time)/5)+1, replace=True)
        
#        campionamento di 5 anni successivi con considerare l'autocorrelazione
        ref_time=[]
        for i in r_time:
            if i <= dataset.time.size-5:
                ref_time.append([i,i+1,i+2,i+3,i+4])
            else:
                j=dataset.time.size
                ref_time.append([j-1, j-2, j-3, j-4, j-5])
        
        new_time = [j for sub in ref_time for j in sub][:len(dataset.time)]
        
        a = dataset.isel(time = new_time)
        a.coords['time'] = dataset['time']
        
#   seleziona a caso N membri del dataset
        a = a.isel(member = np.random.choice(a.member.size, size=len(dataset.member), replace=True))
        
        ancillary_a.append(a.mean('member'))
        
    ancillary_a = xr.concat(ancillary_a, dim='bootstrap')
    
    return ancillary_a

###############################################################################
def standard_bootstrap(dset, k):
    """"""
    ancillary_a=[]
    for i in range(k):        
        a = dset.isel(member = np.random.choice(dset.member.size, size=c.CTL_NUMBER, replace=True))
        ancillary_a.append(a.mean('member'))
    ancillary_a = xr.concat(ancillary_a, dim='bootstrap')   
    return ancillary_a     
    