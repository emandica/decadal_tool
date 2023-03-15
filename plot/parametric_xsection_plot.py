#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 18:21:18 2023

@author: ema
"""

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


import constants as c

import generic_routines as cr

#%%
def parametric_xsection_plot(ds, ds_sign, levels, title=None, sign=0.95):
    
    #plot significance
    if sign == 0.95:
        min_v = 0.025
        max_v = 0.975
    elif sign == 0.90:
        min_v = 0.05
        max_v = 0.95
    elif sign == 0.80:
        min_v = 0.10
        max_v = 0.90
        
    ds = cr.convert_longitudes(ds, 'lon')
    ds_sign = cr.convert_longitudes(ds_sign, 'lon')
    
    #ds = ds.where((ds_sign <= 0.1) | (ds_sign >= 0.9))
    sig = ds.where((ds_sign <= min_v)|(ds_sign >= max_v))
    
#%%plot fields
    plt.rcParams['font.size'] = 18
    fig = plt.figure(figsize=[12,8])
    ax = fig.add_subplot(111)
    
    p = ds[c.VAR].plot(ax=ax, levels=levels,
                   add_labels=False, add_colorbar=False,
                   extend='both',
                   )
    
    _ = sig[c.VAR].plot.contourf(ax=ax, hatches='/',
                                 add_labels=False, add_colorbar=False,
                                 alpha=0
                                 )
    
    # add separate colorbar
    cb = plt.colorbar(p, ticks=levels, shrink=1, extend='both')
    cb.set_label(c.UNITS,fontsize=18)
    cb.ax.tick_params(labelsize=18)

    #Drow gridlines and adjust labels
    ax.grid()
    ax.set_yticks([5000, 10000, 20000, 50000, 85000 ])
    ax.set_ylim(5000,85000)
    ax.invert_yaxis()
    #labels
    ax.set_title(title,fontsize=18)
    ax.set_xlabel('lon')
    ax.set_ylabel('Pa')