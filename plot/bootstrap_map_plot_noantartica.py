#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 15:19:51 2023

@author: ema
"""
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


import constants as c

def bootstrap_map_plot_no_antartica(ds, ds_sign, levels, title=None, sign=0.95):
    
    #plot significance
    if sign == 0.95:
        min_v = 0
        max_v = -1
    elif sign == 0.90:
        min_v = 1
        max_v = -2
    elif sign == 0.80:
        min_v = 2
        max_v = -3
        
    #ds = ds.where((ds[c.VAR] <= ds_sign[c.VAR][2,:,:]) | (ds[c.VAR] >= ds_sign[c.VAR][-3,:,:]))
    neg = np.where(ds[c.VAR] <= ds_sign[c.VAR][min_v,:,:])
    pos = np.where(ds[c.VAR] >= ds_sign[c.VAR][max_v,:,:])
    lons, lats = np.meshgrid(ds.lon, ds.lat) 
    
#%%plot fields
    fig = plt.figure(figsize=[12,8])
    ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
    
    p = ds[c.VAR].plot(ax=ax, levels=levels, transform=ccrs.PlateCarree(),
                   add_labels=False, add_colorbar=False,
                   extend='both'
                   )
    
    _ = ax.scatter(lons[neg], lats[neg], marker = '.', s = 1, c = 'k', alpha = 0.2, transform = ccrs.PlateCarree())
    _ = ax.scatter(lons[pos], lats[pos], marker = '.', s = 1, c = 'k', alpha = 0.2, transform = ccrs.PlateCarree())
    
    #add coastlines
    ax.coastlines()
    
    ax.set_ylim([-60,90])
    
    # add separate colorbar
    cb = plt.colorbar(p, ticks=levels, shrink=0.6, extend='both')
    cb.set_label(c.UNITS,fontsize=18)
    cb.ax.tick_params(labelsize=18)

    #Drow gridlines and adjust labels
    gl = p.axes.gridlines(crs=ccrs.PlateCarree(), draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False
    gl.ylocator = mticker.FixedLocator([-80, -60, -40, -20, 0, 20, 40, 60, 80])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 18, 'color': 'black'}
    gl.ylabel_style = {'size': 18, 'color': 'black'}
    
    #labels
    ax.set_title(title,fontsize=18)