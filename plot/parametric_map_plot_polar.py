#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 15:25:43 2023

@author: ema
"""

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.path as mpath

import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


import constants as c

def parametric_map_plot_polar(ds, ds_sign, levels, title=None, sign=0.95):
    
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
        
    #ds = ds.where((ds_sign <= 0.1) | (ds_sign >= 0.9))
    neg = np.where((ds_sign <= min_v)|(ds_sign >= max_v))
    lons, lats = np.meshgrid(ds.lon, ds.lat)    
    
# to use for correlation, to comment for bias
    #neg = (neg[1],neg[0])
    
#%%plot fields
    fig = plt.figure(figsize=[12,8])
    
    ax = fig.add_subplot(111, projection=ccrs.Orthographic(central_longitude=0,central_latitude=90))
    
    p = ds[c.VAR].plot(ax=ax, levels=levels, transform=ccrs.PlateCarree(),
                   add_labels=False, add_colorbar=False,
                   extend='both',
                   )
    
    _ = ax.scatter(lons[neg], lats[neg], marker = '.', s = 1, c = 'k',
                   alpha = 0.2, transform = ccrs.PlateCarree())
    
    #add coastlines
    ax.coastlines()
    ax.set_global()
    
    ax.set_extent([0,360,30,90], ccrs.PlateCarree())
    
    # add separate colorbar
    cb = plt.colorbar(p, ticks=levels, shrink=1, extend='both')
    
    cb.set_label(c.UNITS,fontsize=18)
    cb.ax.tick_params(labelsize=18)
    
    # Compute a circle in axes coordinates, which we can use as a boundary
# for the map. We can pan/zoom as much as we like - the boundary will be
# permanently circular.
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)
    
    ax.set_boundary(circle, transform=ax.transAxes)

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