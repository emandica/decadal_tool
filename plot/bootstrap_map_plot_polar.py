#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 15:19:51 2023

@author: ema
"""
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.path as mpath

import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


import constants as c

def bootstrap_map_plot_polar(ds, ds_sign, levels, title=None, sign=0.95):
    #ax.set_ylim([-2.5,2.5])
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
        
    ds = ds.where((ds[c.VAR] <= ds_sign[c.VAR][2,:,:]) | (ds[c.VAR] >= ds_sign[c.VAR][-3,:,:]))
    #neg = np.where(ds[c.VAR] <= ds_sign[c.VAR][min_v,:,:])
    #pos = np.where(ds[c.VAR] >= ds_sign[c.VAR][max_v,:,:])
    neg = np.where(ds[c.VAR] <= ds_sign[c.VAR][min_v,:,:],-1,0)
    pos = np.where(ds[c.VAR] >= ds_sign[c.VAR][max_v,:,:],1,0)
    
    lons, lats = np.meshgrid(ds.lon, ds.lat) 
    
#%%plot fields
    fig = plt.figure(figsize=[12,8])
    ax = fig.add_subplot(111, projection=ccrs.Orthographic(central_longitude=0,central_latitude=90))
    
    p = ds[c.VAR].plot(ax=ax, levels=levels, transform=ccrs.PlateCarree(),
                   add_labels=False, add_colorbar=False,
                   extend='both', 
                   #cmap='bwr',
                   )
    
    _ = ax.contourf(lons,lats,neg,levels=[-2,-1,1,2],hatches=["...","","..."],
                    transform = ccrs.PlateCarree(),alpha=0)
    _ = ax.contourf(lons,lats,pos,levels=[-1.1,-0.5,0.5,1.1],hatches=["...","","..."],
                    transform = ccrs.PlateCarree(),alpha=0)
    #_ = ax.scatter(lons[neg], lats[neg], marker = '.', s = 1, c = 'k', alpha = 0.2, transform = ccrs.PlateCarree())
    #_ = ax.scatter(lons[pos], lats[pos], marker = '.', s = 1, c = 'k', alpha = 0.2, transform = ccrs.PlateCarree())
    
    #add coastlines
    ax.coastlines()
    
    ax.set_extent([0,360,30,90], ccrs.PlateCarree())
    
    # add separate colorbar
    cb = plt.colorbar(p, ticks=levels, shrink=1, extend='both')
    cb.set_label(c.UNITS,fontsize=18)
    
    #cb.set_label(c.UNITS,fontsize=18)
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