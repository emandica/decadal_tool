#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 16:21:06 2023

@author: ema
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

k = np.array((0.5,0.458,0.456,0.351,0.381,0.396,0.390,0.5,0.5,0.375,0.5,0.5,0.5,0.419,0.5,0.5,0.438,0.448,0.5,0.5,0.5))
look_up=np.array((0.90,0.85,0.90,0.90,0.90,0.99,0.70,0,0.50,0.90,0.10,0,0.60,0,0,0.5,0.5,0.90,0.90,0.60))

a52o_land_file = '/mnt/d/00-dataset/02-confess/03-confess_data/a52o/land/fc00_lead_0.nc'
dset = xr.open_dataset(a52o_land_file,chunks={'time':'auto'})

slm = dset['LSM'].mean('time')

cvh = xr.where(slm>0.5, dset['CVH'],np.nan)
cvl = xr.where(slm>0.5, dset['CVL'],np.nan)
lai_hv = xr.where(slm>0.5, dset['LAI_HV'],np.nan)
lai_lv = xr.where(slm>0.5, dset['LAI_LV'],np.nan)
tvh = dset['TVH'].round().astype(int).transpose('lat','lon','time')
tvl = dset['TVL'].round().astype(int).transpose('lat','lon','time')

cvh = cvh*(1-np.exp(-k[tvh]*lai_hv))
cvl = cvl*(1-np.exp(-k[tvl]*lai_lv))

#####
a2eq_land_file = '/mnt/d/00-dataset/ICMCL_dec_vegsens_198108201012_regular.nc'
dset_a2eq = xr.open_dataset(a2eq_land_file,chunks={'time':'auto'})

#a2eq_cvh = xr.where(slm>0.5, dset_a2eq['cvh'], np.nan)
a2eq_cvh = dset_a2eq['cvh']
#a2eq_cvl = xr.where(slm>0.5, dset_a2eq['cvl'],np.nan)
a2eq_cvl = dset_a2eq['cvl']
#a2eq_lai_hv = xr.where(slm>0.5, dset_a2eq['lai_hv'],np.nan)
a2eq_lai_hv = dset_a2eq['lai_hv']
#a2eq_lai_lv = xr.where(slm>0.5, dset_a2eq['lai_lv'],np.nan)
a2eq_lai_lv = dset_a2eq['lai_lv']

a2eq_cvh = a2eq_cvh*(1-np.exp(-0.5*a2eq_lai_hv))
a2eq_cvl = a2eq_cvl*(1-np.exp(-0.5*a2eq_lai_lv))

#####
a1ua_land_file = '/mnt/d/00-dataset/icmcl_hist_regular.nc'
dset_a1ua = xr.open_dataset(a1ua_land_file,chunks={'time':'auto'})

#a1ua_cvh = xr.where(slm>0.5, dset_a1ua['cvh'],np.nan)
a1ua_cvh = dset_a1ua['cvh']
#a1ua_cvl = xr.where(slm>0.5, dset_a1ua['cvl'],np.nan)
a1ua_cvl = dset_a1ua['cvl']
#a1ua_lai_hv = xr.where(slm>0.5, dset_a1ua['lai_hv'],np.nan)
a1ua_lai_hv = dset_a1ua['lai_hv']
#a1ua_lai_lv = xr.where(slm>0.5, dset_a1ua['lai_lv'],np.nan)
a1ua_lai_lv = dset_a1ua['lai_lv']
a1ua_tvh = dset_a1ua['tvh'].round().astype(int)
a1ua_tvl = dset_a1ua['tvl'].round().astype(int)

a1ua_cvh = a1ua_cvh#*look_up[a1ua_tvh]
a1ua_cvl = a1ua_cvl#*look_up[a1ua_tvl]


weights = np.cos(np.deg2rad(dset.lat))
weights.name = "weights"

###############################################################################
### CVH ###
#%%
levels=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
title='a52o, CVH, clim'

fig = plt.figure(figsize=[12,8])
ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
p=cvh.sel(time=slice('1993', '2006')).mean('time').plot(levels=levels,transform=ccrs.PlateCarree(),
               add_labels=False, add_colorbar=False,cmap='Greens')

ax.coastlines()
cb = plt.colorbar(p, ticks=levels, shrink=0.6)
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
ax.set_title(title,fontsize=18)

#%%
title='a2eq, CVH, clim'

fig = plt.figure(figsize=[12,8])
ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
p=a2eq_cvh.sel(time=slice('1993', '2006')).mean('time').plot(levels=levels,transform=ccrs.PlateCarree(),
               add_labels=False, add_colorbar=False,cmap='Greens')

ax.coastlines()
cb = plt.colorbar(p, ticks=levels, shrink=0.6)
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
ax.set_title(title,fontsize=18)

#%%
title='a1ua, CVH, clim'

fig = plt.figure(figsize=[12,8])
ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
p=a1ua_cvh.sel(time=slice('1993', '2006')).mean('time').plot(levels=levels,transform=ccrs.PlateCarree(),
               add_labels=False, add_colorbar=False,cmap='Greens')

ax.coastlines()
cb = plt.colorbar(p, ticks=levels, shrink=0.6)
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
ax.set_title(title,fontsize=18)

###############################################################################
### CVL ###
#%%
levels=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
title='a52o, CVL, clim'

fig = plt.figure(figsize=[12,8])
ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
p=cvl.sel(time=slice('1993', '2006')).mean('time').plot(levels=levels,transform=ccrs.PlateCarree(),
               add_labels=False, add_colorbar=False,cmap='Greens')

ax.coastlines()
cb = plt.colorbar(p, ticks=levels, shrink=0.6)
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
ax.set_title(title,fontsize=18)

#%%
title='a2eq, CVL, clim'

fig = plt.figure(figsize=[12,8])
ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
p=a2eq_cvl.sel(time=slice('1993', '2006')).mean('time').plot(levels=levels,transform=ccrs.PlateCarree(),
               add_labels=False, add_colorbar=False,cmap='Greens')

ax.coastlines()
cb = plt.colorbar(p, ticks=levels, shrink=0.6)
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
ax.set_title(title,fontsize=18)

#%%
title='a1ua, CVL, clim'

fig = plt.figure(figsize=[12,8])
ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
p=a1ua_cvl.sel(time=slice('1993', '2006')).mean('time').plot(levels=levels,transform=ccrs.PlateCarree(),
               add_labels=False, add_colorbar=False,cmap='Greens')

ax.coastlines()
cb = plt.colorbar(p, ticks=levels, shrink=0.6)
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
ax.set_title(title,fontsize=18)

###############################################################################
### LAI_HV ###
#%%
levels=[0,1,2,3,4,5,6,7]
title='a52o, LAI_HV, clim'

fig = plt.figure(figsize=[12,8])
ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
p=lai_hv.sel(time=slice('1993', '2006')).mean('time').plot(levels=levels,transform=ccrs.PlateCarree(),
               add_labels=False, add_colorbar=False,cmap='Greens',extend='max')

ax.coastlines()
cb = plt.colorbar(p, ticks=levels, shrink=0.6,extend='max')
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
ax.set_title(title,fontsize=18)

#%%
title='a2eq, LAI_HV, clim'

fig = plt.figure(figsize=[12,8])
ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
p=a2eq_lai_hv.sel(time=slice('1993', '2006')).mean('time').plot(levels=levels,transform=ccrs.PlateCarree(),
               add_labels=False, add_colorbar=False,cmap='Greens',extend='max')

ax.coastlines()
cb = plt.colorbar(p, ticks=levels, shrink=0.6,extend='max')
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
ax.set_title(title,fontsize=18)

#%%
title='a1ua, LAI_HV, clim'

fig = plt.figure(figsize=[12,8])
ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
p=a1ua_lai_hv.sel(time=slice('1993', '2006')).mean('time').plot(levels=levels,transform=ccrs.PlateCarree(),
               add_labels=False, add_colorbar=False,cmap='Greens',extend='max')

ax.coastlines()
cb = plt.colorbar(p, ticks=levels, shrink=0.6,extend='max')
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
ax.set_title(title,fontsize=18)


#%%
title = 'global average lai_hv'

fig = plt.figure(figsize=[12,8])
ax = fig.add_subplot(111)

lai_hv.sel(time=slice('1993', '2020')).weighted(weights).mean('lat').mean('lon').plot(label='a52o')
a2eq_lai_hv.sel(time=slice('1993', '2020')).weighted(weights).mean('lat').mean('lon').plot(label='a2eq')
a1ua_lai_hv.sel(time=slice('1993', '2020')).weighted(weights).mean('lat').mean('lon').plot(label='a1ua')
ax.legend()
ax.set_title(title,fontsize=18)

#%%
title = 'tropical average lai_hv'

fig = plt.figure(figsize=[12,8])
ax = fig.add_subplot(111)

lai_hv.sel(time=slice('1993', '2020')).sel(lat=slice(-15,15)).weighted(weights).mean('lat').mean('lon').plot(label='a52o')
a2eq_lai_hv.sel(time=slice('1993', '2020')).sel(lat=slice(-15,15)).weighted(weights).mean('lat').mean('lon').plot(label='a2eq')
a1ua_lai_hv.sel(time=slice('1993', '2020')).sel(lat=slice(-15,15)).weighted(weights).mean('lat').mean('lon').plot(label='a1ua')
ax.legend()
ax.set_title(title,fontsize=18)

###############################################################################
### LAI_LV ###
#%%
levels=[0,1,2,3,4,5,6,7]
title='a52o, LAI_LV, clim'

fig = plt.figure(figsize=[12,8])
ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
p=lai_lv.sel(time=slice('1993', '2006')).mean('time').plot(levels=levels,transform=ccrs.PlateCarree(),
               add_labels=False, add_colorbar=False,cmap='Greens',extend='max')
slm = dset['LSM']

cvh = xr.where(slm>0.5, dset['CVH'],np.nan)
cvl = xr.where(slm>0.5, dset['CVL'],np.nan)
lai_hv = xr.where(slm>0.5, dset['LAI_HV'],np.nan)
lai_lv = xr.where(slm>0.5, dset['LAI_LV'],np.nan)
tvh = dset['TVH'].round().astype(int)
tvl = dset['TVL'].round().astype(int)

cvh = cvh*(1-np.exp(-k[tvh]*lai_hv))
cvl = cvl*(1-np.exp(-k[tvl]*lai_lv))
cb = 1-cvh-cvl
ax.coastlines()
cb = plt.colorbar(p, ticks=levels, shrink=0.6,extend='max')
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
ax.set_title(title,fontsize=18)

#%%
title='a2eq, LAI_LV, clim'

fig = plt.figure(figsize=[12,8])
ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
p=a2eq_lai_lv.sel(time=slice('1993', '2006')).mean('time').plot(levels=levels,transform=ccrs.PlateCarree(),
               add_labels=False, add_colorbar=False,cmap='Greens',extend='max')

ax.coastlines()
cb = plt.colorbar(p, ticks=levels, shrink=0.6,extend='max')
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
ax.set_title(title,fontsize=18)

#%%
title='a1ua, LAI_LV, clim'

fig = plt.figure(figsize=[12,8])
ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
p=a1ua_lai_lv.sel(time=slice('1993', '2006')).mean('time').plot(levels=levels,transform=ccrs.PlateCarree(),
               add_labels=False, add_colorbar=False,cmap='Greens',extend='max')

ax.coastlines()
cb = plt.colorbar(p, ticks=levels, shrink=0.6,extend='max')
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
ax.set_title(title,fontsize=18)
