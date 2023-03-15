#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 15:12:03 2023

@author: ema
"""
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

import constants as c
import data_agg as da
import level_selection as ls

season = 'Clim' #args.SEASON
lead_exp = [2,3,4] #args.LEAD_LIST

#%%data aggregation
print('starting data aggregation')

da.aggr_datasets(lead_exp, season)
print('data aggregation: OK')

print('starting mean bias')
ls.level_sel(lead_exp, season)


ctl = xr.open_dataset(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc',chunks={'lon':'auto','lat':'auto'})
sens = xr.open_dataset(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc',chunks={'lon':'auto','lat':'auto'})
ref = xr.open_dataset(c.RUN_DIR+c.REF+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc',chunks={'lon':'auto','lat':'auto'})

"""
temperature anomaly
"""
ctl=ctl-ctl.sel(time=slice(ctl['time'][0+lead_exp[0]], ctl['time'][(-1-lead_exp[0])])).mean('time')
sens=sens-sens.sel(time=slice(ctl['time'][0+lead_exp[0]],  ctl['time'][(-1-lead_exp[0])])).mean('time')
ref = ref-ref.sel(time=slice(ctl['time'][0+lead_exp[0]],  ctl['time'][(-1-lead_exp[0])])).mean('time')

"""
global mean
"""
lat_min = 50
lat_max = 60

lon_min = -14
lon_max = 0


LatIndexer, LonIndexer = 'lat', 'lon'
weights = np.cos(np.deg2rad(ctl[c.VAR].lat))
weights.name = "weights"

ctl = ctl[c.VAR].sel(**{LatIndexer: slice(lat_min,lat_max) , LonIndexer: slice(lon_min,lon_max)}).weighted(weights).mean(LatIndexer).mean(LonIndexer)
sens = sens[c.VAR].sel(**{LatIndexer: slice(lat_min,lat_max),LonIndexer:slice(lon_min,lon_max)}).weighted(weights).mean(LatIndexer).mean(LonIndexer)
ref = ref[c.VAR].sel(**{LatIndexer: slice(lat_min,lat_max),LonIndexer:slice(lon_min,lon_max)}).weighted(weights).mean(LatIndexer).mean(LonIndexer)


m_ctl = ctl.mean('member')
s_ctl =ctl.std('member')
corr_ctl=xr.corr(m_ctl, ref).to_numpy()

m_sens = sens.mean('member')
s_sens =sens.std('member')
corr_sens=xr.corr(m_sens, ref).to_numpy()
#%%
fig = plt.figure(figsize=[12,6])
ax = fig.add_subplot(111)
ax.fill_between(ctl.time,m_ctl-s_ctl,m_ctl+s_ctl,alpha=0.5)
m_ctl.plot.line(x='time',c='b',label=c.NAME_CTL+' r='+str(corr_ctl.round(2)))
ax.fill_between(sens.time,m_sens-s_sens,m_sens+s_sens,alpha=0.5)
m_sens.plot.line(x='time',c='r',label=c.NAME_SENS+' r='+str(corr_sens.round(2)))
ref.plot.line(x='time',c='k',label='era5')
ax.legend()
#ax.set_xlim(ctl.time[0],ctl.time[-1])
ax.set_title('area_mean_lat_'+str(lat_min)+'_'+str(lat_max)+'_lon_'+str(lon_min)+'_'+str(lon_max)+'_'+c.VAR+'_lead_'+str(lead_exp) + '_' + season + '_s'+str(c.T_START)+'_' + c.REF)

#%%
plt.savefig(c.OUT_DIR + 'area_mean_lat_'+str(lat_min)+'_'+str(lat_max)+'_lon_'+str(lon_min)+'_'+str(lon_max)+'_'+c.VAR+'_lead_'+str(lead_exp) + '_' + season + '_s'+str(c.T_START)+'_' + c.REF + '.jpg', dpi=300, bbox_inches='tight')    