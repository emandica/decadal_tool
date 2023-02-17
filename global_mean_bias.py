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

ctl = xr.open_dataset(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc')
sens = xr.open_dataset(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc')
ref = xr.open_dataset(c.RUN_DIR+c.REF+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc')

"""
temperature anomaly
"""
ctl=ctl-ctl.mean('time')
sens=sens-sens.mean('time')
ref = ref-ref.mean('time')

"""
global mean
"""
LatIndexer, LonIndexer = 'lat', 'lon'
weights = np.cos(np.deg2rad(ctl[c.VAR].lat))
weights.name = "weights"

ctl = ctl[c.VAR].sel(**{LatIndexer: slice(-90,90)}).weighted(weights).mean(LatIndexer).mean(LonIndexer)
sens = sens[c.VAR].sel(**{LatIndexer: slice(-90,90)}).weighted(weights).mean(LatIndexer).mean(LonIndexer)
ref = ref[c.VAR].sel(**{LatIndexer: slice(-90,90)}).weighted(weights).mean(LatIndexer).mean(LonIndexer)

m_ctl = ctl.mean('member')
s_ctl =ctl.std('member')
m_sens = sens.mean('member')
s_sens =sens.std('member')
#%%
fig = plt.figure(figsize=[12,8])
ax = fig.add_subplot(111)
ax.fill_between(ctl.time,m_ctl-s_ctl,m_ctl+s_ctl,alpha=0.5)
m_ctl.plot.line(x='time',c='b',label=c.NAME_CTL)
ax.fill_between(sens.time,m_sens-s_sens,m_sens+s_sens,alpha=0.5)
m_sens.plot.line(x='time',c='r',label=c.NAME_SENS)
ref.plot.line(x='time',c='k')
ax.legend()
