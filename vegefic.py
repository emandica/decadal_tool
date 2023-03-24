#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 15:18:38 2023

@author: ema
"""

import xarray as xr
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

import generic_routines as cr
import b_bootstrap as b
import constants as c
import plot as pl
import run_bootstrap as rb
import generic_routines as ncr

directory = '/mnt/d/00-dataset/02-confess/'
var = 't2m'
units='K'
season = 'DJF'
lead_exp=[0]

file_ctrl = directory + c.NAME_CTL + '/' + var + '/msmm_' + c.NAME_CTL + '-' + var + '_' + season + '_sm_1x1.nc'
file_sens = directory + c.NAME_SENS + '/' + var + '/msmm_' + c.NAME_SENS + '-' + var + '_' + season + '_sm_1x1.nc'
ctrl = xr.open_dataset(file_ctrl, chunks={'lon':'auto','lat':'auto'})
sens = xr.open_dataset(file_sens, chunks={'lon':'auto','lat':'auto'})
ctrl = ctrl.rename({'t2m': 'tas','number':'member'})
sens = sens.rename({'t2m': 'tas','number':'member'})

ref = ncr.open_era5(c.DIRECTORY_ERA, c.FILE_ERA, c.VAR)
ref = ncr.dataset_season(ref[c.VAR], season)
ref= ref.sel(time=slice('1983','2015'))
ref['time']=ctrl['time']

ctrl_b=ctrl.mean('time')
sens_b=sens.mean('time')

#%%bootstrap
rb.bootstrap(ctrl_b,sens_b,lead_exp,season,rep=1000)

ctrl_b = xr.open_dataset(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_bootstrap_standard.nc')
sens_b = xr.open_dataset(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_bootstrap_standard.nc')

delta = sens_b-ctrl_b

sign_delta = delta.quantile([0.025,0.05,0.10,0.90,0.95,0.975], dim='bootstrap')

ctrl=ctrl.mean('member')
sens=sens.mean('member')

_, p_ctrl = stats.ttest_ind(ctrl[c.VAR].to_numpy(), ref.to_numpy(), axis=0, equal_var=False)    
p_ctrl = xr.DataArray(p_ctrl, coords=ctrl[c.VAR].mean('time').coords)
ctrl = ctrl.mean('time')

_, p_sens = stats.ttest_ind(sens[c.VAR].to_numpy(), ref.to_numpy(), axis=0, equal_var=False)    
p_sens = xr.DataArray(p_sens, coords=sens[c.VAR].mean('time').coords)
sens = sens.mean('time')

ref=ref.mean('time')

delta = sens - ctrl    
ctrl = ctrl-ref
sens=sens-ref

#%%tas
if c.VAR == 'tas':
    levels1=[-5,-4,-3,-2,-1,-0.5,0.5,1,2,3,4,5]
    levels2=[-2.5,-2,-1.5,-1,-0.5,-0.1,0.1,0.5,1,1.5,2,2.5]
    delta=xr.where(delta.apply(np.fabs)<0.1, np.NaN , delta)
#%%psl
elif c.VAR == 'psl':
    levels1=[-500,-400,-300,-200,-100,0,100,200,300,400,500]
    levels2=[-120,-100,-80,-60,-40,-20,-10,10,20,40,60,80,100,120]
    #delta=xr.where(delta.apply(np.fabs)<10, np.NaN , delta)
#%%ua        
elif c.VAR == 'ua':
    levels1=[-2,-1.5,-1,-0.5,0,0.5,1,1.5,2]
    levels2=[-0.5,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5]

title=c.NAME_CTL+'_MBIAS_'+c.VAR+'_lead_' + season + '_s'+str(c.T_START)+ '_' + c.REF    
pl.parametric_map_plot_no_antartica(ctrl,p_ctrl,levels=levels1, title=title, sign=0.95)
plt.savefig(c.OUT_DIR + c.NAME_CTL +'_MBIAS_'+var+'_'+ season + '_s'+str(c.T_START)+ '_' + c.REF + '_'+str(c.SIGN)+'_no_antartica.jpg', dpi=300, bbox_inches='tight')

pl.parametric_map_plot_polar(ctrl,p_ctrl,levels=levels1,title=title, sign=0.95)
plt.savefig(c.OUT_DIR + c.NAME_CTL+'_MBIAS_'+var+'_'+ season + '_s'+str(c.T_START)+ '_' + c.REF + '_'+str(c.SIGN)+'_polar.jpg', dpi=300, bbox_inches='tight')

#%%
title=c.NAME_SENS+'_MBIAS_'+c.VAR+'_lead_'+ season + '_s'+str(c.T_START)+ '_' + c.REF
pl.parametric_map_plot_no_antartica(sens,p_sens,levels=levels1,title=title,sign=c.SIGN)
plt.savefig(c.OUT_DIR + c.NAME_SENS+'_MBIAS_'+var+'_' + season + '_s'+str(c.T_START)+ '_' + c.REF + '_'+str(c.SIGN)+'_no_antartica.jpg', dpi=300, bbox_inches='tight')
    
pl.parametric_map_plot_polar(sens,p_sens,levels=levels1,title=title,sign=c.SIGN)
plt.savefig(c.OUT_DIR + c.NAME_SENS+'_MBIAS_'+var+'_' + season + '_s'+str(c.T_START)+ '_' + c.REF + '_'+str(c.SIGN)+'_polar.jpg', dpi=300, bbox_inches='tight')

#%%
title=c.NAME_SENS+'-'+c.NAME_CTL+'_MBIAS_'+c.VAR+'_'+ season + '_s'+str(c.T_START)+ '_' + c.REF
pl.bootstrap_map_plot_no_antartica(delta,sign_delta,levels=levels2,title=title, sign = c.SIGN)
plt.savefig(c.OUT_DIR + c.NAME_SENS +'-'+c.NAME_CTL+'_MBIAS_'+var+'_'+ season + '_s'+str(c.T_START)+'_' + c.REF + '_'+str(c.SIGN)+'_no_antartica.jpg', dpi=300, bbox_inches='tight')    
    
pl.bootstrap_map_plot_polar(delta,sign_delta,levels=levels2,title=title, sign = c.SIGN)
plt.savefig(c.OUT_DIR +c.NAME_SENS+'-'+c.NAME_CTL+'_MBIAS_'+var+'_'+season + '_s'+str(c.T_START)+'_' + c.REF + '_'+str(c.SIGN)+'_polar.jpg', dpi=300, bbox_inches='tight')
