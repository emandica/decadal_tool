#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 17:34:39 2023

@author: ema
"""

import xarray as xr
import numpy as np
import xskillscore as xs
from scipy import stats
import matplotlib.pyplot as plt

import generic_routines as cr
import b_bootstrap as b
import constants as c
import plot as pl
import run_bootstrap as rb
import generic_routines as ncr
#%%
directory = '/mnt/d/00-dataset/02-confess/'
var = 'u'
units='m/s'
season = 'DJF'
lead_exp=[0]

file_ctrl = c.FILE_DIR + c.NAME_CTL + '/' + var + '850/msmm_' + c.NAME_CTL + '-' + var + '850_' + season + '_sm_1x1.nc'
file_sens = c.FILE_DIR + c.NAME_SENS + '/' + var + '850/msmm_' + c.NAME_SENS + '-' + var + '850_' + season + '_sm_1x1.nc'
ctrl = xr.open_dataset(file_ctrl, chunks={'lon':'auto','lat':'auto'})
sens = xr.open_dataset(file_sens, chunks={'lon':'auto','lat':'auto'})
ctrl = ctrl.rename({var: c.VAR,'number':'member'})
sens = sens.rename({var: c.VAR,'number':'member'})

ref = ncr.open_era5(c.DIRECTORY_ERA, c.FILE_ERA, c.VAR)
ref = ref.sel(plev = c.PLEV*100).squeeze()
ref = ncr.dataset_season(ref[c.VAR], season)
ref= ref.sel(time=slice('1983','2015'))
ref['time']=ctrl['time']
ref=ref.chunk(dict(time=-1))

ctrl = ctrl - ctrl.mean('time')
sens = sens - sens.mean('time')
ref = ref -ref.mean('time')

#%%bootstrap
rb.bootstrap(ctrl,sens,lead_exp,season,rep=1000)

ctrl_b = xr.open_dataset(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_bootstrap_standard.nc')
sens_b = xr.open_dataset(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_bootstrap_standard.nc')

delta1_corr = xr.corr(ctrl_b[c.VAR],ref,dim='time')
delta2_corr = xr.corr(sens_b[c.VAR],ref,dim='time')

delta = delta2_corr - delta1_corr
sign_delta = delta.quantile([0.025,0.05,0.10,0.90,0.95,0.975], dim='bootstrap')

#%%
ctrl=ctrl.mean('member')
sens=sens.mean('member')
ctrl_r = xr.corr(ctrl[c.VAR],ref,dim='time')
sens_r = xr.corr(sens[c.VAR],ref,dim='time')
sign_ctrl = xs.pearson_r_p_value(ctrl[c.VAR], ref, dim='time',keep_attrs=True)
sign_sens = xs.pearson_r_p_value(sens[c.VAR], ref, dim='time',keep_attrs=True)

sign_ctrl = sign_ctrl.transpose('lat','lon')
sign_sens = sign_sens.transpose('lat','lon')

delta=sens_r-ctrl_r

#%%    
levels1=[-0.4,-0.3,-0.2,-0.1,-0.05,0.05,0.1,0.2,0.3,0.4]
levels2=[-0.4,-0.3,-0.2,-0.1,-0.05,0.05,0.1,0.2,0.3,0.4]
    
title=c.NAME_CTL+'_ACC_'+c.VAR+'_lead_' + season + '_s'+str(c.T_START)+ '_' + c.REF    
pl.parametric_map_plot_no_antartica(ctrl_r,sign_ctrl,levels=levels1, title=title, sign=0.95)
plt.savefig(c.OUT_DIR + c.NAME_CTL +'_ACC_'+var+'_'+ season + '_s'+str(c.T_START)+ '_' + c.REF + '_'+str(c.SIGN)+'_no_antartica.jpg', dpi=300, bbox_inches='tight')

pl.parametric_map_plot_polar(ctrl_r,sign_ctrl,levels=levels1,title=title, sign=0.95)
plt.savefig(c.OUT_DIR + c.NAME_CTL+'_ACC_'+var+'_'+ season + '_s'+str(c.T_START)+ '_' + c.REF + '_'+str(c.SIGN)+'_polar.jpg', dpi=300, bbox_inches='tight')

#%%
title=c.NAME_SENS+'_ACC_'+c.VAR+'_lead_'+ season + '_s'+str(c.T_START)+ '_' + c.REF
pl.parametric_map_plot_no_antartica(sens_r,sign_sens,levels=levels1,title=title,sign=c.SIGN)
plt.savefig(c.OUT_DIR + c.NAME_SENS+'_ACC_'+var+'_' + season + '_s'+str(c.T_START)+ '_' + c.REF + '_'+str(c.SIGN)+'_no_antartica.jpg', dpi=300, bbox_inches='tight')
    
pl.parametric_map_plot_polar(sens_r,sign_sens,levels=levels1,title=title,sign=c.SIGN)
plt.savefig(c.OUT_DIR + c.NAME_SENS+'_ACC_'+var+'_' + season + '_s'+str(c.T_START)+ '_' + c.REF + '_'+str(c.SIGN)+'_polar.jpg', dpi=300, bbox_inches='tight')

#%%
title=c.NAME_SENS+'-'+c.NAME_CTL+'_ACC_'+c.VAR+'_'+ season + '_s'+str(c.T_START)+ '_' + c.REF
pl.bootstrap_map_plot_no_antartica(delta,sign_delta,levels=levels2,title=title, sign = c.SIGN)
plt.savefig(c.OUT_DIR + c.NAME_SENS +'-'+c.NAME_CTL+'_ACC_'+var+'_'+ season + '_s'+str(c.T_START)+'_' + c.REF + '_'+str(c.SIGN)+'_no_antartica.jpg', dpi=300, bbox_inches='tight')    
    
pl.bootstrap_map_plot_polar(delta,sign_delta,levels=levels2,title=title, sign = c.SIGN)
plt.savefig(c.OUT_DIR +c.NAME_SENS+'-'+c.NAME_CTL+'_ACC_'+var+'_'+season + '_s'+str(c.T_START)+'_' + c.REF + '_'+str(c.SIGN)+'_polar.jpg', dpi=300, bbox_inches='tight')

    
