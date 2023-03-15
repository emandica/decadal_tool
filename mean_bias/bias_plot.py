#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 27 14:16:12 2023

@author: ema
"""

import xarray as xr

import numpy as np

import matplotlib.pyplot as plt

import constants as c

import plot as pl

def bias_plot(lead_exp, season):
    ctl = xr.open_dataset(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_bootstrap_mbias_'+c.REF+'.nc')
    sens = xr.open_dataset(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_bootstrap_mbias_'+c.REF+'.nc')

#%%
    """
    bias difference
    """
    delta = sens - ctl
#%%
    """
    quantile calculation
    """
    sign_delta = delta.quantile([0.025,0.05,0.10,0.90,0.95,0.975], dim='bootstrap')
    
#%%
    """
    open mean bias
    """
    ctl = xr.open_dataset(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_MBIAS.nc')
    sens = xr.open_dataset(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_MBIAS.nc')
    p_ctl = xr.open_dataarray(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_bias_p_'+c.REF+'.nc')
    p_sens = xr.open_dataarray(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_bias_p_'+c.REF+'.nc')
    ref = xr.open_dataset(c.RUN_DIR+c.REF+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc')    
    ref = ref.mean('time')
    delta = sens - ctl    
    ctl = ctl-ref
    sens=sens-ref
#%% 
    #%%tas
    if c.VAR == 'tas':
        levels1=[-5,-4,-3,-2,-1,0,1,2,3,4,5]
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
    
    title=c.NAME_CTL+'_MBIAS_'+c.VAR+'_lead_'+str(lead_exp) + '_' + season + '_s'+str(c.T_START)+ '_' + c.REF
    pl.parametric_map_plot(ctl,p_ctl,levels=levels1,title=title, sign=0.95)
    plt.savefig(c.OUT_DIR + c.NAME_CTL+'_MBIAS_'+c.VAR+'_lead_'+str(lead_exp) + '_' + season + '_s'+str(c.T_START)+ '_' + c.REF + '_'+str(c.SIGN)+'.jpg', dpi=300, bbox_inches='tight')

    pl.parametric_map_plot_no_antartica(ctl,p_ctl,levels=levels1,title=title, sign=0.95)
    plt.savefig(c.OUT_DIR + c.NAME_CTL+'_MBIAS_'+c.VAR+'_lead_'+str(lead_exp) + '_' + season + '_s'+str(c.T_START)+ '_' + c.REF + '_'+str(c.SIGN)+'_no_antartica.jpg', dpi=300, bbox_inches='tight')
    
    pl.parametric_map_plot_polar(ctl,p_ctl,levels=levels1,title=title, sign=0.95)
    plt.savefig(c.OUT_DIR + c.NAME_CTL+'_MBIAS_'+c.VAR+'_lead_'+str(lead_exp) + '_' + season + '_s'+str(c.T_START)+ '_' + c.REF + '_'+str(c.SIGN)+'_polar.jpg', dpi=300, bbox_inches='tight')
    
    #%%
    title=c.NAME_SENS+'_MBIAS_'+c.VAR+'_lead_'+str(lead_exp) + '_' + season + '_s'+str(c.T_START)+ '_' + c.REF
    pl.parametric_map_plot(sens,p_sens,levels=levels1,title=title,sign=c.SIGN)
    plt.savefig(c.OUT_DIR + c.NAME_SENS+'_MBIAS_'+c.VAR+'_lead_'+str(lead_exp) + '_' + season + '_s'+str(c.T_START)+ '_' + c.REF + '_'+str(c.SIGN)+'.jpg', dpi=300, bbox_inches='tight')
    
    pl.parametric_map_plot_no_antartica(sens,p_sens,levels=levels1,title=title,sign=c.SIGN)
    plt.savefig(c.OUT_DIR + c.NAME_SENS+'_MBIAS_'+c.VAR+'_lead_'+str(lead_exp) + '_' + season + '_s'+str(c.T_START)+ '_' + c.REF + '_'+str(c.SIGN)+'_no_antartica.jpg', dpi=300, bbox_inches='tight')
    
    pl.parametric_map_plot_polar(sens,p_sens,levels=levels1,title=title,sign=c.SIGN)
    plt.savefig(c.OUT_DIR + c.NAME_SENS+'_MBIAS_'+c.VAR+'_lead_'+str(lead_exp) + '_' + season + '_s'+str(c.T_START)+ '_' + c.REF + '_'+str(c.SIGN)+'_polar.jpg', dpi=300, bbox_inches='tight')
    
    #%%
    title=c.NAME_SENS+'-'+c.NAME_CTL+'_MBIAS_'+c.VAR+'_lead_'+str(lead_exp) + '_' + season + '_s'+str(c.T_START)+ '_' + c.REF
    pl.bootstrap_map_plot(delta,sign_delta,levels=levels2,title=title, sign = c.SIGN)
    plt.savefig(c.OUT_DIR + c.NAME_SENS+'-'+c.NAME_CTL+'_MBIAS_'+c.VAR+'_lead_'+str(lead_exp) + '_' + season + '_s'+str(c.T_START)+'_' + c.REF + '_'+str(c.SIGN)+'.jpg', dpi=300, bbox_inches='tight')

    pl.bootstrap_map_plot_no_antartica(delta,sign_delta,levels=levels2,title=title, sign = c.SIGN)
    plt.savefig(c.OUT_DIR + c.NAME_SENS+'-'+c.NAME_CTL+'_MBIAS_'+c.VAR+'_lead_'+str(lead_exp) + '_' + season + '_s'+str(c.T_START)+'_' + c.REF + '_'+str(c.SIGN)+'_no_antartica.jpg', dpi=300, bbox_inches='tight')    
    
    pl.bootstrap_map_plot_polar(delta,sign_delta,levels=levels2,title=title, sign = c.SIGN)
    plt.savefig(c.OUT_DIR + c.NAME_SENS+'-'+c.NAME_CTL+'_MBIAS_'+c.VAR+'_lead_'+str(lead_exp) + '_' + season + '_s'+str(c.T_START)+'_' + c.REF + '_'+str(c.SIGN)+'_polar.jpg', dpi=300, bbox_inches='tight')
