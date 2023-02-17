#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 13:31:50 2023

@author: ema
"""

import xarray as xr

import matplotlib.pyplot as plt

import constants as c

import plot as pl

def msss_plot(lead_exp, season):
    ctl = xr.open_dataset(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_bootstrap_mse_'+c.REF+'.nc')
    sens = xr.open_dataset(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_bootstrap_mse_'+c.REF+'.nc')

#%%
    """
    mse difference
    """
    delta = sens - ctl
#%%
    """
    quantile calculation
    """
    sign_delta = delta.quantile([0.025,0.05,0.10,0.90,0.95,0.975], dim='bootstrap')
#%%
    """
    open mean mse
    """
    ctl = xr.open_dataset(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_mse.nc')
    sens = xr.open_dataset(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_mse.nc')
    ref = xr.open_dataset(c.RUN_DIR+c.REF+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_mse.nc')
    p_ctl = xr.open_dataarray(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_mse_p_'+c.REF+'.nc')
    p_sens = xr.open_dataarray(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_mse_p_'+c.REF+'.nc')  

    ctl = ctl.mean('member')
    sens = sens.mean('member')

    delta = sens - ctl    

    levels1=[-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30]
    levels2=[-11,-9,-7,-5,-3,-1,0,1,3,5,7,9,11]
    
    title=c.NAME_CTL+'_mse_'+c.VAR+'_lead_'+str(lead_exp) + '_' + season + '_s'+str(c.T_START)+ '_' + c.REF
    pl.parametric_map_plot(ctl,p_ctl,levels=levels1,title=title, sign=0.95)
    plt.savefig(c.OUT_DIR + c.NAME_CTL+'_MSE_'+c.VAR+'_lead_'+str(lead_exp) + '_' + season + '_s'+str(c.T_START)+ '_' + c.REF + '.jpg', dpi=300, bbox_inches='tight')
    
    title=c.NAME_SENS+'_mse_'+c.VAR+'_lead_'+str(lead_exp) + '_' + season + '_s'+str(c.T_START)+ '_' + c.REF
    pl.parametric_map_plot(sens, p_sens, levels=levels1, title=title, sign=0.95)
    plt.savefig(c.OUT_DIR + c.NAME_SENS+'_MSE_'+c.VAR+'_lead_'+str(lead_exp) + '_' + season + '_s'+str(c.T_START)+ '_' + c.REF + '.jpg', dpi=300, bbox_inches='tight')
    
    title=c.NAME_SENS+'-'+c.NAME_CTL+'_MSE_'+c.VAR+'_lead_'+str(lead_exp) + '_' + season + '_s'+str(c.T_START)+ '_' + c.REF
    pl.bootstrap_map_plot(delta,sign_delta,levels=levels2,title=title, sign = 0.95)
    plt.savefig(c.OUT_DIR + c.NAME_SENS+'-'+c.NAME_CTL+'_MSE_'+c.VAR+'_lead_'+str(lead_exp) + '_' + season + '_s'+str(c.T_START)+'_' + c.REF + '.jpg', dpi=300, bbox_inches='tight')    