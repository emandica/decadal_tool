#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 11:02:31 2023

@author: ema
"""

import xarray as xr

import matplotlib.pyplot as plt

import constants as c

import generic_routines as cr

import plot as pl

def correlation_plot(lead_exp, season):
    
    if c.BOOTSTRAP=='block':
        ctl = xr.open_dataset(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_bootstrap_correlation_'+c.REF+'.nc', chunks={'lon':'auto','lat':'auto'})
        sens = xr.open_dataset(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_bootstrap_correlation_'+c.REF+'.nc', chunks={'lon':'auto','lat':'auto'})

#%%
        """
        correlation difference
        """
        delta = sens - ctl

#%%
        """
        quantile calculation
        """
        sign_ctl = ctl.quantile([0.025,0.05,0.10,0.90,0.95,0.975], dim='bootstrap')
        sign_sens = sens.quantile([0.025,0.05,0.10,0.90,0.95,0.975], dim='bootstrap')
        sign_delta = delta.quantile([0.025,0.05,0.10,0.90,0.95,0.975], dim='bootstrap')
#%%
        """
        open correlations
        """
        ctl = xr.open_dataset(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_ACC.nc', chunks={'lon':'auto','lat':'auto'})
        sens = xr.open_dataset(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_ACC.nc', chunks={'lon':'auto','lat':'auto'})
        delta = sens - ctl

###############################################################################
#%%
    elif c.BOOTSTRAP=='standard':
        delta1 = xr.open_dataset(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_bootstrap_standard_correlation_'+c.REF+'.nc', chunks={'lon':'auto','lat':'auto'})
        delta2 = xr.open_dataset(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_bootstrap_standard_correlation_'+c.REF+'.nc', chunks={'lon':'auto','lat':'auto'})
        delta = delta1 - delta2
#%%
        """
        quantile calculation
        """
        sign_delta = delta.quantile([0.025,0.05,0.10,0.90,0.95,0.975], dim='bootstrap')
#%%
        """
        open correlations
        """
        ctl = xr.open_dataset(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_ACC.nc', chunks={'lon':'auto','lat':'auto'})
        sign_ctl = xr.open_dataset(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_correlation_p_'+c.REF+'.nc', chunks={'lon':'auto','lat':'auto'})
        
        sens = xr.open_dataset(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_ACC.nc', chunks={'lon':'auto','lat':'auto'})
        sign_sens = xr.open_dataset(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_correlation_p_'+c.REF+'.nc', chunks={'lon':'auto','lat':'auto'})
        
        delta=sens-ctl
#%%    
    levels1=[-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,0,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
    #levels2=[-0.5,-0.4,-0.3,-0.2,0,0.2,0.3,0.4,0.5]
    levels2=[-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,0,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
    
    title=c.NAME_CTL+'_ACC_'+c.VAR+'_lead_'+str(lead_exp) + '_' + season + '_s'+str(c.T_START)+ '_' + c.REF
    pl.parametric_map_plot(ctl,sign_ctl[c.VAR],levels=levels1, title=title, sign=0.95)
    plt.savefig(c.OUT_DIR + c.NAME_CTL+'_ACC_'+c.VAR+'_lead_'+str(lead_exp) + '_' + season + '_s'+str(c.T_START)+ '_' + c.REF + '_no_antartica.jpg', dpi=300, bbox_inches='tight')
    
    title=c.NAME_SENS+'_ACC_'+c.VAR+'_lead_'+str(lead_exp) + '_' + season + '_s'+str(c.T_START)+ '_' + c.REF
    pl.parametric_map_plot(sens,sign_sens[c.VAR],levels=levels1,title=title,sign=0.95)
    plt.savefig(c.OUT_DIR + c.NAME_SENS+'_ACC_'+c.VAR+'_lead_'+str(lead_exp) + '_' + season + '_s'+str(c.T_START)+ '_' + c.REF + '_no_antartica.jpg', dpi=300, bbox_inches='tight')
    
    title=c.NAME_SENS+'-'+c.NAME_CTL+'_ACC_'+c.VAR+'_lead_'+str(lead_exp) + '_' + season + '_s'+str(c.T_START)+ '_' + c.REF
    pl.bootstrap_map_plot(delta,sign_delta,levels=levels2,title=title, sign = 0.95)
    plt.savefig(c.OUT_DIR + c.NAME_SENS+'-'+c.NAME_CTL+'_ACC_'+c.VAR+'_lead_'+str(lead_exp) + '_' + season + '_s'+str(c.T_START)+'_' + c.REF + '_no_antartica.jpg', dpi=300, bbox_inches='tight')