#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 17:38:34 2023

@author: ema
"""

import xarray as xr

import matplotlib.pyplot as plt

import constants as c

import plot as pl
from quantile import quantile_calculation

#%%
def xsection_plot(lead_exp, season):
#%%
    """
    quantile calculation
    """
    sign_delta = quantile_calculation(lead_exp, season)
    
    """
    open xsection
    """
    ctl = xr.open_dataset(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_xsection.nc')
    ctl = ctl.mean('member')
    sens = xr.open_dataset(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_xsection.nc')
    sens = sens.mean('member')
    p_ctl = xr.open_dataarray(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_xsection_p_'+c.REF+'.nc')
    p_sens = xr.open_dataarray(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_xsection_p_'+c.REF+'.nc')
    
    delta = sens - ctl  
    
#%%    
    #levels=[-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5]
    #levels=[-5,-4,-3,-2,-1,0,1,2,3,4,5]
    levels1=[-10,-8,-6,-4,-2,-1.5,-1,-.5,0,.5,1,1.5,2,4,6,8,10]
    levels2=[-1.0,-.8,-.6,-.4,-.2,-.15,-.1,-.05,0,.05,.1,.15,.2,.4,.6,.8,1]
#%%    
    title=c.NAME_CTL+'_xsection_'+c.VAR+'_lead_'+str(lead_exp) + '_' + season + '_s'+str(c.T_START)+ '_' + c.REF
    pl.parametric_xsection_plot(ctl,p_ctl,levels=levels1,title=title, sign=0.95)
    plt.savefig(c.OUT_DIR + c.NAME_CTL+'_xsection_'+c.VAR+'_lead_'+str(lead_exp) + '_' + season + '_s'+str(c.T_START)+ '_' + c.REF + '.jpg', dpi=300, bbox_inches='tight')
#%%
    title=c.NAME_SENS+'_xsection_'+c.VAR+'_lead_'+str(lead_exp) + '_' + season + '_s'+str(c.T_START)+ '_' + c.REF
    pl.parametric_xsection_plot(sens,p_sens,levels=levels1,title=title, sign=0.95)
    plt.savefig(c.OUT_DIR + c.NAME_SENS+'_xsection_'+c.VAR+'_lead_'+str(lead_exp) + '_' + season + '_s'+str(c.T_START)+ '_' + c.REF + '.jpg', dpi=300, bbox_inches='tight')
#%%
    title=c.NAME_SENS+'-'+c.NAME_CTL+'_xsection_'+c.VAR+'_lead_'+str(lead_exp) + '_' + season + '_s'+str(c.T_START)+ '_' + c.REF
    pl.bootstrap_xsection_plot(delta,sign_delta,levels=levels2,title=title, sign = 0.95)
    plt.savefig(c.OUT_DIR + c.NAME_SENS+'-'+c.NAME_CTL+'_xsection_'+c.VAR+'_lead_'+str(lead_exp) + '_' + season + '_s'+str(c.T_START)+'_' + c.REF + '.jpg', dpi=300, bbox_inches='tight')        