#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 16:14:41 2022

@author: ema
"""

import xarray as xr
import numpy as np

###############################################################################
def where_try(tv):
    k = [0, 0.150, 0.020, 2.000, 2.000, 2.000, 2.000, 0.100, 0.013, 0.050,
            0.150, 0.050, 0.0013, 0.050, 0, 0, 0.100, 0.100, 2.000, 0.500, 0]
    n_tv = xr.where(tv==0,k[0],tv)
    n_tv = xr.where(n_tv==1,k[1],n_tv)
    n_tv = xr.where(n_tv==2,k[2],n_tv)
    n_tv = xr.where(n_tv==3,k[3],n_tv)
    n_tv = xr.where(n_tv==4,k[4],n_tv)
    n_tv = xr.where(n_tv==5,k[5],n_tv)
    n_tv = xr.where(n_tv==6,k[6],n_tv)
    n_tv = xr.where(n_tv==7,k[7],n_tv)
    n_tv = xr.where(n_tv==8,k[8],n_tv)
    n_tv = xr.where(n_tv==9,k[9],n_tv)
    n_tv = xr.where(n_tv==10,k[10],n_tv)
    n_tv = xr.where(n_tv==11,k[11],n_tv)
    n_tv = xr.where(n_tv==12,k[12],n_tv)
    n_tv = xr.where(n_tv==13,k[13],n_tv)
    n_tv = xr.where(n_tv==14,k[14],n_tv)
    n_tv = xr.where(n_tv==15,k[15],n_tv)
    n_tv = xr.where(n_tv==16,k[16],n_tv)
    n_tv = xr.where(n_tv==17,k[17],n_tv)
    n_tv = xr.where(n_tv==18,k[18],n_tv)
    n_tv = xr.where(n_tv==19,k[19],n_tv)
    n_tv = xr.where(n_tv==20,k[20],n_tv)
    return n_tv
###############################################################################
def snow_coverage(snow_density, snow_depth):
    DCR = 0.10

    csn = (1000*snow_depth/snow_density)/DCR
    csn = xr.where(csn>=1, 1, csn)

    return csn


def max_intercept(cb, ch, lai_hv, cl, lai_lv, wl):
    WLMAX = 0.0002

    wlm = WLMAX*(cb+ch*lai_hv+cl*lai_lv)
    mi = min(1, wl/wlm)
    return mi


def intercept_coverage(csn, mi=0):
    c3 = (1-csn)*mi
    return c3

###############################################################################
def low_veg_cover(csn, cvl, mi=0):
    c4 = (1 - csn)*(1 - mi)*cvl
    return c4


def snow_on_bare_cover(csn, cvh):
    c5 = csn*(1 - cvh)
    return c5


def hi_veg_cover(csn, cvh, mi=0):
    c6 = (1 - csn)*(1 - mi)*cvh
    return c6


def snow_on_hi(csn, cvh):
    c7 = csn*cvh
    return c7


def bare_soil_cover(csn, cvl, cvh, mi=0):
    c8 = (1 - csn)*(1 - mi)*(1-cvl-cvh)
    return c8

#%%
m_rl = [0, 0.150, 0.020, 2.000, 2.000, 2.000, 2.000, 0.100, 0.013, 0.050,
        0.150, 0.050, 0.0013, 0.050, 0, 0, 0.100, 0.100, 2.000, 0.500, 0]
#h_rl = [0, 0.015, 0.002, 2.000, 2.000, 2.000, 2.000, 0.010, 0.0013, 0.005,
#        0.015, 0.005, 0.00013, 0.005, 0, 0, 0.010, 0.010, 2.000, 0.050, 0]
k = np.array((0.5,0.458,0.456,0.351,0.381,0.396,0.390,0.5,0.5,0.375,0.5,0.5,0.5,0.419,0.5,0.5,0.438,0.448,0.5,0.5,0.5))

a52o_land_file = '/mnt/d/00-dataset/02-confess/03-confess_data/a52o/land/fc01_lead_0.nc'

season = 'DJF'

dset = xr.open_dataset(a52o_land_file)

cvh = dset['CVH']  # .groupby('time.season').mean(dim='time')
cvl = dset['CVL']  # .groupby('time.season').mean(dim='time')
lai_hv = dset['LAI_HV']  # .groupby('time.season').mean(dim='time')
lai_lv = dset['LAI_LV']  # .groupby('time.season').mean(dim='time')
tvh = dset['TVH'].round().astype(int)  # .groupby('time.season').mean(dim='time')
tvl = dset['TVL'].round().astype(int)  # .groupby('time.season').mean(dim='time')

fsr = dset['FSR'] # .groupby('time.season').mean(dim='time') #forecast surface roughness
rsn = dset['RSN']  # .groupby('time.season').mean(dim='time') #snow density
sd = dset['SD']  # .groupby('time.season').mean(dim='time') #snow depth

lsm = dset['LSM']
#%%############################################################################
###COVERS CALCULATION##########################################################
###############################################################################
### vegetation cover
cvh = cvh*(1-np.exp(-k[tvh]*lai_hv))
cvl = cvl*(1-np.exp(-k[tvl]*lai_lv))
cb = 1-cvh-cvl

###snow cover
csn = snow_coverage(rsn, sd)  

###interception
#c1 = intercept_coverage()

#%%############################################################################
###TILE COVER##################################################################
###############################################################################
###3)wet skin

###4)low veg
cl = low_veg_cover(csn, cvl, mi=0)

###5)snow on low and bare
sob = snow_on_bare_cover(csn, cvh)

###6)high veg
ch = hi_veg_cover(csn, cvh, mi=0)

###7)snow high
soh = snow_on_hi(csn, cvh)

###8)bare soil
bs = bare_soil_cover(csn, cvl, cvh, mi=0)
bs = xr.where(lsm>0.5, bs, 0)

#%%############################################################################
###ROUGHNESS LENGTH############################################################
###############################################################################
BLEND_C = 10
###3)wet skin
#ws = np.zeros((1,180,360))
#for step in range(1):
#    for i in range(180):
#        for j in range(360):
#            ws[step,i,j] = cvl[step,i,j]*m_rl[tvh[step,i,j].to_numpy()]+cvh[step,i,j]*m_rl[tvl[step,i,j].to_numpy()]+(1-cvl[step,i,j]-cvh[step,i,j])*m_rl[0]

###4)low veg
tvl=tvl.to_dataset(name='roughness_length')
rl_low = tvl.map(where_try)
n_rl_low = cl/(np.log(BLEND_C/rl_low)**2)

###5)exposed snow
rl_es = sob/(np.log(BLEND_C/m_rl[12])**2)
rl_es=rl_es.to_dataset(name='roughness_length')

###6)high veg
tvh=tvh.to_dataset(name='roughness_length')
rl_hi = tvh.map(where_try)
n_rl_hi = ch/(np.log(BLEND_C/rl_hi)**2)

###7) sheltered snow
n_ss = soh/(np.log(BLEND_C/rl_hi)**2)

###8) bare soil
rl_bs = bs/(np.log(BLEND_C/m_rl[8])**2)
rl_bs=rl_bs.to_dataset(name='roughness_length')

#%%############################################################################
###BLENDING###
###############################################################################
rl_almost_final = (n_rl_low + rl_es + n_rl_hi + n_ss + rl_bs)

rl_final = BLEND_C*np.exp(-1/np.sqrt(rl_almost_final))
