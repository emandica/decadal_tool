#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 20:43:24 2022

@author: ema
"""
###constants definitiontion
###exoperiments data
NAME_SENS =
NAME_CTL =
T_START =
SENS_NUMBER =
CTL_NUMBER =
T_END =
VAR =
UNITS =

###calculation information
no_member = #member to exclude

###flag for reference
REF =

if REF == 'era':
    """era5"""
    ERA_VAR =
    DIRECTORY_ERA =
    FILE_ERA = ERA_VAR+""

elif REF == 'hadcrut':
    """Observations"""
    DIRECTORY_OBS =
    FILE_OBS =

###Directories
FILE_DIR =

RUN_DIR =

OUT_DIR =

#LAND_FILE = '/mnt/d/dataset/02-confess/03-confess_data/a52o/land/fc00_lead_0_r.nc'

###############################################################################
###flags for diagnostics###
###############################################################################
#mean bias
M_BIAS = True

###############################################################################
#correlation
CORRELATION = False

###############################################################################
#mean squerad skill score
MSSS = False

###############################################################################
#bootstrap method
BOOTSTRAP = 'standard'
