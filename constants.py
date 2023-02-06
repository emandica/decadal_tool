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
VAR =
UNITS =
DIM =
PLEV =


###calculation information
no_member = [0] #member to exclude


###flag for reference
REF = 'era'

if REF == 'era':
    """era5"""
    ERA_VAR =
    DIRECTORY_ERA =
    FILE_ERA =

elif REF == 'hadcrut':
    """Observations"""
    DIRECTORY_OBS =
    FILE_OBS =

###Directories
FILE_DIR =

RUN_DIR =

OUT_DIR =

###############################################################################
###flags for diagnostics###
###############################################################################
#mean bias
M_BIAS = False

##############################################################################
#crossection
XSECT = False

###############################################################################
#correlation
CORRELATION = False

###############################################################################
#mean squerad skill score
MSSS = False

###############################################################################
#bootstrap method
BOOTSTRAP = 'standard'
