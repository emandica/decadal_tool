#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 16:58:41 2023

@author: ema
"""
import argparse

import run_xsec as rx
import run_mbias as rm
import run_correlation as rc
import run_msss as rmsss

import constants as c

#%%
# entry point for the program
if __name__ == '__main__':
###############################################################################
###constants definitiontion
    CLI=argparse.ArgumentParser()
    CLI.add_argument(
        "--LEAD_LIST",
        nargs="*",
        type=int,
        default=[],
        )

    CLI.add_argument(
            "-SEASON",
            #nargs="*",
            type=str,
            default=[],
            )

    args=CLI.parse_args()

    SEASON = args.SEASON
    lead_exp = args.LEAD_LIST

#%%Mean bias
    if c.M_BIAS:
        rm.run_mbias(lead_exp, SEASON)

#%%Correlation
    if c.CORRELATION:
        rc.run_corr(lead_exp, SEASON)

#%%MSSS
    if c.MSSS:
        rmsss.run_msss(lead_exp, SEASON)

#%%xsection
    if c.XSECT:
        rx.run_xsec(lead_exp, SEASON)
