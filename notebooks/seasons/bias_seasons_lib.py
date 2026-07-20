"""Plotting delle mappe di bias stagionali.

La funzione sta qui (e non nel notebook) perche' il plotting gira in processi
'spawn' freschi (ProcessPoolExecutor con max_tasks_per_child=1): ogni figura usa
un processo che esce subito dopo, cosi' il sistema operativo recupera tutta la
memoria (map_plot lascia memoria a ogni figura e ne' plt.close ne' gc la
recuperano). I processi spawn pero' non possono importare funzioni definite nel
notebook: devono stare in un modulo importabile come questo.
"""
import os
import sys

# rende config.py (in notebooks/) importabile anche da questa sottocartella
_cfg = os.path.dirname(os.path.abspath(__file__))
while _cfg != os.path.dirname(_cfg):
    if os.path.exists(os.path.join(_cfg, "config.py")):
        if _cfg not in sys.path:
            sys.path.insert(0, _cfg)
        break
    _cfg = os.path.dirname(_cfg)

import logging

import xarray as xr
import matplotlib
matplotlib.use("Agg")          # backend headless: niente ritenzione di figure
import matplotlib.pyplot as plt

from config import FIG_DIR
import albedo_functions as af

LEVELS = [-0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5]
DPI = 300   # il buffer di rendering scala col quadrato del dpi; 300 -> 4800x2700 px


def _load(save_path, fname, v):
    """Carica un array (piccolo) chiudendo subito il file."""
    with xr.open_dataset(f"{save_path}/{fname}") as ds:
        return ds[v].load()


def _plot_one(args):
    """Adattatore per Pool.imap_unordered (che passa un solo argomento: una tupla)."""
    return plot_bias_season(*args)


def plot_bias_season(exp_ctrl, exp_sens, var, y1, y2, season, save_path):
    """Mappe di bias (control / sensitivity) per un range di lead year e una stagione."""
    lead = f"{y1}-{y2}"
    try:
        bias_ctrl = _load(save_path, f"{exp_ctrl}_{var}_lead_{lead}_bias_{season}.nc", var)
        bias_sens = _load(save_path, f"{exp_sens}_{var}_lead_{lead}_bias_{season}.nc", var)
        p_ctrl = _load(save_path, f"{exp_ctrl}_{var}_lead_{lead}_bias_p_{season}.nc", "p")
        p_sens = _load(save_path, f"{exp_sens}_{var}_lead_{lead}_bias_p_{season}.nc", "p")

        af.map_plot(bias_ctrl, p_ctrl, levels=LEVELS,
                    title=f"a) DCPP-CTRL {season} {lead}", cmap="BrBG", antartica=False)
        plt.savefig(f"{FIG_DIR}/{exp_ctrl}_{var}_bias_{season}_{lead}.png",
                    dpi=DPI, bbox_inches="tight")
        plt.close("all")

        af.map_plot(bias_sens, p_sens, levels=LEVELS,
                    title=f"b) DCPP-SENS {season} {lead}", cmap="BrBG", antartica=False)
        plt.savefig(f"{FIG_DIR}/{exp_sens}_{var}_bias_{season}_{lead}.png",
                    dpi=DPI, bbox_inches="tight")
        plt.close("all")
        return f"{season} {lead} ok"
    except Exception as e:
        logging.exception(f"Error {season} {lead}: {e}")
        return f"{season} {lead} ERRORE: {e}"
