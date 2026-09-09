"""Bias globale annuale (04-BIAS.ipynb): per combinazione di lead year (y1-y2).

Le funzioni stanno qui (e non nel notebook) perche' il calcolo gira in processi
'spawn' freschi (multiprocessing.Pool con maxtasksperchild=1): stesso motivo e
stessa soluzione gia' applicata a notebooks/seasons/bias_leadyears_lib.py - il
ProcessPoolExecutor con worker riusati (usato qui in precedenza) si e'
rivelato fragile su questo cluster (BrokenProcessPool: "A process in the
process pool was terminated abruptly"), probabile netCDF4/HDF5 su NFS con
processi forkati. I processi spawn non vedono le funzioni definite nel
notebook: devono stare in un modulo importabile come questo.
"""
import os
import sys

# rende config.py (in notebooks/) importabile
_cfg = os.path.dirname(os.path.abspath(__file__))
while _cfg != os.path.dirname(_cfg):
    if os.path.exists(os.path.join(_cfg, "config.py")):
        if _cfg not in sys.path:
            sys.path.insert(0, _cfg)
        break
    _cfg = os.path.dirname(_cfg)

# HDF5/netCDF4 su NFS + processi possono soffrire di problemi di locking;
# disabilita il file locking HDF5 (fix standard per NFS). Va fatto PRIMA di
# importare xarray/netCDF4.
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

import pandas as pd
import xarray as xr
import xskillscore as xs
import dask
from scipy import stats

from config import POST_DATA, WORK_DIR


def bootstrap_quantile_bias(sens, ctrl, ref, iterations=1000):
    """Identica alla funzione originale nel notebook (era la cella 3): logica
    non toccata, solo spostata qui perche' serve dentro process_lead_years."""
    with dask.config.set(**{"array.slicing.split_large_chunks": True}):
        f = xr.concat([sens, ctrl], dim="member").chunk({"member": -1})
        f_ra = xs.resampling.resample_iterations(f, iterations, "member", replace=True).mean("member").squeeze().compute()
        f_rb = xs.resampling.resample_iterations(f, iterations, "member", replace=True).mean("member").squeeze().compute()
        bias_f_ra = f_ra.mean("time") - ref
        bias_f_rb = f_rb.mean("time") - ref
        delta_bias = bias_f_ra - bias_f_rb
        sig_delta = delta_bias.chunk(dict(iteration=-1)).quantile([0.025, 0.05, 0.10, 0.90, 0.95, 0.975], dim="iteration")
    return sig_delta


def process_lead_years(exp_ctrl, exp_sens, var, era_var, y1, y2, save_path):
    """Identica alla funzione originale nel notebook (era la cella 4), stessa
    logica bit-per-bit. Unico cambiamento strutturale: 'DATA_PATH'/'OBS_PATH'/
    'era_var' erano variabili GLOBALI del notebook (invisibili a un processo
    spawn, che parte da un interprete pulito) - qui 'era_var' e' un parametro
    esplicito, 'DATA_PATH'/'OBS_PATH' sono sostituiti dai loro valori diretti
    (POST_DATA/WORK_DIR, cosi' erano definiti nel notebook originale: DATA_PATH
    = POST_DATA, OBS_PATH = WORK_DIR). Nessuna modifica di logica scientifica.

    Gli errori non vengono piu' loggati su file (logging.exception): i
    processi spawn non condividono l'handler di logging del notebook, quindi
    quei messaggi non arrivavano mai a destinazione. Si ritorna invece una
    stringa di esito, stampata dal notebook - stesso schema gia' in uso in
    bias_leadyears_lib.py e in tutta la pipeline cover_regression.

    NOTA (comportamento preesistente, non modificato qui): nell'originale,
    dset_ctrl_1999/dset_sens_1999/obs_1999 venivano calcolati (slice dal
    1999) ma MAI USATI - bias e bootstrap sotto operano sulla serie temporale
    COMPLETA. Se non era intenzionale andrebbe deciso a parte; qui il
    comportamento e' preservato identico all'originale."""
    lead = f"{y1}-{y2}"
    lead_number = y2 - y1 + 1
    try:
        dset_ctrl_path = POST_DATA / exp_ctrl / "1x1" / var / f"{exp_ctrl}_{var}_Amon_EC-Earth3_dcppA-hindcast_lead_{lead}_1x1_ensemble_rad.nc"
        dset_ctrl = xr.open_dataset(dset_ctrl_path)
        dset_ctrl["time"] = pd.to_datetime(dset_ctrl["time"].values).normalize().to_period("M").start_time

        dset_sens_path = POST_DATA / exp_sens / "1x1" / var / f"{exp_sens}_{var}_Amon_EC-Earth3_dcppA-hindcast_lead_{lead}_1x1_ensemble_rad.nc"
        dset_sens = xr.open_dataset(dset_sens_path)
        dset_sens["time"] = pd.to_datetime(dset_sens["time"].values).normalize().to_period("M").start_time

        obs_path = WORK_DIR / f"ERA5_{era_var}_1x1_{lead_number}year.nc"
        obs = xr.open_dataset(obs_path)
        if lead_number in [2, 4]:
            obs["time"] = pd.to_datetime(obs["time"].values) - pd.DateOffset(months=1)
        obs["time"] = pd.to_datetime(obs["time"].values).normalize().to_period("M").start_time
        obs = obs.rename({f"{era_var}": f"{var}"})
        obs = obs.sel(time=slice(dset_ctrl.time[0], dset_ctrl.time[-1]))

        bias_ctrl_time_series = dset_ctrl[var].mean("member") - obs[var]
        bias_ctrl = bias_ctrl_time_series.mean("time")

        bias_sens_time_series = dset_sens[var].mean("member") - obs[var]
        bias_sens = bias_sens_time_series.mean("time")

        _, ctrl_p_value = stats.ttest_1samp(bias_ctrl_time_series, popmean=0, axis=0)
        _, sens_p_value = stats.ttest_1samp(bias_sens_time_series, popmean=0, axis=0)

        ctrl_p_value = xr.DataArray(data=ctrl_p_value, coords=bias_ctrl.coords, dims=["lat", "lon"], name="p")
        sens_p_value = xr.DataArray(data=sens_p_value, coords=bias_sens.coords, dims=["lat", "lon"], name="p")

        delta = bootstrap_quantile_bias(dset_ctrl, dset_sens, obs.mean("time"))

        ctrl_outfile = f"{save_path}/{exp_ctrl}_{var}_lead_{lead}_bias.nc"
        sens_outfile = f"{save_path}/{exp_sens}_{var}_lead_{lead}_bias.nc"
        ctrl_p_outfile = f"{save_path}/{exp_ctrl}_{var}_lead_{lead}_bias_p.nc"
        sens_p_outfile = f"{save_path}/{exp_sens}_{var}_lead_{lead}_bias_p.nc"
        delta_outfile = f"{save_path}/delta_{var}_lead_{lead}_bias_quantile.nc"

        bias_ctrl.to_dataset(name=var).to_netcdf(ctrl_outfile)
        bias_sens.to_dataset(name=var).to_netcdf(sens_outfile)
        ctrl_p_value.to_netcdf(ctrl_p_outfile)
        sens_p_value.to_netcdf(sens_p_outfile)
        delta.to_netcdf(delta_outfile)

        return f"{lead} ok"
    except Exception as e:
        return f"{lead} ERRORE: {type(e).__name__}: {e}"


def run_one(args):
    """Adattatore per Pool.imap_unordered: un solo argomento (tupla)."""
    exp_ctrl, exp_sens, var, era_var, y1, y2, save_path = args
    return process_lead_years(exp_ctrl, exp_sens, var, era_var, y1, y2, save_path)
