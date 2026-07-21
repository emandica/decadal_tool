"""Bias globale/land/ocean stagionale: lead singolo (0..4) e combinazioni di
lead year (range contigui y1-y2).

Le funzioni stanno qui (e non nel notebook) perche' il calcolo gira in processi
'spawn' freschi (multiprocessing.Pool con maxtasksperchild=1): ogni task usa un
processo che esce subito dopo, cosi' il SO recupera tutta la memoria e non si
accumula stato tra un task e il successivo (il ProcessPoolExecutor con worker
riusati si e' rivelato fragile, per entrambi i percorsi: il pool crashava dopo
alcuni task senza una causa Python catturabile, probabile netCDF4/HDF5 su NFS
con processi forkati). I processi spawn non vedono le funzioni definite nel
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

# HDF5/netCDF4 su NFS + processi possono soffrire di problemi di locking;
# disabilita il file locking HDF5 (fix standard per NFS). Va fatto PRIMA di
# importare xarray/netCDF4.
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

import numpy as np
import pandas as pd
import xarray as xr

from config import POST_DATA, WORK_DIR
import albedo_functions as af


def global_average(dset):
    lat_dim = "lat" if "lat" in dset.dims else "latitude"
    lon_dim = "lon" if "lon" in dset.dims else "longitude"
    weights = np.cos(np.deg2rad(dset[lat_dim]))
    weights.name = "weights"
    return dset.weighted(weights).mean((lat_dim, lon_dim))


def _compute_and_save(var, era_var, exp_ctrl, exp_sens, dset_ctrl, dset_sens, era, tag, season, save_path):
    """Bias scalare per membro (mean(ctrl) - mean(obs)), globale/land/ocean.
    Nomi di output: '{exp}_{var}_{tag}_1x1_global_average_all_members[_land|_ocean]_{season}_1999.nc'."""
    outputs = {}
    for suff, mask in (("", None), ("_land", "land"), ("_ocean", "ocean")):
        if mask is None:
            c, s, e = dset_ctrl[var], dset_sens[var], era[era_var]
        else:
            c = af.land_seas_mask(dset_ctrl[var], mask)
            s = af.land_seas_mask(dset_sens[var], mask)
            e = af.land_seas_mask(era[era_var], mask)
        cg = global_average(c).to_dataset(name=var)
        sg = global_average(s).to_dataset(name=var)
        eg = global_average(e).to_dataset(name=var)

        bias_ctrl = cg.mean("time", skipna=True) - eg.mean("time", skipna=True)
        bias_sens = sg.mean("time", skipna=True) - eg.mean("time", skipna=True)

        n_members = dset_ctrl.sizes["member"]
        member_ids = xr.DataArray(np.arange(1, n_members + 1), dims="member", name="member")
        bias_ctrl["member"] = member_ids
        bias_sens["member"] = member_ids

        ctrl_out = f"{save_path}/{exp_ctrl}_{var}_{tag}_1x1_global_average_all_members{suff}_{season}_1999.nc"
        sens_out = f"{save_path}/{exp_sens}_{var}_{tag}_1x1_global_average_all_members{suff}_{season}_1999.nc"
        bias_ctrl.to_netcdf(ctrl_out)
        bias_sens.to_netcdf(sens_out)
        outputs[suff or "global"] = ctrl_out

    return outputs


def _bias_leadyears(var, era_var, exp_ctrl, exp_sens, y1, y2, season, save_path):
    """Per un range di lead year (combinazione, y1-y2): legge i file multi-anno
    ensemble e calcola il bias. Salva con tag 'leadyears_{y1}-{y2}'."""
    lead = f"{y1}-{y2}"
    lead_number = y2 - y1 + 1

    dset_ctrl = xr.open_dataset(
        POST_DATA / exp_ctrl / "1x1" / var /
        f"{exp_ctrl}_{var}_Amon_EC-Earth3_dcppA-hindcast_lead_{lead}_1x1_ensemble_m{season}_rad.nc")
    dset_sens = xr.open_dataset(
        POST_DATA / exp_sens / "1x1" / var /
        f"{exp_sens}_{var}_Amon_EC-Earth3_dcppA-hindcast_lead_{lead}_1x1_ensemble_m{season}_rad.nc")
    era = xr.open_dataset(WORK_DIR / f"ERA5_{era_var}_1x1_{lead_number}{season}.nc")
    # NON rinominare: si accede a era[era_var] col nome originale ('2t').

    # Allineamento per ANNO solare: per le finestre pari il modello etichetta la
    # media al centro (es. giugno), l'obs resta a dicembre -> 0 tempi in comune.
    dset_ctrl = dset_ctrl.assign_coords(time=pd.to_datetime(dset_ctrl["time"].values).year)
    dset_sens = dset_sens.assign_coords(time=pd.to_datetime(dset_sens["time"].values).year)
    era = era.assign_coords(time=pd.to_datetime(era["time"].values).year)
    dset_ctrl = dset_ctrl.sel(time=slice(1999, None))
    dset_sens = dset_sens.sel(time=slice(1999, None))
    era = era.sel(time=slice(1999, None))

    return _compute_and_save(var, era_var, exp_ctrl, exp_sens, dset_ctrl, dset_sens, era,
                             f"leadyears_{lead}", season, save_path)


def _bias_single_lead(var, era_var, exp_ctrl, exp_sens, lead, season, save_path):
    """Per un singolo lead time (0..4): legge i file annuali ensemble e calcola
    il bias. Salva con tag 'lead_{lead}' (stessi nomi di sempre, nessun prefisso)."""
    ctrl_path = f"{save_path}/{exp_ctrl}/{exp_ctrl}_{var}_Amon_EC-Earth3_dcppA-hindcast_*_gr_1x1_lead{lead}_rad_{season}.nc"
    dset_ctrl = xr.open_mfdataset(ctrl_path, concat_dim="member", combine="nested", chunks="auto").load()
    dset_ctrl["time"] = pd.to_datetime(dset_ctrl["time"].values).to_period("M").to_timestamp()

    sens_path = f"{save_path}/{exp_sens}/{exp_sens}_{var}_Amon_EC-Earth3_dcppA-hindcast_*_gr_1x1_lead{lead}_rad_{season}.nc"
    dset_sens = xr.open_mfdataset(sens_path, concat_dim="member", combine="nested", chunks="auto").load()
    dset_sens["time"] = pd.to_datetime(dset_sens["time"].values).to_period("M").to_timestamp()

    era = xr.open_dataset(f"{WORK_DIR}/ERA5_{era_var}_1x1_{season}.nc", chunks="auto").load()
    era = era.sel(time=slice(dset_ctrl.time[0], dset_ctrl.time[-1]))

    start_date = "1999-09-01"
    dset_ctrl = dset_ctrl.sel(time=slice(start_date, None))
    dset_sens = dset_sens.sel(time=slice(start_date, None))
    era = era.sel(time=slice(start_date, None))

    return _compute_and_save(var, era_var, exp_ctrl, exp_sens, dset_ctrl, dset_sens, era,
                             f"lead_{lead}", season, save_path)


def run_one(args):
    """Adattatore per Pool.imap_unordered (combinazioni di lead year): un solo
    argomento (tupla), ritorna una stringa di esito invece di loggare su file (i
    processi spawn non condividono l'handler di logging del notebook)."""
    var, era_var, exp_ctrl, exp_sens, y1, y2, season, save_path = args
    lead = f"{y1}-{y2}"
    try:
        outputs = _bias_leadyears(var, era_var, exp_ctrl, exp_sens, y1, y2, season, save_path)
        ok = all(os.path.exists(p) for p in outputs.values())
        return f"{season} {lead} {'ok' if ok else 'ATTENZIONE: file non trovato dopo il salvataggio'}"
    except Exception as e:
        return f"{season} {lead} ERRORE: {type(e).__name__}: {e}"


def run_one_single(args):
    """Adattatore per Pool.imap_unordered (lead singolo)."""
    var, era_var, exp_ctrl, exp_sens, lead, season, save_path = args
    try:
        outputs = _bias_single_lead(var, era_var, exp_ctrl, exp_sens, lead, season, save_path)
        ok = all(os.path.exists(p) for p in outputs.values())
        return f"{season} lead{lead} {'ok' if ok else 'ATTENZIONE: file non trovato dopo il salvataggio'}"
    except Exception as e:
        return f"{season} lead{lead} ERRORE: {type(e).__name__}: {e}"
