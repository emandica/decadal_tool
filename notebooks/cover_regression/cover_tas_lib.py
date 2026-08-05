"""Regressione tra vegetation cover (cvh/cvl) e temperatura, per combinazione di
lead year (annuale, come 07-covariance_regression.ipynb).

Due disegni distinti (vedi le due funzioni run_one_adapted / run_one_literal):

ADATTATO (01): delta_tas (SENS-CTRL) vs delta_cover (SENS-CTRL). Non richiede
osservazioni indipendenti di cover: dato che SENS e' forzato con le osservazioni
di vegetation cover (per costruzione cover_SENS ~ obs), delta_cover = cover_SENS
- cover_CTRL rappresenta gia' la correzione rispetto al modello di vegetazione
dinamico di CTRL. Analisi non circolare.

LETTERALE (02): replica esatta dello schema 06/06b/07 (skill-vs-obs per SENS e
CTRL, poi delta), usando cover_SENS come proxy delle osservazioni. Per SENS lo
skill sara' quasi banale (SENS ~ obs per costruzione): incluso comunque su
richiesta esplicita, come termine di confronto con l'analisi adattata.

Le funzioni di lettura dati (tas, cover) sono CONDIVISE tra i due disegni.

Compatibilita' cover/tas: i file di cover (cvh/cvl) sono una serie temporale
CONTINUA (1993-2019, nessuna dimensione member, non organizzata per lead year),
mentre tas e' un hindcast DCPP (member x lead-year, 'time' = un valore per anno
di inizializzazione). Allineamento per ANNO SOLARE (stessa tecnica gia'
validata per il bias stagionale): si converte 'time' in anno e si intersecano i
due dataset sugli anni in comune.

Esecuzione in processi 'spawn' freschi (multiprocessing.Pool con
maxtasksperchild=1): lezione appresa dal debug del bias stagionale in questa
stessa sessione (ProcessPoolExecutor con worker riusati e' fragile, netCDF4/
HDF5 su NFS).
"""
import os
import sys

_cfg = os.path.dirname(os.path.abspath(__file__))
while _cfg != os.path.dirname(_cfg):
    if os.path.exists(os.path.join(_cfg, "config.py")):
        if _cfg not in sys.path:
            sys.path.insert(0, _cfg)
        break
    _cfg = os.path.dirname(_cfg)

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

from config import POST_DATA, WORK_DIR, ERA5_ROOT, FIG_DIR
import albedo_functions as af

# Box di riferimento (Siberia), stessi valori di 07-covariance_regression.ipynb
# e 08-covariance_plot.ipynb (analisi albedo-temperatura).
LAT_MIN, LAT_MAX, LON_MIN, LON_MAX = 54, 70, 88, 110

LEADS = [(y1, y2) for y1 in range(5) for y2 in range(y1 + 1, 5)]


def _slope_pvalue(a, b):
    """Come af.compute_slope_and_pvalue, ma senza il riferimento a 'x' non
    definito presente nell'originale (bug preesistente in albedo_functions.py,
    non toccato qui perche' condiviso con altro codice)."""
    mask = ~np.isnan(a) & ~np.isnan(b)
    if np.sum(mask) < 3 or np.all(a[mask] == a[mask][0]):
        return np.nan, np.nan
    slope, intercept, r, p, std_err = stats.linregress(a[mask], b[mask])
    return slope, p


def _load_tas_ensemble(exp, lead):
    """Media d'ensemble di tas per un lead-year combo (annuale), con 'time'
    convertito in anno solare. Stesso file gia' usato da 04-BIAS/06b."""
    ds = xr.open_dataset(
        POST_DATA / exp / "1x1" / "tas" /
        f"{exp}_tas_Amon_EC-Earth3_dcppA-hindcast_lead_{lead}_1x1_ensemble_rad.nc")
    em = ds["tas"].mean("member")
    em = em.assign_coords(time=pd.to_datetime(em["time"].values).year)
    return em


def _load_cover(exp, var):
    """Serie temporale continua di vegetation cover (cvh o cvl), 'time' in anno
    solare. Nomi file identici a quelli gia' usati in Fig1/Fig2 (final_figure)."""
    if exp == "a1ua":
        fname = f"{exp}_effective_{var}_1x1.nc"
    else:
        fname = f"{exp}_effective_{var}_199311-201910_1x1.nc"
    da = xr.open_dataset(POST_DATA / fname)[var]
    da = da.assign_coords(time=pd.to_datetime(da["time"].values).year)
    # media annuale: la serie e' mensile, il confronto con tas e' per anno
    da = da.groupby("time").mean("time")
    return da


def _era5_obs(era_var, lead_number):
    """ERA5 alla risoluzione del lead-year combo (stesso file usato da 06b/04-BIAS)."""
    obs = xr.open_dataset(WORK_DIR / f"ERA5_{era_var}_1x1_{lead_number}year.nc")
    obs = obs.rename({era_var: "tas"})["tas"].assign_coords(
        time=pd.to_datetime(obs["time"].values).year)
    return obs


def _scatter_plot(box_x, box_y, y_pred, p, r, title, xlabel, ylabel, rho=None, p_spearman=None):
    """Come af.lr_plot, ma senza i due difetti che la rendono inadatta alla
    cover (i cui delta sono ~100x piu' piccoli di quelli dell'albedo):
    (1) limiti degli assi FISSI a [-1,2] (qui: automatici, con margine);
    (2) etichette 'time' convertite in anno (int) passate a pd.to_datetime,
    che le interpreta come nanosecondi dall'epoca Unix -> tutte '1970'
    (qui: l'anno e' gia' un intero, si annota direttamente).

    Se rho/p_spearman sono forniti, li aggiunge in legenda accanto a
    slope/r/p (Pearson): la relazione potrebbe essere monotona ma non
    lineare, e Spearman lo cattura mentre linregress no."""
    fig, ax = plt.subplots(figsize=[10, 8])
    ax.scatter(box_x, box_y, label="Dati", c="blue", alpha=0.7)
    label = f"Pearson: p={p:.2f}, r={r:.2f}"
    if rho is not None:
        label += f"\nSpearman: rho={rho:.2f}, p={p_spearman:.2f}"
    ax.plot(box_x, y_pred, color="red", linewidth=2, label=label)

    if hasattr(box_x, "time"):
        for xi, yi, year in zip(box_x.values, box_y.values, box_x["time"].values):
            ax.annotate(int(year), (xi, yi), fontsize=10, color="black")

    xpad = 0.1 * (box_x.values.max() - box_x.values.min() or 1)
    ypad = 0.1 * (box_y.values.max() - box_y.values.min() or 1)
    ax.set_xlim(box_x.values.min() - xpad, box_x.values.max() + xpad)
    ax.set_ylim(box_y.values.min() - ypad, box_y.values.max() + ypad)

    ax.set_title(title, fontsize=18)
    ax.set_xlabel(xlabel, fontsize=16)
    ax.set_ylabel(ylabel, fontsize=16)
    ax.tick_params(axis="both", labelsize=12)
    ax.grid(True, linestyle="--", alpha=0.5)
    ax.legend(fontsize=12, loc="best")
    plt.tight_layout()


def _map_and_box_regression(delta_x, delta_y, title, png_map, png_scatter, nc_out,
                            xlabel="delta cover", ylabel="delta tas"):
    """Regressione per pixel (mappa, af.map_plot) + regressione scalare sul box
    Siberia (scatter, _scatter_plot). delta_x, delta_y: (time=anno, lat, lon),
    gia' allineati sugli stessi anni."""
    slope_map, p_map = xr.apply_ufunc(
        _slope_pvalue, delta_x, delta_y,
        input_core_dims=[["time"], ["time"]],
        vectorize=True, output_dtypes=[float, float], output_core_dims=[[], []],
    )
    af.map_plot(slope_map, p_map, levels=[-2, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 2],
               title=title, cmap="bwr")
    plt.savefig(png_map, dpi=300, bbox_inches="tight")
    plt.close("all")

    box_x = af.domain_selection(delta_x, LAT_MIN, LAT_MAX, LON_MIN, LON_MAX)
    box_y = af.domain_selection(delta_y, LAT_MIN, LAT_MAX, LON_MIN, LON_MAX)
    box_x, box_y = xr.align(box_x, box_y, join="inner")
    slope, intercept, r, p, std_err = stats.linregress(box_x.values, box_y.values)
    y_pred = slope * box_x.values + intercept

    _scatter_plot(box_x, box_y, y_pred, p, r, title, xlabel, ylabel)
    plt.savefig(png_scatter, dpi=300, bbox_inches="tight")
    plt.close("all")

    xr.Dataset({
        "delta_x": box_x, "delta_y": box_y,
        "slope": xr.DataArray(slope), "intercept": xr.DataArray(intercept),
        "r": xr.DataArray(r), "p": xr.DataArray(p), "std_err": xr.DataArray(std_err),
        "y_pred": xr.DataArray(y_pred, dims=["time"], coords={"time": box_x.time}),
    }).to_netcdf(nc_out)


def _slope_pvalue_full(a, b):
    """Come _slope_pvalue, ma calcola ANCHE la correlazione di rango di
    Spearman accanto alla regressione lineare di Pearson: la relazione tra
    delta_cover e delta_skill_tas potrebbe essere monotona ma non lineare,
    e Spearman la cattura mentre scipy.stats.linregress no."""
    mask = ~np.isnan(a) & ~np.isnan(b)
    if np.sum(mask) < 3 or np.all(a[mask] == a[mask][0]):
        return np.nan, np.nan, np.nan, np.nan
    slope, intercept, r, p, std_err = stats.linregress(a[mask], b[mask])
    rho, p_spear = stats.spearmanr(a[mask], b[mask])
    return slope, p, rho, p_spear


def _map_and_box_regression_hybrid(delta_x, delta_y, title, png_base, nc_out, xlabel, ylabel):
    """Come _map_and_box_regression, ma con DUE mappe (Pearson slope e Spearman
    rho, ciascuna con la propria significativita') e lo scatter con entrambe le
    statistiche in legenda. Usata solo da run_one_hybrid (notebook 03)."""
    slope_map, p_map, rho_map, p_spear_map = xr.apply_ufunc(
        _slope_pvalue_full, delta_x, delta_y,
        input_core_dims=[["time"], ["time"]],
        vectorize=True, output_dtypes=[float, float, float, float],
        output_core_dims=[[], [], [], []],
    )

    af.map_plot(slope_map, p_map, levels=[-2, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 2],
               title=f"{title} (Pearson, slope)", cmap="bwr")
    plt.savefig(f"{png_base}_map_pearson.png", dpi=300, bbox_inches="tight")
    plt.close("all")

    af.map_plot(rho_map, p_spear_map, levels=[-1, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 1],
               title=f"{title} (Spearman, rho)", cmap="PuOr")
    plt.savefig(f"{png_base}_map_spearman.png", dpi=300, bbox_inches="tight")
    plt.close("all")

    box_x = af.domain_selection(delta_x, LAT_MIN, LAT_MAX, LON_MIN, LON_MAX)
    box_y = af.domain_selection(delta_y, LAT_MIN, LAT_MAX, LON_MIN, LON_MAX)
    box_x, box_y = xr.align(box_x, box_y, join="inner")
    slope, intercept, r, p, std_err = stats.linregress(box_x.values, box_y.values)
    rho, p_spear = stats.spearmanr(box_x.values, box_y.values)
    y_pred = slope * box_x.values + intercept

    _scatter_plot(box_x, box_y, y_pred, p, r, title, xlabel, ylabel, rho=rho, p_spearman=p_spear)
    plt.savefig(f"{png_base}_scatter.png", dpi=300, bbox_inches="tight")
    plt.close("all")

    xr.Dataset({
        "delta_x": box_x, "delta_y": box_y,
        "slope": xr.DataArray(slope), "intercept": xr.DataArray(intercept),
        "r": xr.DataArray(r), "p": xr.DataArray(p), "std_err": xr.DataArray(std_err),
        "rho_spearman": xr.DataArray(rho), "p_spearman": xr.DataArray(p_spear),
        "y_pred": xr.DataArray(y_pred, dims=["time"], coords={"time": box_x.time}),
    }).to_netcdf(nc_out)


def run_one_hybrid(args):
    """Notebook 3 (disegno ibrido): prende il pezzo valido di ciascun metodo.

    X = delta_cover (come Notebook 1/adattato: cover_SENS - cover_CTRL). Non
    circolare: SENS e' forzato con le osservazioni, quindi delta_cover misura
    genuinamente l'entita' della correzione rispetto al modello di vegetazione
    dinamico di CTRL.

    Y = delta_skill_tas (come Notebook 2/letterale, ma SOLO il lato tas, dove
    ERA5 e' un riferimento osservativo vero e indipendente). E' un genuino
    miglioramento di skill (non solo una differenza di anomalie come nel
    Notebook 1), e non soffre della circolarita' che affligge lo skill della
    cover nel Notebook 2 (qui la cover non viene mai confrontata con "se
    stessa": si usa solo la sua differenza SENS-CTRL).
    """
    exp_ctrl, exp_sens, var, era_var, y1, y2, save_path = args
    lead = f"{y1}-{y2}"
    lead_number = y2 - y1 + 1
    try:
        # --- Y: delta skill tas (vs ERA5, genuino miglioramento di skill) ---
        tas_ctrl = _load_tas_ensemble(exp_ctrl, lead)
        tas_sens = _load_tas_ensemble(exp_sens, lead)
        obs_tas = _era5_obs(era_var, lead_number)
        tas_ctrl, obs_tas_c = xr.align(tas_ctrl, obs_tas, join="inner")
        tas_sens, obs_tas_s = xr.align(tas_sens, obs_tas, join="inner")

        anom_ctrl = tas_ctrl - tas_ctrl.mean("time")
        anom_sens = tas_sens - tas_sens.mean("time")
        anom_obs_c = obs_tas_c - obs_tas_c.mean("time")
        anom_obs_s = obs_tas_s - obs_tas_s.mean("time")
        skill_tas_ctrl = (anom_ctrl * anom_obs_c) / (tas_ctrl.std("time") * obs_tas_c.std("time"))
        skill_tas_sens = (anom_sens * anom_obs_s) / (tas_sens.std("time") * obs_tas_s.std("time"))
        skill_tas_ctrl, skill_tas_sens = xr.align(skill_tas_ctrl, skill_tas_sens, join="inner")
        delta_skill_tas = skill_tas_sens - skill_tas_ctrl

        # --- X: delta cover (non circolare, come Notebook 1) ---
        cov_ctrl = _load_cover(exp_ctrl, var)
        cov_sens = _load_cover(exp_sens, var)
        cov_ctrl, cov_sens = xr.align(cov_ctrl, cov_sens, join="inner")
        cov_anom_ctrl = cov_ctrl - cov_ctrl.mean("time")
        cov_anom_sens = cov_sens - cov_sens.mean("time")
        delta_cover = cov_anom_sens - cov_anom_ctrl

        delta_skill_tas, delta_cover = xr.align(delta_skill_tas, delta_cover, join="inner")
        if delta_skill_tas.sizes.get("time", 0) < 3:
            return f"{var} {lead} SALTATO: solo {delta_skill_tas.sizes.get('time', 0)} anni in comune"

        title = f"hybrid_delta_skill_tas_vs_delta_{var}_{lead}"
        _map_and_box_regression_hybrid(
            delta_cover, delta_skill_tas, title,
            f"{save_path}/{title}",
            f"{POST_DATA}/hybrid_regression_{var}_{LAT_MIN}_{LAT_MAX}_{LON_MIN}_{LON_MAX}_{lead}.nc",
            xlabel=f"delta {var} (SENS-CTRL)", ylabel="delta skill tas (SENS-CTRL, vs ERA5)",
        )
        return f"{var} {lead} ok ({delta_skill_tas.sizes['time']} anni)"
    except Exception as e:
        return f"{var} {lead} ERRORE: {type(e).__name__}: {e}"


def run_one_adapted(args):
    """Notebook 1: delta_tas (SENS-CTRL) vs delta_cover (SENS-CTRL)."""
    exp_ctrl, exp_sens, var, y1, y2, save_path = args
    lead = f"{y1}-{y2}"
    try:
        tas_ctrl = _load_tas_ensemble(exp_ctrl, lead)
        tas_sens = _load_tas_ensemble(exp_sens, lead)
        anom_ctrl = tas_ctrl - tas_ctrl.mean("time")
        anom_sens = tas_sens - tas_sens.mean("time")
        delta_tas = anom_sens - anom_ctrl  # auto-allineato sugli anni in comune

        cov_ctrl = _load_cover(exp_ctrl, var)
        cov_sens = _load_cover(exp_sens, var)
        cov_ctrl, cov_sens = xr.align(cov_ctrl, cov_sens, join="inner")
        cov_anom_ctrl = cov_ctrl - cov_ctrl.mean("time")
        cov_anom_sens = cov_sens - cov_sens.mean("time")
        delta_cover = cov_anom_sens - cov_anom_ctrl

        delta_tas, delta_cover = xr.align(delta_tas, delta_cover, join="inner")
        if delta_tas.sizes.get("time", 0) < 3:
            return f"{var} {lead} SALTATO: solo {delta_tas.sizes.get('time', 0)} anni in comune"

        title = f"delta_tas_vs_delta_{var}_{lead}"
        _map_and_box_regression(
            delta_cover, delta_tas, title,
            f"{save_path}/{title}_map.png", f"{save_path}/{title}_scatter.png",
            f"{POST_DATA}/adapted_regression_{var}_{LAT_MIN}_{LAT_MAX}_{LON_MIN}_{LON_MAX}_{lead}.nc",
            xlabel=f"delta {var} (SENS-CTRL)", ylabel="delta tas (SENS-CTRL, K)",
        )
        return f"{var} {lead} ok ({delta_tas.sizes['time']} anni)"
    except Exception as e:
        return f"{var} {lead} ERRORE: {type(e).__name__}: {e}"


def run_one_literal(args):
    """Notebook 2: replica 06/06b/07, skill-vs-obs (obs cover = SENS)."""
    exp_ctrl, exp_sens, var, era_var, y1, y2, save_path = args
    lead = f"{y1}-{y2}"
    lead_number = y2 - y1 + 1
    try:
        # --- skill-vs-obs per la temperatura (obs = ERA5), come 06b ---
        tas_ctrl = _load_tas_ensemble(exp_ctrl, lead)
        tas_sens = _load_tas_ensemble(exp_sens, lead)
        obs_tas = _era5_obs(era_var, lead_number)
        tas_ctrl, obs_tas_c = xr.align(tas_ctrl, obs_tas, join="inner")
        tas_sens, obs_tas_s = xr.align(tas_sens, obs_tas, join="inner")

        anom_ctrl = tas_ctrl - tas_ctrl.mean("time")
        anom_sens = tas_sens - tas_sens.mean("time")
        anom_obs_c = obs_tas_c - obs_tas_c.mean("time")
        anom_obs_s = obs_tas_s - obs_tas_s.mean("time")
        skill_tas_ctrl = (anom_ctrl * anom_obs_c) / (tas_ctrl.std("time") * obs_tas_c.std("time"))
        skill_tas_sens = (anom_sens * anom_obs_s) / (tas_sens.std("time") * obs_tas_s.std("time"))
        skill_tas_ctrl, skill_tas_sens = xr.align(skill_tas_ctrl, skill_tas_sens, join="inner")
        delta_tas = skill_tas_sens - skill_tas_ctrl

        # --- skill-vs-obs per la cover (obs proxy = cover_SENS) ---
        cov_ctrl = _load_cover(exp_ctrl, var)
        cov_sens = _load_cover(exp_sens, var)  # anche usato come "obs"
        cov_ctrl, cov_obs = xr.align(cov_ctrl, cov_sens, join="inner")
        anom_cov_ctrl = cov_ctrl - cov_ctrl.mean("time")
        anom_cov_obs = cov_obs - cov_obs.mean("time")
        skill_cov_ctrl = (anom_cov_ctrl * anom_cov_obs) / (cov_ctrl.std("time") * cov_obs.std("time"))
        # SENS vs se stesso: (anomalia_obs/std_obs)^2, non genuino skill previsivo
        # (misura solo la struttura di varianza delle obs, non l'accordo con un
        # dato indipendente) - circolare per costruzione, non banalmente = 1.
        skill_cov_sens = (anom_cov_obs * anom_cov_obs) / (cov_obs.std("time") * cov_obs.std("time"))
        skill_cov_ctrl, skill_cov_sens = xr.align(skill_cov_ctrl, skill_cov_sens, join="inner")
        delta_cover = skill_cov_sens - skill_cov_ctrl

        delta_tas, delta_cover = xr.align(delta_tas, delta_cover, join="inner")
        if delta_tas.sizes.get("time", 0) < 3:
            return f"{var} {lead} SALTATO: solo {delta_tas.sizes.get('time', 0)} anni in comune"

        title = f"literal_delta_tas_vs_delta_{var}_{lead}"
        _map_and_box_regression(
            delta_cover, delta_tas, title,
            f"{save_path}/{title}_map.png", f"{save_path}/{title}_scatter.png",
            f"{POST_DATA}/literal_regression_{var}_{LAT_MIN}_{LAT_MAX}_{LON_MIN}_{LON_MAX}_{lead}.nc",
            xlabel=f"delta skill {var} (SENS-CTRL)", ylabel="delta skill tas (SENS-CTRL)",
        )
        return f"{var} {lead} ok ({delta_tas.sizes['time']} anni)"
    except Exception as e:
        return f"{var} {lead} ERRORE: {type(e).__name__}: {e}"
