"""Regressione cover-temperatura, versione STAGIONALE (DJF/MAM/JJA/SON), per
combinazione di lead year. Analogo stagionale di cover_regression/cover_tas_lib.py
(annuale) - vedi quel modulo per la spiegazione dei due disegni (adattato/letterale).

Differenza principale rispetto alla versione annuale: la cover (cvh/cvl) e' una
serie MENSILE continua senza file stagionali pre-aggregati (a differenza di tas,
per cui si riusano i file DCPP gia' aggregati per stagione, come in Fig3/Fig4
stagionali). La media stagionale della cover viene quindi calcolata qui, con lo
stesso schema di af.DJF_seasonal_mean (resample QS-DEC) generalizzato alle 4
stagioni.

ATTENZIONE - rischio di disallineamento per DJF: dicembre appartiene all'anno
solare precedente rispetto a gennaio-febbraio. La convenzione di etichettatura
usata qui (QS-DEC, label = inizio periodo = dicembre) potrebbe non coincidere
con quella dei file tas DCPP gia' pre-aggregati per stagione. Un disallineamento
di un anno per DJF non darebbe errore, solo un risultato silenziosamente
sbagliato. Verificare con debug_years() prima di lanciare il calcolo completo
(vedi cella dedicata nei notebook).

Stessa esecuzione robusta (processi spawn freschi) adottata in tutta questa
sessione dopo i ripetuti crash di ProcessPoolExecutor con worker riusati.
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

from config import POST_DATA, WORK_DIR, FIG_DIR
import albedo_functions as af

LAT_MIN, LAT_MAX, LON_MIN, LON_MAX = 54, 70, 88, 110
LEADS = [(y1, y2) for y1 in range(5) for y2 in range(y1 + 1, 5)]
SEASONS = ["DJF", "MAM", "JJA", "SON"]
SEASON_START_MONTH = {"DJF": 12, "MAM": 3, "JJA": 6, "SON": 9}


def _slope_pvalue(a, b):
    """Come af.compute_slope_and_pvalue, ma senza il riferimento a 'x' non
    definito presente nell'originale (vedi cover_tas_lib.py, versione annuale)."""
    mask = ~np.isnan(a) & ~np.isnan(b)
    if np.sum(mask) < 3 or np.all(a[mask] == a[mask][0]):
        return np.nan, np.nan
    slope, intercept, r, p, std_err = stats.linregress(a[mask], b[mask])
    return slope, p


def _mask_low_variance(da, pct=10, threshold=None):
    """Maschera (NaN) i pixel dove la variabilita' temporale di da e' troppo
    bassa. Vedi cover_tas_lib.py (versione annuale) per la motivazione e per
    debug_cover_variance (diagnostica usata per scegliere threshold=1e-3,
    confermato sui dati reali: il 75* percentile di std(delta_cover) e'
    ~zero - aree non vegetate - con un salto di ~4 ordini di grandezza prima
    del segnale vero al 90* percentile).

    threshold: se fornito, ha precedenza sul percentile."""
    std = da.std("time")
    if threshold is None:
        threshold = np.nanpercentile(std.values, pct)
    return da.where(std > threshold)


def _load_tas_ensemble_season(exp, lead, season):
    """Media d'ensemble di tas per un lead-year combo, stagionale (file DCPP
    gia' pre-aggregato per stagione, stesso file usato da 04-BIAS_seasons)."""
    ds = xr.open_dataset(
        POST_DATA / exp / "1x1" / "tas" /
        f"{exp}_tas_Amon_EC-Earth3_dcppA-hindcast_lead_{lead}_1x1_ensemble_m{season}_rad.nc")
    em = ds["tas"].mean("member")
    em = em.assign_coords(time=pd.to_datetime(em["time"].values).year)
    return em


def _seasonal_cover(exp, var, season):
    """Media stagionale della cover a partire dalla serie mensile continua.
    Stesso schema di af.DJF_seasonal_mean (resample QS-DEC), generalizzato
    alle 4 stagioni. 'time' diventa l'anno solare del PRIMO mese della stagione
    (dicembre per DJF: la label precede gennaio-febbraio dello stesso inverno)."""
    if exp == "a1ua":
        fname = f"{exp}_effective_{var}_1x1.nc"
    else:
        fname = f"{exp}_effective_{var}_199311-201910_1x1.nc"
    da = xr.open_dataset(POST_DATA / fname)[var]
    resampled = da.resample(time="QS-DEC").mean("time")
    start_month = SEASON_START_MONTH[season]
    seasonal = resampled.sel(time=resampled["time"].dt.month == start_month)
    seasonal = seasonal.assign_coords(time=pd.to_datetime(seasonal["time"].values).year)
    return seasonal


def _era5_obs_season(era_var, lead_number, season):
    """ERA5 stagionale (stesso file usato da 04-BIAS_seasons/06b)."""
    obs = xr.open_dataset(WORK_DIR / f"ERA5_{era_var}_1x1_{lead_number}{season}.nc")
    obs = obs.rename({era_var: "tas"})["tas"].assign_coords(
        time=pd.to_datetime(obs["time"].values).year)
    return obs


def debug_years(exp_ctrl, exp_sens, var, y1, y2, season):
    """Da chiamare PRIMA del calcolo completo: stampa gli anni delle due serie
    (cover stagionale e tas stagionale) per un singolo caso, cosi' si puo'
    verificare a occhio se la coppia rappresenta davvero la stessa stagione
    (rischio di disallineamento di un anno, specifico a DJF - vedi docstring
    del modulo)."""
    lead = f"{y1}-{y2}"
    tas = _load_tas_ensemble_season(exp_ctrl, lead, season)
    cov = _seasonal_cover(exp_ctrl, var, season)
    print(f"tas  {season} {lead}: anni = {sorted(tas['time'].values.tolist())}")
    print(f"cover {season}     : anni = {sorted(cov['time'].values.tolist())}")
    common = sorted(set(tas["time"].values.tolist()) & set(cov["time"].values.tolist()))
    print(f"anni in comune: {common} ({len(common)} totali)")
    if season == "DJF":
        print("Verifica manuale: l'anno X di 'cover DJF' e l'anno X di 'tas DJF' "
              "devono corrispondere allo STESSO inverno (dicembre anno X + "
              "gennaio/febbraio anno X+1). Se il progetto etichetta DJF con "
              "l'anno di GENNAIO invece che di DICEMBRE, qui c'e' uno sfasamento "
              "di un anno da correggere (vedi SEASON_START_MONTH nel modulo).")


def _scatter_plot(box_x, box_y, y_pred, p, r, title, xlabel, ylabel):
    """Vedi cover_tas_lib.py (versione annuale) per la motivazione: af.lr_plot
    ha assi/limiti hardcoded per il caso albedo, inadatti alla cover."""
    fig, ax = plt.subplots(figsize=[10, 8])
    ax.scatter(box_x, box_y, label="Dati", c="blue", alpha=0.7)
    ax.plot(box_x, y_pred, color="red", linewidth=2,
            label=f"Retta di regressione (p={p:.2f}, r={r:.2f})")

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
    slope_map, p_map = xr.apply_ufunc(
        _slope_pvalue, delta_x, delta_y,
        input_core_dims=[["time"], ["time"]],
        vectorize=True, output_dtypes=[float, float], output_core_dims=[[], []],
    )
    af.map_plot(slope_map, p_map, levels=[-2, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 2],
               title=title, cmap="bwr", sign=0.90)
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


def run_one_adapted_season(args):
    """Notebook 01 stagionale: delta_tas (SENS-CTRL) vs delta_cover (SENS-CTRL)."""
    exp_ctrl, exp_sens, var, season, y1, y2, save_path = args
    lead = f"{y1}-{y2}"
    try:
        tas_ctrl = _load_tas_ensemble_season(exp_ctrl, lead, season)
        tas_sens = _load_tas_ensemble_season(exp_sens, lead, season)
        anom_ctrl = tas_ctrl - tas_ctrl.mean("time")
        anom_sens = tas_sens - tas_sens.mean("time")
        delta_tas = anom_sens - anom_ctrl

        cov_ctrl = _seasonal_cover(exp_ctrl, var, season)
        cov_sens = _seasonal_cover(exp_sens, var, season)
        cov_ctrl, cov_sens = xr.align(cov_ctrl, cov_sens, join="inner")
        cov_anom_ctrl = cov_ctrl - cov_ctrl.mean("time")
        cov_anom_sens = cov_sens - cov_sens.mean("time")
        delta_cover = cov_anom_sens - cov_anom_ctrl
        delta_cover = _mask_low_variance(delta_cover, threshold=1e-3)  # vedi versione annuale

        delta_tas, delta_cover = xr.align(delta_tas, delta_cover, join="inner")
        if delta_tas.sizes.get("time", 0) < 3:
            return f"{season} {var} {lead} SALTATO: solo {delta_tas.sizes.get('time', 0)} anni in comune"

        title = f"delta_tas_vs_delta_{var}_{season}_{lead}"
        _map_and_box_regression(
            delta_cover, delta_tas, title,
            f"{save_path}/{title}_map.png", f"{save_path}/{title}_scatter.png",
            f"{POST_DATA}/adapted_regression_{var}_{season}_{LAT_MIN}_{LAT_MAX}_{LON_MIN}_{LON_MAX}_{lead}.nc",
            xlabel=f"delta {var} (SENS-CTRL)", ylabel="delta tas (SENS-CTRL, K)",
        )
        return f"{season} {var} {lead} ok ({delta_tas.sizes['time']} anni)"
    except Exception as e:
        return f"{season} {var} {lead} ERRORE: {type(e).__name__}: {e}"


def run_one_literal_season(args):
    """Notebook 02 stagionale: replica skill-vs-obs (obs cover = SENS)."""
    exp_ctrl, exp_sens, var, era_var, season, y1, y2, save_path = args
    lead = f"{y1}-{y2}"
    lead_number = y2 - y1 + 1
    try:
        tas_ctrl = _load_tas_ensemble_season(exp_ctrl, lead, season)
        tas_sens = _load_tas_ensemble_season(exp_sens, lead, season)
        obs_tas = _era5_obs_season(era_var, lead_number, season)
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

        cov_ctrl = _seasonal_cover(exp_ctrl, var, season)
        cov_sens = _seasonal_cover(exp_sens, var, season)  # anche usato come "obs"
        cov_ctrl, cov_obs = xr.align(cov_ctrl, cov_sens, join="inner")
        anom_cov_ctrl = cov_ctrl - cov_ctrl.mean("time")
        anom_cov_obs = cov_obs - cov_obs.mean("time")
        skill_cov_ctrl = (anom_cov_ctrl * anom_cov_obs) / (cov_ctrl.std("time") * cov_obs.std("time"))
        skill_cov_sens = (anom_cov_obs * anom_cov_obs) / (cov_obs.std("time") * cov_obs.std("time"))
        skill_cov_ctrl, skill_cov_sens = xr.align(skill_cov_ctrl, skill_cov_sens, join="inner")
        delta_cover = skill_cov_sens - skill_cov_ctrl

        delta_tas, delta_cover = xr.align(delta_tas, delta_cover, join="inner")
        if delta_tas.sizes.get("time", 0) < 3:
            return f"{season} {var} {lead} SALTATO: solo {delta_tas.sizes.get('time', 0)} anni in comune"

        title = f"literal_delta_tas_vs_delta_{var}_{season}_{lead}"
        _map_and_box_regression(
            delta_cover, delta_tas, title,
            f"{save_path}/{title}_map.png", f"{save_path}/{title}_scatter.png",
            f"{POST_DATA}/literal_regression_{var}_{season}_{LAT_MIN}_{LAT_MAX}_{LON_MIN}_{LON_MAX}_{lead}.nc",
            xlabel=f"delta skill {var} (SENS-CTRL)", ylabel="delta skill tas (SENS-CTRL)",
        )
        return f"{season} {var} {lead} ok ({delta_tas.sizes['time']} anni)"
    except Exception as e:
        return f"{season} {var} {lead} ERRORE: {type(e).__name__}: {e}"
