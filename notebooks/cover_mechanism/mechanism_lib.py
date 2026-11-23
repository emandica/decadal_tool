"""Passo 1 dell'analisi vegetazione -> albedo -> temperatura: differenze di
STATO MEDIO stagionale SENS - CTRL (il meccanismo, senza skill e senza
osservazioni).

Disegno degli esperimenti: in SENS (a52o) LAI, cover e tipi di vegetazione
sono osservati e prescritti come boundary condition, anche negli anni
previsti; in CTRL (a1ua) la vegetazione viene da un run storico di
LPJ-GUESS. La differenza di vegetazione SENS - CTRL e' quindi soprattutto una
differenza di stato medio (i bias di LPJ-GUESS) e di trend: le anomalie
interannuali di CTRL non sono sincronizzate con quelle osservate. Le
regressioni di cover_regression/ tolgono la media temporale e lavorano solo
sulle anomalie, cioe' sulla parte piu' piccola. Qui si guarda la parte piu'
grande, lo stato medio, stagione per stagione e lead per lead, e la si mette
in relazione con le differenze di albedo, neve e temperatura del modello:

  - DJF/MAM: piu' vegetazione alta copre la neve -> albedo piu' basso -> piu' caldo;
  - JJA: piu' LAI -> piu' evapotraspirazione -> piu' freddo (senza i flussi,
    qui si vede solo l'effetto finale su tas).

Calcolo (notebook 01-mean_state_delta):
  run_one_model  per (variabile del modello, stagione, lead singolo): media
                 d'ensemble per anno di SENS e CTRL, differenza per anno d(t),
                 media nel tempo e p-value col bootstrap a blocchi del repo
                 (af.block_bootstrap_delta_pvalue, statistic='bias': anni a
                 blocchi di 5, gli stessi per i due esperimenti, membri
                 ricampionati separatamente - come le mappe di bias, Fig. 4).
  run_one_cover  per (variabile di vegetazione, stagione): serie stagionali per
                 anno di SENS e CTRL dalla serie mensile continua, con la stessa
                 etichetta dei file del modello (anno del primo mese della
                 stagione: dicembre per DJF, come il resample 'QS-DEC' di
                 03-decadal_season_aggregation).

Figure (notebook 01-Fig_mean_state_delta): mappe delle differenze, relazioni
tra pixel per classi di delta x (con intervalli bootstrap), coerenza tra lead
nel box Siberia e nelle terre boreali.

Le funzioni di calcolo girano in processi 'spawn' freschi
(multiprocessing.Pool con maxtasksperchild=1), come nel resto del repo: per
questo stanno in un modulo importabile e non nel notebook.
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

# HDF5/netCDF4 su NFS + processi: disabilita il file locking (prima di xarray)
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

import re
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt   # niente matplotlib.use(): nel notebook resta il backend inline; le figure si salvano e chiudono
import matplotlib.image as mpimg
from scipy import stats

from config import BC_DATA, CONFESS_DATA, POST_DATA, WORK_DIR
import albedo_functions as af

SEASONS = ["DJF", "MAM", "JJA", "SON"]
LEADS_SINGLE = [0, 1, 2, 3, 4]
START = "1999-09-01"                 # come il bias stagionale per lead singolo
MODEL_VARS = ["tas", "alb", "snd"]
COVER_VARS = ["cvh", "cvl", "lai_hv", "lai_lv"]
COVER_KINDS = set(COVER_VARS) | {"cvt"}        # cvt = cvh + cvl
BOX = (54, 70, 88, 110)              # box Siberia (LAT_MIN, LAT_MAX, LON_MIN, LON_MAX), come cover_regression
SEASON_START_MONTH = {"DJF": 12, "MAM": 3, "JJA": 6, "SON": 9}
UNITS = {"tas": "K", "alb": "-", "snd": "m", "cvh": "frazione", "cvl": "frazione",
         "cvt": "frazione", "lai_hv": "m2/m2", "lai_lv": "m2/m2"}
SEASON_COLORS = {"DJF": "tab:blue", "MAM": "tab:green", "JJA": "tab:red", "SON": "tab:orange"}
STRATA_COLORS = ["tab:blue", "tab:orange", "tab:green", "tab:purple"]


# ---------------------------------------------------------------------------
# Griglia e tempo
# ---------------------------------------------------------------------------
def _fix_grid(da):
    """Nomi lat/lon, longitudini in [0, 360), coordinate arrotondate e ordinate.

    Se le longitudini sono centrate a mezzo grado (visto per il file 'snd' di
    CTRL: -179.5..179.5 invece di 0..359), le riporta sui gradi interi come
    fa cover_tas_lib._load_snow_ensemble (errore di posizione max 0.5 gradi).
    Le latitudini non si toccano: se fossero diverse, _common_grid lo segnala."""
    ren = {k: v for k, v in (("latitude", "lat"), ("longitude", "lon")) if k in da.dims}
    if ren:
        da = da.rename(ren)
    lon = np.asarray(da["lon"].values, dtype=float)
    frac = lon % 1
    if not np.allclose(frac, 0, atol=1e-6) and np.allclose(frac, 0.5, atol=1e-6):
        lon = lon - 0.5
    da = da.assign_coords(lon=np.round(lon % 360, 4),
                          lat=np.round(np.asarray(da["lat"].values, dtype=float), 4))
    return da.sortby("lat").sortby("lon")


def _common_grid(a, b, what):
    """Restringe a e b ai punti di griglia in comune. Errore esplicito se i
    punti in comune sono meno del 90% (griglie davvero diverse): un
    disallineamento silenzioso darebbe relazioni vuote o sbagliate."""
    lat = np.intersect1d(a["lat"].values, b["lat"].values)
    lon = np.intersect1d(a["lon"].values, b["lon"].values)
    frac = min(lat.size / min(a.sizes["lat"], b.sizes["lat"]),
               lon.size / min(a.sizes["lon"], b.sizes["lon"]))
    if frac < 0.9:
        def rng(x, d):
            v = x[d].values
            return f"{d} {v.min():g}..{v.max():g} ({v.size})"
        raise ValueError(f"{what}: griglie diverse - {rng(a, 'lat')}, {rng(a, 'lon')} contro "
                         f"{rng(b, 'lat')}, {rng(b, 'lon')}; punti in comune {lat.size}x{lon.size}")
    return a.sel(lat=lat, lon=lon), b.sel(lat=lat, lon=lon)


def _to_years(da):
    """'time' -> anno (intero). Se due tempi cadono nello stesso anno (visto
    sul cluster per un file 'snd'), li media come fa cover_tas_lib."""
    years = da["time"].dt.year.values
    da = da.assign_coords(time=years)
    if np.unique(years).size < years.size:
        da = da.groupby("time").mean("time")
    return da


# ---------------------------------------------------------------------------
# Variabili del modello (tas, alb, snd): file stagionali per lead singolo
# ---------------------------------------------------------------------------
def model_files(exp, var, season, lead):
    """File per membro prodotti da 03-decadal_season_aggregation
    (create_lead_time_series): WORK_DIR/<exp>/<exp>_<var>_<realm>_..._<membro>_gr_1x1_lead<k>_rad_<stagione>.nc.
    Il realm (Amon, LImon...) e' lasciato libero."""
    pattern = f"{exp}_{var}_*_EC-Earth3_dcppA-hindcast_*_gr_1x1_lead{lead}_rad_{season}.nc"
    return sorted((Path(WORK_DIR) / exp).glob(pattern))


def load_model_season(exp, var, season, lead, start=START):
    """(member, time=anno, lat, lon) per un lead singolo e una stagione.
    Ogni 'time' e' un anno di inizializzazione; l'anno e' quello del primo
    mese della stagione (dicembre per DJF)."""
    files = model_files(exp, var, season, lead)
    if not files:
        raise FileNotFoundError(f"nessun file {exp} {var} {season} lead{lead} in {Path(WORK_DIR) / exp}")
    ds = xr.open_mfdataset(files, concat_dim="member", combine="nested", chunks="auto")
    da = ds[var].load()
    ds.close()
    da = da.sel(time=slice(start, None))
    da = _to_years(da)
    da = da.where(np.isfinite(da))   # l'albedo da rsus/rsds puo' avere inf con rsds = 0
    return _fix_grid(da).transpose("member", "time", "lat", "lon")


def compute_model_delta(exp_ctrl, exp_sens, var, season, lead, out_dir, n_iter=1000, block=5, start=START):
    """Differenza di stato medio SENS - CTRL per (var, stagione, lead) e suo
    p-value; salva <out_dir>/delta_<var>_<stagione>_lead<k>_1999.nc con:
      d          (time, lat, lon)  media d'ensemble SENS - CTRL per anno
      delta      (lat, lon)        media nel tempo di d
      p          (lat, lon)        p-value bilatero (bootstrap a blocchi del repo)
      clim_ctrl, clim_sens (lat, lon)  climatologie (servono per la maschera di neve)."""
    ctrl = load_model_season(exp_ctrl, var, season, lead, start)
    sens = load_model_season(exp_sens, var, season, lead, start)
    ctrl, sens = _common_grid(ctrl, sens, f"{var} {season} lead{lead} CTRL/SENS")
    years = np.intersect1d(ctrl["time"].values, sens["time"].values)
    if years.size < 3:
        raise ValueError(f"solo {years.size} anni in comune tra CTRL e SENS")
    ctrl, sens = ctrl.sel(time=years), sens.sel(time=years)

    p, delta = af.block_bootstrap_delta_pvalue(sens, ctrl, None, statistic="bias",
                                               n_iter=n_iter, block=block)
    d = (sens.mean("member") - ctrl.mean("member")).astype("float32")
    ds = xr.Dataset({
        "d": d.transpose("time", "lat", "lon"),
        "delta": delta.astype("float32"),
        "p": p,
        "clim_ctrl": ctrl.mean(("member", "time")).astype("float32"),
        "clim_sens": sens.mean(("member", "time")).astype("float32"),
    })
    ds.attrs.update(var=var, season=season, lead=int(lead), exp_ctrl=exp_ctrl, exp_sens=exp_sens,
                    n_members_ctrl=int(ctrl.sizes["member"]), n_members_sens=int(sens.sizes["member"]),
                    years=" ".join(str(int(y)) for y in years), n_iter=int(n_iter), block=int(block),
                    test="af.block_bootstrap_delta_pvalue(statistic='bias')")
    out = Path(out_dir) / f"delta_{var}_{season}_lead{lead}_1999.nc"
    ds.to_netcdf(out)
    return out


def run_one_model(args):
    """Adattatore per Pool.imap_unordered: ritorna una stringa di esito."""
    exp_ctrl, exp_sens, var, season, lead, out_dir, n_iter = args
    try:
        out = compute_model_delta(exp_ctrl, exp_sens, var, season, lead, out_dir, n_iter=n_iter)
        n = xr.open_dataset(out).sizes["time"]
        return f"{var} {season} lead{lead} ok ({n} anni) -> {out.name}"
    except Exception as e:
        return f"{var} {season} lead{lead} ERRORE: {type(e).__name__}: {e}"


# ---------------------------------------------------------------------------
# Vegetazione prescritta (cvh, cvl, lai_hv, lai_lv): serie mensili continue
# ---------------------------------------------------------------------------
def _cover_source(exp, var):
    """(file, nome della variabile nel file). Stessi file di cover_regression/
    (cover efficace) e di 02-vegetation_trend (LAI: per CTRL dal file di
    boundary condition icmcl, cioe' la vegetazione LPJ-GUESS)."""
    if var in ("cvh", "cvl"):
        if exp == "a1ua":
            return Path(POST_DATA) / f"{exp}_effective_{var}_1x1.nc", var
        return Path(POST_DATA) / f"{exp}_effective_{var}_199311-201910_1x1.nc", var
    if var in ("lai_hv", "lai_lv"):
        if exp == "a1ua":
            return Path(BC_DATA) / exp / f"icmcl_{exp}_regular_1x1.nc", var
        return (Path(CONFESS_DATA) / exp / "original_resolution" / var /
                f"{exp}_{var}_199311-201910_1x1.nc", var.upper())
    raise ValueError(f"variabile di vegetazione sconosciuta: {var!r}")


def _monthly_time(da, path):
    """Riporta ogni media mensile al primo giorno del mese a cui si riferisce.

    I file di SENS (nome con il periodo, es. ..._199311-201910_...) hanno 312
    mesi, novembre 1993 - ottobre 2019, ma datati al primo giorno del mese
    SUCCESSIVO (1993-12-01 ... 2019-11-01, la fine dell'intervallo di media;
    verificato sul cluster il 2026-11-23). I file di CTRL sono datati a meta'
    mese. Se il nome dichiara il periodo AAAAMM-AAAAMM, lo si usa come
    riferimento: date coerenti col periodo -> si tiene il mese; date spostate
    esattamente di un mese in avanti -> si riportano indietro di un mese;
    altrimenti errore. Senza periodo nel nome si tiene il mese della data.
    Ritorna (DataArray, True se le date sono state spostate)."""
    months = pd.DatetimeIndex(pd.to_datetime(da["time"].values)).to_period("M")
    if not months.is_unique or not months.is_monotonic_increasing:
        raise ValueError(f"{Path(path).name}: mesi ripetuti o non ordinati")
    shifted = False
    m = re.search(r"_(\d{4})(\d{2})-(\d{4})(\d{2})", Path(path).name)
    if m:
        first = pd.Period(f"{m.group(1)}-{m.group(2)}", freq="M")
        last = pd.Period(f"{m.group(3)}-{m.group(4)}", freq="M")
        if months[0] == first + 1 and months[-1] == last + 1:
            months, shifted = months - 1, True
        elif not (months[0] == first and months[-1] == last):
            raise ValueError(f"{Path(path).name}: le date vanno da {months[0]} a {months[-1]}, "
                             f"il nome dichiara {first}-{last} (ne' uguali ne' spostate di un mese)")
    return da.assign_coords(time=months.to_timestamp()), shifted


def seasonal_series(da, season):
    """Media stagionale per anno da una serie mensile. Etichetta = anno del
    primo mese della stagione (dicembre per DJF), come i file del modello
    (resample 'QS-DEC'). Solo stagioni complete (3 mesi)."""
    res = da.resample(time="QS-DEC").mean("time")
    n = xr.ones_like(da["time"], dtype=float).resample(time="QS-DEC").sum()
    keep = (res["time"].dt.month == SEASON_START_MONTH[season]).values & (n.values == 3)
    return res.isel(time=np.flatnonzero(keep))


def compute_cover_series(exp_ctrl, exp_sens, var, season, out_dir):
    """Salva <out_dir>/cover_<var>_<stagione>.nc con le serie stagionali per
    anno di CTRL e SENS (stessa griglia, anni in comune). Gli anni da usare
    si scelgono dopo, uguali a quelli del lead del modello con cui si confronta."""
    series, shifted = {}, []
    for role, exp in (("ctrl", exp_ctrl), ("sens", exp_sens)):
        path, name = _cover_source(exp, var)
        with xr.open_dataset(path) as ds:
            da = ds[name].squeeze(drop=True).load()
        da, moved = _monthly_time(da, path)   # date di SENS spostate di un mese: vedi _monthly_time
        if moved:
            shifted.append(exp)
        series[role] = _fix_grid(_to_years(seasonal_series(da, season)))
    ctrl, sens = _common_grid(series["ctrl"], series["sens"], f"{var} {season} CTRL/SENS")
    years = np.intersect1d(ctrl["time"].values, sens["time"].values)
    ds = xr.Dataset({"ctrl": ctrl.sel(time=years).astype("float32"),
                     "sens": sens.sel(time=years).astype("float32")})
    ds.attrs.update(var=var, season=season, exp_ctrl=exp_ctrl, exp_sens=exp_sens,
                    date_spostate_indietro_di_un_mese=" ".join(shifted) or "nessuno")
    out = Path(out_dir) / f"cover_{var}_{season}.nc"
    ds.to_netcdf(out)
    return out


def run_one_cover(args):
    """Adattatore per Pool.imap_unordered: ritorna una stringa di esito."""
    exp_ctrl, exp_sens, var, season, out_dir = args
    try:
        out = compute_cover_series(exp_ctrl, exp_sens, var, season, out_dir)
        with xr.open_dataset(out) as ds:
            t, moved = ds["time"].values, ds.attrs["date_spostate_indietro_di_un_mese"]
        note = "" if moved == "nessuno" else f", date di {moved} indietro di un mese"
        return f"{var} {season} ok ({t.size} anni, {t.min()}-{t.max()}{note}) -> {out.name}"
    except Exception as e:
        return f"{var} {season} ERRORE: {type(e).__name__}: {e}"


def check_inputs(exp_ctrl, exp_sens, model_vars=MODEL_VARS, cover_vars=COVER_VARS,
                 seasons=SEASONS, leads=LEADS_SINGLE):
    """Da eseguire prima del calcolo: quanti file per membro ci sono per ogni
    variabile/stagione/lead, e se esistono i file della vegetazione."""
    rows = []
    for var in model_vars:
        for season in seasons:
            row = {"variabile": var, "stagione": season}
            for exp in (exp_ctrl, exp_sens):
                row[exp] = " ".join(str(len(model_files(exp, var, season, k))) for k in leads)
            rows.append(row)
    df = pd.DataFrame(rows)
    print(f"File per membro in WORK_DIR/<exp>/ (stagionali, lead singolo), per lead {list(leads)}:")
    print(df.to_string(index=False))
    print("\nVegetazione prescritta:")
    for var in cover_vars:
        for exp in (exp_ctrl, exp_sens):
            path, _ = _cover_source(exp, var)
            print(f"  {exp} {var:7s} {'ok      ' if path.exists() else 'MANCANTE'} {path}")
    return df


# ---------------------------------------------------------------------------
# Lettura dei risultati e maschere (notebook delle figure)
# ---------------------------------------------------------------------------
def open_model_delta(out_dir, var, season, lead):
    path = Path(out_dir) / f"delta_{var}_{season}_lead{lead}_1999.nc"
    if not path.exists():
        raise FileNotFoundError(f"manca {path.name}: eseguire prima 01-mean_state_delta")
    with xr.open_dataset(path) as ds:
        return ds.load()


def cover_delta(out_dir, var, season, years):
    """Delta di vegetazione SENS - CTRL per anno, sugli anni richiesti (quelli
    disponibili). var='cvt' = cvh + cvl."""
    if var == "cvt":
        return cover_delta(out_dir, "cvh", season, years) + cover_delta(out_dir, "cvl", season, years)
    path = Path(out_dir) / f"cover_{var}_{season}.nc"
    if not path.exists():
        raise FileNotFoundError(f"manca {path.name}: eseguire prima 01-mean_state_delta")
    with xr.open_dataset(path) as ds:
        d = (ds["sens"] - ds["ctrl"]).load()
    return d.sel(time=np.intersect1d(years, d["time"].values))


def per_year_delta(out_dir, var, season, lead, years=None):
    """Delta SENS - CTRL per anno (time, lat, lon) per una variabile del
    modello o di vegetazione. Per la vegetazione, gli anni sono quelli del
    lead del modello (years, o se None quelli di tas per quel lead)."""
    if var in COVER_KINDS:
        if years is None:
            years = open_model_delta(out_dir, "tas", season, lead)["time"].values
        return cover_delta(out_dir, var, season, years)
    return open_model_delta(out_dir, var, season, lead)["d"]


def paired(dx, dy, what):
    """Stessi anni e stessa griglia per dx e dy."""
    years = np.intersect1d(dx["time"].values, dy["time"].values)
    if years.size < 3:
        raise ValueError(f"{what}: solo {years.size} anni in comune")
    return _common_grid(dx.sel(time=years), dy.sel(time=years), what)


def land_mask(template):
    """Terre (lsm > 0.5, come af.land_seas_mask) a nord di 60S."""
    m = af.land_seas_mask(xr.ones_like(template, dtype=float), "land").notnull()
    m = m.reindex_like(template, fill_value=False)
    return m & (template["lat"] > -60)


def box_mask(template, box=BOX):
    lat0, lat1, lon0, lon1 = box
    return ((template["lat"] >= lat0) & (template["lat"] <= lat1) &
            (template["lon"] >= lon0) & (template["lon"] <= lon1)).broadcast_like(template)


def snow_mask(out_dir, season, lead, template, thr):
    """Pixel con neve: snd medio (media di CTRL e SENS) > thr. None se il
    file di snd per quella stagione/lead non c'e'."""
    try:
        ds = open_model_delta(out_dir, "snd", season, lead)
    except FileNotFoundError:
        return None
    clim = 0.5 * (ds["clim_ctrl"] + ds["clim_sens"])
    return (clim > thr).reindex_like(template, fill_value=False)


def strata(kind, out_dir, season, lead, template, snow_thr):
    """Sottoinsiemi di pixel su cui calcolare la relazione separatamente:
    'neve' (con / senza neve) oppure 'lat' (fasce di latitudine)."""
    if kind == "neve":
        snow = snow_mask(out_dir, season, lead, template, snow_thr)
        if snow is None:
            return {"tutte (snd mancante)": xr.ones_like(template, dtype=bool)}
        return {f"con neve (snd>{snow_thr:g} m)": snow, "senza neve": ~snow}
    if kind == "lat":
        lat = template["lat"].broadcast_like(template)
        alat = np.abs(lat)
        return {"tropici (<23.5)": alat < 23.5,
                "medie lat. (23.5-50)": (alat >= 23.5) & (alat < 50),
                "boreali (>50N)": lat >= 50}
    raise ValueError(f"strati sconosciuti: {kind!r} (usare 'neve' o 'lat')")


# ---------------------------------------------------------------------------
# Relazione tra pixel per classi di delta x
# ---------------------------------------------------------------------------
def _block_indices(rng, n_t, block):
    """Anni ricampionati a blocchi di `block` anni consecutivi, CIRCOLARI
    (Politis e Romano 1992): a differenza dei blocchi mobili di
    af._block_bootstrap_delta, il primo e l'ultimo anno entrano nei campioni
    con la stessa probabilita' degli altri, e la media bootstrap resta
    centrata sulla media del campione anche se c'e' un trend (con 15-20 anni
    e blocchi di 5 i blocchi mobili sottopesano molto gli estremi)."""
    block = min(block, n_t)
    n_blocks = int(np.ceil(n_t / block))
    starts = rng.integers(0, n_t, n_blocks)
    return ((starts[:, None] + np.arange(block)) % n_t).ravel()[:n_t]


def _weighted_fit(x, y, w):
    """Pendenza e correlazione di Pearson pesate (pesi cos(lat))."""
    sw = w.sum()
    mx, my = (w * x).sum() / sw, (w * y).sum() / sw
    cxy = (w * (x - mx) * (y - my)).sum() / sw
    vx = (w * (x - mx) ** 2).sum() / sw
    vy = (w * (y - my) ** 2).sum() / sw
    if vx <= 0 or vy <= 0:
        return np.nan, np.nan
    return cxy / vx, cxy / np.sqrt(vx * vy)


def binned_relation(dx, dy, mask, n_bins=10, n_boot=1000, block=5, seed=0, min_count=5):
    """Relazione tra pixel delle differenze medie: x = media nel tempo di dx,
    y = media nel tempo di dy, sui pixel di mask.

    - Classi di x a quantili (stesso numero di pixel per classe); in ogni
      classe la media di x e di y pesata con cos(lat).
    - Pendenza e r di Pearson pesati, rho di Spearman.
    - Incertezza: bootstrap a blocchi circolari di anni (blocchi di 5 come nel
      resto del repo) su y, con x FISSO alla media su tutti gli anni, come in
      una regressione a x fisso. Ricampionare anche x, quando x e' una
      variabile del modello col suo rumore (es. l'albedo), schiaccerebbe le
      pendenze bootstrap verso 0 (meno anni distinti -> x piu' rumoroso) e
      l'intervallo non conterrebbe piu' la stima. Lo stesso rumore attenua
      anche la pendenza stimata: dove r e' basso la pendenza va presa con
      cautela. L'intervallo dice quanto la relazione dipende dagli anni e
      dalla variabilita' interna, non tiene conto della correlazione spaziale.
    dx, dy: (time, lat, lon) con stessi anni e griglia; mask: (lat, lon) bool."""
    m = np.asarray(mask.transpose("lat", "lon").values, dtype=bool).ravel()
    n_t = dy.sizes["time"]
    X = np.asarray(dx.transpose("time", "lat", "lon").values, dtype=float).reshape(n_t, -1)[:, m]
    Y = np.asarray(dy.transpose("time", "lat", "lon").values, dtype=float).reshape(n_t, -1)[:, m]
    lat2d = np.broadcast_to(dy["lat"].values[:, None], (dy.sizes["lat"], dy.sizes["lon"]))
    w = np.cos(np.deg2rad(lat2d)).ravel()[m]

    with np.errstate(invalid="ignore"), _quiet():
        xm, ym = np.nanmean(X, 0), np.nanmean(Y, 0)
    ok = np.isfinite(xm) & np.isfinite(ym)
    X, Y, w, xm, ym = X[:, ok], Y[:, ok], w[ok], xm[ok], ym[ok]
    out = {"n": int(ok.sum()), "x": xm, "y": ym}
    if out["n"] < max(3 * min_count, 10) or np.ptp(xm) == 0:
        return None

    edges = np.unique(np.quantile(xm, np.linspace(0, 1, n_bins + 1)))
    if edges.size < 3:
        return None
    nb = edges.size - 1
    idx = np.clip(np.searchsorted(edges, xm, side="right") - 1, 0, nb - 1)
    counts = np.bincount(idx, minlength=nb)
    sw = np.bincount(idx, weights=w, minlength=nb)

    def bin_mean(v):
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.bincount(idx, weights=w * v, minlength=nb) / sw

    slope, r = _weighted_fit(xm, ym, w)
    rng = np.random.default_rng(seed)
    boot_y = np.empty((n_boot, nb))
    boot_slope = np.empty(n_boot)
    for i in range(n_boot):
        t = _block_indices(rng, n_t, block)
        with np.errstate(invalid="ignore"), _quiet():
            yb = np.nanmean(Y[t], 0)
        yb = np.where(np.isfinite(yb), yb, ym)   # pixel senza dati negli anni estratti
        boot_y[i] = bin_mean(yb)
        boot_slope[i] = _weighted_fit(xm, yb, w)[0]

    keep = counts >= min_count
    out.update(
        bin_x=bin_mean(xm)[keep], bin_y=bin_mean(ym)[keep],
        bin_lo=np.nanpercentile(boot_y, 2.5, axis=0)[keep],
        bin_hi=np.nanpercentile(boot_y, 97.5, axis=0)[keep],
        bin_n=counts[keep],
        slope=slope, slope_lo=np.nanpercentile(boot_slope, 2.5),
        slope_hi=np.nanpercentile(boot_slope, 97.5),
        r=r, rho=stats.spearmanr(xm, ym)[0],
    )
    return out


class _quiet:
    """Silenzia i RuntimeWarning di nanmean su colonne tutte NaN."""
    def __enter__(self):
        import warnings
        self._w = warnings.catch_warnings()
        self._w.__enter__()
        warnings.simplefilter("ignore", RuntimeWarning)

    def __exit__(self, *exc):
        return self._w.__exit__(*exc)


def _label(v):
    return f"Δ{v} SENS−CTRL [{UNITS.get(v, '')}]"


def plot_relation(out_dir, xvar, yvar, kind, lead, png, seasons=SEASONS, box=BOX,
                  snow_thr=0.01, n_bins=10, n_bins_box=5, n_boot=1000, max_points=3000):
    """Una figura per (x, y): righe = dominio (terre globali, box Siberia),
    colonne = stagioni. In ogni pannello la nuvola dei pixel (grigia), le
    medie per classe di x con l'intervallo bootstrap al 95% per ciascuno
    strato, e in legenda pendenza [IC 95%], r e numero di pixel.
    Ritorna la tabella delle statistiche (un rigo per pannello e strato)."""
    domains = ("Terre (>60S)", f"Siberia {box[0]}-{box[1]}N {box[2]}-{box[3]}E")
    fig, axes = plt.subplots(2, len(seasons), figsize=(5.2 * len(seasons), 9.5), squeeze=False)
    rows = []
    rng = np.random.default_rng(0)
    for j, season in enumerate(seasons):
        try:
            dy = per_year_delta(out_dir, yvar, season, lead)
            dx = per_year_delta(out_dir, xvar, season, lead, years=dy["time"].values)
            dx, dy = paired(dx, dy, f"{xvar}/{yvar} {season} lead{lead}")
        except (FileNotFoundError, ValueError) as e:
            for i in range(2):
                axes[i, j].text(0.5, 0.5, f"{season}: {e}", ha="center", va="center",
                                wrap=True, fontsize=8, transform=axes[i, j].transAxes)
            continue
        template = dy.isel(time=0, drop=True)
        land = land_mask(template)
        dom_masks = (land, land & box_mask(template, box))
        st = strata(kind, out_dir, season, lead, template, snow_thr)
        for i, (dname, dmask) in enumerate(zip(domains, dom_masks)):
            ax = axes[i, j]
            for k, (sname, smask) in enumerate(st.items()):
                color = STRATA_COLORS[k % len(STRATA_COLORS)]
                res = binned_relation(dx, dy, dmask & smask, n_bins=n_bins if i == 0 else n_bins_box,
                                      n_boot=n_boot)
                if res is None:
                    continue
                sel = rng.permutation(res["n"])[:max_points]
                ax.scatter(res["x"][sel], res["y"][sel], s=3, color=color, alpha=0.12, lw=0)
                ax.errorbar(res["bin_x"], res["bin_y"],
                            yerr=[res["bin_y"] - res["bin_lo"], res["bin_hi"] - res["bin_y"]],
                            color=color, marker="o", ms=4, capsize=2, lw=1.2,
                            label=(f"{sname}: b={res['slope']:.3g} [{res['slope_lo']:.3g}, "
                                   f"{res['slope_hi']:.3g}], r={res['r']:.2f}, n={res['n']}"))
                rows.append({"x": xvar, "y": yvar, "lead": lead, "stagione": season,
                             "dominio": dname, "strato": sname, "n_pixel": res["n"],
                             "pendenza": res["slope"], "IC_2.5": res["slope_lo"],
                             "IC_97.5": res["slope_hi"], "r_pearson": res["r"],
                             "rho_spearman": res["rho"]})
            ax.axhline(0, color="k", lw=0.6)
            ax.axvline(0, color="k", lw=0.6)
            ax.set_title(f"{dname} - {season}", fontsize=11)
            ax.set_xlabel(_label(xvar), fontsize=9)
            ax.set_ylabel(_label(yvar), fontsize=9)
            if ax.get_legend_handles_labels()[0]:
                ax.legend(fontsize=7, loc="best")
    fig.suptitle(f"Δ{yvar} contro Δ{xvar}, stato medio, lead {lead}", fontsize=14)
    fig.tight_layout()
    fig.savefig(png, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Coerenza tra lead e mappe
# ---------------------------------------------------------------------------
def _area_mean(d, mask):
    w = np.cos(np.deg2rad(d["lat"]))
    return d.where(mask).weighted(w).mean(("lat", "lon"))


def _mean_ci(series, n_boot, block, seed=0):
    v = np.asarray(series, dtype=float)
    v = v[np.isfinite(v)]
    if v.size < 3:
        return np.nan, np.nan, np.nan
    rng = np.random.default_rng(seed)
    boot = np.array([v[_block_indices(rng, v.size, block)].mean() for _ in range(n_boot)])
    return v.mean(), np.percentile(boot, 2.5), np.percentile(boot, 97.5)


def plot_lead_consistency(out_dir, png, variables=("alb", "snd", "tas"), seasons=SEASONS,
                          leads=LEADS_SINGLE, box=BOX, n_boot=1000, block=5):
    """Media d'area della differenza SENS - CTRL in funzione del lead, per
    stagione, nel box Siberia e nelle terre boreali (>50N), con intervallo
    bootstrap al 95% sugli anni. Se il meccanismo e' la vegetazione
    prescritta, la differenza deve essere simile a tutti i lead."""
    domains = {f"Siberia {box[0]}-{box[1]}N {box[2]}-{box[3]}E": lambda t: land_mask(t) & box_mask(t, box),
               "Terre boreali (>50N)": lambda t: land_mask(t) & (t["lat"] >= 50)}
    fig, axes = plt.subplots(len(variables), len(domains), figsize=(6.5 * len(domains), 3.4 * len(variables)),
                             squeeze=False)
    rows = []
    offs = np.linspace(-0.2, 0.2, len(seasons))
    for i, var in enumerate(variables):
        for j, (dname, make_mask) in enumerate(domains.items()):
            ax = axes[i, j]
            for off, season in zip(offs, seasons):
                stat = []
                for lead in leads:
                    try:
                        d = open_model_delta(out_dir, var, season, lead)["d"]
                    except FileNotFoundError:
                        stat.append((np.nan, np.nan, np.nan))
                        continue
                    s = _area_mean(d, make_mask(d.isel(time=0, drop=True)))
                    stat.append(_mean_ci(s.values, n_boot, block))
                    rows.append({"variabile": var, "dominio": dname, "stagione": season, "lead": lead,
                                 "delta": stat[-1][0], "IC_2.5": stat[-1][1], "IC_97.5": stat[-1][2]})
                mean, lo, hi = (np.array(v) for v in zip(*stat))
                ax.errorbar(np.asarray(leads) + off, mean, yerr=[mean - lo, hi - mean], marker="o",
                            ms=4, capsize=2, color=SEASON_COLORS.get(season), label=season)
            ax.axhline(0, color="k", lw=0.6)
            ax.set_xticks(list(leads))
            ax.set_xlabel("lead (anni)")
            ax.set_ylabel(_label(var), fontsize=9)
            ax.set_title(dname, fontsize=11)
            if i == 0 and j == 0:
                ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(png, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return pd.DataFrame(rows)


def auto_levels(data, n=10):
    """Livelli simmetrici attorno a 0 dal 98esimo percentile di |data|
    (ignorando gli zeri esatti, es. pixel senza vegetazione), arrotondati a
    una cifra significativa."""
    v = np.abs(np.asarray(data.values, dtype=float))
    v = v[np.isfinite(v) & (v > 0)]
    if v.size == 0:
        return [-1.0, -0.5, 0.0, 0.5, 1.0]
    vmax = float(np.percentile(v, 98))
    mag = 10 ** np.floor(np.log10(vmax))
    vmax = float(np.ceil(vmax / mag) * mag)
    edges = np.linspace(0, vmax, n // 2 + 1)[1:]
    return [float(-e) for e in edges[::-1]] + [0.0] + [float(e) for e in edges]


def plot_maps(out_dir, season, lead, variables, fig_dir, sign=0.95):
    """Una mappa per variabile (af.map_plot, come le altre figure): differenza
    di stato medio SENS - CTRL. Per le variabili del modello la 'x' marca i
    punti NON significativi (p >= 1 - sign, bootstrap a blocchi); la
    vegetazione e' prescritta, quindi senza test. Ritorna i PNG."""
    pngs = []
    for var in variables:
        try:
            if var in COVER_KINDS:
                delta = per_year_delta(out_dir, var, season, lead).mean("time")
                p = xr.zeros_like(delta)
                cmap = "BrBG"
            else:
                ds = open_model_delta(out_dir, var, season, lead)
                delta, p, cmap = ds["delta"], ds["p"], "bwr"
        except (FileNotFoundError, ValueError) as e:
            print(f"  {var} {season} lead{lead}: saltata ({e})")
            continue
        af.map_plot(delta, p, levels=auto_levels(delta), antartica=False, type="pvalue", sign=sign,
                    cmap=cmap, title=f"Δ{var} SENS−CTRL {season}, lead {lead} [{UNITS.get(var, '')}]")
        png = Path(fig_dir) / f"mech_delta_{var}_{season}_lead{lead}_1999.png"
        plt.savefig(png, dpi=200, bbox_inches="tight")
        plt.close("all")
        pngs.append(png)
    return pngs


def composite(pngs, out, ncols=4, title=None):
    """Affianca i PNG in una griglia (come le figure composite del repo)."""
    if not pngs:
        return None
    nrows = int(np.ceil(len(pngs) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(7 * ncols, 4.2 * nrows), squeeze=False)
    for ax in axes.flat:
        ax.axis("off")
    for ax, p in zip(axes.flat, pngs):
        ax.imshow(mpimg.imread(p))
    if title:
        fig.suptitle(title, fontsize=16)
    fig.tight_layout()
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return out
