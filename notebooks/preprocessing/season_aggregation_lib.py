"""Aggregazione stagionale (DJF, MAM, JJA, SON) dei file mensili delle
previsioni decadali, per lead singolo e per combinazione di lead year.

Portato da decadal_tool/decadal/00-process_2025/03-decadal_season_aggregation
(originale, fuori dal repo), con tre differenze:

1. Date: ogni file mensile viene aperto con af.fix_monthly_time. I file non
   CMOR di SENS (es. snd_s1999-r1i1p1f1_199911-200010_1x1.nc, medie mensili
   degli output IFS) sono datati al primo giorno del mese DOPO; l'originale
   raggruppava i mesi per data, quindi per quei file le stagioni erano
   sfasate di un mese (il "JJA" della neve di SENS era mag-lug e il "DJF"
   solo dic-gen). I file CMOR, datati a meta' mese, restano invariati.
2. Le medie su piu' lead year usano solo le inizializzazioni con TUTTI gli
   anni richiesti (l'originale prendeva i file disponibili per posizione: se
   ne mancava uno, la finestra scivolava sugli anni dopo).
3. Esecuzione in processi spawn freschi (multiprocessing.Pool,
   maxtasksperchild=1) come nel resto del repo, invece di ProcessPoolExecutor.
4. Ogni file mensile deve contenere la variabile richiesta (var), altrimenti
   errore; nei file in uscita c'e' solo quella variabile. L'originale
   stagionalizzava qualunque contenuto: nella cartella alb del cluster ci
   sono file con 'bb_albedo' e valori non fisici per l'albedo del modello.

Per il resto stessi passaggi e stessi nomi di file dell'originale, cosi' i
notebook che li leggono (seasons/, cover_mechanism/) non cambiano. In
RUN_PATH = WORK_DIR/<exp>:
  <exp>_<var>_<realm>_EC-Earth3_dcppA-hindcast_s<i>-<m>_gr_<a>11-<a+1>10_1x1_<stag>_rad.nc
      media stagionale per inizializzazione i, membro m, anno di previsione a
  <exp>_<var>_<realm>_EC-Earth3_dcppA-hindcast_<m>_gr_1x1_lead<k>_rad_<stag>.nc
      serie sulle inizializzazioni per lead singolo k
  <exp>_<var>_<realm>_EC-Earth3_dcppA-hindcast_s<i>-<m>_lead_<y1>-<y2>_1x1_m<stag>_rad.nc
      media sui lead year y1..y2 per inizializzazione
In OUT_PATH = POST_DATA/<exp>/1x1/<var>:
  <exp>_<var>_<realm>_EC-Earth3_dcppA-hindcast_<m>_lead_<y1>-<y2>_1x1_m<stag>_rad.nc
  <exp>_<var>_<realm>_EC-Earth3_dcppA-hindcast_lead_<y1>-<y2>_1x1_ensemble_m<stag>_rad.nc

Stagioni come l'originale: ogni file mensile copre nov(a) - ott(a+1);
novembre si scarta, poi DJF = dic-feb, MAM, JJA e SON = set-ott (2 soli mesi:
novembre sta nel file dell'anno dopo). Etichetta = primo giorno della
stagione (resample trimestrale 'QS-DEC'), uguale a quella dell'originale.
"""
import os
import re
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

from pathlib import Path

import pandas as pd
import xarray as xr

from config import CONFESS_DATA, POST_DATA, WORK_DIR
import albedo_functions as af

SEASONS = ["DJF", "MAM", "JJA", "SON"]
SEASON_OF_START_MONTH = {12: "DJF", 3: "MAM", 6: "JJA", 9: "SON"}
MONTHS_PER_SEASON = {"DJF": 3, "MAM": 3, "JJA": 3, "SON": 2}   # SON senza novembre, come l'originale
N_LEADS = 5
TIME_ENCODING = {"time": {"units": "days since 1994-06-01T12:00:00", "dtype": "float64"}}  # come l'originale


def paths(exp, var):
    """(cartella dei file mensili, RUN_PATH, OUT_PATH), come l'originale."""
    return (Path(CONFESS_DATA) / exp / "1x1" / "monthly_start_dates" / var,
            Path(WORK_DIR) / exp,
            Path(POST_DATA) / exp / "1x1" / var)


def stem(exp, var, realm):
    return f"{exp}_{var}_{realm}_EC-Earth3_dcppA-hindcast"


def monthly_file(exp, var, start, member, year):
    """File mensile dell'anno di previsione nov(year) - ott(year+1). Trova sia
    il nome CMOR ({var}_<realm>_EC-Earth3_dcppA-hindcast_s<i>-<m>_gr_<a>11-<a+1>10_1x1.nc)
    sia quello non CMOR ({var}_s<i>-<m>_<a>11-<a+1>10_1x1.nc, es. la neve di
    SENS). None se manca, errore se i candidati sono piu' di uno."""
    data, _, _ = paths(exp, var)
    found = sorted(data.glob(f"{var}_*s{start}-{member}_*{year}11-{year + 1}10_1x1.nc"))
    if len(found) > 1:
        raise ValueError(f"piu' file per s{start} {member} {year}: {[f.name for f in found]}")
    return found[0] if found else None


def _check_variable(ds, var, path):
    """Il file deve contenere la variabile 'var'. Errore altrimenti: nella
    cartella monthly_start_dates/alb del cluster ci sono file con la
    variabile 'bb_albedo' (nome di GLASS) e un andamento non fisico per
    l'albedo del modello (2026-11-23); senza questo controllo verrebbero
    stagionalizzati e sovrascriverebbero i file buoni."""
    if var not in ds.data_vars:
        raise ValueError(f"{Path(path).name}: manca la variabile '{var}' "
                         f"(il file contiene {list(ds.data_vars)}): non e' il file atteso, fermarsi e controllare")


def seasonal_means(path, var):
    """Medie stagionali della variabile 'var' di un file mensile nov-ott:
    (dict stagione -> Dataset con un solo tempo, True se le date sono state
    spostate)."""
    with xr.open_dataset(path) as ds:
        _check_variable(ds, var, path)
        ds = ds[[var]].load()      # solo la variabile (senza limiti temporali e di griglia)
    ds, shifted = af.fix_monthly_time(ds, path, return_shift=True, verbose=False)
    months = ds["time"].dt.month.values
    if ds.sizes["time"] != 12 or months[0] != 11:
        raise ValueError(f"{Path(path).name}: attesi 12 mesi da novembre, trovati {ds.sizes['time']} "
                         f"da {str(ds['time'].values[0])[:7]}")
    ds = ds.isel(time=slice(1, None))            # via novembre, come l'originale
    means = ds.resample(time="QS-DEC").mean("time", keep_attrs=True)
    counts = ds["time"].resample(time="QS-DEC").count()
    out = {}
    for t, n in zip(means["time"].values, counts.values):
        season = SEASON_OF_START_MONTH[pd.Timestamp(t).month]
        if n != MONTHS_PER_SEASON[season]:
            raise ValueError(f"{Path(path).name}: {season} con {n} mesi invece di {MONTHS_PER_SEASON[season]}")
        out[season] = means.sel(time=[t])
    return out, shifted


def existing_realms(exp, var):
    """Realm nei nomi dei file stagionali per lead singolo gia' presenti in
    RUN_PATH (per sovrascriverli invece di creare doppioni)."""
    _, run, _ = paths(exp, var)
    pat = re.compile(rf"^{exp}_{var}_(\w+?)_EC-Earth3_dcppA-hindcast_r\d+i\d+p\d+f\d+_gr_1x1_lead\d_rad_"
                     rf"(DJF|MAM|JJA|SON)\.nc$")
    return sorted({m.group(1) for f in run.glob(f"{exp}_{var}_*_lead*_rad_*.nc") if (m := pat.match(f.name))})


def check_inputs(exp, var, start_years, members):
    """Da eseguire prima: quanti file mensili ci sono per membro, se le date
    di un file d'esempio vanno corrette, e il realm dei file gia' prodotti."""
    expected = len(start_years) * N_LEADS
    counts, example = [], None
    for m in members:
        n = 0
        for i in start_years:
            for k in range(N_LEADS):
                f = monthly_file(exp, var, i, m, i + k)
                if f is not None:
                    n += 1
                    example = example or f
        counts.append(f"{m}:{n}")
    print(f"{exp} {var}: file mensili per membro (attesi {expected}):")
    print("   " + "  ".join(counts))
    if example is not None:
        with xr.open_dataset(example) as ds:
            _check_variable(ds, var, example)   # si ferma qui, prima di calcolare
            t = ds["time"].values
        _, shifted = af.fix_monthly_time(xr.Dataset(coords={"time": t}), example, return_shift=True, verbose=False)
        print(f"   esempio {example.name}: date {str(t[0])[:10]} ... {str(t[-1])[:10]} -> "
              f"{'datate al mese dopo, verranno riportate indietro di un mese' if shifted else 'gia giuste'}")
    realms = existing_realms(exp, var)
    print(f"   file stagionali gia' presenti in {paths(exp, var)[1]}: realm {realms or 'nessuno'}")
    return realms


def run_member(args):
    """Tutti i passaggi che dipendono da un solo membro: medie stagionali per
    inizializzazione e anno, serie per lead singolo, medie sui lead year e
    serie per membro. Ritorna una stringa di esito (processi spawn)."""
    exp, var, realm, member, start_years = args
    _, run, out = paths(exp, var)
    run.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)
    s = stem(exp, var, realm)
    try:
        per_year, missing, n_files, n_shift = {}, [], 0, 0
        for i in start_years:
            for k in range(N_LEADS):
                a = i + k
                f = monthly_file(exp, var, i, member, a)
                if f is None:
                    missing.append(f"s{i}/{a}")
                    continue
                means, shifted = seasonal_means(f, var)
                n_files += 1
                n_shift += shifted
                for season, ds in means.items():
                    ds.to_netcdf(run / f"{s}_s{i}-{member}_gr_{a}11-{a + 1}10_1x1_{season}_rad.nc")
                    per_year[(i, k, season)] = ds

        for season in SEASONS:
            for k in range(N_LEADS):
                parts = [per_year[(i, k, season)] for i in start_years if (i, k, season) in per_year]
                if parts:
                    xr.concat(parts, "time").to_netcdf(run / f"{s}_{member}_gr_1x1_lead{k}_rad_{season}.nc")

            for y1 in range(N_LEADS):
                for y2 in range(y1 + 1, N_LEADS):
                    series = []
                    for i in start_years:
                        parts = [per_year.get((i, k, season)) for k in range(y1, y2 + 1)]
                        if any(p is None for p in parts):
                            continue          # solo inizializzazioni con tutti gli anni y1..y2
                        cat = xr.concat(parts, "time")
                        mean = (cat.mean("time", keep_attrs=True).expand_dims("time")
                                .assign_coords(time=[cat["time"].mean().values]))
                        mean.to_netcdf(run / f"{s}_s{i}-{member}_lead_{y1}-{y2}_1x1_m{season}_rad.nc")
                        series.append(mean)
                    if series:
                        xr.concat(series, "time").to_netcdf(
                            out / f"{s}_{member}_lead_{y1}-{y2}_1x1_m{season}_rad.nc", encoding=TIME_ENCODING)

        msg = f"{member} ok: {n_files} file mensili, "
        msg += f"date di {n_shift} riportate indietro di un mese" if n_shift else "date gia' giuste"
        if missing:
            msg += f"; mancano {len(missing)} anni ({', '.join(missing[:4])}{' ...' if len(missing) > 4 else ''})"
        return msg
    except Exception as e:
        return f"{member} ERRORE: {type(e).__name__}: {e}"


def run_ensemble(args):
    """Per una stagione: unisce i membri di ogni combinazione di lead year
    (file ..._ensemble_m<stag>_rad.nc). Ritorna una stringa di esito."""
    exp, var, realm, season = args
    _, _, out = paths(exp, var)
    s = stem(exp, var, realm)
    try:
        done = []
        for y1 in range(N_LEADS):
            for y2 in range(y1 + 1, N_LEADS):
                files = sorted(out.glob(f"{s}_r*_lead_{y1}-{y2}_1x1_m{season}_rad.nc"))
                if not files:
                    continue
                parts = []
                for f in files:
                    with xr.open_dataset(f) as ds:
                        parts.append(ds.load())
                xr.concat(parts, "member").to_netcdf(out / f"{s}_lead_{y1}-{y2}_1x1_ensemble_m{season}_rad.nc")
                done.append(f"{y1}-{y2}({len(files)})")
        return f"{season} ok: {' '.join(done) if done else 'nessun file per membro'}"
    except Exception as e:
        return f"{season} ERRORE: {type(e).__name__}: {e}"
