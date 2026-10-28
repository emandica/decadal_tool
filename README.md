# decadal_tool

Analisi di **esperimenti di previsione decadale** (progetto CONFESS, EC-Earth3)
e produzione delle figure per la relativa pubblicazione.

Il confronto è tra due set di esperimenti DCPP:

| Sigla  | Ruolo                | Etichetta nelle figure |
|--------|----------------------|------------------------|
| `a1ua` | esperimento *control*     | DCPP-CTRL |
| `a52o` | esperimento *sensitivity* | DCPP-SENS |

## Struttura del repository

```
decadal_tool/
├── environment.yml        # ambiente conda per eseguire i notebook
├── notebooks/             # pipeline di calcolo e plotting
│   ├── config.py             # percorsi centralizzati (radice sovrascrivibile via env)
│   ├── albedo_functions.py   # funzioni condivise (griglie, medie, maschere, bootstrap, plotting mappe)
│   ├── bias_lib.py           # mappe di bias per combinazione di lead year (04-BIAS, 04-Fig4)
│   ├── bias_leadyears_lib.py # bias globale per combinazione di lead year (03-Fig3)
│   ├── NN-<nome>.ipynb       # notebook di CALCOLO   -> salva .nc intermedi
│   ├── NN-FigN_<nome>.ipynb  # notebook di PLOTTING  -> genera le figure
│   ├── seasons/              # varianti stagionali (DJF/MAM/JJA/SON), con i loro *_lib.py
│   ├── cover_regression/     # regressioni delta cover / albedo / neve vs delta skill (+ seasons/)
│   └── altre_figure/         # figure supplementari
└── figures/               # figure finali della pubblicazione (PNG)
```

I moduli `*_lib.py` contengono le funzioni eseguite in parallelo: i notebook le
lanciano in processi `spawn` freschi (`multiprocessing.Pool(maxtasksperchild=1)`),
che non vedono le funzioni definite nel notebook.

## Pipeline

Ogni figura nasce da una **coppia** di notebook:

1. **Calcolo** — `NN-<nome>.ipynb` legge i dati grezzi degli esperimenti,
   calcola la statistica (deviazione standard interannuale, trend, bias, ACC…)
   e salva risultati intermedi in formato NetCDF.
2. **Plotting** — `NN-FigN_<nome>.ipynb` rilegge i NetCDF intermedi, produce i
   singoli pannelli (control / sensitivity / differenza) e li assembla nella
   figura composita finale salvata in `figures/`.

## Mappa notebook → figura

| Figura | Notebook di plotting | Contenuto |
|--------|----------------------|-----------|
| Fig. 1 | `01-Fig1_lai_std` · `01-Fig1_effective_lai_std` · `01-Fig1_cv_std` | Deviazione standard interannuale di LAI / LAI efficace / vegetation cover |
| Fig. 2 | `02-Fig2_lai_trend` · `02-Fig2_effective_lai_trend` · `02-Fig2_cv_trend` | Trend delle stesse variabili |
| Fig. 3 | `03-Fig3_global_bias_tas` | Bias globale della temperatura (tas) |
| Fig. 4 | `04-Fig4_BIAS_plot` (calcolo: `04-BIAS`) | Mappe di bias |
| Fig. 5 | `05-Serie_temporale_anomalie` · `05-Fig5_global_acc` (calcolo: `05-global_acc`) | Serie temporali delle anomalie / ACC di media d'area per lead time e per combinazioni di lead year (con significatività bootstrap) |
| Mappe ACC | `ACC_plot` (calcolo: `yearly_ACC_calculation` · `seasonal_ACC_calculation`) | Mappe spaziali di ACC (anomaly correlation) con significatività bootstrap; `seasonal_ACC_calculation` disegna anche le mappe con correzione FDR |
| Regressioni | `cover_regression/01`–`07` (calcolo e figure nello stesso notebook) | Delta skill (tas, albedo, neve) vs delta cover (cvh, cvl, cvt = cvh + cvl), albedo vs neve: mappe di slope / r per pixel e scatter sul box Siberia; figure in una sottocartella per notebook |
| Suppl. | `altre_figure/`, `seasons/` | Cover, LAI e varianti stagionali |

> I notebook seguono in gran parte la convenzione `NN-<nome>` (calcolo) + `NN-FigN_<nome>`
> (plotting); alcuni ausiliari (es. `ACC_plot`, `yearly_ACC_calculation`) hanno nomi propri.
> Ogni notebook scrive il proprio log in `<nome-notebook>.log` (non versionato).

## Significatività della differenza SENS − CTRL

Le differenze di ACC e di bias tra i due esperimenti sono testate col
**bootstrap a blocchi** di Goddard et al. (2013), nella ricetta di Smith et al.
(2019, Methods): anni ricampionati con reinserimento a blocchi di 5, gli stessi
per SENS, CTRL e osservazione; membri ricampionati con reinserimento per
ciascun anno, separatamente per SENS e CTRL; 1000 iterazioni; test bilaterale
della differenza nulla. Su dati sintetici senza vera differenza è
leggermente conservativo (circa 1–2,5% di falsi positivi a p < 0.05).

Ogni notebook di calcolo interessato (`yearly_ACC_calculation`,
`seasonal_ACC_calculation`, `05-global_acc`, `04-BIAS`, `seasons/04-BIAS_seasons`)
ha nella cella dei parametri:

```python
BOOTSTRAP = 'block'   # oppure 'pooled'
```

`'pooled'` usa il vecchio ricampionamento (membri di SENS e CTRL mescolati in un
unico pool, anni fissi), conservato come alternativa ma mal calibrato (12–20% di
falsi positivi). Le funzioni stanno in `albedo_functions.py` (`block_bootstrap_*`
e le vecchie `bootstrap_*`) e salvano i file nello stesso formato, quindi i
notebook di plotting non dipendono dal metodo scelto.

## Ambiente

```bash
conda env create -f environment.yml
conda activate decadal
jupyter lab
```

## Dati e percorsi

Tutti i percorsi dei dati (compresa la maschera terra-mare usata da
`land_seas_mask`) sono centralizzati in [`notebooks/config.py`](notebooks/config.py)
e derivano da un'unica radice. Sul cluster di calcolo la radice è `/confess/dicarlo`
(default); per eseguire altrove basta impostare una variabile d'ambiente prima di
avviare Jupyter:

```bash
export DECADAL_DATA_ROOT=/percorso/locale/ai/dati
jupyter lab
```

I dataset grezzi (esperimenti, rianalisi ERA5, boundary conditions) **non** sono
inclusi nel repository per dimensione: risiedono sul cluster sotto la radice sopra.
