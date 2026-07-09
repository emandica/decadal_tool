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
│   ├── albedo_functions.py   # funzioni condivise (griglie, medie, plotting mappe)
│   ├── NN-<nome>.ipynb       # notebook di CALCOLO   -> salva .nc intermedi
│   ├── NN-FigN_<nome>.ipynb  # notebook di PLOTTING  -> genera le figure
│   ├── seasons/              # varianti stagionali (DJF/MAM/JJA/SON)
│   └── altre_figure/         # figure supplementari
└── figures/               # figure finali della pubblicazione (PNG)
```

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
| Fig. 5 | `05-Serie_temporale_anomalie` · `05-Fig5_global_acc` | Serie temporali delle anomalie / ACC globale |
| Suppl. | `altre_figure/`, `seasons/` | Cover, LAI e varianti stagionali |

## Ambiente

```bash
conda env create -f environment.yml
conda activate decadal
jupyter lab
```

## Dati e percorsi

Tutti i percorsi dei dati sono centralizzati in [`notebooks/config.py`](notebooks/config.py)
e derivano da un'unica radice. Sul cluster di calcolo la radice è `/confess/dicarlo`
(default); per eseguire altrove basta impostare una variabile d'ambiente prima di
avviare Jupyter:

```bash
export DECADAL_DATA_ROOT=/percorso/locale/ai/dati
jupyter lab
```

I dataset grezzi (esperimenti, rianalisi ERA5, boundary conditions) **non** sono
inclusi nel repository per dimensione: risiedono sul cluster sotto la radice sopra.
