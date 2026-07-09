"""Percorsi centralizzati per la pipeline decadale.

Tutti i percorsi derivano da un'unica radice, sovrascrivibile con la variabile
d'ambiente DECADAL_DATA_ROOT. Sul cluster di calcolo la radice e' /confess/dicarlo;
su un'altra macchina basta esportare la variabile prima di lanciare Jupyter::

    export DECADAL_DATA_ROOT=/percorso/locale/ai/dati

Nei notebook:  from config import POST_DATA, WORK_DIR, FIG_DIR, ...
"""
import os
from pathlib import Path

# Radice dei dati (default: cluster CONFESS). Sovrascrivibile via ambiente.
DATA_ROOT = Path(os.environ.get("DECADAL_DATA_ROOT", "/confess/dicarlo"))

# --- Dati in ingresso ---
CONFESS_DATA = DATA_ROOT / "00-dataset" / "02-confess" / "01-confess_data"   # dati grezzi degli esperimenti
BC_DATA      = DATA_ROOT / "00-dataset" / "02-confess" / "03-bc"             # boundary conditions (land/sea)
ERA5_ROOT    = DATA_ROOT / "00-dataset" / "01-reanalysis" / "01-era5"        # rianalisi ERA5
POST_DATA    = DATA_ROOT / "02-confess" / "05-confess-post-data"             # dati post-processati

# --- Dati calcolati e output ---
WORK_DIR     = DATA_ROOT / "cartella_ordinata"        # .nc intermedi prodotti dai notebook di calcolo
FIG_DIR      = DATA_ROOT / "figure-decadale"          # PNG dei notebook di plotting
FIG_DIR_2025 = DATA_ROOT / "02-confess" / "figures_2025"   # output figure (revisione 2025)
