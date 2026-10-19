# MRMP-N: net photosynthesis of the Atlantic Forest, Cerrado and Caatinga of Bahia, Brazil (2001 to 2025)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22883915.svg)](https://doi.org/10.5281/zenodo.22883915)

Manuscripts, model code and monthly database of the study, which is also the subject of the first author's master's dissertation (PPGCTA, UFSB/IFBA). The paper below is in preparation (not yet submitted):

> Oliveira, Y.A.S., Benfica, N.S., Zanchi, F.B. *Regularized polynomial regression of MODIS net photosynthesis across a tropical hydroclimatic gradient: Atlantic Forest, Cerrado and Caatinga in northeastern Brazil (2001-2025).* In preparation.

The MRMP-N is a Multiple Polynomial Regression Model of order N estimated by Ridge regression. It predicts monthly MODIS net photosynthesis (PSN, MOD17A2HGF) from evapotranspiration, precipitation, land surface temperature, a water availability index and two harmonic components of the annual cycle, with polynomial degree and predictor set selected by repeated cross-validation and verified by year-grouped validation, chronological blocks, permutation importance and null models that preserve the temporal structure of the series.

## Repository layout

| Folder | Content |
|---|---|
| `manuscript/` | LaTeX sources in the Elsevier CAS template (`main.tex` English, `main_pt.tex` Portuguese; `*_dc.tex` two-column variants), bibliography, figures and compiled PDFs in `pdf/`. |
| `manuscript2/` | Drafts of the companion paper on the ENSO analysis (Portuguese). `main2d_pt.tex` is the current version: it leads with the contrasting response patterns of the three biomes along the hydroclimatic gradient, adds the study-area figure, the water-accumulation windows, the seasonal stratification with its direct contrast test, a conceptual synthesis figure and a section on monitoring implications, and keeps the sampling-unit analysis as a robustness section, with `suplementar_pt.tex` holding the full episode-level table. `main2_pt.tex`, `main2b_pt.tex` and `main2c_pt.tex` are the earlier framings, kept for reference. |
| `data/` | `psn_bahia_monthly_2001_2025.csv`, the monthly database (297 periods, 3 biomes, English headers), its data dictionary, the NOAA ONI table and the data licence (CC BY 4.0). |
| `model/` | Python scripts of the model and of all validations, the database files read by the scripts, the results of the 2001 to 2025 run (`resultados_2001_2025/`) and the dissertation figures (`figuras_dissertacao/`). See `model/README.md`. |
| `model/dissertacao/` | Dissertation version of the ENSO anomaly analysis (`Analise_ENSO_Anomalias_PSN_dissertacao.py`), which produces Tables A11 and A12 of the dissertation; the ENSO paper uses `model/Bootstrap_Blocos_ENSO.py` and `model/Estratificacao_Sazonal_ENSO.py` (same effects, confidence intervals and p-values from different resampling routines). See its `README.md`. |
| `model/sensibilidade/` | Scripts of the complementary sensitivity analyses cited in the dissertation: seasonal baselines and stricter temporal schemes, nested selection, Y-randomization with 999 permutations, year-partition sensitivity of the GroupKFold, extended alpha grid and year-by-year bias, refits without the MODIS-derived predictors. See its `README.md`. |
| `model/sensibilidade_qc/` | Sensitivity of the land surface temperature to the MOD11A2 QC_Day filter. See `LEIA-ME.md`. |

## Quick start

```bash
pip install -r requirements.txt
cd model
BIOMA_ATIVO=TODOS python3 Modelo_PSN.py      # fits and validates the model for the three biomes
python3 Robustez_Selecao_PSN.py              # selection under temporal schemes, permutation importance, temporal nulls
python3 Selecao_Variaveis_PSN.py             # exhaustive search of the 10 predictor combinations (slow)
cd dissertacao && python3 Analise_ENSO_Anomalias_PSN_dissertacao.py && cd ..   # dissertation Tables A11 and A12
cd sensibilidade && python3 t01_diagnosticos.py && cd ..                         # one of the sensitivity analyses
```

`bash model/reproduzir_tudo.sh` runs the whole chain.

To compile the manuscript:

```bash
cd manuscript
pdflatex main && bibtex main && pdflatex main && pdflatex main
```

## Data sources

MODIS Collection 6.1 products MOD17A2HGF (PSN), MOD16A2GF (EV, PET), MOD11A2 (LST) and MCD64A1 (burned area); GPM IMERG Final Run V07 (precipitation); NOAA/CPC Oceanic Niño Index. All products were processed in Google Earth Engine and aggregated by biome using the IBGE 1:250,000 biome map clipped to Bahia State. Details are given in the manuscript and in `data/data_dictionary.md`.

## Licences

Code: MIT (see `LICENSE`). Data: CC BY 4.0 (see `data/LICENSE`). Manuscript text and figures: © the authors; all rights reserved until publication.

## Funding

This study was financed in part by the Coordenação de Aperfeiçoamento de Pessoal de Nível Superior - Brasil (CAPES) - Finance Code 001.

## How to cite

This repository is permanently archived on Zenodo.

- **Concept DOI (cite this one):** [10.5281/zenodo.22883915](https://doi.org/10.5281/zenodo.22883915), which always resolves to the most recent version.
- Each release is also archived under its own version DOI, listed on the Zenodo record under *Versions*.

Machine-readable metadata is in `CITATION.cff` (GitHub renders a "Cite this repository" button from it) and in `.zenodo.json`.
