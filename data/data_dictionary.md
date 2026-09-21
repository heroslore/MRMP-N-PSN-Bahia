# Data dictionary: `psn_bahia_monthly_2001_2025.csv`

One row per period and 297 rows (January 2001 to September 2025). Each period is a fixed
window of four consecutive 8-day MODIS composites (about 32 days) assigned to the calendar
month that contains most of its days; see the manuscript, Materials and methods.

| Column | Description | Unit | Source |
|---|---|---|---|
| `year`, `month` | Reference calendar year and month of the window | - | - |
| `ONI` | Oceanic Niño Index of the 3-month season centred on the month | °C | NOAA/CPC (`oni_noaa_cpc.txt`) |
| `ENSO_phase` | El Niño / La Niña / Neutral, official NOAA criterion (ONI ≥ +0.5 or ≤ −0.5 °C for ≥ 5 consecutive overlapping seasons) | - | derived (`model/enso_noaa.py`) |
| `PSN_<biome>` | Net photosynthesis, sum of the four composites | gC m⁻² per period | MOD17A2HGF C6.1 |
| `EV_<biome>` | Actual evapotranspiration, sum | mm per period | MOD16A2GF C6.1 |
| `PET_<biome>` | Potential evapotranspiration, sum | mm per period | MOD16A2GF C6.1 |
| `WAI_<biome>` | Water availability index = EV / PET | dimensionless | derived |
| `PRE_<biome>` | Accumulated precipitation | mm per period | GPM IMERG Final Run V07 |
| `TST_<biome>` | Daytime land surface temperature, mean | °C | MOD11A2 C6.1 |
| `BURN_<biome>` | Burned area (burned pixels × 25 ha) | ha | MCD64A1 C6.1 |

Biome suffixes: `MA` = Atlantic Forest (Mata Atlântica), `CE` = Cerrado, `CA` = Caatinga.
Spatial means were taken over the valid pixels of each biome, delimited by the IBGE 1:250,000
biome map clipped to the Bahia State boundary, in Google Earth Engine.

`model/base_final_2001_2025_plan1_excel_ptbr.csv` is the raw extraction (Portuguese headers,
`;` separator, `,` decimal) and `model/Dados_base_nova_2001_2025.xlsx` is the file read by the
scripts; the CSV in this folder is the same data with English headers and full precision.
