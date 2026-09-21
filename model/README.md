# Model code (MRMP-N)

Scripts of the MRMP-N and of all the validations reported in the manuscript and in the dissertation (PPGCTA, UFSB/IFBA). This folder is self-contained: the scripts read `Dados_base_nova_2001_2025.xlsx` from the same directory and write to `resultados_2001_2025/`.

## Arquivos

| Arquivo | O que é |
|---|---|
| `Modelo_PSN.py` | Script principal: Ridge polinomial (grau 2), GroupKFold por ano, TimeSeriesSplit, Ljung-Box, Y-randomization. |
| `Analise_Biomas_PSN.py` | Análise ENSO: Kruskal-Wallis por fase, Pearson ONI×variáveis, OLS ONI→variáveis, boxplots e dispersões. |
| `base_final_2001_2025_plan1_excel_ptbr.csv` | Base bruta (2001–2025, separador `;`, decimal `,`). |
| `Dados_base_nova_2001_2025.xlsx` | Base convertida (gerada automaticamente pelos scripts a partir do CSV), já com a coluna `Enso`. |
| `enso_noaa.py` | Classificação oficial das fases ENSO (NOAA/CPC: ONI ≥ \|0,5\| por ≥ 5 trimestres consecutivos). Usado pelos dois scripts. |
| `oni_noaa_cpc.txt` | Tabela ONI completa da NOAA (1950–presente), fonte das fases nas bordas do período. Atualizar com `curl -o oni_noaa_cpc.txt https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt`. |
| `Selecao_Variaveis_PSN.py` | Busca exaustiva das 10 combinações C(5,3) de variáveis ambientais por bioma, mesmo pipeline do modelo. |
| `consolidar_resultados.py` | Junta logs e CSVs em `numeros_extra.json` / `enso_resumo.json` (usados pelo script do Word). |
| `Analise_ENSO_Anomalias_PSN.py` | ENSO em anomalias mensais (fase oficial): Tabela 6, mediação, defasagem, Figura 14, Tabelas A4/A5. |
| `Analise_Interanual_PSN.py` | Variabilidade interanual (6.6): soma anual 2001–2024, CV, Sen/Mann-Kendall, Figura 15, Tabela 7. |
| `Sensibilidade_Winsor_PSN.py` | Sensibilidade ao percentil de winsorização (Tabela A6). |
| `Robustez_Selecao_PSN.py` | Atalho que roda só o bloco de robustez do `Modelo_PSN.py` (função `rodar_robustez`): seleção de grau e de variáveis sob GroupKFold/TimeSeriesSplit, estabilidade partição a partição, importância por permutação por bloco temporal e nulos temporais (Tabelas A7–A10; `resultados_2001_2025/robustez/`). |
| `coletar_resultados.py` | Lê os logs de uma rodada completa e sequencial e grava resumo_geral, VIF, graus, JSONs e figuras por bioma. |
| `reproduzir_tudo.sh` | Roda a cadeia inteira, da base bruta ao Word/PDF. |
| `Figuras_Dissertacao.py` | Monta as figuras da dissertação em `figuras_dissertacao/` a partir de `resultados_2001_2025/`. |
| `resultados_2001_2025/` | Resultados da rodada com a base 2001–2025 (ver `RESULTADOS.md`). |

## Como rodar

```bash
cd model
pip install pandas numpy scipy statsmodels scikit-learn matplotlib openpyxl
BIOMA_ATIVO=TODOS python3 Modelo_PSN.py     # MA, CE e CA em paralelo
python3 Analise_Biomas_PSN.py               # análise ENSO
```

Para regenerar a dissertação revisada a partir dos resultados:

```bash
python3 Selecao_Variaveis_PSN.py      # opcional, lento (~15 min)
python3 consolidar_resultados.py
python3 Ablacao_Sazonalidade_PSN.py     # Tabela A3 e R² do ciclo anual
python3 Analise_ENSO_Anomalias_PSN.py  # seção 6.5 (anomalias, mediação, defasagem)
python3 Figuras_Dissertacao.py
python3 banca/aplicar_revisao_docx.py # gera banca/Trabalho_revisado_2001_2025.docx
banca/render_e_paginas.sh             # PDF via LibreOffice + páginas das listas (rodar o script de novo depois)
```

As saídas vão para `saidas_figuras/` (uma subpasta por bioma no modelo).
Essa pasta não é versionada aqui e nada dela é publicado no site do
Controle-Gás: o workflow de Pages copia apenas `index.html`, `manifest.json`
e os ícones.

## Reprodução completa (o que a banca ou a orientação precisa saber)

O código oficial do modelo é `Modelo_PSN.py` desta pasta, com `GRAU_MODELO = 2`,
Mata Atlântica = EV + TST + WAI, Cerrado = EV + PRE + WAI, Caatinga = EV + PRE + TST
(mais SAZsin e SAZcos), winsorização no percentil 3 só no treino, RepeatedKFold 5 × 30,
GroupKFold por ano, TimeSeriesSplit, Ljung-Box, VIF, Y-randomization com 100 permutações e, ao final, o bloco de
robustez (seleção de grau e de variáveis sob validação temporal, estabilidade da seleção, importância por
permutação e nulos temporais: Tabelas A7–A10).
As saídas têm nomes fixos (cada rodada sobrescreve; `NUMERAR_SAIDAS = True` volta à numeração).
Os toggles aceitam variáveis de ambiente: `BIOMA_ATIVO`, `GRAU_MODELO`, `TESTAR_GRAUS`,
`RODAR_YRANDOMIZATION`, `RODAR_ROBUSTEZ`, `RODAR_EM_PARALELO`.

Cadeia completa, da base bruta ao PDF (≈ 1 h; `bash reproduzir_tudo.sh` faz tudo):

| Etapa | Script | Produz |
|---|---|---|
| 0 | `../npp_modis/npp_modis_gee.py` (Google Earth Engine) | `base_final_2001_2025_plan1_excel_ptbr.csv` (MODIS 6.1 + IMERG V07) |
| 1 | `Modelo_PSN.py` (sequencial, graus 1–5, Y-rand) | `saidas_figuras/`, log com todas as métricas |
| 2 | `Analise_Biomas_PSN.py` | ENSO com fases oficiais, dispersão, boxplots |
| 3 | `coletar_resultados.py` | `resumo_geral.csv`, `vif_por_bioma.csv`, `selecao_grau.csv`, `numeros_extra.json`, `enso_resumo.json` |
| 4 | `Selecao_Variaveis_PSN.py`, `Ablacao_Sazonalidade_PSN.py` | Tabelas A2 e A3 |
| 5 | `Analise_ENSO_Anomalias_PSN.py`, `Analise_Interanual_PSN.py` | Seções 6.5 e 6.6 (Tabelas 6, 7, A4, A5; Figuras 14 e 15) |
| 5b | `Sensibilidade_Winsor_PSN.py` | Tabela A6 (sensibilidade à winsorização); as Tabelas A7–A10 saem do próprio `Modelo_PSN.py` (bloco de robustez, copiado pelo coletor) |
| 6 | `Figuras_Dissertacao.py` | todas as figuras em `figuras_dissertacao/` |
| 7 | `banca/aplicar_revisao_docx.py` + `banca/render_e_paginas.sh` (duas vezes) | `banca/Trabalho_revisado_2001_2025.docx` e `.pdf` |

Cada número do texto, das tabelas e das legendas é preenchido pelo script do Word a partir
desses arquivos; nada é digitado à mão. Uma rodada sequencial completa feita em 21/09/2026
reproduziu as Tabelas 2 e 3 exatamente (ver `resultados_2001_2025/RESULTADOS.md`, seção 11).

## Versão citada na dissertação

O texto cita o repositório pela URL fixa `https://github.com/heroslore/Controle-Gas` (pastas
`dissertacao_PSN` e `npp_modis`). A versão exata que gerou o documento entregue à banca é o
commit `9f6227b` da branch `claude/tender-hawking-f81psp` (arquivo
`banca/Trabalho_revisado_2001_2025.docx`). Para congelá-la com um nome estável, sem mudar o
link citado, crie a tag anotada e publique-a (o ambiente automatizado só pode enviar a branch):

```bash
git fetch origin claude/tender-hawking-f81psp
git tag -a dissertacao-v10 9f6227b -m "Dissertação MRMP-N: versão 10 entregue à banca"
git push origin dissertacao-v10
```

Em seguida, no GitHub, transforme a tag em *release* anexando o `.docx` e o `.pdf`: a banca
passa a ter um ponto de acesso único e imutável, e a branch pode continuar a evoluir.
