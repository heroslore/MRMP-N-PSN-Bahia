# Model code (MRMP-N)

Scripts of the MRMP-N and of all the validations reported in the manuscript and in the dissertation (PPGCTA, UFSB/IFBA). This folder is self-contained: the scripts read `Dados_base_nova_2001_2025.xlsx` from the same directory and write to `resultados_2001_2025/`.

## Arquivos

| Arquivo | O que é |
|---|---|
| `Modelo_PSN.py` | Script principal: Ridge polinomial (grau 2), GroupKFold por ano, TimeSeriesSplit, Ljung-Box, Y-randomization. |
| `Analise_Biomas_PSN.py` | Análise ENSO: Kruskal-Wallis por fase, Pearson ONI×variáveis, OLS ONI→variáveis, boxplots e dispersões. |
| `base_final_2001_2025_plan1_excel_ptbr.csv` | Base bruta (2001 a 2025, separador `;`, decimal `,`). |
| `Dados_base_nova_2001_2025.xlsx` | Base convertida (gerada automaticamente pelos scripts a partir do CSV), já com a coluna `Enso`. |
| `enso_noaa.py` | Classificação oficial das fases ENSO (NOAA/CPC: ONI ≥ \|0,5\| por ≥ 5 trimestres consecutivos). Usado pelos dois scripts. |
| `oni_noaa_cpc.txt` | Tabela ONI completa da NOAA (1950 até o presente), fonte das fases nas bordas do período. Atualizar com `curl -o oni_noaa_cpc.txt https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt`. |
| `Selecao_Variaveis_PSN.py` | Busca exaustiva das 10 combinações C(5,3) de variáveis ambientais por bioma, mesmo pipeline do modelo. |
| `consolidar_resultados.py` | Versão anterior do coletor: junta logs e CSVs em `numeros_extra.json` e `enso_resumo.json` e incorpora o log de uma rodada separada da Mata Atlântica (`--log-ma`). Substituído pelo `coletar_resultados.py` na cadeia atual. |
| `Analise_ENSO_Anomalias_PSN.py` | ENSO em anomalias mensais (fase oficial): Tabela 6, mediação, defasagem, Figura 14, Tabelas A4/A5. |
| `Analise_Interanual_PSN.py` | Variabilidade interanual (6.6): soma anual 2001 a 2024, CV, Sen/Mann-Kendall, Figura 15, Tabela 7. |
| `Sensibilidade_Winsor_PSN.py` | Sensibilidade ao percentil de winsorização (Tabela A6). |
| `Robustez_Selecao_PSN.py` | Atalho que roda só o bloco de robustez do `Modelo_PSN.py` (função `rodar_robustez`): seleção de grau e de variáveis sob GroupKFold/TimeSeriesSplit, estabilidade partição a partição, importância por permutação por bloco temporal e nulos temporais (Tabelas A7 a A10; `resultados_2001_2025/robustez/`). |
| `coletar_resultados.py` | Lê os logs de uma rodada completa e sequencial e grava resumo_geral, VIF, graus, JSONs e figuras por bioma. |
| `reproduzir_tudo.sh` | Roda a cadeia inteira, da base bruta às tabelas e figuras. |
| `Figuras_Dissertacao.py` | Monta as figuras da dissertação em `figuras_dissertacao/` a partir de `resultados_2001_2025/`; com `FIG_IDIOMA=EN`, grava as Figs. 2 a 5 do artigo do modelo em inglês em `figuras_artigo1/en/` (em PT, copia as mesmas figuras para `figuras_artigo1/pt/`). |
| `resultados_2001_2025/` | Resultados da rodada com a base 2001 a 2025 (ver `RESULTADOS.md`). |
| `dissertacao/` | Versão da dissertação da análise do ENSO (`Analise_ENSO_Anomalias_PSN_dissertacao.py`), que gera as Tabelas A11 e A12 da dissertação; ver `dissertacao/README.md`. |
| `sensibilidade/` | Scripts das análises complementares de sensibilidade citadas na dissertação (`auditoria_comum.py`, `t01` a `t06`); ver `sensibilidade/README.md`. |
| `sensibilidade_qc/` | Sensibilidade da TST ao filtro QC_Day do MOD11A2; ver `sensibilidade_qc/LEIA-ME.md`. |
| `Figura_Area_Estudo.py` | Figura 1 do artigo do ENSO (mapa dos biomas e climatologia de precipitação), em `figuras_artigo2/`. |


## Como rodar

```bash
cd model
pip install -r ../requirements.txt
BIOMA_ATIVO=TODOS python3 Modelo_PSN.py     # MA, CE e CA em paralelo
python3 Analise_Biomas_PSN.py               # análise ENSO
```

O `Modelo_PSN.py` grava as saídas em `saidas_figuras/` (uma subpasta por bioma e `robustez/`), pasta que não é versionada.
Quem copia essas saídas para `resultados_2001_2025/` é o `coletar_resultados.py`, a partir do log de uma rodada sequencial
(`RODAR_EM_PARALELO=0`).

## Reprodução completa

O código do modelo é `Modelo_PSN.py` desta pasta, com `GRAU_MODELO = 2`,
Mata Atlântica = EV + TST + WAI, Cerrado = EV + PRE + WAI, Caatinga = EV + PRE + TST
(mais SAZsin e SAZcos), winsorização no percentil 3 só no treino, RepeatedKFold 5 × 30,
GroupKFold por ano, TimeSeriesSplit, Ljung-Box, VIF, Y-randomization com 100 permutações e, ao final, o bloco de
robustez (seleção de grau e de variáveis sob validação temporal, estabilidade da seleção, importância por
permutação e nulos temporais: Tabelas A7 a A10).
As saídas têm nomes fixos (cada rodada sobrescreve; `NUMERAR_SAIDAS = True` volta à numeração).
Os toggles aceitam variáveis de ambiente: `BIOMA_ATIVO`, `GRAU_MODELO`, `TESTAR_GRAUS`,
`RODAR_YRANDOMIZATION`, `RODAR_ROBUSTEZ`, `RODAR_EM_PARALELO`.

`bash reproduzir_tudo.sh` roda a cadeia inteira a partir da base bruta (cerca de 1 h de CPU; o modelo e as análises
complementares são as partes lentas):

| Etapa | Script | Produz |
|---|---|---|
| 1 | `Modelo_PSN.py` (sequencial, graus 1 a 5, Y-randomization, robustez) | `saidas_figuras/` e `resultados_2001_2025/log_modelo_psn_completo.txt` |
| 2 | `Analise_Biomas_PSN.py` | ENSO com fases oficiais, dispersão, boxplots; `resultados_2001_2025/log_analise_biomas.txt` |
| 3 | `coletar_resultados.py` | `resumo_geral.csv`, `vif_por_bioma.csv`, `selecao_grau.csv`, `numeros_extra.json`, `enso_resumo.json`, figuras e coeficientes por bioma, `robustez/` |
| 4 | `Selecao_Variaveis_PSN.py`, `Ablacao_Sazonalidade_PSN.py` | Tabelas A2 e A3 |
| 5 | `Analise_ENSO_Anomalias_PSN.py`, `Analise_Interanual_PSN.py` | Seções 6.5 e 6.6 (Tabelas 6, 7, A4, A5; Figuras 14 e 15) |
| 6 | `Sensibilidade_Winsor_PSN.py` | Tabela A6 (sensibilidade à winsorização) |
| 7 | `Figuras_Dissertacao.py` (e `FIG_IDIOMA=EN`) | figuras da dissertação em `figuras_dissertacao/` e Figs. 2 a 5 do artigo do modelo em `figuras_artigo1/` |
| 8 | `Bootstrap_Blocos_ENSO.py`, `Teste_Pareado_Corrigido.py`, `Figuras_Artigo2.py`, `Figura_Area_Estudo.py` | unidade amostral do ENSO com Benjamini-Hochberg, teste pareado com correção de Nadeau-Bengio, figuras do artigo do ENSO em `figuras_artigo2/` |
| 9 | `dissertacao/Analise_ENSO_Anomalias_PSN_dissertacao.py` | Tabelas A11 e A12 da dissertação, em `dissertacao/resultados/` |
| 10 | `sensibilidade/t01` a `t06` | análises complementares de sensibilidade, em `sensibilidade/resultados/` |
| 11 | `sensibilidade_qc/Recalcular_297.py` | sensibilidade da TST ao filtro QC_Day na série de 297 meses |

Outros scripts que geram arquivos versionados e não fazem parte do `reproduzir_tudo.sh`:
`Estratificacao_Sazonal_ENSO.py` e `Memoria_Hidrica.py` (CSVs em `resultados_2001_2025/enso_anomalias/` e
Figuras 3 e 4 de `figuras_artigo2/`), `Figura_Conceitual.py` (Figura 5 do artigo do ENSO), `Recolorir_Mapa.py`
(mapa dos biomas em paleta acessível, a partir de `../manuscript/figs/fig1_map.png`) e `Robustez_Selecao_PSN.py`
(só o bloco de robustez; console em `resultados_2001_2025/log_robustez.txt`).

A base bruta `base_final_2001_2025_plan1_excel_ptbr.csv` foi extraída no Google Earth Engine (MODIS 6.1 e IMERG V07);
a rotina de extração não faz parte deste repositório, de modo que a série mensal é o ponto de partida da reprodução.

As figuras e os resultados de `resultados_2001_2025/`, `figuras_dissertacao/` e `figuras_artigo2/` foram regenerados com o
código desta versão (os valores numéricos são os mesmos da rodada sequencial de 21/09/2026 descrita em
`resultados_2001_2025/RESULTADOS.md`; mudam só os rótulos das figuras). Os logs `log_modelo_psn.txt` e
`log_modelo_psn_MA_EV_TST_WAI.txt` são registros de rodadas anteriores e foram mantidos como estavam, exceto pelo prefixo
dos caminhos de pasta, substituído por `.`.
