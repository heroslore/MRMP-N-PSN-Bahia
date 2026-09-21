# Resultados — base 2001–2025 (n = 297)

Rodada em 12/09/2026 com `Modelo_PSN.py`. **Atualização (mesmo dia):** a busca exaustiva de combinações (`Selecao_Variaveis_PSN.py`, seção 9) mostrou que na Mata Atlântica o conjunto EV + TST + WAI supera EV + PRE + TST; a MA foi rodada de novo com esse conjunto e as tabelas abaixo já refletem isso (CE e CA inalterados). Rodada original com `Modelo_PSN.py` (BIOMA_ATIVO=TODOS, grau 2,
RODAR_YRANDOMIZATION=True, TESTAR_GRAUS=False) e `Analise_Biomas_PSN.py`.
Base: `base_final_2001_2025_plan1_excel_ptbr.csv`, 300 linhas; removidas
out/nov/dez de 2025 (precipitação ainda não publicada) → 297 linhas.
Logs completos em `log_modelo_psn.txt` e `log_analise_biomas.txt`.

## 1. Desempenho do modelo (RepeatedKFold, 150 rodadas)

| Bioma | R² treino | R² teste | IC95% teste | Gap (pp) | RMSE | MAE | GroupKFold ano | TimeSeriesSplit | Queda TS (pp) |
|---|---|---|---|---|---|---|---|---|---|
| MA (EV+TST+WAI) | 78,84 | **73,65** | 60,06–83,40 | 5,19 | 7,97 | 6,29 | 73,73 | 66,98 | 6,66 |
| CE | 97,25 | **96,66** | 95,40–97,71 | 0,59 | 6,74 | 5,43 | 96,38 | 95,73 | 0,93 |
| CA | 95,19 | **94,16** | 90,79–96,25 | 1,03 | 6,69 | 5,08 | 93,81 | 91,20 | 2,97 |

CE e CA conferem com a tabela do handoff (96,66 / 94,16). MA com o conjunto antigo (EV+PRE+TST) dava 70,96 (gap 5,14), igual ao handoff; com EV+TST+WAI passou a 73,65. R² out-of-fold da Figura 3: MA 73,85 %, CE 96,88 %, CA 94,65 %.

## 2. Y-randomization (100 permutações, validação cruzada) — NOVO

| Bioma | R² CV original | R² CV permutados (média ± dp) | Δ (pp) | p empírico | Status |
|---|---|---|---|---|---|
| MA | 73,81 % | -10,61 ± 3,09 % | 84,42 | 0,0099 | APROVADO |
| CE | 96,66 % | −10,88 ± 2,64 % | 107,54 | 0,0099 | APROVADO |
| CA | 94,16 % | −13,04 ± 3,07 % | 107,20 | 0,0099 | APROVADO |

p = 0,0099 = 1/(100+1): o modelo original superou todas as 100 permutações.
Critério operacional de 60 pp atingido nos três biomas. Figuras:
`MA/yrandomization_MA.png`, `CE/…`, `CA/…`.

## 3. VIF completo por bioma — NOVO (`vif_por_bioma.csv`)

| Variável | MA (EV+TST+WAI) | CE | CA |
|---|---|---|---|
| EV | 2,65 | **36,99** | 5,02 |
| PRE | — | 4,02 | 1,68 |
| TST | 4,23 | — | 3,90 |
| WAI | 3,54 | **31,61** | — |
| saz_sin | 1,65 | 4,49 | 2,63 |
| saz_cos | 3,92 | 3,59 | 3,76 |

MA e CA: todos abaixo de 5 (CA EV = 5,02, no limite); na MA, EV e WAI **não** são colineares (VIF 2,65 e 3,54), ao contrário do Cerrado. CE: EV_CE e WAI_CE
com colinearidade grave, como já descrito no handoff.

## 4. Diagnóstico dos resíduos

| Bioma | Shapiro-Wilk p | Normais | Ljung-Box(12) p | ACF lag1 | Alpha Ridge |
|---|---|---|---|---|---|
| MA | 0,6299 | Sim | 0,0000 | +0,505 | 0,71 |
| CE | 0,4232 | Sim | 0,0052 | +0,201 | 0,10 |
| CA | 0,0015 | Não | 0,0000 | +0,420 | 0,10 |

Autocorrelação residual nos três biomas, idêntica ao achado do handoff.
CA é o único com resíduos não normais (caudas mais pesadas no QQ-plot). Com EV+TST+WAI, a ACF(1) da MA subiu para 0,50 (era 0,30 com EV+PRE+TST).

## 5. Coeficientes Ridge (Figura 5) — NOVO

Arquivos: `MA/coeficientes_ridge_MA.csv`, `CE/…`, `CA/…` (20 termos cada:
5 lineares, 5 quadráticos, 10 interações). Importância relativa dos termos
lineares (|coef| / soma |coef| lineares):

| Termo | MA | CE | CA |
|---|---|---|---|
| EV | +22,79 | +56,48 | +33,12 |
| PRE | — | −12,18 | −4,54 |
| TST | −11,11 | — | −4,08 |
| WAI | −16,66 | −6,57 | — |
| saz_sin | −8,88 | −7,73 | −5,94 |
| saz_cos | −8,84 | −11,62 | −6,49 |

Importância relativa (soma de |coef| de todos os termos que envolvem a variável, normalizada a 100 %), usada na Figura 8 da dissertação: MA — EV 26.9 %, saz_cos 20.9 %, WAI 19.0 %, saz_sin 18.8 %, TST 14.5 %; CE — EV 37.4 %, saz_sin 19.5 %, WAI 17.9 %, saz_cos 16.1 %, PRE 9.1 %; CA — EV 42.8 %, saz_cos 18.2 %, TST 15.2 %, saz_sin 13.6 %, PRE 10.2 %.

Maiores interações: CE `WAI_CE×saz_sin` +16,50 e `EV_CE×WAI_CE` −13,04
(reflexo da colinearidade EV×WAI); MA `saz_sin×saz_cos` +8,85 e `EV_MA×saz_cos` +3,97; CA
`EV_CA×saz_cos` +5,52.

## 6. Análise ENSO (classificação OFICIAL NOAA/CPC) — NOVO

Fase por mês definida pela regra operacional da NOAA: ONI >= +0,5 (El Niño)
ou <= −0,5 (La Niña) por pelo menos 5 trimestres móveis consecutivos, usando
a tabela completa oni.ascii.txt (cópia em `oni_noaa_cpc.txt`) para as bordas.
O ONI do CSV é idêntico ao da NOAA nos 300 meses. Em relação ao limiar
simples, só **um mês mudou**: outubro/2016 (ONI −0,51, sequência de 1 mês)
passou de La Niña para Neutro. Distribuição: Neutro 149, El Niño 76,
La Niña 72 meses.

Kruskal-Wallis por fase (médias La Niña | Neutro | El Niño):

| Var | MA | CE | CA |
|---|---|---|---|
| PSN | 135,3 \| 136,0 \| 128,6 — p=0,065 ns | 104,2 \| 93,8 \| 92,9 — p=0,094 ns | 94,2 \| 87,0 \| 79,7 — **p=0,010 \*** |
| EV | p=0,006 ** | p=0,0003 *** | p=0,0015 ** |
| PRE | p=0,349 ns | p<0,0001 *** | p=0,0003 *** |
| TST | p<0,0001 *** | p=0,003 ** | p<0,0001 *** |

Pearson ONI × PSN (não depende da fase): MA r=−0,240 (p<0,0001), CE
r=−0,071 (ns), CA r=−0,155 (p=0,007). ONI × TST positivo em MA (0,175 **)
e CA (0,134 *).

OLS ONI → variáveis: R² máximo 3,05 % (TST_MA); EV_MA 2,00 %; todos os
demais abaixo de 2 %. Sustenta a afirmação de que o ONI explica menos de
~3 % da variância das variáveis climáticas.

Regressão simples PSN × variáveis (R² ajustado): EV explica 74,8 % (CE) e
76,2 % (CA) mas só 12,5 % em MA; WAI 81,5 % em CA; TST 59,2 % (CE) e 54,2 %
(CA), 11,2 % em MA. Área queimada (log) negativa nos três biomas.

Observação: na tabela ONI atual da NOAA (período-base atualizado), o evento
frio de 2016-17 não alcança 5 trimestres com ONI <= −0,5 (set −0,42,
out −0,51, nov −0,49), por isso não conta como La Niña aqui. Se a
dissertação citar a lista histórica de episódios da NOAA, vale registrar
essa nuance.

## 7. Observações sobre os scripts enviados

- O `Modelo_PSN.py` enviado ainda **numera** as saídas (`_1`, `_2`…) via
  `proximo_nome`, e está com `TESTAR_GRAUS=False`. O handoff descreve uma
  versão com nomes fixos e `TESTAR_GRAUS=True`. Os resultados não mudam,
  mas a versão no repositório não é a "final" descrita no handoff.
- Corrigido: os dois scripts agora gravam o mesmo
  `Dados_base_nova_2001_2025.xlsx`, com a coluna `Enso` oficial (módulo
  `enso_noaa.py`).

## 8. Seleção do grau polinomial (1 a 5) — NOVO

Rodada separada com `TESTAR_GRAUS=True` (150 partições por grau). R² de
teste médio e diferença treino–teste em pp:

| Grau | Termos | MA | CE | CA |
|---|---|---|---|---|
| 1 | 5 | 65,2 (1,5) | 90,1 (1,2) | 92,9 (0,6) |
| **2** | 20 | **73,6 (5,2)** | **96,7 (0,6)** | **94,2 (1,0)** |
| 3 | 55 | 74,4 (6,9) | 96,4 (1,0) | 93,4 (2,4) |
| 4 | 125 | 68,4 (12,2) | 96,1 (1,6) | 63,8 (31,4) |
| 5 | 251 | 65,1 (16,8) | 93,7 (3,9) | 52,3 (42,8) |

Grau 2 é o ótimo nos três biomas pelo critério R² teste − 0,5 × gap (MA com EV+TST+WAI: o grau 3 tem R² teste 0,8 pp maior, mas gap 1,7 pp maior). Valores de MA já com EV+TST+WAI. Figuras:
`MA/selecao_grau_polinomial_MA.png`, `CE/…`, `CA/…`.

## 9. Busca exaustiva de combinações de variáveis (C(5,3) = 10 por bioma) — NOVO

`Selecao_Variaveis_PSN.py`, mesmo pipeline do modelo (RepeatedKFold 5×30, grau 2). Arquivo: `selecao_variaveis.csv`. Três melhores por bioma (R² teste, gap):

- **MA**: EV + TST + WAI (73,6 %, 5,2 pp); EV + PRE + TST (71,0 %, 5,1 pp); EV + PRE + WAI (69,0 %, 5,4 pp)
- **CE**: EV + PRE + WAI (96,7 %, 0,6 pp); EV + PRE + TST (96,3 %, 0,7 pp); EV + PRE + BURNlog (95,8 %, 1,2 pp)
- **CA**: EV + PRE + TST (94,2 %, 1,0 pp); EV + PRE + WAI (94,0 %, 1,1 pp); PRE + TST + WAI (93,9 %, 1,1 pp)

Conclusão: CE (EV+PRE+WAI) e CA (EV+PRE+TST) confirmam os conjuntos da qualificação. Na **MA o conjunto ótimo mudou** para EV+TST+WAI (+2,7 pp de R² de teste, mesmo gap). O `Modelo_PSN.py` foi atualizado para esse conjunto e a MA foi rodada de novo (`log_modelo_psn_MA_EV_TST_WAI.txt`). Para voltar ao conjunto antigo, basta trocar a linha `'x'` da MA em `config_biomas`.

## 10. ENSO em anomalias mensais, mediação e defasagem — NOVO (21/09/2026)

Script `Analise_ENSO_Anomalias_PSN.py` (fase oficial NOAA da coluna `Enso`: 149 neutros, 76 El Niño, 72 La Niña).
Saídas em `enso_anomalias/` e `enso_anomalias.json`; Tabela 6, Figura 14, Tabelas A4 (ampliada) e A5 da dissertação.

- Anomalia = valor − média do mês do calendário; Δ em % da média; testes (Kruskal-Wallis, Mann-Whitney,
  Fligner-Killeen, qui-quadrado nos extremos P10/P90) sobre as anomalias nas unidades originais.
- PSN: La Niña +11,7 % no Cerrado (p < 0,001) e +10,4 % na Caatinga (p = 0,002); El Niño −5,4 % na Mata
  Atlântica (p = 0,010), com 22 % dos meses abaixo do P10 contra 5 % nos neutros (qui-quadrado p < 0,001).
- Mediação pelo modelo (Ridge, alpha médio, série completa): La Niña → +8,0 % (CE, obs. +11,7) e +9,8 % (CA,
  obs. +10,4), via EV; El Niño → −3,1 % na MA (obs. −5,4), via TST.
- Defasagem (Spearman ONI × anomalia PSN): CA máximo no mês 0 (ρ = −0,19), significativo até 4 meses; MA máximo
  aos 2 meses (−0,16), até 5; CE significativo de 1 a 12 meses (máximo −0,16 aos 10). Compósitos na Figura 14.
- Diferenças em relação à versão 6 do docx (outra sessão): fases 149/76/72 em vez de 151/76/70 (jan–fev/2001
  são La Niña pela tabela completa da NOAA); Δ mudam ≤ 0,7 pp; importâncias com/sem harmônicos recalculadas com
  o pipeline do modelo (Tabela A4).

## 11. Prova de reprodutibilidade (rodada completa e sequencial, 21/09/2026)

`BIOMA_ATIVO=TODOS RODAR_EM_PARALELO=0 TESTAR_GRAUS=1 RODAR_YRANDOMIZATION=1 python3 Modelo_PSN.py`
(log em `log_modelo_psn_completo.txt`) seguido de `coletar_resultados.py`. Comparação com os
arquivos que alimentam a dissertação (`resumo_geral.csv`, `selecao_grau.csv`, `vif_por_bioma.csv`,
`numeros_extra.json`): **diferença zero** em todas as métricas das Tabelas 2 e 3, nos 15 pares
grau × bioma, nos VIF, nos alfas médios, no R² fora da amostra e nas importâncias da Figura 8
(única diferença: `vif_max` do Cerrado gravado com 2 casas, 36,99 contra 36,990854).
Ou seja: a cadeia CSV bruto → `Modelo_PSN.py` → `coletar_resultados.py` → Word regenera a
dissertação sem nenhum número digitado à mão.

## 12. Robustez da seleção, importância por permutação e nulos temporais (`Robustez_Selecao_PSN.py`)

Rodado em 22/09/2026 em resposta à revisão do artigo (JSAES). Mesmo pipeline e semente do modelo.
Saídas em `robustez/`; Tabelas A7–A10 da dissertação.

**Grau polinomial por esquema (R² de teste, %)**

| bioma | grau | GroupKFold | TimeSeriesSplit |
|---|---|---|---|
| Caatinga | 1 | 92.78 | 91.49 |
| Caatinga | 2 | 93.81 | 91.2 |
| Caatinga | 3 | 92.76 | 90.4 |
| Caatinga | 4 | 91.32 | 86.99 |
| Caatinga | 5 | 85.16 | 86.66 |
| Cerrado | 1 | 89.99 | 88.46 |
| Cerrado | 2 | 96.38 | 95.73 |
| Cerrado | 3 | 96.17 | 93.68 |
| Cerrado | 4 | 96.14 | 93.39 |
| Cerrado | 5 | 95.77 | 92.24 |
| Mata Atlântica | 1 | 65.45 | 58.7 |
| Mata Atlântica | 2 | 73.73 | 66.98 |
| Mata Atlântica | 3 | 74.51 | 60.97 |
| Mata Atlântica | 4 | 67.38 | 56.32 |
| Mata Atlântica | 5 | 67.19 | 53.9 |

Só duas exceções ao grau 2: grau 3 na Mata Atlântica sob GroupKFold (+0,8 pp) e grau 1 na Caatinga sob TimeSeriesSplit (+0,3 pp).

**Combinações: duas melhores por esquema**

| bioma | esquema | rank | variaveis | r2_teste |
|---|---|---|---|---|
| Caatinga | GroupKFold | 1 | EV + PRE + TST | 93.81 |
| Caatinga | GroupKFold | 2 | EV + PRE + WAI | 93.75 |
| Caatinga | RepeatedKFold | 1 | EV + PRE + TST | 94.16 |
| Caatinga | RepeatedKFold | 2 | EV + PRE + WAI | 94.01 |
| Caatinga | TimeSeriesSplit | 1 | EV + PRE + WAI | 91.92 |
| Caatinga | TimeSeriesSplit | 2 | EV + PRE + TST | 91.2 |
| Cerrado | GroupKFold | 1 | EV + PRE + WAI | 96.38 |
| Cerrado | GroupKFold | 2 | EV + PRE + TST | 95.95 |
| Cerrado | RepeatedKFold | 1 | EV + PRE + WAI | 96.66 |
| Cerrado | RepeatedKFold | 2 | EV + PRE + TST | 96.27 |
| Cerrado | TimeSeriesSplit | 1 | EV + PRE + WAI | 95.73 |
| Cerrado | TimeSeriesSplit | 2 | EV + PRE + TST | 95.07 |
| Mata Atlântica | GroupKFold | 1 | EV + TST + WAI | 73.73 |
| Mata Atlântica | GroupKFold | 2 | EV + PRE + TST | 70.54 |
| Mata Atlântica | RepeatedKFold | 1 | EV + TST + WAI | 73.65 |
| Mata Atlântica | RepeatedKFold | 2 | EV + PRE + TST | 70.96 |
| Mata Atlântica | TimeSeriesSplit | 1 | EV + TST + WAI | 66.98 |
| Mata Atlântica | TimeSeriesSplit | 2 | EV + PRE + TST | 63.47 |

**Estabilidade partição a partição (RepeatedKFold): frequência com que cada combinação foi a melhor**

| bioma | variaveis | r2_teste | freq_melhor_pct |
|---|---|---|---|
| Mata Atlântica | EV + TST + WAI | 73.6 | 66.0 |
| Mata Atlântica | EV + PRE + TST | 71.0 | 19.3 |
| Mata Atlântica | EV + PRE + WAI | 69.0 | 12.0 |
| Mata Atlântica | EV + WAI + BURNlog | 66.4 | 0.7 |
| Mata Atlântica | EV + PRE + BURNlog | 66.0 | 1.3 |
| Mata Atlântica | EV + TST + BURNlog | 58.4 | 0.7 |
| Mata Atlântica | PRE + TST + WAI | 55.7 | 0.0 |
| Mata Atlântica | PRE + TST + BURNlog | 47.7 | 0.0 |
| Mata Atlântica | PRE + WAI + BURNlog | 47.6 | 0.0 |
| Mata Atlântica | TST + WAI + BURNlog | 41.4 | 0.0 |
| Cerrado | EV + PRE + WAI | 96.7 | 82.7 |
| Cerrado | EV + PRE + TST | 96.3 | 14.0 |
| Cerrado | EV + PRE + BURNlog | 95.8 | 2.7 |
| Cerrado | EV + TST + WAI | 95.2 | 0.7 |
| Cerrado | EV + WAI + BURNlog | 94.2 | 0.0 |
| Cerrado | PRE + TST + WAI | 93.8 | 0.0 |
| Cerrado | EV + TST + BURNlog | 93.3 | 0.0 |
| Cerrado | PRE + WAI + BURNlog | 92.6 | 0.0 |
| Cerrado | TST + WAI + BURNlog | 88.8 | 0.0 |
| Cerrado | PRE + TST + BURNlog | 71.3 | 0.0 |
| Caatinga | EV + PRE + TST | 94.2 | 36.0 |
| Caatinga | EV + PRE + WAI | 94.0 | 22.0 |
| Caatinga | PRE + TST + WAI | 93.9 | 28.7 |
| Caatinga | EV + PRE + BURNlog | 93.7 | 4.7 |
| Caatinga | EV + TST + WAI | 93.0 | 4.0 |
| Caatinga | PRE + WAI + BURNlog | 92.8 | 2.0 |
| Caatinga | EV + WAI + BURNlog | 92.6 | 1.3 |
| Caatinga | EV + TST + BURNlog | 92.5 | 0.7 |
| Caatinga | TST + WAI + BURNlog | 89.7 | 0.7 |
| Caatinga | PRE + TST + BURNlog | 77.9 | 0.0 |

**Diferença pareada 1ª − 2ª (150 partições)**

| bioma | primeira | segunda | dif_media_pp | ic95_inf | ic95_sup | prop_particoes_primeira_maior | p_wilcoxon |
|---|---|---|---|---|---|---|---|
| Mata Atlântica | EV + TST + WAI | EV + PRE + TST | 2.6885 | 2.0986 | 3.2956 | 74.0 | 0.0 |
| Cerrado | EV + PRE + WAI | EV + PRE + TST | 0.3903 | 0.3339 | 0.4484 | 86.0 | 0.0 |
| Caatinga | EV + PRE + TST | EV + PRE + WAI | 0.1486 | 0.0529 | 0.2443 | 62.0 | 0.0023 |

**Importância por permutação: parcela de cada variável na queda total do R² (%)**

| bioma | variavel | GroupKFold | TimeSeriesSplit |
|---|---|---|---|
| Caatinga | EV | 80.4 | 76.9 |
| Caatinga | PRE | 2.4 | 3.8 |
| Caatinga | SAZcos | 8.9 | 4.7 |
| Caatinga | SAZsin | 4.6 | 3.2 |
| Caatinga | TST | 3.7 | 11.4 |
| Cerrado | EV | 85.3 | 77.7 |
| Cerrado | PRE | 4.1 | 3.7 |
| Cerrado | SAZcos | 4.4 | 6.1 |
| Cerrado | SAZsin | 3.6 | 5.0 |
| Cerrado | WAI | 2.6 | 7.5 |
| Mata Atlântica | EV | 41.9 | 38.2 |
| Mata Atlântica | SAZcos | 16.6 | 16.7 |
| Mata Atlântica | SAZsin | 11.0 | 15.2 |
| Mata Atlântica | TST | 9.4 | 10.0 |
| Mata Atlântica | WAI | 21.1 | 19.9 |

**Nulos que preservam a estrutura temporal (R² de validação cruzada, %)**

| bioma | teste | n | r2_original | r2_nulo_media | r2_nulo_dp | r2_nulo_max | p_empirico |
|---|---|---|---|---|---|---|---|
| Mata Atlântica | Deslocamento circular (todos os 296) | 296 | 73.81 | -2.978 | 6.114 | 52.435 | 0.003 |
| Mata Atlântica | Deslocamento circular múltiplo de 12 meses (calendário preservado) | 24 | 73.81 | -3.233 | 5.03 | 7.28 | 0.04 |
| Mata Atlântica | Permutação de anos inteiros (calendário preservado) | 100 | 73.81 | 4.428 | 2.748 | 11.518 | 0.01 |
| Cerrado | Deslocamento circular (todos os 296) | 296 | 96.659 | 36.295 | 12.864 | 81.772 | 0.003 |
| Cerrado | Deslocamento circular múltiplo de 12 meses (calendário preservado) | 24 | 96.659 | 35.479 | 10.977 | 58.714 | 0.04 |
| Cerrado | Permutação de anos inteiros (calendário preservado) | 100 | 96.659 | 62.379 | 1.444 | 67.001 | 0.01 |
| Caatinga | Deslocamento circular (todos os 296) | 296 | 94.162 | 23.915 | 8.879 | 72.523 | 0.003 |
| Caatinga | Deslocamento circular múltiplo de 12 meses (calendário preservado) | 24 | 94.162 | 22.464 | 7.45 | 36.472 | 0.04 |
| Caatinga | Permutação de anos inteiros (calendário preservado) | 100 | 94.162 | 40.205 | 2.323 | 44.739 | 0.01 |
