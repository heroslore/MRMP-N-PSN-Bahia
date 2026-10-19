# Análises complementares de sensibilidade

Scripts das análises de sensibilidade citadas na dissertação (seção 5.5, Validação do Modelo, e resultados do capítulo 6). Elas não substituem as estimativas principais do `Modelo_PSN.py`: delimitam a leitura do desempenho. Todas usam o mesmo pipeline de ajuste do modelo, com a winsorização (percentil 3 da PSN), a padronização, a expansão polinomial e a escolha do α (GridSearchCV, grade [0,1; 1; 10; 50; 100], 5 folds) restritas ao treino de cada partição, e os preditores de `config_biomas` (MA: EV, TST, WAI; CE: EV, PRE, WAI; CA: EV, PRE, TST; todos com SAZsin e SAZcos).

## Como rodar

```bash
cd model/sensibilidade
python3 t01_diagnosticos.py               # cerca de 30 s
python3 t04_groupkfold_sensibilidade.py   # cerca de 30 s
python3 t05_alfa_redundancia_deriva.py    # cerca de 30 s
python3 t06_circularidade_modis.py        # cerca de 2 min
python3 t02_selecao_aninhada.py           # cerca de 15 min
python3 t03_yrandomization.py             # cerca de 15 min
```

A base é lida de `model/Dados_base_nova_2001_2025.xlsx`; as saídas vão para `model/sensibilidade/resultados/`. Os scripts paralelos usam até 4 processos (variável `AUDIT_N_JOBS`; o `t04` usa 10 e o `t01`, 3). Sementes fixas (`random_state=42` nas partições e `default_rng` com semente nas permutações). Os tempos acima são de uma máquina com 12 núcleos.

## Scripts

| Script | O que calcula | Saídas | Onde aparece na dissertação |
|---|---|---|---|
| `auditoria_comum.py` | Funções comuns: leitura da base, pipeline de ajuste (`ajustar`, `prever`, `r2_cv`) e conjuntos de preditores. | | |
| `t01_diagnosticos.py` | Sob RepeatedKFold 5 × 30, GroupKFold por ano, TimeSeriesSplit, TimeSeriesSplit com intervalo de 12 meses e GroupKFold com exclusão de 3 e 6 meses vizinhos: R² do modelo, de um Ridge só com os harmônicos e da climatologia mensal do treino; redução do erro quadrático em relação à climatologia; frequência do α nas bordas da grade; ACF(1) e Ljung-Box dos resíduos dentro e fora da amostra; posto e número de condição da matriz polinomial; VIF. | `t01_esquemas_e_linhas_de_base.csv`, `t01_diagnosticos.json` | Desempenho preditivo do modelo por bioma: redução do erro em relação à climatologia (RepeatedKFold: 90,6% no Cerrado, 89,9% na Caatinga, 72,0% na Mata Atlântica); R² só com os harmônicos (63,4%, 43,5% e 9,5%); TimeSeriesSplit com intervalo de 12 meses (Mata Atlântica 65,7%) e GroupKFold com exclusão de 6 meses (Mata Atlântica 65,3%). |
| `t02_selecao_aninhada.py` | Otimismo da seleção: em cada partição externa (RepeatedKFold 5 × 6, GroupKFold por ano com a partição fixa e TimeSeriesSplit), a combinação de variáveis e o grau são escolhidos de novo só com o treino, por um RepeatedKFold interno 5 × 3, e o R² de teste é comparado ao do conjunto e do grau selecionados. | `t02_aninhada_particoes.csv`, `t02_aninhada_resumo.csv` | Desempenho preditivo do modelo por bioma: otimismo de −0,1 a 1,6 ponto percentual, nulo no Cerrado e máximo na Mata Atlântica sob TimeSeriesSplit (67,0% com a seleção fixa contra 65,4% com a aninhada). |
| `t03_yrandomization.py` | (a) Y-randomization com 999 permutações (RepeatedKFold 5 × 3, α fixo), que mostra que 1/101 = 0,0099 é o piso do p empírico com as 100 permutações do modelo; (b) Y-randomization com a seleção das 10 combinações dentro do procedimento, com permutação de meses e de anos inteiros. | `t03a_*.csv`, `t03b_*.csv` | Complementa a leitura do p = 0,0099 da Y-randomization. |
| `t04_groupkfold_sensibilidade.py` | R² do GroupKFold por ano sob outras partições dos anos: GroupKFold do scikit-learn instalado, a mesma regra com ordenação estável (a partição fixa do modelo), 100 partições aleatórias (`shuffle=True`, sementes 0 a 99), um ano de fora (LOYO) e 5 blocos contíguos de 5 anos. | `t04_groupkfold_todos.csv`, `t04_groupkfold_shuffle_resumo.csv` | Ferramentas Computacionais e Reprodutibilidade e Desempenho preditivo do modelo por bioma: na Mata Atlântica, 73,7% na partição adotada contra 69,6% na partição do GroupKFold do scikit-learn 1.8.0 com numpy 2.4.2; em 100 partições aleatórias, 68,1% a 73,5% (percentis 2,5 e 97,5; média de 71,1%), com o valor adotado igualado ou superado em 3% delas. |
| `t05_alfa_redundancia_deriva.py` | (a) grade ampliada do α (`logspace(-4, 3, 15)`); (b) ajuste sem o termo saz_cos², redundante porque sen² + cos² = 1; (c) viés fora da amostra por ano, com cada ano deixado fora do ajuste, e anomalia anual da TST em relação à climatologia mensal de 2001 a 2019. | `t05a_grade_alfa_e_redundancia.csv`, `t05b_importancia_parametrizacao.json`, `t05c_vies_por_ano_LOYO.csv` | Diagnóstico estatístico do modelo: anomalias da TST em 2024 (−1,19 °C na Mata Atlântica, −1,71 °C no Cerrado, −2,62 °C na Caatinga) e viés por ano na Mata Atlântica (positivo de 2001 a 2009, negativo de 2015 a 2024 e +3,67% em 2025). |
| `t06_circularidade_modis.py` | R² de teste (RepeatedKFold 5 × 30, GroupKFold por ano com a partição fixa e TimeSeriesSplit) de modelos reajustados com o conjunto selecionado, sem a EV e o WAI (preditores derivados do MODIS), só com a EV e os harmônicos, só com a PRE e os harmônicos e só com os harmônicos; viés médio fora da amostra com e sem a winsorização. | `t06_circularidade_modis.csv`, `t06b_vies_winsorizacao.csv` | Seção 5.4 e discussão da dependência entre a PSN e a EV do MODIS. |

## Observação sobre a partição dos anos

O `t01` e o esquema `oficial_atual` do `t04` usam o `GroupKFold(n_splits=5)` do scikit-learn instalado, cujo desempate entre os 24 anos com 12 meses depende da versão do numpy e da CPU. Com scikit-learn 1.8.0 e numpy 2.4.2, essa partição difere da partição fixa do `Modelo_PSN.py` (por isso o GroupKFold com exclusão de 6 meses, 65,3% na Mata Atlântica, é descrito na dissertação como "obtido em outra partição dos anos"). Em outro ambiente, essas linhas podem mudar. O `t02`, o `t05` e o `t06` usam a partição com ordenação estável, igual à do modelo, e não dependem do ambiente.

Versões usadas: Python 3.14.3, numpy 2.4.2, scikit-learn 1.8.0, scipy 1.17.1, statsmodels 0.14.6, pandas 3.0.1.
