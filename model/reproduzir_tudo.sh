#!/bin/bash
# Reprodução completa da dissertação a partir da base bruta (~1 h de CPU; o modelo é a parte lenta).
# Base: base_final_2001_2025_plan1_excel_ptbr.csv.
# ATENCAO: a rotina de extracao no Google Earth Engine que gera essa base
# (npp_modis_gee.py) ainda NAO integra este repositorio. Sem ela, a serie mensal
# nao e reproduzivel a partir dos produtos originais, apenas as analises abaixo.
set -e; cd "$(dirname "$0")"
pip install -q pandas numpy scipy statsmodels scikit-learn matplotlib openpyxl python-docx pillow lxml
export BIOMA_ATIVO=TODOS RODAR_EM_PARALELO=0 TESTAR_GRAUS=1 RODAR_YRANDOMIZATION=1 RODAR_ROBUSTEZ=1 GRAU_MODELO=2
python3 Modelo_PSN.py            | tee resultados_2001_2025/log_modelo_psn_completo.txt   # modelo, VIF, resíduos, Y-rand, graus 1-5, robustez (A7-A10)
python3 Analise_Biomas_PSN.py    | tee resultados_2001_2025/log_analise_biomas.txt         # ENSO (fases oficiais), dispersão, boxplots
python3 coletar_resultados.py                                                            # resumo_geral, VIF, graus, JSONs, figuras por bioma
python3 Selecao_Variaveis_PSN.py                                                         # 10 combinações por bioma (Tabela A2)
python3 Ablacao_Sazonalidade_PSN.py                                                      # Tabela A3 e ciclo anual
python3 Analise_ENSO_Anomalias_PSN.py                                                    # 6.5: Tabela 6, Figura 14, Tabelas A4/A5
python3 Analise_Interanual_PSN.py                                                        # 6.6: Figura 15, Tabela 7
python3 Sensibilidade_Winsor_PSN.py                                                      # Tabela A6
python3 Figuras_Dissertacao.py                                                           # figuras da dissertação (PT) e Figs. 2 a 5 do artigo 1 em figuras_artigo1/pt/
FIG_IDIOMA=EN python3 Figuras_Dissertacao.py                                             # Figs. 2 a 5 do artigo 1 em inglês, em figuras_artigo1/en/
python3 Bootstrap_Blocos_ENSO.py                                                         # artigo 2: mês vs. episódio como unidade amostral, com BH
python3 Teste_Pareado_Corrigido.py                                                       # artigo 1: teste pareado com correção de Nadeau-Bengio
python3 Figuras_Artigo2.py                                                               # artigo 2: figuras 1 a 3
python3 Figura_Area_Estudo.py                                                            # artigo 2: figura 1 (área de estudo)
(cd dissertacao && python3 Analise_ENSO_Anomalias_PSN_dissertacao.py)                    # dissertação: Tabelas A11 e A12 (versão da dissertação do ENSO)
(cd sensibilidade && for s in t01_diagnosticos t04_groupkfold_sensibilidade t05_alfa_redundancia_deriva t06_circularidade_modis t02_selecao_aninhada t03_yrandomization; do python3 $s.py; done)   # análises complementares (~35 min)
(cd sensibilidade_qc && python3 Recalcular_297.py)                                       # sensibilidade da TST ao QC_Day, 297 meses
echo "Pronto: resultados em resultados_2001_2025/, dissertacao/resultados/, sensibilidade/resultados/ e sensibilidade_qc/; figuras em figuras_dissertacao/, figuras_artigo1/ e figuras_artigo2/"
