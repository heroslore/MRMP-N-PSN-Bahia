#!/bin/bash
# Reprodução completa da dissertação a partir da base bruta (~1 h de CPU; o modelo é a parte lenta).
# Base: base_final_2001_2025_plan1_excel_ptbr.csv (gerada por ../npp_modis/npp_modis_gee.py no Google Earth Engine).
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
python3 Figuras_Dissertacao.py                                                           # figuras da dissertação
python3 banca/aplicar_revisao_docx.py && bash banca/render_e_paginas.sh                  # Word + PDF (1ª passagem: páginas)
python3 banca/aplicar_revisao_docx.py && bash banca/render_e_paginas.sh                  # 2ª passagem: listas e sumário com páginas
cp /tmp/claude-0/-home-user-Controle-Gas/*/scratchpad/docx/render/rev.pdf banca/Trabalho_revisado_2001_2025.pdf 2>/dev/null || true
echo "Pronto: banca/Trabalho_revisado_2001_2025.docx"
