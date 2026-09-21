# -*- coding: utf-8 -*-
"""
Atalho para rodar SOMENTE o bloco de robustez do Modelo_PSN.py (Tabelas A7 a A10), sem refazer o
modelo principal. O código está em Modelo_PSN.rodar_robustez (seção "ROBUSTEZ DA SELEÇÃO E
COMPLEMENTOS DE VALIDAÇÃO"); o modelo completo já o executa ao final quando RODAR_ROBUSTEZ=1.

Usa o conjunto de variáveis de config_biomas e o alfa médio de cada bioma gravado em
resultados_2001_2025/numeros_extra.json (o mesmo usado no Y-randomization).
Saída: resultados_2001_2025/robustez/*.csv e robustez.json
"""
import os, json
import Modelo_PSN as M

if __name__ == '__main__':
    RES = os.path.join(M.BASE_DIR, 'resultados_2001_2025')
    numx = json.load(open(os.path.join(RES, 'numeros_extra.json'), encoding='utf-8'))
    alphas = {b: numx[b]['alpha_medio'] for b in ['MA', 'CE', 'CA']}
    M.rodar_robustez(M.dados_total_base, M.config_biomas, alphas, ['MA', 'CE', 'CA'], RES,
                     max_workers=min(4, os.cpu_count() or 1))
