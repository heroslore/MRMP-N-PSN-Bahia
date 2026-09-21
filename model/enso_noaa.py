# -*- coding: utf-8 -*-
"""
Classificação OFICIAL das fases ENSO (NOAA/CPC) a partir do ONI.

Definição operacional da NOAA (Climate Prediction Center):
  - El Niño: ONI >= +0,5 °C por pelo menos 5 trimestres móveis consecutivos
    (DJF, JFM, FMA, ...).
  - La Niña: ONI <= -0,5 °C por pelo menos 5 trimestres móveis consecutivos.
  - Caso contrário: Neutro.

Um mês isolado com ONI >= 0,5 NÃO é El Niño se não fizer parte de uma
sequência de 5 ou mais. Por isso o critério de limiar simples (usado antes)
superestima os meses de El Niño/La Niña perto das transições de fase.

Cada trimestre móvel é atribuído ao seu mês central (DJF -> janeiro,
JFM -> fevereiro, ..., NDJ -> dezembro), que é exatamente a convenção da
coluna ONI da base (conferido: os 300 valores 2001-2025 do CSV são iguais
aos de oni.ascii.txt da NOAA).

Fonte: https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt
(cópia local: oni_noaa_cpc.txt). A tabela completa (1950-presente) é usada
para que as sequências que começam antes de 2001 ou continuam depois do
último mês da base sejam contadas corretamente nas bordas.
"""

import os
import numpy as np
import pandas as pd

NOME_ONI_NOAA = 'oni_noaa_cpc.txt'
LIMIAR        = 0.5
MIN_TRIMESTRES = 5

_TRIMESTRES = ['DJF', 'JFM', 'FMA', 'MAM', 'AMJ', 'MJJ',
               'JJA', 'JAS', 'ASO', 'SON', 'OND', 'NDJ']


def _fases_por_sequencia(oni):
    """
    Recebe um vetor de ONI em ordem cronológica, sem lacunas, e devolve um
    vetor de fases ('El Niño', 'La Niña', 'Neutro') aplicando a regra dos
    5 trimestres consecutivos.
    """
    oni = np.asarray(oni, dtype=float)
    fase = np.full(len(oni), 'Neutro', dtype=object)
    for rotulo, cond in (('El Niño', oni >= LIMIAR), ('La Niña', oni <= -LIMIAR)):
        inicio = None
        for i, ok in enumerate(list(cond) + [False]):   # sentinela para fechar a última sequência
            if ok and inicio is None:
                inicio = i
            elif not ok and inicio is not None:
                if i - inicio >= MIN_TRIMESTRES:
                    fase[inicio:i] = rotulo
                inicio = None
    return fase


def carregar_oni_noaa(pasta):
    """Lê a tabela oni.ascii.txt da NOAA e devolve DataFrame (ANO, MÊS, ONI)."""
    caminho = os.path.join(pasta, NOME_ONI_NOAA)
    if not os.path.exists(caminho):
        return None
    tab = pd.read_csv(caminho, sep=r'\s+')
    tab['MÊS'] = tab['SEAS'].map({s: i + 1 for i, s in enumerate(_TRIMESTRES)})
    tab = tab.rename(columns={'YR': 'ANO', 'ANOM': 'ONI'})[['ANO', 'MÊS', 'ONI']]
    return tab.sort_values(['ANO', 'MÊS']).reset_index(drop=True)


def classificar_fase_enso_noaa(df, pasta, col_ano='ANO', col_mes='MÊS', col_oni='ONI',
                               verboso=True):
    """
    Devolve uma Series (alinhada ao índice de df) com a fase ENSO oficial de
    cada mês. Usa a tabela completa da NOAA (oni_noaa_cpc.txt) se ela
    existir na pasta; senão aplica a regra só sobre o ONI da própria base
    (com aviso: as sequências das bordas podem ficar truncadas).
    """
    tab = carregar_oni_noaa(pasta)
    if tab is None:
        if verboso:
            print(f"AVISO: '{NOME_ONI_NOAA}' não encontrado — aplicando a regra dos 5 "
                  f"trimestres apenas sobre o ONI da base (bordas podem ficar truncadas).")
        base = df[[col_ano, col_mes, col_oni]].rename(
            columns={col_ano: 'ANO', col_mes: 'MÊS', col_oni: 'ONI'})
        tab = base.sort_values(['ANO', 'MÊS']).reset_index(drop=True)
    tab = tab.copy()
    tab['Enso'] = _fases_por_sequencia(tab['ONI'].values)

    chave = df[[col_ano, col_mes]].rename(columns={col_ano: 'ANO', col_mes: 'MÊS'})
    fase = chave.merge(tab[['ANO', 'MÊS', 'Enso']], on=['ANO', 'MÊS'], how='left')['Enso']
    fase.index = df.index
    faltando = fase.isna().sum()
    if faltando:
        raise ValueError(f"{faltando} mês(es) da base sem ONI na tabela da NOAA — "
                         f"atualize {NOME_ONI_NOAA}.")

    if verboso:
        ultimo = tab.iloc[-1]
        print(f"Fase ENSO: classificação oficial NOAA (ONI >= |0,5| por >= {MIN_TRIMESTRES} "
              f"trimestres consecutivos), tabela até {_TRIMESTRES[int(ultimo['MÊS']) - 1]} "
              f"{int(ultimo['ANO'])}.")
    return fase


def classificar_fase_enso_limiar(oni):
    """Critério antigo (limiar simples, sem persistência) — mantido só para comparação."""
    oni = np.asarray(oni, dtype=float)
    return np.select([oni >= LIMIAR, oni <= -LIMIAR], ['El Niño', 'La Niña'], default='Neutro')
