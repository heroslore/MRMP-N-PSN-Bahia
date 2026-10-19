# -*- coding: utf-8 -*-
"""
Sensibilidade, script 06: quanto do R² vem de preditores derivados do MODIS (EV e WAI do MOD16A2GF)
na previsão da PSN (MOD17A2HGF), ponto da circularidade PSN-EV.
Comando: python t06_circularidade_modis.py   (saída em resultados/t06_*.csv)

Pipeline oficial (auditoria_comum.ajustar: winsor 3% no treino, scaler, grau 2, GridSearchCV
do alfa), mesmas partições dos scripts 01 e 02: RKF 5x30 (random_state=42), GroupKFold por ano com a
atribuição gulosa estável (reproduz o log oficial) e TimeSeriesSplit(5).
Conjuntos por bioma (todos com saz_sin e saz_cos):
  oficial        conjunto de config_biomas
  so_harm        só harmônicos
  EV             EV
  sem_EV_WAI     conjunto oficial sem EV e sem WAI (o que sobra: PRE e/ou TST)
  PRE            só precipitação (IMERG, fonte não MODIS)
Parte (b): viés médio fora da amostra (obs - pred) do conjunto oficial com e sem a
winsorização inferior de 3% no treino, sob GroupKFold por ano (gulosa estável) e RKF 5x30
(OOF por repetição), para medir se o corte só da cauda inferior desloca as predições.
Processo único (sem paralelismo).
"""
import os, time
import numpy as np, pandas as pd
from sklearn.model_selection import RepeatedKFold, TimeSeriesSplit
from auditoria_comum import carregar_base, r2_cv, ajustar, prever, X_OFICIAL, PASTA
from sklearn.metrics import r2_score
from t04_groupkfold_sensibilidade import gulosa_estavel

OUT = os.path.join(PASTA, 'resultados'); os.makedirs(OUT, exist_ok=True)
H = ['saz_sin', 'saz_cos']

if __name__ == '__main__':
    t0 = time.time(); d = carregar_base(); n = len(d); Z = np.zeros(n)
    esq = {'RKF5x30': list(RepeatedKFold(n_splits=5, n_repeats=30, random_state=42).split(Z)),
           'GroupKFold_ano': gulosa_estavel(d['ANO'].values),
           'TimeSeriesSplit': list(TimeSeriesSplit(5).split(Z))}
    lin = []
    for b in ['MA', 'CE', 'CA']:
        y = d[f'PSN_{b}']
        conj = {'oficial': X_OFICIAL[b], 'so_harm': H, 'EV': [f'EV_{b}'] + H,
                'sem_EV_WAI': [c for c in X_OFICIAL[b] if not c.startswith(('EV_', 'WAI_'))],
                'PRE': [f'PRE_{b}'] + H}
        for nome, cols in conj.items():
            for e, sp in esq.items():
                te, tr, _ = r2_cv(d[cols], y, sp, grau=2)
                lin.append(dict(bioma=b, conjunto=nome, variaveis='+'.join(cols), esquema=e, r2_teste=te, r2_treino=tr))
    R = pd.DataFrame(lin); R.to_csv(os.path.join(OUT, 't06_circularidade_modis.csv'), index=False)
    P = R.pivot_table(index=['bioma', 'conjunto'], columns='esquema', values='r2_teste', sort=False)
    pd.set_option('display.width', 250)
    print(P.round(2).to_string())
    vb = []
    for b in ['MA', 'CE', 'CA']:
        y = d[f'PSN_{b}'].values; x = d[X_OFICIAL[b]]
        for e in ['RKF5x30', 'GroupKFold_ano']:
            for pw in (3, 0):
                res, r2s = [], []
                for tr, te in esq[e]:
                    sc, pf, m, _, _ = ajustar(x.iloc[tr], y[tr], 2, pct_winsor=pw)
                    p = prever(sc, pf, m, x.iloc[te]); res.append(y[te] - p); r2s.append(r2_score(y[te], p) * 100)
                r = np.concatenate(res)
                vb.append(dict(bioma=b, esquema=e, pct_winsor=pw, vies_medio=r.mean(), vies_pct_media_psn=100 * r.mean() / y.mean(),
                               frac_residuo_negativo=(r < 0).mean() * 100, r2_teste_media_particoes=np.mean(r2s)))
    VB = pd.DataFrame(vb); VB.to_csv(os.path.join(OUT, 't06b_vies_winsorizacao.csv'), index=False)
    print(); print(VB.round(3).to_string(index=False))
    print(f'duração: {time.time() - t0:.1f} s')
