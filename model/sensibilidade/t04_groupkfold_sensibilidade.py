# -*- coding: utf-8 -*-
"""
Sensibilidade, script 04: dependência do R² GroupKFold (por ano) em relação à composição dos folds.
Comando: python t04_groupkfold_sensibilidade.py   (saída em resultados/t04_*.csv)

Motivo: o script 01 (sklearn 1.8.0, numpy 2.4.2) obteve R² GroupKFold diferente do log oficial,
com RepeatedKFold e TimeSeriesSplit idênticos. GroupKFold sem shuffle distribui grupos
de tamanho igual (24 anos com 12 meses) na ordem de np.argsort (não estável), de modo que
o empate entre anos é resolvido de forma dependente de implementação.

Esquemas:
  oficial_atual   GroupKFold(5) desta instalação
  argsort_estavel mesma regra gulosa, com argsort estável (kind='stable')
  shuffle_s       GroupKFold(5, shuffle=True, random_state=s), s = 0..99
  LOYO            deixa um ano de fora (25 partições)
  blocos5         5 blocos contíguos de 5 anos (2001-05, 2006-10, ...)
"""
import os, time
import numpy as np, pandas as pd
from concurrent.futures import ProcessPoolExecutor
from sklearn.model_selection import GroupKFold, LeaveOneGroupOut
from sklearn.metrics import r2_score
from auditoria_comum import carregar_base, ajustar, prever, X_OFICIAL, PASTA

OUT = os.path.join(PASTA, 'resultados'); os.makedirs(OUT, exist_ok=True)
OFICIAL_LOG = {'MA': [64.5, 81.2, 74.1, 77.6, 71.3], 'CE': [96.0, 97.4, 96.9, 95.1, 96.6], 'CA': [91.4, 95.6, 94.9, 91.5, 95.5]}


def gulosa_estavel(anos, k=5):
    ug, gi = np.unique(anos, return_inverse=True); cnt = np.bincount(gi)
    idx = np.argsort(cnt, kind='stable')[::-1]; carga = np.zeros(k); g2f = np.zeros(len(ug))
    for j in idx:
        f = np.argmin(carga); carga[f] += cnt[j]; g2f[j] = f
    fold = g2f[gi]
    return [(np.where(fold != f)[0], np.where(fold == f)[0]) for f in range(k)]


def avaliar(args):
    b, nome, splits = args
    d = carregar_base(); y = d[f'PSN_{b}'].values; x = d[X_OFICIAL[b]]; n = len(y)
    r2f = []; pred = np.full(n, np.nan)
    for tr, te in splits:
        sc, pf, m, _, _ = ajustar(x.iloc[tr], y[tr], 2); p = prever(sc, pf, m, x.iloc[te]); pred[te] = p
        r2f.append(r2_score(y[te], p) * 100)
    anos_teste = [sorted(set(d['ANO'].values[te].tolist())) for _, te in splits]
    return dict(bioma=b, esquema=nome, r2_media_folds=np.mean(r2f), r2_agregado=r2_score(y, pred) * 100,
                folds=' '.join(f'{v:.1f}' for v in r2f), anos_teste=str(anos_teste) if len(splits) <= 5 else '')


if __name__ == '__main__':
    t0 = time.time(); d = carregar_base(); anos = d['ANO'].values; n = len(d); Z = np.zeros(n)
    esq = {'oficial_atual': list(GroupKFold(5).split(Z, groups=anos)),
           'argsort_estavel': gulosa_estavel(anos),
           'LOYO': list(LeaveOneGroupOut().split(Z, groups=anos))}
    bl = (anos - 2001) // 5
    esq['blocos5'] = [(np.where(bl != k)[0], np.where(bl == k)[0]) for k in range(5)]
    for s in range(100):
        esq[f'shuffle_{s}'] = list(GroupKFold(5, shuffle=True, random_state=s).split(Z, groups=anos))
    tarefas = [(b, k, v) for b in ['MA', 'CE', 'CA'] for k, v in esq.items()]
    with ProcessPoolExecutor(max_workers=10) as ex:
        res = pd.DataFrame(list(ex.map(avaliar, tarefas, chunksize=10)))
    res.to_csv(os.path.join(OUT, 't04_groupkfold_todos.csv'), index=False)
    pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 200)
    fixos = res[~res.esquema.str.startswith('shuffle_')]
    print(fixos.round(2).to_string(index=False))
    sh = res[res.esquema.str.startswith('shuffle_')].groupby('bioma')['r2_media_folds']
    resumo = sh.describe(percentiles=[0.025, 0.5, 0.975]).round(2)
    resumo['oficial_log'] = [np.mean(OFICIAL_LOG[b]) for b in resumo.index]
    resumo.to_csv(os.path.join(OUT, 't04_groupkfold_shuffle_resumo.csv'))
    print('\nGroupKFold(shuffle=True), 100 sementes: R² médio dos 5 folds\n' + resumo.to_string())
    print('\nfolds do log oficial:', OFICIAL_LOG)
    print(f'duração: {time.time() - t0:.1f} s')
