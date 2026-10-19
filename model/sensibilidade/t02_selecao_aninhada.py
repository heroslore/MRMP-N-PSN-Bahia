# -*- coding: utf-8 -*-
"""
Sensibilidade, script 02: otimismo da seleção de variáveis e de grau (validação aninhada).
Comando: python t02_selecao_aninhada.py   (saída em resultados/t02_*.csv)

No oficial, o conjunto de variáveis (Selecao_Variaveis_PSN.py, 10 combinações C(5,3) +
harmônicos, grau 2) e o grau (TESTAR_GRAUS, critério R²teste - 0,5 x gap, graus 1 a 5)
foram escolhidos com o MESMO RepeatedKFold 5x30 que depois mede o desempenho.
Aqui a seleção é refeita DENTRO de cada partição externa, só com o treino:
  passo 1: combinação com maior R² médio num RKF interno 5x3 (grau 2);
  passo 2: grau 1..5 da combinação escolhida pelo critério oficial no mesmo RKF interno;
  ajuste no treino externo e avaliação no teste externo.
Compara com o conjunto e o grau oficiais avaliados nas MESMAS partições externas.
Partições externas: RKF 5x6 (random_state=42; 30 partições, subconjunto das 150 oficiais),
GroupKFold por ano (atribuição gulosa estável, igual à do log oficial) e TimeSeriesSplit(5).
O RKF interno usa 5x3 (não 5x30) por custo; isso é uma simplificação declarada.
"""
import os, time, itertools
import numpy as np, pandas as pd
from concurrent.futures import ProcessPoolExecutor
from sklearn.model_selection import RepeatedKFold, TimeSeriesSplit
from sklearn.metrics import r2_score
from auditoria_comum import N_JOBS, carregar_base, ajustar, prever, colunas_combo, X_OFICIAL, CANDIDATAS, PASTA
from t04_groupkfold_sensibilidade import gulosa_estavel

OUT = os.path.join(PASTA, 'resultados'); os.makedirs(OUT, exist_ok=True)
COMBOS = list(itertools.combinations(CANDIDATAS, 3))
COMBO_OF = {'MA': ('EV', 'TST', 'WAI'), 'CE': ('EV', 'PRE', 'WAI'), 'CA': ('EV', 'PRE', 'TST')}


def cv_interna(x, y, grau):
    te, tr = [], []
    for a, b in RepeatedKFold(n_splits=5, n_repeats=3, random_state=42).split(x):
        sc, pf, m, _, r2tr = ajustar(x.iloc[a], y[a], grau)
        te.append(r2_score(y[b], prever(sc, pf, m, x.iloc[b])) * 100); tr.append(r2tr)
    return np.mean(te), np.mean(tr)


def tarefa(args):
    b, esquema, k, tr, te = args
    d = carregar_base(); y = d[f'PSN_{b}'].values
    dtr = d.iloc[tr].reset_index(drop=True); ytr = y[tr]
    # passo 1: combinação (grau 2)
    sc1 = {c: cv_interna(dtr[colunas_combo(b, c)], ytr, 2)[0] for c in COMBOS}
    c_best = max(sc1, key=sc1.get)
    # passo 2: grau
    sc2 = {}
    for g in range(1, 6):
        r2te, r2tr = cv_interna(dtr[colunas_combo(b, c_best)], ytr, g)
        sc2[g] = r2te - 0.5 * (r2tr - r2te)
    g_best = max(sc2, key=sc2.get)
    def aval(c, g):
        cols = colunas_combo(b, c); sc, pf, m, a, _ = ajustar(d[cols].iloc[tr], ytr, g)
        p = prever(sc, pf, m, d[cols].iloc[te]); return r2_score(y[te], p) * 100, p
    r2_nest, p_nest = aval(c_best, g_best)
    r2_of, p_of = aval(COMBO_OF[b], 2)
    r2_c2, _ = aval(c_best, 2)
    return dict(bioma=b, esquema=esquema, particao=k, combo_escolhido='+'.join(c_best), grau_escolhido=g_best,
                r2_aninhado=r2_nest, r2_oficial_fixo=r2_of, r2_combo_aninhado_grau2=r2_c2,
                combo_oficial_venceu=(c_best == COMBO_OF[b]), grau2_venceu=(g_best == 2),
                te=list(map(int, te)), p_nest=list(map(float, p_nest)), p_of=list(map(float, p_of)))


if __name__ == '__main__':
    t0 = time.time(); d = carregar_base(); n = len(d); Z = np.zeros(n)
    esq = {'RKF5x6': list(RepeatedKFold(n_splits=5, n_repeats=6, random_state=42).split(Z)),
           'GroupKFold_ano': gulosa_estavel(d['ANO'].values),
           'TimeSeriesSplit': list(TimeSeriesSplit(5).split(Z))}
    tarefas = [(b, e, k, tr, te) for b in ['MA', 'CE', 'CA'] for e, sp in esq.items() for k, (tr, te) in enumerate(sp)]
    with ProcessPoolExecutor(max_workers=N_JOBS) as ex:
        res = list(ex.map(tarefa, tarefas))
    R = pd.DataFrame(res)
    R.drop(columns=['te', 'p_nest', 'p_of']).to_csv(os.path.join(OUT, 't02_aninhada_particoes.csv'), index=False)
    lin = []
    for (b, e), g in R.groupby(['bioma', 'esquema'], sort=False):
        lin.append(dict(bioma=b, esquema=e, n_particoes=len(g), r2_oficial_fixo=g.r2_oficial_fixo.mean(),
                        r2_aninhado=g.r2_aninhado.mean(), otimismo_pp=g.r2_oficial_fixo.mean() - g.r2_aninhado.mean(),
                        r2_combo_aninhado_grau2=g.r2_combo_aninhado_grau2.mean(),
                        pct_combo_oficial=100 * g.combo_oficial_venceu.mean(), pct_grau2=100 * g.grau2_venceu.mean(),
                        combos=g.combo_escolhido.value_counts().to_dict(), graus=g.grau_escolhido.value_counts().to_dict()))
    S = pd.DataFrame(lin); S.to_csv(os.path.join(OUT, 't02_aninhada_resumo.csv'), index=False)
    pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 120)
    print(S.round(2).to_string(index=False))
    print(f'duração: {time.time() - t0:.1f} s')
