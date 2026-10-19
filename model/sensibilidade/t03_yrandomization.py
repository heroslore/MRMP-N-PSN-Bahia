# -*- coding: utf-8 -*-
"""
Sensibilidade, script 03: Y-randomization.
Comando: python t03_yrandomization.py   (saída em resultados/t03_*.csv)

(a) Resolução do p empírico: o oficial usa 100 permutações e p = (n_maior_igual + 1)/(N + 1)
    (bloco de Y-randomization do Modelo_PSN.py), cujo mínimo possível é 1/101 = 0,0099. Aqui: mesmo pipeline
    do _r2_cv oficial (alfa fixo = alfa médio, winsorização no treino), N = 999 permutações
    i.i.d., RKF 5x3 (não 5x30, por custo). Mostra que 0,0099 é o piso, não a estimativa.
(b) Y-randomization COM a seleção de variáveis dentro do procedimento: para cada permutação,
    o R² da melhor das 10 combinações (grau 2, GridSearchCV do alfa, RKF 5x2), comparado
    com o mesmo máximo nos dados reais. Dois nulos: permutação i.i.d. de meses e permutação
    de anos inteiros (calendário preservado, como nos nulos temporais do Modelo_PSN.py). 100 permutações.
"""
import os, time, itertools
import numpy as np, pandas as pd
from concurrent.futures import ProcessPoolExecutor
from sklearn.model_selection import RepeatedKFold
from sklearn.metrics import r2_score
from auditoria_comum import N_JOBS, carregar_base, ajustar, prever, colunas_combo, X_OFICIAL, CANDIDATAS, PASTA

OUT = os.path.join(PASTA, 'resultados'); os.makedirs(OUT, exist_ok=True)
COMBOS = list(itertools.combinations(CANDIDATAS, 3))
ALPHA_MEDIO = {'MA': 0.706, 'CE': 0.1, 'CA': 0.1}   # recalculado no t01 (alpha_medio_recalc)


def r2_fixo(x, y, splits, alpha):
    r = []
    for a, b in splits:
        sc, pf, m, _, _ = ajustar(x.iloc[a], y[a], 2, alpha=alpha)
        r.append(r2_score(y[b], prever(sc, pf, m, x.iloc[b])) * 100)
    return np.mean(r)


def tarefa_a(args):
    b, lote = args
    d = carregar_base(); y = d[f'PSN_{b}'].values; x = d[X_OFICIAL[b]]
    sp = list(RepeatedKFold(n_splits=5, n_repeats=3, random_state=42).split(x))
    out = []
    for k in lote:
        yp = y.copy() if k < 0 else np.random.default_rng(k).permutation(y)
        out.append(dict(bioma=b, k=k, r2=r2_fixo(x, yp, sp, ALPHA_MEDIO[b])))
    return out


def perm_anos(y, anos, seed):
    completos = [a for a in np.unique(anos) if (anos == a).sum() == 12]
    perm = np.random.default_rng(seed).permutation(completos); yn = y.copy()
    for dest, orig in zip(completos, perm):
        yn[anos == dest] = y[anos == orig]
    return yn


def tarefa_b(args):
    b, tipo, k = args
    d = carregar_base(); y = d[f'PSN_{b}'].values; anos = d['ANO'].values
    if k < 0: yp = y.copy()
    elif tipo == 'iid': yp = np.random.default_rng(5000 + k).permutation(y)
    else: yp = perm_anos(y, anos, 9000 + k)
    sp = list(RepeatedKFold(n_splits=5, n_repeats=2, random_state=42).split(np.zeros(len(y))))
    melhor, nome = -np.inf, None; sc_of = None
    for c in COMBOS:
        cols = colunas_combo(b, c); r = []
        for a, t in sp:
            sc, pf, m, _, _ = ajustar(d[cols].iloc[a], yp[a], 2)
            r.append(r2_score(yp[t], prever(sc, pf, m, d[cols].iloc[t])) * 100)
        r = np.mean(r)
        if set(cols) == set(X_OFICIAL[b]): sc_of = r
        if r > melhor: melhor, nome = r, '+'.join(c)
    return dict(bioma=b, tipo=tipo, k=k, r2_max=melhor, combo_max=nome, r2_combo_oficial=sc_of)


if __name__ == '__main__':
    t0 = time.time(); N_A = 999
    tarefas_a = [(b, list(range(i, min(i + 50, N_A)))) for b in ['MA', 'CE', 'CA'] for i in range(0, N_A, 50)]
    tarefas_a += [(b, [-1]) for b in ['MA', 'CE', 'CA']]
    tarefas_b = [(b, t, k) for b in ['MA', 'CE', 'CA'] for t in ('iid', 'anos') for k in range(-1, 100)]
    with ProcessPoolExecutor(max_workers=N_JOBS) as ex:
        A = pd.DataFrame([r for lst in ex.map(tarefa_a, tarefas_a) for r in lst])
        B = pd.DataFrame(list(ex.map(tarefa_b, tarefas_b, chunksize=4)))
    A.to_csv(os.path.join(OUT, 't03a_perm999_bruto.csv'), index=False)
    B.to_csv(os.path.join(OUT, 't03b_perm_com_selecao_bruto.csv'), index=False)
    la = []
    for b, g in A.groupby('bioma'):
        o = g[g.k < 0].r2.iloc[0]; nul = g[g.k >= 0].r2.values
        la.append(dict(bioma=b, N=len(nul), r2_original=o, r2_nulo_media=nul.mean(), r2_nulo_dp=nul.std(), r2_nulo_max=nul.max(),
                       p_empirico=((nul >= o).sum() + 1) / (len(nul) + 1), z=(o - nul.mean()) / nul.std(),
                       p_piso_com_100=1 / 101))
    LA = pd.DataFrame(la); LA.to_csv(os.path.join(OUT, 't03a_resumo.csv'), index=False)
    lb = []
    for (b, t), g in B.groupby(['bioma', 'tipo']):
        o = g[g.k < 0].iloc[0]; nul = g[g.k >= 0]
        lb.append(dict(bioma=b, nulo=t, N=len(nul), r2_max_original=o.r2_max, combo_original=o.combo_max,
                       nulo_max_media=nul.r2_max.mean(), nulo_max_dp=nul.r2_max.std(), nulo_max_max=nul.r2_max.max(),
                       nulo_combo_oficial_media=nul.r2_combo_oficial.mean(),
                       p_com_selecao=((nul.r2_max >= o.r2_max).sum() + 1) / (len(nul) + 1),
                       delta_com_selecao_pp=o.r2_max - nul.r2_max.mean(),
                       delta_sem_selecao_pp=o.r2_combo_oficial - nul.r2_combo_oficial.mean(),
                       inflacao_nulo_por_selecao_pp=nul.r2_max.mean() - nul.r2_combo_oficial.mean()))
    LB = pd.DataFrame(lb); LB.to_csv(os.path.join(OUT, 't03b_resumo.csv'), index=False)
    pd.set_option('display.width', 250)
    print(LA.round(4).to_string(index=False)); print(); print(LB.round(3).to_string(index=False))
    print(f'duração: {time.time() - t0:.1f} s')
