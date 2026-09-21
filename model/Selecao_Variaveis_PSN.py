# -*- coding: utf-8 -*-
"""
Busca exaustiva das combinações de 3 variáveis ambientais (entre EV, PRE, TST,
WAI e BURN_log) + saz_sin + saz_cos, por bioma, com o MESMO pipeline do
Modelo_PSN.py (winsorização 3% só no treino, StandardScaler, PolynomialFeatures
grau 2, Ridge com GridSearchCV, RepeatedKFold 5x30, random_state=42).
Saída: resultados_2001_2025/selecao_variaveis.csv
"""
import os, itertools, numpy as np, pandas as pd
from concurrent.futures import ProcessPoolExecutor
from sklearn.preprocessing import StandardScaler, PolynomialFeatures
from sklearn.linear_model import Ridge
from sklearn.model_selection import RepeatedKFold, GridSearchCV
from sklearn.metrics import r2_score

BASE = os.path.dirname(os.path.abspath(__file__))
df = pd.read_excel(os.path.join(BASE, 'Dados_base_nova_2001_2025.xlsx'))
df.columns = df.columns.str.strip()
df['saz_sin'] = np.sin(2 * np.pi * df['MÊS'] / 12); df['saz_cos'] = np.cos(2 * np.pi * df['MÊS'] / 12)
for b in ['MA', 'CE', 'CA']: df[f'BURN_{b}_log'] = np.log1p(df[f'BURN_{b}'])

def avaliar(args):
    bioma, combo = args
    x = df[list(combo) + ['saz_sin', 'saz_cos']]; y = df[f'NP_{bioma}']
    rkf = RepeatedKFold(n_splits=5, n_repeats=30, random_state=42)
    tr_l, te_l, al = [], [], []
    for tr, te in rkf.split(x):
        x_tr, x_te, y_tr, y_te = x.iloc[tr], x.iloc[te], y.iloc[tr].copy(), y.iloc[te]
        y_tr = y_tr.clip(lower=np.percentile(y_tr.values, 3))
        sc = StandardScaler(); Xtr = sc.fit_transform(x_tr); Xte = sc.transform(x_te)
        pf = PolynomialFeatures(degree=2, include_bias=False); Xtr = pf.fit_transform(Xtr); Xte = pf.transform(Xte)
        gs = GridSearchCV(Ridge(), param_grid={'alpha': [0.1, 1.0, 10.0, 50.0, 100.0]}, cv=5, scoring='r2')
        gs.fit(Xtr, y_tr.values); m = gs.best_estimator_
        tr_l.append(r2_score(y_tr, m.predict(Xtr)) * 100); te_l.append(r2_score(y_te, m.predict(Xte)) * 100)
        al.append(gs.best_params_['alpha'])
    return dict(bioma=bioma, variaveis=' + '.join(v.replace(f'_{bioma}', '') for v in combo),
                r2_treino=np.mean(tr_l), r2_teste=np.mean(te_l), dp_teste=np.std(te_l),
                gap_pp=np.mean(tr_l) - np.mean(te_l), alpha_medio=np.mean(al))

if __name__ == '__main__':
    tarefas = []
    for b in ['MA', 'CE', 'CA']:
        cand = [f'EV_{b}', f'PRE_{b}', f'TST_{b}', f'WAI_{b}', f'BURN_{b}_log']
        tarefas += [(b, c) for c in itertools.combinations(cand, 3)]
    with ProcessPoolExecutor(max_workers=4) as ex:
        res = list(ex.map(avaliar, tarefas))
    out = pd.DataFrame(res).sort_values(['bioma', 'r2_teste'], ascending=[True, False])
    out.to_csv(os.path.join(BASE, 'resultados_2001_2025', 'selecao_variaveis.csv'), index=False)
    print(out.round(2).to_string(index=False))
