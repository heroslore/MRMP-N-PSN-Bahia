# -*- coding: utf-8 -*-
"""Sensibilidade do MRMP-N ao percentil de winsorização da PSN no treino (0, 1, 3 e 5%), mesmo pipeline
do Modelo_PSN.py (RepeatedKFold 5 x 30, GridSearchCV do alfa). Saída: resultados_2001_2025/sensibilidade_winsor.csv"""
import os, json, numpy as np, pandas as pd
from concurrent.futures import ProcessPoolExecutor
from sklearn.preprocessing import StandardScaler, PolynomialFeatures
from sklearn.linear_model import Ridge
from sklearn.model_selection import RepeatedKFold, GridSearchCV
from sklearn.metrics import r2_score, mean_squared_error
BASE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(BASE, 'resultados_2001_2025')
df = pd.read_excel(os.path.join(BASE, 'Dados_base_nova_2001_2025.xlsx')); df.columns = df.columns.str.strip()
df['saz_sin'] = np.sin(2 * np.pi * df['MÊS'] / 12); df['saz_cos'] = np.cos(2 * np.pi * df['MÊS'] / 12)
X_VARS = {b: json.load(open(os.path.join(RES, 'numeros_extra.json'), encoding='utf-8'))[b]['x'] for b in ['MA', 'CE', 'CA']}
NOME = {'MA': 'Mata Atlântica', 'CE': 'Cerrado', 'CA': 'Caatinga'}

def rodar(args):
    b, pct = args
    cols = [f'{v}_{b}' for v in X_VARS[b]] + ['saz_sin', 'saz_cos']; x, y = df[cols], df[f'NP_{b}']
    tr_l, te_l, rm_l = [], [], []
    for tr, te in RepeatedKFold(n_splits=5, n_repeats=30, random_state=42).split(x):
        y_tr = y.iloc[tr].clip(lower=np.percentile(y.iloc[tr], pct)) if pct > 0 else y.iloc[tr]; y_te = y.iloc[te]
        sc = StandardScaler(); pf = PolynomialFeatures(degree=2, include_bias=False)
        Xtr = pf.fit_transform(sc.fit_transform(x.iloc[tr])); Xte = pf.transform(sc.transform(x.iloc[te]))
        m = GridSearchCV(Ridge(), {'alpha': [0.1, 1.0, 10.0, 50.0, 100.0]}, cv=5, scoring='r2').fit(Xtr, y_tr.values).best_estimator_
        p = m.predict(Xte); tr_l.append(r2_score(y_tr, m.predict(Xtr)) * 100); te_l.append(r2_score(y_te, p) * 100); rm_l.append(np.sqrt(mean_squared_error(y_te, p)))
    return dict(bioma=NOME[b], percentil=pct, r2_treino=round(np.mean(tr_l), 1), r2_teste=round(np.mean(te_l), 1), gap_pp=round(np.mean(tr_l) - np.mean(te_l), 1), rmse=round(np.mean(rm_l), 2))

if __name__ == '__main__':
    with ProcessPoolExecutor(max_workers=4) as ex:
        res = list(ex.map(rodar, [(b, p) for b in ['MA', 'CE', 'CA'] for p in (0, 1, 3, 5)]))
    out = pd.DataFrame(res); out.to_csv(os.path.join(RES, 'sensibilidade_winsor.csv'), index=False); print(out.to_string(index=False))
