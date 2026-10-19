# -*- coding: utf-8 -*-
"""
Análises complementares de sensibilidade do MRMP-N: funções comuns aos scripts.

Replica o pipeline do Modelo_PSN.py:
  winsorização inferior no percentil 3 da PSN de treino (só treino)
  -> StandardScaler (fit no treino) -> PolynomialFeatures(grau, include_bias=False)
  -> GridSearchCV(Ridge(), alpha em [0.1, 1, 10, 50, 100], cv=5, scoring='r2').
Os preditores por bioma são os de config_biomas no Modelo_PSN.py.

A base é lida de model/Dados_base_nova_2001_2025.xlsx (pasta acima desta); nada
aqui altera os dados. As saídas vão para model/sensibilidade/resultados/.
"""
import os
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler, PolynomialFeatures
from sklearn.linear_model import Ridge
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import r2_score

PASTA = os.path.dirname(os.path.abspath(__file__))        # model/sensibilidade/
BASE = os.path.dirname(PASTA)                              # model/
ALPHAS = [0.1, 1.0, 10.0, 50.0, 100.0]
SEMENTE = 42
# número de processos dos scripts paralelos (padrão 4; ajuste pela variável AUDIT_N_JOBS)
N_JOBS = int(os.environ.get('AUDIT_N_JOBS', '4'))
X_OFICIAL = {  # config_biomas do Modelo_PSN.py
    'MA': ['EV_MA', 'TST_MA', 'WAI_MA', 'saz_sin', 'saz_cos'],
    'CE': ['EV_CE', 'PRE_CE', 'WAI_CE', 'saz_sin', 'saz_cos'],
    'CA': ['EV_CA', 'PRE_CA', 'TST_CA', 'saz_sin', 'saz_cos'],
}
CANDIDATAS = ['EV', 'PRE', 'TST', 'WAI', 'BURNlog']   # candidatas da seleção de variáveis


def carregar_base():
    d = pd.read_excel(os.path.join(BASE, 'Dados_base_nova_2001_2025.xlsx'))
    d.columns = d.columns.str.strip()
    d = d.rename(columns={'NP_MA': 'PSN_MA', 'NP_CE': 'PSN_CE', 'NP_CA': 'PSN_CA'})
    d['saz_sin'] = np.sin(2 * np.pi * d['MÊS'] / 12)
    d['saz_cos'] = np.cos(2 * np.pi * d['MÊS'] / 12)
    for b in ['MA', 'CE', 'CA']:
        d[f'BURN_{b}_log'] = np.log1p(d[f'BURN_{b}'])
    t = d['ANO'] * 12 + d['MÊS']
    assert t.is_monotonic_increasing and (t.diff().dropna() == 1).all(), 'base fora de ordem ou com lacunas'
    return d


def colunas_combo(bioma, combo):
    return [f'{v}_{bioma}' if v != 'BURNlog' else f'BURN_{bioma}_log' for v in combo] + ['saz_sin', 'saz_cos']


def ajustar(x_tr, y_tr, grau=2, alpha=None, pct_winsor=3):
    """Pipeline oficial. Devolve (sc, pf, modelo, alpha, r2_treino_winsor_%)."""
    y_tr = np.asarray(y_tr, float).ravel()
    if pct_winsor:
        y_tr = np.clip(y_tr, np.percentile(y_tr, pct_winsor), None)
    sc = StandardScaler()
    pf = PolynomialFeatures(degree=grau, include_bias=False)
    Xtr = pf.fit_transform(sc.fit_transform(x_tr))
    if alpha is None:
        # error_score='raise' para que nenhum ajuste falhe em silêncio (nan)
        gs = GridSearchCV(Ridge(), {'alpha': ALPHAS}, cv=5, scoring='r2', error_score='raise').fit(Xtr, y_tr)
        m, a = gs.best_estimator_, gs.best_params_['alpha']
    else:
        m, a = Ridge(alpha=alpha).fit(Xtr, y_tr), alpha
    return sc, pf, m, a, r2_score(y_tr, m.predict(Xtr)) * 100


def prever(sc, pf, m, x):
    return m.predict(pf.transform(sc.transform(x)))


def r2_cv(x, y, splits, grau=2, alpha=None):
    """R² de teste médio (%) por partição e R² de treino médio, pipeline oficial."""
    te, tr = [], []
    for a, b in splits:
        sc, pf, m, _, r2tr = ajustar(x.iloc[a], y.iloc[a], grau, alpha)
        te.append(r2_score(y.iloc[b], prever(sc, pf, m, x.iloc[b])) * 100)
        tr.append(r2tr)
    return float(np.mean(te)), float(np.mean(tr)), te
