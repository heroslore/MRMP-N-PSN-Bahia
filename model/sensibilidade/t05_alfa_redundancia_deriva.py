# -*- coding: utf-8 -*-
"""
Sensibilidade, script 05: grade do alfa, termo redundante dos harmônicos e erro por ano.
Comando: python t05_alfa_redundancia_deriva.py   (saída em resultados/t05_*.csv)

(a) Grade ampliada do alfa: logspace(-4, 3, 15) no GridSearchCV (o oficial usa
    [0.1, 1, 10, 50, 100] e o script 01 mostrou alfa = 0,1, a borda inferior, em 100% das
    partições do CE e da CA). R² RKF 5x30, GroupKFold (gulosa estável) e TSS.
(b) Redundância exata: sin² + cos² = 1 torna a matriz de grau 2 de posto incompleto
    (script 01). Ajuste no RKF 5x30 sem a coluna saz_cos² (base equivalente do mesmo espaço)
    e comparação das somas de |coef| por variável (a "importância" da Figura 8) entre as
    duas parametrizações, na série completa com o alfa médio oficial.
(c) Viés fora da amostra por ano: média do resíduo (obs - pred) por ano no
    LeaveOneGroupOut (um ano de fora) e anomalia anual da TST (MOD11A2, satélite Terra)
    em relação à climatologia mensal de 2001-2019 (descritivo; hipótese de deriva).
"""
import os, time, json
import numpy as np, pandas as pd
from concurrent.futures import ProcessPoolExecutor
from sklearn.model_selection import RepeatedKFold, TimeSeriesSplit, LeaveOneGroupOut, GridSearchCV
from sklearn.preprocessing import StandardScaler, PolynomialFeatures
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score
from auditoria_comum import N_JOBS, carregar_base, X_OFICIAL, PASTA
from t04_groupkfold_sensibilidade import gulosa_estavel

OUT = os.path.join(PASTA, 'resultados'); os.makedirs(OUT, exist_ok=True)
GRADE_AMPLA = list(np.logspace(-4, 3, 15))
ALPHA_MEDIO = {'MA': 0.706, 'CE': 0.1, 'CA': 0.1}


def feats(x_tr, x_te, sem_cos2=False):
    sc = StandardScaler(); pf = PolynomialFeatures(2, include_bias=False)
    A = pf.fit_transform(sc.fit_transform(x_tr)); B = pf.transform(sc.transform(x_te))
    nomes = pf.get_feature_names_out(list(x_tr.columns))
    if sem_cos2:
        keep = nomes != 'saz_cos^2'; A, B, nomes = A[:, keep], B[:, keep], nomes[keep]
    return A, B, nomes


def cv(b, esquema, grade, sem_cos2):
    d = carregar_base(); y = d[f'PSN_{b}'].values; x = d[X_OFICIAL[b]]; n = len(y)
    sp = {'RKF5x30': list(RepeatedKFold(n_splits=5, n_repeats=30, random_state=42).split(x)),
          'GroupKFold_ano': gulosa_estavel(d['ANO'].values),
          'TimeSeriesSplit': list(TimeSeriesSplit(5).split(x))}[esquema]
    r2, al = [], []
    for tr, te in sp:
        ytr = np.clip(y[tr], np.percentile(y[tr], 3), None)
        A, B, _ = feats(x.iloc[tr], x.iloc[te], sem_cos2)
        gs = GridSearchCV(Ridge(), {'alpha': grade}, cv=5, scoring='r2', error_score='raise').fit(A, ytr)
        r2.append(r2_score(y[te], gs.best_estimator_.predict(B)) * 100); al.append(gs.best_params_['alpha'])
    al = np.array(al)
    return dict(bioma=b, esquema=esquema, grade='ampla' if len(grade) > 5 else 'oficial', sem_cos2=sem_cos2,
                r2_teste=np.mean(r2), alpha_mediana=float(np.median(al)), alpha_min=float(al.min()), alpha_max=float(al.max()),
                pct_na_borda_inferior=100 * np.mean(al == min(grade)))


def tarefa(args):
    return cv(*args)


def importancia(b, sem_cos2):
    d = carregar_base(); y = d[f'PSN_{b}'].values; x = d[X_OFICIAL[b]]
    A, _, nomes = feats(x, x, sem_cos2)
    m = Ridge(alpha=ALPHA_MEDIO[b]).fit(A, y); acc = {}
    for nm, c in zip(nomes, m.coef_):
        for v in set(t.split('^')[0] for t in nm.split(' ')):
            acc[v] = acc.get(v, 0) + abs(c)
    s = sum(acc.values())
    return {k: round(100 * v / s, 1) for k, v in acc.items()}, float(r2_score(y, m.predict(A)) * 100)


def vies_por_ano(b):
    d = carregar_base(); y = d[f'PSN_{b}'].values; x = d[X_OFICIAL[b]]; anos = d['ANO'].values
    pred = np.full(len(y), np.nan)
    for tr, te in LeaveOneGroupOut().split(x, groups=anos):
        ytr = np.clip(y[tr], np.percentile(y[tr], 3), None)
        A, B, _ = feats(x.iloc[tr], x.iloc[te])
        gs = GridSearchCV(Ridge(), {'alpha': [0.1, 1.0, 10.0, 50.0, 100.0]}, cv=5, scoring='r2', error_score='raise').fit(A, ytr)
        pred[te] = gs.best_estimator_.predict(B)
    r = pd.DataFrame({'ANO': anos, 'res': y - pred, 'obs': y})
    g = r.groupby('ANO').agg(vies_medio=('res', 'mean'), obs_medio=('obs', 'mean'))
    g['vies_pct'] = 100 * g.vies_medio / g.obs_medio
    tst = d[f'TST_{b}']; base = d['ANO'] <= 2019
    clim = d[base].groupby('MÊS')[f'TST_{b}'].mean()
    d['anom_tst'] = tst - d['MÊS'].map(clim)
    g['anom_TST_C'] = d.groupby('ANO')['anom_tst'].mean()
    g['bioma'] = b
    return g.reset_index()


if __name__ == '__main__':
    t0 = time.time()
    tarefas = [(b, e, g, s) for b in ['MA', 'CE', 'CA'] for e in ['RKF5x30', 'GroupKFold_ano', 'TimeSeriesSplit']
               for g, s in [([0.1, 1.0, 10.0, 50.0, 100.0], False), (GRADE_AMPLA, False), ([0.1, 1.0, 10.0, 50.0, 100.0], True)]]
    with ProcessPoolExecutor(max_workers=N_JOBS) as ex:
        A = pd.DataFrame(list(ex.map(tarefa, tarefas)))
        V = pd.concat(list(ex.map(vies_por_ano, ['MA', 'CE', 'CA'])))
    A.to_csv(os.path.join(OUT, 't05a_grade_alfa_e_redundancia.csv'), index=False)
    V.to_csv(os.path.join(OUT, 't05c_vies_por_ano_LOYO.csv'), index=False)
    imp = {b: {'com_cos2': importancia(b, False), 'sem_cos2': importancia(b, True)} for b in ['MA', 'CE', 'CA']}
    json.dump(imp, open(os.path.join(OUT, 't05b_importancia_parametrizacao.json'), 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
    pd.set_option('display.width', 250)
    print(A.round(3).to_string(index=False)); print()
    print(json.dumps(imp, ensure_ascii=False, indent=1)); print()
    print(V[V.ANO >= 2015].round(2).to_string(index=False))
    print(f'duração: {time.time() - t0:.1f} s')
