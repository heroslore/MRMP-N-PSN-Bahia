# -*- coding: utf-8 -*-
"""
Sensibilidade, script 01: diagnósticos de métrica, validação, colinearidade e resíduos.
Comando: python t01_diagnosticos.py   (saída em resultados/t01_*.csv e no console)

Itens:
 A  RKF 5x30 oficial recalculado: R² médio por partição x R² agregado (OOF por repetição);
    RMSE médio por partição x RMSE agregado; frequência do alfa escolhido (bordas da grade).
 B  Linhas de base sazonais sob RKF, GroupKFold(ano) e TimeSeriesSplit:
    climatologia mensal do treino; Ridge só com harmônicos; modelo oficial.
    Skill sobre a climatologia: SS = 1 - SSE_modelo / SSE_climatologia (agregado no teste).
 C  Variantes temporais mais estritas: TimeSeriesSplit com gap=12; GroupKFold com
    zona de exclusão (h-block) de 3 e 6 meses em torno dos anos de teste; composição dos folds.
 D  Resíduos: ACF(1) e Ljung-Box(12) dos resíduos no ajuste completo (como o oficial)
    x resíduos fora da amostra (KFold OOF, GroupKFold OOF).
 E  Posto da matriz polinomial (identidade sin²+cos²=1), número de condição, VIF.
 F  Winsorização: meses afetados por partição; R² de treino contra y winsorizado x bruto.
"""
import os, json, time
import numpy as np, pandas as pd
from concurrent.futures import ProcessPoolExecutor
from sklearn.model_selection import RepeatedKFold, GroupKFold, TimeSeriesSplit, KFold
from sklearn.preprocessing import StandardScaler, PolynomialFeatures
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score
from statsmodels.stats.outliers_influence import variance_inflation_factor
from statsmodels.tools import add_constant
from statsmodels.tsa.stattools import acf
from statsmodels.stats.diagnostic import acorr_ljungbox
from auditoria_comum import carregar_base, ajustar, prever, X_OFICIAL, ALPHAS, PASTA

OUT = os.path.join(PASTA, 'resultados'); os.makedirs(OUT, exist_ok=True)


def esquemas(d):
    n = len(d); anos = d['ANO'].values
    return {
        'RKF5x30': list(RepeatedKFold(n_splits=5, n_repeats=30, random_state=42).split(np.zeros(n))),
        'GroupKFold_ano': list(GroupKFold(n_splits=5).split(np.zeros(n), groups=anos)),
        'TimeSeriesSplit': list(TimeSeriesSplit(n_splits=5).split(np.zeros(n))),
        'TSS_gap12': list(TimeSeriesSplit(n_splits=5, gap=12).split(np.zeros(n))),
        'GKF_hblock3': hblock(GroupKFold(n_splits=5).split(np.zeros(n), groups=anos), 3),
        'GKF_hblock6': hblock(GroupKFold(n_splits=5).split(np.zeros(n), groups=anos), 6),
    }


def hblock(splits, h):
    out = []
    for tr, te in splits:
        te_set = np.asarray(te)
        dist = np.min(np.abs(np.asarray(tr)[:, None] - te_set[None, :]), axis=1)
        out.append((np.asarray(tr)[dist > h], te))
    return out


def clim_pred(d, tr, te, y):
    mes = d['MÊS'].values
    c = pd.Series(y[tr]).groupby(mes[tr]).mean()
    return c.reindex(mes[te]).values


def avaliar_bioma(b):
    t0 = time.time()
    d = carregar_base(); y = d[f'PSN_{b}'].values; x = d[X_OFICIAL[b]]; xh = d[['saz_sin', 'saz_cos']]
    n = len(d); res = {'bioma': b}
    E = esquemas(d)
    linhas = []
    oof = {}
    for nome, splits in E.items():
        r2m, r2h, r2c, ssm, rm, ma, al, r2tr_w, r2tr_b, nwin = [], [], [], [], [], [], [], [], [], []
        sse_m = sse_c = sse_h = sst = 0.0
        pred_oof = np.full(n, np.nan)
        for tr, te in splits:
            sc, pf, m, a, r2tr = ajustar(x.iloc[tr], y[tr], 2)
            p = prever(sc, pf, m, x.iloc[te]); pred_oof[te] = p
            sch, pfh, mh, _, _ = ajustar(xh.iloc[tr], y[tr], 2)
            ph = prever(sch, pfh, mh, xh.iloc[te])
            pc = clim_pred(d, tr, te, y)
            yt = y[te]
            r2m.append(r2_score(yt, p) * 100); r2h.append(r2_score(yt, ph) * 100); r2c.append(r2_score(yt, pc) * 100)
            ssm.append((1 - np.sum((yt - p) ** 2) / np.sum((yt - pc) ** 2)) * 100)
            rm.append(np.sqrt(np.mean((yt - p) ** 2))); ma.append(np.mean(np.abs(yt - p))); al.append(a)
            sse_m += np.sum((yt - p) ** 2); sse_c += np.sum((yt - pc) ** 2); sse_h += np.sum((yt - ph) ** 2)
            sst += np.sum((yt - yt.mean()) ** 2)
            lim = np.percentile(y[tr], 3); nwin.append(int((y[tr] < lim).sum()))
            r2tr_w.append(r2tr); r2tr_b.append(r2_score(y[tr], prever(sc, pf, m, x.iloc[tr])) * 100)
        oof[nome] = pred_oof
        al = np.array(al)
        linhas.append(dict(bioma=b, esquema=nome, n_particoes=len(splits),
                           r2_modelo_media_particoes=np.mean(r2m), r2_harmonicos_media=np.mean(r2h),
                           r2_climatologia_media=np.mean(r2c), skill_vs_clim_media_particoes=np.mean(ssm),
                           skill_vs_clim_agregado=(1 - sse_m / sse_c) * 100,
                           skill_harm_vs_clim_agregado=(1 - sse_h / sse_c) * 100,
                           rmse_media_particoes=np.mean(rm), mae_media_particoes=np.mean(ma),
                           pct_alpha_0_1=100 * np.mean(al == 0.1), pct_alpha_100=100 * np.mean(al == 100.0),
                           alpha_media=al.mean(), alpha_mediana=np.median(al),
                           alpha_freq=json.dumps({str(k): int((al == k).sum()) for k in ALPHAS}),
                           n_winsorizados_medio=np.mean(nwin), r2_treino_vs_y_winsor=np.mean(r2tr_w),
                           r2_treino_vs_y_bruto=np.mean(r2tr_b),
                           treino_medio=np.mean([len(a) for a, _ in splits]), teste_medio=np.mean([len(c) for _, c in splits])))
    # A: agregado por repetição no RKF (cada repetição cobre os 297 meses uma vez)
    rk = E['RKF5x30']; r2_rep, rmse_rep = [], []
    for r in range(30):
        pr = np.full(n, np.nan)
        for tr, te in rk[5 * r:5 * r + 5]:
            sc, pf, m, _, _ = ajustar(x.iloc[tr], y[tr], 2); pr[te] = prever(sc, pf, m, x.iloc[te])
        r2_rep.append(r2_score(y, pr) * 100); rmse_rep.append(np.sqrt(np.mean((y - pr) ** 2)))
    res['RKF_r2_agregado_media_reps'] = float(np.mean(r2_rep)); res['RKF_rmse_agregado_media_reps'] = float(np.mean(rmse_rep))
    # composição dos folds temporais
    anos = d['ANO'].values
    res['GKF_anos_teste'] = [sorted(set(anos[te].tolist())) for _, te in E['GroupKFold_ano']]
    res['TSS_periodos_teste'] = [f"{d.ANO[te[0]]}-{d['MÊS'][te[0]]:02d} a {d.ANO[te[-1]]}-{d['MÊS'][te[-1]]:02d} (treino {len(tr)})" for tr, te in E['TimeSeriesSplit']]
    # D: resíduos in-sample (como o oficial: alpha = média aritmética dos alfas RKF, y bruto) x fora da amostra
    alpha_medio = float(pd.DataFrame(linhas).query("esquema=='RKF5x30'").alpha_media.iloc[0])
    sc = StandardScaler(); pf = PolynomialFeatures(2, include_bias=False); Xp = pf.fit_transform(sc.fit_transform(x))
    res_in = y - Ridge(alpha=alpha_medio).fit(Xp, y).predict(Xp)
    def diag(r):
        r = r[~np.isnan(r)]; return dict(acf1=float(acf(r, nlags=1, fft=False)[1]), lb12_p=float(acorr_ljungbox(r, lags=[12], return_df=True)['lb_pvalue'].iloc[0]), n=len(r))
    kf = np.full(n, np.nan)
    for tr, te in KFold(5, shuffle=True, random_state=42).split(x):
        sc2, pf2, m2, _, _ = ajustar(x.iloc[tr], y[tr], 2); kf[te] = prever(sc2, pf2, m2, x.iloc[te])
    res['residuos'] = {'in_sample_oficial': diag(res_in), 'oof_KFold5': diag(y - kf), 'oof_GroupKFold': diag(y - oof['GroupKFold_ano'])}
    res['alpha_medio_recalc'] = alpha_medio
    # E: posto e condição da matriz grau 2 (base inteira, só diagnóstico)
    for g in (2, 3):
        Xg = PolynomialFeatures(g, include_bias=True).fit_transform(StandardScaler().fit_transform(x))
        sv = np.linalg.svd(Xg, compute_uv=False)
        res[f'grau{g}_colunas_com_intercepto'] = int(Xg.shape[1]); res[f'grau{g}_posto'] = int(np.linalg.matrix_rank(Xg))
        res[f'grau{g}_menor_valor_singular'] = float(sv.min()); res[f'grau{g}_condicao'] = float(sv.max() / sv.min())
    s, c = d['saz_sin'].values, d['saz_cos'].values
    res['max_desvio_sin2_mais_cos2'] = float(np.max(np.abs(s ** 2 + c ** 2 - 1)))
    Xv = add_constant(d[X_OFICIAL[b]])
    res['VIF'] = {v: float(variance_inflation_factor(Xv.values, i + 1)) for i, v in enumerate(X_OFICIAL[b])}
    res['corr_preditores'] = d[X_OFICIAL[b]].corr().round(3).to_dict()
    res['n'] = n; res['duracao_s'] = round(time.time() - t0, 1)
    return linhas, res


if __name__ == '__main__':
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=3) as ex:
        saidas = list(ex.map(avaliar_bioma, ['MA', 'CE', 'CA']))
    tab = pd.DataFrame([l for ls, _ in saidas for l in ls])
    tab.to_csv(os.path.join(OUT, 't01_esquemas_e_linhas_de_base.csv'), index=False)
    extra = {r['bioma']: r for _, r in saidas}
    json.dump(extra, open(os.path.join(OUT, 't01_diagnosticos.json'), 'w', encoding='utf-8'), ensure_ascii=False, indent=1, default=float)
    pd.set_option('display.width', 250); pd.set_option('display.max_columns', 40)
    print(tab.drop(columns=['alpha_freq']).round(2).to_string(index=False))
    print(tab[['bioma', 'esquema', 'alpha_freq']].to_string(index=False))
    for b, r in extra.items():
        print(f"\n### {b}")
        for k, v in r.items():
            if k not in ('bioma', 'corr_preditores'):
                print(f"  {k}: {v}")
        print('  corr:', json.dumps(r['corr_preditores']))
    print(f"\nduração total: {time.time() - t0:.1f} s")
