# -*- coding: utf-8 -*-
"""
Análise do ENSO em anomalias mensais (seção 6.5 da dissertação), com a fase ENSO OFICIAL
da NOAA já gravada na coluna `Enso` de Dados_base_nova_2001_2025.xlsx (módulo enso_noaa.py,
que usa a tabela completa do ONI para as bordas da série).

Reproduz, com o modelo da dissertação (grau 2, Ridge, alpha médio), as análises pedidas
pela orientação e feitas originalmente em npp_modis/modelo/analise_enso_sazonalidade.py:
  1. Sazonalidade: importância das variáveis com e sem SAZsin/SAZcos e correlação das
     anomalias (Tabela A4 ampliada).
  2. Efeito das fases sobre as anomalias de PSN e dos preditores: posição (Kruskal-Wallis,
     Mann-Whitney), dispersão (Fligner-Killeen) e extremos (qui-quadrado, % < P10) (Tabela 6).
  3. Mediação: quanto da resposta da PSN o modelo reproduz via cada preditor.
  4. Defasagem: Spearman ONI(t) x anomalia PSN(t+L) e compósitos por fase (Figura 14).
  5. Valores exatos dos boxplots brutos por fase (Tabela A5).
Saídas: resultados_2001_2025/enso_anomalias/*.csv, enso_anomalias.json e
        figuras_dissertacao/fig14_compositos_enso.png
"""
import os, json, numpy as np, pandas as pd
from scipy import stats
from sklearn.preprocessing import StandardScaler, PolynomialFeatures
from sklearn.linear_model import Ridge
from sklearn.model_selection import GridSearchCV
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(BASE, 'resultados_2001_2025')
OUT = os.path.join(RES, 'enso_anomalias'); os.makedirs(OUT, exist_ok=True)
FIG = os.path.join(BASE, 'figuras_dissertacao'); os.makedirs(FIG, exist_ok=True)
NUMX = json.load(open(os.path.join(RES, 'numeros_extra.json'), encoding='utf-8'))
NOME = {'MA': 'Mata Atlântica', 'CE': 'Cerrado', 'CA': 'Caatinga'}
FASES = ['El Niño', 'Neutro', 'La Niña']; ATIVAS = ['El Niño', 'La Niña']
SAZ = ['saz_sin', 'saz_cos']

d = pd.read_excel(os.path.join(BASE, 'Dados_base_nova_2001_2025.xlsx')); d.columns = d.columns.str.strip()
d['saz_sin'] = np.sin(2 * np.pi * d['MÊS'] / 12); d['saz_cos'] = np.cos(2 * np.pi * d['MÊS'] / 12)
fase = d['Enso'].astype(str).values
assert set(fase) <= set(FASES), set(fase)
n_fase = {f: int((fase == f).sum()) for f in FASES}
J = {'n_fase': n_fase}

# distribuição das fases pelo calendário
cal = {}
for f in FASES:
    m = d.loc[fase == f, 'MÊS'].value_counts(normalize=True).sort_index() * 100
    cal[f] = {int(k): round(float(x), 1) for k, x in m.items()}
J['fase_por_mes_pct'] = cal
J['fase_nov_fev'] = {f: [min(cal[f].get(m, 0) for m in (11, 12, 1, 2)), max(cal[f].get(m, 0) for m in (11, 12, 1, 2))] for f in ATIVAS}
J['fase_mai_jul'] = {f: [min(cal[f].get(m, 0) for m in (5, 6, 7)), max(cal[f].get(m, 0) for m in (5, 6, 7))] for f in ATIVAS}

def anom(col):
    clim = d.groupby('MÊS')[col].transform('mean')
    return d[col] - clim, 100 * (d[col] - clim) / clim

def box(v):
    v = np.asarray(v, float); q1, med, q3 = np.percentile(v, [25, 50, 75]); iqr = q3 - q1
    return dict(n=len(v), minimo=v.min(), Q1=q1, mediana=med, media=v.mean(), Q3=q3, maximo=v.max(),
                outliers=int(((v < q1 - 1.5 * iqr) | (v > q3 + 1.5 * iqr)).sum()))

# ------------------------------------------------------------------ Tabela 6 e A5
t6, a5, bruto = [], [], []
for b in ['MA', 'CE', 'CA']:
    for var in ['PSN', 'EV', 'PRE', 'TST', 'WAI']:
        col = f'NP_{b}' if var == 'PSN' else f'{var}_{b}'
        a_abs, a_pct = anom(col)
        g = {f: a_abs[fase == f].values for f in FASES}; gp = {f: a_pct[fase == f].values for f in FASES}
        kw = stats.kruskal(*g.values()); eps2 = max((kw.statistic - 2) / (len(d) - 3), 0)
        p10, p90 = np.percentile(a_abs, [10, 90])
        tab = np.array([[np.sum(g[f] < p10), np.sum(g[f] > p90), np.sum((g[f] >= p10) & (g[f] <= p90))] for f in FASES])
        row = dict(bioma=NOME[b], variavel=var,
                   delta_EN=gp['El Niño'].mean() - gp['Neutro'].mean(), delta_LN=gp['La Niña'].mean() - gp['Neutro'].mean(),
                   p_mw_EN=stats.mannwhitneyu(g['El Niño'], g['Neutro']).pvalue, p_mw_LN=stats.mannwhitneyu(g['La Niña'], g['Neutro']).pvalue,
                   p_kw=kw.pvalue, eps2=eps2, p_fligner=stats.fligner(*g.values()).pvalue, p_chi2=stats.chi2_contingency(tab)[1],
                   **{f'pct_abaixo_P10_{k}': 100 * np.mean(g[f] < p10) for k, f in zip(['EN', 'N', 'LN'], FASES)},
                   **{f'pct_acima_P90_{k}': 100 * np.mean(g[f] > p90) for k, f in zip(['EN', 'N', 'LN'], FASES)},
                   media_abs_EN=g['El Niño'].mean() - g['Neutro'].mean(), media_abs_LN=g['La Niña'].mean() - g['Neutro'].mean())
        t6.append(row)
    for var in ['PSN', 'EV', 'PRE', 'TST', 'WAI']:
        col = f'NP_{b}' if var == 'PSN' else f'{var}_{b}'
        for f in FASES:
            a5.append(dict(bioma=NOME[b], variavel=var, fase=f, **box(d.loc[fase == f, col])))
        bruto.append(dict(bioma=NOME[b], variavel=var, **{f'media_{k}': d.loc[fase == f, col].mean() for k, f in zip(['EN', 'N', 'LN'], FASES)}))
T6 = pd.DataFrame(t6); A5 = pd.DataFrame(a5); BR = pd.DataFrame(bruto)
# Correção de Benjamini-Hochberg (FDR 5%) por família de testes: Mann-Whitney (30), Kruskal (15), Fligner (15), qui-quadrado (15)
def bh(pvals, q=0.05):
    p = np.asarray(pvals, float); n = len(p); order = np.argsort(p); ranked = p[order]
    thr = q * (np.arange(1, n + 1) / n); ok = ranked <= thr; k = np.where(ok)[0].max() + 1 if ok.any() else 0
    sig = np.zeros(n, bool); sig[order[:k]] = True; return sig
mw = bh(np.concatenate([T6.p_mw_EN.values, T6.p_mw_LN.values])); T6['bh_mw_EN'] = mw[:len(T6)]; T6['bh_mw_LN'] = mw[len(T6):]
T6['bh_kw'] = bh(T6.p_kw.values); T6['bh_fligner'] = bh(T6.p_fligner.values); T6['bh_chi2'] = bh(T6.p_chi2.values)
T6.round(4).to_csv(os.path.join(OUT, 'tabela6_anomalias_por_fase.csv'), index=False)
A5.round(3).to_csv(os.path.join(OUT, 'tabelaA5_boxplots_por_fase.csv'), index=False)
BR.round(2).to_csv(os.path.join(OUT, 'medias_brutas_por_fase.csv'), index=False)
J['tabela6'] = {f"{r.bioma}|{r.variavel}": {k: (float(v) if isinstance(v, (int, float, np.floating)) else v) for k, v in r._asdict().items() if k != 'Index'} for r in T6.itertuples()}
J['a5'] = {f"{r.bioma}|{r.variavel}|{r.fase}": {k: float(v) for k, v in r._asdict().items() if k not in ('Index', 'bioma', 'variavel', 'fase')} for r in A5.itertuples()}
J['bruto'] = {f"{r.bioma}|{r.variavel}": {k: float(v) for k, v in r._asdict().items() if k.startswith('media')} for r in BR.itertuples()}
J['eps2_range'] = [float(T6.eps2.min()), float(T6.eps2.max())]
# TST da MA em °C (anomalia absoluta média sob El Niño)
J['tst_MA_EN_graus'] = float(T6[(T6.bioma == 'Mata Atlântica') & (T6.variavel == 'TST')].media_abs_EN.iloc[0])

# ------------------------------------------------------------------ modelo (mesmo da Figura 9) e mediação
def ajuste_completo(b, cols):
    X = d[cols].values; y = d[f'NP_{b}'].values
    sc = StandardScaler(); pf = PolynomialFeatures(degree=2, include_bias=False)
    Xf = pf.fit_transform(sc.fit_transform(X)); m = Ridge(alpha=NUMX[b]['alpha_medio']).fit(Xf, y)
    return lambda Xn: m.predict(pf.transform(sc.transform(Xn)))

med = []
for b in ['MA', 'CE', 'CA']:
    xs = [f'{v}_{b}' for v in NUMX[b]['x']]; pred = ajuste_completo(b, xs + SAZ)
    clim = d.groupby('MÊS')[xs].mean().loc[d['MÊS']].values
    base = np.column_stack([clim, d[SAZ].values]); psn_clim = pred(base)
    apsn = anom(f'NP_{b}')[1]
    for f in ATIVAS:
        obs = apsn[fase == f].mean() - apsn[fase == 'Neutro'].mean(); tot = 0.0
        for j, x in enumerate(xs):
            a_abs = anom(x)[0]; dv = a_abs[fase == f].mean() - a_abs[fase == 'Neutro'].mean()
            Xc = base.copy(); Xc[:, j] += dv; e = 100 * (pred(Xc) - psn_clim).mean() / psn_clim.mean()
            med.append(dict(bioma=NOME[b], fase=f, preditor=NUMX[b]['x'][j], anomalia_preditor=dv, efeito_PSN_pct=e))
        Xc = base.copy()
        for j, x in enumerate(xs):
            a_abs = anom(x)[0]; Xc[:, j] += a_abs[fase == f].mean() - a_abs[fase == 'Neutro'].mean()
        e = 100 * (pred(Xc) - psn_clim).mean() / psn_clim.mean()
        med.append(dict(bioma=NOME[b], fase=f, preditor='TODOS', anomalia_preditor=np.nan, efeito_PSN_pct=e, observado_pct=obs))
MED = pd.DataFrame(med); MED.round(3).to_csv(os.path.join(OUT, 'mediacao_via_modelo.csv'), index=False)
J['mediacao'] = {f"{r.bioma}|{r.fase}|{r.preditor}": dict(efeito=float(r.efeito_PSN_pct), observado=(None if pd.isna(r.observado_pct) else float(r.observado_pct))) for r in MED.itertuples()}

# ------------------------------------------------------------------ defasagem
lag_rows, comp = [], []
LAGS_C = [0, 1, 2, 3, 4, 6, 9, 12]
for b in ['MA', 'CE', 'CA']:
    apsn = anom(f'NP_{b}')[1].values
    for L in range(13):
        r, pv = stats.spearmanr(d.ONI.values[:len(d) - L] if L else d.ONI.values, apsn[L:])
        lag_rows.append(dict(bioma=NOME[b], lag=L, rho=r, p=pv))
    neutro = apsn[fase == 'Neutro']
    for f in ATIVAS:
        idx = np.where(fase == f)[0]
        for L in LAGS_C:
            ii = idx + L; ii = ii[ii < len(d)]; vals = apsn[ii]
            rng = np.random.default_rng(42); boot = [rng.choice(vals, len(vals)).mean() for _ in range(2000)]
            comp.append(dict(bioma=NOME[b], fase=f, lag=L, n=len(vals), media=vals.mean(), ic_inf=np.percentile(boot, 2.5), ic_sup=np.percentile(boot, 97.5),
                             p_vs_neutro=stats.mannwhitneyu(vals, neutro).pvalue))
LAG = pd.DataFrame(lag_rows); CP = pd.DataFrame(comp)
LAG.round(4).to_csv(os.path.join(OUT, 'spearman_ONI_x_PSN_lags.csv'), index=False); CP.round(3).to_csv(os.path.join(OUT, 'compositos_PSN_por_fase.csv'), index=False)
J['lag'] = {NOME[b]: [dict(lag=int(r.lag), rho=float(r.rho), p=float(r.p)) for r in LAG[LAG.bioma == NOME[b]].itertuples()] for b in ['MA', 'CE', 'CA']}
J['comp'] = {f"{r.bioma}|{r.fase}|{r.lag}": dict(media=float(r.media), ic=[float(r.ic_inf), float(r.ic_sup)], p=float(r.p_vs_neutro), n=int(r.n)) for r in CP.itertuples()}

# eventos (>= 5 meses consecutivos na fase)
ev = []; i = 0
while i < len(d):
    if fase[i] != 'Neutro':
        j = i
        while j + 1 < len(d) and fase[j + 1] == fase[i]: j += 1
        if j - i + 1 >= 5:
            row = dict(fase=fase[i], inicio=f"{int(d.ANO[i])}-{int(d['MÊS'][i]):02d}", fim=f"{int(d.ANO[j])}-{int(d['MÊS'][j]):02d}", meses=j - i + 1,
                       oni_pico=float(d.ONI.iloc[i:j + 1].abs().max()))
            for b in ['MA', 'CE', 'CA']:
                a = anom(f'NP_{b}')[1]; row[f'durante_{b}'] = float(a.iloc[i:j + 1].mean()); row[f'depois3m_{b}'] = float(a.iloc[j + 1:j + 4].mean()) if j + 1 < len(d) else np.nan
            ev.append(row)
        i = j + 1
    else: i += 1
EV = pd.DataFrame(ev); EV.round(2).to_csv(os.path.join(OUT, 'eventos_ENSO_e_PSN.csv'), index=False)
J['eventos'] = EV.to_dict('records')

# ------------------------------------------------------------------ Figura 14 (3 linhas, uma por bioma)
COR = {'El Niño': '#D62728', 'La Niña': '#1F77B4'}
plt.rcParams.update({'font.size': 12, 'axes.titlesize': 14, 'axes.labelsize': 12.5, 'xtick.labelsize': 11, 'ytick.labelsize': 11, 'legend.fontsize': 11.5})
fig, axes = plt.subplots(3, 1, figsize=(9.5, 12.5), sharex=True)
for k, b in enumerate(['MA', 'CE', 'CA']):
    ax = axes[k]; ax.axhline(0, color='#555', lw=1)
    for f in ATIVAS:
        q = CP[(CP.bioma == NOME[b]) & (CP.fase == f)]
        ax.errorbar(q.lag + (-0.12 if f == 'El Niño' else 0.12), q.media, yerr=[q.media - q.ic_inf, q.ic_sup - q.media], fmt='-', color=COR[f], lw=1.8, capsize=4, zorder=2)
        sig = q.p_vs_neutro < 0.05
        ax.scatter(q.lag[sig] + (-0.12 if f == 'El Niño' else 0.12), q.media[sig], s=70, color=COR[f], edgecolors='white', zorder=3, label=f'{f} (p < 0,05)')
        ax.scatter(q.lag[~sig] + (-0.12 if f == 'El Niño' else 0.12), q.media[~sig], s=70, facecolors='white', edgecolors=COR[f], lw=1.8, zorder=3, label=f'{f} (n.s.)')
    ax.set_title(f'({"abc"[k]}) {NOME[b]}', fontweight='bold', loc='left'); ax.set_ylabel('Anomalia média da PSN (%)')
    ax.set_xticks(LAGS_C); ax.grid(alpha=0.3); ax.spines[['top', 'right']].set_visible(False)
    if k == 0: ax.legend(ncol=2, frameon=False, loc='upper right')
axes[-1].set_xlabel('Meses após o mês em fase ENSO (defasagem)')
fig.tight_layout(); fig.savefig(os.path.join(FIG, 'fig14_compositos_enso.png'), dpi=300, bbox_inches='tight'); plt.close(fig)

# ------------------------------------------------------------------ sazonalidade: importância com/sem harmônicos e correlações
def importancia(b, cols):
    """Mesmo procedimento da Figura 8 / coeficientes_ridge_*.csv do Modelo_PSN.py: Ridge com o alfa médio
    dos folds, ajustado à série completa, sem winsorização; importância = soma de |coef| por variável."""
    X = d[cols].values; y = d[f'NP_{b}'].values
    sc = StandardScaler(); pf = PolynomialFeatures(degree=2, include_bias=False); Xf = pf.fit_transform(sc.fit_transform(X))
    m = Ridge(alpha=NUMX[b]['alpha_medio']).fit(Xf, y)
    names = pf.get_feature_names_out(cols); acc = {}
    for nm, c in zip(names, m.coef_):
        for var in set(x.replace(f'_{b}', '').replace('^2', '') for x in nm.split(' ')):
            acc[var] = acc.get(var, 0) + abs(c)
    s = sum(acc.values()); return {k: 100 * v / s for k, v in acc.items()}
SAZR = pd.read_csv(os.path.join(RES, 'sazonalidade_variaveis.csv'))
a4 = []
for b in ['MA', 'CE', 'CA']:
    xs = [f'{v}_{b}' for v in NUMX[b]['x']]; i_sem = importancia(b, xs); i_com = importancia(b, xs + SAZ)
    apsn_abs = anom(f'NP_{b}')[0]
    for v in NUMX[b]['x']:
        a4.append(dict(bioma=NOME[b], variavel=v, r2_ciclo=float(SAZR[(SAZR.bioma == NOME[b]) & (SAZR.variavel == v)].r2.iloc[0]),
                       r_bruto=stats.pearsonr(d[f'{v}_{b}'], d[f'NP_{b}'])[0], r_anom=stats.pearsonr(anom(f'{v}_{b}')[0], apsn_abs)[0],
                       imp_sem=i_sem[v], imp_com=i_com[v], imp_fig8=NUMX['importancia'][b].get(v)))
    for s_ in ('saz_sin', 'saz_cos'):
        a4.append(dict(bioma=NOME[b], variavel=s_.replace('saz_', 'SAZ'), r2_ciclo=np.nan, r_bruto=np.nan, r_anom=np.nan, imp_sem=np.nan, imp_com=i_com[s_], imp_fig8=NUMX['importancia'][b].get(s_)))
A4 = pd.DataFrame(a4); A4.round(3).to_csv(os.path.join(OUT, 'tabelaA4_sazonalidade_importancia.csv'), index=False)
J['a4'] = {f"{r.bioma}|{r.variavel}": {k: (None if pd.isna(v) else float(v)) for k, v in r._asdict().items() if k not in ('Index', 'bioma', 'variavel')} for r in A4.itertuples()}
json.dump(J, open(os.path.join(RES, 'enso_anomalias.json'), 'w', encoding='utf-8'), ensure_ascii=False, indent=1, default=float)

pd.set_option('display.width', 250)
print('fases:', n_fase)
print(T6[['bioma', 'variavel', 'delta_EN', 'delta_LN', 'p_mw_EN', 'p_mw_LN', 'p_kw', 'p_fligner', 'p_chi2', 'pct_abaixo_P10_EN', 'pct_abaixo_P10_N', 'pct_abaixo_P10_LN']].round(3).to_string(index=False))
print(MED.round(2).to_string(index=False))
print(A4.round(2).to_string(index=False))
print(EV.round(1).to_string(index=False))
print('eps2', J['eps2_range'], 'TST MA EN °C', round(J['tst_MA_EN_graus'], 2))
