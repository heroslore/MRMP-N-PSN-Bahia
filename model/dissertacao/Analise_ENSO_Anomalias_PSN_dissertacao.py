# -*- coding: utf-8 -*-
"""
Análise do ENSO em anomalias mensais (seção 6.5 da dissertação), com a fase ENSO OFICIAL
da NOAA já gravada na coluna `Enso` de Dados_base_nova_2001_2025.xlsx (módulo enso_noaa.py,
que usa a tabela completa do ONI para as bordas da série).

Reproduz, com o modelo da dissertação (grau 2, Ridge, alpha médio), as análises pedidas
pela orientação e feitas originalmente numa rotina exploratória anterior:
  1. Sazonalidade: importância das variáveis com e sem SAZsin/SAZcos e correlação das
     anomalias (Tabela A4 ampliada).
  2. Efeito das fases sobre as anomalias de PSN e dos preditores: posição (Kruskal-Wallis,
     Mann-Whitney), dispersão (Fligner-Killeen) e extremos (qui-quadrado, % < P10) (Tabela 6).
  2b. Unidade amostral: a mesma diferença (anomalia média da fase ativa menos a dos meses
     neutros) com intervalo de confiança de 95% por bootstrap de 10.000 reamostragens em dois
     esquemas: mês como unidade (i.i.d.) e episódio como unidade (blocos, sorteando sequências
     contíguas inteiras de meses na mesma fase), com correção de Benjamini-Hochberg sobre os
     valores-p do esquema por episódio. Tabela A11.
  2c. Estratificação sazonal: a mesma diferença restrita ao trimestre climatologicamente mais
     chuvoso e ao mais seco de cada bioma, com o episódio como unidade e contraste direto entre
     os dois trimestres pareado por episódio. Tabela A12.
  2d. Escala de integração hídrica: correlação de Spearman entre a anomalia de PSN e a anomalia
     da chuva acumulada em janelas de 1 a 12 períodos. Tabela A13.
  3. Mediação: quanto da resposta da PSN o modelo reproduz via cada preditor.
  4. Defasagem: Spearman ONI(t) x anomalia PSN(t+L) e compósitos por fase (Figura 14).
  5. Valores exatos dos boxplots brutos por fase (Tabela A5).
Versão usada na dissertação. Roda a partir de model/dissertacao/ e lê a base e as saídas
do modelo em model/ (Dados_base_nova_2001_2025.xlsx, resultados_2001_2025/numeros_extra.json
e resultados_2001_2025/sazonalidade_variaveis.csv).
Saídas (nesta pasta, sem sobrescrever as do artigo do ENSO): resultados/enso_anomalias/*.csv,
        resultados/enso_anomalias.json e figuras/fig14_compositos_enso.png
"""
import os, json, numpy as np, pandas as pd
from scipy import stats
from sklearn.preprocessing import StandardScaler, PolynomialFeatures
from sklearn.linear_model import Ridge
from sklearn.model_selection import GridSearchCV
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

AQUI = os.path.dirname(os.path.abspath(__file__))           # model/dissertacao/
BASE = os.path.dirname(AQUI); RES = os.path.join(BASE, 'resultados_2001_2025')   # entradas em model/
SAIDA = os.path.join(AQUI, 'resultados')
OUT = os.path.join(SAIDA, 'enso_anomalias'); os.makedirs(OUT, exist_ok=True)
FIG = os.path.join(AQUI, 'figuras'); os.makedirs(FIG, exist_ok=True)
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

# ------------------------------------------------------------------ Tabela A11: unidade amostral (mês x episódio)
# Os meses de uma mesma fase não são observações independentes: pertencem a um número pequeno de
# episódios do ENSO, dentro dos quais a série é fortemente autocorrelacionada. A estatística de
# interesse (Δ = anomalia percentual média da fase ativa − a dos meses neutros) é a mesma nos dois
# esquemas; o que muda é a reamostragem. O bootstrap em blocos (Künsch, 1989; Politis e Romano,
# 1994) sorteia episódios inteiros, preservando a dependência interna de cada um.
EPISODIO = np.zeros(len(d), int)
for i in range(len(d)):
    EPISODIO[i] = 1 if i == 0 else EPISODIO[i - 1] + int(fase[i] != fase[i - 1])
_eps = {f: [np.where((EPISODIO == e))[0] for e in np.unique(EPISODIO[fase == f])] for f in FASES}
J['episodios'] = {f: dict(n_meses=n_fase[f], n_sequencias=len(_eps[f]),
                          duracao_mediana=float(np.median([len(i) for i in _eps[f]])),
                          duracao_min=int(min(len(i) for i in _eps[f])), duracao_max=int(max(len(i) for i in _eps[f]))) for f in FASES}
J['episodios']['total_sequencias'] = int(EPISODIO.max())
N_BOOT = 10000

def _ic_boot(vals_ativa, vals_neutro, blocos_ativa=None, blocos_neutro=None, semente=42):
    """IC95% e p bilateral da diferença de médias, por bootstrap i.i.d. (blocos=None) ou em blocos."""
    rng = np.random.default_rng(semente); difs = np.empty(N_BOOT)
    for k in range(N_BOOT):
        if blocos_ativa is None:
            a = rng.choice(vals_ativa, size=len(vals_ativa), replace=True)
            n = rng.choice(vals_neutro, size=len(vals_neutro), replace=True)
        else:
            ia = rng.integers(0, len(blocos_ativa), len(blocos_ativa)); a = np.concatenate([blocos_ativa[j] for j in ia])
            inn = rng.integers(0, len(blocos_neutro), len(blocos_neutro)); n = np.concatenate([blocos_neutro[j] for j in inn])
        difs[k] = a.mean() - n.mean()
    lo, hi = np.percentile(difs, [2.5, 97.5])
    p = min(1.0, 2 * min((difs <= 0).mean(), (difs >= 0).mean()))
    return lo, hi, p, difs.std()

a11 = []
for b in ['MA', 'CE', 'CA']:
    for var in ['PSN', 'EV', 'PRE', 'TST', 'WAI']:
        col = f'NP_{b}' if var == 'PSN' else f'{var}_{b}'
        _, a_pct = anom(col); ap = a_pct.values
        blocos = {f: [ap[i] for i in _eps[f]] for f in FASES}
        for f in ATIVAS:
            delta = ap[fase == f].mean() - ap[fase == 'Neutro'].mean()
            lo_m, hi_m, p_m, se_m = _ic_boot(ap[fase == f], ap[fase == 'Neutro'])
            lo_e, hi_e, p_e, se_e = _ic_boot(None, None, blocos[f], blocos['Neutro'])
            a11.append(dict(bioma=NOME[b], variavel=var, fase=f, delta_pp=delta,
                            mes_ic_inf=lo_m, mes_ic_sup=hi_m, mes_p=p_m, mes_largura=hi_m - lo_m,
                            ep_ic_inf=lo_e, ep_ic_sup=hi_e, ep_p=p_e, ep_largura=hi_e - lo_e,
                            razao_largura=(hi_e - lo_e) / (hi_m - lo_m),
                            sig_mes=bool(lo_m * hi_m > 0), sig_episodio=bool(lo_e * hi_e > 0)))
A11 = pd.DataFrame(a11)
A11.round(4).to_csv(os.path.join(OUT, 'tabelaA11_unidade_amostral.csv'), index=False)
_sm, _se = int(A11.sig_mes.sum()), int(A11.sig_episodio.sum())
J['unidade_amostral'] = dict(
    n_testes=len(A11), sig_mes=_sm, sig_episodio=_se, perderam=_sm - _se,
    razao_largura_mediana=float(A11.razao_largura.median()), razao_largura_max=float(A11.razao_largura.max()),
    razao_largura_max_caso=A11.loc[A11.razao_largura.idxmax(), ['bioma', 'variavel', 'fase']].to_dict(),
    resistem=[{k: (float(v) if isinstance(v, (int, float, np.floating)) and k not in ('bioma', 'variavel', 'fase') else v)
               for k, v in r._asdict().items() if k != 'Index'} for r in A11[A11.sig_episodio].itertuples()],
    perdidos=[{k: (float(v) if isinstance(v, (int, float, np.floating)) and k not in ('bioma', 'variavel', 'fase') else v)
               for k, v in r._asdict().items() if k != 'Index'} for r in A11[A11.sig_mes & ~A11.sig_episodio].itertuples()],
    tabela={f"{r.bioma}|{r.variavel}|{r.fase}": {k: (float(v) if not isinstance(v, (str, bool)) else v)
            for k, v in r._asdict().items() if k not in ('Index', 'bioma', 'variavel', 'fase')} for r in A11.itertuples()})
print("\n===== UNIDADE AMOSTRAL (Tabela A11) =====")
print(f"sequências contíguas: {J['episodios']['total_sequencias']} " +
      "; ".join(f"{f}: {J['episodios'][f]['n_meses']} meses em {J['episodios'][f]['n_sequencias']} sequências" for f in FASES))
print(f"significativos: {_sm}/30 com o mês como unidade, {_se}/30 com o episódio como unidade")
print(A11[['bioma', 'variavel', 'fase', 'delta_pp', 'mes_ic_inf', 'mes_ic_sup', 'mes_p', 'ep_ic_inf', 'ep_ic_sup', 'ep_p', 'razao_largura']].round(3).to_string(index=False))

# correção de Benjamini-Hochberg sobre os valores-p do esquema por episódio (família dos 30 testes)
A11['bh_ep'] = bh(A11.ep_p.values)
_ord = np.argsort(A11.ep_p.values); _n = len(A11)
_q = np.minimum.accumulate((A11.ep_p.values[_ord] * _n / np.arange(1, _n + 1))[::-1])[::-1]
A11['q_ep'] = np.empty(_n); A11.loc[A11.index[_ord], 'q_ep'] = _q
A11.round(4).to_csv(os.path.join(OUT, 'tabelaA11_unidade_amostral.csv'), index=False)
J['unidade_amostral']['sig_episodio_bh'] = int(A11.bh_ep.sum())
J['unidade_amostral']['q_ep_min'] = float(A11.q_ep.min())
J['unidade_amostral']['q_resistem'] = {f"{r.bioma}|{r.variavel}|{r.fase}": float(r.q_ep) for r in A11[A11.sig_episodio].itertuples()}
for _r in J['unidade_amostral']['resistem']:
    _r['q_ep'] = float(A11[(A11.bioma == _r['bioma']) & (A11.variavel == _r['variavel']) & (A11.fase == _r['fase'])].q_ep.iloc[0])
print(f"após Benjamini-Hochberg sobre os p por episódio: {int(A11.bh_ep.sum())}/30 significativos "
      f"(menor q = {A11.q_ep.min():.3f})")

# ------------------------------------------------------------------ Tabela A12: estratificação sazonal
# O canal hídrico só pode operar quando há chuva a ser modulada. A mesma diferença entre fases foi
# recalculada dentro do trimestre civil climatologicamente mais chuvoso e do mais seco de cada bioma,
# com o episódio como unidade; o contraste entre os dois trimestres é pareado, reamostrando os mesmos
# episódios simultaneamente nos dois recortes.
def _trimestres(b):
    clim = d.groupby('MÊS')[f'PRE_{b}'].mean()
    soma = {m: sum(clim[((m - 1 + k) % 12) + 1] for k in range(3)) for m in range(1, 13)}
    ini_w, ini_s = max(soma, key=soma.get), min(soma, key=soma.get)
    tri = lambda m: [((m - 1 + k) % 12) + 1 for k in range(3)]
    return tri(ini_w), tri(ini_s), soma[ini_w], soma[ini_s], float(clim.sum())

a12 = []
for b in ['MA', 'CE', 'CA']:
    wet, dry, mm_w, mm_s, mm_ano = _trimestres(b)
    _, ap_s = anom(f'NP_{b}'); ap = ap_s.values
    for f in ATIVAS:
        est = {}
        for rot, meses in (('chuvoso', wet), ('seco', dry)):
            sel = d['MÊS'].isin(meses).values
            bl_a = [ap[i[np.isin(i, np.where(sel)[0])]] for i in _eps[f]]
            bl_n = [ap[i[np.isin(i, np.where(sel)[0])]] for i in _eps['Neutro']]
            bl_a = [x for x in bl_a if len(x)]; bl_n = [x for x in bl_n if len(x)]
            delta = np.concatenate(bl_a).mean() - np.concatenate(bl_n).mean()
            lo, hi, p, _ = _ic_boot(None, None, bl_a, bl_n)
            est[rot] = dict(delta=delta, lo=lo, hi=hi, p=p, n_meses=int(sum(len(x) for x in bl_a)), n_ep=len(bl_a))
        # contraste pareado entre trimestres: mesmos episódios reamostrados nos dois recortes
        selw, seld = d['MÊS'].isin(wet).values, d['MÊS'].isin(dry).values
        A_w = [ap[i[np.isin(i, np.where(selw)[0])]] for i in _eps[f]]; A_d = [ap[i[np.isin(i, np.where(seld)[0])]] for i in _eps[f]]
        N_w = [ap[i[np.isin(i, np.where(selw)[0])]] for i in _eps['Neutro']]; N_d = [ap[i[np.isin(i, np.where(seld)[0])]] for i in _eps['Neutro']]
        # o mesmo sorteio de episódios vale para os dois trimestres (contraste pareado por episódio);
        # cada episódio contribui com os meses que possui em cada recorte, e as reamostragens em que
        # algum dos quatro grupos fica vazio são descartadas
        rng = np.random.default_rng(42); difs = np.full(N_BOOT, np.nan)
        for k in range(N_BOOT):
            ia = rng.integers(0, len(A_w), len(A_w)); inn = rng.integers(0, len(N_w), len(N_w))
            gr = [[A_w[j] for j in ia if len(A_w[j])], [A_d[j] for j in ia if len(A_d[j])],
                  [N_w[j] for j in inn if len(N_w[j])], [N_d[j] for j in inn if len(N_d[j])]]
            if any(len(g) == 0 for g in gr): continue
            aw, ad, nw, nd = (np.concatenate(g) for g in gr)
            difs[k] = (aw.mean() - nw.mean()) - (ad.mean() - nd.mean())
        difs = difs[~np.isnan(difs)]
        p_c = min(1.0, 2 * min((difs <= 0).mean(), (difs >= 0).mean()))
        a12.append(dict(bioma=NOME[b], fase=f, tri_chuvoso='-'.join(str(m) for m in wet), mm_chuvoso=mm_w,
                        tri_seco='-'.join(str(m) for m in dry), mm_seco=mm_s, mm_ano=mm_ano,
                        **{f'{k2}_{r2_}': est[r2_][k2] for r2_ in ('chuvoso', 'seco') for k2 in ('delta', 'lo', 'hi', 'p', 'n_meses', 'n_ep')},
                        n_boot_validos=len(difs), contraste=est['chuvoso']['delta'] - est['seco']['delta'],
                        contraste_lo=np.percentile(difs, 2.5), contraste_sup=np.percentile(difs, 97.5), contraste_p=p_c))
A12 = pd.DataFrame(a12)
A12['bh_chuvoso'] = bh(A12.p_chuvoso.values); A12['bh_seco'] = bh(A12.p_seco.values)
A12.round(4).to_csv(os.path.join(OUT, 'tabelaA12_estratificacao_sazonal.csv'), index=False)
J['sazonal'] = dict(
    trimestres={NOME[b]: dict(zip(('chuvoso', 'seco', 'mm_chuvoso', 'mm_seco', 'mm_ano'), _trimestres(b))) for b in ['MA', 'CE', 'CA']},
    tabela={f"{r.bioma}|{r.fase}": {k: (float(v) if not isinstance(v, (str, bool, np.bool_)) else (bool(v) if isinstance(v, np.bool_) else v))
            for k, v in r._asdict().items() if k not in ('Index', 'bioma', 'fase')} for r in A12.itertuples()},
    sig_chuvoso=[f"{r.bioma}|{r.fase}" for r in A12[A12.bh_chuvoso].itertuples()],
    sig_seco=[f"{r.bioma}|{r.fase}" for r in A12[A12.bh_seco].itertuples()])
print("\n===== ESTRATIFICAÇÃO SAZONAL (Tabela A12) =====")
print(A12[['bioma', 'fase', 'tri_chuvoso', 'delta_chuvoso', 'p_chuvoso', 'bh_chuvoso', 'delta_seco', 'p_seco', 'contraste', 'contraste_p']].round(3).to_string(index=False))

# ------------------------------------------------------------------ Tabela A13: escala de integração hídrica
# Quanto tempo de chuva a produtividade de cada bioma integra: correlação entre a anomalia de PSN e a
# anomalia da chuva acumulada nos L períodos terminados no próprio período, para L de 1 a 12.
a13 = []
for b in ['MA', 'CE', 'CA']:
    _, ap_s = anom(f'NP_{b}'); psn = ap_s.values; pre = d[f'PRE_{b}'].values
    for L in range(1, 13):
        acc = pd.Series(pre).rolling(L).sum()
        clim = acc.groupby(d['MÊS'].values).transform('mean')
        a_acc = 100 * (acc - clim) / clim; ok = ~a_acc.isna()
        rho = stats.spearmanr(a_acc[ok], psn[ok.values])
        a13.append(dict(bioma=NOME[b], janela=L, rho=rho.statistic, p=rho.pvalue, n=int(ok.sum())))
A13 = pd.DataFrame(a13)
A13.round(4).to_csv(os.path.join(OUT, 'tabelaA13_chuva_acumulada.csv'), index=False)
J['chuva_acumulada'] = {}
for b in ['MA', 'CE', 'CA']:
    sub = A13[A13.bioma == NOME[b]]; imax = sub.rho.idxmax()
    J['chuva_acumulada'][b] = dict(rho_L1=float(sub[sub.janela == 1].rho.iloc[0]), janela_max=int(sub.loc[imax, 'janela']),
                                   rho_max=float(sub.loc[imax, 'rho']), rho_L12=float(sub[sub.janela == 12].rho.iloc[0]),
                                   retencao_L12_pct=float(100 * sub[sub.janela == 12].rho.iloc[0] / sub.loc[imax, 'rho']),
                                   curva=[float(x) for x in sub.sort_values('janela').rho.values])
print("\n===== CHUVA ACUMULADA (Tabela A13) =====")
print(A13.pivot_table(index='janela', columns='bioma', values='rho').round(3).to_string())

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
json.dump(J, open(os.path.join(SAIDA, 'enso_anomalias.json'), 'w', encoding='utf-8'), ensure_ascii=False, indent=1, default=float)

pd.set_option('display.width', 250)
print('fases:', n_fase)
print(T6[['bioma', 'variavel', 'delta_EN', 'delta_LN', 'p_mw_EN', 'p_mw_LN', 'p_kw', 'p_fligner', 'p_chi2', 'pct_abaixo_P10_EN', 'pct_abaixo_P10_N', 'pct_abaixo_P10_LN']].round(3).to_string(index=False))
print(MED.round(2).to_string(index=False))
print(A4.round(2).to_string(index=False))
print(EV.round(1).to_string(index=False))
print('eps2', J['eps2_range'], 'TST MA EN °C', round(J['tst_MA_EN_graus'], 2))
