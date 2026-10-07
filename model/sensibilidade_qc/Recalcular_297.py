# -*- coding: utf-8 -*-
"""
Recalcula a analise de sensibilidade ao QC_Day sobre a base exata dos artigos.

A rodada de extracao cobriu 300 periodos (ate dezembro de 2025) e usou uma
classificacao ENSO propria. A base dos manuscritos termina em setembro de 2025
(297 periodos), porque a NASA encerrou a versao 07 do IMERG nesse mes, e a
precipitacao limita a serie; e sua classificacao ENSO tem 76 meses de El Nino,
149 neutros e 72 de La Nina.

Este script trunca as tres versoes da TST em setembro de 2025 e refaz todos os
indicadores com a climatologia e a classificacao ENSO da propria base, de modo
que os numeros sejam diretamente citaveis nos manuscritos.

Saidas: tst_qc_resumo_297.csv, tst_qc_enso_297.csv
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AQUI = os.path.join(BASE, 'sensibilidade_qc')
SIG = {'Mata Atlântica': 'MA', 'Cerrado': 'CE', 'Caatinga': 'CA'}
VERS = {'tst_sem_filtro': 'sem filtro', 'tst_erro_3k': 'erro ≤ 3 K', 'tst_boa': 'boa qualidade'}
NB, SEED = 10000, 42

S = pd.read_csv(os.path.join(AQUI, 'tst_qc_series_mensais.csv'))
S = S[(S.ano < 2025) | (S.mes <= 9)].copy()
B = pd.read_excel(os.path.join(BASE, 'Dados_base_nova_2001_2025.xlsx'))
B.columns = B.columns.str.strip()

# classificacao ENSO e blocos (episodios) da base dos artigos
B = B.sort_values(['ANO', 'MÊS']).reset_index(drop=True)
fase = B['Enso'].astype(str).values
bloco = np.zeros(len(B), int); k = 0
for i in range(1, len(B)):
    if fase[i] != fase[i - 1]:
        k += 1
    bloco[i] = k
chave = {(a, m): i for i, (a, m) in enumerate(zip(B['ANO'], B['MÊS']))}

def anom_pct(v, meses):
    s = pd.Series(v); clim = s.groupby(meses).transform('mean')
    return (100 * (s - clim) / clim).values

def bh(p, q=0.05):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); ps = p[o]
    sig = np.zeros(n, bool)
    kk = np.where(ps <= (np.arange(1, n + 1) / n) * q)[0]
    if len(kk):
        sig[o[:kk.max() + 1]] = True
    adj = np.minimum.accumulate((ps * n / np.arange(1, n + 1))[::-1])[::-1]
    qv = np.empty(n); qv[o] = np.minimum(adj, 1.0)
    return sig, qv

def trimestres(sg):
    pre = B[f'PRE_{sg}'].groupby(B['MÊS']).mean()
    som = {m: sum(pre[(m + i - 1) % 12 + 1] for i in range(3)) for m in range(1, 13)}
    mm = lambda t: [(t + i - 1) % 12 + 1 for i in range(3)]
    return mm(max(som, key=som.get)), mm(min(som, key=som.get))

# ------------------------------------------------------------------ resumo
linhas = []
for nome, sg in SIG.items():
    d = S[S.bioma == nome].sort_values(['ano', 'mes']).reset_index(drop=True)
    chuv, seco = trimestres(sg)
    base_v = d['tst_sem_filtro'].values
    for col in ['tst_erro_3k', 'tst_boa']:
        alvo = d[col].values
        desc = 100 - d['retido_pct_erro_3k' if col == 'tst_erro_3k' else 'retido_pct_boa'].values
        for rot, sel in [('ano todo', np.ones(len(d), bool)),
                         ('trimestre chuvoso', d['mes'].isin(chuv).values),
                         ('trimestre seco', d['mes'].isin(seco).values)]:
            dif = (alvo - base_v)[sel]
            i = int(np.argmax(np.abs(dif)))
            r = stats.pearsonr(base_v[sel], alvo[sel])[0]
            an_b = anom_pct(base_v, d['mes'].values)[sel]
            rho_an = stats.spearmanr(an_b, desc[sel])[0] if np.std(desc[sel]) > 0 else np.nan
            rho_br = stats.spearmanr(base_v[sel], desc[sel])[0] if np.std(desc[sel]) > 0 else np.nan
            linhas.append(dict(bioma=nome, versao=VERS[col], periodo=rot, n_meses=int(sel.sum()),
                               dif_media_C=round(float(dif.mean()), 4),
                               dif_mediana_C=round(float(np.median(dif)), 4),
                               rmse_C=round(float(np.sqrt((dif ** 2).mean())), 4),
                               r_pearson=round(float(r), 4),
                               maior_dif_abs_C=round(float(np.abs(dif).max()), 4),
                               mes_maior_dif=f"{d['ano'].values[sel][i]}-{d['mes'].values[sel][i]:02d}",
                               pct_descartado_medio=round(float(desc[sel].mean()), 4),
                               pct_descartado_max=round(float(desc[sel].max()), 4),
                               rho_anomalia_vs_descarte=None if np.isnan(rho_an) else round(float(rho_an), 4),
                               rho_bruto_vs_descarte=None if np.isnan(rho_br) else round(float(rho_br), 4)))
R = pd.DataFrame(linhas)
R.to_csv(os.path.join(AQUI, 'tst_qc_resumo_297.csv'), index=False, encoding='utf-8')

# ------------------------------------------------------------------ ENSO
rng = np.random.default_rng(SEED)
ens = []
for nome, sg in SIG.items():
    d = S[S.bioma == nome].sort_values(['ano', 'mes']).reset_index(drop=True)
    idx = np.array([chave[(a, m)] for a, m in zip(d['ano'], d['mes'])])
    f_i = fase[idx]; b_i = bloco[idx]
    for col, vrot in VERS.items():
        a_abs = d[col].values - pd.Series(d[col].values).groupby(d['mes'].values).transform('mean').values
        a_pct = anom_pct(d[col].values, d['mes'].values)
        p10, p90 = np.percentile(a_abs, [10, 90])
        neutro = a_abs[f_i == 'Neutro']
        for ph in ['El Niño', 'La Niña']:
            g = a_abs[f_i == ph]
            p_mw = stats.mannwhitneyu(g, neutro).pvalue
            # bootstrap em blocos por episodio, sobre a anomalia percentual
            bf = [a_pct[(b_i == u)] for u in sorted(set(b_i[f_i == ph]))]
            bn = [a_pct[(b_i == u)] for u in sorted(set(b_i[f_i == 'Neutro']))]
            obs = np.concatenate(bf).mean() - np.concatenate(bn).mean()
            dist = np.empty(NB)
            for j in range(NB):
                rf = np.concatenate([bf[i] for i in rng.integers(0, len(bf), len(bf))])
                rn = np.concatenate([bn[i] for i in rng.integers(0, len(bn), len(bn))])
                dist[j] = rf.mean() - rn.mean()
            lo, hi = np.percentile(dist, [2.5, 97.5])
            p_blk = min(max(2 * min((dist <= 0).mean(), (dist >= 0).mean()), 1 / NB), 1.0)
            ens.append(dict(bioma=nome, versao=vrot, fase=ph, n=int((f_i == ph).sum()),
                            n_episodios=len(bf),
                            anomalia_media_C=round(float(g.mean() - neutro.mean()), 4),
                            anomalia_media_pp=round(float(obs), 4),
                            p_mannwhitney=float('%.4g' % p_mw),
                            ic_blk_inf=round(float(lo), 3), ic_blk_sup=round(float(hi), 3),
                            p_bootstrap_blocos=round(float(p_blk), 4),
                            pct_acima_P90=round(100 * float((g > p90).mean()), 2),
                            pct_neutro_acima_P90=round(100 * float((neutro > p90).mean()), 2)))
E = pd.DataFrame(ens)
for v in VERS.values():
    m = E.versao == v
    sig, qv = bh(E.loc[m, 'p_mannwhitney'].values)
    E.loc[m, 'q_mw'] = np.round(qv, 4); E.loc[m, 'bh_mw'] = sig
E.to_csv(os.path.join(AQUI, 'tst_qc_enso_297.csv'), index=False, encoding='utf-8')

pd.set_option('display.width', 260); pd.set_option('display.max_columns', 40)
print('=== RESUMO (297 períodos) ===')
print(R[R.versao == 'boa qualidade'][['bioma','periodo','dif_media_C','rmse_C','r_pearson',
      'maior_dif_abs_C','mes_maior_dif','pct_descartado_medio','pct_descartado_max',
      'rho_anomalia_vs_descarte']].to_string(index=False))
print()
print('erro ≤ 3 K — descarte máximo em qualquer bioma/recorte:',
      R[R.versao == 'erro ≤ 3 K'].pct_descartado_max.max())
print()
print('=== ENSO (297 períodos, classificação da base: 76/149/72) ===')
print(E[['bioma','versao','fase','n','n_episodios','anomalia_media_C','anomalia_media_pp',
         'p_mannwhitney','q_mw','bh_mw','ic_blk_inf','ic_blk_sup','p_bootstrap_blocos',
         'pct_acima_P90','pct_neutro_acima_P90']].to_string(index=False))
