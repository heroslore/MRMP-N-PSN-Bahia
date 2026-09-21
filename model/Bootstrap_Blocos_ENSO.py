# -*- coding: utf-8 -*-
"""
Unidade amostral em compositos ENSO: mes vs. episodio (Artigo 2).

Os testes usuais (Mann-Whitney, Kruskal-Wallis) e o bootstrap i.i.d. tratam cada
mes como uma observacao independente. Nao sao: os 76 meses de El Nino pertencem
a 8 episodios e os 72 de La Nina a 9 sequencias contiguas, e a serie dentro de
um mesmo episodio e fortemente autocorrelacionada. O efetivo e o numero de
episodios, nao o de meses.

Este script calcula, para a MESMA estatistica -- a diferenca entre a anomalia
percentual media da fase e a dos meses neutros --, dois intervalos de confianca:

  (a) bootstrap i.i.d., sorteando MESES com reposicao;
  (b) bootstrap em blocos, sorteando EPISODIOS inteiros com reposicao,
      o que preserva a dependencia interna de cada episodio.

Saidas: resultados_2001_2025/enso_anomalias/bootstrap_unidade_amostral.csv
        resultados_2001_2025/enso_anomalias/bootstrap_unidade_amostral_meta.json
"""
import os, json
import numpy as np
import pandas as pd

BASE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(BASE, 'resultados_2001_2025')
OUT = os.path.join(RES, 'enso_anomalias')
os.makedirs(OUT, exist_ok=True)

NOME = {'MA': 'Mata Atlântica', 'CE': 'Cerrado', 'CA': 'Caatinga'}
FASES = ['El Niño', 'Neutro', 'La Niña']
ATIVAS = ['El Niño', 'La Niña']
B = 10000
SEED = 42

d = pd.read_excel(os.path.join(BASE, 'Dados_base_nova_2001_2025.xlsx'))
d.columns = d.columns.str.strip()
fase = d['Enso'].astype(str).values
assert set(fase) <= set(FASES), set(fase)

# blocos = corridas contiguas de meses na mesma fase
bloco_id = np.zeros(len(d), dtype=int)
k = 0
for i in range(1, len(d)):
    if fase[i] != fase[i - 1]:
        k += 1
    bloco_id[i] = k
n_meses = {f: int((fase == f).sum()) for f in FASES}
n_blocos = {f: len(set(bloco_id[fase == f])) for f in FASES}
dur = {f: [int(len(bloco_id[bloco_id == b])) for b in sorted(set(bloco_id[fase == f]))] for f in FASES}

def anom_pct(col):
    clim = d.groupby('MÊS')[col].transform('mean')
    return (100 * (d[col] - clim) / clim).values

def blocos_da_fase(a, f):
    return [a[bloco_id == b] for b in sorted(set(bloco_id[fase == f]))]

def ic_p(dist, B):
    lo, hi = np.percentile(dist, [2.5, 97.5])
    p = 2 * min((dist <= 0).mean(), (dist >= 0).mean())
    return float(lo), float(hi), float(min(max(p, 1 / B), 1.0))

rng = np.random.default_rng(SEED)
linhas = []
for b in ['MA', 'CE', 'CA']:
    for var in ['PSN', 'EV', 'PRE', 'TST', 'WAI']:
        col = f'NP_{b}' if var == 'PSN' else f'{var}_{b}'
        a = anom_pct(col)
        mN = a[fase == 'Neutro']
        bl_N = blocos_da_fase(a, 'Neutro')
        for f in ATIVAS:
            mF = a[fase == f]
            bl_F = blocos_da_fase(a, f)
            obs = mF.mean() - mN.mean()

            # (a) i.i.d. sobre meses
            d_iid = np.empty(B)
            for j in range(B):
                d_iid[j] = rng.choice(mF, len(mF)).mean() - rng.choice(mN, len(mN)).mean()
            lo_i, hi_i, p_i = ic_p(d_iid, B)

            # (b) blocos = episodios
            d_blk = np.empty(B)
            for j in range(B):
                rF = np.concatenate([bl_F[i] for i in rng.integers(0, len(bl_F), len(bl_F))])
                rN = np.concatenate([bl_N[i] for i in rng.integers(0, len(bl_N), len(bl_N))])
                d_blk[j] = rF.mean() - rN.mean()
            lo_b, hi_b, p_b = ic_p(d_blk, B)

            linhas.append(dict(
                bioma=NOME[b], variavel=var, fase=f,
                n_meses=len(mF), n_blocos=len(bl_F),
                delta_pct=round(float(obs), 3),
                iid_inf=round(lo_i, 3), iid_sup=round(hi_i, 3), p_iid=round(p_i, 4),
                blk_inf=round(lo_b, 3), blk_sup=round(hi_b, 3), p_blk=round(p_b, 4),
                largura_iid=round(hi_i - lo_i, 3), largura_blk=round(hi_b - lo_b, 3),
                razao_largura=round((hi_b - lo_b) / (hi_i - lo_i), 2),
                sig_iid=bool(lo_i > 0 or hi_i < 0), sig_blk=bool(lo_b > 0 or hi_b < 0)))

T = pd.DataFrame(linhas)
T.to_csv(os.path.join(OUT, 'bootstrap_unidade_amostral.csv'), index=False, encoding='utf-8')
meta = dict(n_meses=n_meses, n_blocos=n_blocos, duracao_blocos=dur, B=B, seed=SEED,
            razao_largura_mediana=float(T.razao_largura.median()),
            razao_largura_min=float(T.razao_largura.min()),
            razao_largura_max=float(T.razao_largura.max()),
            n_sig_iid=int(T.sig_iid.sum()), n_sig_blk=int(T.sig_blk.sum()), n_testes=len(T))
json.dump(meta, open(os.path.join(OUT, 'bootstrap_unidade_amostral_meta.json'), 'w', encoding='utf-8'),
          ensure_ascii=False, indent=1)

print('meses por fase :', n_meses)
print('blocos por fase:', n_blocos)
print('duracao dos blocos ativos:', {f: dur[f] for f in ATIVAS})
print()
print(T[T.variavel == 'PSN'][['bioma', 'fase', 'delta_pct', 'iid_inf', 'iid_sup', 'p_iid',
                              'blk_inf', 'blk_sup', 'p_blk', 'razao_largura']].to_string(index=False))
print()
print('significativos: i.i.d. %d/%d | blocos %d/%d' % (meta['n_sig_iid'], len(T), meta['n_sig_blk'], len(T)))
print('alargamento do IC: mediana %.2fx (min %.2f, max %.2f)' %
      (meta['razao_largura_mediana'], meta['razao_largura_min'], meta['razao_largura_max']))
