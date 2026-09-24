# -*- coding: utf-8 -*-
"""
Confinamento sazonal da associacao entre as fases do ENSO e a PSN (Artigo 2).

Motivacao fisica, definida antes de olhar o desfecho: o ENSO modula a estacao
chuvosa, e em sistemas limitados pela agua e nessa estacao que um excedente
hidrico pode se converter em produtividade. Agregar o ano inteiro mistura os
meses em que esse canal opera com os meses em que ele nao opera.

A analise repete o teste por fase da Secao 2.4, com o EPISODIO como unidade
(bootstrap em blocos, como na Secao 2.8), restrito a dois trimestres definidos
pela climatologia de precipitacao de cada bioma:
  - trimestre chuvoso: os tres meses civis consecutivos de maior precipitacao media;
  - trimestre seco:    os tres meses civis consecutivos de menor precipitacao media.
O trimestre seco funciona como controle negativo: se a associacao decorre do
canal hidrico da estacao chuvosa, ela deve estar ausente nele.

Correcao para comparacoes multiplas por Benjamini-Hochberg em duas familias:
  (i)  os 6 testes de cada trimestre (3 biomas x 2 fases), que e a familia primaria;
  (ii) os 12 testes dos dois trimestres, como analise de sensibilidade.

Saidas: resultados_2001_2025/enso_anomalias/estratificacao_sazonal.csv
        figuras_artigo2/fig3_sazonal.png
"""
import os
import numpy as np
import pandas as pd

BASE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE, 'resultados_2001_2025', 'enso_anomalias')
FIG = os.path.join(BASE, 'figuras_artigo2')
os.makedirs(OUT, exist_ok=True); os.makedirs(FIG, exist_ok=True)

NOME = {'MA': 'Mata Atlântica', 'CE': 'Cerrado', 'CA': 'Caatinga'}
ATIVAS = ['El Niño', 'La Niña']
B = 20000
SEED = 42

d = pd.read_excel(os.path.join(BASE, 'Dados_base_nova_2001_2025.xlsx'))
d.columns = d.columns.str.strip()
fase = d['Enso'].astype(str).values

# blocos = corridas contiguas de meses na mesma fase oficial (o episodio)
bloco = np.zeros(len(d), int); k = 0
for i in range(1, len(d)):
    if fase[i] != fase[i - 1]:
        k += 1
    bloco[i] = k

def anom_pct(col):
    clim = d.groupby('MÊS')[col].transform('mean')
    return (100 * (d[col] - clim) / clim).values

def trimestres(b):
    """Trimestres civis consecutivos de maior e de menor precipitacao media."""
    pre = d[f'PRE_{b}'].groupby(d['MÊS']).mean()
    soma = {m: sum(pre[(m + i - 1) % 12 + 1] for i in range(3)) for m in range(1, 13)}
    chuv = max(soma, key=soma.get); seco = min(soma, key=soma.get)
    meses = lambda t: [(t + i - 1) % 12 + 1 for i in range(3)]
    return ('chuvoso', meses(chuv)), ('seco', meses(seco))

def bh(p, q=0.05):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); ps = p[o]
    sig = np.zeros(n, bool)
    kk = np.where(ps <= (np.arange(1, n + 1) / n) * q)[0]
    if len(kk):
        sig[o[:kk.max() + 1]] = True
    adj = np.minimum.accumulate((ps * n / np.arange(1, n + 1))[::-1])[::-1]
    qv = np.empty(n); qv[o] = np.minimum(adj, 1.0)
    return sig, qv

rng = np.random.default_rng(SEED)
linhas = []
for b in ['MA', 'CE', 'CA']:
    a = anom_pct(f'NP_{b}')
    for rot, meses in trimestres(b):
        sel = d['MÊS'].isin(meses).values
        bl_N = [a[(bloco == u) & sel] for u in sorted(set(bloco[fase == 'Neutro']))
                if ((bloco == u) & sel).sum() > 0]
        for f in ATIVAS:
            bl_F = [a[(bloco == u) & sel] for u in sorted(set(bloco[fase == f]))
                    if ((bloco == u) & sel).sum() > 0]
            obs = np.concatenate(bl_F).mean() - np.concatenate(bl_N).mean()
            dist = np.empty(B)
            for j in range(B):
                rF = np.concatenate([bl_F[i] for i in rng.integers(0, len(bl_F), len(bl_F))])
                rN = np.concatenate([bl_N[i] for i in rng.integers(0, len(bl_N), len(bl_N))])
                dist[j] = rF.mean() - rN.mean()
            lo, hi = np.percentile(dist, [2.5, 97.5])
            p = min(max(2 * min((dist <= 0).mean(), (dist >= 0).mean()), 1 / B), 1.0)
            linhas.append(dict(bioma=NOME[b], trimestre=rot, meses='-'.join(map(str, meses)),
                               fase=f, n_meses=int(sel.sum() & 0xFFFFFFFF) if False else int(np.sum(sel & (fase == f))),
                               n_episodios=len(bl_F), delta_pct=round(float(obs), 2),
                               ic_inf=round(float(lo), 2), ic_sup=round(float(hi), 2),
                               p=round(float(p), 4)))

T = pd.DataFrame(linhas)
# (i) familia primaria: os 6 testes de cada trimestre
for rot in ['chuvoso', 'seco']:
    m = T.trimestre == rot
    sig, qv = bh(T.loc[m, 'p'].values)
    T.loc[m, 'q_trimestre'] = np.round(qv, 4); T.loc[m, 'bh_trimestre'] = sig
# (ii) sensibilidade: os 12 testes
sig12, q12 = bh(T['p'].values)
T['q_12testes'] = np.round(q12, 4); T['bh_12testes'] = sig12

T.to_csv(os.path.join(OUT, 'estratificacao_sazonal.csv'), index=False, encoding='utf-8')
pd.set_option('display.width', 220)
print(T.to_string(index=False))
print('\nsobrevivem ao BH na familia do trimestre:', int(T.bh_trimestre.sum()), 'de', len(T))
print('sobrevivem ao BH na familia dos 12 testes :', int(T.bh_12testes.sum()), 'de', len(T))

# ------------------------------------------------------------------ figura
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

COR_FASE = {'El Niño': '#CC3311', 'La Niña': '#0077BB'}
TINTA, TINTA2, GRADE = '#1a1a1a', '#555555', '#d8d8d8'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
                     'axes.edgecolor': '#999999', 'axes.linewidth': 0.8,
                     'xtick.color': TINTA2, 'ytick.color': TINTA2,
                     'axes.labelcolor': TINTA, 'text.color': TINTA})
BIOMAS = ['Mata Atlântica', 'Cerrado', 'Caatinga']
fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.0), sharex=True)
for ax, rot, tit in zip(axes, ['chuvoso', 'seco'], ['Trimestre chuvoso', 'Trimestre seco']):
    linhas_y = [(b, f) for b in BIOMAS for f in ATIVAS]
    ypos = {k: len(linhas_y) - 1 - i for i, k in enumerate(linhas_y)}
    for (b, f), y in ypos.items():
        r = T[(T.bioma == b) & (T.trimestre == rot) & (T.fase == f)].iloc[0]
        c = COR_FASE[f]
        ax.plot([r.ic_inf, r.ic_sup], [y] * 2, color=c, lw=3.0, solid_capstyle='round', zorder=2)
        ax.plot([r.delta_pct], [y], 's', color=c, ms=7.5, mec='white', mew=1.1, zorder=3)
        if r.bh_trimestre:
            ax.plot([r.delta_pct], [y + 0.34], marker='*', color=c, ms=10, zorder=4)
    ax.axvline(0, color='#777777', lw=1.1, zorder=1)
    ax.set_yticks(list(ypos.values()))
    ax.set_yticklabels([f'{b}\n{f}' for b, f in ypos] if rot == 'chuvoso' else [''] * len(ypos), fontsize=9.5)
    ax.set_title(tit, fontsize=11.5, pad=8)
    ax.set_xlabel('anomalia de PSN em relação\naos meses neutros (%)', fontsize=10)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    ax.grid(axis='x', color=GRADE, lw=0.7, zorder=0); ax.set_axisbelow(True)
    ax.set_ylim(-0.7, len(linhas_y) - 0.3)
fig.legend(handles=[Line2D([], [], color=COR_FASE['El Niño'], lw=3, label='El Niño'),
                    Line2D([], [], color=COR_FASE['La Niña'], lw=3, label='La Niña'),
                    Line2D([], [], color='#777777', marker='*', ls='none', ms=10,
                           label='sobrevive à correção de Benjamini-Hochberg')],
           loc='lower center', ncol=3, frameon=False, fontsize=9.5, bbox_to_anchor=(0.5, -0.04))
fig.tight_layout(rect=[0, 0.05, 1, 1])
fig.savefig(os.path.join(FIG, 'fig3_sazonal.png'), dpi=300, bbox_inches='tight')
print('\nfigura em', os.path.join(FIG, 'fig3_sazonal.png'))
