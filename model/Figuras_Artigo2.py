# -*- coding: utf-8 -*-
"""
Figuras do Artigo 2 (ENSO / unidade amostral).

  fig1_episodios.png   anomalia mensal de PSN com os episodios ENSO sombreados
  fig2_unidade.png     delta% por fase: IC do bootstrap i.i.d. (mes) vs em blocos (episodio)
  fig3_defasagem.png   Spearman ONI(t) x anomalia PSN(t+L)

Paleta verificada para daltonismo (Okabe-Ito adaptada): o par verde/laranja da
dissertacao tem separacao quase nula sob protanopia e foi substituido.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

BASE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(BASE, 'resultados_2001_2025')
ENS = os.path.join(RES, 'enso_anomalias')
OUT = os.path.join(BASE, 'figuras_artigo2')
os.makedirs(OUT, exist_ok=True)

COR_BIOMA = {'Mata Atlântica': '#009E73', 'Cerrado': '#E69F00', 'Caatinga': '#CC3311'}
COR_FASE = {'El Niño': '#CC3311', 'La Niña': '#0077BB'}
TINTA, TINTA2, GRADE = '#1a1a1a', '#555555', '#d8d8d8'
BIOMAS = ['Mata Atlântica', 'Cerrado', 'Caatinga']
SIGLA = {'Mata Atlântica': 'MA', 'Cerrado': 'CE', 'Caatinga': 'CA'}

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
                     'axes.edgecolor': '#999999', 'axes.linewidth': 0.8,
                     'xtick.color': TINTA2, 'ytick.color': TINTA2,
                     'axes.labelcolor': TINTA, 'text.color': TINTA})

def limpar(ax, grade_y=True):
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    if grade_y:
        ax.grid(axis='y', color=GRADE, lw=0.7, zorder=0); ax.set_axisbelow(True)

d = pd.read_excel(os.path.join(BASE, 'Dados_base_nova_2001_2025.xlsx'))
d.columns = d.columns.str.strip()
d['data'] = pd.to_datetime(dict(year=d['ANO'], month=d['MÊS'], day=1))
fase = d['Enso'].astype(str).values

def anom_pct(col):
    clim = d.groupby('MÊS')[col].transform('mean')
    return (100 * (d[col] - clim) / clim).values

# corridas contiguas por fase, para sombrear
corridas = []
ini = 0
for i in range(1, len(d) + 1):
    if i == len(d) or fase[i] != fase[ini]:
        corridas.append((fase[ini], d['data'].iloc[ini], d['data'].iloc[i - 1], i - ini))
        ini = i

# ------------------------------------------------------------------ Figura 1
fig, axes = plt.subplots(3, 1, figsize=(9.2, 5.9), sharex=True)
for ax, b in zip(axes, BIOMAS):
    a = anom_pct(f'NP_{SIGLA[b]}')
    for f, t0, t1, n in corridas:
        if f in COR_FASE:
            ax.axvspan(t0, t1, color=COR_FASE[f], alpha=0.14, lw=0, zorder=0)
    ax.axhline(0, color='#999999', lw=0.9, zorder=1)
    ax.plot(d['data'], a, color=TINTA, lw=1.3, zorder=3)
    limpar(ax)
    ax.set_ylabel('anomalia de PSN (%)', fontsize=10)
    ax.text(0.006, 0.93, b, transform=ax.transAxes, fontsize=11.5,
            fontweight='bold', va='top', color=COR_BIOMA[b])
axes[-1].set_xlabel('')
axes[0].legend(handles=[Patch(facecolor=COR_FASE['El Niño'], alpha=0.35, label='El Niño'),
                        Patch(facecolor=COR_FASE['La Niña'], alpha=0.35, label='La Niña')],
               loc='upper right', frameon=False, ncol=2, fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(OUT, 'fig1_episodios.png'), dpi=300, bbox_inches='tight')
plt.close(fig)

# ------------------------------------------------------------------ Figura 2
T = pd.read_csv(os.path.join(ENS, 'bootstrap_unidade_amostral.csv'))
P = T[T.variavel == 'PSN'].copy()
fig, ax = plt.subplots(figsize=(8.6, 4.6))
linhas = [(b, f) for b in BIOMAS for f in ['El Niño', 'La Niña']]
ypos = {k: len(linhas) - 1 - i for i, k in enumerate(linhas)}
for (b, f), y in ypos.items():
    r = P[(P.bioma == b) & (P.fase == f)].iloc[0]
    c = COR_FASE[f]
    ax.plot([r.iid_inf, r.iid_sup], [y + 0.17] * 2, color=c, lw=1.6, alpha=0.5,
            solid_capstyle='round', zorder=2)
    ax.plot([r.delta_pct], [y + 0.17], 'o', mfc='white', mec=c, mew=1.6, ms=7, zorder=3)
    ax.plot([r.blk_inf, r.blk_sup], [y - 0.17] * 2, color=c, lw=3.2,
            solid_capstyle='round', zorder=2)
    ax.plot([r.delta_pct], [y - 0.17], 's', color=c, ms=7.5, mec='white', mew=1.1, zorder=3)
    ax.text(25.6, y, f'{r.n_blocos} ep.', fontsize=9.5, color=TINTA2, va='center', ha='right')
ax.axvline(0, color='#777777', lw=1.1, zorder=1)
ax.set_yticks(list(ypos.values()))
ax.set_yticklabels([f'{b}\n{f}' for b, f in ypos], fontsize=10)
for t, (b, f) in zip(ax.get_yticklabels(), ypos):
    t.set_color(TINTA)
ax.set_xlabel('anomalia de PSN em relação aos meses neutros (%)')
ax.set_xlim(-30, 27)
ax.set_ylim(-0.75, len(linhas) - 0.25)
limpar(ax, grade_y=False)
ax.grid(axis='x', color=GRADE, lw=0.7, zorder=0); ax.set_axisbelow(True)
ax.legend(handles=[Line2D([], [], color='#777777', lw=1.6, alpha=0.6, marker='o',
                          mfc='white', mec='#777777', ms=7, label='mês como unidade (i.i.d.)'),
                   Line2D([], [], color='#777777', lw=3.2, marker='s', ms=7.5,
                          mec='white', label='episódio como unidade (blocos)')],
          loc='lower left', frameon=False, fontsize=9.8)
fig.tight_layout()
fig.savefig(os.path.join(OUT, 'fig2_unidade.png'), dpi=300, bbox_inches='tight')
plt.close(fig)

# ------------------------------------------------------------------ Figura 3
L = pd.read_csv(os.path.join(ENS, 'spearman_ONI_x_PSN_lags.csv'))
fig, ax = plt.subplots(figsize=(8.6, 4.2))
for b in BIOMAS:
    q = L[L.bioma == b].sort_values('lag')
    c = COR_BIOMA[b]
    ax.plot(q.lag, q.rho, '-', color=c, lw=2.2, zorder=3)
    sig = q[q.p < 0.05]
    ax.plot(sig.lag, sig.rho, 'o', color=c, ms=7, mec='white', mew=1.2, zorder=4)
    ns = q[q.p >= 0.05]
    ax.plot(ns.lag, ns.rho, 'o', mfc='white', mec=c, mew=1.6, ms=6, zorder=4)
    ult = q.iloc[-1]
    ax.text(ult.lag + 0.22, ult.rho, b, color=c, fontsize=10.5, va='center', fontweight='bold')
ax.axhline(0, color='#999999', lw=0.9, zorder=1)
ax.set_xlabel('defasagem $L$ (meses) entre o ONI e a anomalia de PSN')
ax.set_ylabel('ρ de Spearman')
ax.set_xlim(-0.4, 16.4)
ax.set_xticks(range(0, 13, 2))
limpar(ax)
ax.legend(handles=[Line2D([], [], color='#777777', marker='o', ms=7, ls='none',
                          mec='white', mew=1.2, label='p < 0,05'),
                   Line2D([], [], color='#777777', marker='o', ms=6, ls='none',
                          mfc='white', mec='#777777', mew=1.6, label='p ≥ 0,05')],
          loc='lower right', frameon=False, fontsize=9.8, ncol=2)
fig.tight_layout()
fig.savefig(os.path.join(OUT, 'fig3_defasagem.png'), dpi=300, bbox_inches='tight')
plt.close(fig)

print('figuras em', OUT)
for f in sorted(os.listdir(OUT)):
    print(' ', f, os.path.getsize(os.path.join(OUT, f)) // 1024, 'KB')
