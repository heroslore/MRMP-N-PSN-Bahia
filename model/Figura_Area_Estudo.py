# -*- coding: utf-8 -*-
"""
Figura da area de estudo (Artigo 2): mapa dos biomas e climatologia de precipitacao.

O painel (a) reaproveita o mapa dos tres biomas na Bahia usado no artigo
companheiro, recolorido para a paleta acessivel das demais figuras
(Recolorir_Mapa.py). O painel (b) mostra a precipitacao media de cada mes civil por
bioma, calculada da propria serie (GPM IMERG Final V07, 2001-2025), com o
trimestre climatologicamente mais chuvoso sombreado. Ele torna visiveis as
duas dimensoes do gradiente: o volume anual e a concentracao sazonal.

Saida: figuras_artigo2/fig1_area_estudo.png
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from PIL import Image

BASE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(BASE, 'figuras_artigo2')
MAPA = os.path.join(BASE, 'figuras_artigo2', 'mapa_biomas_acessivel.png')
os.makedirs(FIG, exist_ok=True)

NOME = {'MA': 'Mata Atlântica', 'CE': 'Cerrado', 'CA': 'Caatinga'}
COR = {'Mata Atlântica': '#009E73', 'Cerrado': '#E69F00', 'Caatinga': '#CC3311'}
TINTA, TINTA2, GRADE = '#1a1a1a', '#555555', '#d8d8d8'
MESES = ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']

d = pd.read_excel(os.path.join(BASE, 'Dados_base_nova_2001_2025.xlsx'))
d.columns = d.columns.str.strip()

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
                     'axes.edgecolor': '#999999', 'axes.linewidth': 0.8,
                     'xtick.color': TINTA2, 'ytick.color': TINTA2,
                     'axes.labelcolor': TINTA, 'text.color': TINTA})

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.0, 4.6),
                               gridspec_kw={'width_ratios': [1.0, 1.25]})

# (a) mapa
ax1.imshow(np.asarray(Image.open(MAPA).convert('RGB')))
ax1.set_axis_off()
ax1.set_title('(a) Biomas da Bahia', fontsize=11.5, pad=8)

# (b) climatologia de precipitacao
anual = {}
for b in ['MA', 'CE', 'CA']:
    clim = d.groupby('MÊS')[f'PRE_{b}'].mean()
    anual[NOME[b]] = clim.sum()
    ax2.plot(range(1, 13), clim.values, '-o', color=COR[NOME[b]], lw=2.2, ms=5.5,
             mec='white', mew=1.0, zorder=3, label=NOME[b])
# trimestre chuvoso comum aos tres biomas (nov-dez-jan)
for x0, x1 in [(0.5, 1.5), (10.5, 12.5)]:
    ax2.axvspan(x0, x1, color='#0077BB', alpha=0.08, lw=0, zorder=0)
ax2.text(12.0, ax2.get_ylim()[1] * 0.96, 'trimestre\nchuvoso', fontsize=9,
         color=TINTA2, ha='right', va='top', linespacing=1.3)
ax2.set_xticks(range(1, 13)); ax2.set_xticklabels(MESES)
ax2.set_xlim(0.4, 12.6)
ax2.set_ylabel('precipitação média do período (mm)')
ax2.set_title('(b) Climatologia de precipitação, 2001–2025', fontsize=11.5, pad=8)
ax2.spines['top'].set_visible(False); ax2.spines['right'].set_visible(False)
ax2.grid(axis='y', color=GRADE, lw=0.7, zorder=0); ax2.set_axisbelow(True)
ax2.legend(handles=[Line2D([], [], color=COR[n], lw=2.4, marker='o', ms=5.5, mec='white',
                           label=f'{n} ({anual[n]:.0f} mm nos 12 períodos)') for n in NOME.values()],
           loc='upper center', frameon=False, fontsize=9.5)
fig.tight_layout()
fig.savefig(os.path.join(FIG, 'fig1_area_estudo.png'), dpi=300, bbox_inches='tight')
print('precipitação anual climatológica (mm):', {k: round(v) for k, v in anual.items()})
print('figura em', os.path.join(FIG, 'fig1_area_estudo.png'))
