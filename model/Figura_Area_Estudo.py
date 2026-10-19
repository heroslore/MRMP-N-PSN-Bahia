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
from PIL import Image, ImageDraw, ImageFont
import matplotlib.font_manager as fm

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
# A linha de fonte embutida no raster do mapa e reescrita com virgula entre "IBGE" e
# "Biomas (2019)". Somente a faixa dessa linha (linhas 955 a 977 do
# raster 1261x984, x de 180 a 1080, so texto sobre fundo branco) e apagada e o
# texto e reescrito em DejaVu Sans, mesma cor cinza, mesma posicao, com virgula.
# Nenhum pixel do mapa (areas, contornos, legenda, eixos) e alterado.
mapa_img = Image.open(MAPA).convert('RGB')
_arr = np.asarray(mapa_img).copy()
assert _arr.shape == (984, 1261, 3)
assert _arr[955:978, :180].min() > 240 and _arr[955:978, 1080:].min() > 240
_arr[955:978, 180:1080] = 255
mapa_img = Image.fromarray(_arr)
FONTE_TXT = ('Fonte: IBGE, Biomas (2019) e limite estadual (2020). Datum SIRGAS 2000 '
             '(EPSG:4674). Contorno da América do Sul: Natural Earth.')
ImageDraw.Draw(mapa_img).text(
    (193, 970), FONTE_TXT, fill=(87, 87, 87), anchor='ls',
    font=ImageFont.truetype(fm.findfont('DejaVu Sans'), 12.6))
ax1.imshow(np.asarray(mapa_img))
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
# rotulo na base da faixa N-D, regiao sem curvas
ax2.text(11.5, 4.0, 'trimestre\nchuvoso', fontsize=9,
         color=TINTA2, ha='center', va='bottom', linespacing=1.3)
ax2.set_xticks(range(1, 13)); ax2.set_xticklabels(MESES)
ax2.set_xlim(0.4, 12.6)
ax2.set_ylabel('precipitação média do período (mm)')
ax2.set_title('(b) Climatologia de precipitação, 2001 a 2025', fontsize=11.5, pad=8)
ax2.spines['top'].set_visible(False); ax2.spines['right'].set_visible(False)
ax2.grid(axis='y', color=GRADE, lw=0.7, zorder=0); ax2.set_axisbelow(True)
ax2.legend(handles=[Line2D([], [], color=COR[n], lw=2.4, marker='o', ms=5.5, mec='white',
                           label=f'{n} ({anual[n]:.0f} mm nos 12 períodos)') for n in NOME.values()],
           loc='upper center', frameon=False, fontsize=9.5)
fig.tight_layout()
fig.savefig(os.path.join(FIG, 'fig1_area_estudo.png'), dpi=300, bbox_inches='tight')
print('precipitação anual climatológica (mm):', {k: round(v) for k, v in anual.items()})
print('figura em', os.path.join(FIG, 'fig1_area_estudo.png'))
# registro dos valores plotados
_clim = pd.DataFrame({b: d.groupby('MÊS')[f'PRE_{b}'].mean() for b in ['MA', 'CE', 'CA']})
_clim.to_csv(os.path.join(FIG, 'valores_plotados_fig1.csv'), float_format='%.6f')
print('climatologia mensal (mm):')
print(_clim.round(3).to_string())
print('anual exato (mm):', {k: round(v, 3) for k, v in anual.items()})
