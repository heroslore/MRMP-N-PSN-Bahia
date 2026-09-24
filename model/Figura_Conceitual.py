# -*- coding: utf-8 -*-
"""
Esquema conceitual dos tres padroes de associacao (Artigo 2).

Cada linha resume, para um bioma, a fase do ENSO com que a PSN se associou, a
anomalia climatica observada nessa fase e a anomalia de produtividade que a
acompanha, com a assinatura temporal correspondente. As setas representam
padroes compativeis com os resultados, nao relacoes causais demonstradas.

Todos os valores sao os ja reportados no texto; nenhum e recalculado aqui.

Saida: figuras_artigo2/fig5_conceitual.png
"""
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

BASE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(BASE, 'figuras_artigo2')
os.makedirs(FIG, exist_ok=True)

COR = {'Mata Atlântica': '#009E73', 'Cerrado': '#E69F00', 'Caatinga': '#CC3311'}
TINTA, TINTA2 = '#1a1a1a', '#555555'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'text.color': TINTA})

LINHAS = [
    ('Mata Atlântica', 'El Niño',
     'temperatura da superfície\n$+3{,}1$ p.p. (≈ $+0{,}9$ °C)',
     'PSN $-5{,}4$ p.p.\nmeses extremos: 4,7 % → 22,4 %',
     'associação máxima\ncom 2 meses de defasagem'),
    ('Cerrado', 'La Niña',
     'evapotranspiração $+10{,}2$ p.p.\ndisponibilidade hídrica $+14{,}0$ p.p.',
     'PSN $+11{,}7$ p.p.\n(incerto por episódio: $p = 0{,}081$)',
     'associação sustentada\nao longo de 12 meses'),
    ('Caatinga', 'La Niña',
     'evapotranspiração $+10{,}9$ p.p.\ndisponibilidade hídrica $+12{,}5$ p.p.',
     'PSN $+10{,}4$ p.p.\ntrimestre chuvoso: $+25{,}0$ p.p.',
     'associação máxima em fase,\ndecaindo em poucos meses'),
]

fig, ax = plt.subplots(figsize=(11.2, 5.4))
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis('off')

X = [2, 23, 47, 73]           # inicio de cada coluna
W = [18, 21, 23, 25]
TOP = 84
ALT = 20
GAP = 6.5

for j, t in enumerate(['Bioma e fase', 'Anomalia climática observada',
                       'Anomalia de produtividade', 'Assinatura temporal']):
    ax.text(X[j] + W[j] / 2, TOP + 6, t, ha='center', va='bottom',
            fontsize=10.5, fontweight='bold', color=TINTA2)

for i, (bioma, fase, clima, psn, tempo) in enumerate(LINHAS):
    y = TOP - i * (ALT + GAP) - ALT
    c = COR[bioma]
    # coluna 1: bioma e fase
    ax.add_patch(FancyBboxPatch((X[0], y), W[0], ALT, boxstyle='round,pad=0.6,rounding_size=1.5',
                                facecolor=c, edgecolor='none', alpha=0.16, zorder=1))
    ax.text(X[0] + W[0] / 2, y + ALT * 0.62, bioma, ha='center', va='center',
            fontsize=11, fontweight='bold', color=c)
    ax.text(X[0] + W[0] / 2, y + ALT * 0.26, fase, ha='center', va='center',
            fontsize=10.5, color=TINTA)
    # colunas 2 a 4
    for j, txt in [(1, clima), (2, psn), (3, tempo)]:
        ax.add_patch(FancyBboxPatch((X[j], y), W[j], ALT, boxstyle='round,pad=0.6,rounding_size=1.5',
                                    facecolor='#f4f4f4' if j < 3 else 'white',
                                    edgecolor='#d0d0d0' if j == 3 else 'none', lw=0.9, zorder=1))
        ax.text(X[j] + W[j] / 2, y + ALT / 2, txt, ha='center', va='center',
                fontsize=9.5, color=TINTA, linespacing=1.5)
    # setas
    for j in range(3):
        ax.add_patch(FancyArrowPatch((X[j] + W[j] + 0.6, y + ALT / 2), (X[j + 1] - 0.6, y + ALT / 2),
                                     arrowstyle='-|>', mutation_scale=13, lw=1.6,
                                     color=c, zorder=3))

ax.text(50, 2.5, 'As setas indicam padrões de associação compatíveis com os resultados, '
                 'não relações causais demonstradas.',
        ha='center', va='bottom', fontsize=9.2, style='italic', color=TINTA2)

fig.tight_layout()
fig.savefig(os.path.join(FIG, 'fig5_conceitual.png'), dpi=300, bbox_inches='tight')
print('figura em', os.path.join(FIG, 'fig5_conceitual.png'))
