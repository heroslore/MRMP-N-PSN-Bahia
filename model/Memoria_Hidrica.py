# -*- coding: utf-8 -*-
"""
Janela de integracao hidrica de cada bioma (Artigo 2).

A defasagem entre o ONI e a anomalia de PSN descreve o tempo de resposta ao
forcante remoto, mas com correlacoes fracas. Aqui a pergunta e outra e mais
proxima da vegetacao: sobre que janela de tempo a chuva acumulada se associa
melhor a anomalia de produtividade de cada bioma? A janela de maior associacao
e uma medida empirica de quanto tempo o sistema integra a agua recebida.

Para cada janela w de 1 a 12 periodos calcula-se a chuva acumulada nos w
periodos terminados em t, sua anomalia em relacao a climatologia da mesma
janela, e a correlacao de Spearman com a anomalia de PSN em t.

Os meses sao tratados como independentes, como na analise de defasagem; os
coeficientes descrevem a forma da curva e nao sustentam teste individual.

Saidas: resultados_2001_2025/enso_anomalias/memoria_hidrica.csv
        figuras_artigo2/fig3_resposta_temporal.png
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

BASE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE, 'resultados_2001_2025', 'enso_anomalias')
FIG = os.path.join(BASE, 'figuras_artigo2')
os.makedirs(FIG, exist_ok=True)

NOME = {'MA': 'Mata Atlântica', 'CE': 'Cerrado', 'CA': 'Caatinga'}
JANELAS = list(range(1, 13))

d = pd.read_excel(os.path.join(BASE, 'Dados_base_nova_2001_2025.xlsx'))
d.columns = d.columns.str.strip()

def anom_pct(col):
    clim = d.groupby('MÊS')[col].transform('mean')
    return (100 * (d[col] - clim) / clim).values

def acumulado(v, w):
    return np.array([v[max(0, i - w + 1):i + 1].sum() for i in range(len(v))])

linhas = []
for b in ['MA', 'CE', 'CA']:
    psn = anom_pct(f'NP_{b}')
    pre = d[f'PRE_{b}'].values
    clim = d.groupby('MÊS')[f'PRE_{b}'].transform('mean').values
    for w in JANELAS:
        a_pre = 100 * (acumulado(pre, w) - acumulado(clim, w)) / acumulado(clim, w)
        rho, p = stats.spearmanr(a_pre[w:], psn[w:])
        linhas.append(dict(bioma=NOME[b], janela_meses=w, rho=round(float(rho), 3),
                           p=float('%.3g' % p), n=len(psn) - w))

M = pd.DataFrame(linhas)
M.to_csv(os.path.join(OUT, 'memoria_hidrica.csv'), index=False, encoding='utf-8')
piv = M.pivot(index='bioma', columns='janela_meses', values='rho')
pd.set_option('display.width', 200)
print(piv.to_string())
print('\njanela de maior |rho| por bioma:')
for b in NOME.values():
    s = M[M.bioma == b]
    print(f"  {b:<16} {int(s.loc[s.rho.abs().idxmax(), 'janela_meses']):>2} meses  (rho = {s.rho.max():+.2f})")

# ------------------------------------------------------------------ figura
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

COR = {'Mata Atlântica': '#009E73', 'Cerrado': '#E69F00', 'Caatinga': '#CC3311'}
TINTA, TINTA2, GRADE = '#1a1a1a', '#555555', '#d8d8d8'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
                     'axes.edgecolor': '#999999', 'axes.linewidth': 0.8,
                     'xtick.color': TINTA2, 'ytick.color': TINTA2,
                     'axes.labelcolor': TINTA, 'text.color': TINTA})

L = pd.read_csv(os.path.join(OUT, 'spearman_ONI_x_PSN_lags.csv'))
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.2, 4.2))

for b in NOME.values():
    q = L[L.bioma == b].sort_values('lag'); c = COR[b]
    ax1.plot(q.lag, q.rho, '-', color=c, lw=2.2, zorder=3)
    sig = q[q.p < 0.05]; ns = q[q.p >= 0.05]
    ax1.plot(sig.lag, sig.rho, 'o', color=c, ms=6.5, mec='white', mew=1.2, zorder=4)
    ax1.plot(ns.lag, ns.rho, 'o', mfc='white', mec=c, mew=1.5, ms=5.5, zorder=4)
ax1.axhline(0, color='#999999', lw=0.9, zorder=1)
ax1.set_xlabel('defasagem $L$ (meses) entre o ONI\ne a anomalia de PSN', fontsize=10)
ax1.set_ylabel('ρ de Spearman')
ax1.set_title('(a) Resposta ao forçante remoto', fontsize=11.5, pad=8)
ax1.set_xticks(range(0, 13, 3)); ax1.set_ylim(-0.72, 0.72)

for b in NOME.values():
    q = M[M.bioma == b].sort_values('janela_meses'); c = COR[b]
    ax2.plot(q.janela_meses, q.rho, '-o', color=c, lw=2.2, ms=6, mec='white', mew=1.1, zorder=3)
    i = q.rho.abs().idxmax()
    ax2.plot([M.loc[i, 'janela_meses']], [M.loc[i, 'rho']], marker='o', color=c,
             ms=12, mfc='none', mew=1.8, zorder=4)
ax2.axhline(0, color='#999999', lw=0.9, zorder=1)
ax2.set_xlabel('janela de acumulação da chuva\nanterior (meses)', fontsize=10)
ax2.set_title('(b) Associação com a água acumulada', fontsize=11.5, pad=8)
ax2.set_xticks(range(0, 13, 3)); ax2.set_ylim(-0.72, 0.72)

for ax in (ax1, ax2):
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    ax.grid(axis='y', color=GRADE, lw=0.7, zorder=0); ax.set_axisbelow(True)

fig.legend(handles=[Line2D([], [], color=COR[b], lw=2.4, marker='o', ms=6,
                           mec='white', label=b) for b in NOME.values()],
           loc='lower center', ncol=3, frameon=False, fontsize=10, bbox_to_anchor=(0.5, -0.05))
fig.tight_layout(rect=[0, 0.07, 1, 1])
fig.savefig(os.path.join(FIG, 'fig3_resposta_temporal.png'), dpi=300, bbox_inches='tight')
print('\nfigura em', os.path.join(FIG, 'fig3_resposta_temporal.png'))
