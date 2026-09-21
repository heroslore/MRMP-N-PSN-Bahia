# -*- coding: utf-8 -*-
"""
Variabilidade interanual da PSN (subseção 6.6 da dissertação): soma anual da PSN por bioma
(anos completos 2001–2024; 2025 = janeiro a setembro, parcial), coeficiente de variação,
tendência monotônica (estimador de Sen + teste de Mann-Kendall), piores/melhores anos,
comparação com o NPP anual do MOD17A3HGF (npp_modis/resultados/npp_anual_2001_2025.csv)
e fase ENSO dominante de cada ano (>= 6 meses na fase, coluna `Enso` oficial da base).
Saídas: resultados_2001_2025/interanual/*.csv, interanual.json e figuras_dissertacao/fig15_interanual.png
"""
import os, json, numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(BASE, 'resultados_2001_2025')
OUT = os.path.join(RES, 'interanual'); os.makedirs(OUT, exist_ok=True)
FIG = os.path.join(BASE, 'figuras_dissertacao'); os.makedirs(FIG, exist_ok=True)
NPP_CSV = os.path.join(BASE, '..', 'npp_modis', 'resultados', 'npp_anual_2001_2025.csv')
NOME = {'MA': 'Mata Atlântica', 'CE': 'Cerrado', 'CA': 'Caatinga'}; COR = {'MA': '#2CA02C', 'CE': '#FF7F0E', 'CA': '#D62728'}
NPPCOL = {'MA': 'FMA_NPP', 'CE': 'Cerrado_NPP', 'CA': 'Caatinga_NPP'}
ANO_INI, ANO_FIM = 2001, 2024

d = pd.read_excel(os.path.join(BASE, 'Dados_base_nova_2001_2025.xlsx')); d.columns = d.columns.str.strip()
npp = pd.read_csv(NPP_CSV) if os.path.exists(NPP_CSV) else None

def mann_kendall(y):
    y = np.asarray(y, float); n = len(y)
    s = sum(np.sign(y[j] - y[i]) for i in range(n) for j in range(i + 1, n))
    var = n * (n - 1) * (2 * n + 5) / 18
    z = (s - np.sign(s)) / np.sqrt(var) if s != 0 else 0.0
    return 2 * (1 - stats.norm.cdf(abs(z))), z

# fase dominante do ano
fase_ano = {}
for ano, g in d.groupby('ANO'):
    vc = g['Enso'].value_counts(); f = vc.idxmax()
    fase_ano[int(ano)] = f if (vc.max() >= 6 and f != 'Neutro') else 'Neutro'

J = {'ano_ini': ANO_INI, 'ano_fim': ANO_FIM, 'fase_ano': fase_ano, 'biomas': {}}
rows, anual = [], []
for b in ['MA', 'CE', 'CA']:
    an = d.groupby('ANO')[f'NP_{b}'].sum(); full = an.loc[ANO_INI:ANO_FIM]; anos = full.index.values.astype(float)
    slope, intercept, lo, hi = stats.theilslopes(full.values, anos)
    p_mk, z_mk = mann_kendall(full.values)
    tend_pct = 100 * slope * (ANO_FIM - ANO_INI) / full.mean()
    n_meses_25 = int((d.ANO == 2025).sum()); m25 = int(d.loc[d.ANO == 2025, 'MÊS'].max())
    parcial25 = float(d.loc[d.ANO == 2025, f'NP_{b}'].sum())
    ref25 = float(d[(d.ANO <= ANO_FIM) & (d['MÊS'] <= m25)].groupby('ANO')[f'NP_{b}'].sum().mean())
    r_npp = p_npp = np.nan
    if npp is not None:
        m = npp[(npp.ano >= ANO_INI) & (npp.ano <= ANO_FIM)].set_index('ano')[NPPCOL[b]].reindex(full.index)
        ok = m.notna(); r_npp, p_npp = stats.pearsonr(full[ok], m[ok])
    piores = [int(a) for a in full.nsmallest(3).index]; melhores = [int(a) for a in full.nlargest(3).index]
    info = dict(media=float(full.mean()), dp=float(full.std(ddof=1)), cv=float(100 * full.std(ddof=1) / full.mean()),
                sen_g_ano=float(slope), sen_pct_periodo=float(tend_pct), sen_ic95_pct=[float(100 * lo * (ANO_FIM - ANO_INI) / full.mean()), float(100 * hi * (ANO_FIM - ANO_INI) / full.mean())],
                p_mk=float(p_mk), piores=piores, melhores=melhores, min_ano=int(full.idxmin()), min_val=float(full.min()), max_ano=int(full.idxmax()), max_val=float(full.max()),
                parcial_2025=parcial25, meses_2025=m25, ref_2025=ref25, dif_2025_pct=float(100 * (parcial25 / ref25 - 1)), r_npp=float(r_npp), p_npp=float(p_npp),
                npp_media=(float(m[ok].mean()) if npp is not None else None))
    J['biomas'][b] = info
    rows.append(dict(bioma=NOME[b], **{k: v_ for k, v_ in info.items() if not isinstance(v_, list)}, piores=', '.join(map(str, piores)), melhores=', '.join(map(str, melhores))))
    for a, val in an.items():
        anual.append(dict(bioma=NOME[b], ano=int(a), psn_anual=float(val), parcial=(int(a) > ANO_FIM), fase_dominante=fase_ano.get(int(a)),
                          npp_anual=(float(npp.set_index('ano')[NPPCOL[b]].get(int(a), np.nan)) if npp is not None else np.nan)))
T7 = pd.DataFrame(rows); T7.round(3).to_csv(os.path.join(OUT, 'tabela7_interanual.csv'), index=False)
AN = pd.DataFrame(anual); AN.round(2).to_csv(os.path.join(OUT, 'psn_anual_por_bioma.csv'), index=False)
json.dump(J, open(os.path.join(RES, 'interanual.json'), 'w', encoding='utf-8'), ensure_ascii=False, indent=1)

# ---------------- Figura 15: três painéis, legenda única no rodapé
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
plt.rcParams.update({'font.size': 12, 'axes.titlesize': 13.5, 'axes.labelsize': 12.5, 'xtick.labelsize': 11, 'ytick.labelsize': 11, 'legend.fontsize': 11})
fig, axes = plt.subplots(3, 1, figsize=(10, 12.5), sharex=True)
MESES = ['jan', 'fev', 'mar', 'abr', 'mai', 'jun', 'jul', 'ago', 'set', 'out', 'nov', 'dez']
for k, b in enumerate(['MA', 'CE', 'CA']):
    ax = axes[k]; q = AN[AN.bioma == NOME[b]]; full = q[~q.parcial]; par = q[q.parcial]; info = J['biomas'][b]
    for a, f in fase_ano.items():
        if f != 'Neutro': ax.axvspan(a - 0.5, a + 0.5, color='#D62728' if f == 'El Niño' else '#1F77B4', alpha=0.10, lw=0)
    ax.plot(full.ano, full.psn_anual, '-o', color=COR[b], lw=1.8, ms=6)
    xs = np.array([ANO_INI, ANO_FIM]); ax.plot(xs, info['media'] + info['sen_g_ano'] * (xs - full.ano.mean()), '--', color='k', lw=1.5)
    if len(par):
        ax.plot(par.ano, par.psn_anual, 'o', mfc='white', mec=COR[b], mew=2, ms=8)
        ax.annotate(f"2025\n(jan–{MESES[info['meses_2025'] - 1]})", (float(par.ano.iloc[0]), float(par.psn_anual.iloc[0])), textcoords='offset points', xytext=(-6, 8), ha='right', fontsize=9.5, color='#444')
    for a in info['piores']:
        ax.annotate(str(a), (a, float(full[full.ano == a].psn_anual.iloc[0])), textcoords='offset points', xytext=(0, -15), ha='center', fontsize=9.5, color='#444')
    sinal = '−' if info['sen_pct_periodo'] < 0 else '+'
    ax.set_title(f"({'abc'[k]}) {NOME[b]}: tendência de Sen {sinal}{abs(info['sen_pct_periodo']):.1f}% em {ANO_INI}–{ANO_FIM} "
                 f"(Mann-Kendall p = {info['p_mk']:.2f}); CV = {info['cv']:.1f}%".replace('.', ','), fontweight='bold', loc='left')
    ax.set_ylabel('PSN anual (g C·m⁻²·ano⁻¹)'); ax.grid(alpha=0.3); ax.spines[['top', 'right']].set_visible(False)
    ax.margins(y=0.12)
axes[-1].set_xlabel('Ano'); axes[-1].set_xticks(range(2001, 2026, 2))
handles = [Line2D([], [], color='#555', marker='o', lw=1.8, label='PSN anual (soma dos 12 meses)'),
           Line2D([], [], color='k', ls='--', lw=1.5, label='tendência de Sen (2001–2024)'),
           Line2D([], [], marker='o', mfc='white', mec='#555', mew=2, ls='', ms=8, label='2025 parcial (fora da tendência)'),
           Patch(color='#D62728', alpha=0.25, label='ano com El Niño dominante (≥ 6 meses)'),
           Patch(color='#1F77B4', alpha=0.25, label='ano com La Niña dominante (≥ 6 meses)')]
fig.legend(handles=handles, loc='lower center', ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.005))
fig.tight_layout(rect=(0, 0.05, 1, 1)); fig.savefig(os.path.join(FIG, 'fig15_interanual.png'), dpi=300, bbox_inches='tight'); plt.close(fig)
pd.set_option('display.width', 250); print(T7.round(2).to_string(index=False)); print('fases dominantes:', {a: f for a, f in fase_ano.items() if f != 'Neutro'})
