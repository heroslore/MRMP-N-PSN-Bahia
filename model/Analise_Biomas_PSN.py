# -*- coding: utf-8 -*-
"""
Análise do ENSO sobre a Fotossíntese Líquida (PSN) e variáveis ambientais.
Bahia | Mata Atlântica, Cerrado e Caatinga | 2001-2025

Gera as figuras de distribuição por fase ENSO (uma por bioma, com os painéis
PSN, EV, PRE e TST), o teste de Kruskal-Wallis por fase, a correlação de Pearson
ONI x variáveis, a regressão linear simples ONI -> variáveis climáticas, e a
dispersão PSN x variáveis ambientais.

-------------------------------------------------------------------------------
BASE DE DADOS DE ENTRADA
-------------------------------------------------------------------------------
Basta colocar o CSV bruto atualizado (base_final_2001_2025_plan1_excel_ptbr.csv)
na mesma pasta do script e rodar — o script detecta o CSV sozinho, renomeia as
colunas para o padrão interno, remove as linhas com precipitação ainda não
publicada pela NASA, calcula a fase ENSO (La Niña / Neutro / El Niño) a partir
do ONI, e salva o Excel (Dados_base_nova_2001_2025.xlsx) que o script usa —
o MESMO arquivo que o Modelo_PSN.py também lê, então rodar qualquer um dos
dois primeiro já prepara a base para o outro. Se apagar o CSV e deixar só o
.xlsx pronto, o script usa ele direto.

A fase ENSO segue a definição OFICIAL da NOAA/CPC (módulo enso_noaa.py):
El Niño / La Niña só quando o ONI fica >= +0,5 / <= -0,5 por pelo menos 5
trimestres móveis consecutivos; caso contrário, Neutro. A tabela completa da
NOAA (oni_noaa_cpc.txt) é usada para contar corretamente as sequências nas
bordas do período. O sombreamento do gráfico de série temporal continua
mostrando o limiar simples de +/-0,5 apenas como referência visual.
"""

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # backend não-interativo
import matplotlib.pyplot as plt
from scipy.stats import kruskal, pearsonr
from scipy import stats
import statsmodels.api as sm
from enso_noaa import classificar_fase_enso_noaa

BASE_DIR   = os.path.dirname(os.path.abspath(__file__))
PASTA_SAIDA = os.path.join(BASE_DIR, 'saidas_figuras')
os.makedirs(PASTA_SAIDA, exist_ok=True)

def caminho(nome):
    return os.path.join(BASE_DIR, nome)

def saida(nome):
    """Nome fixo dentro de saidas_figuras/ — cada execução sobrescreve a anterior."""
    return os.path.join(PASTA_SAIDA, nome)

# =============================================================================
# PREPARAÇÃO DA BASE A PARTIR DO CSV BRUTO (opcional — só roda se o CSV existir)
# Mesma lógica e o MESMO arquivo de saída usados pelo Modelo_PSN.py.
# =============================================================================

NOME_CSV_BRUTO   = 'base_final_2001_2025_plan1_excel_ptbr.csv'
NOME_XLSX_MODELO = 'Dados_base_nova_2001_2025.xlsx'

RENOME_COLUNAS_CSV = {
    'ano':                   'ANO',
    'mes':                   'MÊS',
    'ONI':                   'ONI',
    'FMA_PSN':               'NP_MA',
    'Cerrado_PSN':           'NP_CE',
    'Caatinga_PSN':          'NP_CA',
    'FMA_IDA':               'WAI_MA',
    'Cerrado_IDA':           'WAI_CE',
    'Caatinga_IDA':          'WAI_CA',
    'FMA_Evap':              'EV_MA',
    'Cerrado_Evap':          'EV_CE',
    'Caatinga_Evap':         'EV_CA',
    'FMA_Precip':            'PRE_MA',
    'Cerrado_Precip':        'PRE_CE',
    'Caatinga_Precip':       'PRE_CA',
    'FMA_Temp':              'TST_MA',
    'Cerrado_Temp':          'TST_CE',
    'Caatinga_Temp':         'TST_CA',
    'FMA_AreaQueimada':      'BURN_MA',
    'Cerrado_AreaQueimada':  'BURN_CE',
    'Caatinga_AreaQueimada': 'BURN_CA',
    # FMA_PET, Cerrado_PET, Caatinga_PET não são usadas — mantidas e ignoradas.
}

COLUNAS_PRECIP = ['PRE_MA', 'PRE_CE', 'PRE_CA']


def preparar_base_a_partir_do_csv():
    """
    Lê o CSV bruto, renomeia as colunas para o padrão do modelo, remove as
    linhas com precipitação ausente, calcula a fase ENSO oficial (NOAA) a
    partir do ONI, e salva o Excel que o script usa. Retorna True se o CSV foi encontrado e
    processado, False caso contrário.
    """
    caminho_csv = caminho(NOME_CSV_BRUTO)
    if not os.path.exists(caminho_csv):
        return False

    df = pd.read_csv(caminho_csv, sep=';', decimal=',')
    df = df.rename(columns=RENOME_COLUNAS_CSV)

    antes = len(df)
    linhas_removidas = df[df[COLUNAS_PRECIP].isna().any(axis=1)][['ANO', 'MÊS']]
    df = df.dropna(subset=COLUNAS_PRECIP).reset_index(drop=True)
    depois = len(df)

    df['Enso'] = classificar_fase_enso_noaa(df, BASE_DIR)

    print(f"\n===== PREPARAÇÃO DA BASE (a partir do CSV) =====")
    print(f"Arquivo lido: {NOME_CSV_BRUTO}")
    print(f"Linhas originais: {antes}")
    if depois < antes:
        print("Linhas removidas (precipitação ausente — dado ainda não publicado pela NASA):")
        print(linhas_removidas.to_string(index=False))
    print(f"Linhas finais usadas: {depois}")
    print("Distribuição de fases ENSO (classificação oficial NOAA):")
    print(df['Enso'].value_counts().to_string())

    df.to_excel(caminho(NOME_XLSX_MODELO), index=False)
    print(f"Base atualizada salva em {NOME_XLSX_MODELO}\n")
    return True


csv_processado = preparar_base_a_partir_do_csv()
if not csv_processado:
    print(f"\n(CSV '{NOME_CSV_BRUTO}' não encontrado na pasta — "
          f"usando '{NOME_XLSX_MODELO}' existente diretamente.)\n")

# =============================================================================
# CARREGAMENTO E PREPARO
# =============================================================================

dados_total = pd.read_excel(caminho(NOME_XLSX_MODELO))
dados_total.columns = dados_total.columns.str.strip()

# Padroniza a variável resposta como PSN (a base interna usa o prefixo NP_).
dados_total = dados_total.rename(columns={'NP_MA': 'PSN_MA',
                                          'NP_CE': 'PSN_CE',
                                          'NP_CA': 'PSN_CA'})

# Se a base já veio com uma coluna 'Enso' própria (ex.: planilha antiga
# Dados_Benfica_.xlsx usada manualmente), ela é respeitada; senão, a fase é
# calculada pela regra oficial da NOAA (ver enso_noaa.py).
if 'Enso' not in dados_total.columns:
    dados_total['Enso'] = classificar_fase_enso_noaa(dados_total, BASE_DIR)

dados_total['DATA'] = pd.to_datetime(
    dados_total['ANO'].astype(str) + '-' + dados_total['MÊS'].astype(str) + '-01')

biomas = {'Mata Atlântica': 'PSN_MA', 'Cerrado': 'PSN_CE', 'Caatinga': 'PSN_CA'}
sigla_bioma = {'Mata Atlântica': 'MA', 'Cerrado': 'CE', 'Caatinga': 'CA'}

cores_biomas = {'Mata Atlântica': '#185FA5', 'Cerrado': '#3B6D11', 'Caatinga': '#BA7517'}

# Ordem e cores das fases (mesma ordem das figuras: La Niña, Neutro, El Niño)
fases       = ['La Niña', 'Neutro', 'El Niño']
cores_fases = {'La Niña': '#5B9BD5', 'Neutro': '#A6A6A6', 'El Niño': '#ED7D6E'}

def rotulo_sig(p):
    if p < 0.001: return '***'
    if p < 0.01:  return '**'
    if p < 0.05:  return '*'
    return 'ns'

# =============================================================================
# DESCRITIVO: PICOS E ANO DE MAIOR PSN
# =============================================================================

print("\n===== PICO HISTÓRICO MENSAL DE PSN POR BIOMA =====")
for nome, col in biomas.items():
    lin_max = dados_total.loc[dados_total[col].idxmax()]
    lin_min = dados_total.loc[dados_total[col].idxmin()]
    print(f"\n--- {nome} ---")
    print(f"  MAIOR PSN : {lin_max[col]:.1f} | {int(lin_max['MÊS'])}/{int(lin_max['ANO'])} "
          f"| ONI={lin_max['ONI']:.2f} | Fase={lin_max['Enso']}")
    print(f"  MENOR PSN : {lin_min[col]:.1f} | {int(lin_min['MÊS'])}/{int(lin_min['ANO'])} "
          f"| ONI={lin_min['ONI']:.2f} | Fase={lin_min['Enso']}")
    print(f"  Média     : {dados_total[col].mean():.1f} +/- {dados_total[col].std():.1f}")

df_anual = dados_total.groupby('ANO').agg(
    {'PSN_MA': 'mean', 'PSN_CE': 'mean', 'PSN_CA': 'mean', 'ONI': 'mean'}).reset_index()

print("\n===== ANO DE MAIOR PSN MÉDIA ANUAL POR BIOMA =====")
for nome, col in biomas.items():
    lin = df_anual.loc[df_anual[col].idxmax()]
    print(f"{nome}: ANO {int(lin['ANO'])} | PSN={lin[col]:.1f} | ONI médio={lin['ONI']:.2f}")

# =============================================================================
# SÉRIE TEMPORAL DE PSN E ONI
# =============================================================================

fig_st, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 8), sharex=True,
                                   gridspec_kw={'height_ratios': [3, 1]})
for nome, col in biomas.items():
    ax1.plot(dados_total['DATA'], dados_total[col], label=nome,
             color=cores_biomas[nome], linewidth=1.2, alpha=0.85)
    im = dados_total[col].idxmax()
    ax1.scatter(dados_total['DATA'][im], dados_total[col][im],
                color=cores_biomas[nome], zorder=5, s=80, marker='*')
ax1.set_ylabel('PSN (gC/m²/mês)')
ax1.set_title(f'Série Temporal de PSN - 3 Biomas ({int(dados_total["ANO"].min())}-{int(dados_total["ANO"].max())})')
ax1.legend(loc='upper right'); ax1.grid(alpha=0.3)

ax2.fill_between(dados_total['DATA'], dados_total['ONI'],
                 where=(dados_total['ONI'] > 0.5), color='#D85A30', alpha=0.6, label='El Niño')
ax2.fill_between(dados_total['DATA'], dados_total['ONI'],
                 where=(dados_total['ONI'] < -0.5), color='#185FA5', alpha=0.6, label='La Niña')
ax2.plot(dados_total['DATA'], dados_total['ONI'], color='gray', linewidth=0.8)
ax2.axhline(0, color='black', linewidth=0.5, linestyle='--')
ax2.set_ylabel('ONI'); ax2.legend(loc='upper right', fontsize=8); ax2.grid(alpha=0.3)

plt.tight_layout()
plt.savefig(saida('serie_temporal_biomas.png'), dpi=200)
plt.close(fig_st)
print("Figura salva: serie_temporal_biomas.png")

# =============================================================================
# KRUSKAL-WALLIS POR FASE ENSO E CORRELAÇÃO DE PEARSON (ONI x VARIÁVEIS)
# =============================================================================

variaveis = {
    'PSN': {'MA': 'PSN_MA', 'CE': 'PSN_CE', 'CA': 'PSN_CA'},
    'EV':  {'MA': 'EV_MA',  'CE': 'EV_CE',  'CA': 'EV_CA'},
    'PRE': {'MA': 'PRE_MA', 'CE': 'PRE_CE', 'CA': 'PRE_CA'},
    'TST': {'MA': 'TST_MA', 'CE': 'TST_CE', 'CA': 'TST_CA'},
}

print("\n===== KRUSKAL-WALLIS POR FASE ENSO =====")
for var, cols in variaveis.items():
    print(f"\n[ {var} ]")
    for nome, sig in sigla_bioma.items():
        grupos = [dados_total[dados_total['Enso'] == f][cols[sig]].dropna().values for f in fases]
        medias = [g.mean() for g in grupos]
        stat, p = kruskal(*grupos)
        print(f"  {nome}: La Niña={medias[0]:.1f} | Neutro={medias[1]:.1f} | El Niño={medias[2]:.1f} "
              f"| H={stat:.2f} | p={p:.4f} {rotulo_sig(p)}")

print("\n===== CORRELAÇÃO DE PEARSON: ONI x VARIÁVEIS =====")
for var, cols in variaveis.items():
    print(f"\n[ {var} ]")
    for nome, sig in sigla_bioma.items():
        par = dados_total[['ONI', cols[sig]]].dropna()
        r, p = pearsonr(par['ONI'], par[cols[sig]])
        print(f"  {nome}: r={r:.3f} | p={p:.4f} {rotulo_sig(p)} "
              f"({'positiva' if r > 0 else 'negativa'})")

# =============================================================================
# FIGURAS POR BIOMA: DISTRIBUIÇÃO POR FASE ENSO (PSN, EV, PRE, TST)
# Uma figura por bioma, reproduzindo as Figuras 8, 9 e 10 da dissertação.
# =============================================================================

paineis = [
    ('PSN', r'PSN (gC$\cdot$m$^{-2}\cdot$mês$^{-1}$)'),
    ('EV',  r'EV (mm$\cdot$mês$^{-1}$)'),
    ('PRE', r'PRE (mm$\cdot$mês$^{-1}$)'),
    ('TST', r'TST (°C)'),
]

for nome, sig in sigla_bioma.items():
    fig, axes = plt.subplots(1, 4, figsize=(16, 5.2))
    for k, (var, ylab) in enumerate(paineis):
        ax  = axes[k]
        col = variaveis[var][sig]
        grupos = [dados_total[dados_total['Enso'] == f][col].dropna().values for f in fases]

        bp = ax.boxplot(grupos, patch_artist=True,
                        medianprops=dict(color='black', linewidth=2))
        for patch, f in zip(bp['boxes'], fases):
            patch.set_facecolor(cores_fases[f])
            patch.set_alpha(0.85)

        stat, p = kruskal(*grupos)
        ax.text(0.97, 0.97, ('KW: p < 0,001' if p < 0.001 else f'KW: p = {p:.3f}'.replace('.', ',')) + f' {rotulo_sig(p)}',
                transform=ax.transAxes, ha='right', va='top', fontsize=12,
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                          edgecolor='0.7', alpha=0.9))

        ax.set_xticks([1, 2, 3])
        ax.set_xticklabels(fases, fontsize=12)
        ax.tick_params(axis='y', labelsize=11)
        ax.set_title(var, fontsize=15, fontweight='bold')
        ax.set_ylabel(ylab, fontsize=12)
        ax.grid(alpha=0.3, axis='y')

    fig.suptitle(f'Distribuição por Fase ENSO - {nome}', fontsize=17, fontweight='bold')
    plt.tight_layout()
    plt.savefig(saida(f'boxplot_enso_{sig}.png'), dpi=600, bbox_inches='tight')
    plt.close(fig)
    print(f"Figura salva: boxplot_enso_{sig}.png")

# =============================================================================
# REGRESSÃO LINEAR SIMPLES: ONI -> EV, PRE E TST
# Sustenta a afirmação de que o ONI explica menos de 2% da variância das
# variáveis climáticas (R² do OLS). Aqui a estatística F é válida (OLS, não Ridge).
# =============================================================================

dependentes = {
    'EV':  {'MA': 'EV_MA',  'CE': 'EV_CE',  'CA': 'EV_CA'},
    'PRE': {'MA': 'PRE_MA', 'CE': 'PRE_CE', 'CA': 'PRE_CA'},
    'TST': {'MA': 'TST_MA', 'CE': 'TST_CE', 'CA': 'TST_CA'},
}

print("\n===== REGRESSÃO SIMPLES: ONI -> VARIÁVEIS AMBIENTAIS =====")
for var, cols in dependentes.items():
    print(f"\n[ {var} ]")
    for nome, sig in sigla_bioma.items():
        reg = dados_total[['ONI', cols[sig]]].dropna()
        modelo = sm.OLS(reg[cols[sig]], sm.add_constant(reg[['ONI']])).fit()
        print(f"  {nome}: {var} = {modelo.params['const']:.2f} + {modelo.params['ONI']:.4f} x ONI "
              f"| R²={modelo.rsquared:.4f} ({modelo.rsquared*100:.2f}%) "
              f"| F={modelo.fvalue:.2f} | p={modelo.f_pvalue:.4f} {rotulo_sig(modelo.f_pvalue)}")

# =============================================================================
# DISPERSÃO PSN x VARIÁVEIS AMBIENTAIS (estilo Benfica et al., 2022)
# =============================================================================

for suf in ['MA', 'CE', 'CA']:
    dados_total[f'BURN_{suf}_log'] = np.log1p(dados_total[f'BURN_{suf}'])

var_disp = {
    'Precipitação (mm)':      {'MA': 'PRE_MA', 'CE': 'PRE_CE', 'CA': 'PRE_CA'},
    'Evapotranspiração (mm)': {'MA': 'EV_MA',  'CE': 'EV_CE',  'CA': 'EV_CA'},
    'WAI':                    {'MA': 'WAI_MA', 'CE': 'WAI_CE', 'CA': 'WAI_CA'},
    'Temperatura (°C)':       {'MA': 'TST_MA', 'CE': 'TST_CE', 'CA': 'TST_CA'},
    'Área queimada log(ha)':  {'MA': 'BURN_MA_log', 'CE': 'BURN_CE_log', 'CA': 'BURN_CA_log'},
}
psn_cols   = {'MA': 'PSN_MA', 'CE': 'PSN_CE', 'CA': 'PSN_CA'}
nomes_bio  = {'MA': 'Mata Atlântica', 'CE': 'Cerrado', 'CA': 'Caatinga'}
cores_disp = {'MA': '#185FA5', 'CE': '#2E5A0B', 'CA': '#A3620F'}

n_vars = len(var_disp)
fig, axes = plt.subplots(n_vars, 3, figsize=(14, n_vars * 3.5), constrained_layout=True)

print(f"\n===== REGRESSÃO SIMPLES PSN x VARIÁVEIS (n={len(dados_total)}) =====")
print(f"{'Variável':<26} {'Bioma':<18} {'R²adj':>6} {'Slope':>10} {'p-valor':>12} Signif.")
print("-" * 82)

for row, (var, cols) in enumerate(var_disp.items()):
    for c, bio in enumerate(['MA', 'CE', 'CA']):
        ax = axes[row, c]
        x = dados_total[cols[bio]].values
        y = dados_total[psn_cols[bio]].values
        mask = ~(np.isnan(x) | np.isnan(y))
        x, y = x[mask], y[mask]

        slope, interc, r, p, se = stats.linregress(x, y)
        r2_adj = 1 - (1 - r**2) * (len(x) - 1) / (len(x) - 2)

        ax.scatter(x, y, color=cores_disp[bio], alpha=0.55, s=16, edgecolors='white', linewidths=0.3)
        xl = np.linspace(x.min(), x.max(), 200)
        ax.plot(xl, slope * xl + interc, color='#C00000', linewidth=2.3)

        p_str = "p < 0,001" if p < 0.001 else f"p = {p:.3f}".replace('.', ',')
        ax.set_title(f"R²aj = {r2_adj:.2f}   b = {slope:.2f}   {p_str}".replace('.', ','), fontsize=11.5, pad=6)
        if row == 0:
            ax.text(0.5, 1.24, nomes_bio[bio], transform=ax.transAxes, ha='center',
                    fontsize=17, fontweight='bold', color=cores_disp[bio])
        ax.set_xlabel(var, fontsize=11)
        ax.set_ylabel(r'PSN (gC$\cdot$m$^{-2}\cdot$mês$^{-1}$)' if c == 0 else '', fontsize=11)
        ax.tick_params(labelsize=9.5)
        ax.grid(alpha=0.15)

        sig = rotulo_sig(p)
        print(f"{var:<26} {nomes_bio[bio]:<18} {r2_adj:>6.3f} {slope:>10.3f} {p:>12.4e} {sig}")

_meses_abrev = ['Jan', 'Fev', 'Mar', 'Abr', 'Mai', 'Jun', 'Jul', 'Ago', 'Set', 'Out', 'Nov', 'Dez']
_ini = dados_total.sort_values('DATA').iloc[0]; _fim = dados_total.sort_values('DATA').iloc[-1]
fig.suptitle(
    'Regressão Linear Simples: PSN x Variáveis Ambientais - 3 Biomas\n'
    f'Dados mensais (n={len(dados_total)}) | '
    f'{_meses_abrev[int(_ini["MÊS"]) - 1]}/{int(_ini["ANO"])}-{_meses_abrev[int(_fim["MÊS"]) - 1]}/{int(_fim["ANO"])}',
    fontsize=15, fontweight='bold', y=1.04)
plt.savefig(saida('dispersao_psn_variaveis_biomas.png'), dpi=300, bbox_inches='tight')
plt.close(fig)
print("\nFigura salva: dispersao_psn_variaveis_biomas.png")
