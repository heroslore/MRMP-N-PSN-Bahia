# -*- coding: utf-8 -*-
"""
Monta as figuras da dissertação a partir das saídas curadas em
resultados_2001_2025/ (nomes fixos). Gera em figuras_dissertacao/:

  fig01_psn_conceito.png        (nova) GPP -> PSN -> NPP
  fig_validacao_esquema.png     (nova) RepeatedKFold / GroupKFold / TimeSeriesSplit
  fig03_obs_vs_pred.png         3 painéis (MA, CE, CA)
  fig04_selecao_grau.png        3 linhas (MA, CE, CA)
  fig05_importancia.png         importância relativa (gerada dos coeficientes)
  fig06_residuos.png            3 linhas (MA, CE, CA)
  fig07_yrandomization.png      3 painéis
  fig08/09/10_enso_*.png        boxplots por fase ENSO (cópias)
  fig_dispersao.png             PSN x variáveis (cópia)
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
from PIL import Image, ImageDraw, ImageFont

BASE = os.path.dirname(os.path.abspath(__file__))
RES  = os.path.join(BASE, 'resultados_2001_2025')
OUT  = os.path.join(BASE, 'figuras_dissertacao')
os.makedirs(OUT, exist_ok=True)

BIOMAS = [('MA', 'Mata Atlântica', '#2CA02C'), ('CE', 'Cerrado', '#FF7F0E'), ('CA', 'Caatinga', '#D62728')]
LARGURA_FINAL = 4800   # px (~ 30 cm a 400 dpi; o Word reduz para a largura da página)

def _font(px):
    for f in ['/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf',
              matplotlib.font_manager.findfont('DejaVu Sans:bold')]:
        try: return ImageFont.truetype(f, px)
        except Exception: pass
    return ImageFont.load_default()

def compor(arquivos, saida, orientacao='h', rotulos=None, largura=LARGURA_FINAL, fundo='white', pad=40):
    """Junta PNGs lado a lado ('h') ou empilhados ('v'), com letras (a), (b), (c)."""
    ims = [Image.open(a).convert('RGB') for a in arquivos]
    if orientacao == 'h':
        h = min(im.height for im in ims)
        ims = [im.resize((int(im.width * h / im.height), h), Image.LANCZOS) for im in ims]
        W = sum(im.width for im in ims) + pad * (len(ims) - 1); H = h
        tela = Image.new('RGB', (W, H), fundo); x = 0
        for im in ims:
            tela.paste(im, (x, 0)); x += im.width + pad
        pos = []; x = 0
        for im in ims: pos.append((x, 0)); x += im.width + pad
    else:
        w = min(im.width for im in ims)
        ims = [im.resize((w, int(im.height * w / im.width)), Image.LANCZOS) for im in ims]
        W = w; H = sum(im.height for im in ims) + pad * (len(ims) - 1)
        tela = Image.new('RGB', (W, H), fundo); y = 0; pos = []
        for im in ims:
            tela.paste(im, (0, y)); pos.append((0, y)); y += im.height + pad
    if rotulos:
        d = ImageDraw.Draw(tela); f = _font(int(H * 0.035) if orientacao == 'h' else int(W * 0.018))
        for (x, y), r in zip(pos, rotulos):
            d.text((x + 10, y + 5), r, fill='black', font=f)
    if tela.width > largura:
        tela = tela.resize((largura, int(tela.height * largura / tela.width)), Image.LANCZOS)
    tela.save(saida, dpi=(300, 300), optimize=True)
    print('salvo', os.path.basename(saida), tela.size)

# ---------------------------------------------------------------- Fig 3, 4, 6, 7, 8-10
compor([os.path.join(RES, b, f'yrandomization_{b}.png') for b, _, _ in BIOMAS],
       os.path.join(OUT, 'fig07_yrandomization.png'), 'h', ['(a)', '(b)', '(c)'])
for n, (b, _, _) in zip([8, 9, 10], BIOMAS):
    im = Image.open(os.path.join(RES, f'boxplot_enso_{b}.png')).convert('RGB')
    im = im.resize((LARGURA_FINAL, int(im.height * LARGURA_FINAL / im.width)), Image.LANCZOS)
    im.save(os.path.join(OUT, f'fig{n:02d}_enso_{b}.png'), dpi=(300, 300), optimize=True)
    print('salvo', f'fig{n:02d}_enso_{b}.png', im.size)
Image.open(os.path.join(RES, 'dispersao_psn_variaveis_biomas.png')).convert('RGB').save(
    os.path.join(OUT, 'fig_dispersao.png'), dpi=(200, 200), optimize=True)
print('salvo fig_dispersao.png')

# ---------------------------------------------------------------- Fig 5: importância relativa
plt.rcParams.update({'font.size': 13, 'font.family': 'DejaVu Sans'})
NOMES = {'EV': 'EV', 'PRE': 'PRE', 'TST': 'TST', 'WAI': 'WAI', 'saz_sin': 'SAZsin', 'saz_cos': 'SAZcos'}
fig, axes = plt.subplots(1, 3, figsize=(16, 5.4))
letras = ['(a)', '(b)', '(c)']
for ax, (b, nome, cor), letra in zip(axes, BIOMAS, letras):
    c = pd.read_csv(os.path.join(RES, b, f'coeficientes_ridge_{b}.csv'))
    imp = {}
    for _, r in c.iterrows():
        termo = r['feature'].replace('^2', '')
        vars_ = set(v.replace(f'_{b}', '') for v in termo.split(' '))
        for v in vars_:
            imp[v] = imp.get(v, 0.0) + abs(r['coef'])
    s = pd.Series(imp); s = (s / s.sum() * 100).sort_values()
    ax.barh([NOMES.get(k, k) for k in s.index], s.values, color=cor, height=0.62)
    for i, v in enumerate(s.values):
        ax.text(v + 0.6, i, f'{v:.1f}%', va='center', fontsize=12.5, fontweight='bold')
    ax.set_title(f'{letra} {nome}', fontweight='bold', fontsize=15)
    ax.tick_params(labelsize=12)
    ax.set_xlabel('Importância relativa (%)', fontweight='bold')
    ax.set_xlim(0, s.max() * 1.22)
    ax.grid(axis='x', alpha=0.3); ax.spines[['top', 'right']].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(OUT, 'fig05_importancia.png'), dpi=300, bbox_inches='tight')
plt.close(fig); print('salvo fig05_importancia.png')

# ---------------------------------------------------------------- Fig 4: seleção do grau (a partir do CSV)
sg = pd.read_csv(os.path.join(RES, 'selecao_grau.csv'))
plt.rcParams.update({'font.size': 12, 'font.family': 'DejaVu Sans'})
fig, axes = plt.subplots(2, 3, figsize=(16, 10.5), gridspec_kw={'height_ratios': [3, 2.3]})
for j, (b, nome, cor) in enumerate(BIOMAS):
    d = sg[sg['bioma'] == b].sort_values('grau')
    ax = axes[0, j]
    piso = np.floor(min(d['r2_teste'].iloc[:3].min() - d['dp_teste'].iloc[:3].max(), d['r2_treino'].iloc[0]) / 10) * 10 - 10
    piso = max(piso, 0)
    ax.plot(d['grau'], d['r2_treino'], '-o', color='#555555', lw=2, ms=8, label='R² treino')
    ax.errorbar(d['grau'], d['r2_teste'], yerr=d['dp_teste'], fmt='-s', color=cor, lw=2.2, ms=9,
                capsize=4, label='R² teste (± 1 dp)')
    for g_, r2_, dp_ in zip(d['grau'], d['r2_teste'], d['dp_teste']):
        if r2_ >= piso:
            ax.annotate(f'{r2_:.1f}%'.replace('.', ','), (g_, r2_), xytext=(0, -16), textcoords='offset points',
                        ha='center', va='top', fontsize=11, fontweight='bold', color=cor,
                        bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='none', alpha=0.85))
        else:
            ax.annotate(f'{r2_:.0f}%', (g_, piso + (100 - piso) * 0.42), ha='center', va='bottom', fontsize=12,
                        color=cor, fontweight='bold')
            ax.plot([g_], [piso + 0.8], marker='v', color=cor, ms=9, clip_on=False)
    ax.axvline(2, color='#999999', ls=':', lw=1.5, zorder=0)
    ax.set_ylim(piso, 100); ax.set_xticks([1, 2, 3, 4, 5]); ax.set_xlim(0.6, 5.4)
    ax.set_title(f'({"abc"[j]}) {nome}', fontweight='bold', fontsize=15)
    ax.set_ylabel('R² (%)' if j == 0 else ''); ax.grid(alpha=0.3)
    if j == 0: ax.legend(loc='lower right', fontsize=11, framealpha=0.95)
    ax = axes[1, j]
    cores_b = ['#2CA02C' if g_ == 2 else ('#FF7F0E' if v_ <= 10 else '#9E9E9E') for g_, v_ in zip(d['grau'], d['gap_pp'])]
    ax.bar(d['grau'], d['gap_pp'], color=cores_b, width=0.6)
    ax.axhline(10, color='#D62728', ls='--', lw=1.3)
    ax.set_yscale('log'); ax.set_ylim(0.4, 1000)
    ax.set_yticks([1, 10, 100]); ax.set_yticklabels(['1', '10', '100'])
    for g_, v_, t_ in zip(d['grau'], d['gap_pp'], d['termos']):
        ax.text(g_, v_ * 1.18, f'{v_:.1f}'.replace('.', ','), ha='center', va='bottom', fontsize=11.5,
                fontweight='bold' if g_ == 2 else 'normal')
        ax.text(g_, v_ * 2.3, f'({int(t_)} termos)', ha='center', va='bottom', fontsize=9.5, style='italic', color='#333')
    ax.set_xticks([1, 2, 3, 4, 5]); ax.set_xlim(0.6, 5.4)
    ax.set_xlabel('Grau polinomial', fontsize=13)
    ax.set_ylabel('Diferença treino − teste (pp, escala log)' if j == 0 else '', fontsize=12)
    ax.text(5.35, 520, 'referência 10 pp', ha='right', va='bottom', fontsize=10.5, color='#D62728')
    ax.grid(axis='y', alpha=0.3)
fig.suptitle('Seleção do grau polinomial (1 a 5) — com nº de termos e R² de teste anotados', fontsize=15, fontweight='bold', y=0.995)
fig.text(0.5, 0.005, 'Grau 2 (verde) escolhido como grau comum: no Cerrado e na Caatinga já supera o grau 3 em R² de teste; na Mata Atlântica,\n'
         'o grau 3 ganha +0,8 pp às custas de mais que o dobro de termos (55 vs. 20) e maior diferença treino−teste.',
         ha='center', va='bottom', fontsize=11.5, style='italic', color='#333')
plt.tight_layout(rect=(0, 0.045, 1, 0.97))
plt.savefig(os.path.join(OUT, 'fig04_selecao_grau.png'), dpi=300, bbox_inches='tight'); plt.close(fig)
print('salvo fig04_selecao_grau.png (gerada)')

# ---------------------------------------------------------------- Fig nova: conceito da PSN
fig, ax = plt.subplots(figsize=(15, 5.8)); ax.set_xlim(0, 15); ax.set_ylim(0, 5.8); ax.axis('off')
def caixa(x, y, w, h, titulo, sub, cor, tcor='white'):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.02,rounding_size=0.15', fc=cor, ec='none'))
    ax.text(x + w / 2, y + h * 0.70, titulo, ha='center', va='center', fontsize=17, fontweight='bold', color=tcor)
    ax.text(x + w / 2, y + h * 0.30, sub, ha='center', va='center', fontsize=10.5, color=tcor, linespacing=1.3)
def seta(x0, x1, y, texto):
    ax.add_patch(FancyArrowPatch((x0, y), (x1, y), arrowstyle='-|>', mutation_scale=24, lw=2.4, color='#333333'))
    ax.text((x0 + x1) / 2, y + 0.22, texto, ha='center', va='bottom', fontsize=10, color='#8B0000', linespacing=1.25)
caixa(0.2, 1.7, 2.9, 1.9, 'GPP', 'Produtividade Primária\nBruta: todo o carbono\nfixado pela fotossíntese', '#1B5E20')
seta(3.2, 5.95, 2.65, '− respiração de\nmanutenção de folhas\ne raízes finas')
caixa(6.05, 1.7, 2.9, 1.9, 'PSN', 'Fotossíntese Líquida:\nsaldo de carbono em\n8 dias (MOD17A2H)', '#2E7D32')
seta(9.05, 11.8, 2.65, '− respiração de\nmanutenção do lenho\ne raízes grossas\n− respiração de\ncrescimento')
caixa(11.9, 1.7, 2.9, 1.9, 'NPP', 'Produtividade Primária\nLíquida: agregação\nanual (MOD17A3)', '#558B2F')
ax.text(0.2, 1.05, 'Atmosfera (CO₂)  →  fotossíntese  →  GPP  →  PSN  →  NPP: cada etapa desconta uma parcela de respiração.',
        fontsize=11.5, color='#333')
ax.text(0.2, 0.5, 'Variável-resposta deste estudo: PSN mensal, por bioma (compostos de 8 dias do MOD17A2H agregados ao mês).',
        fontsize=11.5, color='#333')
ax.text(0.2, 5.2, 'Fluxos de carbono estimados pelo algoritmo MODIS/MOD17', fontsize=14, fontweight='bold')
plt.savefig(os.path.join(OUT, 'fig01_psn_conceito.png'), dpi=300, bbox_inches='tight'); plt.close(fig)
print('salvo fig01_psn_conceito.png')

# ---------------------------------------------------------------- Fig nova: esquema de validação
n = 297
fig, axes = plt.subplots(3, 1, figsize=(13, 8.6), gridspec_kw={'height_ratios': [5.2, 4, 4]})
TREINO, TESTE = '#BDBDBD', '#D85A30'
ax = axes[0]
rng = np.random.RandomState(42)
for rep in range(3):
    perm = rng.permutation(n); folds = np.array_split(perm, 5)
    for k in range(5):
        y = rep * 6 + k
        ax.broken_barh([(0, n)], (y, 0.8), color=TREINO)
        idx = np.sort(folds[k]); blocos = []; ini = idx[0]; prev = idx[0]
        for i in idx[1:]:
            if i != prev + 1: blocos.append((ini, prev - ini + 1)); ini = i
            prev = i
        blocos.append((ini, prev - ini + 1))
        ax.broken_barh(blocos, (y, 0.8), color=TESTE)
        ax.text(-3, y + 0.4, f'rep. {rep + 1}, fold {k + 1}', ha='right', va='center', fontsize=8)
ax.text(n / 2, 18.3, '⋯ e assim até a repetição 30, totalizando 150 partições', ha='center', va='center',
        fontsize=9.5, style='italic')
ax.set_ylim(-0.5, 19.2); ax.invert_yaxis(); ax.set_xlim(0, n)
ax.set_yticks([]); ax.set_xlabel('Índice da observação mensal (jan/2001 → set/2025), n = 297', fontsize=9)
ax.set_title('(a) RepeatedKFold, 5 folds × 30 repetições: em cada partição, 80% treino (cinza) e 20% teste (laranja),\n'
             'sorteados sem respeitar a ordem temporal; cada mês é teste 30 vezes', fontsize=10.5, loc='left')
ax.spines[['top', 'right', 'left']].set_visible(False)
anos_lista = np.arange(2001, 2026); grupos = np.array_split(anos_lista, 5)
meses_por_ano = {a: 12 for a in anos_lista}; meses_por_ano[2025] = 9
ini_ano = {}; acc = 0
for a in anos_lista: ini_ano[a] = acc; acc += meses_por_ano[a]
ax = axes[1]
for k, g in enumerate(grupos):
    ax.broken_barh([(0, n)], (k, 0.8), color=TREINO)
    ax.broken_barh([(ini_ano[a], meses_por_ano[a]) for a in g], (k, 0.8), color=TESTE)
    ax.text(-3, k + 0.4, f'fold {k + 1}', ha='right', va='center', fontsize=8)
    ax.text(ini_ano[g[0]] + sum(meses_por_ano[a] for a in g) / 2, k + 0.4, f'{g[0]}–{g[-1]}',
            ha='center', va='center', fontsize=8, color='white', fontweight='bold')
ax.set_ylim(-0.5, 5.3); ax.invert_yaxis(); ax.set_xlim(0, n); ax.set_yticks([])
ax.set_xticks([ini_ano[a] for a in anos_lista[::4]]); ax.set_xticklabels(anos_lista[::4], fontsize=8)
ax.set_title('(b) GroupKFold por ano: anos inteiros ficam fora do treino (5 grupos de 5 anos)', fontsize=10.5, loc='left')
ax.spines[['top', 'right', 'left']].set_visible(False)
ax = axes[2]
from sklearn.model_selection import TimeSeriesSplit
for k, (tr, te) in enumerate(TimeSeriesSplit(n_splits=5).split(np.arange(n))):
    ax.broken_barh([(tr[0], len(tr))], (k, 0.8), color=TREINO)
    ax.broken_barh([(te[0], len(te))], (k, 0.8), color=TESTE)
    ax.text(-3, k + 0.4, f'fold {k + 1}', ha='right', va='center', fontsize=8)
    ax.text(te[0] + len(te) / 2, k + 0.4, f'teste: {len(te)} meses', ha='center', va='center', fontsize=8, color='white', fontweight='bold')
    ax.text(tr[0] + len(tr) / 2, k + 0.4, f'treino: {len(tr)} meses', ha='center', va='center', fontsize=8, color='#333')
ax.set_ylim(-0.5, 5.3); ax.invert_yaxis(); ax.set_xlim(0, n); ax.set_yticks([])
ax.set_xticks([ini_ano[a] for a in anos_lista[::4]]); ax.set_xticklabels(anos_lista[::4], fontsize=8)
ax.set_title('(c) TimeSeriesSplit com janela expansível: treina no passado (cinza) e testa no bloco seguinte (laranja)', fontsize=10.5, loc='left')
ax.spines[['top', 'right', 'left']].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(OUT, 'fig_validacao_esquema.png'), dpi=300, bbox_inches='tight'); plt.close(fig)
print('salvo fig_validacao_esquema.png')

# ---------------------------------------------------------------- Fig: fluxo metodológico (refeito, base 2001-2025)
fig, ax = plt.subplots(figsize=(14, 8.6)); ax.set_xlim(0, 14); ax.set_ylim(0, 8.6); ax.axis('off')
def bloco(x, y, w, h, num, titulo, linhas, cor_t, cor_c):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.02,rounding_size=0.12', fc=cor_c, ec='#555', lw=0.8))
    ax.add_patch(FancyBboxPatch((x, y + h - 0.62), w, 0.62, boxstyle='round,pad=0.02,rounding_size=0.12', fc=cor_t, ec='none'))
    ax.add_patch(plt.Circle((x + 0.32, y + h - 0.31), 0.19, fc='white', ec='none'))
    ax.text(x + 0.32, y + h - 0.31, str(num), ha='center', va='center', fontsize=9.5, fontweight='bold', color=cor_t)
    ax.text(x + w / 2 + 0.15, y + h - 0.31, titulo, ha='center', va='center', fontsize=11, fontweight='bold', color='white')
    ax.text(x + w / 2, y + (h - 0.62) / 2, linhas, ha='center', va='center', fontsize=8.6, color='#222', linespacing=1.35)
def seta_h(x0, x1, y): ax.add_patch(FancyArrowPatch((x0, y), (x1, y), arrowstyle='-|>', mutation_scale=16, lw=1.6, color='#333'))
W, H, G = 3.05, 1.55, 0.35; xs = [0.3 + i * (W + G) for i in range(4)]
y1, y2, y3 = 6.6, 4.0, 1.3
azul, azul_c = '#1F3B5C', '#D6E2F0'; verde, verde_c = '#1B5E3A', '#D5EADF'; marrom, marrom_c = '#5C4632', '#E8DED2'; roxo, roxo_c = '#3B2A5C', '#DED6EC'
bloco(xs[0], y1, W, H, 1, 'Base de dados', 'MODIS (PSN, EV, TST, BURN)\nIMERG (PRE) · ONI (NOAA)\n2001–2025, n = 297 meses', azul, azul_c)
bloco(xs[1], y1, W, H, 2, 'Pré-processamento', 'Agregação mensal por bioma\nlog(1 + área queimada)\nexclusão da PET (r = 0,73–0,80)', azul, azul_c)
bloco(xs[2], y1, W, H, 3, 'Eng. de variáveis', 'Sazonalidade harmônica\n(SAZsin, SAZcos) fixa\nem todos os modelos', verde, verde_c)
bloco(xs[3], y1, W, H, 4, 'Seleção de variáveis', '10 combinações C(5,3)\nde EV, PRE, TST, WAI, BURN\nmaior R² de teste', verde, verde_c)
for i in range(3): seta_h(xs[i] + W, xs[i + 1], y1 + H / 2)
ax.add_patch(FancyBboxPatch((0.12, y2 - 0.35), 13.76, H + 1.05, boxstyle='round,pad=0.02,rounding_size=0.15', fc='#FBF4E8', ec='#B08A50', lw=1.3, ls='--'))
ax.text(7, y2 + H + 0.42, 'Executado dentro de cada partição da validação cruzada (sem vazamento de informação)', ha='center', va='center', fontsize=10, style='italic', color='#7A5A20', fontweight='bold')
bloco(xs[0], y2, W, H, 5, 'Winsorização', 'Percentil 3 da PSN\ncalculado só no treino', marrom, marrom_c)
bloco(xs[1], y2, W, H, 6, 'Padronização', 'StandardScaler\nμ = 0, σ = 1 (ajustado no treino)', marrom, marrom_c)
bloco(xs[2], y2, W, H, 7, 'Expansão polinomial', 'Grau 2 (20 termos),\nescolhido entre 1 e 5', marrom, marrom_c)
bloco(xs[3], y2, W, H, 8, 'Regressão Ridge', 'Penalização L2\nα por GridSearchCV (5 folds)', marrom, marrom_c)
for i in range(3): seta_h(xs[i] + W, xs[i + 1], y2 + H / 2)
bloco(xs[0], y3, W, H, 9, 'Validação cruzada', 'RepeatedKFold 5 × 30\nR², RMSE, MAE (média ± dp)', roxo, roxo_c)
bloco(xs[1], y3, W, H, 10, 'Validação temporal', 'GroupKFold por ano\nTimeSeriesSplit expansível', roxo, roxo_c)
bloco(xs[2], y3, W, H, 11, 'Diagnóstico', 'VIF · Shapiro-Wilk\nLjung-Box · ACF', roxo, roxo_c)
bloco(xs[3], y3, W, H, 12, 'Y-randomization', '100 permutações da PSN\nΔ (pp) e p-valor empírico', roxo, roxo_c)
for i in range(3): seta_h(xs[i] + W, xs[i + 1], y3 + H / 2)
# setas verticais entre linhas (última caixa da linha -> primeira da seguinte)
for ya, yb in [(y1, y2), (y2, y3)]:
    xm = xs[3] + W / 2; ym = (ya + yb + H) / 2
    ax.plot([xm, xm, xs[0] + W / 2, xs[0] + W / 2], [ya, ym + 0.02, ym + 0.02, yb + H + 0.12], color='#333', lw=1.6)
    ax.add_patch(FancyArrowPatch((xs[0] + W / 2, yb + H + 0.3), (xs[0] + W / 2, yb + H + 0.02), arrowstyle='-|>', mutation_scale=16, lw=1.6, color='#333'))
plt.savefig(os.path.join(OUT, 'fig_fluxo.png'), dpi=300, bbox_inches='tight'); plt.close(fig)
print('salvo fig_fluxo.png')

# ---------------------------------------------------------------- Fig 6 (obs vs pred, OOF) e Fig 9 (resíduos) gerados dos dados
import json
from scipy import stats
from sklearn.preprocessing import StandardScaler, PolynomialFeatures
from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold, GridSearchCV
from sklearn.metrics import r2_score
from statsmodels.stats.diagnostic import acorr_ljungbox
from statsmodels.tsa.stattools import acf
NUMX = json.load(open(os.path.join(RES, 'numeros_extra.json'), encoding='utf-8'))
df = pd.read_excel(os.path.join(BASE, 'Dados_base_nova_2001_2025.xlsx')); df.columns = df.columns.str.strip()
df['saz_sin'] = np.sin(2 * np.pi * df['MÊS'] / 12); df['saz_cos'] = np.cos(2 * np.pi * df['MÊS'] / 12)
plt.rcParams.update({'font.size': 13, 'axes.titlesize': 15, 'axes.labelsize': 13, 'xtick.labelsize': 11.5,
                     'ytick.labelsize': 11.5, 'legend.fontsize': 11.5})
UN = r'gC$\cdot$m$^{-2}\cdot$mês$^{-1}$'
fig6, ax6 = plt.subplots(1, 3, figsize=(17, 5.8))
fig9, ax9 = plt.subplots(3, 3, figsize=(17, 15))
for k, (b, nome, cor) in enumerate(BIOMAS):
    xcols = [f'{v_}_{b}' for v_ in NUMX[b]['x']] + ['saz_sin', 'saz_cos']
    X = df[xcols]; y = df[f'NP_{b}']
    # OOF (KFold 5, shuffle, rs=42), winsorização 3% só no treino, GridSearch alpha
    pred = np.zeros(len(y))
    for tr, te in KFold(n_splits=5, shuffle=True, random_state=42).split(X):
        ytr = y.iloc[tr].clip(lower=np.percentile(y.iloc[tr], 3))
        sc = StandardScaler(); pf = PolynomialFeatures(degree=2, include_bias=False)
        Xtr = pf.fit_transform(sc.fit_transform(X.iloc[tr])); Xte = pf.transform(sc.transform(X.iloc[te]))
        gs = GridSearchCV(Ridge(), {'alpha': [0.1, 1.0, 10.0, 50.0, 100.0]}, cv=5, scoring='r2').fit(Xtr, ytr)
        pred[te] = gs.best_estimator_.predict(Xte)
    r2oof = r2_score(y, pred) * 100
    ax = ax6[k]
    lim = [min(y.min(), pred.min()) - 5, max(y.max(), pred.max()) + 5]
    ax.plot(lim, lim, 'k--', lw=1.3, label='Predição perfeita (1:1)')
    ax.scatter(y, pred, s=26, color=cor, alpha=0.6, edgecolors='none')
    ax.set_xlim(lim); ax.set_ylim(lim); ax.set_aspect('equal', adjustable='box')
    ax.set_title(f'({"abc"[k]}) {nome}', fontweight='bold')
    ax.set_xlabel(f'PSN observada ({UN})', fontweight='bold'); ax.set_ylabel(f'PSN predita ({UN})' if k == 0 else '', fontweight='bold')
    ax.text(0.04, 0.95, f'R² fora da amostra = {r2oof:.1f}%'.replace('.', ','), transform=ax.transAxes, va='top',
            fontsize=13, fontweight='bold', bbox=dict(boxstyle='round,pad=0.35', fc='white', ec=cor, lw=2))
    ax.legend(loc='lower right'); ax.grid(alpha=0.3); ax.spines[['top', 'right']].set_visible(False)
    # diagnóstico: ajuste completo, y não winsorizado, alpha médio
    sc = StandardScaler(); pf = PolynomialFeatures(degree=2, include_bias=False)
    Xf = pf.fit_transform(sc.fit_transform(X)); m = Ridge(alpha=NUMX[b]['alpha_medio']).fit(Xf, y)
    fitted = m.predict(Xf); res = y.values - fitted
    W, pW = stats.shapiro(res); lb = acorr_ljungbox(res, lags=[12], return_df=True); ac = acf(res, nlags=3)
    a0, a1, a2 = ax9[k]
    a0.scatter(fitted, res, s=22, color=cor, alpha=0.6, edgecolors='none'); a0.axhline(0, color='k', ls='--', lw=1.2)
    a0.set_title('Resíduos vs valores ajustados', fontweight='bold'); a0.set_xlabel(f'Valores ajustados ({UN})'); a0.set_ylabel(f'Resíduos ({UN})')
    (osm, osr), (slope, inter, _) = stats.probplot(res, dist='norm')
    a1.scatter(osm, osr, s=22, color=cor, alpha=0.7, edgecolors='none'); a1.plot(osm, slope * np.array(osm) + inter, 'k-', lw=1.6)
    a1.set_title('QQ-plot dos resíduos', fontweight='bold'); a1.set_xlabel('Quantis teóricos'); a1.set_ylabel('Valores ordenados')
    a2.hist(res, bins=28, color=cor, alpha=0.85, edgecolor='white')
    a2.set_title('Distribuição dos resíduos', fontweight='bold'); a2.set_xlabel(f'Resíduo ({UN})'); a2.set_ylabel('Frequência')
    ptxt = lambda p: 'p < 0,001' if p < 0.001 else f'p = {p:.3f}'.replace('.', ',')
    a2.text(0.03, 0.96, f'Shapiro-Wilk: W = {W:.3f}, {ptxt(pW)}\nLjung-Box(12): {ptxt(float(lb["lb_pvalue"].iloc[0]))}\nACF(1) = {ac[1]:+.2f}'.replace('.', ',').replace('W = 0,', 'W = 0,'),
            transform=a2.transAxes, va='top', fontsize=12, bbox=dict(boxstyle='round,pad=0.35', fc='white', ec=cor, lw=2))
    for a in (a0, a1, a2): a.grid(alpha=0.3); a.spines[['top', 'right']].set_visible(False)
    a0.text(-0.32, 0.5, nome, transform=a0.transAxes, rotation=90, va='center', ha='center', fontsize=17, fontweight='bold', color=cor)
    print(f'{b}: R² OOF {r2oof:.2f} | Shapiro p {pW:.4f} | LB p {float(lb["lb_pvalue"].iloc[0]):.4f} | ACF1 {ac[1]:.3f}')
fig6.tight_layout(); fig6.savefig(os.path.join(OUT, 'fig03_obs_vs_pred.png'), dpi=300, bbox_inches='tight'); plt.close(fig6)
fig9.tight_layout(); fig9.savefig(os.path.join(OUT, 'fig06_residuos.png'), dpi=300, bbox_inches='tight'); plt.close(fig9)
print('salvo fig03_obs_vs_pred.png e fig06_residuos.png (gerados)')
