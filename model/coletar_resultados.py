# -*- coding: utf-8 -*-
"""
Coleta os resultados de uma rodada COMPLETA e SEQUENCIAL do Modelo_PSN.py
(BIOMA_ATIVO=TODOS RODAR_EM_PARALELO=0 TESTAR_GRAUS=1 RODAR_YRANDOMIZATION=1, log salvo em
resultados_2001_2025/log_modelo_psn_completo.txt) e do Analise_Biomas_PSN.py (log_analise_biomas.txt),
e grava em resultados_2001_2025/: resumo_geral.csv, vif_por_bioma.csv, selecao_grau.csv,
numeros_extra.json, enso_resumo.json e as figuras/coeficientes por bioma.

Uso: python3 coletar_resultados.py [--log LOG_MODELO] [--log-enso LOG_ENSO] [--saida PASTA]
(por padrão lê/grava em resultados_2001_2025/; use --saida para conferir sem sobrescrever).
"""
import os, re, json, shutil, argparse
import pandas as pd, numpy as np
BASE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(BASE, 'resultados_2001_2025'); SF = os.path.join(BASE, 'saidas_figuras')
ap = argparse.ArgumentParser(); ap.add_argument('--log', default=os.path.join(RES, 'log_modelo_psn_completo.txt'))
ap.add_argument('--log-enso', default=os.path.join(RES, 'log_analise_biomas.txt')); ap.add_argument('--saida', default=RES); a = ap.parse_args()
OUT = a.saida; os.makedirs(OUT, exist_ok=True)
NOME = {'MA': 'Mata Atlântica', 'CE': 'Cerrado', 'CA': 'Caatinga'}; SIG = {v: k for k, v in NOME.items()}

def num(pat, txt, grp=1, cast=float):
    m = re.search(pat, txt); return cast(m.group(grp)) if m else None

log = open(a.log, encoding='utf-8').read(); log = re.sub(r'Linhas de (Treino|Teste) \[[^\]]*\]', '', log)
# ---- seções por bioma (rodada sequencial)
starts = [(m.start(), m.group(1)) for m in re.finditer(r'===== RODANDO MODELO RIDGE - (MATA ATLÂNTICA|CERRADO|CAATINGA) =====', log)]
sec = {}
for i, (pos, nm) in enumerate(starts):
    fim = starts[i + 1][0] if i + 1 < len(starts) else len(log)
    b = {'MATA ATLÂNTICA': 'MA', 'CERRADO': 'CE', 'CAATINGA': 'CA'}[nm]; sec[b] = log[pos:fim]
assert len(sec) == 3, f'log com {len(sec)} biomas; rode com BIOMA_ATIVO=TODOS e RODAR_EM_PARALELO=0'

rows, vifs, graus, extra = [], [], [], {}
for b in ['MA', 'CE', 'CA']:
    t = sec[b]
    r = dict(bioma=b, nome=NOME[b],
        r2_treino_medio=num(r'R2 Treino médio: ([\d.]+)', t), r2_treino_dp=num(r'R2 Treino desvio padrão: ([\d.]+)', t),
        r2aj_treino_medio=num(r'R2 Ajustado Treino médio: ([\d.]+)', t), r2_teste_medio=num(r'R2 Teste médio: ([\d.]+)', t),
        r2_teste_dp=num(r'R2 Teste desvio padrão: ([\d.]+)', t),
        ic95_lower=num(r'Intervalo empírico 95% R² Teste: \[([\d.]+)% - ([\d.]+)%\]', t, 1), ic95_upper=num(r'Intervalo empírico 95% R² Teste: \[([\d.]+)% - ([\d.]+)%\]', t, 2),
        rmse_teste_medio=num(r'RMSE Teste médio: ([\d.]+)', t), mae_teste_medio=num(r'MAE Teste médio: ([\d.]+)', t),
        r2_groupkfold_ano=num(r'GroupKFold por ano\s*: ([\d.]+)%', t), r2_timeseriessplit=num(r'TimeSeriesSplit expansível: ([\d.]+)%', t),
        queda_groupkfold_pp=num(r'Queda RKF:GroupKFold\s*: ([-+\d.]+) pp', t), queda_timeseries_pp=num(r'Queda RKF:TimeSeries\s*: ([-+\d.]+) pp', t),
        shapiro_p=num(r'Shapiro-Wilk: stat=[\d.]+, p=([\d.]+)', t), ljungbox_p=num(r'Ljung-Box\(12\): p=([\d.]+)', t), acf_lag1=num(r'ACF lag1=([-+\d.]+)', t),
        yrand_diferenca_pp=num(r'Diferença \(Δ CV\)\s*: ([\d.]+) pp', t), yrand_p_empirico=num(r'p-valor empírico unilateral\s*: ([\d.]+)', t))
    r['gap_overfitting_pp'] = round(r['r2_treino_medio'] - r['r2_teste_medio'], 2)
    r['residuos_normais'] = 'Sim' if r['shapiro_p'] > 0.05 else 'Não'
    r['yrand_status'] = 'APROVADO' if (r['yrand_diferenca_pp'] or 0) >= 60 else 'VERIFICAR'
    r['arquivo_coef'] = f'coeficientes_ridge_{b}.csv'
    blk = re.search(r'===== VIF =====\n(.*?)\nAlpha', t, flags=re.S).group(1)
    vv = re.findall(r'\d+\s+(\S+)\s+([\d.]+)', blk); r['vif_max'] = max(float(x) for _, x in vv)
    vifs += [dict(bioma=b, variavel=v_, VIF=round(float(x), 3)) for v_, x in vv]
    rows.append(r)
    tab = re.search(r' Grau \| Termos.*?\n=+\n(.*?)\n=+\n', t, flags=re.S)
    if tab:
        for ln in tab.group(1).split('\n'):
            m = re.match(r'\s*(\d)\s*\|\s*(\d+)\s*\|\s*([-\d.]+)%\s*\|\s*([-\d.]+)%\s*\|\s*([-\d.]+)%\s*\|\s*([-\d.]+)p\.p\s*\|\s*([-\d.]+)\s*\|\s*\[([-\d.]+)% - ([-\d.]+)%\]', ln)
            if m: graus.append([b] + [float(x) for x in m.groups()])
    yr = re.search(r'modelo original\s*: ([-\d.]+)%\nR² CV médio - modelos permutados\s*: ([-\d.]+)% ± ([\d.]+)', t)
    extra[b] = dict(x=None, mae_dp=num(r'MAE Teste desvio padrão: ([\d.]+)', t), shapiro_w=num(r'Shapiro-Wilk: stat=([\d.]+)', t),
                    yr_orig=float(yr.group(1)) if yr else None, yr_perm_media=float(yr.group(2)) if yr else None, yr_perm_dp=float(yr.group(3)) if yr else None,
                    oof_r2=num(r'R² out-of-fold \(Figura 3\): ([\d.]+)%', t), alpha_medio=num(r'Alpha médio utilizado na análise de resíduos: ([\d.]+)', t))
    # figuras e coeficientes da rodada (nomes fixos)
    os.makedirs(os.path.join(OUT, b), exist_ok=True)
    for f in [f'coeficientes_ridge_{b}.csv', f'analise_residuos_{b}.png', f'obs_vs_pred_{b}.png', f'yrandomization_{b}.png', f'selecao_grau_polinomial_{b}.png']:
        src = os.path.join(SF, b, f)
        if os.path.exists(src): shutil.copy(src, os.path.join(OUT, b, f))
    coef = pd.read_csv(os.path.join(OUT, b, f'coeficientes_ridge_{b}.csv'))
    extra[b]['x'] = [f.replace(f'_{b}', '') for f in coef['feature'].iloc[:3]]
pd.DataFrame(rows).to_csv(os.path.join(OUT, 'resumo_geral.csv'), index=False)
pd.DataFrame(vifs).to_csv(os.path.join(OUT, 'vif_por_bioma.csv'), index=False)
if graus:
    pd.DataFrame(graus, columns=['bioma', 'grau', 'termos', 'r2_treino', 'r2_teste', 'r2_teste_dp', 'gap_pp', 'rmse', 'ic_lo', 'ic_hi']).to_csv(os.path.join(OUT, 'selecao_grau.csv'), index=False)
imp = {}
for b in ['MA', 'CE', 'CA']:
    c = pd.read_csv(os.path.join(OUT, b, f'coeficientes_ridge_{b}.csv')); acc = {}
    for _, rr in c.iterrows():
        for var in set(x.replace(f'_{b}', '') for x in rr['feature'].replace('^2', '').split(' ')):
            acc[var] = acc.get(var, 0) + abs(rr['coef'])
    s_ = sum(acc.values()); imp[b] = {k: round(v_ / s_ * 100, 1) for k, v_ in acc.items()}
extra['importancia'] = imp
json.dump(extra, open(os.path.join(OUT, 'numeros_extra.json'), 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
# ---- ENSO (Analise_Biomas_PSN.py): mesmo parser do consolidar_resultados.py
if os.path.exists(a.log_enso):
    la = open(a.log_enso, encoding='utf-8').read(); nm = {'Mata Atlântica': 'MA', 'Cerrado': 'CE', 'Caatinga': 'CA'}
    kw = {b: {} for b in nm.values()}
    s1 = re.search(r'===== KRUSKAL-WALLIS POR FASE ENSO =====(.*?)===== CORRELAÇÃO', la, flags=re.S).group(1)
    for var, blk in re.findall(r'\[ (\w+) \]\n(.*?)(?=\n\[|\Z)', s1, flags=re.S):
        for m in re.finditer(r'(Mata Atlântica|Cerrado|Caatinga): La Niña=([\d.]+) \| Neutro=([\d.]+) \| El Niño=([\d.]+) \| H=([\d.]+) \| p=([\d.]+)', blk):
            b = nm[m.group(1)]; kw[b][var] = float(m.group(6)); kw[b][f'{var}_medias'] = [float(m.group(2)), float(m.group(3)), float(m.group(4))]
    pe, pep = {}, {}
    s2 = re.search(r'===== CORRELAÇÃO DE PEARSON: ONI x VARIÁVEIS =====\n\n\[ PSN \]\n(.*?)\n\[', la, flags=re.S).group(1)
    for m in re.finditer(r'(Mata Atlântica|Cerrado|Caatinga): r=([-\d.]+) \| p=([\d.]+)', s2): pe[nm[m.group(1)]] = float(m.group(2)); pep[nm[m.group(1)]] = float(m.group(3))
    ols = [float(x) for x in re.findall(r'R²=([\d.]+) \(', la)]; disp = {b: {} for b in nm.values()}
    for m in re.finditer(r'^(Precipitação \(mm\)|Evapotranspiração \(mm\)|WAI|Temperatura \(°C\)|Área queimada log\(ha\))\s+(Mata Atlântica|Cerrado|Caatinga)\s+([\d.]+)', la, flags=re.M):
        var = {'Precipitação (mm)': 'PRE', 'Evapotranspiração (mm)': 'EV', 'WAI': 'WAI', 'Temperatura (°C)': 'TST', 'Área queimada log(ha)': 'BURN'}[m.group(1)]
        disp[nm[m.group(2)]][var] = round(float(m.group(3)) * 100, 1)
    for b in disp: disp[b]['max'] = max(disp[b].values())
    fases = dict(re.findall(r'^(Neutro|El Niño|La Niña)\s+(\d+)$', la, flags=re.M))
    enso = dict(kw=kw, pearson=pe, pearson_p=pep, ols_max=round(max(ols) * 100, 2), disp=disp, n_neutro=int(fases['Neutro']), n_elnino=int(fases['El Niño']), n_lanina=int(fases['La Niña']))
    json.dump(enso, open(os.path.join(OUT, 'enso_resumo.json'), 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
    for f in ['boxplot_enso_MA.png', 'boxplot_enso_CE.png', 'boxplot_enso_CA.png', 'dispersao_psn_variaveis_biomas.png', 'serie_temporal_biomas.png']:
        if os.path.exists(os.path.join(SF, f)): shutil.copy(os.path.join(SF, f), os.path.join(OUT, f))
rob_src = os.path.join(SF, 'robustez')
if os.path.isdir(rob_src):                       # bloco de robustez do Modelo_PSN.py (Tabelas A7-A10)
    os.makedirs(os.path.join(OUT, 'robustez'), exist_ok=True)
    for f in os.listdir(rob_src): shutil.copy(os.path.join(rob_src, f), os.path.join(OUT, 'robustez', f))
    print('robustez copiada de saidas_figuras/robustez')
print(pd.DataFrame(rows)[['bioma', 'r2_teste_medio', 'gap_overfitting_pp', 'rmse_teste_medio', 'mae_teste_medio', 'r2_groupkfold_ano', 'r2_timeseriessplit', 'vif_max', 'shapiro_p', 'acf_lag1', 'yrand_diferenca_pp']].to_string(index=False))
print('gravado em', OUT)
