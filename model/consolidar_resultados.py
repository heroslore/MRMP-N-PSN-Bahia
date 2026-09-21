# -*- coding: utf-8 -*-
"""Consolida números dos logs/CSVs em resultados_2001_2025/ (numeros_extra.json, enso_resumo.json)
e incorpora a rodada da Mata Atlântica com EV+TST+WAI (log separado)."""
import os, re, json, shutil, sys
import pandas as pd, numpy as np
BASE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(BASE, 'resultados_2001_2025')
SCR = '/tmp/claude-0/-home-user-Controle-Gas/b644115f-efef-5303-8a82-3af14e08b70e/scratchpad'
LOG_MA = os.path.join(SCR, 'modelo_MA_wai.log')
rg = pd.read_csv(os.path.join(RES, 'resumo_geral.csv'))

def num(pat, txt, grp=1, cast=float):
    m = re.search(pat, txt); return cast(m.group(grp)) if m else None

# ---------- 1. MA com EV+TST+WAI: atualiza resumo_geral, VIF, grau, figuras
if os.path.exists(LOG_MA):
    t = open(LOG_MA, encoding='utf-8').read()
    t = re.sub(r'Linhas de (Treino|Teste) \[[^\]]*\]', '', t)
    row = rg.index[rg.bioma == 'MA'][0]
    upd = {
        'r2_treino_medio': num(r'R2 Treino médio: ([\d.]+)', t), 'r2_treino_dp': num(r'R2 Treino desvio padrão: ([\d.]+)', t),
        'r2aj_treino_medio': num(r'R2 Ajustado Treino médio: ([\d.]+)', t), 'r2_teste_medio': num(r'R2 Teste médio: ([\d.]+)', t),
        'r2_teste_dp': num(r'R2 Teste desvio padrão: ([\d.]+)', t), 'rmse_teste_medio': num(r'RMSE Teste médio: ([\d.]+)', t),
        'mae_teste_medio': num(r'MAE Teste médio: ([\d.]+)', t),
        'ic95_lower': num(r'Intervalo empírico 95% R² Teste: \[([\d.]+)% - ([\d.]+)%\]', t, 1),
        'ic95_upper': num(r'Intervalo empírico 95% R² Teste: \[([\d.]+)% - ([\d.]+)%\]', t, 2),
        'r2_groupkfold_ano': num(r'GroupKFold por ano\s*: ([\d.]+)%', t), 'r2_timeseriessplit': num(r'TimeSeriesSplit expansível: ([\d.]+)%', t),
        'queda_groupkfold_pp': num(r'Queda RKF:GroupKFold\s*: ([-+\d.]+) pp', t), 'queda_timeseries_pp': num(r'Queda RKF:TimeSeries\s*: ([-+\d.]+) pp', t),
        'shapiro_p': num(r'Shapiro-Wilk: stat=[\d.]+, p=([\d.]+)', t), 'ljungbox_p': num(r'Ljung-Box\(12\): p=([\d.]+)', t),
        'acf_lag1': num(r'ACF lag1=([-+\d.]+)', t), 'yrand_diferenca_pp': num(r'Diferença \(Δ CV\)\s*: ([\d.]+) pp', t),
        'yrand_p_empirico': num(r'p-valor empírico unilateral\s*: ([\d.]+)', t), 'yrand_status': 'APROVADO', 'arquivo_coef': 'coeficientes_ridge_MA.csv',
    }
    upd['gap_overfitting_pp'] = round(upd['r2_treino_medio'] - upd['r2_teste_medio'], 2)
    upd['residuos_normais'] = 'Sim' if upd['shapiro_p'] > 0.05 else 'Não'
    vif_blk = re.search(r'===== VIF =====\n(.*?)\nAlpha', t, flags=re.S).group(1)
    vifs = re.findall(r'\d+\s+(\S+)\s+([\d.]+)', vif_blk); upd['vif_max'] = max(float(x) for _, x in vifs)
    for k, val in upd.items(): rg.loc[row, k] = val
    vif = pd.read_csv(os.path.join(RES, 'vif_por_bioma.csv')); vif = vif[vif.bioma != 'MA']
    vif = pd.concat([pd.DataFrame([{'bioma': 'MA', 'variavel': a, 'VIF': round(float(b), 3)} for a, b in vifs]), vif]).reset_index(drop=True)
    vif.to_csv(os.path.join(RES, 'vif_por_bioma.csv'), index=False)
    sg = pd.read_csv(os.path.join(RES, 'selecao_grau.csv')); sg = sg[sg.bioma != 'MA']
    rows = []
    for ln in re.search(r' Grau \| Termos.*?\n=+\n(.*?)\n=+\n', t, flags=re.S).group(1).split('\n'):
        m = re.match(r'\s*(\d)\s*\|\s*(\d+)\s*\|\s*([-\d.]+)%\s*\|\s*([-\d.]+)%\s*\|\s*([-\d.]+)%\s*\|\s*([-\d.]+)p\.p\s*\|\s*([-\d.]+)\s*\|\s*\[([-\d.]+)% - ([-\d.]+)%\]', ln)
        if m: rows.append(['MA'] + [float(x) for x in m.groups()])
    sg = pd.concat([pd.DataFrame(rows, columns=sg.columns), sg]).reset_index(drop=True); sg.to_csv(os.path.join(RES, 'selecao_grau.csv'), index=False)
    SF = os.path.join(BASE, 'saidas_figuras', 'MA')
    def ultimo(prefixo, ext):
        c = sorted([f for f in os.listdir(SF) if f.startswith(prefixo) and f.endswith(ext)], key=lambda f: int(re.search(r'_(\d+)\.', f).group(1)))
        return os.path.join(SF, c[-1])
    for pref, ext, dest in [('coeficientes_ridge_MA', '.csv', 'coeficientes_ridge_MA.csv'), ('analise_residuos_MA', '.png', 'analise_residuos_MA.png'),
                            ('obs_vs_pred_MA', '.png', 'obs_vs_pred_MA.png'), ('yrandomization_MA', '.png', 'yrandomization_MA.png'),
                            ('selecao_grau_polinomial_MA', '.png', 'selecao_grau_polinomial_MA.png')]:
        shutil.copy(ultimo(pref, ext), os.path.join(RES, 'MA', dest))
    open(os.path.join(RES, 'log_modelo_psn_MA_EV_TST_WAI.txt'), 'w', encoding='utf-8').write(re.sub(r'\n{3,}', '\n\n', t))
    rg.to_csv(os.path.join(RES, 'resumo_geral.csv'), index=False)
    print('MA (EV+TST+WAI) incorporada:', {k: upd[k] for k in ['r2_teste_medio', 'gap_overfitting_pp', 'vif_max', 'shapiro_p', 'acf_lag1', 'yrand_diferenca_pp']})

# ---------- 2. numeros_extra.json
rg = rg.set_index('bioma')
logs = {'MA': open(LOG_MA, encoding='utf-8').read() if os.path.exists(LOG_MA) else '',
        'CE': open(os.path.join(RES, 'log_modelo_psn.txt'), encoding='utf-8').read(),
        'CA': open(os.path.join(RES, 'log_modelo_psn.txt'), encoding='utf-8').read()}
extra = {}
for b in ['MA', 'CE', 'CA']:
    t = logs[b]; r2t = rg.loc[b, 'r2_teste_medio']
    # MAE dp: linha após 'MAE Teste médio: X' com X igual ao do bioma
    mae_dp = None
    for m in re.finditer(r'MAE Teste médio: ([\d.]+)\nMAE Teste desvio padrão: ([\d.]+)', t):
        if abs(float(m.group(1)) - rg.loc[b, 'mae_teste_medio']) < 0.011: mae_dp = float(m.group(2))
    shap_w = None
    for m in re.finditer(r'Shapiro-Wilk: stat=([\d.]+), p=([\d.]+)', t):
        if abs(float(m.group(2)) - rg.loc[b, 'shapiro_p']) < 0.0011: shap_w = float(m.group(1))
    yr = None
    for m in re.finditer(r'modelo original\s*: ([-\d.]+)%\nR² CV médio - modelos permutados\s*: ([-\d.]+)% ± ([\d.]+)', t):
        if abs(float(m.group(1)) - r2t) < 1.5: yr = (float(m.group(1)), float(m.group(2)), float(m.group(3)))
    oof = None
    for m in re.finditer(r'R² out-of-fold \(Figura 3\): ([\d.]+)%\s*\|\s*R² RepeatedKFold \(tabelas\): ([\d.]+)%', t):
        if abs(float(m.group(2)) - r2t) < 0.011: oof = float(m.group(1))
    pares = re.findall(r'Alpha médio utilizado na análise de resíduos: ([\d.]+)\nArquivo coeficientes_ridge_(\w\w)', t)
    alpha_medio = [float(a) for a, bb in pares if bb == b][-1]
    coef = pd.read_csv(os.path.join(RES, b, f'coeficientes_ridge_{b}.csv'))
    xv = [f.replace(f'_{b}', '') for f in coef['feature'].iloc[:3]]
    extra[b] = dict(x=xv, mae_dp=mae_dp, shapiro_w=shap_w, yr_orig=yr[0], yr_perm_media=yr[1], yr_perm_dp=yr[2], oof_r2=oof, alpha_medio=alpha_medio)
imp = {}
for b in ['MA', 'CE', 'CA']:
    c = pd.read_csv(os.path.join(RES, b, f'coeficientes_ridge_{b}.csv')); acc = {}
    for _, r in c.iterrows():
        for var in set(x.replace(f'_{b}', '') for x in r['feature'].replace('^2', '').split(' ')):
            acc[var] = acc.get(var, 0) + abs(r['coef'])
    s = sum(acc.values()); imp[b] = {k: round(v_ / s * 100, 1) for k, v_ in acc.items()}
extra['importancia'] = imp
json.dump(extra, open(os.path.join(RES, 'numeros_extra.json'), 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
print(json.dumps(extra, ensure_ascii=False))

# ---------- 3. enso_resumo.json
la = open(os.path.join(RES, 'log_analise_biomas.txt'), encoding='utf-8').read()
nm = {'Mata Atlântica': 'MA', 'Cerrado': 'CE', 'Caatinga': 'CA'}
kw = {b: {} for b in nm.values()}
sec = re.search(r'===== KRUSKAL-WALLIS POR FASE ENSO =====(.*?)===== CORRELAÇÃO', la, flags=re.S).group(1)
for var, blk in re.findall(r'\[ (\w+) \]\n(.*?)(?=\n\[|\Z)', sec, flags=re.S):
    for m in re.finditer(r'(Mata Atlântica|Cerrado|Caatinga): La Niña=([\d.]+) \| Neutro=([\d.]+) \| El Niño=([\d.]+) \| H=([\d.]+) \| p=([\d.]+)', blk):
        b = nm[m.group(1)]; kw[b][var] = float(m.group(6)); kw[b][f'{var}_medias'] = [float(m.group(2)), float(m.group(3)), float(m.group(4))]
pe, pep = {}, {}
sec = re.search(r'===== CORRELAÇÃO DE PEARSON: ONI x VARIÁVEIS =====\n\n\[ PSN \]\n(.*?)\n\[', la, flags=re.S).group(1)
for m in re.finditer(r'(Mata Atlântica|Cerrado|Caatinga): r=([-\d.]+) \| p=([\d.]+)', sec):
    pe[nm[m.group(1)]] = float(m.group(2)); pep[nm[m.group(1)]] = float(m.group(3))
ols = [float(x) for x in re.findall(r'R²=([\d.]+) \(', la)]
disp = {b: {} for b in nm.values()}
for m in re.finditer(r'^(Precipitação \(mm\)|Evapotranspiração \(mm\)|WAI|Temperatura \(°C\)|Área queimada log\(ha\))\s+(Mata Atlântica|Cerrado|Caatinga)\s+([\d.]+)', la, flags=re.M):
    var = {'Precipitação (mm)': 'PRE', 'Evapotranspiração (mm)': 'EV', 'WAI': 'WAI', 'Temperatura (°C)': 'TST', 'Área queimada log(ha)': 'BURN'}[m.group(1)]
    disp[nm[m.group(2)]][var] = round(float(m.group(3)) * 100, 1)
for b in disp: disp[b]['max'] = max(disp[b].values())
fases = dict(re.findall(r'^(Neutro|El Niño|La Niña)\s+(\d+)$', la, flags=re.M))
enso = dict(kw=kw, pearson=pe, pearson_p=pep, ols_max=round(max(ols) * 100, 2), disp=disp,
            n_neutro=int(fases['Neutro']), n_elnino=int(fases['El Niño']), n_lanina=int(fases['La Niña']))
json.dump(enso, open(os.path.join(RES, 'enso_resumo.json'), 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
print(json.dumps(enso, ensure_ascii=False)[:600])
