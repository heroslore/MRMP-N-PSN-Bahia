# -*- coding: utf-8 -*-
"""
Correcao do teste pareado entre as duas melhores combinacoes de preditores (Artigo 1).

A comparacao reportada usa as 150 particoes do RepeatedKFold 5 x 30 como se fossem
150 observacoes independentes. Nao sao: as 30 repeticoes reembaralham os MESMOS 297
meses, de modo que os conjuntos de treino de particoes diferentes se sobrepoem
fortemente. Isso subestima a variancia da diferenca e infla a significancia -- o
mesmo problema que o artigo companheiro trata ao passar do mes para o episodio.

A correcao padrao e o teste t reamostrado corrigido de Nadeau e Bengio (2003), que
substitui a variancia 1/J por (1/J + n2/n1), onde n1 e n2 sao os tamanhos de treino
e de teste de cada particao.

Saida: resultados_2001_2025/robustez/diferenca_pareada_corrigida.csv
"""
import os, sys
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import r2_score

BASE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BASE)
os.environ.setdefault('BIOMA_ATIVO', 'TODOS')
import Modelo_PSN as M

OUT = os.path.join(BASE, 'resultados_2001_2025', 'robustez')
NOME = {'MA': 'Mata Atlântica', 'CE': 'Cerrado', 'CA': 'Caatinga'}

dif = pd.read_csv(os.path.join(OUT, 'diferenca_pareada.csv'))
sig = {v: k for k, v in NOME.items()}

def r2_por_particao(bioma, combo, dd):
    """R2 de teste em cada uma das 150 particoes, identico a _rob_tarefa_combo."""
    x = dd[M._rob_colunas(bioma, combo)]; y = dd[f'PSN_{bioma}']
    anos = dd['ANO'].values; out = []
    for tr, te in M._rob_splits('RepeatedKFold', x, anos):
        sc, pf, m, _ = M._rob_fit(x.iloc[tr], y.iloc[tr].copy(), 2)
        out.append(r2_score(y.iloc[te], M._rob_pred(sc, pf, m, x.iloc[te])) * 100)
    return np.array(out)

linhas = []
for _, r in dif.iterrows():
    b = sig[r['bioma']]
    dd = M._rob_pacote(M.dados_total_base, b)
    c1 = tuple(r['primeira'].split(' + ')); c2 = tuple(r['segunda'].split(' + '))
    d = r2_por_particao(b, c1, dd) - r2_por_particao(b, c2, dd)
    J = len(d); md = d.mean(); s2 = d.var(ddof=1)

    # tamanhos de treino e teste das particoes (5 folds sobre n meses)
    x = dd[M._rob_colunas(b, c1)]
    n1, n2 = np.mean([[len(tr), len(te)] for tr, te in M._rob_splits('RepeatedKFold', x, dd['ANO'].values)], axis=0)

    # (a) ingenuo: particoes tratadas como independentes
    se_ing = np.sqrt(s2 / J)
    t_ing = md / se_ing
    p_ing = 2 * stats.t.sf(abs(t_ing), J - 1)
    p_wil = stats.wilcoxon(d).pvalue

    # (b) Nadeau-Bengio: variancia corrigida pela sobreposicao dos treinos
    se_nb = np.sqrt((1 / J + n2 / n1) * s2)
    t_nb = md / se_nb
    p_nb = 2 * stats.t.sf(abs(t_nb), J - 1)
    ic_nb = md + np.array([-1, 1]) * stats.t.ppf(0.975, J - 1) * se_nb

    linhas.append(dict(
        bioma=r['bioma'], primeira=r['primeira'], segunda=r['segunda'],
        dif_media_pp=round(md, 4), n_particoes=J, n_treino=round(n1), n_teste=round(n2),
        ic95_ingenuo_inf=round(md - 1.96 * se_ing, 4), ic95_ingenuo_sup=round(md + 1.96 * se_ing, 4),
        p_wilcoxon=float('%.3g' % p_wil), p_t_ingenuo=float('%.3g' % p_ing),
        ic95_nb_inf=round(ic_nb[0], 4), ic95_nb_sup=round(ic_nb[1], 4),
        p_nadeau_bengio=round(p_nb, 4),
        razao_erro_padrao=round(se_nb / se_ing, 2),
        prop_particoes_primeira_maior=round(100 * (d > 0).mean(), 1),
        sig_nb=bool(ic_nb[0] > 0 or ic_nb[1] < 0)))

T = pd.DataFrame(linhas)
T.to_csv(os.path.join(OUT, 'diferenca_pareada_corrigida.csv'), index=False, encoding='utf-8')
pd.set_option('display.width', 200)
print(T[['bioma', 'dif_media_pp', 'ic95_ingenuo_inf', 'ic95_ingenuo_sup', 'p_wilcoxon',
         'ic95_nb_inf', 'ic95_nb_sup', 'p_nadeau_bengio', 'razao_erro_padrao', 'sig_nb']].to_string(index=False))
print('\nsalvo em', os.path.join(OUT, 'diferenca_pareada_corrigida.csv'))
