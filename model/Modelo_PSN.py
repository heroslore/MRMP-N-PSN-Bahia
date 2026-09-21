# -*- coding: utf-8 -*-
"""
Modelo de Regressão Ridge Polinomial para estimativa de PSN
Localidade: BA, Mata Atlântica, Cerrado e Caatinga | Período: 2001-2025

Autor:
PPGCTA - Programa de Pós-Graduação em Ciências e Tecnologia Ambiental

-------------------------------------------------------------------------------
BASE DE DADOS DE ENTRADA
-------------------------------------------------------------------------------
Basta colocar o CSV bruto atualizado (base_final_2001_2025_plan1_excel_ptbr.csv)
na mesma pasta do script e rodar — nada mais precisa ser feito manualmente. O
script detecta o CSV sozinho, renomeia as colunas para o padrão interno,
remove as linhas com precipitação ainda não publicada pela NASA, e salva o
Excel (Dados_base_nova_2001_2025.xlsx) que o modelo de fato usa. Se você
apagar o CSV e deixar só o .xlsx já pronto, o script usa ele direto.

-------------------------------------------------------------------------------
COMO ESCOLHER O(S) BIOMA(S) A RODAR
-------------------------------------------------------------------------------
BIOMA_ATIVO = os.environ.get('BIOMA_ATIVO', 'CE')

- Rodar só um bioma: dê o "Run" normal — ele cai no padrão 'CE'. Para rodar
  outro, troque o 'CE' logo abaixo por 'MA' ou 'CA' e rode novamente.
- Rodar os 3 de uma vez: troque o 'CE' abaixo por 'TODOS'. Por padrão os três
  rodam EM PARALELO (RODAR_EM_PARALELO = os.environ.get('RODAR_EM_PARALELO', '1') == '1'   # 0 = biomas em sequência (log ordenado), um processo por núcleo) — cada
  bioma já lê/escreve só na própria subpasta, então não há conflito. Para
  rodar em sequência em vez de paralelo, troque RODAR_EM_PARALELO para False.
- Também dá pra rodar pelo terminal sem editar nada, definindo a variável de
  ambiente antes do comando, ex.: BIOMA_ATIVO=TODOS python3 "Modelo_PSN.py"

-------------------------------------------------------------------------------
ONDE OS ARQUIVOS SÃO SALVOS
-------------------------------------------------------------------------------
Cada bioma tem sua própria subpasta dentro de "saidas_figuras/":

    saidas_figuras/
    ├── MA/
    │   ├── coeficientes_ridge_MA_1.csv
    │   ├── analise_residuos_MA_1.png
    │   ├── obs_vs_pred_MA_1.png
    │   ├── yrandomization_MA_1.png        (se RODAR_YRANDOMIZATION=True)
    │   ├── selecao_grau_polinomial_MA_1.png (se TESTAR_GRAUS=True)
    │   └── ...
    ├── CE/
    │   └── ...
    ├── CA/
    │   └── ...
    └── resumo_geral_1.csv   (só quando roda 'TODOS' — tabela comparativa)

Nenhum arquivo é sobrescrito: a cada nova execução, o número no final do
nome (_1, _2, _3...) avança automaticamente. O histórico de rodadas antigas
fica preservado.

-------------------------------------------------------------------------------
TOGGLES DE VALIDAÇÃO (mantidos do script original)
-------------------------------------------------------------------------------
TESTAR_GRAUS          — compara os graus polinomiais 1 a GRAU_MAXIMO_TESTE.
RODAR_YRANDOMIZATION  — roda as 100 permutações com validação cruzada
                         (Y-randomization fora da amostra, conforme a
                         metodologia descrita no texto da dissertação).
Ambos aplicam-se igualmente a todos os biomas selecionados.
"""

# =============================================================================
# IMPORTAÇÕES
# =============================================================================

import os
import json
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # backend não-interativo — necessário para rodar em paralelo (vários processos)
import matplotlib.pyplot as plt
import scipy.stats as stats

from itertools import combinations
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed

from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
from sklearn.preprocessing import StandardScaler, PolynomialFeatures
from sklearn.linear_model import Ridge
from sklearn.model_selection import RepeatedKFold, GridSearchCV, TimeSeriesSplit, GroupKFold, KFold

from statsmodels.stats.outliers_influence import variance_inflation_factor

# =============================================================================
# DIRETÓRIO BASE E FUNÇÕES AUXILIARES DE CAMINHO
# =============================================================================

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

def caminho(nome):
    return os.path.join(BASE_DIR, nome)

# =============================================================================
# PASTA DE SAÍDA DAS FIGURAS — uma subpasta por bioma, nomes numerados
# A cada execução é garantida a existência de uma pasta 'saidas_figuras/<bioma>'
# na mesma pasta do script. Nenhum arquivo é sobrescrito: o número final do
# nome (_1, _2, _3...) avança sozinho a cada nova rodada.
# =============================================================================

PASTA_SAIDA = os.path.join(BASE_DIR, 'saidas_figuras')
os.makedirs(PASTA_SAIDA, exist_ok=True)


def pasta_do_bioma(bioma):
    """Retorna (e cria, se preciso) a subpasta de saída de um bioma."""
    pasta = os.path.join(PASTA_SAIDA, bioma)
    os.makedirs(pasta, exist_ok=True)
    return pasta


NUMERAR_SAIDAS = False

def proximo_nome(pasta, nome_base, extensao):
    """
    Gera o próximo nome disponível no formato 'nome_base_N.extensao',
    sem nunca sobrescrever um arquivo já existente (N cresce a cada rodada).
    """
    # Nomes fixos (cada rodada sobrescreve a anterior), conforme decidido pelo autor.
    # Para voltar à numeração _1, _2, ..., defina NUMERAR_SAIDAS = True.
    if not NUMERAR_SAIDAS:
        return os.path.join(pasta, f"{nome_base}.{extensao}")
    n = 1
    while os.path.exists(os.path.join(pasta, f"{nome_base}_{n}.{extensao}")):
        n += 1
    return os.path.join(pasta, f"{nome_base}_{n}.{extensao}")


# =============================================================================
# ESTILO GLOBAL DOS GRÁFICOS (fontes ampliadas para leitura em impressão)
# Aplica-se a TODAS as figuras geradas pelo script.
# =============================================================================

plt.rcParams.update({
    'font.size':        13,
    'axes.titlesize':   16,
    'axes.titleweight': 'bold',
    'axes.labelsize':   14,
    'axes.labelweight': 'bold',
    'xtick.labelsize':  12,
    'ytick.labelsize':  12,
    'legend.fontsize':  12,
    'figure.titlesize': 18,
    'figure.titleweight': 'bold',
})

# =============================================================================
# PREPARAÇÃO DA BASE A PARTIR DO CSV BRUTO (opcional — só roda se o CSV existir)
# -----------------------------------------------------------------------------
# Se o arquivo base_final_2001_2025_plan1_excel_ptbr.csv estiver na pasta do
# script, ele é lido, as colunas são renomeadas para o padrão que o modelo
# espera (MÊS, NP_*, WAI_*, EV_*, PRE_*, TST_*, BURN_*), e as linhas com
# precipitação ausente (meses recentes que a NASA ainda não publicou) são
# removidas automaticamente. O resultado é salvo como
# Dados_base_nova_2001_2025.xlsx, sobrescrevendo a versão anterior.
#
# Se o CSV não existir na pasta, o script usa direto o
# Dados_base_nova_2001_2025.xlsx existente, sem precisar do CSV.
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
    # FMA_PET, Cerrado_PET, Caatinga_PET não são usadas pelo modelo —
    # mantidas com o nome original e simplesmente ignoradas.
}

# Fase ENSO oficial (NOAA): mesma coluna 'Enso' que o Analise_Biomas_PSN.py
# grava, para que os dois scripts produzam o MESMO Dados_base_nova_2001_2025.xlsx.
from enso_noaa import classificar_fase_enso_noaa

COLUNAS_PRECIP = ['PRE_MA', 'PRE_CE', 'PRE_CA']


def preparar_base_a_partir_do_csv():
    """
    Lê o CSV bruto, renomeia as colunas para o padrão do modelo, remove as
    linhas com precipitação ausente e salva o Excel que o modelo usa.
    Retorna True se o CSV foi encontrado e processado, False caso contrário.
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
    print(f"Linhas finais usadas no modelo: {depois}")

    df.to_excel(caminho(NOME_XLSX_MODELO), index=False)
    print(f"Base atualizada salva em {NOME_XLSX_MODELO}\n")
    return True


if __name__ == '__main__':
    csv_processado = preparar_base_a_partir_do_csv()
    if not csv_processado:
        print(f"\n(CSV '{NOME_CSV_BRUTO}' não encontrado na pasta — "
              f"usando '{NOME_XLSX_MODELO}' existente diretamente.)\n")

# =============================================================================
# CARREGAMENTO E PRÉ-PROCESSAMENTO DOS DADOS (feito uma única vez)
# =============================================================================

dados_total_base = pd.read_excel(caminho(NOME_XLSX_MODELO))

# Remove espaços invisíveis nos nomes das colunas (problema comum em planilhas Excel)
dados_total_base.columns = dados_total_base.columns.str.strip()

# Padroniza o nome da variável resposta como PSN (Fotossíntese Líquida) por bioma.
# A planilha de origem usa o prefixo NP_; aqui mapeamos para PSN_ sem alterar o arquivo.
dados_total_base = dados_total_base.rename(columns={
    'NP_MA': 'PSN_MA',
    'NP_CE': 'PSN_CE',
    'NP_CA': 'PSN_CA',
})

# Sazonalidade harmônica: decompõe o ciclo anual em duas componentes ortogonais.
# saz_sin e saz_cos descrevem completamente qualquer ciclo periódico de 12 meses
# sem impor linearidade no tempo (diferente de usar o número do mês diretamente).
dados_total_base['saz_sin'] = np.sin(2 * np.pi * dados_total_base['MÊS'] / 12)
dados_total_base['saz_cos'] = np.cos(2 * np.pi * dados_total_base['MÊS'] / 12)

# Transformação log(1+x) em área queimada: normaliza a distribuição fortemente
# assimétrica de BURN_MA, que concentra zeros e possui raros picos muito elevados.
dados_total_base['BURN_MA_log'] = np.log1p(dados_total_base['BURN_MA'])
dados_total_base['BURN_CE_log'] = np.log1p(dados_total_base['BURN_CE'])
dados_total_base['BURN_CA_log'] = np.log1p(dados_total_base['BURN_CA'])

# =============================================================================
# CONFIGURAÇÃO DOS BIOMAS
# =============================================================================

config_biomas = {
    'MA': {
        'nome':     'Mata Atlântica',
        'y':        'PSN_MA',
        # Conjunto ótimo pela busca exaustiva C(5,3) com a base 2001-2025
        # (Selecao_Variaveis_PSN.py): EV + TST + WAI (R² teste 73,7%) supera
        # EV + PRE + TST (71,0%), conjunto ótimo da base 2001-2020.
        'x':        ['EV_MA', 'TST_MA', 'WAI_MA', 'saz_sin', 'saz_cos'],
        'winsor_y': 'PSN_MA',
        'cor':      '#2CA02C',   # verde
    },
    'CE': {
        'nome':     'Cerrado',
        'y':        'PSN_CE',
        'x':        ['EV_CE', 'PRE_CE', 'WAI_CE', 'saz_sin', 'saz_cos'],
        'winsor_y': 'PSN_CE',
        'cor':      '#FF7F0E',   # laranja
    },
    'CA': {
        'nome':     'Caatinga',
        'y':        'PSN_CA',
        'x':        ['EV_CA', 'PRE_CA', 'TST_CA', 'saz_sin', 'saz_cos'],
        'winsor_y': 'PSN_CA',
        'cor':      '#D62728',   # vermelho
    },
}

# =============================================================================
# CONFIGURAÇÃO DO GRAU POLINOMIAL, VALIDAÇÃO E Y-RANDOMIZATION (comum aos biomas)
# GRAU_MODELO : grau usado no modelo principal (validação cruzada completa).
#                Altere este valor para rodar o modelo com o grau desejado.
#
# TESTAR_GRAUS: ativa a comparação automática de graus ao final do script.
#                True  = compara todos os graus de 1 até GRAU_MAXIMO_TESTE.
#                False = pula a comparação (execução mais rápida).
#
# GRAU_MAXIMO_TESTE: até qual grau comparar quando TESTAR_GRAUS = True.
#                     Recomendado: máximo 5 (acima disso há overfitting severo
#                     com cerca de 300 meses e 5 variáveis preditoras).
#
# RODAR_YRANDOMIZATION: True roda as 100 permutações com validação cruzada
#                        (mais lento); False pula a Y-randomization.
# =============================================================================

GRAU_MODELO       = int(os.environ.get('GRAU_MODELO', 2))     # grau do modelo principal (ou variável de ambiente)
TESTAR_GRAUS      = os.environ.get('TESTAR_GRAUS', '0') == '1'   # 1 para comparar graus 1..GRAU_MAXIMO_TESTE
GRAU_MAXIMO_TESTE = 5     # até qual grau comparar (válido se TESTAR_GRAUS=True)

RODAR_YRANDOMIZATION = os.environ.get('RODAR_YRANDOMIZATION', '1') == '1'   # 0 para pular as 100 permutações

# RODAR_ROBUSTEZ: ao final da rodada, executa o bloco "ROBUSTEZ DA SELEÇÃO E
# COMPLEMENTOS DE VALIDAÇÃO" (Tabelas A7 a A10 da dissertação): seleção do grau e
# das combinações de variáveis sob GroupKFold por ano e TimeSeriesSplit, estabilidade
# da seleção partição a partição, importância por permutação por bloco temporal e
# nulos que preservam a estrutura temporal da PSN. Demora ~15 min com 4 núcleos.
RODAR_ROBUSTEZ = os.environ.get('RODAR_ROBUSTEZ', '1') == '1'   # 0 para pular

numero_colunas_agrupamento = 5

# =============================================================================
# BIOMA(S) A RODAR
# 'MA', 'CE', 'CA'  → roda só aquele bioma
# 'TODOS'           → roda os três em sequência (ou em paralelo, ver abaixo)
# =============================================================================

BIOMA_ATIVO = os.environ.get('BIOMA_ATIVO', 'TODOS')  # altere aqui: 'MA', 'CE', 'CA' ou 'TODOS'

# Só faz diferença quando BIOMA_ATIVO = 'TODOS'.
# True  → roda os 3 biomas ao mesmo tempo, em processos separados (usa 3 núcleos).
# False → roda um de cada vez, em sequência (mais lento, usa 1 núcleo por vez).
# Seguro em paralelo porque cada bioma já lê/escreve só na própria subpasta
# (saidas_figuras/MA, saidas_figuras/CE, saidas_figuras/CA).
RODAR_EM_PARALELO = os.environ.get('RODAR_EM_PARALELO', '1') == '1'   # 0 = biomas em sequência (log ordenado)

if BIOMA_ATIVO.upper() == 'TODOS':
    biomas_para_rodar = ['MA', 'CE', 'CA']
else:
    biomas_para_rodar = [BIOMA_ATIVO]


# =============================================================================
# FUNÇÃO PRINCIPAL — roda o modelo completo para UM bioma
# (mesma lógica do script original; caminhos de arquivo e nomes de variável
# adaptados para pasta própria por bioma e numeração incremental)
# =============================================================================

def rodar_modelo_bioma(bioma, dados_total):

    pasta = pasta_do_bioma(bioma)

    cfg       = config_biomas[bioma]
    colunas_y = [cfg['y']]
    colunas_x = cfg['x']

    print(f"\n===== RODANDO MODELO RIDGE - {cfg['nome'].upper()} =====")
    print(f"Variável resposta : {cfg['y']}")
    print(f"Preditores        : {colunas_x}\n")
    # Winsorização da variável resposta: o percentil de corte é calculado dentro de
    # cada fold, apenas sobre o conjunto de treino (ver loop de validação cruzada),
    # de modo que as observações de teste não influenciam o limite aplicado.

    dados_x = dados_total[colunas_x]
    dados_y = dados_total[colunas_y]
    # =============================================================================
    # CONJUNTO DE PREDITORES
    # O conjunto ótimo de 5 preditores por bioma foi definido pela busca exaustiva
    # entre as C(7,5) = 21 combinações das 7 candidatas (EV, PRE, TST, WAI, BURN_log,
    # saz_sin, saz_cos). Aqui já se utiliza esse conjunto final (definido em config_biomas),
    # de modo que a expansão a seguir opera sobre as 5 variáveis selecionadas.
    # =============================================================================

    numero_colunas_agrupamento = 5
    lista_comb_colunas = list(combinations(colunas_x, numero_colunas_agrupamento))
    lista_r2_treino    = []
    lista_r2_teste     = []
    lista_r2aj_treino  = []
    lista_rmse_teste   = []
    lista_mae_teste    = []
    alphas_escolhidos  = []
    for f in ["melhores_combinacoes_20.txt", "melhores_combinacoes_60.txt", "melhores_combinacoes.txt"]:
        open(os.path.join(pasta, f), "w").close()
    # =============================================================================
    # CONFIGURAÇÃO DA VALIDAÇÃO CRUZADA
    # RepeatedKFold com 5 splits e 30 repetições = 150 partições de validação.
    # As partições são parcialmente sobrepostas (não independentes), por isso os
    # percentis dos escores são reportados como intervalo empírico, não como
    # intervalo de confiança no sentido estatístico estrito.
    # =============================================================================

    rkf = RepeatedKFold(
        n_splits=5,
        n_repeats=30,
        random_state=42        # semente fixa para reprodutibilidade
    )

    # Contadores e controle de progresso
    qtd_rodadas_linhas = 150                               # total de folds (5 splits × 30 repeats)
    tot_rodadas        = qtd_rodadas_linhas * len(lista_comb_colunas)
    rodada             = 0
    data_e_hora_atual  = datetime.now()

    # Variáveis para rastrear o melhor resultado por fold
    r2_treino_melhor                = 0
    r2_teste_melhor                 = 0
    melhor_parametos_equacao_treino = ''

    # =============================================================================
    # LOOP PRINCIPAL DE VALIDAÇÃO CRUZADA
    # =============================================================================

    for train_index, test_index in rkf.split(dados_x):

        # Separação treino/teste para X e Y
        x_train  = dados_x.iloc[train_index]
        x_test   = dados_x.iloc[test_index]
        y_treino = dados_y.iloc[train_index]
        y_teste  = dados_y.iloc[test_index]

        # Índices das observações de treino e teste (registrados no log de saída)
        linhas_treino = str(sorted([int(v) for v in x_train.index]))
        linhas_teste  = str(sorted([int(v) for v in x_test.index]))

        melhor_resultado_treino = '\n\nNenhum resultado com R² treino e teste > 60%'

        # Itera sobre o conjunto de preditores selecionado para o bioma.
        for i in lista_comb_colunas:

            x_treino = x_train[list(i)]
            x_teste  = x_test[list(i)]

            # Garante alinhamento correto de índices entre X e Y
            y_treino = dados_y.loc[x_treino.index].copy()
            y_teste  = dados_y.loc[x_teste.index].copy()

            # ---------------------------------------------------------------
            # WINSORIZAÇÃO SEM LEAKAGE: o limite do percentil 3% é calculado
            # EXCLUSIVAMENTE no conjunto de treino e aplicado apenas a ele.
            # O conjunto de teste permanece intacto, garantindo avaliação
            # realista do desempenho fora da amostra.
            # ---------------------------------------------------------------
            limite_winsor = np.percentile(y_treino.values, 3)
            y_treino = y_treino.clip(lower=limite_winsor)

            rodada += 1

            # ---------------------------------------------------------------
            # ETAPA 1: Padronização (StandardScaler)
            # Transforma cada variável para média 0 e desvio-padrão 1.
            # Necessário para que a regularização Ridge trate os coeficientes
            # de forma equitativa, independente da escala original das variáveis.
            # ---------------------------------------------------------------
            scaler          = StandardScaler()
            x_treino_scaled = scaler.fit_transform(x_treino)  # ajusta e transforma treino
            x_teste_scaled  = scaler.transform(x_teste)       # apenas transforma (sem reajuste)

            # ---------------------------------------------------------------
            # ETAPA 2: Expansão Polinomial
            # O grau é definido pela variável GRAU_MODELO no topo do script.
            # Gera termos polinomiais e de interação entre variáveis.
            # ---------------------------------------------------------------
            poly          = PolynomialFeatures(degree=GRAU_MODELO, include_bias=False)
            X_treino_poly = poly.fit_transform(x_treino_scaled)
            X_teste_poly  = poly.transform(x_teste_scaled)

            # ---------------------------------------------------------------
            # ETAPA 3: Regressão Ridge com seleção automática de alpha
            # GridSearchCV testa 5 valores de alpha via validação cruzada interna (cv=5).
            # O alpha controla a força da regularização L2 que penaliza coeficientes
            # grandes, evitando overfitting causado pela expansão polinomial.
            # ---------------------------------------------------------------
            alphas   = [0.1, 1.0, 10.0, 50.0, 100.0]
            ridge_cv = GridSearchCV(Ridge(),
                                    param_grid={'alpha': alphas},
                                    cv=5,
                                    scoring='r2')
            ridge_cv.fit(X_treino_poly, y_treino.values.ravel())

            # Reutiliza o melhor estimador e o re-treina no conjunto completo de treino
            maquina_preditiva = ridge_cv.best_estimator_
            maquina_preditiva.fit(X_treino_poly, y_treino.values.ravel())
            alphas_escolhidos.append(ridge_cv.best_params_['alpha'])

            # ---------------------------------------------------------------
            # PREDIÇÃO E CÁLCULO DO R²
            # ---------------------------------------------------------------
            y_pred   = maquina_preditiva.predict(X_treino_poly)
            r2       = r2_score(y_treino, y_pred) * 100           # R² treino (%)

            y_pred_teste = maquina_preditiva.predict(X_teste_poly)
            r2_teste     = r2_score(y_teste, y_pred_teste) * 100  # R² teste (%)

            # ---------------------------------------------------------------
            # MÉTRICAS COMPLEMENTARES
            # ---------------------------------------------------------------

            # Erros de predição no conjunto de teste
            rmse_teste = np.sqrt(mean_squared_error(y_teste, y_pred_teste))
            mae_teste  = mean_absolute_error(y_teste, y_pred_teste)

            # Tamanhos dos conjuntos e número de preditores após expansão polinomial
            n_tr = len(y_treino)
            n_te = len(y_teste)
            k    = X_treino_poly.shape[1]   # 20 termos para 5 variáveis (grau 2)

            R2_tr = r2 / 100
            R2_te = r2_teste / 100

            # R² ajustado: penaliza pelo número de preditores, evitando inflação do R²
            r2aj_tr = 1 - (1 - R2_tr) * (n_tr - 1) / (n_tr - k - 1)

            # A estatística F e o p-valor não são calculados: a fórmula clássica
            # pressupõe estimação por MQO, com graus de liberdade iguais a k. Sob
            # regularização L2 os graus de liberdade efetivos são menores que k, o
            # que torna F/p-valor inadequados. O desempenho é avaliado por R² de
            # validação cruzada, RMSE, MAE e Y-randomization.

            # Armazena resultados de todos os folds para as estatísticas finais
            lista_r2_treino.append(r2)
            lista_r2_teste.append(r2_teste)
            lista_r2aj_treino.append(r2aj_tr * 100)
            lista_rmse_teste.append(rmse_teste)
            lista_mae_teste.append(mae_teste)

            # ---------------------------------------------------------------
            # FORMATAÇÃO DO RESULTADO E CONTROLE DE TEMPO
            # ---------------------------------------------------------------
            resultado = (
                f'\n#################  Rodada {rodada} / {tot_rodadas}  '
                f'{round((rodada / tot_rodadas) * 100, 2)} %'
                f'\nParâmetros {i}'
                f'\nLinhas de Treino {linhas_treino}'
                f'\nLinhas de Teste {linhas_teste}'
                f'\nR2 de Treino = {round(r2, 2)}'
                f'\nR2 de Teste = {round(r2_teste, 2)}'
            )
            parametos_equacao = (
                f'Coeficientes (pesos): {maquina_preditiva.coef_}'
                f'\nInterceptor (constante): {maquina_preditiva.intercept_}'
            )

            data_e_hora_pos   = datetime.now()
            tempo_restante    = (data_e_hora_pos - data_e_hora_atual) * (tot_rodadas - rodada)
            data_e_hora_atual = datetime.now()

            # Imprime no terminal apenas folds com R² treino e teste > 50%
            if r2 > 50 and r2_teste > 50:
                print(resultado)
                print(f'\n{rodada} / {tot_rodadas}  {round((rodada / tot_rodadas) * 100, 2)} %'
                      f'\nTempo restante {tempo_restante}')

            # ---------------------------------------------------------------
            # SALVAMENTO DOS RESULTADOS EM ARQUIVO
            # melhores_combinacoes_20.txt: todos os folds (sem filtro)
            # melhores_combinacoes_60.txt: folds com R² treino>50% e teste>60%
            # melhores_combinacoes.txt   : melhor resultado por fold
            # ---------------------------------------------------------------
            with open(os.path.join(pasta, 'melhores_combinacoes_20.txt'), 'a') as armazena_resultado:
                armazena_resultado.write(resultado + '\n')

            if r2_teste > 60 and r2 > 50:
                with open(os.path.join(pasta, 'melhores_combinacoes_60.txt'), 'a') as armazena_resultado:
                    armazena_resultado.write(resultado + '\n')
                    armazena_resultado.write(parametos_equacao + '\n')

            # Atualiza o melhor resultado de treino do fold atual
            if r2 > r2_treino_melhor:
                r2_treino_melhor                = r2
                melhor_parametos_equacao_treino = parametos_equacao
                melhor_resultado_treino         = (
                    f'\nRodada {rodada} / {tot_rodadas}  '
                    f'{round((rodada / tot_rodadas) * 100, 2)} %'
                    f'\nParâmetros {i}'
                    f'\nR2 de Treino = {round(r2, 2)}'
                    f'\nR2 de Teste = {round(r2_teste, 2)}'
                )

            if r2_teste > r2_teste_melhor:
                r2_teste_melhor = r2_teste

        # Salva o melhor resultado (treino) encontrado no fold atual
        with open(os.path.join(pasta, 'melhores_combinacoes.txt'), 'a') as armazena_melhor_resultado:
            armazena_melhor_resultado.write(
                f'\n********* Melhor Resultado do Fold ************'
                f'\nLinhas de Treino {linhas_treino}'
                f'\nLinhas de Teste {linhas_teste}'
                f'\n{melhor_resultado_treino}\n'
                f'{melhor_parametos_equacao_treino}\n\n'
            )

        # Reinicia os melhores para o próximo fold
        r2_treino_melhor = 0
        r2_teste_melhor  = 0


    # =============================================================================
    # ESTATÍSTICAS GERAIS (RESUMO DE TODOS OS FOLDS)
    # =============================================================================

    print("\nDiretório atual:", os.getcwd())

    print("\n===== ESTATÍSTICAS GERAIS =====")
    print("R2 Treino médio:",                  round(np.mean(lista_r2_treino),    2))
    print("R2 Treino desvio padrão:",          round(np.std(lista_r2_treino),     2))
    print("R2 Ajustado Treino médio:",         round(np.mean(lista_r2aj_treino),  2))
    print("R2 Ajustado Treino desvio padrão:", round(np.std(lista_r2aj_treino),   2))
    print("R2 Teste médio:",                   round(np.mean(lista_r2_teste),     2))
    print("R2 Teste desvio padrão:",           round(np.std(lista_r2_teste),      2))
    print("RMSE Teste médio:",                 round(np.mean(lista_rmse_teste),   2))
    print("MAE Teste médio:",                  round(np.mean(lista_mae_teste),    2))
    print("MAE Teste desvio padrão:",          round(np.std(lista_mae_teste),     2))

    # Intervalo empírico de 95% do R² de teste (percentis 2,5 e 97,5).
    # Nota: as partições do RepeatedKFold são sobrepostas, então este é um intervalo
    # empírico da distribuição dos escores, não um intervalo de confiança clássico.
    ic_lower = np.percentile(lista_r2_teste, 2.5)
    ic_upper = np.percentile(lista_r2_teste, 97.5)
    print(f"\nIntervalo empírico 95% R² Teste: [{ic_lower:.2f}% - {ic_upper:.2f}%]")

    # =============================================================================
    # VALIDAÇÃO TEMPORAL
    # Como o RepeatedKFold embaralha as observações, meses autocorrelacionados podem
    # ficar divididos entre treino e teste. Para verificar se o desempenho não depende
    # dessa proximidade temporal, o modelo é reavaliado por dois esquemas que preservam
    # a ordem cronológica:
    #   (1) GroupKFold por ANO: anos inteiros ficam fora do treino em cada partição.
    #   (2) TimeSeriesSplit (janela expansível): treina no passado e prevê o futuro.
    # Em ambos, a winsorização é estimada dentro do fold (apenas no treino).
    # =============================================================================

    print("\n===== VALIDAÇÃO TEMPORAL (ROBUSTEZ) =====")

    def _fit_eval_fold(x_tr, x_te, y_tr, y_te):
        """Ajusta o pipeline (winsor no treino + scaler + poly + Ridge/GridSearch)
        e retorna o R² de teste (%). Réplica fiel do pipeline principal."""
        y_tr = y_tr.copy()
        limite = np.percentile(y_tr.values, 3)
        y_tr = y_tr.clip(lower=limite)

        sc = StandardScaler()
        x_tr_sc = sc.fit_transform(x_tr)
        x_te_sc = sc.transform(x_te)

        pf = PolynomialFeatures(degree=GRAU_MODELO, include_bias=False)
        X_tr_p = pf.fit_transform(x_tr_sc)
        X_te_p = pf.transform(x_te_sc)

        gs = GridSearchCV(Ridge(), param_grid={'alpha': [0.1, 1.0, 10.0, 50.0, 100.0]},
                          cv=5, scoring='r2')
        gs.fit(X_tr_p, y_tr.values.ravel())
        m = gs.best_estimator_
        return r2_score(y_te, m.predict(X_te_p)) * 100

    # --- (1) GroupKFold por ANO -------------------------------------------------
    # Requer a coluna ANO na planilha; agrupa observações do mesmo ano.
    anos_grupo = dados_total.loc[dados_x.index, 'ANO'].values
    gkf = GroupKFold(n_splits=5)
    r2_gkf = []
    for tr_idx, te_idx in gkf.split(dados_x, dados_y, groups=anos_grupo):
        r2_gkf.append(_fit_eval_fold(
            dados_x.iloc[tr_idx], dados_x.iloc[te_idx],
            dados_y.iloc[tr_idx], dados_y.iloc[te_idx]))

    # --- (2) TimeSeriesSplit (janela expansível) --------------------------------
    # Pressupõe que as linhas estão em ordem cronológica (verificado nos dados).
    tss = TimeSeriesSplit(n_splits=5)
    r2_tss = []
    for tr_idx, te_idx in tss.split(dados_x):
        r2_tss.append(_fit_eval_fold(
            dados_x.iloc[tr_idx], dados_x.iloc[te_idx],
            dados_y.iloc[tr_idx], dados_y.iloc[te_idx]))

    print(f"R² teste - RepeatedKFold (principal) : {np.mean(lista_r2_teste):.2f}%")
    print(f"R² teste - GroupKFold por ano        : {np.mean(r2_gkf):.2f}%  "
          f"(folds: {', '.join(f'{v:.1f}' for v in r2_gkf)})")
    print(f"R² teste - TimeSeriesSplit expansível: {np.mean(r2_tss):.2f}%  "
          f"(folds: {', '.join(f'{v:.1f}' for v in r2_tss)})")

    gap_gkf = np.mean(lista_r2_teste) - np.mean(r2_gkf)
    gap_tss = np.mean(lista_r2_teste) - np.mean(r2_tss)
    print(f"\nQueda RKF:GroupKFold  : {gap_gkf:+.2f} pp")
    print(f"Queda RKF:TimeSeries  : {gap_tss:+.2f} pp")
    if max(gap_gkf, gap_tss) < 10:
        print("Desempenho ROBUSTO à validação temporal (quedas < 10 pp): "
              "o R² não é artefato de vazamento por autocorrelação.")
    else:
        print("Queda relevante sob validação temporal: parte do desempenho "
              "dependia da estrutura temporal. Interpretar como modelo diagnóstico.")

    # =============================================================================
    # DIAGNÓSTICO DE MULTICOLINEARIDADE - VIF (Variance Inflation Factor)
    # Calculado sobre as variáveis originais (antes da expansão polinomial).
    # Referência: VIF < 5 aceitável | 5-10 moderado | > 10 grave
    # =============================================================================

    from statsmodels.tools import add_constant
    X_vif   = dados_total[colunas_x].dropna()
    X_const = add_constant(X_vif)

    vif_data = pd.DataFrame({
        "Variável": X_vif.columns,
        "VIF": [variance_inflation_factor(X_const.values, i+1)
                for i in range(X_vif.shape[1])]
    })

    print("\n===== VIF =====")
    print(vif_data)

    # =============================================================================
    # ANÁLISE DE RESÍDUOS - MODELO AJUSTADO EM TODOS OS DADOS
    # O modelo é ajustado no conjunto completo (todos os meses da base; 297 na série 2001-2025) apenas para diagnóstico visual.
    # A avaliação preditiva real é feita pela validação cruzada acima.
    # =============================================================================


    X_full = dados_total[colunas_x].values
    y_full = dados_total[colunas_y].values.ravel()

    # Padronização e expansão polinomial no conjunto completo
    scaler_full = StandardScaler()
    X_sc        = scaler_full.fit_transform(X_full)

    poly_full = PolynomialFeatures(degree=GRAU_MODELO, include_bias=False)
    X_poly    = poly_full.fit_transform(X_sc)

    # Ajuste do modelo com alpha médio dos folds para diagnóstico de resíduos
    alpha_medio  = np.mean(alphas_escolhidos)
    print(f"Alpha médio utilizado na análise de resíduos: {alpha_medio:.4f}")
    modelo_final = Ridge(alpha=alpha_medio)
    modelo_final.fit(X_poly, y_full)

    # =============================================================================
    # EXPORTAÇÃO DOS COEFICIENTES DO MODELO FINAL
    # =============================================================================

    nomes_features = poly_full.get_feature_names_out(colunas_x)
    coef_final = pd.DataFrame({
        'feature': nomes_features,
        'coef':    modelo_final.coef_,
        'bioma':   bioma
    })
    caminho_coef = proximo_nome(pasta, f'coeficientes_ridge_{bioma}', 'csv')
    coef_final.to_csv(caminho_coef, index=False)
    print(f"Arquivo {os.path.basename(caminho_coef)} salvo com sucesso!")
    print(coef_final)

    y_pred_full = modelo_final.predict(X_poly)
    residuos    = y_full - y_pred_full

    # Teste de Shapiro-Wilk: verifica normalidade dos resíduos
    # H0: resíduos seguem distribuição normal (p > 0.05: não rejeita H0)
    stat_sw, p_sw = stats.shapiro(residuos)
    print(f"\n===== NORMALIDADE DOS RESÍDUOS =====")
    print(f"Shapiro-Wilk: stat={stat_sw:.4f}, p={p_sw:.4f}")
    print("Resíduos normais" if p_sw > 0.05 else "Resíduos NÃO normais")

    # -----------------------------------------------------------------------------
    # AUTOCORRELAÇÃO DOS RESÍDUOS (ACF + Ljung-Box)
    # Complementa a validação temporal: resíduos SEM autocorrelação sobrante
    # sustentam que o esquema de validação não está vazando por proximidade
    # temporal. ACF(lag1) alto (> ~0,3) ou Ljung-Box p < 0,05 indicam
    # autocorrelação residual relevante.
    # -----------------------------------------------------------------------------
    ljung_box_p = np.nan
    acf_lag1    = np.nan
    try:
        from statsmodels.tsa.stattools import acf as _acf
        from statsmodels.stats.diagnostic import acorr_ljungbox as _ljung
        ac = _acf(residuos, nlags=3, fft=False)
        lb = _ljung(residuos, lags=[12], return_df=True)
        ljung_box_p = lb['lb_pvalue'].values[0]
        acf_lag1    = ac[1]
        print(f"\n===== AUTOCORRELAÇÃO DOS RESÍDUOS =====")
        print(f"ACF lag1={ac[1]:+.3f} | lag2={ac[2]:+.3f} | lag3={ac[3]:+.3f}")
        print(f"Ljung-Box(12): p={lb['lb_pvalue'].values[0]:.4f}")
        if ac[1] < 0.3 and lb['lb_pvalue'].values[0] > 0.05:
            print("Sem autocorrelação residual relevante.")
        else:
            print("Autocorrelação residual presente - limitação a reportar.")
    except ImportError:
        print("(statsmodels não instalado - ACF dos resíduos não calculada)")

    # =============================================================================
    # GRÁFICOS DE DIAGNÓSTICO DOS RESÍDUOS
    # Painel 1: Resíduos vs Valores Ajustados - detecta heterocedasticidade
    # Painel 2: QQ-Plot - avalia aderência à distribuição normal
    # Painel 3: Histograma dos resíduos - visualiza a forma da distribuição
    # =============================================================================

    COR_BIOMA = cfg['cor']   # cor característica do bioma ativo
    UNID = r'(gC$\cdot$m$^{-2}\cdot$mês$^{-1}$)'

    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

    # Painel 1 - Resíduos vs Valores Ajustados
    axes[0].scatter(y_pred_full, residuos, alpha=0.55, s=24,
                    color=COR_BIOMA, edgecolor='none')
    axes[0].axhline(0, color='#333333', linestyle='--', linewidth=1.4)
    axes[0].set_xlabel(f'Valores Ajustados {UNID}')
    axes[0].set_ylabel(f'Resíduos {UNID}')
    axes[0].set_title('Resíduos vs Valores Ajustados')

    # Painel 2 - QQ-Plot (colorido por bioma)
    (osm, osr), (sl, inter, _) = stats.probplot(residuos, dist="norm")
    axes[1].scatter(osm, osr, s=22, color=COR_BIOMA, alpha=0.6, edgecolor='none')
    axes[1].plot(osm, sl*osm + inter, color='#333333', linewidth=1.7)
    axes[1].set_xlabel('Quantis teóricos')
    axes[1].set_ylabel('Valores ordenados')
    axes[1].set_title('QQ-Plot dos Resíduos')

    # Painel 3 - Histograma + caixa estatística (Shapiro, Ljung-Box, ACF)
    axes[2].hist(residuos, bins=30, color=COR_BIOMA, alpha=0.85, edgecolor='white')
    axes[2].set_xlabel(f'Resíduo {UNID}')
    axes[2].set_ylabel('Frequência')
    axes[2].set_title('Distribuição dos Resíduos')
    try:
        from statsmodels.stats.diagnostic import acorr_ljungbox as _lb2
        from statsmodels.tsa.stattools import acf as _acf2
        _lbp = _lb2(residuos, lags=[12], return_df=True)['lb_pvalue'].values[0]
        _ac1 = _acf2(residuos, nlags=1, fft=False)[1]
        _caixa = (f'Shapiro-Wilk\nW = {stat_sw:.4f}\np = {p_sw:.4f}\n'
                  f'Ljung-Box(12)\np = {_lbp:.4f}\nACF(1) = {_ac1:+.2f}')
    except Exception:
        _caixa = f'Shapiro-Wilk\nW = {stat_sw:.4f}\np = {p_sw:.4f}'
    axes[2].text(0.03, 0.97, _caixa, transform=axes[2].transAxes,
                 va='top', ha='left', fontsize=11, fontweight='bold',
                 bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                           edgecolor=COR_BIOMA, linewidth=1.7, alpha=0.95))

    # Rótulo do bioma na LATERAL ESQUERDA (nome), rotacionado e colorido,
    # no lugar do título superior. Facilita a montagem posterior do painel de 3 linhas.
    fig.text(0.012, 0.5, cfg["nome"],
             rotation=90, va='center', ha='center',
             fontsize=18, fontweight='bold', color=COR_BIOMA)

    plt.tight_layout(rect=[0.045, 0, 1, 1])
    caminho_residuos = proximo_nome(pasta, f'analise_residuos_{bioma}', 'png')
    plt.savefig(caminho_residuos, dpi=600, bbox_inches='tight', pad_inches=0.25)
    plt.close(fig)
    print(f"Gráfico salvo em {os.path.basename(caminho_residuos)} (600 dpi)")

    # =============================================================================
    # GRÁFICO OBSERVADO vs PREDITO (Figura 3)
    # Os pontos são predições OUT-OF-FOLD: cada observação é prevista por um modelo
    # que não a viu no treino (KFold de 5 partições). O pipeline sem leakage é
    # replicado em cada fold (winsorização 3% só no treino, padronização, expansão
    # polinomial e Ridge com alpha por GridSearch). Assim, os pontos plotados e o R²
    # anotado provêm do MESMO procedimento fora da amostra, sem a incoerência de
    # plotar predição in-sample e anotar um R² de validação cruzada.
    # =============================================================================

    y_pred_oof = np.full(len(y_full), np.nan)
    kf_oof     = KFold(n_splits=5, shuffle=True, random_state=42)

    for tr_idx, te_idx in kf_oof.split(dados_x):
        x_tr = dados_x.iloc[tr_idx]
        x_te = dados_x.iloc[te_idx]
        y_tr = dados_y.iloc[tr_idx].copy()

        # winsorização sem leakage (percentil 3% calculado só no treino)
        limite = np.percentile(y_tr.values, 3)
        y_tr   = y_tr.clip(lower=limite)

        sc     = StandardScaler()
        xtr_sc = sc.fit_transform(x_tr)
        xte_sc = sc.transform(x_te)

        pf    = PolynomialFeatures(degree=GRAU_MODELO, include_bias=False)
        Xtr_p = pf.fit_transform(xtr_sc)
        Xte_p = pf.transform(xte_sc)

        gs = GridSearchCV(Ridge(), param_grid={'alpha': [0.1, 1.0, 10.0, 50.0, 100.0]},
                          cv=5, scoring='r2')
        gs.fit(Xtr_p, y_tr.values.ravel())
        y_pred_oof[te_idx] = gs.best_estimator_.predict(Xte_p)

    # R² calculado sobre as predições out-of-fold (exatamente os pontos plotados)
    r2_oof = r2_score(y_full, y_pred_oof) * 100

    fig_op, ax_op = plt.subplots(figsize=(7.5, 7.5))
    ax_op.scatter(y_full, y_pred_oof, alpha=0.55, s=32,
                  color=COR_BIOMA, edgecolor='white', linewidth=0.4)

    # Linha 1:1 (predição perfeita)
    lim_min = min(y_full.min(), np.nanmin(y_pred_oof))
    lim_max = max(y_full.max(), np.nanmax(y_pred_oof))
    marg    = 0.05 * (lim_max - lim_min)
    lo, hi  = lim_min - marg, lim_max + marg
    ax_op.plot([lo, hi], [lo, hi], color='#333333', linestyle='--',
               linewidth=1.6, label='Predição perfeita (1:1)')
    ax_op.set_xlim(lo, hi); ax_op.set_ylim(lo, hi)

    ax_op.set_xlabel(f'PSN observada {UNID}', fontsize=15, fontweight='bold')
    ax_op.set_ylabel(f'PSN predita {UNID}',   fontsize=15, fontweight='bold')
    ax_op.set_title(f'Observado vs. Predito - {cfg["nome"]}',
                    fontsize=16, fontweight='bold')
    ax_op.tick_params(labelsize=12)

    # Caixa com o R² out-of-fold (mesmos pontos do gráfico)
    ax_op.text(0.04, 0.96,
               f'R² (out-of-fold) = {r2_oof:.1f}%',
               transform=ax_op.transAxes, va='top', ha='left',
               fontsize=13, fontweight='bold',
               bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                         edgecolor=COR_BIOMA, linewidth=1.7, alpha=0.95))
    ax_op.legend(loc='lower right', fontsize=12, framealpha=0.95)

    plt.tight_layout()
    caminho_obspred = proximo_nome(pasta, f'obs_vs_pred_{bioma}', 'png')
    plt.savefig(caminho_obspred, dpi=600, bbox_inches='tight', pad_inches=0.2)
    plt.close(fig_op)
    print(f"Gráfico salvo em {os.path.basename(caminho_obspred)} (600 dpi)")
    print(f"R² out-of-fold (Figura 3): {r2_oof:.2f}%  |  "
          f"R² RepeatedKFold (tabelas): {np.mean(lista_r2_teste):.2f}%")

    # =============================================================================
    # Y-RANDOMIZATION (MEDIDA FORA DA AMOSTRA)
    # Verifica se o desempenho reflete relação real entre preditores e resposta.
    # A permutação da variável resposta é avaliada pela mesma estrutura de validação
    # cruzada do teste principal (RepeatedKFold 5x30) e pelo mesmo pipeline sem
    # leakage, com uma diferença: o alpha é fixado no valor médio dos folds
    # (alpha_medio), em vez de reselecionado por GridSearch em cada fold. Modelo
    # original e modelos permutados usam esse mesmo alpha fixo, de modo que a
    # comparação é internamente consistente. Medido fora da amostra, o modelo
    # permutado se comporta como modelo nulo (R² próximo de zero ou negativo), sem
    # a inflação por sobreajuste que ocorreria numa medida in-sample.
    # Dois critérios são reportados:
    #   (1) Delta = R²_CV(original) - R²_CV(permutado médio) >= 60 pp, critério
    #       operacional adotado por Guimarães et al. (2024) para modelos QSAR.
    #   (2) p-valor empírico unilateral: proporção de permutações cujo R²_CV iguala
    #       ou supera o do modelo original; significativo quando p < 0,05.
    # =============================================================================

    yrand_diferenca = np.nan
    yrand_p_empirico = np.nan
    yrand_status = 'NAO_RODADO'

    if RODAR_YRANDOMIZATION:
        print("\n===== Y-RANDOMIZATION COM VALIDAÇÃO CRUZADA =====")

        N_PERM_YR = 100

        def _r2_cv(y_vetor):
            """R² de teste médio pela mesma estrutura de validação cruzada do teste
            principal (RepeatedKFold 5x30), replicando o pipeline sem leakage:
            winsorização no percentil 3% calculada apenas no treino de cada fold,
            padronização, expansão polinomial e Ridge com alpha fixado no valor médio
            dos folds (alpha_medio), aplicado igualmente ao modelo original e aos
            permutados."""
            y_serie = pd.Series(y_vetor, index=dados_x.index)
            rkf_yr  = RepeatedKFold(n_splits=5, n_repeats=30, random_state=42)
            r2s = []
            for tr_idx, te_idx in rkf_yr.split(dados_x):
                x_tr = dados_x.iloc[tr_idx]
                x_te = dados_x.iloc[te_idx]
                y_tr = y_serie.iloc[tr_idx].copy()
                y_te = y_serie.iloc[te_idx]

                # winsorização sem leakage (só no treino), idêntica ao modelo principal
                limite = np.percentile(y_tr.values, 3)
                y_tr   = y_tr.clip(lower=limite)

                sc     = StandardScaler()
                xtr_sc = sc.fit_transform(x_tr)
                xte_sc = sc.transform(x_te)

                pf    = PolynomialFeatures(degree=GRAU_MODELO, include_bias=False)
                Xtr_p = pf.fit_transform(xtr_sc)
                Xte_p = pf.transform(xte_sc)

                m = Ridge(alpha=alpha_medio)
                m.fit(Xtr_p, y_tr.values.ravel())
                r2s.append(r2_score(y_te, m.predict(Xte_p)) * 100)
            return np.mean(r2s)

        # R² CV do modelo original (resposta na ordem correta)
        y_orig = dados_total[colunas_y].values.ravel().copy()
        r2_cv_original = _r2_cv(y_orig)

        # R² CV de cada modelo permutado
        lista_r2_cv_perm = []
        np.random.seed(42)
        for perm in range(N_PERM_YR):
            y_perm = y_orig.copy()
            np.random.shuffle(y_perm)
            r2_cv_perm = _r2_cv(y_perm)
            lista_r2_cv_perm.append(r2_cv_perm)
            if (perm + 1) % 10 == 0:
                print(f"  Permutação {perm+1}/{N_PERM_YR} - R² CV médio = {r2_cv_perm:.2f}%")

        r2_perm_medio = np.mean(lista_r2_cv_perm)
        r2_perm_dp    = np.std(lista_r2_cv_perm)
        diferenca     = r2_cv_original - r2_perm_medio

        # p-valor empírico unilateral: fração de permutações que atingiu ou superou
        # o R² CV do modelo original (correção +1/+1 para não subestimar em amostra finita)
        n_maior_igual = np.sum(np.array(lista_r2_cv_perm) >= r2_cv_original)
        p_empirico    = (n_maior_igual + 1) / (N_PERM_YR + 1)

        print(f"\nR² CV médio - modelo original        : {r2_cv_original:.2f}%")
        print(f"R² CV médio - modelos permutados     : {r2_perm_medio:.2f}% ± {r2_perm_dp:.2f}")
        print(f"Diferença (Δ CV)                     : {diferenca:.2f} pp")
        print(f"p-valor empírico unilateral          : {p_empirico:.4f}")

        criterio_delta = diferenca >= 60
        criterio_p     = p_empirico < 0.05

        yrand_diferenca  = diferenca
        yrand_p_empirico = p_empirico
        yrand_status     = 'APROVADO' if (criterio_delta and criterio_p) else 'VERIFICAR'

        if criterio_delta:
            print(f"Critério operacional de 60 pp ATINGIDO - diferença de {diferenca:.2f} pp")
        else:
            print(f"Critério operacional de 60 pp NÃO atingido - diferença de {diferenca:.2f} pp")
        if criterio_p:
            print(f"p-valor < 0,05 - desempenho estatisticamente significativo")
        else:
            print(f"p-valor ≥ 0,05 - desempenho não significativo")

        # -----------------------------------------------------------------------------
        # GRÁFICO DA Y-RANDOMIZATION
        # Fontes ampliadas e cores vivas para leitura em impressão (600 dpi).
        # Paleta: azul-marinho vivo (modelo original), verde vivo (média permutados),
        # cinza médio (distribuição permutada). Título e valores em destaque colorido.
        # -----------------------------------------------------------------------------
        COR_ORIGINAL = '#0B5FA5'   # azul vivo - R² do modelo original
        COR_MEDIA    = '#1DB954'   # verde vivo - média dos permutados
        COR_BARRAS   = '#9AA0A6'   # cinza médio - distribuição permutada

        fig_yr, ax_yr = plt.subplots(figsize=(9.5, 6.0))

        ax_yr.hist(lista_r2_cv_perm, bins=25, color=COR_BARRAS,
                   edgecolor='white', linewidth=0.6,
                   label='R² CV modelos permutados')
        ax_yr.axvline(r2_perm_medio, color=COR_MEDIA, linestyle='--',
                      linewidth=3.0, label=f'Média permutados ({r2_perm_medio:.2f}%)')
        ax_yr.axvline(r2_cv_original, color=COR_ORIGINAL, linestyle='-',
                      linewidth=3.5, label=f'R² CV modelo original ({r2_cv_original:.2f}%)')

        # Rótulo do valor original destacado junto da linha
        ax_yr.annotate(f'{r2_cv_original:.2f}%',
                       xy=(r2_cv_original, ax_yr.get_ylim()[1]*0.92),
                       xytext=(-10, 0), textcoords='offset points',
                       ha='right', va='top', fontsize=16, fontweight='bold',
                       color=COR_ORIGINAL)

        ax_yr.set_xlabel('R² de teste por validação cruzada (%)',
                         fontsize=16, fontweight='bold')
        ax_yr.set_ylabel('Frequência', fontsize=16, fontweight='bold')
        ax_yr.tick_params(axis='both', labelsize=13)

        status   = 'APROVADO' if (criterio_delta and criterio_p) else 'VERIFICAR'
        cor_stat = '#0B5FA5' if status == 'APROVADO' else '#C0392B'

        # Título em duas linhas: nome do bioma e, abaixo, Δ / p / status coloridos
        ax_yr.set_title(f'Y-Randomization (CV) - {cfg["nome"]}',
                        fontsize=18, fontweight='bold', pad=34, color='#202124')
        ax_yr.text(0.5, 1.015,
                   f'Δ = {diferenca:.2f} pp    |    p = {p_empirico:.4f}    |    {status}',
                   transform=ax_yr.transAxes, ha='center', va='bottom',
                   fontsize=15, fontweight='bold', color=cor_stat)

        ax_yr.legend(fontsize=13, framealpha=0.95)
        plt.tight_layout()
        caminho_yr = proximo_nome(pasta, f'yrandomization_{bioma}', 'png')
        plt.savefig(caminho_yr, dpi=600, bbox_inches='tight')
        plt.close(fig_yr)
        print(f"Gráfico salvo em {os.path.basename(caminho_yr)} (600 dpi)")
    else:
        print("\n===== Y-RANDOMIZATION DESATIVADA (RODAR_YRANDOMIZATION = False) =====")
        print("Bloco de Y-randomization pulado. Defina RODAR_YRANDOMIZATION = False para rodar.")

    # =============================================================================
    # SELEÇÃO DO GRAU POLINOMIAL
    # Executado quando TESTAR_GRAUS = True. Compara os graus 1 a GRAU_MAXIMO_TESTE
    # sob a mesma validação cruzada do modelo principal (RepeatedKFold 5x30 = 150
    # folds). Critério de seleção: maior R² de teste associado ao menor gap
    # treino-teste. Guimarães et al. (2024) adotaram o grau 4 para dados moleculares;
    # para os dados ambientais deste estudo o grau é determinado empiricamente.
    # =============================================================================

    if TESTAR_GRAUS:

        print(f"\n===== SELEÇÃO DO GRAU POLINOMIAL (graus 1 a {GRAU_MAXIMO_TESTE}) =====")
        print("Aguarde - rodando 150 folds para cada grau...\n")

        from sklearn.utils import shuffle as sk_shuffle

        resultados_graus = {}

        for degree in range(1, GRAU_MAXIMO_TESTE + 1):

            lst_r2_tr = []
            lst_r2_te = []
            lst_rmse  = []
            lst_mae   = []

            rkf_grau = RepeatedKFold(n_splits=5, n_repeats=30, random_state=42)

            for tr_idx, te_idx in rkf_grau.split(dados_x):
                x_tr = dados_x.iloc[tr_idx]
                x_te = dados_x.iloc[te_idx]
                y_tr = dados_y.iloc[tr_idx]
                y_te = dados_y.iloc[te_idx]

                # Winsorização SEM leakage: percentil 3% calculado apenas no treino,
                # idêntico ao pipeline principal (garante coerência do R² de treino
                # e do gap de overfitting entre a Figura 4 e as tabelas).
                _lim_g = np.percentile(y_tr.values, 3)
                y_tr   = y_tr.clip(lower=_lim_g)

                sc   = StandardScaler()
                xtr_sc = sc.fit_transform(x_tr)
                xte_sc = sc.transform(x_te)

                pf     = PolynomialFeatures(degree=degree, include_bias=False)
                Xtr_p  = pf.fit_transform(xtr_sc)
                Xte_p  = pf.transform(xte_sc)

                rc = GridSearchCV(Ridge(),
                                  param_grid={'alpha': [0.1, 1.0, 10.0, 50.0, 100.0]},
                                  cv=5, scoring='r2')
                rc.fit(Xtr_p, y_tr.values.ravel())
                m = rc.best_estimator_
                m.fit(Xtr_p, y_tr.values.ravel())

                lst_r2_tr.append(r2_score(y_tr, m.predict(Xtr_p)) * 100)
                lst_r2_te.append(r2_score(y_te, m.predict(Xte_p)) * 100)
                lst_rmse.append(np.sqrt(mean_squared_error(y_te, m.predict(Xte_p))))
                lst_mae.append(mean_absolute_error(y_te, m.predict(Xte_p)))

            n_termos = PolynomialFeatures(degree=degree,
                                          include_bias=False).fit_transform(
                                          np.zeros((1, len(colunas_x)))).shape[1]

            resultados_graus[degree] = {
                'r2_treino' : np.mean(lst_r2_tr),
                'r2_teste'  : np.mean(lst_r2_te),
                'std_teste' : np.std(lst_r2_te),
                'overfitting': np.mean(lst_r2_tr) - np.mean(lst_r2_te),
                'rmse'      : np.mean(lst_rmse),
                'mae'       : np.mean(lst_mae),
                'ic_lower'  : np.percentile(lst_r2_te, 2.5),
                'ic_upper'  : np.percentile(lst_r2_te, 97.5),
                'n_termos'  : n_termos,
            }

            r = resultados_graus[degree]
            marca = ' <- MODELO ATUAL' if degree == GRAU_MODELO else ''
            print(f"  Grau {degree} | {n_termos:>3} termos | "
                  f"R²treino={r['r2_treino']:.1f}% | "
                  f"R²teste={r['r2_teste']:.1f}% | "
                  f"Overfit={r['overfitting']:.1f}p.p | "
                  f"RMSE={r['rmse']:.4f}{marca}")

        # Tabela resumo
        print("\n" + "="*90)
        print(f"{'Grau':>5} | {'Termos':>6} | {'R²Treino':>9} | {'R²Teste':>8} | "
              f"{'Std':>6} | {'Overfit':>9} | {'RMSE':>8} | Int. emp. 95% Teste")
        print("="*90)
        melhor_grau = max(resultados_graus,
                          key=lambda d: resultados_graus[d]['r2_teste'] -
                                        0.5 * resultados_graus[d]['overfitting'])
        for d, r in resultados_graus.items():
            marca = ' <- ÓTIMO' if d == melhor_grau else ''
            marca += ' <- ATUAL' if d == GRAU_MODELO and d != melhor_grau else ''
            print(f"{d:>5} | {r['n_termos']:>6} | {r['r2_treino']:>8.1f}% | "
                  f"{r['r2_teste']:>7.1f}% | {r['std_teste']:>5.1f}% | "
                  f"{r['overfitting']:>8.1f}p.p | {r['rmse']:>8.4f} | "
                  f"[{r['ic_lower']:.1f}% - {r['ic_upper']:.1f}%]{marca}")
        print("="*90)
        print(f"\nGrau sugerido pelo critério R²teste − 0.5×overfitting: {melhor_grau}")
        if melhor_grau == GRAU_MODELO:
            print(f"O grau atual ({GRAU_MODELO}) já é o ótimo para estes dados.")
        else:
            print(f"Considere alterar GRAU_MODELO de {GRAU_MODELO} para {melhor_grau} "
                  f"e rodar novamente.")

        # Grafico
        graus      = list(resultados_graus.keys())
        r2_trs     = [resultados_graus[d]['r2_treino']   for d in graus]
        r2_tes     = [resultados_graus[d]['r2_teste']    for d in graus]
        stds       = [resultados_graus[d]['std_teste']   for d in graus]
        overfits   = [resultados_graus[d]['overfitting'] for d in graus]
        rmses      = [resultados_graus[d]['rmse']        for d in graus]
        n_termos_l = [resultados_graus[d]['n_termos']    for d in graus]

        fig2, axes2 = plt.subplots(1, 3, figsize=(17, 5))

        # Painel 1 - R² treino vs teste
        axes2[0].plot(graus, r2_trs, 'o-', color='steelblue',
                      linewidth=2, markersize=7, label='R² Treino')
        axes2[0].errorbar(graus, r2_tes, yerr=stds, fmt='s--', color='#d62728',
                          linewidth=2, markersize=7, capsize=5,
                          label='R² Teste (±1dp)')
        axes2[0].axvline(melhor_grau, color='green', linestyle=':',
                         alpha=0.7, label=f'Grau ótimo ({melhor_grau})')
        axes2[0].axvline(GRAU_MODELO, color='orange', linestyle='--',
                         alpha=0.7, label=f'Grau atual ({GRAU_MODELO})')
        axes2[0].set_xlabel('Grau Polinomial')
        axes2[0].set_ylabel('R² (%)')
        axes2[0].set_title('R² Treino vs Teste por Grau')
        axes2[0].legend(fontsize=8)
        axes2[0].grid(alpha=0.3)
        axes2[0].set_xticks(graus)

        # Painel 2 - Gap de overfitting
        cores_ov = ['#2ca02c' if o < 10 else '#ff7f0e' if o < 20 else '#d62728'
                    for o in overfits]
        bars2 = axes2[1].bar(graus, overfits, color=cores_ov,
                             edgecolor='white', width=0.5)
        for bar, val in zip(bars2, overfits):
            axes2[1].text(bar.get_x() + bar.get_width()/2,
                          bar.get_height() + 0.3,
                          f'{val:.1f}', ha='center', va='bottom', fontsize=9)
        axes2[1].axhline(10, color='orange', linestyle='--',
                         alpha=0.7, label='Limite aceitável (10p.p.)')
        axes2[1].set_xlabel('Grau Polinomial')
        axes2[1].set_ylabel('Treino − Teste (p.p.)')
        axes2[1].set_title('Gap de Overfitting por Grau')
        axes2[1].legend(fontsize=8)
        axes2[1].grid(alpha=0.3)
        axes2[1].set_xticks(graus)

        # Painel 3 - RMSE e número de termos
        ax_twin = axes2[2].twinx()
        axes2[2].plot(graus, rmses, 'o-', color='purple',
                      linewidth=2, markersize=7, label='RMSE')
        ax_twin.bar(graus, n_termos_l, alpha=0.2, color='gray',
                    width=0.4, label='Nº termos')
        axes2[2].set_xlabel('Grau Polinomial')
        axes2[2].set_ylabel('RMSE', color='purple')
        ax_twin.set_ylabel('Nº de termos', color='gray')
        axes2[2].set_title('RMSE e Complexidade por Grau')
        axes2[2].grid(alpha=0.3)
        axes2[2].set_xticks(graus)
        l1, lb1 = axes2[2].get_legend_handles_labels()
        l2, lb2 = ax_twin.get_legend_handles_labels()
        axes2[2].legend(l1 + l2, lb1 + lb2, fontsize=8)

        plt.suptitle(
            f'Seleção do Grau Polinomial Ótimo (testados: 1 a {GRAU_MAXIMO_TESTE})',
            fontsize=12, fontweight='bold')
        plt.tight_layout()
        caminho_grau = proximo_nome(pasta, f'selecao_grau_polinomial_{bioma}', 'png')
        plt.savefig(caminho_grau, dpi=600, bbox_inches='tight')
        plt.close(fig2)
        print(f"Gráfico salvo em {os.path.basename(caminho_grau)} (600 dpi)")

    # =========================================================================
    # RESUMO DESTE BIOMA (usado na tabela comparativa final, quando 'TODOS')
    # =========================================================================
    return {
        'bioma':              bioma,
        'nome':               cfg['nome'],
        'r2_treino_medio':    round(np.mean(lista_r2_treino), 2),
        'r2_treino_dp':       round(np.std(lista_r2_treino), 2),
        'r2aj_treino_medio':  round(np.mean(lista_r2aj_treino), 2),
        'r2_teste_medio':     round(np.mean(lista_r2_teste), 2),
        'r2_teste_dp':        round(np.std(lista_r2_teste), 2),
        'ic95_lower':         round(ic_lower, 2),
        'ic95_upper':         round(ic_upper, 2),
        'gap_overfitting_pp': round(np.mean(lista_r2_treino) - np.mean(lista_r2_teste), 2),
        'rmse_teste_medio':   round(np.mean(lista_rmse_teste), 2),
        'mae_teste_medio':    round(np.mean(lista_mae_teste), 2),
        'r2_groupkfold_ano':  round(np.mean(r2_gkf), 2),
        'r2_timeseriessplit': round(np.mean(r2_tss), 2),
        'queda_groupkfold_pp':   round(gap_gkf, 2),
        'queda_timeseries_pp':   round(gap_tss, 2),
        'vif_max':            round(vif_data['VIF'].max(), 2),
        'shapiro_p':          round(p_sw, 4),
        'residuos_normais':   'Sim' if p_sw > 0.05 else 'Não',
        'ljungbox_p':         round(ljung_box_p, 4) if not np.isnan(ljung_box_p) else None,
        'acf_lag1':           round(acf_lag1, 3) if not np.isnan(acf_lag1) else None,
        'yrand_diferenca_pp': round(yrand_diferenca, 2) if not np.isnan(yrand_diferenca) else None,
        'yrand_p_empirico':   round(yrand_p_empirico, 4) if not np.isnan(yrand_p_empirico) else None,
        'yrand_status':       yrand_status,
        'arquivo_coef':       os.path.basename(caminho_coef),
        'alpha_medio':        round(float(alpha_medio), 4),
    }


# =============================================================================
# ROBUSTEZ DA SELEÇÃO E COMPLEMENTOS DE VALIDAÇÃO (Tabelas A7 a A10)
#
# Responde a quatro perguntas que a validação principal deixa em aberto:
#   (a) o grau 2 e (b) o conjunto de variáveis continuariam a ser escolhidos se a
#       própria seleção usasse esquemas que consideram a estrutura temporal
#       (GroupKFold por ano e TimeSeriesSplit)?
#   (c) quão estável é a escolha do conjunto partição a partição no RepeatedKFold
#       (frequência de vitória; diferença pareada 1ª-2ª com IC bootstrap e Wilcoxon)?
#   (d) que importância cada variável tem fora da amostra, por bloco temporal
#       (importância por permutação, Breiman 2001), independentemente da escala dos
#       coeficientes?
#   (e) o desempenho resiste a nulos que preservam a autocorrelação e o ciclo
#       sazonal da PSN (deslocamento circular; permutação de anos inteiros)?
#
# Mesmo pipeline do modelo (winsorização 3% só no treino, StandardScaler,
# PolynomialFeatures, Ridge com GridSearchCV; random_state = 42). As funções de
# tarefa recebem tudo por argumento para funcionar com multiprocessing em qualquer
# sistema operacional (spawn no Windows). Saídas: saidas_figuras/robustez/.
# =============================================================================

_ALPHAS_GRADE = [0.1, 1.0, 10.0, 50.0, 100.0]
_VARS_CANDIDATAS = ['EV', 'PRE', 'TST', 'WAI', 'BURNlog']

def _rob_colunas(bioma, combo):
    return [f'{v}_{bioma}' if v != 'BURNlog' else f'BURN_{bioma}_log' for v in combo] + ['saz_sin', 'saz_cos']

def _rob_pacote(dados_total, bioma):
    """Subconjunto pequeno da base (uma linha por mês) enviado a cada tarefa."""
    cols = ['ANO', 'MÊS', f'PSN_{bioma}', 'saz_sin', 'saz_cos'] + [f'{v}_{bioma}' for v in ('EV', 'PRE', 'TST', 'WAI')] + [f'BURN_{bioma}_log']
    return dados_total[cols].reset_index(drop=True)

def _rob_fit(x_tr, y_tr, grau=2, alpha=None):
    y_tr = y_tr.clip(lower=np.percentile(y_tr.values, 3))
    sc = StandardScaler(); pf = PolynomialFeatures(degree=grau, include_bias=False)
    Xtr = pf.fit_transform(sc.fit_transform(x_tr))
    if alpha is None:
        m = GridSearchCV(Ridge(), {'alpha': _ALPHAS_GRADE}, cv=5, scoring='r2').fit(Xtr, y_tr.values).best_estimator_
    else:
        m = Ridge(alpha=alpha).fit(Xtr, y_tr.values)
    return sc, pf, m, r2_score(y_tr, m.predict(Xtr)) * 100

def _rob_pred(sc, pf, m, x): return m.predict(pf.transform(sc.transform(x)))

def _rob_splits(esquema, x, anos):
    if esquema == 'GroupKFold':      return list(GroupKFold(n_splits=5).split(x, groups=anos))
    if esquema == 'TimeSeriesSplit': return list(TimeSeriesSplit(n_splits=5).split(x))
    return list(RepeatedKFold(n_splits=5, n_repeats=30, random_state=42).split(x))

def _rob_tarefa_grau(args):
    bioma, nome, x_sel, grau, esq, dd = args
    x, y, anos = dd[_rob_colunas(bioma, x_sel)], dd[f'PSN_{bioma}'], dd['ANO'].values; tr_l, te_l = [], []
    for tr, te in _rob_splits(esq, x, anos):
        sc, pf, m, r2tr = _rob_fit(x.iloc[tr], y.iloc[tr].copy(), grau)
        tr_l.append(r2tr); te_l.append(r2_score(y.iloc[te], _rob_pred(sc, pf, m, x.iloc[te])) * 100)
    return dict(bioma=nome, grau=grau, esquema=esq, r2_treino=np.mean(tr_l), r2_teste=np.mean(te_l), gap_pp=np.mean(tr_l) - np.mean(te_l))

def _rob_tarefa_combo(args):
    bioma, combo, esq, dd = args
    x, y, anos = dd[_rob_colunas(bioma, combo)], dd[f'PSN_{bioma}'], dd['ANO'].values; te_l = []
    for tr, te in _rob_splits(esq, x, anos):
        sc, pf, m, _ = _rob_fit(x.iloc[tr], y.iloc[tr].copy(), 2)
        te_l.append(r2_score(y.iloc[te], _rob_pred(sc, pf, m, x.iloc[te])) * 100)
    return dict(bioma=bioma, variaveis=' + '.join(combo), esquema=esq, r2_teste=np.mean(te_l), folds=te_l)

def _rob_tarefa_perm(args):
    bioma, nome, x_sel, esq, dd = args
    vars_ = list(x_sel) + ['SAZsin', 'SAZcos']; x, y, anos = dd[_rob_colunas(bioma, x_sel)], dd[f'PSN_{bioma}'], dd['ANO'].values
    rng = np.random.default_rng(42); quedas = {v: [] for v in vars_}; razao = {v: [] for v in vars_}; parcela = {v: [] for v in vars_}; r2_base = []
    for tr, te in _rob_splits(esq, x, anos):
        sc, pf, m, _ = _rob_fit(x.iloc[tr], y.iloc[tr].copy(), 2)
        x_te = x.iloc[te].reset_index(drop=True); y_te = y.iloc[te].values
        pb = _rob_pred(sc, pf, m, x_te); r2b = r2_score(y_te, pb) * 100; r2_base.append(r2b); rmse_b = np.sqrt(np.mean((y_te - pb) ** 2))
        q_bloco = {}
        for j, v in enumerate(vars_):
            qs, rz = [], []
            for _ in range(20):
                xp = x_te.copy(); xp.iloc[:, j] = rng.permutation(xp.iloc[:, j].values); pp = _rob_pred(sc, pf, m, xp)
                qs.append(r2b - r2_score(y_te, pp) * 100); rz.append(np.sqrt(np.mean((y_te - pp) ** 2)) / rmse_b)
            q_bloco[v] = np.mean(qs); quedas[v].append(np.mean(qs)); razao[v].append(np.mean(rz))
        tot = sum(max(q, 0) for q in q_bloco.values())
        for v in vars_: parcela[v].append(max(q_bloco[v], 0) / tot * 100 if tot > 0 else np.nan)
    return [dict(bioma=nome, esquema=esq, variavel=v, queda_media_pp=np.mean(quedas[v]), queda_dp_pp=np.std(quedas[v]),
                 razao_rmse_media=np.mean(razao[v]), razao_rmse_dp=np.std(razao[v]), parcela_media_pct=np.nanmean(parcela[v]),
                 parcela_dp_pct=np.nanstd(parcela[v]), r2_base=np.mean(r2_base)) for v in vars_]

def _rob_r2_cv_fixo(bioma, x_sel, alpha, dd, y_vetor):
    x, anos = dd[_rob_colunas(bioma, x_sel)], dd['ANO'].values; y = pd.Series(y_vetor, index=x.index); r2s = []
    for tr, te in _rob_splits('RepeatedKFold', x, anos):
        sc, pf, m, _ = _rob_fit(x.iloc[tr], y.iloc[tr].copy(), 2, alpha=alpha)
        r2s.append(r2_score(y.iloc[te], _rob_pred(sc, pf, m, x.iloc[te])) * 100)
    return np.mean(r2s)

def _rob_tarefa_nulo(args):
    bioma, x_sel, alpha, tipo, k, dd = args
    y = dd[f'PSN_{bioma}'].values.copy(); anos = dd['ANO'].values
    if tipo == 'shift':
        y_n = np.roll(y, k)
    else:                                   # permutação de anos inteiros (anos completos), último ano parcial fixo
        completos = [a for a in np.unique(anos) if (anos == a).sum() == 12]
        rng = np.random.default_rng(1000 + k); perm = rng.permutation(completos); y_n = y.copy()
        for a_dest, a_orig in zip(completos, perm):
            y_n[anos == a_dest] = y[anos == a_orig]
    return dict(bioma=bioma, tipo=tipo, k=int(k), r2=_rob_r2_cv_fixo(bioma, x_sel, alpha, dd, y_n))

def rodar_robustez(dados_total, config, alphas, biomas, pasta_saida, max_workers=4):
    """Executa (a)-(e) para os biomas indicados e grava CSVs + robustez.json em pasta_saida/robustez.
    alphas: dict bioma -> alfa médio dos folds do modelo principal (usado nos nulos, como no Y-randomization)."""
    OUT = os.path.join(pasta_saida, 'robustez'); os.makedirs(OUT, exist_ok=True)
    ESQ_T = ['GroupKFold', 'TimeSeriesSplit']
    nome = {b: config[b]['nome'] for b in biomas}
    x_sel = {b: [c.replace(f'_{b}', '') for c in config[b]['x'] if c not in ('saz_sin', 'saz_cos')] for b in biomas}
    dd = {b: _rob_pacote(dados_total, b) for b in biomas}
    n_meses = len(dados_total)
    print("\n===== ROBUSTEZ DA SELEÇÃO E COMPLEMENTOS DE VALIDAÇÃO (Tabelas A7 a A10) =====")
    with ProcessPoolExecutor(max_workers=max_workers) as ex:
        grau = pd.DataFrame(list(ex.map(_rob_tarefa_grau, [(b, nome[b], x_sel[b], g, e, dd[b]) for b in biomas for g in range(1, 6) for e in ESQ_T])))
        grau.to_csv(os.path.join(OUT, 'grau_temporal.csv'), index=False)
        print("\n(a) R² de teste por grau e esquema temporal:\n" + grau.round(2).to_string(index=False))
        combos = list(combinations(_VARS_CANDIDATAS, 3))
        res = list(ex.map(_rob_tarefa_combo, [(b, c, e, dd[b]) for b in biomas for c in combos for e in ESQ_T + ['RepeatedKFold']]))
        perm = pd.DataFrame([r for lst in ex.map(_rob_tarefa_perm, [(b, nome[b], x_sel[b], e, dd[b]) for b in biomas for e in ESQ_T]) for r in lst])
        perm.to_csv(os.path.join(OUT, 'permutation_importance.csv'), index=False)
        print("\n(d) Importância por permutação fora da amostra:\n" + perm.round(2).to_string(index=False))
        tarefas = ([(b, x_sel[b], alphas[b], 'shift', k, dd[b]) for b in biomas for k in range(1, n_meses)] +
                   [(b, x_sel[b], alphas[b], 'anos', k, dd[b]) for b in biomas for k in range(100)])
        nulos = pd.DataFrame(list(ex.map(_rob_tarefa_nulo, tarefas, chunksize=8)))
    # (b) posição de cada combinação por esquema
    tab = pd.DataFrame([{k: v for k, v in r.items() if k != 'folds'} for r in res])
    tab['rank'] = tab.groupby(['bioma', 'esquema'])['r2_teste'].rank(ascending=False, method='min').astype(int)
    tab['bioma'] = tab['bioma'].map(nome); tab.sort_values(['bioma', 'esquema', 'rank']).to_csv(os.path.join(OUT, 'combos_temporal.csv'), index=False)
    print("\n(b) Combinações sob os esquemas temporais (3 primeiras):\n" + tab[tab['rank'] <= 3].sort_values(['bioma', 'esquema', 'rank']).round(2).to_string(index=False))
    # (c) estabilidade partição a partição no RepeatedKFold
    est_rows, dif_rows = [], []
    for b in biomas:
        rk = [r for r in res if r['bioma'] == b and r['esquema'] == 'RepeatedKFold']
        F = np.array([r['folds'] for r in rk]); nomes = [r['variaveis'] for r in rk]
        vence = np.bincount(F.argmax(axis=0), minlength=len(nomes)) / F.shape[1] * 100
        ordem = np.argsort(-F.mean(axis=1))
        for i in ordem: est_rows.append(dict(bioma=nome[b], variaveis=nomes[i], r2_teste=F[i].mean(), freq_melhor_pct=vence[i]))
        i1, i2 = ordem[0], ordem[1]; d = F[i1] - F[i2]
        rng = np.random.default_rng(42); boots = [rng.choice(d, size=len(d), replace=True).mean() for _ in range(5000)]
        dif_rows.append(dict(bioma=nome[b], primeira=nomes[i1], segunda=nomes[i2], dif_media_pp=d.mean(), ic95_inf=np.percentile(boots, 2.5),
                             ic95_sup=np.percentile(boots, 97.5), prop_particoes_primeira_maior=(d > 0).mean() * 100, p_wilcoxon=stats.wilcoxon(d).pvalue))
    est = pd.DataFrame(est_rows); est.to_csv(os.path.join(OUT, 'estabilidade_selecao.csv'), index=False)
    dif = pd.DataFrame(dif_rows); dif.to_csv(os.path.join(OUT, 'diferenca_pareada.csv'), index=False)
    print("\n(c) Estabilidade da seleção (RepeatedKFold):\n" + est.round(2).to_string(index=False) + "\n" + dif.round(4).to_string(index=False))
    # (e) nulos temporais
    nulos.to_csv(os.path.join(OUT, 'nulos_temporais_bruto.csv'), index=False)
    orig = {b: _rob_r2_cv_fixo(b, x_sel[b], alphas[b], dd[b], dd[b][f'PSN_{b}'].values) for b in biomas}; nul_rows = []
    for b in biomas:
        for tipo, sub in [(f'Deslocamento circular (todos os {n_meses - 1})', nulos[(nulos.bioma == b) & (nulos.tipo == 'shift')]),
                          ('Deslocamento circular múltiplo de 12 meses (calendário preservado)', nulos[(nulos.bioma == b) & (nulos.tipo == 'shift') & (nulos.k % 12 == 0)]),
                          ('Permutação de anos inteiros (calendário preservado)', nulos[(nulos.bioma == b) & (nulos.tipo == 'anos')])]:
            r = sub['r2'].values
            nul_rows.append(dict(bioma=nome[b], teste=tipo, n=len(r), r2_original=orig[b], r2_nulo_media=r.mean(), r2_nulo_dp=r.std(),
                                 r2_nulo_max=r.max(), p_empirico=((r >= orig[b]).sum() + 1) / (len(r) + 1)))
    nul = pd.DataFrame(nul_rows); nul.to_csv(os.path.join(OUT, 'nulos_temporais.csv'), index=False)
    print("\n(e) Nulos que preservam a estrutura temporal:\n" + nul.round(3).to_string(index=False))
    # síntese em JSON (lida pelo script do Word)
    js = {}
    for b in biomas:
        nb = nome[b]; sel = ' + '.join(x_sel[b])
        js[b] = dict(selecionado=sel,
                     grau={e: {int(g): float(grau[(grau.bioma == nb) & (grau.esquema == e) & (grau.grau == g)]['r2_teste'].iloc[0]) for g in range(1, 6)} for e in ESQ_T},
                     melhor_grau={e: int(grau[(grau.bioma == nb) & (grau.esquema == e)].sort_values('r2_teste', ascending=False)['grau'].iloc[0]) for e in ESQ_T},
                     rank_selecionado={e: int(tab[(tab.bioma == nb) & (tab.esquema == e) & (tab.variaveis == sel)]['rank'].iloc[0]) for e in ESQ_T},
                     top3={e: tab[(tab.bioma == nb) & (tab.esquema == e)].sort_values('rank').head(3)[['variaveis', 'r2_teste']].values.tolist() for e in ESQ_T},
                     freq_melhor=float(est[(est.bioma == nb) & (est.variaveis == sel)]['freq_melhor_pct'].iloc[0]),
                     dif=dif[dif.bioma == nb].iloc[0].to_dict(),
                     perm={e: perm[(perm.bioma == nb) & (perm.esquema == e)][['variavel', 'queda_media_pp', 'queda_dp_pp', 'razao_rmse_media', 'razao_rmse_dp', 'parcela_media_pct', 'parcela_dp_pct']].values.tolist() for e in ESQ_T},
                     nulos=nul[nul.bioma == nb].to_dict('records'))
    json.dump(js, open(os.path.join(OUT, 'robustez.json'), 'w', encoding='utf-8'), ensure_ascii=False, indent=1, default=float)
    print(f"\nRobustez concluída: {OUT}")
    return js

# =============================================================================
# EXECUÇÃO — roda 1 bioma, ou os 3 em sequência, ou os 3 em paralelo
#
# O bloco "if __name__ == '__main__':" é obrigatório no Windows quando se usa
# multiprocessing/ProcessPoolExecutor — sem ele, cada processo filho tentaria
# reimportar o script e disparar o Pool de novo, entrando em loop infinito.
# =============================================================================

if __name__ == '__main__':

    resumos = []

    if len(biomas_para_rodar) > 1 and RODAR_EM_PARALELO:
        print(f"\n>>> Rodando {len(biomas_para_rodar)} biomas EM PARALELO "
              f"({', '.join(biomas_para_rodar)}) — um processo por núcleo.\n"
              f">>> Os prints de cada bioma podem aparecer intercalados no console;\n"
              f">>> isso não afeta o resultado — cada bioma calcula de forma independente.\n")

        with ProcessPoolExecutor(max_workers=len(biomas_para_rodar)) as executor:
            futuros = {
                executor.submit(rodar_modelo_bioma, bioma, dados_total_base): bioma
                for bioma in biomas_para_rodar
            }
            for futuro in as_completed(futuros):
                bioma = futuros[futuro]
                try:
                    resumos.append(futuro.result())
                    print(f">>> Bioma {bioma} concluído.")
                except Exception as e:
                    print(f">>> ERRO no bioma {bioma}: {e}")
                    raise

        ordem = {b: i for i, b in enumerate(biomas_para_rodar)}
        resumos.sort(key=lambda r: ordem[r['bioma']])

    else:
        for bioma in biomas_para_rodar:
            resumo = rodar_modelo_bioma(bioma, dados_total_base)
            resumos.append(resumo)

    # =========================================================================
    # TABELA COMPARATIVA FINAL (só faz sentido quando rodou mais de um bioma)
    # =========================================================================
    if len(resumos) > 1:
        df_resumo = pd.DataFrame(resumos)
        print("\n\n===== RESUMO COMPARATIVO — TODOS OS BIOMAS =====")
        print(df_resumo.to_string(index=False))

        caminho_resumo = proximo_nome(PASTA_SAIDA, 'resumo_geral', 'csv')
        df_resumo.to_csv(caminho_resumo, index=False)
        print(f"\nTabela comparativa salva em {os.path.basename(caminho_resumo)} "
              f"(pasta 'saidas_figuras/')")

    # =========================================================================
    # ROBUSTEZ DA SELEÇÃO E COMPLEMENTOS DE VALIDAÇÃO (Tabelas A7 a A10)
    # =========================================================================
    if RODAR_ROBUSTEZ:
        rodar_robustez(dados_total_base, config_biomas, {r['bioma']: r['alpha_medio'] for r in resumos},
                       [r['bioma'] for r in resumos], PASTA_SAIDA, max_workers=min(4, os.cpu_count() or 1))
