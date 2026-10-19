# Versão da dissertação da análise do ENSO

`Analise_ENSO_Anomalias_PSN_dissertacao.py` é a rotina que gerou, na dissertação, a seção 6.5 (Tabela 6, Figura 14, Tabelas A4 e A5) e as Tabelas A11 (unidade amostral: mês ou episódio como unidade do bootstrap) e A12 (estratificação pelo trimestre mais chuvoso e pelo mais seco), além da Tabela A13 (chuva acumulada).

O artigo do ENSO usa as rotinas da pasta acima: `Bootstrap_Blocos_ENSO.py` (unidade amostral) e `Estratificacao_Sazonal_ENSO.py` (estratificação sazonal). As duas versões estimam os mesmos efeitos (as mesmas diferenças médias de anomalia entre a fase ativa e os meses neutros); os intervalos de confiança e os valores de p diferem ligeiramente porque vêm de rotinas de reamostragem diferentes. Para conferir os números da dissertação, use este script; para os do artigo, os da pasta `model/`.

## Como rodar

```bash
cd model/dissertacao
python3 Analise_ENSO_Anomalias_PSN_dissertacao.py
```

Entradas (lidas de `model/`): `Dados_base_nova_2001_2025.xlsx`, `resultados_2001_2025/numeros_extra.json` e `resultados_2001_2025/sazonalidade_variaveis.csv`. Saídas, gravadas nesta pasta para não sobrescrever as do artigo: `resultados/enso_anomalias/*.csv` (entre eles `tabelaA11_unidade_amostral.csv` e `tabelaA12_estratificacao_sazonal.csv`), `resultados/enso_anomalias.json` e `figuras/fig14_compositos_enso.png`. O bootstrap usa 10.000 reamostragens com semente fixa; a execução leva menos de um minuto.

Em relação à versão original, só os caminhos de entrada e saída e dois trechos do cabeçalho de comentário foram alterados; os cálculos são os mesmos.
