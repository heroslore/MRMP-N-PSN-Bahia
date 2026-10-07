# Sensibilidade da TST ao filtro de qualidade QC_Day

Reprocessamento do MOD11A2 (Coleção 6.1) sob três critérios de QC_Day, para
avaliar a limitação declarada nos dois manuscritos. **A série da base principal
não foi alterada**; estes arquivos são apenas a análise de sensibilidade.

Os três critérios:

1. **sem filtro** — como na base atual: apenas pixels sem leitura são excluídos;
2. **erro ≤ 3 K** — QC_Day bits 6–7 ≠ `11`;
3. **boa qualidade** — critério restritivo sobre a qualidade mandatória.

## Resultados

**O critério intermediário é redundante nesta região.** Descartou 0% dos pixels
em todos os biomas e em todos os períodos: na Coleção 6.1 os pixels com nuvem ou
erro acima de 3 K já saem como valor de preenchimento. A série da base principal
corresponde, na prática, a esse critério.

**O critério restritivo troca qualidade por representatividade.** Descarta em
média 71,8% dos pixels da Mata Atlântica (máximo 94,8%), 36,3% da Caatinga e
19,7% do Cerrado, e aquece as séries em 1,36, 0,60 e 0,23 °C, porque os pixels
removidos são sistematicamente mais frios. O aquecimento concentra-se no
trimestre seco: 2,02 °C na Mata Atlântica e 1,23 °C na Caatinga. A maior
diferença pontual é 4,15 °C, em outubro de 2011 na Mata Atlântica. A correlação
de Pearson entre as séries nunca cai abaixo de 0,97.

**Dentro de cada trimestre, meses mais quentes perdem menos pixels**: o rho de
Spearman entre a TST sem filtro e o percentual descartado vai de -0,23 a -0,69.
No ano todo esse sinal desaparece (rho de +0,08 na Mata Atlântica), porque a
sazonalidade domina.

**As conclusões sobre o ENSO não mudam.** A anomalia térmica da Mata Atlântica
em El Niño passa de +0,888 para +0,776 °C e permanece significativa
(p = 0,0007); a frequência de meses no decil mais quente passa de 25,0% para
22,4%. Na Caatinga, a anomalia de La Niña passa de -0,717 para -0,690 °C
(p de 0,012 para 0,005). Nenhum efeito troca de sinal ou de significância.

## Recálculo sobre a base exata dos artigos

A rodada de extração cobriu 300 períodos (até dezembro de 2025) e usou uma
classificação ENSO própria, com 73 meses de La Niña. A base dos manuscritos
termina em setembro de 2025 (297 períodos), porque a precipitação do IMERG V07
foi encerrada nesse mês, e sua classificação tem 76 / 149 / 72.

`Recalcular_297.py` trunca as três versões em setembro de 2025 e refaz todos os
indicadores com a climatologia e a classificação ENSO da própria base. Os
resultados em `tst_qc_resumo_297.csv` e `tst_qc_enso_297.csv` são os citados nos
manuscritos.

**Validação do reprocessamento.** A versão sem filtro reproduz a TST da base dos
artigos com diferença máxima de 0,0005 °C e correlação de 1,000000 nos três
biomas, o que confirma que a extração seguiu o mesmo protocolo. A anomalia da
Mata Atlântica em El Niño sai como +3,07 pontos percentuais, idêntica à da
Tabela 2 do artigo do ENSO.

**Uma qualificação que só aparece no recálculo.** Sob o teste de posição com
correção de Benjamini-Hochberg, a anomalia térmica da Mata Atlântica em El Niño
sobrevive nas três versões (q = 0,0001, 0,0001 e 0,0037). Sob o bootstrap em
blocos por episódio, que é o padrão do artigo do ENSO, ela sobrevive nas duas
primeiras (p = 0,006 e 0,007) mas não na restritiva (p = 0,071), cujo intervalo
passa a tocar o zero. O efeito não troca de sinal nem de magnitude, mas a
afirmação de que resiste a todos os critérios vale para o teste de posição, não
para a inferência por episódio.

## Arquivos

- `tst_qc_series_mensais.csv` — 900 linhas (3 biomas x 300 períodos)
- `tst_qc_resumo.csv` — indicadores por bioma, versão e recorte sazonal
- `tst_qc_enso.csv` — anomalias por fase do ENSO nas três versões
