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

## Duas diferenças em relação à base do artigo

Registradas para quem for comparar os números diretamente:

1. **A série aqui vai até dezembro de 2025** (300 períodos), enquanto a base
   dos manuscritos termina em setembro de 2025 (297). A climatologia mensal de
   outubro, novembro e dezembro difere, portanto, ligeiramente.
2. **A classificação ENSO tem um mês de diferença**: 76 El Niño, 73 La Niña e
   148 neutros, contra 76 / 72 / 149 nos manuscritos.

Nenhuma das duas altera as conclusões, mas elas explicam por que a anomalia da
Mata Atlântica aparece como +3,04 pontos percentuais aqui e +3,07 nos
manuscritos. Uma nova rodada restrita a setembro de 2025 e à mesma tabela de ONI
tornaria os valores diretamente citáveis.

## Arquivos

- `tst_qc_series_mensais.csv` — 900 linhas (3 biomas x 300 períodos)
- `tst_qc_resumo.csv` — indicadores por bioma, versão e recorte sazonal
- `tst_qc_enso.csv` — anomalias por fase do ENSO nas três versões
