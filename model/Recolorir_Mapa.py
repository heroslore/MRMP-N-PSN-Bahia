# -*- coding: utf-8 -*-
"""
Recolore o mapa dos biomas para a paleta acessivel usada nas demais figuras.

O mapa original usa verde #509E52, laranja #E29B3F e amarelo-palido #F3E79B.
Sob protanopia o par verde/laranja tem separacao quase nula, e o amarelo-palido
nao corresponde a cor da Caatinga nos graficos. As tres areas sao preenchimentos
chapados, de modo que a troca pode ser feita por substituicao de cor, com
tolerancia para as bordas suavizadas.

As cores novas sao versoes levemente clareadas da paleta das figuras, para que
um preenchimento de area nao fique pesado, mantendo a matiz e portanto a
identidade entre os paineis.

Saida: figuras_artigo2/mapa_biomas_acessivel.png
"""
import os
import numpy as np
from PIL import Image

BASE = os.path.dirname(os.path.abspath(__file__))
ORIG = os.path.abspath(os.path.join(BASE, '..', 'manuscript', 'figs', 'fig1_map.png'))
SAIDA = os.path.join(BASE, 'figuras_artigo2', 'mapa_biomas_acessivel.png')

TROCAS = {
    (0x50, 0x9E, 0x52): (0x33, 0xB2, 0x95),   # Mata Atlântica  -> verde-azulado
    (0xE2, 0x9B, 0x3F): (0xEB, 0xB2, 0x33),   # Cerrado         -> âmbar
    (0xF3, 0xE7, 0x9B): (0xD6, 0x5C, 0x41),   # Caatinga        -> vermelho
}
TOL = 25.0

im = np.asarray(Image.open(ORIG).convert('RGB')).astype(np.int16)
out = im.copy()
for antiga, nova in TROCAS.items():
    dist = np.sqrt(((im - np.array(antiga)) ** 2).sum(axis=2))
    m = dist < TOL
    out[m] = nova
    print(f'#{antiga[0]:02X}{antiga[1]:02X}{antiga[2]:02X} -> '
          f'#{nova[0]:02X}{nova[1]:02X}{nova[2]:02X}: {int(m.sum()):>7} pixels')

Image.fromarray(out.astype(np.uint8)).save(SAIDA)
print('salvo em', SAIDA)
