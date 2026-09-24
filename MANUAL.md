# Manual do Pipeline de Extração de Dados DFT (Turbomole)

## 1. Visão geral

Este manual descreve, passo a passo, o fluxo completo usado para ir de
uma geometria molecular até as figuras finais (cargas de Bader,
espectro óptico TD-DFT, IR simulado, DOS/pDOS — via Mulliken, Becke ou
SCPA —, correção de gap por GW-BSE, exciton binding, RDF por grupo
funcional, ordem de ligação/SEN): do cálculo DFT no cluster Turbomole,
passando pela extração local com os parsers Python deste repositório,
até a geração das figuras no Veusz.

É o mesmo fluxo para qualquer molécula nova — o que muda entre
moléculas são apenas os valores específicos indicados como
placeholders (`<...>`) ao longo do texto: nomes de arquivo (ex.:
`groups.json`), índices de átomo usados no `$pop mo ... dos atoms ...`
e o número de MOs/estados do TD-DFT. **Nunca cole um placeholder
literal** (ex. `<N_MOs>`) num comando — substitua sempre pelo valor
real antes de rodar.

Os scripts Python usados aqui estão documentados individualmente em
[README.md](README.md); este manual foca no *procedimento*, não na
referência de cada flag.

### 1.1. Estrutura do pipeline

| Estágio | Onde roda | O que faz |
|---|---|---|
| 0 — Preparação de geometria | local | Converte `.xyz` em `coord` (Turbomole) |
| 1 — Cálculos DFT | cluster (PBS) | SCF, Bader, TD-DFT, Hessiana/IR, população |
| 1b — DOS/PDOS (Mulliken) | cluster (PBS) | Rodadas extras de `$pop ... dos`, só quando for gerar esse gráfico |
| 1b-alt — Diagnóstico de esquema de população | cluster, cópia separada | Testa Mulliken vs. Löwdin vs. NBO quando a pDOS por canal `s`/`p`/`d` tem dip negativo |
| 1b-alt2 — PDOS por elemento via Becke | local (Multiwfn + Molden2AIM) | Solução funcional sem densidade negativa por elemento (O/C/H) |
| 1b-alt3 — PDOS por elemento×momento angular via SCPA | local (Multiwfn) | Decomposição s/p por elemento, sem o artefato de sinal do Mulliken |
| 1c — Quebra de simetria de spin | cluster, cópia separada | Diagnóstico opcional de instabilidade de spin |
| 1d — Correção quasipartícula G0W0 | cluster (GW-BSE/TDA) + local | DOS orbital KS vs. G0W0, base do exciton binding |
| 1e — Isosuperfícies HOMO/LUMO e NTO (S1) | cluster (`proper`/`dscf -proper`) + local (render) | Visualização real-space dos orbitais de fronteira e do par hole/particle do S1, p/ investigar caráter de transferência de carga |
| 1f — Espectro Raman harmônico | cluster (PBS, `egrad`+`intense`) | Intensidades Raman projetadas nos modos normais já calculados no Estágio 1 |
| 2 — Extração local | local | Parsers Python → datasets Veusz (`.dat`) |
| 3 — Figuras | local | Scripts `veusz/plot_*.py` → figuras finais |
| 4 — Comparações entre rodadas | local, sob demanda | Compara funcionais/níveis de teoria já extraídos (não gera figura) |

### 1.2. Pré-requisitos

- Acesso ao cluster onde o Turbomole está instalado, com os binários
  `ridft`, `escf`, `aoforce`, `dscf` no `PATH` e uma pasta de cálculo
  já preparada (`control`, `basis`, `mos`/orbitais iniciais etc., via
  `define`).
- O binário `bader` (código do grupo Henkelman) disponível em
  `$PATH` ou em `~/local_bin/bader`.
- Ambiente Python local com as dependências deste repositório
  instaladas — veja [README.md § Instalação](README.md#instalação).
- **Só para os Estágios 1b-alt2/1b-alt3** (PDOS via Multiwfn): binários
  `Multiwfn_noGUI` e `molden2aim.exe` compilados em `~/local_bin/`
  (nunca dentro de uma pasta de cálculo — ver a armadilha do
  `Multiwfnpath` no Estágio 1b-alt2). Citação obrigatória no texto do
  paper: Lu, T.; Chen, F. *J. Comput. Chem.* **2012**, 33, 580–592;
  Lu, T. *J. Chem. Phys.* **2024**, 161, 082503.
- **Só para as figuras específicas de Artepillin C** (estrutura 2D/3D
  anotada com SEN: `plot_bond_order_structure_artepillinC.py`,
  `plot_molecule_vesta_style.py`): ambiente conda `vasp_env` com
  `rdkit`, `ase` e `povray` — a base deste repositório (`.venv` do
  [README.md](README.md)) não os inclui.
- **Só para o Estágio 1e** (renderização de isosuperfícies):
  `matplotlib` e `scikit-image` (ver `requirements.txt`) — rodam no
  Python do sistema mesmo sem GUI/`$DISPLAY`, ao contrário de
  VMD/PyMOL/Multiwfn com GUI, que não estavam disponíveis neste
  ambiente quando o estágio foi documentado.
- **Só para o Estágio 1f** (Raman): binário `intense` do Turbomole,
  além de `egrad` — confirme que ambos estão no `PATH` junto com os
  demais (`module load turbomole`).

---

## 2. Estágio 0 — Preparação da geometria

**Objetivo:** converter a geometria inicial (tipicamente um `.xyz`
exportado de um editor molecular) para o formato `coord` exigido pelo
Turbomole (unidades atômicas/Bohr, símbolos em minúsculo, bloco
`$coord` / `$end`), antes de rodar o `define`.

```bash
python ../scripts/xyz_to_coord.py --xyz molecula.xyz --output coord
```

Por padrão o script recusa sobrescrever um `coord` já existente; use
`--force` se for intencional. Depois de gerado o `coord`, siga o
`define` interativo do Turbomole normalmente (base, funcional,
ocupação etc.) para produzir o `control` inicial do Estágio 1.

---

## 3. Estágio 1 — Cálculos DFT (cluster, via job PBS)

**Objetivo:** produzir todas as saídas brutas do Turbomole que
alimentam a extração do Estágio 2: densidade convergida, cubo para
Bader, excitações TD-DFT, Hessiana/IR e população.

**Procedimento.** Rodar em sequência — pode ser um único job PBS ou
vários, dependendo do tamanho da molécula:

```bash
# 1. SCF base (HSE06 + D3(BJ) + senex)
ridft > ridft.out

# 2. Cubo de densidade para Bader
ridft -proper > proper.out

# 3. Bader charges (fora do Turbomole)
BADER_BIN=$(command -v bader || echo ~/local_bin/bader)
$BADER_BIN td.cub > bader_run.out
tail -5 ACF.dat   # checar NUMBER OF ELECTRONS vs. esperado

# 4. Reconverge densidade apertada antes do TD-DFT
ridft > ridft_tight.out

# 5. TD-DFT (excitações -> JDOS + óptico)
escf > escf.out

# 6. Hessiana / IR
aoforce > aoforce.out
grep -i "imaginary" aoforce.out   # checar ausência de frequência imaginária

# 7. População (Mulliken + SEN/bond-order)
dscf -proper > pop.out
```

**Verificação — checagens de sanidade antes de seguir pro Estágio 2:**

- `grep -A 4 "HOMO-LUMO Separation" proper.out`
- `tail -5 ACF.dat` — elétrons integrados perto do esperado (Z total)
- `awk '$1 ~ /^[0-9]+$/ && $1 > 6 {print $3}' vibspectrum | sort -n | head -5` — sem frequência negativa

Se alguma dessas checagens falhar, resolva antes de prosseguir — os
parsers do Estágio 2 não validam a física dos resultados, só o
formato dos arquivos.

---

## 4. Estágio 1b — DOS/PDOS via Mulliken (opcional)

**Quando usar:** apenas quando for gerar a figura de DOS/pDOS
(`veusz/plot_dos.py`). São rodadas extras de população restrita a
subconjuntos de átomos.

**Pré-condições — como obter os valores dos placeholders:**

- `<N_MOs>` = número de funções de base: `grep -n "nbf(AO)" control`,
  ou `grep -A 3 "SCF-basis functions" ridft.out`.
- `<range_O>`, `<range_C>`, `<range_H>` = intervalos de índice de
  átomo na ordem em que aparecem no `$coord`:
  `grep -A <N+1> "^\$coord" coord | tail -n +2 | head -<N> | cat -n`
  (a primeira linha do `$coord` é o cabeçalho `natoms=...` e **não**
  conta como átomo).

Todos os `<...>` abaixo são placeholders — substitua pelos valores
reais da molécula antes de rodar.

```bash
# Total
sed -i 's/^\$pop.*/$pop mo 1-<N_MOs> dos/' control
dscf -proper > dos_total.out
mv dos dos_total

# Oxigênio (ou o heteroátomo relevante — ver nota abaixo)
sed -i 's/^\$pop.*/$pop mo 1-<N_MOs> dos atoms <range_O>/' control
dscf -proper > dos_oxygen.out
mv dos dos_oxygen

# Carbono
sed -i 's/^\$pop.*/$pop mo 1-<N_MOs> dos atoms <range_C>/' control
dscf -proper > dos_carbon.out
mv dos dos_carbon

# Hidrogênio
sed -i 's/^\$pop.*/$pop mo 1-<N_MOs> dos atoms <range_H>/' control
dscf -proper > dos_hydrogen.out
mv dos dos_hydrogen
```

**Depois de terminar:** restaure o `$pop` original antes de rodar o
Estágio 1 de novo (passo 7, população para Mulliken/SEN):

```bash
sed -i 's/^\$pop.*/$pop mulliken loewdin nbo paboon/' control
```

**Limitação conhecida:** com `$pop mo dos` puro (Mulliken), os canais
`DOS_s`/`DOS_p`/`DOS_d` podem ficar transitoriamente negativos em
regiões de forte hibridização (overlap cruzado distribuído de forma
arbitrária entre canais), mesmo com `TDOS` sempre correto — foi assim
que apareceu o caso do azobenzeno (`TDOS < DOS_s` em ~0.7% dos pontos,
`DOS_s` a +820 / `DOS_p` a -476 no mesmo ponto, cancelando quase
exatamente). Não é erro de cálculo. Se esse artefato aparecer, veja o
Estágio 1b-alt (diagnóstico) e os Estágios 1b-alt2/1b-alt3 (soluções
funcionais via Multiwfn).

---

## 5. Estágio 1b-alt — Diagnóstico de esquema de população da DOS (opcional)

**Quando usar:** se a pDOS por canal `l` (`DOS_s`/`DOS_p`/`DOS_d`) do
Estágio 1b tiver dips negativos suspeitos.

**Motivação:** testar se trocar o esquema de população do `$pop mo
... dos` (Mulliken → Löwdin → NBO) reduz as excursões negativas.
Löwdin (ortogonalização simétrica antes de projetar) tende a
reduzi-las; NBO eliminaria o problema por construção, mas **não está
confirmado que o Turbomole aceita `dos` combinado com `nbo`/`loewdin`
no mesmo `$pop`** — teste primeiro num caso pequeno (`total` isolado).

Fazer numa **cópia separada**, nunca na pasta de produção:

```bash
cp -r <pasta_producao> <pasta_producao>_dos_scheme_test
cd <pasta_producao>_dos_scheme_test
```

1. **Baseline Mulliken** — repete o Estágio 1b (total + carbono +
   hidrogênio + heteroátomo), salvando como `dos_*_mulliken`.
2. **Löwdin:**
   ```bash
   sed -i 's/^\$pop.*/$pop mo 1-<N_MOs> dos loewdin/' control
   dscf -proper > dos_total_loewdin.out
   grep -i "error\|abnormal\|not implement\|unknown" dos_total_loewdin.out
   mv dos dos_total_loewdin
   # repita com "loewdin atoms <range>" para carbono/hidrogênio/heteroátomo
   ```
   Se o `dscf` reclamar de keyword desconhecida, ou o arquivo `dos` não
   for regravado (confira o mtime), este build do Turbomole não aceita
   `loewdin` junto com `dos` — reporte e pule para o passo 3.
3. **NBO** (suporte incerto — teste `total` isolado primeiro):
   ```bash
   sed -i 's/^\$pop.*/$pop mo 1-<N_MOs> dos nbo/' control
   dscf -proper > dos_total_nbo.out
   grep -i "error\|abnormal\|not implement\|unknown" dos_total_nbo.out
   ```
   Se falhar (provável — `nbo` costuma ser análise separada sobre
   `pop.out`, não opção de projeção do `dos`), documente a tentativa e
   siga com Löwdin ou Mulliken.
4. **Restaure o `$pop` de produção** antes de sair da pasta de teste:
   ```bash
   sed -i 's/^\$pop.*/$pop mulliken loewdin nbo paboon/' control
   ```
5. **Compare estatisticamente** — rode `parse_dos.py` uma vez por
   esquema (`--output` distinto por rodada) e conte pontos negativos
   por canal entre `-20` e `10` eV para cada `dos_dataset_<esquema>.dat`
   (`DOS_s`, `DOS_p`, `DOS_d`, e quantos pontos têm `TDOS < DOS_s`).

Se Löwdin reduzir bem os contadores, adote-o no Estágio 1b daqui pra
frente (documente a mudança no header do dataset/figura); se não,
mantenha Mulliken e trate como limitação conhecida no texto — não é
erro de cálculo, só a decomposição por canal `l` que tem esse
artefato. Para eliminar o problema de vez (não só reduzir), use os
Estágios 1b-alt2 (por elemento) ou 1b-alt3 (por elemento × momento
angular).

---

## 6. Estágio 1b-alt2 — PDOS por elemento via Becke (Multiwfn + Molden2AIM) (opcional)

**Status:** resolve de fato a densidade negativa por elemento (O/C/H),
confirmado no artepillin C (nenhum ponto negativo na janela -20 a
10 eV). Não cobre a decomposição por canal orbital (s/p) — para isso,
ver o Estágio 1b-alt3.

**Motivação:** o `$pop mo ... dos` do Turbomole é hardcoded em
Mulliken para projeção de DOS, independente do resto do `control` (ver
Estágio 1b-alt). A solução funcional usa ferramentas externas, a
partir do arquivo Molden que o Turbomole gera (`$last step tm2molden`).

**Setup (uma vez por máquina)** — compile `Multiwfn_noGUI` e
`molden2aim.exe` direto em `~/local_bin/Multiwfn/` (nunca dentro da
pasta de um cálculo específico: o Multiwfn usa um buffer Fortran de
tamanho fixo para ler `Multiwfnpath`, e caminhos longos são truncados
*silenciosamente* — sintoma: `WARNING: "settings.ini" was found
neither in current folder nor in the path defined by "Multiwfnpath"`
mesmo com o arquivo existindo lá):

```bash
mkdir -p ~/local_bin
unzip Multiwfn_*_bin_Linux_noGUI.zip -d /tmp/multiwfn_unzip
mv /tmp/multiwfn_unzip/Multiwfn_*_bin_Linux_noGUI ~/local_bin/Multiwfn
git clone https://github.com/zorkzou/Molden2AIM.git ~/local_bin/Multiwfn/Molden2AIM
cd ~/local_bin/Multiwfn/Molden2AIM/src
gfortran -O3 edflib.f90 edflib-pbe0.f90 molden2aim.f90 -o ~/local_bin/Multiwfn/molden2aim.exe
ln -s ~/local_bin/Multiwfn/Multiwfn_noGUI ~/local_bin/Multiwfn_noGUI
ln -s ~/local_bin/Multiwfn/molden2aim.exe ~/local_bin/molden2aim.exe
```

No `~/.bashrc` (depois do bloco `conda initialize`):

```bash
export PATH="$HOME/local_bin:$PATH"
export Multiwfnpath="$HOME/local_bin/Multiwfn"
export OMP_STACKSIZE=500M
ulimit -s unlimited 2>/dev/null || ulimit -s 65536
```

Verifique com `cd /tmp && Multiwfn_noGUI < /dev/null` (banner, **sem**
o aviso de `settings.ini`) e `molden2aim.exe < /dev/null`.

**Passo 1 — corrigir o Molden do Turbomole:**

```bash
tm2molden   # se ainda não existir
sed -i '1i [Program] turbomole' molden.input
echo 'carsph=1' >> m2a.ini   # força saída em funções esféricas
./molden2aim.exe -i molden.input
```

Confirme que gerou `molden_new.molden` (**não** `molden.wfn` — o WFN
não serve para PDOS, só para partição no espaço real tipo Bader/AIM).

**Passo 2 — rodar o Multiwfn** (em sistemas grandes, >800 funções de
base, evite stack overflow):

```bash
ulimit -s unlimited
export OMP_STACKSIZE=500M
./Multiwfn_noGUI molden_new.molden
```

No menu interativo:

```
10                # Plot total DOS, PDOS, ...
-1                # Define fragments for PDOS
  1 -> 1-<range_O>   # fragmento 1 = oxigênio (ou heteroátomo)
  2 -> <range_C>     # fragmento 2 = carbono
  3 -> <range_H>     # fragmento 3 = hidrogênio
  e                  # exporta config para DOSfrag.txt (reaproveitável)
  0                  # volta
7                 # Set method for PDOS -> escolhe 4 (Becke)
                  # NÃO usa 1 (Mulliken) nem 2 (SCPA) para esta curva --
                  # o próprio Multiwfn avisa que não são robustos para
                  # MOs desocupados, exatamente onde o problema aparece
2                 # ajusta faixa de energia (ex.: -0.74 a 0.37 a.u. == -20/+10 eV)
3                 # ajusta FWHM (ex.: 0.005 a.u.)
0                 # Draw TDOS graph!
```

No menu de pós-processamento: `7` (Toggle showing PDOS curves → Yes),
depois `3` (Export curve and line data) — gera `DOS_curve.txt`
(`Energy (a.u.) | TDOS | PDOS frag.1 | frag.2 | frag.3 | frag.4-10
vazios`) e `orginfo.txt`. Sair do módulo de DOS e reentrar reseta os
fragmentos — use `-1` → `i` (import de `DOSfrag.txt`) para recarregar.

**Passo 3 — converter para dataset Veusz:**

```bash
python ../scripts/parse_becke_pdos.py --curve DOS_curve.txt \
    --output dos_becke_dataset.dat
```

**Passo 4 — plotar:**

```bash
python ~/scripts/veusz/plot_becke_pdos.py --data dos_becke_dataset.dat \
    --orginfo orginfo.txt --output becke_pdos_figure

# ou, com HOMO/LUMO já em mãos (ex. grep -A 4 "HOMO-LUMO Separation" proper.out):
python ~/scripts/veusz/plot_becke_pdos.py --data dos_becke_dataset.dat \
    --homo -5.79375 --lumo -2.10956 --output becke_pdos_figure
```

`--orginfo` deriva HOMO/LUMO sozinho (ocupado mais alto / desocupado
mais baixo pela coluna de ocupação de `orginfo.txt`) — confirmado
batendo com o `HOMO-LUMO Separation` do Turbomole a menos de
0.001 eV.

---

## 7. Estágio 1b-alt3 — PDOS por elemento × momento angular via SCPA (Multiwfn) (opcional)

**Quando usar:** quando além da decomposição por elemento (Estágio
1b-alt2) você também precisa da decomposição s/p por elemento — o que
Becke não resolve (é um método de partição de espaço real, não de
função de base). O método SCPA opera na granularidade de função de
base (via os comandos `l s`/`l p`/`l d` do Multiwfn combinados com
restrição de átomos por `cond`), e garante composições de fragmento
sempre entre 0% e 100% — ao contrário do Mulliken cru, fonte
documentada do artefato de sinal negativo (Estágio 1b).

> **Nota:** o caminho exato do menu interativo do Multiwfn para SCPA
> (equivalente ao `7 → 4` do Estágio 1b-alt2, mas escolhendo o método
> `2`/SCPA e definindo fragmentos por elemento×momento angular) ainda
> não está documentado passo a passo neste repositório — só o formato
> de entrada/saída dos scripts abaixo. Documente aqui assim que
> validado (ver também Multiwfn manual, Seção 2.4).

Dois parsers, mesmo `DOS_curve.txt` de origem:

- `parse_scpa_pdos_angular.py` — fragmentos fixos, na ordem C-s, C-p,
  O-s, O-p, H-s, H-p (a que foi usada para o artepillin C):
  ```bash
  python ../scripts/parse_scpa_pdos_angular.py --curve DOS_curve.txt \
      --output dos_scpa_angular_dataset.dat
  ```
- `parse_dos_scpa.py` — generaliza o anterior para qualquer conjunto de
  fragmentos, na ordem em que foram definidos no menu `-1 Define
  fragments` do Multiwfn:
  ```bash
  # elemento x momento angular (6 fragmentos)
  python ../scripts/parse_dos_scpa.py --curve DOS_curve.txt \
      --labels pDOS_Cs pDOS_Cp pDOS_Os pDOS_Op pDOS_Hs pDOS_Hp \
      --output dos_scpa_angular_dataset.dat

  # só momento angular, molécula inteira (3 fragmentos)
  python ../scripts/parse_dos_scpa.py --curve DOS_curve.txt \
      --labels DOS_s DOS_p DOS_d \
      --output dos_scpa_total_dataset.dat
  ```

**Atenção:** esta decomposição de 6 fragmentos não cobre funções de
polarização `d` (ou maiores) em C/O — a soma
`pDOS_Cs+pDOS_Cp+pDOS_Os+pDOS_Op+pDOS_Hs+pDOS_Hp` não bate exatamente
com `TDOS`; a diferença é a contribuição da camada `d`.

**Plotar:**

```bash
python ~/scripts/veusz/plot_scpa_pdos_angular.py \
    --data dos_scpa_angular_dataset.dat --output scpa_pdos_figure

# com suavização extra (0 desliga; comece em 0.1 eV) e HOMO/LUMO:
python ~/scripts/veusz/plot_scpa_pdos_angular.py \
    --data dos_scpa_angular_dataset.dat --sigma 0.1 \
    --orginfo orginfo.txt --output scpa_pdos_figure
```

**Cuidado com `--sigma`:** mantenha bem abaixo da largura do gap
HOMO-LUMO, ou a cauda suavizada dos picos de fronteira vaza pro gap
(confirmado no export do artepillin C: `--sigma 0.4` fez isso contra
um gap de 3.68 eV; `--sigma 0.1` manteve os picos nítidos e isolados).
Cor por elemento (O vermelho, C azul, H verde, mesma paleta do
`plot_becke_pdos.py`), estilo de linha por momento angular (sólido =
s, tracejado = p); TDOS em preto.

---

## 8. Estágio 1c — Diagnóstico de quebra de simetria de spin (opcional)

**Objetivo:** testar se a solução RHF/RKS de camada fechada é estável
(sem instabilidade de spin), comparando com uma solução tripleto
"desconvergida de volta" para singleto (broken-symmetry).

**Atenção:** fazer sempre numa **cópia separada** da pasta principal,
nunca na pasta de produção:

```bash
cp -r <pasta_producao> <pasta_producao>_broken_symmetry
cd <pasta_producao>_broken_symmetry
```

Esse diagnóstico precisa do `define` interativo — a sintaxe manual de
`$soes`/`flip` escrita direto no `control` **não funciona de forma
confiável** (já foi tentado e falhou silenciosamente mais de uma vez).

**8.1. Configurar ocupação tripleto (via `define`)**

No menu de ocupação — chega lá respondendo `n` a "DO YOU ACCEPT THIS
OCCUPATION" durante o `eht`, ou via `&` a partir do `GENERAL MENU`
numa sessão que reaproveita o `mos` existente:

```
t          # CHOOSE UHF TRIPLET OCCUPATION
*          # salva e segue
```

Confirme `#a=<N/2+1> #b=<N/2-1>` na tela (ex.: 82 alpha / 80 beta para
162 elétrons).

**8.2. Rodar o SCF do tripleto**

```bash
ridft > ridft_triplet.out
grep -i "convergence\|total energy" ridft_triplet.out | tail -5
```

**8.3. Reabrir o `define`, reaproveitando os orbitais do tripleto**

Responda `n` à pergunta de deletar data groups, e na pergunta
`DO YOU WANT TO CHANGE THESE DATA?` aceite o default `mos` existente.
**Não use o `flip`** — não funcionou de forma confiável; a ocupação é
editada diretamente no `control` no próximo passo.

**8.4. Editar manualmente `$alpha shells`/`$beta shells`**

Iguale a ocupação (forçando singleto, `Ms=0`), mantendo os orbitais
espacialmente distintos herdados do tripleto — essa é a técnica que
de fato funcionou:

```bash
grep -n "alpha shells\|beta shells" control
```

Localize as linhas de ocupação logo abaixo de cada cabeçalho (ex.:
`a       1-82   ( 1 )` e `a       1-80   ( 1 )`) e edite para igualar
os dois ao número total de elétrons ocupados no RHF original (ex.:
`1-81` nos dois):

```bash
sed -i '<linha_alpha>s/.*/ a       1-<N_ocup>                                   ( 1 )/' control
sed -i '<linha_beta>s/.*/ a       1-<N_ocup>                                   ( 1 )/' control
```

**8.5. Reconvergir e comparar com a energia RHF original**

```bash
ridft > ridft_broken_symmetry.out
grep -i "convergence\|total energy" ridft_broken_symmetry.out | tail -5
```

**Interpretação:** se a energia final voltar a ser (quase) idêntica à
energia RHF fechada original (diferença da ordem de `1e-5` a `1e-7`
Hartree), **não há instabilidade de spin** — a solução fechada é o
mínimo real. Se convergir para uma energia visivelmente mais baixa, há
uma solução de spin quebrado mais estável, e vale investigar mais.

---

## 9. Estágio 1d — Correção quasipartícula G0W0 (opcional)

**Objetivo:** obter a densidade orbital (DOS não-projetada, peso igual
por orbital) tanto no nível Kohn-Sham (HSE06) quanto corrigido por
G0W0, sobre a mesma grade de energia, para comparação direta da
abertura do gap na região de fronteira — e alimentar o exciton binding
do Estágio 4 (`--gwbse`).

**Pré-condição:** uma rodada GW-BSE/TDA do Turbomole que produza
`qpenergies.dat` (correções de quasipartícula G0W0) e o `escf_bse.out`
correspondente. O procedimento exato foi sondado interativamente no
`define` (não havia nenhum exemplo desse fluxo no repositório antes
disso) e está automatizado em `setup_gw_bse_template.py`.

**9.1. Gerar os diretórios `gw/` e `bse/`**

```bash
python ../scripts/setup_gw_bse_template.py \
    --xyz molecula.xyz --charge 0 --functional hse06 \
    --nstates 10 --output-dir gwbse
```

`--functional` aceita `hse06`/`pbe`/`pbe0` (mesma lógica de
RI/dispersão/grid dos outros scripts de setup); `--nstates` é o número
de estados excitônicos singleto calculados no BSE; `--no-fullspec`
restringe o GW à região próxima do gap (mais barato, mas insuficiente
para a DOS completa desta seção); `--scfconv`/`--grid` seguem o mesmo
padrão dos demais scripts de setup (argumentos no topo do arquivo, via
`build_arg_parser()`).

**Importante — são duas rodadas `escf` separadas, com `control`s
diferentes** (`gw/control` tem `$rigw`; `bse/control` tem `$bse`, sem
`$rigw`) — não existe uma única rodada que faça as duas coisas juntas.

**9.2. Rodar o GW primeiro** (produz `qpenergies.dat`):

```bash
cd gwbse/gw
ridft > ridft.out
escf > escf.out
```

**9.3. Copiar `qpenergies.dat` para `bse/` e rodar o BSE**:

```bash
cp qpenergies.dat ../bse/
cd ../bse
ridft > ridft.out
escf > escf_bse.out
```

O BSE só usa a correção GW se `qpenergies.dat` já estiver na pasta
*antes* do `escf` — copie primeiro, sempre.

**Checagem de sanidade** (confirma que o BSE realmente leu a correção
GW, e não caiu de volta pro TD-DFT puro): grep `escf_bse.out` por
"Dominant contributions" e compare as energias de orbital citadas ali
com as colunas G0W0 (não Kohn-Sham) de `qpenergies.dat` — devem bater.
Se baterem com o valor Kohn-Sham puro em vez do G0W0, algo deu errado
(ex.: `qpenergies.dat` não foi copiado antes do `escf`).

**Gotchas do `define` encontrados durante a sondagem** (documentados
com mais detalhe no docstring de `setup_gw_bse_template.py`):
- No submenu `gw` (dentro do `GENERAL MENU`), o comando `rigw` sozinho
  **não liga nada** — redesenha o mesmo menu em "off" silenciosamente,
  sem erro. Precisa ser `rigw on` explicitamente. Já o `bse` (dentro do
  menu `ex`) liga normalmente como toggle simples — a mesma convenção
  on/off não vale para todo submenu do `define`.
- `fullspec` (dentro do submenu do `rigw`) é o que dá QP de **todos**
  os orbitais, não só do HOMO/LUMO — sem isso a DOS completa desta
  seção não é possível, só o gap.

**Ainda não testado**: `pbe0` como funcional de referência junto com
RI-GW (só `hse06`/`pbe` foram validados de ponta a ponta), sistemas de
camada aberta (UHF), e qualquer molécula além de uma água de 3 átomos
usada para a sondagem — meça custo/escala numa molécula pequena antes
de rodar em produção.

**9.4. Extrair e plotar a DOS comparativa** (KS vs. G0W0):

```bash
python ../scripts/parse_qpenergies_dos.py --qpenergies qpenergies.dat \
    --output-ks dos_ks_dataset.dat --output-qp dos_qp_dataset.dat \
    --broadening 0.136 --npoints 4000 --emin -20 --emax 10
```

É deliberadamente uma densidade orbital não-projetada (mesma convenção
de peso igual usada em `parse_excitations.py --output-density`) — não
envolve projeção atômica, então não tem o artefato de sinal do
Mulliken.

**Plotar** (dois painéis empilhados, mesmo eixo de energia):

```bash
python ~/scripts/veusz/plot_qpenergies_dos.py --ks dos_ks_dataset.dat \
    --qp dos_qp_dataset.dat --output qpdos_figure --emin -20 --emax 10 \
    --sigma 0.01
```

`--sigma` aplica suavização gaussiana extra sobre os dados já
discretos de `qpenergies.dat` (mesma técnica do Estágio 1b-alt3).

---

## 10. Estágio 1e — Isosuperfícies HOMO/LUMO e NTO do S1 (opcional)

**Objetivo:** visualização real-space dos orbitais de fronteira
(HOMO/LUMO canônicos do estado fundamental) e, separadamente, do par
Natural Transition Orbital (NTO) hole/particle do S1 — usado para
investigar caráter de transferência de carga sugerido por outro nível
de teoria (ex. divergência HSE06 vs. CAM-B3LYP no gap/estado S1). Não
confunda os dois: HOMO/LUMO são orbitais canônicos do estado
fundamental; os NTOs são a decomposição SVD da densidade de transição
do S1 — comandos e arquivos `.cub` diferentes, não misture.

**Trabalhe sempre numa cópia da pasta de produção**, nunca na
original — o `dscf -proper`/`proper` regravam `control` (novos data
groups `$pointval`, `$ntos_occ`/`$ntos_vir`):

```bash
cp -r <pasta_producao> <pasta_producao>_orbitals
cd <pasta_producao>_orbitals
```

**10.1. HOMO/LUMO** (sintaxe confirmada, comando leve — roda no login
node sem problema):

```bash
grep -n "\$closed shells" -A 1 control   # confirma o índice do HOMO (topo da faixa ocupada)
sed -i '/^\$end/i $pointval mo <HOMO>-<LUMO> fmt=cub' control
dscf -proper > pointval_homolumo.out
ls *.cub   # confirme os nomes reais gerados (ex. 81a.cub, 82a.cub) -- não presuma
```

`grep -A 4 "HOMO-LUMO Separation" pointval_homolumo.out` mostra as
energias de cada orbital, úteis para o título da figura.

**10.2. NTO do S1 via `proper`** (Turbomole 7.8.1 já tem essa
funcionalidade — não é preciso esperar a 7.9): a partir da pasta com o
`escf.out` (ou `escf_d3.out`) já calculado, rode `proper` interativo:

```
mos          # menu "get MOs, LMOs, or NTOs"
dftnto       # NTOs ground->excited state via escf (não `ntos`, que é p/ ricc2)
1            # número do estado excitado (1 = S1)
1            # definição do NTO: 1 = renormalized excitation part (XX)
```

**Bug confirmado (7.8.1):** a opção de definição `2`
("excitation+deexcitation part (X+Y)(X+Y)") causa **SIGSEGV** em
`cc_ntos.f:129`, reprodutível. Use sempre a opção `1`; não tem
workaround documentado além de evitar a opção 2.

Isso grava `nto_occ`/`nto_vir` (coeficientes em base CAO, um par por
autovalor) e imprime a tabela de contribuição — confirme que o
`Frequency` impresso bate com a energia do estado desejado em
`escf.out`/`escf_d3.out` (garante que "estado 1" é de fato o S1) e que
há um par dominante (`%contrib` alto) antes de prosseguir; se a
excitação for multi-configuracional, mais de um par de NTO importa e a
visualização de um único par não conta a história toda.

**10.3. Gerar os cubos dos NTOs** — `$pointval nto` lê diretamente de
`nto_occ`/`nto_vir` (não precisa de conversão Molden):

```bash
sed -i 's/^\$pointval.*/$pointval nto <indice_dominante> fmt=cub/' control
dscf -proper > pointval_nto.out
ls *.cub   # ex. nto_occ_1.cub (hole), nto_vir_1.cub (particle)
```

**Se `proper` não tiver `dftnto`** (versão mais antiga que a 7.8.1):
não invente workaround — reporte e pare. A alternativa documentada é
gerar o Molden (`tm2molden`) e usar a análise de NTO nativa do Multiwfn
(módulo de análise de excitação) a partir do `molden_new.molden` já
validado no Estágio 1b-alt2.

**10.4. Renderizar as isosuperfícies:** sem GUI/`$DISPLAY` (e sem
VMD/PyMOL/Multiwfn-com-GUI instalados no cluster, só
`Multiwfn_noGUI` — cuja opção `0 Show molecular structure and view
isosurface` é um no-op nesse build, confirmado testando), use
`render_cube_isosurface.py` (matplotlib + scikit-image, 100% headless,
sem instalar nada novo):

```bash
python ../scripts/render_cube_isosurface.py --cube 81a.cub \
    --out homo.png --isovalue 0.03 --title "HOMO (orbital <HOMO>a)"
python ../scripts/render_cube_isosurface.py --cube nto_occ_1.cub \
    --out nto_hole_S1.png --isovalue 0.03 --title "NTO hole, S1"
```

O isovalor (`--isovalue`, em e$^{-1/2}$bohr$^{-3/2}$) é escolha
estética, não física — documente sempre qual valor foi usado (`0.02`–
`0.05` costuma funcionar bem; ajuste se a superfície sair vazia ou
saturada). O script lê o cubo diretamente (formato Gaussian-cube
padrão, `OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z`), desenha
lóbulo positivo/negativo em cores separadas, e desenha átomos/ligações
(cutoff 1.75 Å) só como referência estrutural.

**Se a ferramenta certa não estiver disponível** (ex. ambiente sem
VMD/PyMOL/GUI, como aconteceu aqui): não instale software novo
silenciosamente nem force um workaround sem avisar — reporte o gap e
peça a decisão (script Python local vs. entregar só os `.cub` para
renderização na máquina do usuário vs. instalar algo novo).

---

## 11. Estágio 1f — Espectro Raman harmônico (opcional)

**Objetivo:** intensidades Raman por modo normal, a partir da mesma
Hessiana/modos normais já calculados no `aoforce` do Estágio 1 — só
falta a derivada do tensor de polarizabilidade (o IR usa derivada do
dipolo, que já se tem; Raman usa polarizabilidade).

**Mecanismo confirmado (`DOC/Documentation.pdf` local, Seção 15.2, e
o script `$TURBODIR/scripts/raman`): não é automático a partir só do
`egrad`.** É um procedimento de 3 passos: `aoforce` (frequências/modos
— já existe, não refaz) → `egrad` com `$scfinstab polly` (derivadas
cartesianas de polarizabilidade estática) → **`intense`** (projeta nos
modos normais e escreve as intensidades Raman). O `vibspectrum`
reserva uma coluna `RAMAN` ao lado de `IR` mesmo antes desse fluxo, mas
só com o símbolo de regra de seleção (`YES`/`-`) — sem números até o
`intense` rodar.

**Trabalhe numa cópia da pasta de produção**, nunca na original:

```bash
cp -r <pasta_producao> <pasta_producao>_raman
cd <pasta_producao>_raman
sed -i '/^\$end/i $scfinstab polly' control
grep -nE '^\$disp3|^\$senex|gridsize' control   # confirma que bate com a produção (mesmo nível de teoria)
```

**11.1. Custo — estime antes de rodar em produção completa.** `egrad`
com `polly` é uma resposta analítica (Lagrangiano/Z-vector) sobre
todas as 3N coordenadas cartesianas — pela documentação, "computation
of polarizability derivatives at the computational cost which is only
2–3 higher than for the electronic polarizability itself" — mas isso
ainda é da mesma ordem de grandeza do `aoforce` (que já resolve CPHF
para as 3N coordenadas). Compare com os tempos já medidos na produção
(`aoforce_d3.out`, `escf_d3.out` — `tail -5` de cada, campo
`total wall-time`) antes de decidir rodar interativamente ou via job.

**11.2. Rodar via job PBS** (mesmo template do Estágio 1, `job.pbs` —
troque só o corpo do job):

```bash
# job_egrad_polly.pbs: mesmo cabeçalho #PBS do job.pbs de produção
egrad > egrad_polly.out
qsub job_egrad_polly.pbs
qstat -f <job_id>   # confirme job_state = R e o exec_host
```

Depois de terminar, confirme no log que ele de fato gerou derivadas de
polarizabilidade (procure "polarizability derivative" ou similar) antes
de seguir para o `intense`.

**11.3. Projetar nos modos normais** (`intense` — separado do
`egrad`, não pule este passo):

```bash
intense > intense_raman.out
grep -c "^\s*[0-9]" vibspectrum   # confirma que a coluna RAMAN agora tem números, não só YES/-
```

**Não rode `intense` numa pasta sem os dados de `polly` prontos** —
sem `$scfinstab dynpol`/`polly` no `control` e sem a saída do `egrad`,
ele tenta rodar mesmo assim (não há checagem de pré-condição) e só
deixa um marcador inofensivo (`$actual step intense`) em `control`,
sem gerar Raman de verdade.

**11.4. Extração e figura** — adapte `parse_ir.py` (ou escreva um
`parse_raman.py` análogo) para ler a coluna RAMAN do `vibspectrum` em
vez de IR, e reaproveite o broadening Lorentziano/Gaussiano de
`plot_ir.py` com a mesma FWHM de 10 cm⁻¹ usada no IR, para manter as
duas figuras comparáveis no mesmo paper.

**Armadilhas conhecidas:**
1. Colunas IR e RAMAN do `vibspectrum` são duas intensidades
   diferentes por modo — não some nem escale uma pra virar a outra.
2. Confirme `$disp3`/`$senex`/`gridsize`/`$scfconv` idênticos aos da
   produção do IR já publicado, ou as duas figuras deixam de ser
   diretamente comparáveis no mesmo paper.
3. Se `egrad`/`intense` reclamar de incompatibilidade com `$senex`/
   `$rik` (RIK/troca exata seminumérica nem sempre é suportada em todo
   tipo de resposta linear) — não troque o esquema de troca exata
   sozinho; reporte o erro exato e peça a decisão.

---

## 12. Estágio 2 — Extração local (parsers Python)

**Objetivo:** transformar as saídas brutas do Turbomole/Bader (e, se
aplicável, do Multiwfn/GW-BSE) em datasets no formato "descriptor" do
Veusz, prontos para plotar.

```bash
python ../scripts/parse_bader.py --acf ACF.dat --coord coord --output bader_dataset.dat

python ../scripts/parse_excitations.py --escf escf.out \
    --output-sticks excitation_sticks.dat --output-jdos jdos_dataset.dat \
    --broadening 0.4 --npoints 2000

python ../scripts/parse_ir.py --vibspectrum vibspectrum --output ir_dataset.dat

python ../scripts/parse_mulliken.py --pop pop.out --output mulliken_dataset.dat

python ../scripts/parse_bond_order.py --pop pop.out --output bond_order_dataset.dat

python ../scripts/parse_dos.py --total dos_total --carbon dos_carbon \
    --hydrogen dos_hydrogen --oxygen dos_oxygen --output dos_dataset.dat

python ../scripts/compute_rdf.py --coord coord --groups groups.json \
    --sigma 0.15 --rmax 14.0 --npoints 500 --output rdf_dataset.dat

# Cross-check (opcional, validação)
python ../scripts/compare_charges.py --bader bader_dataset.dat \
    --mulliken mulliken_dataset.dat --output charge_comparison.dat

# Médias por grupo funcional (opcional, pro texto do paper)
python ../scripts/group_average_charges.py --groups groups.json --bader bader_dataset.dat
```

`parse_dos.py` só é necessário se o Estágio 1b foi rodado. **Nota:** a
flag chama-se `--oxygen` mesmo em moléculas sem oxigênio (ex.:
azobenzeno, que usa nitrogênio) — passe o arquivo do heteroátomo real
nessa flag mesmo assim (funciona numericamente; só a coluna de saída
fica rotulada `pDOS_O` em vez de `pDOS_N` — renomeie manualmente antes
de importar no Veusz, ou use `--element4-label` no `plot_dos.py` para
rotular a *figura* corretamente sem tocar no dataset).

A PDOS por elemento vinda daqui é Mulliken — pode ter pontos negativos
(ver Estágio 1b). Para a versão sem esse artefato, gere
`dos_becke_dataset.dat` (Estágio 1b-alt2) ou
`dos_scpa_angular_dataset.dat`/`dos_scpa_total_dataset.dat` (Estágio
1b-alt3) no lugar deste.

**Se os Estágios 1b-alt2/1b-alt3/1d foram rodados**, os parsers
correspondentes também entram aqui:

```bash
# PDOS por elemento via Becke (Estágio 1b-alt2)
python ../scripts/parse_becke_pdos.py --curve DOS_curve.txt \
    --output dos_becke_dataset.dat

# PDOS por elemento x momento angular via SCPA (Estágio 1b-alt3)
python ../scripts/parse_dos_scpa.py --curve DOS_curve.txt \
    --labels pDOS_Cs pDOS_Cp pDOS_Os pDOS_Op pDOS_Hs pDOS_Hp \
    --output dos_scpa_angular_dataset.dat

# DOS orbital KS vs. G0W0 (Estágio 1d)
python ../scripts/parse_qpenergies_dos.py --qpenergies qpenergies.dat \
    --output-ks dos_ks_dataset.dat --output-qp dos_qp_dataset.dat \
    --broadening 0.136 --npoints 4000 --emin -20 --emax 10
```

Cada script aceita `--help` para a lista completa de opções (ver
também [README.md](README.md)).

---

## 13. Estágio 3 — Figuras (Veusz)

**Objetivo:** gerar as figuras finais a partir dos datasets do
Estágio 2.

```bash
python ~/scripts/veusz/plot_rdf.py --data rdf_dataset.dat --output rdf_figure \
    --panel gr_carboxylic_acid "Carboxylic Acid" black a \
    --panel gr_phenolic_ring "Phenolic Ring" red b \
    --panel gr_prenyl_1 "Prenyl 1" green c \
    --panel gr_prenyl_2 "Prenyl 2" blue d
    # --columns N controla o número de colunas da grade (default 2)

python ~/scripts/veusz/plot_dos.py --data dos_dataset.dat \
    --jdos jdos_dataset.dat --output dos_figure --emin -20 --emax 10 \
    --element4-label O --homo-lumo homo_lumo_dataset.dat

# azobenzeno (N no lugar de O), sem overlay de HOMO/LUMO:
python ~/scripts/veusz/plot_dos.py --data dos_dataset.dat \
    --jdos jdos_dataset.dat --output dos_figure --emin -20 --emax 10 \
    --element4-label N

python ~/scripts/veusz/plot_optical.py --sticks excitation_sticks.dat \
    --jdos jdos_dataset.dat --output optical_figure

python ~/scripts/veusz/plot_bader.py --data bader_dataset.dat --output bader_figure

python ~/scripts/veusz/plot_bond_order.py --data bond_order_dataset.dat \
    --labels bond_order_dataset_labels.csv --top 20 --output bond_order_figure

python ~/scripts/veusz/plot_ir.py --data ir_dataset.dat --output ir_figure \
    --fwhm 10 --npoints 3000
```

Os grupos funcionais do `--panel` do `plot_rdf.py` (nome do dataset,
título, cor, letra do painel) são passados na linha de comando — não
é mais necessário editar o script por molécula; ajuste apenas os
argumentos `--panel` para os grupos de `groups.json` da molécula em
questão. `--homo-lumo` no `plot_dos.py` é opcional — sem ele, o painel
TDOS sai sem as linhas de HOMO/LUMO e sem o rótulo de gap.

**Se os Estágios 1b-alt2/1b-alt3/1d foram rodados**, as figuras
correspondentes:

```bash
# PDOS por elemento via Becke (Estágio 1b-alt2)
python ~/scripts/veusz/plot_becke_pdos.py --data dos_becke_dataset.dat \
    --orginfo orginfo.txt --output becke_pdos_figure

# PDOS por elemento x momento angular via SCPA (Estágio 1b-alt3)
python ~/scripts/veusz/plot_scpa_pdos_angular.py \
    --data dos_scpa_angular_dataset.dat --sigma 0.1 \
    --orginfo orginfo.txt --output scpa_pdos_figure

# DOS orbital KS vs. G0W0 (Estágio 1d)
python ~/scripts/veusz/plot_qpenergies_dos.py --ks dos_ks_dataset.dat \
    --qp dos_qp_dataset.dat --output qpdos_figure --sigma 0.01
```

### Variantes específicas de molécula (Artepillin C)

Fora do fluxo genérico acima, existem versões hardcoded para a atomagem
de Artepillin C — **não reaproveite para outra molécula ou uma
estrutura renumerada** sem revisar `bond_order_dataset_labels.csv`
primeiro:

| Script | Ambiente | Função |
|---|---|---|
| `veusz/plot_bond_order_artepillinC.py` | base (`.venv`) | Como `plot_bond_order.py`, mas colorido por categoria química (anel aromático / carbonila / C=C / outros), com as ligações citadas no texto do artigo |
| `plot_bond_order_structure_artepillinC.py` | conda `vasp_env` (rdkit + matplotlib) | Figura combinada: (a) barras de SEN, (b) estrutura 2D via RDKit com as mesmas ligações coloridas direto na geometria (projeção PCA da geometria DFT relaxada, não o gerador 2D nativo do RDKit) |
| `plot_molecule_vesta_style.py` | conda `vasp_env` (ase + povray) | Renderização 3D estilo VESTA (ball-and-stick) da geometria relaxada via POV-Ray. **Gotcha:** o build `povray` do conda-forge quebra com o bloco de câmera default do ASE (sem `angle` explícito) — precisa de câmera perspectiva explícita (`location` + `angle` + `right`/`up` unitários), e os átomos precisam estar centrados antes de renderizar |

```bash
~/miniconda3/envs/vasp_env/bin/python3 plot_bond_order_structure_artepillinC.py \
    --xyz coord.xyz --labels bond_order_dataset_labels.csv --top 20 \
    --output bond_order_structure_figure

~/miniconda3/envs/vasp_env/bin/python3 plot_molecule_vesta_style.py \
    --xyz coord.xyz --labels bond_order_dataset_labels.csv \
    --output molecule_structure_vesta
```

> `veusz/plot_dos_molecules.py` também está neste diretório mas é uma
> cópia mais antiga de `plot_dos.py` (sem `--element4-label` nem
> `--homo-lumo`) — mantida por enquanto, mas prefira `plot_dos.py`
> para figuras novas.

---

## 14. Estágio 4 — Comparações entre rodadas/níveis de teoria (sob demanda)

Não fazem parte da sequência sequencial 0→3 — não geram figura, geram
tabelas comparativas (CSV) a partir de saídas já extraídas de
**múltiplas** rodadas (funcionais diferentes, ou TD-DFT vs. GW-BSE).

**`compute_exciton_binding.py`** — energia de ligação do exciton,
`E_b = E_gap(HOMO-LUMO) - E_exciton(S1)`. Para TD-DFT, `E_gap` vem do
gap HOMO-LUMO Kohn-Sham do próprio funcional (`ridft`/`dscf`, "HOMO-LUMO
gap:"); para GW-BSE/TDA, `E_gap` vem da coluna G0W0 de `qpenergies.dat`
(Estágio 1d), nunca da coluna Kohn-Sham — GW corrige o gap
monoparticula antes do BSE somar a correção do exciton por cima:

```bash
python ../scripts/compute_exciton_binding.py \
    --tddft "TD-DFT/HSE06":artepillin_C_d3_hse06/ridft_tight_d3.out:artepillin_C_d3_hse06/escf_d3.out \
    --tddft "TD-DFT/B3LYP":artepillin_C_d3_b3lyp/ridft_tight_b3lyp.out:artepillin_C_d3_b3lyp/escf_b3lyp.out \
    --gwbse "GW-BSE/TDA":artepillin_C_d3_hse06_gwbse/qpenergies.dat:artepillin_C_d3_hse06_gwbse/escf_bse.out:81:82 \
    --output exciton_binding_energy.csv
```

(os dois últimos números do `--gwbse` são os rótulos de orbital
HOMO/LUMO em `qpenergies.dat`, ex. `81:82`.)

**`compare_low_states.py`** — tabela longa dos N estados excitados mais
baixos de várias rodadas `escf.out` lado a lado (energia real
`energy_eV`, gap orbital do par dominante `orbital_gap_eV`, e
`binding_energy_eV` = a diferença entre os dois — um binding
*por estado*, calculado contra o par occ/virt dominante daquele
estado, não necessariamente HOMO-LUMO; **não confundir** com o
binding de `compute_exciton_binding.py`, que sempre usa o par
HOMO/LUMO fixo):

```bash
python ../scripts/compare_low_states.py \
    --run hse06-tddft:artepillin_C_d3_hse06/escf_d3.out \
    --run b3lyp-tddft:artepillin_C_d3_b3lyp/escf_b3lyp.out \
    --run camb3lyp-tddft:artepillin_C_d3_camb3lyp/escf_camb3lyp.out \
    --run hse06-bse-tda:artepillin_C_d3_hse06_gwbse/escf_bse.out \
    --nstates 10 \
    --output low_states_comparison.csv
```

Scripts relacionados, já existentes, fora da sequência principal (ver
[README.md](README.md) para a lista completa): `compare_functionals.py`
(excitação dominante entre funcionais), `convergence_table.py`
(convergência numérica), `ir_error_stats.py` (erro de IR vs.
experimento).

---

## 15. Solução de problemas

| Sintoma | Causa provável / correção |
|---|---|
| Figura em branco | Quase sempre um dataset apontando para o arquivo errado, ou nome de coluna que não bate entre o `.dat` e o script de plot. Confira com `head -1 arquivo.dat` (o `descriptor`). |
| `Invalid axis range` | Alguma propriedade Veusz inválida no script (ex.: `autoRange` com valor errado) — remova a linha; o eixo já auto-escala por padrão. |
| `$end` ou data group sumindo do `control` | Edições manuais (`cat >>`, `sed`) ficaram *depois* do `$end`, e não antes. Confira com `grep -n '^\$end' control` vs. `wc -l control`. |
| `WARNING: atom count mismatch` (`parse_bader.py`) | Número de átomos no `coord` diverge do número de linhas em `ACF.dat` — cubo e `coord` são de rodadas diferentes; regere o cubo (`ridft -proper`) a partir do `coord` atual. |
| `WARNING: unknown element(s)` (`parse_bader.py`) | Elemento fora da tabela `ATOMIC_NUMBER` (só cobre H–Ca) — adicione o número atômico no dicionário antes de confiar nas cargas líquidas. |
| `DOS_s`/`DOS_p`/`DOS_d` negativos, `TDOS` correto | Artefato de sinal do Mulliken (Estágio 1b), não é bug — ver Estágios 1b-alt/1b-alt2/1b-alt3. |
| `WARNING: "settings.ini" was found neither in current folder nor in the path defined by "Multiwfnpath"` (Multiwfn) | `$Multiwfnpath` truncado — o Multiwfn usa um buffer Fortran de tamanho fixo; caminhos longos (~116 caracteres já basta) são cortados silenciosamente. Instale em `~/local_bin/Multiwfn`, nunca dentro de uma pasta de cálculo. |
| Multiwfn trava/crasha em moléculas grandes (>800 funções de base) | Stack overflow — rode com `ulimit -s unlimited` e `export OMP_STACKSIZE=500M` antes de `Multiwfn_noGUI`. |
| `povray` reclama `"Viewing angle has to be smaller than 180 degrees"` (`plot_molecule_vesta_style.py`) | Bug do parser do build conda-forge com o bloco de câmera default do ASE (sem `angle` explícito) — use uma câmera perspectiva explícita com `angle` definido. |
| `proper` trava com SIGSEGV em `cc_ntos.f` (Estágio 1e, `dftnto`) | Bug confirmado na definição de NTO `2` ("(X+Y)(X+Y)") no Turbomole 7.8.1 — use sempre a definição `1` ("renormalized XX"). |
| Coluna `RAMAN` do `vibspectrum` só tem `YES`/`-`, sem números (Estágio 1f) | `intense` não rodou (ou rodou sem os dados do `egrad`+`polly` prontos) — confira `$scfinstab polly` no `control` e rode `egrad` antes do `intense`, nessa ordem. |

---

## 16. Referência rápida

Descrição de cada script, dependências e formato dos datasets/`groups.json`:
ver [README.md](README.md).
