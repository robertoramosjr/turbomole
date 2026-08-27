# Manual do Pipeline de Extração de Dados DFT (Turbomole)

## 1. Visão geral

Este manual descreve, passo a passo, o fluxo completo usado para ir de
uma geometria molecular até as figuras finais (cargas de Bader,
espectro óptico TD-DFT, IR simulado, DOS/pDOS, RDF por grupo
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
| 1b — DOS/PDOS | cluster (PBS) | Rodadas extras de `$pop ... dos`, só quando for gerar esse gráfico |
| 1c — Quebra de simetria de spin | cluster, cópia separada | Diagnóstico opcional de instabilidade de spin |
| 2 — Extração local | local | Parsers Python → datasets Veusz (`.dat`) |
| 3 — Figuras | local | Scripts `veusz/plot_*.py` → figuras finais |

### 1.2. Pré-requisitos

- Acesso ao cluster onde o Turbomole está instalado, com os binários
  `ridft`, `escf`, `aoforce`, `dscf` no `PATH` e uma pasta de cálculo
  já preparada (`control`, `basis`, `mos`/orbitais iniciais etc., via
  `define`).
- O binário `bader` (código do grupo Henkelman) disponível em
  `$PATH` ou em `~/local_bin/bader`.
- Ambiente Python local com as dependências deste repositório
  instaladas — veja [README.md § Instalação](README.md#instalação).

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

## 4. Estágio 1b — DOS/PDOS (opcional)

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

# Oxigênio
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

---

## 5. Estágio 1c — Diagnóstico de quebra de simetria de spin (opcional)

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

**5.1. Configurar ocupação tripleto (via `define`)**

No menu de ocupação — chega lá respondendo `n` a "DO YOU ACCEPT THIS
OCCUPATION" durante o `eht`, ou via `&` a partir do `GENERAL MENU`
numa sessão que reaproveita o `mos` existente:

```
t          # CHOOSE UHF TRIPLET OCCUPATION
*          # salva e segue
```

Confirme `#a=<N/2+1> #b=<N/2-1>` na tela (ex.: 82 alpha / 80 beta para
162 elétrons).

**5.2. Rodar o SCF do tripleto**

```bash
ridft > ridft_triplet.out
grep -i "convergence\|total energy" ridft_triplet.out | tail -5
```

**5.3. Reabrir o `define`, reaproveitando os orbitais do tripleto**

Responda `n` à pergunta de deletar data groups, e na pergunta
`DO YOU WANT TO CHANGE THESE DATA?` aceite o default `mos` existente.
**Não use o `flip`** — não funcionou de forma confiável; a ocupação é
editada diretamente no `control` no próximo passo.

**5.4. Editar manualmente `$alpha shells`/`$beta shells`**

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

**5.5. Reconvergir e comparar com a energia RHF original**

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

## 6. Estágio 2 — Extração local (parsers Python)

**Objetivo:** transformar as saídas brutas do Turbomole/Bader em
datasets no formato "descriptor" do Veusz, prontos para plotar.

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

`parse_dos.py` só é necessário se o Estágio 1b foi rodado. Cada script
aceita `--help` para a lista completa de opções (ver também
[README.md](README.md)).

---

## 7. Estágio 3 — Figuras (Veusz)

**Objetivo:** gerar as figuras finais a partir dos datasets do
Estágio 2.

```bash
python ~/scripts/veusz/plot_rdf.py --data rdf_dataset.dat --output rdf_figure \
    --panel gr_carboxylic_acid "Carboxylic Acid" black a \
    --panel gr_phenolic_ring "Phenolic Ring" red b \
    --panel gr_prenyl_1 "Prenyl 1" green c \
    --panel gr_prenyl_2 "Prenyl 2" blue d

python ~/scripts/veusz/plot_dos.py --data dos_dataset.dat --jdos jdos_dataset.dat \
    --output dos_figure --emin -20 --emax 10

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
questão.

> **Pendência conhecida:** falta o script `veusz/plot_optical.py`
> (figura de absorção óptica, a partir de `excitation_sticks.dat` +
> `jdos_dataset.dat` gerados no Estágio 2) — precisa ser escrito antes
> de gerar essa figura. Ver também [README.md § Avisos
> conhecidos](README.md#avisos-conhecidos).

---

## 8. Solução de problemas

| Sintoma | Causa provável / correção |
|---|---|
| Figura em branco | Quase sempre um dataset apontando para o arquivo errado, ou nome de coluna que não bate entre o `.dat` e o script de plot. Confira com `head -1 arquivo.dat` (o `descriptor`). |
| `Invalid axis range` | Alguma propriedade Veusz inválida no script (ex.: `autoRange` com valor errado) — remova a linha; o eixo já auto-escala por padrão. |
| `$end` ou data group sumindo do `control` | Edições manuais (`cat >>`, `sed`) ficaram *depois* do `$end`, e não antes. Confira com `grep -n '^\$end' control` vs. `wc -l control`. |
| `WARNING: atom count mismatch` (`parse_bader.py`) | Número de átomos no `coord` diverge do número de linhas em `ACF.dat` — cubo e `coord` são de rodadas diferentes; regere o cubo (`ridft -proper`) a partir do `coord` atual. |
| `WARNING: unknown element(s)` (`parse_bader.py`) | Elemento fora da tabela `ATOMIC_NUMBER` (só cobre H–Ca) — adicione o número atômico no dicionário antes de confiar nas cargas líquidas. |

---

## 9. Referência rápida

Descrição de cada script, dependências e formato dos datasets/`groups.json`:
ver [README.md](README.md).
