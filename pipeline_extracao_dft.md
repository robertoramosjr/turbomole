# Pipeline de Extração de Dados — Checklist

Repita esse fluxo pra qualquer molécula nova. Ajuste só: nomes de
arquivo específicos da molécula (`groups.json`, índices de átomo no
`$pop mo ... dos atoms ...`), e o número de MOs/estados no TD-DFT.

---

## Estágio 1 — Cálculos DFT (cluster, via job PBS)

Rodar em sequência (pode ser um único job ou vários, dependendo do
tamanho da molécula):

```bash
# 1. SCF base (HSE06 + D3(BJ) + senex)
ridft > ridft.out

# 2. Cubo de densidade para Bader
ridft -proper > proper.out

# 3. Bader charges (fora do Turbomole)
BADER_BIN=$(command -v bader || echo ~/local_bin/bader)
$BADER_BIN td.cub > bader_run.out
tail -5 ACF.dat   # checar NUMBER OF ELECTRONS vs. esperado

# 3b. Remove o $pointval do control -- ridft -proper (passo 2) o deixa
#     lá com o grid explícito (grid1/grid2/grid3/origin) gravado. Se
#     não for removido, TODA chamada futura de `dscf -proper` (passo 7,
#     e cada rodada do Estágio 1b) recalcula a densidade nesse grid
#     inteiro (centenas de milhares a milhões de pontos) -- trava por
#     minutos em vez de rodar em segundos. Já aconteceu (azo_trans).
grep -n '^\$pointval' control   # confirma que existe antes de mexer
sed -i '/^\$pointval/,/^origin/d' control

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

**Checagens de sanidade antes de seguir pro Estágio 2:**
- `grep -A 4 "HOMO-LUMO Separation" proper.out`
- `tail -5 ACF.dat` — elétrons integrados perto do esperado (Z total)
- `awk '$1 ~ /^[0-9]+$/ && $1 > 6 {print $3}' vibspectrum | sort -n | head -5` — sem negativo

---

## Estágio 1b — DOS/PDOS (rodadas extras, só quando for gerar esse gráfico)

`<N_MOs>` = número de funções de base (`grep -n "nbf(AO)" control`, ou
`grep -A 3 "SCF-basis functions" ridft.out`). Os ranges de átomo
(`<range_4>` para o 4º elemento — O, N, etc., conforme a molécula —,
`<range_C>`, `<range_H>`) vêm da ordem no `$coord`
(`grep -A <N+1> "^\$coord" coord | tail -n +2 | head -<N> | cat -n` —
a primeira linha do `$coord` é o cabeçalho `natoms=...`, não conta
como átomo). **Todos são placeholders — substitua pelos valores
reais antes de rodar**, nunca cole `<...>` literal

**Atenção:** o `dscf` regrava o `control` ao final e *remove* o data
group `$pop` (ele é consumido, não é persistente). Isso significa que
um `sed -i 's/^\$pop.*/.../' control` simples só funciona na
*primeira* substituição depois de um `define`/rodada anterior que
ainda tenha deixado `$pop` no arquivo — a partir da segunda chamada
não há mais linha `$pop` pra casar, o `sed` não erra (só não substitui
nada) e o `dscf` seguinte roda **sem nenhuma análise de população**,
sem gerar o arquivo `dos` e sem avisar. Sempre usar a forma abaixo, que
insere a linha quando ela não existe:

```bash
set_pop() {
  if grep -q '^\$pop' control; then
    sed -i "s/^\$pop.*/\$pop $1/" control
  else
    sed -i "/^\$end/i \$pop $1" control
  fi
}

# Total
set_pop "mo 1-<N_MOs> dos"
dscf -proper > dos_total.out
mv dos dos_total

# 4º elemento (O, N, ... conforme a molécula)
set_pop "mo 1-<N_MOs> dos atoms <range_4>"
dscf -proper > dos_element4.out
mv dos dos_element4

# Carbono
set_pop "mo 1-<N_MOs> dos atoms <range_C>"
dscf -proper > dos_carbon.out
mv dos dos_carbon

# Hidrogênio
set_pop "mo 1-<N_MOs> dos atoms <range_H>"
dscf -proper > dos_hydrogen.out
mv dos dos_hydrogen

# Restaura $pop mulliken (pro Estágio 1, passo 7, se for rodar de novo)
set_pop "mulliken loewdin nbo paboon"
```

**Checagem depois de cada rodada:** confirma que o arquivo `dos` foi
mesmo gerado (`ls -la dos`) antes do `mv` — se `$pop` não tiver sido
aplicado, o `dscf` roda normalmente (SCF + `-proper`) mas não cria
`dos`, e o `mv dos dos_total` vai falhar (ou, pior, mover um `dos`
antigo de sobra).

---

## Estágio 1c — Diagnóstico de quebra de simetria de spin (broken-symmetry, opcional)

Teste pra verificar se a solução RHF/RKS fechada é estável (sem
instabilidade de spin). Fazer numa **cópia separada** da pasta
principal, nunca na de produção. Precisa de `define` interativo (a
sintaxe manual de `$soes`/`flip` não funciona bem escrita direto no
`control` — já tentamos e falhou silenciosamente mais de uma vez).

```bash
cp -r <pasta_producao> <pasta_producao>_broken_symmetry
cd <pasta_producao>_broken_symmetry
```

**1. Configurar ocupação tripleto (via `define`)** — no menu de
ocupação (chega lá respondendo `n` a "DO YOU ACCEPT THIS OCCUPATION"
durante o `eht`, ou via `&` a partir do `GENERAL MENU` numa sessão
que reaproveita o `mos` existente):

```
t          # CHOOSE UHF TRIPLET OCCUPATION
*          # salva e segue
```

Confirma `#a=<N/2+1> #b=<N/2-1>` na tela (ex: 82 alpha / 80 beta pra
162 elétrons).

**2. Rodar o SCF do tripleto:**

```bash
ridft > ridft_triplet.out
grep -i "convergence\|total energy" ridft_triplet.out | tail -5
```

**3. Reabrir o `define`, reaproveitando os orbitais do tripleto**
(responde `n` à pergunta de deletar data groups, e na pergunta
`DO YOU WANT TO CHANGE THESE DATA?` aceita o default `mos` existente
— **não precisa do `flip`**, que não funcionou de forma confiável;
vamos editar a ocupação direto no `control` depois).

**4. Editar manualmente `$alpha shells`/`$beta shells`** pra igualar
a ocupação (forçando singleto, `Ms=0`), mantendo os orbitais
espacialmente distintos herdados do tripleto — essa é a técnica que
realmente funcionou, não o comando `flip` interativo:

```bash
grep -n "alpha shells\|beta shells" control
```

Localiza as linhas de ocupação logo abaixo de cada cabeçalho (ex:
`a       1-82   ( 1 )` e `a       1-80   ( 1 )`) e edita pra igualar
os dois ao número total de elétrons ocupados no RHF original
(ex: `1-81` nos dois):

```bash
sed -i '<linha_alpha>s/.*/ a       1-<N_ocup>                                   ( 1 )/' control
sed -i '<linha_beta>s/.*/ a       1-<N_ocup>                                   ( 1 )/' control
```

**5. Reconverge e compara com a energia RHF original:**

```bash
ridft > ridft_broken_symmetry.out
grep -i "convergence\|total energy" ridft_broken_symmetry.out | tail -5
```

Se a energia final voltar a ser (quase) idêntica à energia RHF
fechada original (diferença na ordem de `1e-5` a `1e-7` Hartree),
**não há instabilidade de spin** — a solução fechada é o mínimo
real. Se convergir pra uma energia visivelmente mais baixa, há uma
solução de spin quebrado mais estável, e vale investigar mais.

---

## Estágio 1d — GW / GW-BSE (opcional)

Correção de quasipartícula G0W0 (`qpenergies.dat`) e, em cima dela,
excitons via Bethe-Salpeter (`escf_bse.out`). São **duas rodadas
`escf` separadas, com `control`s diferentes** — não dá para combinar
num único `define`/`escf`. Setup automatizado (sondado via `define`
porque não havia nenhum exemplo desse fluxo no repositório):

```bash
python ../scripts/setup_gw_bse_template.py \
    --xyz molecula.xyz --charge 0 --functional hse06 \
    --nstates 10 --output-dir gwbse

# 1. G0W0 (RI-GW, fullspec = todos os orbitais, não só o gap)
cd gwbse/gw
ridft > ridft.out
escf > escf.out          # gera qpenergies.dat

# 2. BSE/TDA sobre os orbitais já corrigidos por GW
cp qpenergies.dat ../bse/
cd ../bse
ridft > ridft.out
escf > escf_bse.out
```

**Checagem de sanidade**: nas "Dominant contributions" de
`escf_bse.out`, as energias de orbital citadas devem bater com as
colunas G0W0 (não Kohn-Sham) de `qpenergies.dat` — se baterem com o
valor DFT puro, o BSE não pegou a correção GW.

**Ainda não testado**: `pbe0` como funcional de referência junto com
RI-GW, sistemas de camada aberta (UHF), e qualquer molécula real (só
validado numa água de 3 átomos). Ver o docstring de
`setup_gw_bse_template.py` para os detalhes de cada quirk do `define`
(ex.: `rigw` sozinho não liga nada, precisa ser `rigw on`; já `bse`
sozinho liga normalmente — inconsistência do próprio `define`, não bug
do script).

Ver Estágio 1d do `MANUAL.md` para os parsers/plots
(`parse_qpenergies_dos.py`, `plot_qpenergies_dos.py`) que consomem
esse `qpenergies.dat`.

---

## Estágio 2 — Extração local (parsers Python)

```bash
python ../scripts/parse_bader.py --acf ACF.dat --coord coord --output bader_dataset.dat

# HOMO-LUMO gap -- antes só aparecia no terminal/.out do passo 2
# (proper.out), sem registro estruturado nenhum
python ../scripts/parse_homo_lumo.py --proper proper.out --output homo_lumo_dataset.dat

python ../scripts/parse_excitations.py --escf escf.out \
    --output-sticks excitation_sticks.dat --output-jdos jdos_dataset.dat \
    --broadening 0.4 --npoints 2000

python ../scripts/parse_ir.py --vibspectrum vibspectrum --output ir_dataset.dat

python ../scripts/parse_mulliken.py --pop pop.out --output mulliken_dataset.dat

python ../scripts/parse_bond_order.py --pop pop.out --output bond_order_dataset.dat

python ../scripts/parse_dos.py --total dos_total --carbon dos_carbon \
    --hydrogen dos_hydrogen --element4 O dos_oxygen --output dos_dataset.dat
# azobenzeno (N no lugar de O):
#   --element4 N dos_element4

python ../scripts/compute_rdf.py --coord coord --groups groups.json \
    --sigma 0.15 --rmax 14.0 --npoints 500 --output rdf_dataset.dat

# Cross-check (opcional, validação)
python ../scripts/compare_charges.py --bader bader_dataset.dat \
    --mulliken mulliken_dataset.dat --output charge_comparison.dat

# Médias por grupo funcional (opcional, pro texto do paper)
python ../scripts/group_average_charges.py --groups groups.json --bader bader_dataset.dat
```

---

## Estágio 3 — Plots (Veusz)

```bash
python ~/scripts/veusz/plot_rdf.py --data rdf_dataset.dat --output rdf_figure

python ~/scripts/veusz/plot_dos.py --data dos_dataset.dat --jdos jdos_dataset.dat \
    --output dos_figure --emin -20 --emax 10 --element4-label O \
    --homo-lumo homo_lumo_dataset.dat
# azobenzeno: --element4-label N

python ~/scripts/veusz/plot_optical.py --sticks excitation_sticks.dat \
    --jdos jdos_dataset.dat --output optical_figure

python ~/scripts/veusz/plot_bader.py --data bader_dataset.dat --output bader_figure

python ~/scripts/veusz/plot_bond_order.py --data bond_order_dataset.dat \
    --labels bond_order_dataset_labels.csv --top 20 --output bond_order_figure

python ~/scripts/veusz/plot_ir.py --data ir_dataset.dat --output ir_figure \
    --fwhm 10 --npoints 3000
```

**Lembrete**: cada molécula tem sua própria lista `PANELS` (nomes de
grupo) dentro do `plot_rdf.py` — confirma que está usando a versão
certa do script antes de rodar (mantenha arquivos separados por
molécula, tipo `plot_rdf_<molecula>.py`, pra evitar rodar com a lista
errada).

---

## Se algo der errado

- **Figura em branco**: quase sempre é dataset apontando pro arquivo
  errado, ou nome de coluna não bate entre o `.dat` e o script de
  plot (`head -1 arquivo.dat` pra conferir o `descriptor`).
- **`Invalid axis range`**: alguma propriedade Veusz inválida no
  script (ex: `autoRange` com valor errado) — remove a linha, o eixo
  já auto-escala por padrão.
- **`$end` ou data group sumindo**: confirma que edições manuais
  (`cat >>`, `sed`) ficaram *antes* do `$end` no `control`
  (`grep -n '^\$end' control` vs. `wc -l control`).