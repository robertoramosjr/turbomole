# turbomole — scripts de extração e plotagem (pipeline DFT)

Scripts Python usados para processar as saídas de cálculos DFT do
Turbomole (Bader, TD-DFT, IR, DOS, RDF, ordem de ligação/SEN) e gerar
os datasets e figuras (Veusz) usados no trabalho de extração de dados
DFT. Companheiro do [MANUAL.md](MANUAL.md), que descreve o fluxo
completo do cálculo à figura final.

## Instalação

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

O pacote `veusz` traz o módulo de embedding (`veusz.embed`) usado
pelos scripts em `veusz/`. Em algumas distribuições Linux a instalação
via `pip` tem problemas com a dependência de Qt/PyQt — nesse caso,
prefira o pacote do sistema (ex.: `apt install veusz` /
`dnf install veusz`) e garanta que o `python3` usado enxerga os pacotes
do sistema (ou crie o venv com `--system-site-packages`).

Os scripts em `parse_*.py`, `compute_rdf.py`, `xyz_to_coord.py` etc.
(fora de `veusz/`) só precisam de `numpy` — funcionam mesmo sem Veusz
instalado, o que é útil para rodar a extração num nó do cluster que
não tem Qt.

## Estrutura

### Estágio 0 — preparação de geometria

| Script | Função |
|---|---|
| `xyz_to_coord.py` | Converte um `.xyz` padrão (Å) em `coord` no formato Turbomole (Bohr, símbolos em minúsculo, bloco `$coord`/`$end`). Primeiro passo antes do `define`. |

### Estágio 2 — extração local (parsers)

| Script | Função |
|---|---|
| `parse_bader.py` | Lê `ACF.dat` (Bader/Henkelman) + `coord`, calcula cargas líquidas de Bader (referência all-electron: `Z - população`). |
| `parse_excitations.py` | Lê `escf.out` (TD-DFT/RPA), gera espectro de sticks, envelope alargado ponderado por força do oscilador (JDOS óptico) e densidade de excitações de peso igual (painel de DOS). |
| `parse_ir.py` | Lê `vibspectrum` (saída do `aoforce`), descarta os 6 primeiros modos (translação/rotação), gera dataset de IR. |
| `parse_mulliken.py` | Lê a tabela de cargas de Mulliken em `pop.out` (`$pop mulliken`, via `-proper`). |
| `parse_bond_order.py` | Lê pares de "shared electron number" (SEN) do `pop.out` — análogo Mulliken ao COHP/ICOHP do VASP+LOBSTER. |
| `parse_dos.py` | Combina as rodadas de `$pop ... dos` (total + projetada por elemento, incluindo subcamadas s/p/d) num único dataset. |
| `compute_rdf.py` | RDF (função de distribuição radial) alargada por gaussiana, por grupo funcional, a partir de `coord` + `groups.json`. |

### Estágio 2 (opcional) — validação cruzada e agregação

| Script | Função |
|---|---|
| `compare_charges.py` | Compara cargas de Bader vs. Mulliken átomo a átomo (correlação, concordância de sinal, maiores discrepâncias). |
| `group_average_charges.py` | Média e faixa das cargas líquidas de Bader por grupo funcional (`groups.json`), para o texto do artigo. |

### Estágio 3 — figuras (Veusz)

| Script | Função |
|---|---|
| `veusz/plot_bader.py` | Gráfico de barras das cargas líquidas de Bader por átomo. |
| `veusz/plot_dos.py` | Grade 2×2: TDOS (com s/p/d), perfil tipo JDOS, pDOS por elemento. |
| `veusz/plot_ir.py` | Espectro de IR simulado: envelope lorentziano + sticks brutos sobrepostos. |
| `veusz/plot_bond_order.py` | Barras dos N pares átomo-átomo com maior SEN. |
| `veusz/plot_rdf.py` | Grade de painéis de RDF, um por grupo funcional passado via `--panel`. |

> **Faltando:** o manual (Estágio 3) também referencia um
> `veusz/plot_optical.py` para a figura de absorção óptica
> (`excitation_sticks.dat` + `jdos_dataset.dat`), que ainda não existe
> neste diretório. Precisa ser escrito antes de rodar essa etapa —
> ver aviso equivalente no [MANUAL.md](MANUAL.md).

### Validação e comparação adicionais (fora da sequência principal)

Não fazem parte do fluxo sequencial dos Estágios 1–3; são usados sob
demanda para seções específicas do artigo (robustez numérica,
comparação de funcionais, erro de IR vs. experimento):

| Script | Função |
|---|---|
| `compare_functionals.py` | Extrai a excitação TD-DFT dominante de vários `escf.out` (um por funcional) e monta tabela comparativa (ex.: HSE06 vs. híbrido global vs. híbrido de longo alcance). |
| `convergence_table.py` | Tabela de convergência numérica (energia SCF + excitação dominante) entre rodadas que variam grade de integração, threshold de SCF, base etc., com a mesma geometria. |
| `ir_error_stats.py` | MAE, RMSE e desvio absoluto máximo entre bandas de IR calculadas e experimentais, a partir de um CSV de pares já atribuídos manualmente. |

## Formato dos datasets

Os scripts de extração escrevem arquivos `.dat` no formato de
"descriptor" do Veusz (primeira linha `descriptor nome1,nome2,...`,
colunas separadas por espaço nas linhas seguintes) — importáveis
diretamente no Veusz ou pelos scripts `veusz/plot_*.py`, que os leem
via `veusz.embed`.

`groups.json` (usado por `compute_rdf.py` e `group_average_charges.py`)
mapeia nome do grupo funcional para os índices 1-based dos átomos no
`coord`, por exemplo:

```json
{
  "carboxylic_acid": [16, 21, 22, 2, 3, 31, 45, 46],
  "phenolic_ring": [4, 5, 6, 9, 10, 11, 1, 27, 28, 32]
}
```

## Uso

Cada script tem `--help` com a lista completa de flags. Para o fluxo
completo, do cálculo Turbomole à figura final (incluindo os comandos
de cada etapa e os checks de sanidade entre elas), veja o
[MANUAL.md](MANUAL.md).

## Avisos conhecidos

- `veusz/plot_optical.py` (figura de absorção óptica) está referenciado
  no manual mas não existe neste repositório ainda.
- `veusz/plot_bond_order.py` recebe `--data` mas não usa o arquivo
  diretamente (mantido só por consistência de interface com os outros
  `plot_*.py`) — os dados vêm de `--labels`.
- A tabela `ATOMIC_NUMBER` em `parse_bader.py` cobre só H–Ca; moléculas
  com elementos fora dessa faixa precisam da tabela estendida antes de
  confiar nas cargas líquidas (o script avisa e grava `nan` nesse
  caso).
