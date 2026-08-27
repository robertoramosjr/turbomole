# turbomole — scripts de extração e plotagem (pipeline DFT)

Scripts Python usados para processar as saídas de cálculos DFT do
Turbomole (Bader, TD-DFT, IR, DOS/PDOS — Mulliken, Becke ou SCPA —,
correção GW-BSE, exciton binding, RDF, ordem de ligação/SEN) e gerar
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

Duas dependências **não** cobertas por este `requirements.txt`, só
necessárias para partes específicas do pipeline (ver
[MANUAL.md § Pré-requisitos](MANUAL.md#12-pré-requisitos)):

- **Multiwfn + Molden2AIM** (binários compilados à parte, em
  `~/local_bin/`) — só para as rotas de PDOS via Becke/SCPA
  (`parse_becke_pdos.py`, `parse_dos_scpa.py`,
  `parse_scpa_pdos_angular.py`).
- **Ambiente conda `vasp_env`** (`rdkit`, `ase`, `povray`) — só para as
  figuras de estrutura molecular específicas de Artepillin C
  (`plot_bond_order_structure_artepillinC.py`,
  `plot_molecule_vesta_style.py`).

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
| `parse_dos.py` | Combina as rodadas de `$pop ... dos` (total + projetada por elemento, incluindo subcamadas s/p/d) num único dataset. PDOS por elemento aqui é Mulliken — pode ter canais negativos (ver PDOS via Multiwfn abaixo). A flag chama-se `--oxygen` mesmo para o heteroátomo de moléculas sem oxigênio (ex. azobenzeno/N) — funciona numericamente, só o rótulo da coluna de saída fica `pDOS_O`. |
| `compute_rdf.py` | RDF (função de distribuição radial) alargada por gaussiana, por grupo funcional, a partir de `coord` + `groups.json`. |

### Estágio 2 — PDOS alternativas via Multiwfn (opcional, sem o artefato de sinal do Mulliken)

Entrada: `DOS_curve.txt` (e `orginfo.txt`) exportados do Multiwfn a
partir do Molden corrigido pelo Molden2AIM — ver
[MANUAL.md, Estágios 1b-alt2/1b-alt3](MANUAL.md#6-estágio-1b-alt2--pdos-por-elemento-via-becke-multiwfn--molden2aim-opcional).

| Script | Função |
|---|---|
| `parse_becke_pdos.py` | PDOS por elemento (O/C/H) via partição Becke (espaço real) — garante contribuições não-negativas; não resolve momento angular. |
| `parse_scpa_pdos_angular.py` | PDOS por elemento × momento angular (s/p) via partição SCPA (função de base) — 6 fragmentos fixos na ordem C-s, C-p, O-s, O-p, H-s, H-p. |
| `parse_dos_scpa.py` | Generalização de `parse_scpa_pdos_angular.py`: mesmo método SCPA, qualquer conjunto de fragmentos via `--labels`, na ordem definida no Multiwfn. |

### Estágio 2 — DOS via correção quasipartícula G0W0 (opcional)

| Script | Função |
|---|---|
| `parse_qpenergies_dos.py` | Lê `qpenergies.dat` (correções G0W0 do Turbomole) e gera duas densidades orbitais alargadas por gaussiana (não-projetadas, peso igual por orbital), uma no nível Kohn-Sham e outra G0W0, na mesma grade de energia. |

### Estágio 2 (opcional) — validação cruzada e agregação

| Script | Função |
|---|---|
| `compare_charges.py` | Compara cargas de Bader vs. Mulliken átomo a átomo (correlação, concordância de sinal, maiores discrepâncias). |
| `group_average_charges.py` | Média e faixa das cargas líquidas de Bader por grupo funcional (`groups.json`), para o texto do artigo. |

### Estágio 3 — figuras (Veusz)

| Script | Função |
|---|---|
| `veusz/plot_bader.py` | Gráfico de barras das cargas líquidas de Bader por átomo. |
| `veusz/plot_dos.py` | Grade 2×2: TDOS (com s/p/d), perfil tipo JDOS, pDOS por elemento (`--element4-label` define o rótulo do 4º elemento; `--homo-lumo` sobrepõe HOMO/LUMO + gap, opcional). |
| `veusz/plot_optical.py` | Espectro de absorção óptica TD-DFT: envelope gaussiano alargado sobreposto ao espectro de sticks bruto. |
| `veusz/plot_ir.py` | Espectro de IR simulado: envelope lorentziano + sticks brutos sobrepostos. |
| `veusz/plot_bond_order.py` | Barras dos N pares átomo-átomo com maior SEN. |
| `veusz/plot_rdf.py` | Grade de painéis de RDF, um por grupo funcional passado via `--panel` (`--columns` controla o número de colunas). |
| `veusz/plot_becke_pdos.py` | Plota `dos_becke_dataset.dat` (TDOS + pDOS-O/C/H sobrepostos), com HOMO/LUMO marcados via `--orginfo` (deriva de `orginfo.txt` do Multiwfn) ou `--homo`/`--lumo` diretos. |
| `veusz/plot_scpa_pdos_angular.py` | Plota a PDOS SCPA por elemento × momento angular (cor por elemento, estilo de linha sólido/tracejado por s/p); `--sigma` aplica suavização gaussiana extra — manter bem abaixo da largura do gap HOMO-LUMO. |
| `veusz/plot_qpenergies_dos.py` | Empilha as densidades orbitais Kohn-Sham e G0W0 de `parse_qpenergies_dos.py` em dois painéis com o mesmo eixo de energia, para comparar a abertura do gap sob a correção GW. |

### Estágio 3 — variantes específicas de Artepillin C

Hardcoded para a atomagem de Artepillin C — não reaproveitar para
outra molécula sem revisar `bond_order_dataset_labels.csv` primeiro:

| Script | Ambiente | Função |
|---|---|---|
| `veusz/plot_bond_order_artepillinC.py` | base (`.venv`) | Como `plot_bond_order.py`, mas colorido por categoria química (anel aromático / carbonila / C=C / outros), com rótulos nas ligações citadas no artigo. |
| `plot_bond_order_structure_artepillinC.py` | conda `vasp_env` (rdkit) | Figura combinada: barras de SEN + estrutura 2D via RDKit com as mesmas ligações coloridas direto na geometria (projeção PCA da geometria DFT relaxada). |
| `plot_molecule_vesta_style.py` | conda `vasp_env` (ase + povray) | Renderização 3D estilo VESTA (ball-and-stick) da geometria relaxada via POV-Ray. |

> `veusz/plot_dos_molecules.py` é uma cópia mais antiga de
> `plot_dos.py` (sem `--element4-label` nem `--homo-lumo`) — mantida
> por enquanto, mas prefira `plot_dos.py` para figuras novas.

### Estágio 4 — comparações entre rodadas/níveis de teoria (sob demanda)

Não fazem parte do fluxo sequencial dos Estágios 1–3; são usados sob
demanda para seções específicas do artigo (robustez numérica,
comparação de funcionais/níveis de teoria, erro de IR vs.
experimento):

| Script | Função |
|---|---|
| `compare_functionals.py` | Extrai a excitação TD-DFT dominante de vários `escf.out` (um por funcional) e monta tabela comparativa (ex.: HSE06 vs. híbrido global vs. híbrido de longo alcance). |
| `compare_low_states.py` | Extrai os N estados excitados mais baixos de várias rodadas `escf.out`/`escf_bse.out` (TD-DFT com funcionais diferentes, ou GW-BSE/TDA) numa única tabela longa, com energia real, gap orbital do par dominante e o binding por estado entre os dois. |
| `compute_exciton_binding.py` | Energia de ligação do exciton, `E_gap(HOMO-LUMO) - E_exciton(S1)`, por nível de teoria — gap Kohn-Sham para TD-DFT, gap G0W0 (de `qpenergies.dat`) para GW-BSE/TDA. |
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

- `veusz/plot_bond_order.py` recebe `--data` mas não usa o arquivo
  diretamente (mantido só por consistência de interface com os outros
  `plot_*.py`) — os dados vêm de `--labels`.
- A tabela `ATOMIC_NUMBER` em `parse_bader.py` cobre só H–Ca; moléculas
  com elementos fora dessa faixa precisam da tabela estendida antes de
  confiar nas cargas líquidas (o script avisa e grava `nan` nesse
  caso).
- `parse_dos.py` mantém a flag `--oxygen` mesmo em moléculas sem
  oxigênio — ver nota na tabela do Estágio 2 acima.
- `veusz/plot_dos_molecules.py` é uma cópia mais antiga de
  `veusz/plot_dos.py`, mantida no repositório mas não recomendada para
  figuras novas (ver Estágio 3).
- O caminho exato do menu interativo do Multiwfn para a rota SCPA
  (`parse_scpa_pdos_angular.py`/`parse_dos_scpa.py`) ainda não está
  documentado passo a passo — só o formato de entrada/saída dos
  scripts. Ver [MANUAL.md, Estágio 1b-alt3](MANUAL.md#7-estágio-1b-alt3--pdos-por-elemento--momento-angular-via-scpa-multiwfn-opcional).
