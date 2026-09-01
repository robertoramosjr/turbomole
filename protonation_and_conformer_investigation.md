# Investigação de Estados de Protonação e População Conformacional (Paper 2)

## 1. Visão geral — o que é isto e por que existe

O Paper 1 caracterizou a artepilina C usando **um único confôrmero de
referência** (HSE06-D3(BJ)/def2-TZVP, via Turbomole) — a molécula
neutra, numa única geometria. Isso é uma simplificação: em solução, a
qualquer temperatura acima do zero absoluto, uma molécula flexível não
existe como uma única geometria fixa, e sim como uma **população de
confôrmeros** em equilíbrio térmico, cada um contribuindo para as
propriedades observadas (espectro, energia, etc.) proporcionalmente ao
seu peso de Boltzmann. Além disso, a artepilina C tem dois grupos
ácidos (carboxila e fenol) cujo estado de protonação muda com o pH —
então "a molécula" na verdade significa "uma de quatro espécies
químicas possíveis", dependendo do pH do meio.

O Paper 2 substitui essa simplificação por algo mais realista: para
cada estado de protonação relevante, gerar a **população real de
confôrmeros** (não só o mínimo global) e usá-la para calcular
propriedades ponderadas por população — em particular, espectros de
IV (infravermelho) e UV-Vis diretamente comparáveis aos espectros
experimentais já citados no Paper 1.

**O obstáculo é puramente computacional.** Uma busca conformacional
tipicamente gera 20–300+ estruturas por estado de protonação. Refazer
o cálculo completo de produção do Paper 1 (HSE06/def2-TZVP + Hessiana)
em cada uma delas é inviável — só a Hessiana de produção de **um**
confôrmero já levou ~1h no Paper 1. Este documento descreve o **funil**
construído para resolver isso: uma sequência de etapas que vai
ficando mais cara e mais precisa, descartando estruturas
irrelevantes a cada passo, até sobrar só um punhado de confôrmeros que
merecem o tratamento completo de produção.

### 1.1. Como ler este documento

Cada seção de estágio segue o mesmo formato: **o que faz**, **por que
esta escolha e não outra** (a parte que interessa mesmo para quem não
vai rodar o comando, só entender a decisão), **qual script usar** e
**quando usar**. Trechos marcados como "⚠️ Pegadinha" documentam um erro
real que já aconteceu nesta investigação e como foi corrigido — vale
ler antes de repetir o mesmo passo, porque alguns desses erros
produziram resultados numericamente plausíveis mas errados (não
travaram com uma mensagem de erro óbvia).

### 1.2. Convenção de generalidade dos scripts

Todo script novo escrito para esta investigação segue a mesma
convenção do restante do repositório (ver [MANUAL.md](MANUAL.md) /
[README.md](README.md)): **nenhum é hardcoded para a artepilina C**.
Carga, número de átomos, geometria de entrada, funcional, base, nível
de energia — tudo é `argparse`, nada fixo no código. Isso significa
que **o funil inteiro pode ser reaplicado a qualquer outra molécula ou
conjunto de moléculas** sem editar uma linha de código-fonte, só
trocando os argumentos de linha de comando (geometria de entrada,
cargas dos estados de protonação, etc.). A única parte específica da
artepilina C é a *química* das decisões (quais grupos são
protonáveis, quantos estados existem) — o *código* é genérico.

Exceção conhecida (herdada do resto do repositório, não desta
investigação): os scripts `plot_bond_order_structure_artepillinC.py` e
similares em `veusz/`, que têm a numeração de átomos da artepilina C
hardcoded — não fazem parte deste funil.

---

## 2. Os 4 estados de protonação

A artepilina C é um ácido diprótico (carboxila + fenol), com 4
microestados possíveis:

| Microestado | Carboxila | Fenol | Carga total | Prioridade |
|---|---|---|---|---|
| Neutro | -COOH | -OH | 0 | Alta — já existia como confôrmero único no Paper 1, agora ganha ensemble completo |
| Monoânion, tautômero A (carboxilato) | -COO⁻ | -OH | −1 | Alta — dominante em pH fisiológico |
| Diânion | -COO⁻ | -O⁻ | −2 | Média |
| Monoânion, tautômero B (fenolato) | -COOH | -O⁻ | −1 | Baixa — população desprezível esperada em qualquer pH relevante (diferença de várias unidades entre os dois pKa), mas incluído como checagem de consistência |

Nenhum grupo da molécula é básico o suficiente para ser protonado na
janela de pH aquoso (0–14) — **não existe estado catiônico** nesta
investigação.

**Regra que percorre todo o funil**: nunca comparar energia absoluta
entre microestados de carga diferente sem uma correção de referência
apropriada (ex. referência a H⁺/H₃O⁺ em solução). O ranking/filtro de
cada etapa roda **dentro** de cada estado de protonação, nunca entre
eles.

Decisão tomada em 2026-08-28: **todos os 4 estados rodam o funil
completo**, sem exceção — inclusive o neutro, que precisa de uma busca
CREST nova (o confôrmero único do Paper 1 não é reaproveitado como
"ensemble", já que o objetivo aqui é ter uma população real e
comparável aos outros 3 estados).

---

## 3. Pré-requisitos / ambiente

- Ambiente conda dedicado `paper2_funnel` (criado 2026-08-27, via
  `conda create -n paper2_funnel -c conda-forge crest xtb scikit-learn
  numpy scipy pandas matplotlib`) — crest 3.0.2, xtb 6.7.1,
  scikit-learn 1.9.0. **Gotcha de PATH**: `~/software/TURBOMOLE/bin/...`
  vem antes do `bin/` do ambiente conda no `PATH` mesmo depois de
  `conda activate` — uma chamada nua a `xtb` resolve para a versão
  mais antiga (6.6.1) empacotada junto do Turbomole, não a do ambiente
  (6.7.1). Não afeta o CREST em si (ele usa tblite internamente, não
  chama o binário `xtb`), mas use o caminho absoluto
  (`~/.conda/envs/paper2_funnel/bin/xtb`) se algum script novo chamar
  `xtb` diretamente.
- `pexpect` (`pip install pexpect` dentro do ambiente) — necessário
  para os scripts que automatizam o `define` do Turbomole via um
  pseudo-terminal real (ver Estágios 3 e 4; `stdin` via pipe simples
  se comporta de forma diferente/enganosa com o `define`, não use
  heredoc para isso).
- TURBOMOLE já instalado em `~/software/TURBOMOLE`, no `PATH`.
- Cluster GridUNESP via **SLURM** (`sbatch`/`squeue`/`sinfo`) — **não**
  PBS, apesar de arquivos antigos (`job.pbs`) no repositório sugerirem
  o contrário (`qsub`/`qstat` nem estão instalados nesta máquina).
  Partições: `short` (1 dia), `medium` (7 dias), `long` (30 dias),
  `gpu`. Conta `laser` (única, não precisa passar `--account`
  explicitamente).

---

## 4. Estágio 0 — Geometria de partida por estado de protonação

**Objetivo:** a partir da geometria neutra de referência (Paper 1),
gerar a geometria de partida (carga + átomos removidos) de cada
estado de protonação.

**Como funciona:** em vez de hardcodar índices de átomo, o script
detecta os grupos hidroxila automaticamente por geometria (um O ligado
a exatamente um C e um H; se esse C também tiver um segundo O a
distância de carbonila, classifica como "carboxílico", senão "outro" —
ex. fenólico). Isso foi validado contra uma checagem manual de
distância de ligação antes de confiar no resultado.

**Script:** `build_protonation_state.py`

```bash
# inspeciona os grupos hidroxila detectados, não escreve nada
python ~/work_turbomole/scripts/build_protonation_state.py \
    --xyz structure.xyz --list-oh

# monoânion tautômero A (carboxilato)
python ~/work_turbomole/scripts/build_protonation_state.py \
    --xyz structure.xyz --remove carboxylic \
    --output-dir artepillin_c_monoanion_a_d3_hse06

# monoânion tautômero B (fenolato)
python ~/work_turbomole/scripts/build_protonation_state.py \
    --xyz structure.xyz --remove other \
    --output-dir artepillin_c_monoanion_b_d3_hse06

# diânion
python ~/work_turbomole/scripts/build_protonation_state.py \
    --xyz structure.xyz --remove carboxylic other \
    --output-dir artepillin_c_dianion_d3_hse06
```

Cada pasta de saída recebe `structure.xyz`, `coord` (Turbomole, via
`xyz_to_coord.py`) e `CHARGE.txt` com a carga total resultante.

**Convenção de nome de pasta**: `artepillin_c_<estado>_d3_<funcional>`
— o funcional faz parte do nome porque múltiplos funcionais (ex. HSE06
e outro em consideração) podem eventualmente reotimizar a mesma
molécula; a geometria de partida em si não depende do funcional, só o
nome da pasta em que ela é colocada.

**Validação feita:** GFN2-xTB single-point em cada estado, conferindo
carga total e convergência de SCF antes de prosseguir.

---

## 5. Estágio 1 — Busca conformacional (CREST/GFN2-xTB)

**Objetivo:** gerar o ensemble bruto de confôrmeros de cada estado de
protonação.

**Por que CREST/GFN2-xTB:** é o padrão de fato para busca
conformacional em química computacional — o algoritmo iMTD-GC
(metadinâmica com viés iterativo + cruzamento genético de estruturas)
explora a superfície de energia potencial de forma muito mais
eficiente que uma busca exaustiva de ângulos diedros, e o nível
semi-empírico GFN2-xTB é rápido o bastante para rodar centenas de
otimizações sem ser proibitivo.

**Script:** `job_crest_search.sh` (SLURM)

```bash
sbatch ~/work_turbomole/scripts/job_crest_search.sh \
    /caminho/para/artepillin_c_<estado>_d3_hse06 [janela_energia_kcal]
```

Lê a carga de `CHARGE.txt` na pasta alvo, roda
`crest coord --gfn2 --chrg <carga> -T <threads> -ewin <janela>`
(padrão iMTD-GC, janela padrão 6 kcal/mol).

**⚠️ Pegadinha corrigida — oversubscription de threads:** a primeira
versão do script setava `OMP_NUM_THREADS=$SLURM_NTASKS` **junto** com
o `-T $SLURM_NTASKS` do próprio CREST. O CREST já paraleliza
internamente entre as ~14 metadinâmicas simultâneas de cada iteração;
com `OMP_NUM_THREADS` também setado, cada uma dessas metadinâmicas
tentava *ela mesma* usar todos os threads do OpenBLAS, gerando uma
explosão de threads competindo pelos mesmos núcleos. Sintoma: logs de
500 MB–1 GB só de avisos repetidos do OpenBLAS ("Detect OpenMP Loop"),
e um desempenho muito mais lento que um benchmark de referência
single-thread. Correção: `OMP_NUM_THREADS=1` +
`OPENBLAS_NUM_THREADS=1`, deixando só o `-T` do CREST controlar a
concorrência entre metadinâmicas.

**Resultado obtido (todos os 4 estados, 2026-08-28):**

| Estado | Confôrmeros únicos (janela 6 kcal/mol) | População do mais estável | Wall-time (16 threads) |
|---|---|---|---|
| Neutro | 337 | 8.5% | 56 min |
| Monoânion A (carboxilato) | 322 | 4.8% | 1h05 |
| Monoânion B (fenolato) | 37 | 36.9% | 1h17 |
| Diânion | 21 | 46.6% | 54 min |

Note que fenolato e diânion saíram com ensembles bem menores e mais
concentrados num único confôrmero — plausível quimicamente (sítios
aniônicos tendem a travar a geometria via eletrostática/ligação de
H intramolecular), mas vale conferir visualmente antes de dar como
certo, não é intuitivo à primeira vista.

**Saída usada nas próximas etapas:** `crest_conformers.xyz` (geometrias
do ensemble), `crest.energies` (energia relativa GFN2-xTB, kcal/mol) e
`cre_members` (degenerescência de cada confôrmero — ver Estágio 2).

---

## 6. Estágio 2 — Camada A, parte 1: filtro de população de Boltzmann

**Objetivo:** descartar, dentro de cada estado de protonação, os
confôrmeros cuja população de Boltzmann é irrelevante à temperatura
ambiente — não vale gastar DFT caro num confôrmero que praticamente
não existe em equilíbrio térmico.

**Script:** `filter_boltzmann_population.py`

```bash
python ~/work_turbomole/scripts/filter_boltzmann_population.py \
    --energies crest.energies --conformers crest_conformers.xyz \
    --members cre_members \
    --output-dir layer_a [--cutoff 0.95] [--temp 298.15]
```

Mantém, em ordem de energia, o menor conjunto de confôrmeros cuja
população cumulativa atinge o corte (padrão 95% a 298.15 K), e escreve
`population_filtered.xyz` (sobreviventes) + `population_table.csv`
(tabela completa, incluindo os descartados, para transparência).

**⚠️ Pegadinha corrigida — degenerescência de rotâmeros:** o CREST já
deduplica automaticamente rotâmeros quase-idênticos (ex. uma rotação
de metila) num único confôrmero representante em
`crest_conformers.xyz`/`crest.energies` — mas isso apaga a informação
de *quantos* rotâmeros aquele confôrmero representa. Tratar cada linha
de `crest.energies` como um microestado de peso igual **subestima**
sistematicamente a população dos confôrmeros que tinham mais rotâmeros
equivalentes: no estado neutro, o mínimo global representa 58
rotâmeros quase-idênticos — sem pesar por isso, sua população calculada
saía 1.2% em vez dos 8.5% corretos (valor que o próprio CREST já
reporta internamente). Correção: peso de Boltzmann de cada confôrmero
multiplicado pela sua degenerescência, lida da coluna 1 do arquivo
`cre_members` do CREST. Verificado batendo com o valor interno do
CREST a menos de 0.02 pontos percentuais em todos os 4 estados.

**Resultado obtido (corte de 95% cumulativo, 298.15 K):**

| Estado | Sobreviventes / total | População do dominante |
|---|---|---|
| Neutro | 170 / 337 | 8.5% |
| Monoânion A | 143 / 322 | 4.8% |
| Monoânion B | 11 / 37 | 36.9% |
| Diânion | 3 / 21 | 46.6% |

---

## 7. Estágio 3 — Camada A, parte 2: refino em DFT barato (r2SCAN-3c)

**Objetivo:** reotimizar geometria + energia dos sobreviventes do
Estágio 2 num nível de DFT muito mais confiável que GFN2-xTB, mas
ainda tratável para centenas de estruturas — e usar essa energia para
reranquear a população dentro de cada estado.

**Por que r2SCAN-3c:** é um "método composto" (r2SCAN meta-GGA +
base def2-mTZVP + dispersão D4 + correção de superposição de base gCP,
tudo com parâmetros pré-otimizados) — dá geometrias e energias
próximas de DFT híbrido de qualidade a uma fração do custo, sem
precisar montar manualmente dispersão/correção de base como seria
necessário com uma GGA simples tipo PBE-D3/def2-SVP.

**Scripts (nesta ordem):**

```bash
# 1. Um template por estado (define roda uma vez por carga, não por confôrmero —
#    a atribuição de base não depende da geometria)
python ~/work_turbomole/scripts/setup_r2scan3c_template.py \
    --xyz structure.xyz --charge <carga> --output-dir layer_a/template

# 2. Clona o template pra cada confôrmero sobrevivente, trocando só o coord
python ~/work_turbomole/scripts/fanout_conformer_jobs.py \
    --template layer_a/template --ensemble layer_a/population_filtered.xyz \
    --output-dir layer_a/conformer_jobs

# 3. Lista combinada de todos os estados + submissão em array SLURM
#    (ver Estágio 5 do repositório-mãe para o padrão de array; aqui,
#    um exemplo mínimo)
sbatch --array=1-N%30 ~/work_turbomole/scripts/job_r2scan3c_optimize.sh \
    <arquivo_com_uma_pasta_por_linha>

# 4. Depois que o array termina: extrai a energia final de cada pasta
#    e reranqueia a população da Camada A usando DFT em vez de GFN2-xTB
python ~/work_turbomole/scripts/rerank_layer_a_dft.py \
    --population-table layer_a/population_table.csv \
    --conformer-jobs layer_a/conformer_jobs \
    --output layer_a/population_table_dft.csv
```

**⚠️ Pegadinhas corrigidas ao automatizar o `define`:**

1. **Dois prompts em branco separados no início** do `define` (um para
   "ler defaults de outro control?", outro para o título) — tratar
   como um só desalinha todo o resto do roteiro (sintoma visto:
   "NO ATOMS, NO MOLECULE, NOTHING!" mesmo com um `coord` válido).
2. **Métodos "-3c" não têm a base atribuída automaticamente** ao
   escolher o funcional no menu `dft` — fica com a base genérica
   padrão (`def-SV(P)`, pequena demais) até você mesmo atribuir
   `def2-mTZVP` explicitamente com `b all def2-mTZVP` **antes** de
   entrar no menu `dft`. Em compensação, a dispersão D4 e a correção
   gCP **são** aplicadas automaticamente pelo `ridft`/`jobex` ao
   detectar o funcional `r2scan-3c` — só a base precisa de ajuda.
3. **`$rij` precisa ser adicionado manualmente** ao `control` (o
   `define` não liga sozinho nesta sequência) — e, se for inserido à
   mão num `control` já existente, tem que ficar **antes** do `$end`,
   nunca depois (o Turbomole ignora tudo depois do `$end`).
4. **Critério de convergência geométrica**: os grupos prenila da
   artepilina C são flexíveis o bastante para que o `jobex` não
   declare convergência formal mesmo depois de 35+ ciclos — mas a
   energia total já está estável na 6ª casa decimal (Hartree) desde
   por volta do ciclo 15–17, cinco ordens de grandeza abaixo da escala
   que importa para peso de Boltzmann (~1 kcal/mol). Decisão: usar um
   **teto fixo de ciclos** (`jobex -c 10`) como critério de parada
   nesta etapa de triagem, em vez de esperar convergência formal.
   **Isso não se aplica** à reotimização final de produção (Estágio 8),
   que mantém os critérios rígidos originais do Paper 1.
5. **Reranking por DFT precisa da mesma correção de degenerescência**
   do Estágio 2 — `rerank_layer_a_dft.py` usa a coluna `degeneracy` que
   `filter_boltzmann_population.py` já carrega de `cre_members`. Foi
   um quase-erro real: a primeira versão deste script esqueceu a
   ponderação e teria reintroduzido o mesmo viés corrigido no Estágio
   2, só que na reclassificação por DFT.

**Resultado obtido:** dos 327 confôrmeros no total, 326 convergiram
sem intervenção; 1 (monoânion A, confôrmero #125) precisou de
`$scfiterlimit` maior e amortecimento inicial (`$scfdamp start`) mais
forte para o SCF convergir — corrigido manualmente, sem impacto no
restante. **O confôrmero de maior população mudou** entre GFN2-xTB e
r2SCAN-3c em 3 dos 4 estados (neutro, monoânion A, monoânion B) — só o
diânion manteve o mesmo confôrmero dominante. Isso é o resultado
esperado (é exatamente por isso que existe esta etapa de refino em
DFT), não um sinal de erro.

---

## 8. Estágio 4 — Camada A, parte 3: espectros de IV e UV-Vis

**Objetivo:** gerar o espectro de infravermelho (via Hessiana,
`aoforce`) e o espectro de UV-Vis (via TD-DFT, `escf`) de cada
sobrevivente, para depois compor o **espectro ponderado por
população** de cada estado de protonação — o entregável principal
desta camada, comparável direto aos espectros experimentais já citados
no Paper 1.

### 8.1. Por que não r2SCAN-3c também aqui

Seria natural tentar reaproveitar o mesmo r2SCAN-3c do Estágio 3. Não
dá: confirmado rodando de verdade (não só lendo a documentação) que
tanto `aoforce` quanto `escf` **abortam** neste Turbomole com o erro
`Invalid value of nfun in <mgga_r2>!` — uma limitação real de
implementação para derivadas analíticas (Hessiana e resposta de
TD-DFT) com funcionais meta-GGA como o r2SCAN, não uma questão de
parâmetro. (Um teste inicial mais curto, que parou antes do `aoforce`
chegar na parte que falha, sugeriu erroneamente que "funcionava" — a
lição prática é não declarar uma capacidade confirmada a partir de um
teste truncado.)

**Solução adotada**: **PBE0/def2-SVP** (híbrido global + base padrão,
com RI-K, já que RI-J puro não converge em tempo viável para um
híbrido) rodado como single-point sobre a geometria já otimizada em
r2SCAN-3c — só a etapa de propriedade espectroscópica muda de
funcional/base, a geometria e o ranking de energia continuam em
r2SCAN-3c.

**Scripts:**

```bash
# 1. Template PBE0/def2-SVP por estado (RI-J + RI-K + bloco TD-DFT)
python ~/work_turbomole/scripts/setup_pbe0_svp_template.py \
    --xyz structure.xyz --charge <carga> --nstates 10 \
    --output-dir layer_a/pbe0_template

# 2. Clona o template, usando a geometria JÁ OTIMIZADA em r2SCAN-3c
#    (o Turbomole atualiza 'coord' no lugar durante o jobex, então não
#    precisa reconverter nada)
python ~/work_turbomole/scripts/fanout_pbe0_jobs.py \
    --template layer_a/pbe0_template \
    --r2scan3c-jobs layer_a/conformer_jobs \
    --output-dir layer_a/pbe0_jobs

# 3. Array SLURM: ridft (reconverge SCF) -> aoforce (Hessiana/IV) -> escf (TD-DFT/UV-Vis)
sbatch --array=1-N%25 ~/work_turbomole/scripts/job_pbe0_ir_uvvis.sh \
    <arquivo_com_uma_pasta_por_linha>
```

**⚠️ Pegadinhas corrigidas:**

- **Não reaproveite a base `def2-mTZVP` do r2SCAN-3c para o PBE0/RI-K**
  — a atribuição automática de base auxiliar RI-JK falha silenciosamente
  para essa base (o `ridft` recusa com "Problem reading basis set(s)"):
  `def2-mTZVP` é uma base sob medida para a receita "-3c", sem
  cobertura completa na biblioteca RI-JK do Turbomole. `def2-SVP` é
  uma base padrão, com cobertura completa, e funciona sem drama.
- **PBE0 precisa de RI-K** para ser tratável — confirmado que RI-J puro
  não passou da 2ª iteração de SCF em 4 minutos.
- O `enable_rik` já existente em `prepare_functional_benchmark.py` usa
  uma sequência de `stdin` de tamanho fixo para o `define`, que
  pressupõe um número específico de confirmações de "apagar data group
  órfão?" — não generaliza (falhou aqui, porque nosso `control` tinha
  um número diferente de grupos órfãos). `setup_pbe0_svp_template.py`
  usa uma abordagem mais robusta: keep enviando linha em branco até
  reconhecer um marcador de menu esperado (ex. "GENERAL MENU"), em vez
  de contar quantos "enters" mandar às cegas.

### 8.2. Custo real medido — e por que isso importa

Medido diretamente num confôrmero (8 núcleos, mesma máquina/partição
usada no restante do funil):

| Etapa | Tempo (wall-clock) |
|---|---|
| `ridft` (reconverge SCF) | 9min31s |
| `aoforce` (Hessiana/IV) | **3h00min42s** |
| `escf` (TD-DFT/UV-Vis, 10 estados) | 16min58s |
| **Total por confôrmero** | **~3h27min** |

Isso é **~3x mais caro** que a própria Hessiana de produção do Paper 1
em HSE06/def2-TZVP (~1h), apesar do PBE0/def2-SVP ser nominalmente o
nível "barato" desta etapa — o custo vem da troca exata (RI-K) do
híbrido PBE0, não do tamanho da base. **Não assuma que um nível
"mais barato" no nome é automaticamente mais rápido na prática** —
meça antes de escalar.

Para os 327 confôrmeros dos 4 estados juntos, isso soma **~9.000
núcleo-horas**, e mesmo com 25 tarefas simultâneas no array SLURM,
**~37–40 horas de wall-time** (quase 2 dias corridos).

### 8.3. A alternativa considerada e descartada — e quando ela faz sentido

Antes de aceitar esse custo, foi avaliada a alternativa de **não**
rodar IR/UV-Vis nos 327 confôrmeros, e sim só num subconjunto pequeno
(os mesmos ~15–30 por estado que a Camada B, descrita no Estágio 6,
já ia rotular para o aprendizado ativo), usando um modelo de ML para
**prever** a propriedade espectroscópica do resto do ensemble e montar
o espectro ponderado a partir daí.

**Essa alternativa foi rejeitada para esta rodada** — decisão explícita
de priorizar acurácia total sobre economia de tempo de cluster,
aceitando as ~2 dias de wall-time do cálculo completo em DFT. Mas ela
continua sendo a estratégia certa para casos onde o ensemble é grande
demais para o orçamento disponível (por exemplo, uma molécula muito
mais flexível que a artepilina C, ou um estudo com muitas moléculas em
paralelo) — daí valer documentar como fazer direito, para quando for
necessário.

**Importante: o modelo de ML certo aqui não é o mesmo da Camada B**
(Estágio 6). O GPR com kernel Matérn 5/2 usado lá foi desenhado e
validado para **ranking** (Spearman ρ, Recall@10) — acertar *qual*
confôrmero tem o maior desvio de propriedade, não o **valor absoluto**
da propriedade em si. Usar esse mesmo modelo para gerar os pontos
reais de um espectro simulado (posição da banda, intensidade) seria
usá-lo fora do que foi validado a fazer: um modelo pode ranquear
perfeitamente e ainda errar o valor absoluto o bastante para distorcer
a forma da banda.

Para essa finalidade — prever o **valor** de uma propriedade
espectroscópica (frequência de um modo vibracional, energia de
excitação, força de oscilador) com acurácia suficiente para reconstruir
um espectro —, a literatura de ML para espectroscopia (delta-learning
sobre observáveis espectroscópicos, ex. Ramakrishnan *et al.* 2015 para
a ideia geral de Δ-ML; trabalhos dos grupos Dral e Barbatti
especificamente sobre espectros IV/UV-Vis) recomenda um desenho
diferente:

- **Modelo**: Kernel Ridge Regression (KRR) ou GPR com kernel RBF/Matérn
  com ARD (*Automatic Relevance Determination* — um comprimento de
  escala por *feature*, não um só) em vez do Matérn 5/2 isotrópico da
  Camada B. Alternativa robusta se o conjunto rotulado crescer além de
  ~50–100 pontos: Random Forest ou Gradient Boosting Regressor sobre
  descritores estruturais (distâncias/ângulos relevantes, cargas) —
  mais tolerante a uma superfície de resposta pouco suave que GPR/KRR.
- **Métrica de validação**: MAE e RMSE por validação cruzada
  (leave-one-out ou k-fold) contra o valor de DFT — **não**
  Spearman/Recall@10. O objetivo aqui é acurácia absoluta, não
  ordenação relativa.
- **Delta-learning continua fazendo sentido** (prever a diferença entre
  o valor caro e um proxy barato, não o valor absoluto do zero) — só a
  métrica de sucesso e a escolha de modelo/kernel mudam em relação à
  Camada B.
- **Reporte a incerteza**: se o modelo tiver variância preditiva
  calibrada (GPR/KRR com kernel bayesiano), propague essa incerteza
  para a média ponderada por população do espectro final (ex. banda de
  incerteza em torno da curva), em vez de tratar a previsão como exata.

Se esta rota for adotada no futuro, o script de treino/validação **não
deve reaproveitar** o mesmo código de ranking da Camada B — precisa de
um script próprio, validado por MAE/RMSE, específico para esse fim.
Não implementado nesta rodada (decisão do usuário, 2026-08-29:
priorizar acurácia total).

### 8.4. Montagem do espectro ponderado por população

**Objetivo:** combinar os espectros individuais dos sobreviventes de
cada estado (Estágio 4.3) num único envelope contínuo, ponderado pela
população de cada confôrmero (a mesma população DFT do Estágio 3.4) —
este sim é o entregável final da Camada A.

**Script:** `build_population_weighted_spectrum.py`

```bash
python ~/work_turbomole/scripts/build_population_weighted_spectrum.py \
    --kind ir \
    --population-table layer_a/population_table_dft.csv \
    --jobs-dir layer_a/pbe0_jobs \
    --output layer_a/ir_population_weighted.dat

python ~/work_turbomole/scripts/build_population_weighted_spectrum.py \
    --kind uvvis \
    --population-table layer_a/population_table_dft.csv \
    --jobs-dir layer_a/pbe0_jobs \
    --output layer_a/uvvis_population_weighted.dat
```

Cada confôrmero contribui com seu próprio espectro em sticks (IV:
número de onda + intensidade; UV-Vis: energia + força de oscilador),
alargado individualmente (Lorentziana pro IV — FWHM 10 cm⁻¹, mesmo
padrão de `plot_ir.py`; Gaussiana pro UV-Vis — σ 0.4 eV, mesmo padrão
de `parse_excitations.py`) e multiplicado pelo peso de população
daquele confôrmero antes de somar no envelope final.

**⚠️ Pegadinha real, sistemática, não um caso isolado — modos
vibracionais imaginários:** checando os 327 confôrmeros, **100% deles**
mostram modos vibracionais imaginários genuínos (não ruído numérico
perto de zero — chegam a -273 cm⁻¹), tipicamente ~6 por confôrmero.
Causa: a geometria foi otimizada em r2SCAN-3c com convergência
afrouxada (Estágio 3, teto de 10 ciclos), mas a Hessiana foi calculada
num nível diferente (PBE0/def2-SVP) que nunca fez parte dessa
otimização — a estrutura não é um ponto estacionário de verdade na
superfície do PBE0. Nos dados observados, os modos de
translação/rotação genuínos saem bem perto de 0.00 cm⁻¹, e vêm
**depois** dos modos imaginários na ordenação crescente de frequência
do Turbomole — ou seja, a convenção já existente em `parse_ir.py` de
"descartar os 6 primeiros modos por posição" **não** separa
corretamente translação/rotação de modo imaginário aqui.
`build_population_weighted_spectrum.py` filtra por **valor**, não por
posição: só mantém modos com número de onda acima de
`--min-wavenumber` (padrão 10 cm⁻¹, calibrado pelos dados desta
molécula — ajuste se a menor vibração real esperada for mais baixa que
isso).

**Decisão do usuário (2026-08-29)**: descartar os modos imaginários
(não modelar como banda em frequência negativa — não corresponde a
absorção física real) e seguir em frente, documentando como limitação
conhecida da etapa de triagem. **Isso não deveria se repetir no
espectro de produção final** (Estágio 6), que reotimiza o top-N
selecionado com os critérios de convergência rígidos originais do
Paper 1 num único nível de teoria consistente (geometria e Hessiana no
mesmo funcional/base) — se aparecer de novo lá, é motivo pra investigar
a fundo antes de aceitar o resultado, não pra aplicar o mesmo
descarte automaticamente.

**⚠️ Pegadinha corrigida — código de saída do Turbomole não é
confiável para detectar falha:** ao rodar `build_population_weighted_spectrum.py`
sobre os 327 confôrmeros, 1 (de 327) não tinha `vibspectrum` —
investigando, o `ridft` daquele confôrmero imprimiu **"ATTENTION:
ridft did not converge!"** mas ainda assim terminou com "ridft ended
normally" e código de saída `0`; o `aoforce`/`escf` rodados em cima
dessa SCF quebrada também imprimiram "ended abnormally" internamente
enquanto **também** retornavam código de saída `0`. Um script de job
que só checa `$?` depois de cada binário do Turbomole (como
`job_pbe0_ir_uvvis.sh` e `job_r2scan3c_optimize.sh`) pode dar como
sucesso um cálculo que na verdade falhou — **sempre confira também o
texto de saída** (`grep -l "ended abnormally" */aoforce.out` etc.) numa
rodada em lote, não confie só no código de saída agregado do SLURM.
Correção usada: mesmo remédio do Estágio 3 (`$scfiterlimit` maior,
`$scfdamp start` mais forte — coincidiu de novo com um gap HOMO-LUMO
quase nulo, ~0.05 eV, sintoma típico de oscilação de SCF entre
orbitais quase degenerados). Depois de um `aoforce` abortado, o
Turbomole deixa um marcador `$actual step` órfão no `control` que
bloqueia uma nova tentativa direta ("CONTRL dead = actual step") — rode
`actual -r` na pasta antes de tentar de novo.

---

## 9. Estágio 5 (em andamento desde 2026-08-31) — Camada B: aprendizado ativo

- **Dois surrogates independentes de ML**, rodando sobre o ensemble já
  filtrado pela Camada A (não o ensemble bruto do CREST — decisão
  2026-08-28, ver justificativa detalhada no histórico do projeto):
  - **B1** — desvio da banda carbonila: alvo caro = frequência
    harmônica ν(C=O) via `aoforce` completo no nível de produção;
    proxy barato a testar entre GFN2-xTB e r2SCAN-3c (teste de
    correlação antes de comprometer orçamento de Hessianas).
  - **B2** — caráter da excitação S1: alvo caro = TD-DFT em
    HSE06/def2-TZVP (mesma janela de 60 estados do Paper 1); proxy
    barato = sTDA (não GFN2-xTB — proxy errado para excitações).
- **Modelo**: GPR, kernel Matérn 5/2, delta-learning — decisão já
  tomada, não trocar sem aviso prévio (ver regras de parada abaixo).
- **Métrica de sucesso**: Spearman ρ e Recall@10 contra um subconjunto
  rotulado por DFT held-out — **não** MAE/RMSE (ver contraste com o
  Estágio 8.3: aqui o objetivo é ranking, lá seria acurácia absoluta —
  são propósitos diferentes, não a mesma tarefa com nomes diferentes).
- **Aprendizado ativo**: usa a incerteza do próprio GPR para escolher
  quais confôrmeros rotular a seguir, repete até a métrica de ranking
  estabilizar — critério de parada exato ainda não definido, precisa
  ser mostrado ao usuário (curva de convergência) antes de ser fixado.

### 9.1. Preparação de ambiente (concluída 2026-08-31)

**sTDA (proxy barato do B2) precisou ser compilado do zero.** O
repositório original (`grimme-lab/stda`) foi descontinuado e sucedido
pelo `std2` (mesmo grupo Grimme, mesmo método sTDA/sTD-DFT, só
rebatizado) — sem pacote conda/pip pronto, só código-fonte Fortran.

```bash
# ferramentas de build (uma vez, dentro do ambiente paper2_funnel)
pip install meson ninja

# a lib openblas do ambiente só tinha o .so versionado — o meson
# precisa do nome sem versão pra achar a biblioteca
ln -sf ~/.conda/envs/paper2_funnel/lib/libopenblas.so.0 \
       ~/.conda/envs/paper2_funnel/lib/libopenblas.so

git clone https://github.com/Theoretical-Chemistry-Group-UCLouvain/std2.git ~/software/std2
cd ~/software/std2
export LIBRARY_PATH="$CONDA_PREFIX/lib:$LIBRARY_PATH"
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"
meson setup _build -Dla_backend=openblas   # baixa o libcint sozinho
meson compile -C _build

ln -sf ~/software/std2/_build/std2 ~/local_bin/std2
```

Suíte de testes embutida (`tests/run_tests.sh`) validada: os testes de
energia de excitação sTDA (`EXCI_TEST1/2/3`, o que o B2 realmente usa)
passam limpo; os testes de hiperpolarizabilidade/absorção de dois
fótons falham numericamente — não usados aqui, mas não confie nesse
build para essas propriedades sem investigar mais.

**⚠️ Nota**: `~/local_bin` não está no `PATH` desta máquina (apesar do
`MANUAL.md` sugerir que devia estar, pro Multiwfn) — use o caminho
absoluto (`~/local_bin/std2` ou `~/software/std2/_build/std2`) em
qualquer script novo, não assuma que `std2` resolve sozinho.

### 9.2. Seleção do subconjunto inicial a rotular (concluído 2026-08-31)

**Script:** `select_active_learning_subset.py`

**Estratégia** (dois critérios combinados, não só um):
1. **Representatividade** — os `--n-population-seed` confôrmeros de
   maior população (padrão 3) entram sempre. São os que mais pesam na
   propriedade final ponderada, então o surrogate nunca deveria estar
   extrapolando para eles.
2. **Diversidade** — o resto das vagas (até `--n-target`, padrão 20,
   dentro da faixa 15–30 decidida para o projeto) é preenchido por
   amostragem de ponto mais distante (*farthest-point sampling*) sobre
   RMSD estrutural (alinhamento de Kabsch), a partir dos já
   selecionados: a cada passo, adiciona o confôrmero cujo vizinho mais
   próximo já selecionado está mais longe. Isso importa porque um GPR
   treinado só no aglomerado de menor energia não tem informação de
   como a correção de delta-learning se comporta no resto do espaço
   conformacional — e a incerteza que o aprendizado ativo usa pra
   escolher o próximo rótulo fica pouco confiável fora desse
   aglomerado.

Se o ensemble filtrado da Camada A já for menor ou igual ao alvo
(casos do fenolato e diânion aqui), seleciona todo mundo, sem
precisar do algoritmo.

```bash
python ~/work_turbomole/scripts/select_active_learning_subset.py \
    --population-table layer_a/population_table_dft.csv \
    --conformer-jobs layer_a/conformer_jobs \
    --n-target 20 --n-population-seed 3 \
    --output layer_b/initial_labeled_subset.csv
```

**Resultado obtido:**

| Estado | Selecionados / disponíveis | População coberta |
|---|---|---|
| Neutro | 20 / 170 (3 população + 17 diversidade) | 26.5% |
| Monoânion A | 20 / 143 (3 população + 17 diversidade) | 26.6% |
| Monoânion B | 11 / 11 (todos) | 100% |
| Diânion | 3 / 3 (todos) | 100% |

### 9.3. Rotulagem de produção (em andamento desde 2026-08-31)

**Objetivo:** gerar o "alvo caro" de verdade pro B1 (Hessiana) e B2
(TD-DFT) — HSE06-D3(BJ)/def2-TZVP, exatamente os parâmetros de
produção do Paper 1, em single-point sobre a geometria já otimizada em
r2SCAN-3c (sem reotimizar geometria nesse nível — isso fica reservado
pro Estágio 6, só nos confôrmeros finais selecionados).

**Scripts:**

```bash
# 1. Template por estado (HSE06/def2-TZVP, RI-J só — não precisa de
#    RI-K como o PBE0 da Camada A, HSE06 é híbrido blindado)
python ~/work_turbomole/scripts/setup_hse06_tzvp_template.py \
    --xyz structure.xyz --charge <carga> --nstates 60 \
    --output-dir layer_b/hse06_template

# 2. Clona o template só pro subconjunto selecionado (9.2), não pra
#    todo o ensemble da Camada A
python ~/work_turbomole/scripts/fanout_labeled_subset_jobs.py \
    --template layer_b/hse06_template \
    --subset layer_b/initial_labeled_subset.csv \
    --r2scan3c-jobs layer_a/conformer_jobs \
    --output-dir layer_b/labeled_jobs

# 3. Array SLURM: ridft -> aoforce -> escf
sbatch --array=1-54%20 ~/work_turbomole/scripts/job_hse06_labeling.sh \
    <arquivo_com_uma_pasta_por_linha>
```

**⚠️ Pegadinha corrigida — `auxbasis` não deve ser exigido no
template:** ao copiar o padrão do fan-out da Camada A (PBE0/RI-K, que
gera `auxbasis` explicitamente via `define`), o primeiro rascunho desse
script recusava rodar porque o template HSE06 não tinha `auxbasis` —
mas isso é esperado e não é erro: o Turbomole gera esse arquivo
sozinho, na primeira chamada de `ridft`/`aoforce`, quando `$rij` está
ligado mas ainda não existe um arquivo correspondente (confirmado
batendo com o comportamento real do r2SCAN-3c, cujo template também
nunca teve `auxbasis` pronto). Só métodos que passam pelo menu
`rijk`/`jkbas` do `define` (RI-K, ex. PBE0 da Camada A) geram o arquivo
de antemão — não generalize essa exigência pra todo template.

Job 5056026 submetido 2026-08-31 (54 confôrmeros, `--array=1-54%20`,
6h de teto por tarefa, sem estimativa formal de custo antes de
submeter — decisão explícita do usuário para agilizar).

**⚠️ 6h não foi suficiente.** Checando no dia seguinte: 40 das 54
tarefas morreram por `TIMEOUT` (não por erro — bateram o teto de 6h
com o `aoforce` ainda no meio da Hessiana). O `ridft` sozinho já levou
41min num dos casos investigados (bem mais que a referência de ~35s
que eu tinha, que era de um cenário diferente) — HSE06/def2-TZVP com
60 estados de TD-DFT em 8 núcleos é mais caro do que o esperado.
Correção: `actual -r` em cada uma das 40 pastas (limpa o marcador
`$actual step` órfão deixado pelo `SIGTERM` do timeout — mesma
pegadinha da Seção 8.4, mesmo remédio), resubmetidas com 20h de teto
(`job_hse06_labeling_retry.sh`, job 5059270, `--array=1-40%15`),
aproveitando o checkpoint `restarthess` que o `aoforce` deixou salvo em
vez de recomeçar do zero.

**As outras 14 (índices 41-54) também não terminaram** — checando de
novo pouco depois, todas ainda estavam `RUNNING`, 7 delas a ~20min do
próprio teto de 6h. Cancelei antes que estourassem (sem sentido deixar
rodar até o fim quando o padrão já mostrado pelas outras 40 indicava
que não iam terminar a tempo mesmo), limpei os marcadores, e
resubmeti com um script melhorado
(`job_hse06_labeling_retry2.sh`): **16 núcleos em vez de 8** (nós do
cluster têm 56 núcleos cada, sobra espaço) e variáveis de ambiente
OpenMP explícitas (`OMP_NUM_THREADS` igual ao `PARNODES` — nunca maior,
senão reproduz o bug de oversubscription do CREST/OpenBLAS visto antes
—, `OMP_STACKSIZE=4G`, `OMP_PLACES=cores`, `OMP_PROC_BIND=close`) pros
binários `*_omp` do Turbomole aproveitarem melhor os núcleos. Como o
lote de 40 (job 5059270) ainda estava quase todo pendente (36/40 sem
nem ter começado), cancelei e resubmeti ele também com o script
melhorado, em vez de deixar rodando com a configuração mais lenta —
perda de progresso desprezível. Estado final: dois arrays com o
script de 16 núcleos — job 5059893 (14 tarefas) e job 5059916 (40
tarefas), 54 no total.

### 9.4. Proxy barato do B1 — decisão forçada, sem teste de gate-zero

O desenho original previa testar GFN2-xTB vs. r2SCAN-3c como baseline
do delta-learning do B1, escolhendo o que melhor correlacionar com a
frequência real. **Não dá pra fazer esse teste**: r2SCAN-3c não
calcula frequência vibracional nenhuma neste Turbomole (mesma
limitação de meta-GGA da Seção 8.1 — `aoforce` aborta com `Invalid
value of nfun`), então só sobra um candidato viável.

**Decisão (2026-08-31)**: usar GFN2-xTB diretamente como baseline do
B1, sem teste comparativo — não é uma escolha ativa entre opções, é
a única opção que sobrou depois da limitação do r2SCAN-3c. Documentado
aqui para não parecer uma etapa pulada por acidente.

### 9.5. Próximos passos (ainda não implementados)

1. ~~Selecionar o subconjunto inicial a rotular~~ ✅ (9.2)
2. 🔄 Gerar inputs de produção e rodar (9.3) — em andamento, 40/54
   precisaram de resubmissão por timeout (ver pegadinha na 9.3).
3. ~~Teste de gate-zero do B1~~ ✅ decidido sem teste (9.4) — GFN2-xTB
   direto.
4. Rodar GFN2-xTB single-point em cada confôrmero do subconjunto
   rotulado (a mesma frequência que já vem do CREST/GFN2-xTB da
   Camada A não serve — precisa da frequência vibracional, não só a
   energia; ainda não calculada).
5. Implementar o pipeline de GPR/Matérn-5/2 + aprendizado ativo (código
   novo, ainda não escrito).

**Nota de operação**: `batch_gfn2xtb_co_stretch.py` também foi rodado
por engano direto na máquina de login antes de ir pra fila — mesmo
erro do Estágio 9.3, mesmo remédio: `job_gfn2xtb_co_stretch.sh`
(job 5060530, 5h de teto, partição `short`).

---

## 10. Estágio 6 (futuro, não implementado) — Seleção final e produção

- Camada A entrega diretamente os sobreviventes do filtro de
  população — não precisa de corte adicional.
- Camada B entrega **duas listas top-N separadas** (uma por surrogate,
  B1 e B2 podem apontar para confôrmeros diferentes — não force uma
  lista única).
- Só os confôrmeros finalmente selecionados (N entre 5–10 por estado)
  vão para reotimização de produção **usando exatamente os mesmos
  parâmetros de convergência do Paper 1** (HSE06-D3(BJ)/def2-TZVP,
  grid m4, `$senex`, `$rij`, `$scfconv 8`, `$denconv 1d-7`,
  `threchange 1e-7`, `thrmaxgrad 5e-4`, `thrrmsgrad 2e-4`,
  `thrmaxdispl 5e-4`, `thrrmsdispl 2e-4`) — para manter os dois papers
  comparáveis. **Não afrouxar esses critérios** "pra ficar mais
  rápido" sem aviso explícito, mesmo que o resto do funil use critérios
  mais soltos (Estágio 3.4) — a comparabilidade com o Paper 1 depende
  disso.

---

## 11. Tabela de referência rápida

| Estágio | Script | Quando usar |
|---|---|---|
| 0 | `build_protonation_state.py` | Uma vez por estado de protonação novo (detecta hidroxilas, remove H, gera carga) |
| 1 | `job_crest_search.sh` | Uma vez por estado, busca conformacional completa (SLURM) |
| 2 | `filter_boltzmann_population.py` | Depois do CREST, filtra por população de Boltzmann (com degenerescência) |
| 3.1 | `setup_r2scan3c_template.py` | Uma vez por estado (carga), monta o template de DFT barato |
| 3.2 | `fanout_conformer_jobs.py` | Clona o template pra cada confôrmero sobrevivente do filtro |
| 3.3 | `job_r2scan3c_optimize.sh` | Array SLURM, otimização r2SCAN-3c de cada confôrmero |
| 3.4 | `rerank_layer_a_dft.py` | Depois do array, reranqueia população por energia DFT |
| 4.1 | `setup_pbe0_svp_template.py` | Uma vez por estado, monta o template de IV/UV-Vis (PBE0/def2-SVP) |
| 4.2 | `fanout_pbe0_jobs.py` | Clona o template usando a geometria já otimizada em r2SCAN-3c |
| 4.3 | `job_pbe0_ir_uvvis.sh` | Array SLURM: SCF + Hessiana + TD-DFT de cada confôrmero |
| 4.4 | `build_population_weighted_spectrum.py` | Depois do array, monta o espectro IV/UV-Vis ponderado por população — entregável final da Camada A |
| 5 | *(não implementado)* | Aprendizado ativo B1/B2 sobre o ensemble filtrado pela Camada A |
| 6 | *(não implementado)* | Reotimização de produção do top-N final, parâmetros idênticos ao Paper 1 |

---

## 12. Solução de problemas / pegadinhas conhecidas

| Sintoma | Causa provável / correção |
|---|---|
| Logs do CREST de centenas de MB, cheios de "OpenBLAS Warning: Detect OpenMP Loop" | `OMP_NUM_THREADS` setado junto com o `-T` do CREST — oversubscription de threads. Use `OMP_NUM_THREADS=1` + `OPENBLAS_NUM_THREADS=1`, deixe só o `-T` do CREST controlar concorrência. |
| População de Boltzmann de um confôrmero parece baixa demais comparado ao que o CREST reporta internamente | Esqueceu de pesar por degenerescência (`cre_members`) — cada linha de `crest.energies` já é um confôrmero deduplicado, não um microestado de peso 1. |
| `define` diz "NO ATOMS, NO MOLECULE, NOTHING!" mesmo com um `coord` válido | Contagem errada de prompts em branco no início do `define` (são dois separados: "ler defaults de outro control?" e "título", não um só). |
| `ridft`/`jobex` recusa rodar com "Option $rij not found!" | `$rij` não foi adicionado ao `control`, ou foi adicionado depois do `$end` (Turbomole ignora tudo depois disso). |
| `ridft` recusa com "Problem reading basis set(s)" ao tentar RI-K | Base não-padrão (ex. `def2-mTZVP` do r2SCAN-3c) sem cobertura completa na biblioteca RI-JK — use uma base padrão como `def2-SVP` para etapas que precisem de RI-K. |
| `aoforce`/`escf` aborta com "Invalid value of nfun in \<mgga_r2\>!" | Limitação real de implementação: derivadas analíticas (Hessiana, resposta TD-DFT) não estão disponíveis para funcionais meta-GGA (r2SCAN, r2SCAN-3c) neste Turbomole. Troque para um híbrido/GGA convencional (ex. PBE0) só para essa etapa, mantendo a geometria já otimizada. |
| `jobex` não declara convergência mesmo depois de muitos ciclos, mas a energia já está estável | Provável modo geométrico muito plano (ex. rotação de grupo flexível) — considere um teto fixo de ciclos como critério de parada para etapas de triagem (não para produção final). |
| `scontrol update jobid=... TimeLimit=...` recusa com "Access/permission denied" | Usuário comum não pode estender o walltime de um job já rodando neste cluster — precisa cancelar e resubmeter com mais tempo (ideally usando o `crest.restart`/checkpoint equivalente, se existir, para não perder o progresso). |
| `squeue` mostra só "1 PENDING" pra um array job com muito mais tarefas pendentes | `squeue` comprime ranges de tarefas pendentes contíguas numa linha só por padrão — use `squeue -r` (expande) para uma contagem confiável. |
| `vibspectrum` mostra vários modos com número de onda bem negativo (não perto de zero) | Modo vibracional imaginário genuíno — geometria não é ponto estacionário no nível em que a Hessiana foi calculada (comum quando geometria e Hessiana vêm de níveis de teoria diferentes, ex. Estágio 4). Filtre por valor (número de onda mínimo), não por posição — a convenção "primeiros 6 modos" pressupõe que não há modo imaginário antes da translação/rotação. |
| Array SLURM reporta sucesso (`sacct` mostra exit `0:0`) mas um arquivo de saída esperado (ex. `vibspectrum`) não existe | O Turbomole nem sempre retorna código de saída diferente de zero quando um passo interno falha ("ridft did not converge", "force ended abnormally" etc. podem coexistir com exit code 0) — confira o texto de saída (`grep -l "ended abnormally"`), não só o código de saída, ao validar um lote grande. |
| `ridft`/`aoforce` recusa reiniciar com "CONTRL dead = actual step" | Marcador `$actual step` órfão deixado por uma falha anterior (ex. `force` abortado) — rode `actual -r` na pasta antes de tentar de novo. |

---

## 13. Status desta investigação (atualizado 2026-08-31)

- ✅ Ambiente configurado e validado (Seção 3).
- ✅ Geometrias de partida dos 4 estados de protonação geradas (Estágio 0).
- ✅ Busca CREST completa nos 4 estados (Estágio 1).
- ✅ Filtro de população de Boltzmann nos 4 estados, com correção de
  degenerescência (Estágio 2).
- ✅ Refino r2SCAN-3c dos 327 sobreviventes + reranking por DFT
  (Estágio 3).
- ✅ Espectros IV/UV-Vis em PBE0/def2-SVP dos 327 sobreviventes
  (Estágio 4.1–4.3) — array SLURM 5042994 terminado, 326/327 sem
  intervenção, 1 corrigido manualmente (falha silenciosa de SCF não
  sinalizada pelo código de saída — ver Seção 8.4/12).
- ✅ Espectro ponderado por população (IV e UV-Vis) montado para os 4
  estados (Estágio 4.4), com a ressalva documentada dos modos
  imaginários (Seção 8.4). **Camada A está completa** para os 4
  estados de protonação.
- ⬜ Não iniciado: aprendizado ativo B1/B2 (Estágio 5).
- ⬜ Não iniciado: seleção final e reotimização de produção (Estágio 6).

**Nota de operação (2026-08-31)**: retentativas manuais de cálculo
único (ex. corrigir 1 confôrmero problemático) devem sempre ir para a
fila SLURM (`sbatch`), mesmo sendo "só um"— não rodar direto na
máquina de login/interativa, mesmo em background. Foi um erro cometido
nesta investigação (corrigido depois de aviso do administrador do
cluster).
