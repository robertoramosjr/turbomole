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
4. ✅ Frequências GFN2-xTB do subconjunto rotulado (9.6).
5. Implementar o pipeline de GPR/Matérn-5/2 + aprendizado ativo (código
   novo, ainda não escrito) — **próximo passo real**.

### 9.6. Identificação do modo C=O e dados de treino do B1 (concluído 2026-09-01)

**Problema**: com ~130-140 modos vibracionais por confôrmero, pegar "a
banda mais intensa numa janela de frequência" não é confiável — em
particular nos estados aniônicos, a carboxila vira carboxilato e o
padrão de bandas na região muda (mais de um modo forte perto um do
outro).

**Solução**: `identify_co_stretch_mode.py` projeta o deslocamento
relativo (O menos C) de cada modo no vetor da própria ligação C=O
(mesma detecção geométrica do Estágio 0, sem índice hardcoded) — o
modo com maior projeção é o estiramento C=O. Testado e validado:
sempre um vencedor claro, com margem de 6-7x sobre o segundo colocado.

Isso precisa dos vetores de deslocamento normal, não só
frequência+intensidade — formato Gaussian-98 (`g98.out`), que tanto o
`xtb --hess` quanto o script `aoforce2g98` do próprio Turbomole (
`~/software/TURBOMOLE/scripts/aoforce2g98 aoforce.out > arquivo_g98.out`)
produzem, então o mesmo parser serve pros dois lados do delta-learning.

**Scripts:**

```bash
# proxy barato (GFN2-xTB), local, rápido -- mas vai pra fila mesmo assim
# (pegadinha de operação: rodar direto na máquina de login, mesmo pra
# algo "rápido", já rendeu um aviso do administrador nesta investigação)
sbatch ~/work_turbomole/scripts/job_gfn2xtb_co_stretch.sh

# alvo caro (HSE06), a partir das Hessianas de produção já rodadas (9.3)
python ~/work_turbomole/scripts/batch_hse06_co_stretch.py \
    --subset layer_b/initial_labeled_subset.csv \
    --labeled-jobs layer_b/labeled_jobs \
    --r2scan3c-jobs layer_a/conformer_jobs \
    --gfn2xtb-table layer_b/gfn2xtb_co_stretch.csv \
    --output layer_b/b1_training_data.csv
```

**Resultado — a suposição do delta-learning se confirma com dado
real** (não só assumida, como a "armadilha conhecida" original deste
projeto pedia pra testar antes de investir no loop inteiro):

| Estado | Correção média (HSE06 − GFN2-xTB) | Espalhamento entre confôrmeros |
|---|---|---|
| Neutro | +28.22 cm⁻¹ | 2.40 cm⁻¹ |
| Monoânion A | −79.92 cm⁻¹ | 3.46 cm⁻¹ |
| Monoânion B | +34.68 cm⁻¹ | 1.51 cm⁻¹ |
| Diânion | −65.40 cm⁻¹ | 2.42 cm⁻¹ |

Espalhamento pequeno frente à média em todos os estados → a correção é
quase um deslocamento constante, fácil de aprender com poucos pontos.
Cada estado tem um deslocamento sistemático bem diferente dos outros —
reforça que o modelo tem que ser treinado separado por estado (já era
o desenho, agora com evidência).

### 9.7. GPR do B1 — "iteração 0" (bootstrap, sem aprendizado ativo ainda)

**Scripts:** `compute_b1_features.py` (distância O···H),
`compute_b1_dihedral_feature.py` (diedro de conjugação O=C-Cα=Cβ do
braço acrílico — confirmado por geometria: Cα-C(carbonila) 1.465 Å,
Cα=Cβ 1.353 Å), `train_b1_gpr.py` (GPR Matérn 5/2 + validação cruzada
leave-one-out, apropriada pro tamanho pequeno da amostra — N=11 a 20).

**⚠️ Pegadinha corrigida — ângulo bruto tem descontinuidade:** a
primeira tentativa usou o diedro bruto (`oc_dihedral_deg`, −180° a
+180°) como feature direto, e isso piorou o modelo — dois confôrmeros
quase idênticos (ex. +179.9° e −179.93°, ambos quase totalmente
planares) ficam parecendo "máximo diferentes" pro GPR nessa
representação, por causa do salto na borda ±180°. Corrigido usando só
`planarity_deviation_deg` (dobrado pra faixa 0–90°, sem
descontinuidade) — a mesma informação de conjugação, sem o artefato.

**Resultado (só a distância O···H + a torção de planaridade, sem
aprendizado ativo ainda):**

| Estado | Spearman ρ (LOO) | MAE vs. média (baseline) |
|---|---|---|
| Neutro | −0.276 (não significativo, p=0.24) | empata com a média |
| Monoânion A | **−0.973** (muito forte, p<0.001) | pior que a média |
| Monoânion B | **0.773** (significativo, p=0.005) | melhor que a média |
| Diânion | N=3, insuficiente pra validar | — |

Achado interessante: tanto a distância O···H quanto o diedro de
conjugação variam muito pouco entre os confôrmeros do subconjunto
rotulado (diedro sempre <5° de desvio da planaridade em todo mundo) —
plausível quimicamente (torcer o braço acrílico conjugado custa
energia, então o CREST não encontra confôrmeros de baixa energia
torcidos ali; a diversidade conformacional real da molécula vem mais
das cadeias prenila do que dessa região). Isso limita quanto essas
duas features conseguem explicar sozinhas, principalmente no estado
neutro (onde a própria correção já varia muito pouco, só 2.4 cm⁻¹ no
total — pode ser que não sobre sinal real pra aprender além da média).

**Nota sobre a métrica**: o Spearman ρ do monoânion A é forte mesmo com
MAE pior que a média — isso não é contraditório, porque a métrica de
sucesso definida para o B1 é ranking (Spearman/Recall@10), não erro
absoluto (ver Estágio 9). Um modelo pode acertar a ordem relativa e
ainda assim comprimir os valores absolutos por causa da regularização
do GPR com poucos pontos.

Isso ainda é só o bootstrap ("iteração 0") — nenhum aprendizado ativo
de verdade rodou ainda (nenhum novo confôrmero foi rotulado com base na
incerteza do modelo). Próximo passo real: estender o proxy GFN2-xTB +
features pro ensemble completo da Camada A (não só o subconjunto
rotulado), pra poder ranquear todo mundo e começar o loop de fato —
ainda não feito, precisa de decisão do usuário sobre escopo/custo antes
de rodar mais confôrmeros.

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

### 9.8. Auditoria do surrogate B1 (2026-09-01)

O usuário pediu, com razão, um diagnóstico rigoroso antes de confiar no
resultado da Seção 9.7.

**Achado matemático confirmado empiricamente**: um controle "prever a
média dos outros pontos de treino" sob validação leave-one-out dá
**ρ de Spearman = −1.000 exato, sempre**, pra qualquer conjunto de
dados — consequência direta da álgebra do LOO (a previsão vira função
linear decrescente do valor verdadeiro), não coincidência. Confirmado
nos 4 estados via `diagnose_b1_surrogate.py`.

**Com um kernel GPR devidamente regularizado (ARD — length-scale por
feature — + `ConstantKernel` explícito pra variância de sinal +
`WhiteKernel` pro ruído), o Q² saiu negativo nos 4 estados**:

| Estado | ρ controle | ρ modelo (ARD) | Q² |
|---|---|---|---|
| Neutro | −1.000 | −0.007 | −0.376 |
| Monoânion A | −1.000 | −0.114 (era −0.973 com o kernel isotrópico simplificado da Seção 9.7) | −0.190 |
| Monoânion B | −1.000 | +0.509 | −0.061 |
| Diânion (N=3) | −1.000 | −1.000 (idêntico ao controle) | −1.250 |

O diânion mostra o artefato na forma mais pura. **O ρ=−0.973 do
monoânion A reportado na Seção 9.7 não se sustenta** — era o mesmo
artefato de reversão à média, manifestado através de um kernel mal
regularizado, não sinal real.

**Auditoria de features** (`audit_b1_features.py`): bond order
(`parse_bond_order.py`) e carga Bader (`parse_bader.py`) do
átomo/ligação específico **não existem** pra nenhum dos 54 confôrmeros
rotulados (os jobs de produção só rodaram `ridft+aoforce+escf`, sem
etapa de população). Mas a distância O···H **sozinha** já tem
correlação de Pearson significativa no monoânion B (r=0.839, p=0.001)
— o GPR de 2 features não está aproveitando esse sinal, possível
super-parametrização com N=11. Neutro e monoânion A não mostram
correlação individual significativa em nenhuma das duas features
atuais.

**Piso de ruído numérico** (`job_noise_floor_test.sh`, grid m4→m5,
scfconv 8→9, denconv 1d-7→1d-8, mesma geometria): desvio de
0.05–0.24 cm⁻¹ entre as duas convergências — 10-20x menor que o
espalhamento real de Δ entre confôrmeros (1.5–3.5 cm⁻¹). **O sinal que
o B1 tenta aprender não é ruído numérico disfarçado** — existe de
verdade, acima do piso. A falha do GPR aponta pra features
insuficientes (ou super-parametrização), não pro alvo ser inaprendível
em princípio.

**Status**: aguardando decisão do usuário sobre adicionar bond
order/Bader ao vetor de features (exigiria uma etapa barata de
população — `dscf -proper`/`ridft -proper`+`bader` sobre o `mos` já
convergido, não uma Hessiana nova) antes de tentar o modelo conjunto
multi-estado (explicitamente adiado a pedido do usuário) ou o loop de
aprendizado ativo de verdade.

### 9.9. Modelo de produção do B1 para o monoânion B — FECHADO (2026-09-09)

Decisão final, não candidato — não reabrir com GPR a menos que N cresça
substancialmente além dos 11 pontos atuais.

**Modelo**: regressão linear simples (OLS, 2 parâmetros) só com a
distância O···H, ajustada nos 11 pontos completos (não a média dos
folds LOO):

```
delta_cm1 = 1392.13 * oh_distance_angstrom + (-3114.30)
```

(faixa de `oh_distance_angstrom` nos dados: 2.2615–2.2623 Å; desvio
residual do ajuste completo: 0.283 cm⁻¹)

**Validação**: Q²_LOO = 0.594 (positivo, robusto — validado sob o
mesmo protocolo leave-one-out usado no GPR). Correlação de Pearson
estável sob remoção de qualquer ponto individual (r entre 0.787 e
0.888, nenhum salto grande). Nenhum ponto excede o limiar de Cook's
distance 4/n=0.364 (máximo observado: 0.267).

**Registrado explicitamente — não usar**:
- **Combinação das 2 features** (O···H + planaridade): Q²=0.489, PIOR
  que a distância O···H sozinha (0.594) — a feature de planaridade
  dilui o sinal, não soma.
- **Planaridade isolada**: Cook's distance do ponto índice 2 = 1.44,
  acima do limiar 4/n — o pouco sinal aparente (r=−0.306, já não
  significativo) depende demais de um único confôrmero pra ser
  confiável.

**Motivo de não usar GPR aqui**: um OLS de 2 parâmetros já supera o
GPR ARD de 4+ hiperparâmetros (Q²=−0.061) por uma margem grande — a
falha do GPR foi super-parametrização pra N=11, não falta de sinal.
Ver Seção 9.8 pro diagnóstico completo e `test_b1_linear_baseline.py`
pro script que produziu esse resultado.

### 9.10. Etapa de população (bond order + Bader) — em andamento (2026-09-09)

**⚠️ Correção da estimativa de custo da Seção 9.8**: a referência usada
antes (`proper_d3.out`, ~10s) era do passo errado — a geração real do
cubo de densidade pro Bader é outro arquivo
(`badergrid/proper_badergrid.out`, ~1min38s). Custo real por
confôrmero: `dscf -proper` (bond order, sem `$pointval`) ~16s +
`ridft -proper` (cubo) ~3min + `bader` ~25s ≈ 3.7 min. Pros 40
confôrmeros aprovados: **~10-20 núcleo-horas** (corrigido pra cima em
relação à estimativa original, ainda trivial em escala).

**⚠️ Pegadinha nova — `$pop` e `$pointval` não podem estar juntos no
mesmo `control` ao rodar os dois passos em sequência**: confirmado que,
com os dois presentes, tanto `dscf -proper` quanto `ridft -proper`
tentam a avaliação de grade (cara, ~3min) — rodar os dois assim duplica
o trabalho caro em vez de separar o barato (população) do caro (cubo).
O Paper 1 evitava isso mantendo o `control` do `badergrid/` separado do
principal. Aqui, `job_population_analysis.sh` roda `dscf -proper`
primeiro (control sem `$pointval`), só depois injeta o bloco
`$pointval` (calculado por confôrmero — 5.0 bohr de folga além da
caixa delimitadora atômica, ~0.15 bohr de espaçamento, valores
reproduzidos a partir da grade real do Paper 1) e roda `ridft -proper`.

**Scripts**: `setup_population_job.py` (monta o diretório por
confôrmero, calcula a grade), `job_population_analysis.sh` (array
SLURM: dscf → injeta pointval → ridft → bader).

**Escopo**: 40 confôrmeros aprovados (neutro 20 + monoânion A 20,
reaproveitando `mos` já convergido) — array SLURM 5103802.

### 9.11. Protocolo de seleção de features/modelo — automatizado (2026-09-09)

**Script**: `select_b1_feature_model.py` — implementa a ordem exigida
(Pearson por feature isolada → OLS 1-feature sob LOO → combinação só
se superar a melhor isolada → GPR só sob critérios explícitos, nunca
como primeira tentativa). Validado reproduzindo exatamente a decisão
já fechada do monoânion B (Seção 9.9) a partir do zero, e confirmando
que neutro/monoânion A ainda não têm nenhuma feature candidata
(Q²≤0) com as 2 features geométricas atuais — esperado, é por isso que
a etapa de população da Seção 9.10 está rodando.

Assim que os dados de bond order/Bader saírem, rodar este script com
o vetor expandido é o próximo passo — ainda não feito.

### 9.12. Resultado do vetor expandido — monoânion A fechado, neutro ainda sem sinal (2026-09-09)

**Script**: `compute_b1_population_features.py` extrai SEN (bond order,
via `pop.out`) e carga Bader (via `ACF.dat`) especificamente do
O/C da carbonila (mesma detecção geométrica dos outros scripts, índices
variam por confôrmero/estado — ex. O3-C22 no neutro, O2-C22 no
monoânion A, já que a numeração muda com a desprotonação).

**Monoânion A — modelo fechado**: a carga Bader do oxigênio carbonílico
é candidata real (r=0.606, p=0.0046, Q²_LOO=+0.277 — todas as outras 4
features, incluindo as 2 antigas, continuam Q²≤0 pra esse estado):

```
delta_cm1 = 4.7223 * bader_charge_carbonyl_o + (-73.9249)
```

(faixa de `bader_charge_carbonyl_o`: −1.603 a −1.155 e; desvio residual
do ajuste completo: 0.756 cm⁻¹). Faz sentido quimicamente: a densidade
eletrônica no O carbonílico deveria se relacionar direto com o quanto
a ligação C=O muda de caráter entre confôrmeros.

**Neutro — ainda sem sinal**, mesmo com as 5 features testadas (2
geométricas + bond order + 2 cargas Bader). Todas com Q²≤0. Não é um
resultado forçado — o protocolo (`select_b1_feature_model.py`,
Seção 9.11) não recomenda GPR aqui, e não há feature candidata pra
fechar um modelo de produção como foi feito pro monoânion A/B. Possível
explicação: o alvo do neutro já tem o menor espalhamento absoluto entre
os 4 estados (2.4 cm⁻¹ total) — pode não sobrar sinal aprendível acima
do que já foi descartado como ruído numérico (Seção 9.8, piso de
0.05-0.24 cm⁻¹) com as features tentadas até agora.

**Status do B1, por estado, ao final desta rodada**:

| Estado | Modelo | Q²_LOO |
|---|---|---|
| Neutro | nenhum ainda — sem feature candidata | — |
| Monoânion A | OLS 1-feature (carga Bader do O carbonílico) | +0.277 |
| Monoânion B | OLS 1-feature (distância O···H) | +0.594 |
| Diânion | não avaliado (N=3, fora de escopo) | — |

### 9.13. Calibração por permutação — monoânion A confirmado, neutro confirmado sem sinal (2026-09-09)

**Script**: `permutation_test_b1_features.py` — calibra o risco de
comparação múltipla do protocolo "testa 5 features, reporta a melhor"
(Seção 9.11). Embaralha o alvo Δ (mantendo as features intactas),
repete o procedimento completo 1000 vezes, constrói a distribuição
nula de "melhor Q² em 5 tentativas, só por acaso". Usa a fórmula
fechada (exata, via matriz hat/PRESS) do LOO pra regressão de 1
feature, validada batendo com o valor já obtido por refit explícito
antes de confiar nas 1000 repetições.

**⚠️ Bug próprio corrigido**: a primeira versão do script tinha um
limiar de interpretação (percentil ≤60 = "normal") arbitrário demais —
classificou erroneamente o resultado do neutro (percentil 72.4%, bem
dentro do corpo comum da nula) como "anômalo". Corrigido: só o percentil
baixo (<10, pior que quase todo ruído puro) seria de fato uma anomalia
a investigar; qualquer coisa no meio do corpo da distribuição nula é
esperada e não significa nada de errado nos dados.

**Resultado**:

| Estado | Q² observado | Percentil na nula | p-valor empírico | Limiar honesto (percentil 95) |
|---|---|---|---|---|
| Monoânion A | +0.277 | **98.3%** | 0.017 | 0.146 |
| Neutro | −0.011 | 72.4% | 0.276 | 0.139 |

**Monoânion A confirmado**: Q²=0.277 supera o limiar de percentil 95
da nula — não é explicado por acaso mesmo testando 5 features e
reportando só a melhor. Reforça a decisão da Seção 9.12 (modelo com
carga Bader do O carbonílico).

**Neutro confirmado sem sinal, de verdade**: o "zero" (nenhuma das 5
features passou de Q²=0) cai dentro do corpo comum esperado da
distribuição nula (percentil 72.4%, nem extremo alto nem baixo) — ou
seja, é exatamente o que se esperaria ver mesmo que não houvesse
nenhum sinal real pra aprender. Isso não é um artefato de features
ruins ou processo malfeito — reforça que o estado neutro genuinamente
não tem, até agora, um modelo B1 recuperável com o vetor de features
tentado.

### 9.14. B1 fechado — status final (2026-09-09)

**Neutro — fechado como offset, não descartado**: sem feature
recuperável (Seção 9.13 confirmou por permutação que a ausência de
sinal é real, não falta de tentativa), o estado neutro fecha com

```
delta_neutro = 28.222 ± 0.676 cm-1  (N=20, média ± desvio padrão)
```

Esse ± deve ser propagado como incerteza na etapa de especiação
(Ω_i(pH)), não descartado — o estado continua participando do funil
normalmente, só não tem uma correção dependente de confôrmero, e sim
um deslocamento médio com barra de erro.

**Monoânion A — faixa de validade documentada**: o modelo
`delta_cm1 = 4.7223 * bader_charge_carbonyl_o + (-73.9249)` (Seção
9.12) só vale dentro da faixa de treino:

```
bader_charge_carbonyl_o ∈ [-1.6034, -1.1553] e
```

Fora dessa faixa é extrapolação — usar `apply_b1_model.py`, que
sinaliza cada previsão como `interpolation` ou `EXTRAPOLATION` em vez
de devolver um número com a mesma confiança implícita. Dado que
Q²=0.277 já é moderado dentro da faixa treinada, uma previsão
extrapolada deve ser tratada como pouco confiável ou não usada — não
uma diferença cosmética de rótulo.

**Protocolo congelado pra reuso** (B2, calibração de pKa): ver
[feature_selection_protocol.md](feature_selection_protocol.md) — as 4
etapas completas (Pearson individual → OLS-1-feature sob LOO →
combinação só se justificada → calibração por permutação obrigatória),
mais a regra de faixa de validade pós-fechamento e o critério explícito
de quando (não) usar GPR.

**Tabela final do B1:**

| Estado | N | Modelo | Q²_LOO | Validação |
|---|---|---|---|---|
| Neutro | 20 | offset (28.222 ± 0.676 cm⁻¹) | — | permutação confirma ausência real de sinal (percentil 72.4%, p=0.276) |
| Monoânion A | 20 | Bader(O carbonílico) | 0.277 | permutação: percentil 98.3%, p=0.017 |
| Monoânion B | 11 | distância O···H | 0.594 | Cook's D (máx 0.267, abaixo do limiar) + leave-one-out do r (0.787–0.888) |
| Diânion | 3 | fora de escopo | — | N insuficiente pra qualquer validação |

**B1 está fechado.** Próximo: B2 (caráter da excitação S1/CT) — ainda
não iniciado, nenhum dado de sTDA ou TD-DFT de produção rotulado até
agora.

### 9.15. B2 — extração dos alvos e do proxy inicial (2026-09-09)

Alvo caro: S1 (energia + força de oscilador) já disponível de graça a
partir dos mesmos jobs de rotulagem de produção HSE06/def2-TZVP do B1
(`escf` já rodou com 60 estados singlete, igual ao Paper 1) —
`batch_extract_s1_hse06.py` extrai o estado de menor energia por
conformero, sem custo de DFT adicional. Rodado para os 54 confôrmeros
rotulados nos 4 estados.

Proxy barato: sTDA via `std2` sobre a função de onda r2SCAN-3c já
computada na Camada A (também sem custo adicional), convertida por
`tm2molden` (dois prompts em branco, não um — pegadinha nova) em vez
da antiga rota `wfn.xtb` do xtb (`-xtb` do std2), que não funciona mais
com as versões modernas do xtb (6.6/6.7 não geram mais esse arquivo).
`batch_stda_s1.py` roda esse pipeline; `-ax` controla a fração de
exchange exato assumida pela referência (0.0 = apropriado pra uma
referência meta-GGA pura como r2SCAN-3c).

Achado inicial: energia de S1 tem um delta bem comportado entre sTDA e
HSE06 em todos os 4 estados. Força de oscilador, não — a razão
f_HSE06/f_sTDA varia de 1.5x a 8000x dependendo do estado, claramente
fora de controle.

### 9.16. Diagnóstico da força de oscilador — não é o parâmetro `ax` (2026-09-09)

Hipótese inicial: `-ax 0.0` é inadequado porque HSE06 tem 25% de
exchange exato e o sTDA está rodando com uma referência que assume 0%.
`diagnose_stda_oscillator.py` separou erro ABSOLUTO |f_HSE06 - f_sTDA|
de erro RELATIVO (razão) — a razão explode trivialmente quando o
denominador (f_sTDA) é pequeno, então uma razão extrema não implica
necessariamente um erro absoluto grande. Resultado com `ax=0.0`:
monoânion A (ambos os métodos concordam que S1 é escuro) e diânion (a
maioria das razões extremas É artefato de denominador pequeno) estão
OK; neutro e monoânion B têm erro absoluto genuinamente grande
(médias 0.378 e 0.791, numa escala de f que vai de 0 a ~1) — não é
artefato de denominador.

Rerodando com `-ax 0.25` (igual à fração de exchange exato do HSE06)
para todos os 54 confôrmeros: **o erro absoluto não muda em nenhum dos
4 estados** (neutro 0.3778→0.3772, monoânion B 0.7911→0.7925,
monoânion A e diânion também inalterados dentro do ruído). `ax` não é
o parâmetro responsável.

**Causa raiz encontrada**: inspecionando o espectro completo do sTDA
(não só o estado 1) para confôrmeros representativos do neutro e do
monoânion B, em TODO o grid `ax` de 0.0 a 0.30 (7 valores testados), o
estado sTDA de menor energia é sistematicamente um estado ESCURO
espúrio (f~0.0001-0.02), e o estado realmente brilhante (que bate com
a transição dominante do S1 real do HSE06, mesmo par orbital
HOMO→LUMO, intensidade comparável) fica deslocado pro rank 2-4. Isso é
robusto em toda a faixa de `ax` testada — não é um artefato do
parâmetro de exchange, é uma limitação estrutural do sTDA (aproximação
monopolar) pra esse sistema: ele empurra pra baixo um estado escuro
que na referência real (HSE06 TD-DFT completo) não fica abaixo do
estado brilhante. Comparar "S1 por energia" entre os dois métodos está
comparando estados eletrônicos DIFERENTES pro neutro e monoânion B.

Consequência prática: nenhum ajuste de `ax`, nem um grid completo em
log-espaço, resolveria isso — a inversão de rank nunca desaparece em
nenhum ponto do grid testado. A via de correção tem que ser outra:
identificar o estado sTDA certo por caráter, não por posição de
energia nem (sozinho) por intensidade.

### 9.17. Seleção de estado por caráter orbital, não por energia ou intensidade (2026-09-09)

Tentativa 1 (maior f dentro de uma janela de energia): funciona melhor
que "menor energia", mas tem dois problemas próprios, achados por
`analyze_stda_state_selection.py`:

- **Sem limite de janela**, a busca por maior f é sequestrada em 2/20
  confôrmeros do neutro por um estado espúrio de energia bem alta
  (~6 eV, HSE06 real ~3.9 eV) com f grande por acaso e caráter orbital
  errado. Limitando a busca aos primeiros 10 estados por energia (não
  precisava de mais estados calculados, precisava de MENOS no
  critério de busca), isso desaparece — rank fica em 3-4 (neutro) e
  sempre 2 (monoânion B).
- **Contraste de f é fraco mesmo corrigido**: mediana 2.5x (neutro) e
  5.3x (monoânion B), nunca passa de ~10x — há quase sempre um segundo
  estado de intensidade comparável por perto. Apesar disso, o estado
  de maior f (na janela de 10) bate com o caráter HOMO→LUMO do HSE06
  em 18/20 (neutro) e 11/11 (monoânion B) — as exceções são empates
  próximos (contraste ~1.26x), não erros grosseiros.

Critério final adotado, em `batch_stda_s1_character_matched.py`:
extrai o índice HOMO local de cada rodada do sTDA a partir da linha
`oMOs in TDA: N` do stdout do `std2` (o sTDA trunca e renumera
localmente o espaço ativo de orbitais por rodada — índices absolutos
não são comparáveis entre valores de `ax`, nem com a numeração do
Turbomole), calcula o par ocupado→virtual dominante do S1 real do
HSE06 relativo ao SEU PRÓPRIO HOMO (via a linha `number of occupied
orbitals :` do `escf.out`), e escolhe, dentro da janela de 10 estados,
o estado sTDA cujo par dominante tem o MESMO offset relativo. Cai pra
"maior f na janela" só se nenhum estado bater esse offset (não
aconteceu nenhuma vez nos 54 confôrmeros rotulados testados — 0
fallbacks).

Achado lateral relevante: pra monoânion A e diânion, o S1 real do
HSE06 já é escuro por natureza (f~0.0003-0.0006) — esses dois estados
NUNCA tiveram problema de seleção de estado (o "estado 1 por energia"
original já concordava bem, daí o erro absoluto pequeno visto desde a
Seção 9.16). Rodar o critério de "maior f" nesses dois estados dá
100% de "mismatch" contra o caráter — não porque o critério de caráter
esteja errado, mas porque "maior f" está procurando um estado
brilhante que simplesmente não é o S1 real ali. **O critério de maior
f / caráter orbital só se aplica a neutro e monoânion B** (S1
brilhante); monoânion A e diânion continuam com a extração original
por energia mais baixa.

### 9.18. B2 — força de oscilador: resultado do delta-learning com rótulos corrigidos (2026-09-09)

Com os rótulos sTDA corrigidos por caráter orbital
(`build_b2_osc_delta_target.py` junta HSE06 + sTDA casado, calcula
`delta_energy_eV` e `delta_log10_f` com piso de 1e-4 antes do log,
igual à convenção de "quase-zero" já usada no diagnóstico), rodando o
protocolo completo (`select_b1_feature_model.py` e
`permutation_test_b1_features.py`, generalizados com `--target-column`
pra reuso fora do B1) sobre as MESMAS features geométricas já
calculadas pro B1 (distância O···H, diedro/planaridade de conjugação,
ordem de ligação SEN e carga de Bader do carbonilo):

| Estado | N | Alvo | Melhor feature | Q²_LOO | Permutação |
|---|---|---|---|---|---|
| Neutro | 20 | ΔE(S1) | Bader(O carbonílico) | +0.007 | percentil 65.3%, p=0.347 — sem sinal real |
| Neutro | 20 | Δlog₁₀(f) | nenhuma | ≤0 | percentil 55.4%, p=0.446 — sem sinal, consistente |
| Monoânion A | 20 | ΔE(S1) | nenhuma | ≤0 | percentil 53.7%, p=0.463 — sem sinal |
| Monoânion A | 20 | Δlog₁₀(f) | nenhuma | ≤0 | percentil 37.1%, p=0.629 — sem sinal |
| **Monoânion B** | 11 | ΔE(S1) | **distância O···H** | **+0.596** | **percentil 99.2%, p=0.008 — CONFIRMADO** |
| **Monoânion B** | 11 | Δlog₁₀(f) | **distância O···H** | **+0.412** | **percentil 98.9%, p=0.011 — CONFIRMADO** |
| Diânion | 3 | — | — | — | N insuficiente (protocolo exige N≥4) |

Achado notável: a MESMA feature que fechou o B1 do monoânion B
(distância O···H) também prevê tanto a energia quanto a força de
oscilador do S1 nesse estado — plausível fisicamente, já que essa
distância está ligada a uma ligação de hidrogênio intramolecular que
afeta a conjugação do sistema todo, não só o modo C=O. Neutro e
monoânion A fecham como "sem modelo recuperável" com as features
atuais, igual ao padrão já visto no B1 do neutro — calibrado por
permutação, não por falta de tentativa.

### 9.19. Caracterização de confiança pra produção — só onde o critério se aplica (2026-09-09)

Pros confôrmeros SEM referência HSE06 (o resto do ensemble completo da
Camada A, além do subconjunto rotulado), não há como validar caráter
orbital — só o contraste de f (sTDA sozinho) fica disponível como
sinal de confiança, e só faz sentido pra neutro/monoânion B (Seção
9.17 — monoânion A/diânion não usam esse critério). Da contagem da
Camada A completa: neutro tem 170 confôrmeros (150 ainda não
rotulados), monoânion B tem exatamente 11 — ou seja, JÁ 100% rotulado,
nada novo a caracterizar ali.

`characterize_stda_selection_confidence.py` deriva um limiar a partir
do maior contraste visto entre os 2 casos conhecidos de mismatch de
caráter no neutro (Seção 9.17): confôrmeros de produção com contraste
≤ esse limiar são marcados `uncertain_selection`. Rodado nos 170
confôrmeros do neutro via SLURM (fila `short`, ~4s/confôrmero medido
em amostra, ~12 min total — não rodado no login node).

**Resultado (170/170 processados)**: limiar derivado = 1.26x (maior
contraste entre os 2 mismatches conhecidos; nenhum dos 18 acertos
conhecidos cai nessa faixa — separação limpa nos dados rotulados).
Distribuição de contraste no ensemble completo: mínimo 1.05x, mediana
2.44x, máximo 5.72x — consistente com a amostra rotulada, sem cauda
inesperada. **12/170 (7.1%) marcados `uncertain_selection`**
(job_0002, 0011, 0012, 0070, 0072, 0095, 0097, 0100, 0152, 0153, 0154,
0162) — esses precisam de tratamento especial (incerteza inflada ou
checagem manual antes do ranking top-N final) quando o surrogate do
B2 for aplicado ao ensemble completo. Monoânion B: ensemble completo
(11) já é 100% igual ao subconjunto rotulado — nada de novo a marcar
ali. Monoânion A/diânion: critério não se aplica (Seção 9.17).

**Aplicar o surrogate treinado (Seção 9.18) aos 327 confôrmeros ainda
não foi feito** — próxima etapa, agora que o limiar de confiança está
definido.

### 9.20. Antes de aplicar aos 327: comparação antes/depois, margem do limiar, checagem dos incertos (2026-09-10)

**Antes/depois da correção de seleção de estado.** Reconstruindo o
alvo de delta com os rótulos ANTIGOS (sTDA por rank de energia,
mesmo `ax=0.25`, sem casamento de caráter) e rodando o mesmo
protocolo, para isolar o efeito da correção de seleção de estado (não
confundir com o efeito de mudar `ax`):

| Estado | Alvo | Q²_LOO ANTES (rank) | Q²_LOO DEPOIS (caráter) | Mudança |
|---|---|---|---|---|
| Neutro | ΔE(S1) | nenhuma feature (todas ≤0) | +0.007 (não sobrevive permutação) | sem mudança prática |
| Neutro | Δlog₁₀(f) | co_bond_order_sen: **+0.203**, permutação percentil 95.5%, p=0.045 (marginal) | nenhuma feature (todas ≤0), percentil 55.4%, p=0.446 | **o "sinal" ANTES desaparece DEPOIS** — era o rótulo errado (mistura de estado) correlacionando por acaso com uma feature, não física real |
| Monoânion B | ΔE(S1) | oh_distance: +0.092 | oh_distance: **+0.596** | correção aumenta o Q² em ~6.5x |
| Monoânion B | Δlog₁₀(f) | nenhuma feature (todas ≤0) | oh_distance: **+0.412**, permutação percentil 98.9%, p=0.011 | **sinal só existe DEPOIS da correção** — com rótulos errados, a força de oscilador do monoânion B não tinha nenhum sinal recuperável |

Conclusão: a correção de seleção de estado não é cosmética — ela
inverte a conclusão em 2 dos 4 casos (cria sinal real onde não havia
nenhum, e remove um "sinal" espúrio que só existia por causa do
rótulo errado).

**Margem de segurança do limiar de contraste (170 confôrmeros do
neutro, Seção 9.19):**

| Faixa | N | % |
|---|---|---|
| <1.26x (descartado) | 12 | 7.1% |
| 1.26–1.5x (intermediário) | 12 | 7.1% |
| 1.5–2.0x | 16 | 9.4% |
| >2.0x (confiável) | 130 | 76.5% |

O limiar de 1.26x veio de só 2 casos conhecidos — evidência fina, como
o Roberto apontou. Proposta (a decidir): tratar a faixa intermediária
1.26–2.0x (28/170, 16.5%) com incerteza inflada no ranking final, em
vez de uma fronteira rígida em 1.26x — os 12 casos abaixo de 1.26x
already viraram objeto de checagem direta de caráter (próximo item) em
vez de descarte automático.

**Checagem de caráter orbital dos 12 `uncertain_selection` (sem rodar
HSE06 novo)** — `check_uncertain_conformer_character.py` compara o
offset (occ−HOMO, virt−LUMO) do estado escolhido (maior f na janela)
contra o padrão de falha conhecido (offset (−2,0), visto nos 2
mismatches confirmados do treino) e verifica se existe um competidor
(0,0) (HOMO→LUMO real) por perto:

- **2/12 já são do subconjunto rotulado** (job_0011, job_0072) — são
  literalmente os 2 casos usados pra calibrar o limiar, não incógnitas
  novas. Já têm o valor certo extraído por
  `batch_stda_s1_character_matched.py`; não entram no problema de
  "confôrmero sem rótulo".
- **8/10 novos (job_0002, 0095, 0097, 0100, 0152, 0153, 0154, 0162):
  LIKELY_CORRECT** — o estado escolhido já tem offset (0,0)
  (HOMO→LUMO real); o contraste baixo é contra outro estado
  quimicamente diferente, não o mesmo padrão de falha do treino.
  Provavelmente seleção correta apesar da margem fina.
- **1/10 (job_0070): PROBABLE_MISMATCH** — offset (−2,0) igual ao
  padrão de falha conhecido, com um competidor (0,0) genuíno logo ao
  lado (f=0.294 vs f=0.309 escolhido, quase empate) — mesma assinatura
  dos 2 casos de treino.
- **1/10 (job_0012): mesmo padrão não-canônico do job_0011** — números
  quase idênticos ao job_0011 (f, E) mesmo expandindo a janela pra 8
  estados, nenhum estado (0,0) aparece em nenhum dos dois; ambos têm
  em vez disso um estado (−1,0) fraco (f~0.056, E~2.36eV). Forte
  evidência circunstancial de que job_0012 é um confôrmero quase-
  duplicado do job_0011 (mesmo tipo de geometria), herdando o mesmo
  caráter não-canônico do S1 real (um dos 2/20 casos do treino onde o
  próprio HSE06 não é HOMO→LUMO limpo) — não uma detecção nova e
  independente de erro.

**Resumo prático**: dos 12 flagados, só **1 (job_0070)** tem evidência
de ser um erro de seleção genuíno e novo; **1 (job_0012)** é
provavelmente um caso do mesmo tipo não-canônico já visto no treino
(precisaria da mesma extração especial, não da extração padrão por
maior-f); os outros 10 (8 novos + 2 já rotulados) estão resolvidos ou
provavelmente corretos apesar do contraste baixo.

### 9.21. B2 força de oscilador/energia — FECHADO e aplicado aos 327 (2026-09-10)

**Monoânion B, ΔE(S1) — calibração completa**: Q²=0.5956 (distância
O···H), permutação = percentil 99.2%, p=0.008 — CONFIRMADO, mesmo
padrão de calibração já usado no B1 e no Δlog₁₀(f) do próprio B2
(percentil 98.9%, p=0.011). **Os dois alvos do monoânion B (energia e
força de oscilador) estão fechados e confirmados.**

**Neutro — fechado como offset nos dois alvos**, mesmo princípio do
B1 (Seção 9.14): nenhuma feature sobrevive à calibração por permutação
em nenhum dos dois alvos (Seção 9.20). Offsets (N=20, média ± desvio
padrão amostral):

```
delta_energy_eV(neutro)  = 1.3206 ± 0.2140 eV
delta_log10_f(neutro)    = 0.2066 ± 0.4305 (em log10(f))
```

Essas incertezas devem ser propagadas na etapa de composição do
espectro ponderado por população, não descartadas — mesmo tratamento
já dado ao offset do B1 do neutro.

**Nota importante sobre a validade do modelo de monoânion B**: a faixa
de treino da distância O···H é extremamente estreita
(2.2615-2.2623 Å, span de 0.0008 Å) — a inclinação ajustada é
consequentemente muito íngreme (−51.6 eV/Å pra energia, −110.4 pra
log₁₀f). Isso não invalida a calibração (Q² e permutação já dão conta
do ajuste dentro da faixa), mas TORNA CRÍTICO nunca extrapolar esse
modelo além da faixa treinada. Na prática isso não é um problema para
o ranking atual: o ensemble completo do monoânion B na Camada A tem
exatamente 11 confôrmeros — o mesmo número já rotulado com HSE06 real.
**Não há confôrmero novo do monoânion B pra prever** — os 11 usam o
valor real do HSE06 diretamente (`method=hse06_ground_truth`), o
modelo de regressão fica documentado pra reuso futuro (se a busca
CREST desse estado for expandida) mas não é aplicado agora.

**Tratamento especial aplicado (Seção 9.20)**:
- `job_0070`: estado forçado por caráter (offset (0,0), E=2.841 eV,
  f=0.2942) em vez do de maior f (E=2.815 eV, f=0.3092) — assinatura
  de falha confirmada.
- `job_0012`: tratado como a mesma classe não-canônica do `job_0011`
  (offset (−1,0), E=2.355 eV, f=0.0563) em vez do de maior f.

**Aplicação aos 327 confôrmeros** — `apply_b2_osc_surrogate.py`,
rodado pro neutro (170) e monoânion B (11); monoânion A e diânion
ficam fora do escopo desta etapa (S1 real é escuro, o critério
maior-f/oscilador não se aplica, e nenhum modelo de delta sobreviveu à
calibração pra eles — Seção 9.18):

| Estado | N total | HSE06 real | Estado forçado (especial) | Surrogate confiável (contraste>2.0x) | Surrogate incerteza inflada (1.26–2.0x) | Excluído (≤1.26x, sem tratamento) |
|---|---|---|---|---|---|---|
| Neutro | 170 | 20 | 2 (job_0070, job_0012) | 118 | 22 | 8 |
| Monoânion B | 11 | 11 | 0 | 0 | 0 | 0 |

Os 8 excluídos do neutro (job_0002, 0095, 0097, 0100, 0152, 0153,
0154, 0162) são os mesmos que a Seção 9.20 classificou como
LIKELY_CORRECT (caráter já bate com HOMO→LUMO real apesar do
contraste baixo) — ficam fora do ranking por decisão explícita
(exceto os com tratamento especial confirmado), não porque haja
evidência concreta de erro; podem ser revisitados manualmente depois.

Saída final por confôrmero em
`layer_b/b2_osc_predictions_full_ensemble.csv` (ambos os estados):
energia e força de oscilador previstas, `energy_std`/`log10_f_std`
propagados, `method` (`hse06_ground_truth` / `surrogate_forced_state`
/ `surrogate_max_f`), `special_treatment`, `f_contrast_ratio`,
`uncertain_inflated`, `excluded`.

### 9.22. Reconstrução do espectro ponderado e lote de confirmação DFT (2026-09-10)

**Tarefa 1 — população conformacional**: já existe e já está em
produção — `layer_a/population_table_dft.csv` (coluna
`r2scan3c_population`, indexada por `job_dir`), a mesma tabela que
`build_population_weighted_spectrum.py` já usa pros espectros de
IV/UV-Vis do B1. Nada novo precisou ser gerado.

**Tarefa 2 — ranking de contribuição espectral** (população × f,
`rank_b2_spectral_contribution.py`): neutro tem 8.8% da população
total em confôrmeros `excluded` (sem f utilizável). Achado lateral
relevante: vários confôrmeros com `method=surrogate_max_f` no topo do
ranking do neutro têm f previsto **> 1.0** (ex. job_0151: f=1.078,
job_0032: f=1.032) — fisicamente implausível pra uma transição de
elétron único. É consequência direta de somar um offset FIXO em
log-espaço (Δlog₁₀f médio = +0.207, ou seja, ×1.61 em escala linear)
a confôrmeros cujo f de sTDA já estava no topo da faixa vista no
treino (~0.28-0.67) — o offset extrapola em ESPAÇO DE SAÍDA mesmo sem
extrapolar em espaço de feature (não há feature aqui, é um offset
fixo). Isso reforça diretamente a prioridade da Tarefa 3: são
exatamente esses confôrmeros de maior peso e sem rótulo real que mais
precisam de confirmação.

**Tarefa 3 — lote de confirmação DFT do neutro** (10 confôrmeros,
`select_b2_dft_confirmation_batch.py`, exclui os 8 `excluded` +
job_0070 + job_0012):

| job_dir | população | f previsto | motivo |
|---|---|---|---|
| job_0101 | 0.0409 | 1.017 | maior peso, offset puro |
| job_0151 | 0.0234 | 1.078 | maior peso, offset puro |
| job_0023 | 0.0229 | 1.008 | maior peso, offset puro |
| job_0116 | 0.0234 | 0.943 | maior peso, offset puro |
| job_0032 | 0.0198 | 1.032 | maior peso, offset puro |
| job_0096 | 0.0188 | 0.557 | incerteza inflada, maior peso |
| job_0010 | 0.0147 | 0.513 | incerteza inflada, maior peso |
| job_0038 | 0.0115 | 0.643 | incerteza inflada, maior peso |
| job_0079 | 0.0004 | 0.953 | diversidade forçada (RMSD 2.74 Å) |
| job_0018 | 0.0222 | 0.537 | diversidade forçada (RMSD 2.61 Å) |

Os 2 conformeros de diversidade forçada ficaram a >2.6 Å (átomos
pesados, Kabsch) do lote já selecionado — os 8 primeiros (maior
peso + incerteza inflada) formam um cluster conformacional bem mais
compacto entre si, confirmando que a diversidade forçada não era
redundante.

Monoânion B: nenhum lote — 11/11 já são HSE06 real (Seção 9.21).

**Tarefa 4 — estimativa de custo (não submetido)**: medido nos 54
confôrmeros já rotulados, `escf` (TD-DFT, 60 estados singlete,
HSE06/def2-TZVP, 16 núcleos) custa em média **~7.5 núcleo-horas**
(7h15m-7h43m de CPU em 3 amostras), mais **~0.27 núcleo-horas** do
`ridft` (SCF inicial) que precede -- nenhum cálculo de Hessian
(`aoforce`, que sozinho custa ~27 núcleo-horas) é necessário aqui, só
o TD-DFT. Geometria e `mos` de partida reaproveitados do r2SCAN-3c já
convergido (Camada A), sem otimização nova.

```
10 confôrmeros × ~7.8 núcleo-horas/confôrmero ≈ 78 núcleo-horas total
```

Escala "trivial" comparado às rotulagens anteriores (~5500
núcleo-horas os lotes grandes de Camada A). Tempo de parede esperado
~30-35 min por confôrmero com 16 núcleos, rodando em paralelo como
array SLURM (mesma configuração de `job_hse06_labeling_retry2.sh`) —
o lote inteiro roda em paralelo, não em série.

**Nada disso foi submetido.** Aguardando aprovação do lote antes de
rodar qualquer HSE06 novo.

### 9.23. Lote de confirmação DFT rodado — resultado preocupante pro offset do neutro (2026-09-10)

Aprovado e rodado (10 confôrmeros, `escf`-only, sem otimização
geométrica — reaproveitando a geometria r2SCAN-3c, como nos 54 já
rotulados; sem bloco `$optimize`, então o pedido de aumentar pra 800
passos não se aplicou aqui, fica registrado pro Estágio 6 futuro).
Pegadinha nova de infraestrutura: a lista de jobs pro array SLURM
precisa estar no filesystem compartilhado do projeto, não em `/tmp`
(`/tmp` é local por nó neste cluster — os nós de computação não
enxergam o `/tmp` do login node, causando falha silenciosa "no valid
directory" em toda tarefa do array até a lista ser movida). Custo real
ficou dentro do estimado, 10/10 tarefas `ridft`+`escf` exit 0, sem
texto de falha.

**Erro de processo próprio**: `job_0018` (escolhido como "diversidade
forçada" por `select_b2_dft_confirmation_batch.py`) já era um dos 20
confôrmeros originalmente rotulados — o filtro de diversidade não
excluiu os já-rotulados (`hse06_ground_truth`) do pool de candidatos,
só excluiu os `excluded`/tratados manualmente. Rodou de novo por
engano (~7.8 núcleo-horas desperdiçadas, resultado bateu exato com o
valor já conhecido, como esperado). Script corrigido pra excluir
`hse06_ground_truth` do pool de diversidade também.

**Comparação real vs. previsto (9 confôrmeros genuinamente novos)**:

| job_dir | E real | E previsto | ΔE | f real | f previsto | razão | confiança |
|---|---|---|---|---|---|---|---|
| job_0101 | 4.012 | 4.057 | +0.045 | 0.656 | 1.017 | 1.55x | confiável |
| job_0151 | 4.003 | 4.049 | +0.045 | 0.665 | 1.078 | 1.62x | confiável |
| job_0023 | 3.966 | 4.063 | +0.096 | **0.063** | **1.008** | **16.0x** | **confiável (contraste 3.57x)** |
| job_0116 | 3.993 | 4.059 | +0.066 | 0.604 | 0.943 | 1.56x | confiável |
| job_0032 | 3.926 | 4.168 | +0.241 | 0.315 | 1.032 | 3.28x | confiável |
| job_0096 | 3.870 | 4.163 | +0.293 | 0.335 | 0.557 | 1.66x | incerteza inflada |
| job_0010 | 3.916 | 4.134 | +0.218 | 0.419 | 0.513 | 1.22x | incerteza inflada |
| job_0038 | 3.895 | 4.176 | +0.281 | 0.282 | 0.643 | 2.28x | incerteza inflada |
| job_0079 | 4.006 | 4.048 | +0.042 | 0.616 | 0.953 | 1.55x | diversidade |

Energia: viés sistemático, todas as 9 previsões SUPERESTIMAM (ΔE
médio +0.15 eV, máx +0.293 eV) — maior que o esperado pelo desvio
padrão do offset (0.214 eV) sozinho; o offset calibrado em N=20 não
está reproduzindo tão bem nesse lote novo quanto a validação LOO
sugeria.

**Achado mais sério — `job_0023`**: estava no grupo de **maior peso
populacional, contraste alto (3.57x, bem acima do limiar de 2.0x de
"confiável")** — exatamente o tipo de confôrmero que o sistema de
confiança deveria proteger. O sTDA (mesmo com seleção de estado por
maior-f dentro da janela, contraste ok) encontrou um estado brilhante
(f~0.63) que o HSE06 real NÃO confirma — o S1 real é quase escuro
(f=0.063), erro de 16x. **Isso mostra que o contraste de f não
protege contra esse tipo de erro** — contraste alto só garante que o
sTDA tem um estado claramente dominante DENTRO DO PRÓPRIO sTDA, não
que esse estado corresponda à transição real do HSE06. É um modo de
falha diferente do que a Seção 9.17-9.20 mapeou (lá, dois estados sTDA
competiam; aqui, o sTDA está confiante mas errado).

**Consequência prática, a decidir com o Roberto**: a confiança
depositada no critério de contraste pra força de oscilador do neutro
é mais frágil do que a calibração LOO/permutação sozinha sugeria. Os
118 confôrmeros "confiáveis" (contraste>2.0x, Seção 9.21) podem conter
outros casos como o job_0023 sem que o sistema atual tenha como
detectar. Recomendação a avaliar: tratar a intensidade do espectro do
neutro (não a posição de energia, que errou bem menos) com incerteza
qualitativa muito maior do que a Seção 9.21 propunha, ou rodar mais
confirmações antes de finalizar a Camada B do neutro.

### 9.24. Separando viés de seleção de viés real, e uma nova assinatura de falha (2026-09-10)

**Tarefa 1 — simulação de viés por seleção de extremos**
(`simulate_offset_selection_bias.py`): modela o pior caso (ordem de
seleção = ordem do resíduo real, mais favorável possível à hipótese
"é só seleção") -- 5000 repetições, ruído ~N(0, 0.214):

| Grupo | Pool | N selecionado | Viés médio observado | Nulo do pior caso (média ± dp) | Percentil observado |
|---|---|---|---|---|---|
| Maior peso, offset puro | 118 | 5 | +0.0986 eV | 0.447 ± 0.051 | **0.0** (bem abaixo) |
| Incerteza inflada, maior peso | 22 | 3 | +0.2638 eV | 0.325 ± 0.078 | **22.1** (abaixo da mediana) |

**O viés observado fica ABAIXO do que a seleção por extremos, no pior
caso, poderia produzir sozinha** -- ou seja, a magnitude do viés
observado (+0.10 a +0.26 eV) é totalmente compatível com efeito de
seleção (real ou parcial), sem precisar invocar nenhum viés adicional
do offset além do já conhecido. Reforça isso: o único confôrmero do
lote que NÃO foi selecionado por valor extremo (`job_0079`, escolhido
só por diversidade geométrica) teve o menor ΔE de todos (+0.042 eV) --
consistente com "quanto mais a seleção mirou no valor previsto, maior
o viés aparente". **Conclusão da Tarefa 1: o viés de energia visto no
lote NÃO é evidência de que o offset esteja quebrado -- é
majoritariamente (possivelmente inteiramente) efeito de seleção.**
Corrige o tom de alarme da Seção 9.23 quanto à ENERGIA especificamente
(a força de oscilador é outra história, ver abaixo).

**Tarefa 2 — caráter orbital do job_0023**: HSE06 real -- transição
dominante occ=80, virt=82 (HOMO=81), offset **(-1,+1)**, peso 78%,
f=0.063 (quase escuro). sTDA escolhido (maior f, contraste 3.6x) --
offset **(0,0)** (HOMO→LUMO canônico limpo, coeficiente 0.91),
f=0.626 (brilhante). **Padrão NOVO, diferente dos já vistos**:
job_0070 tinha o sTDA escolhendo o estado ERRADO (não-canônico) com um
competidor certo (canônico) por perto, contraste baixo, capturável
pelo critério de contraste. Aqui é o oposto: o sTDA está CONFIANTE
(contraste alto) no estado CANÔNICO certo-pela-própria-lógica-do-sTDA,
mas o HSE06 real prefere um estado não-canônico bem mais fraco e
distante em ranking dentro do próprio sTDA -- um modo de falha que o
contraste não pode, por construção, detectar (mede confiança interna
do sTDA, não correspondência com o HSE06 real).

**Tarefa 3 — varredura de risco nos 112 restantes** (`screen_confident_conformers_for_hidden_state_risk.py`,
critério: existe ALGUM estado não-canônico com f≥0.02 em qualquer
lugar da janela de 10 estados?): **112/112 (100%) flagados** --
critério sem nenhum poder discriminante. Estados não-canônicos fracos
(f~0.02-0.2) são ubíquos no espectro de baixa energia desse sistema
(plausível quimicamente: múltiplos grupos doadores parcialmente
conjugados geram várias excitações locais quase-degeneradas), não uma
anomalia rara exclusiva do job_0023. **Resultado inconclusivo**: com
apenas 1 caso confirmado desse tipo de falha, não há informação
suficiente pra calibrar um limiar que separe risco real de
contaminação eletrônica normal. Não dá pra pré-filtrar os 117 por essa
via.

**Tarefa 4 — lote de controle por amostragem aleatória** (não
submetido), `select_random_population_weighted_batch.py`: sorteio
sem reposição, ponderado por população (não pelo valor previsto),
dentro do bucket "confiável" (contraste>2.0x, ainda não confirmado,
N=112 elegíveis):

| job_dir | população | f previsto |
|---|---|---|
| job_0131 | 0.0053 | 0.807 |
| job_0120 | 0.0150 | 0.791 |
| job_0065 | 0.0257 | 0.788 |
| job_0106 | 0.0022 | 1.058 |
| job_0087 | 0.0013 | 0.955 |
| job_0138 | 0.0055 | 0.859 |
| job_0057 | 0.0033 | 1.055 |
| job_0081 | 0.0051 | 1.060 |

**Custo estimado**: mesmo nível (`escf`-only, ~7.8 núcleo-horas/confôrmero,
sem otimização geométrica) -- **8 × 7.8 ≈ 62 núcleo-horas total**,
mesma escala trivial de antes. Esse lote mede a taxa de erro real
"típica" do bucket confiável, sem o viés de seleção da Tarefa 1 --
complementa (não substitui) os lotes por extremo de valor, que
continuam relevantes pra limitar o pior caso nos confôrmeros de maior
peso. **Nada submetido.**

### 9.25. Lote de controle rodado, bug crítico no casamento de caráter encontrado e corrigido, recalibração de magnitude confirmada (2026-09-11)

**Lote de controle (8 aleatórios, ponderados por população) rodado**:
job_0057, 0065, 0081, 0087, 0106, 0120, 0131, 0138, todos `escf`+`ridft`
exit 0, sem falha de texto. Resultado bruto: **8/8 superestimam f**,
sem exceção, razão 1.62-2.63x -- direcional e limpo, não explicável
por regressão à média (esse lote não foi selecionado por valor
extremo).

**Bug crítico encontrado ao montar o pool limpo para a regressão**: o
algoritmo de casamento de caráter (`batch_stda_s1_character_matched.py`)
escolhia, entre estados sTDA com a mesma assinatura orbital do HSE06,
o de MAIOR COEFICIENTE isolado -- mas quando o caráter verdadeiro está
fragmentado entre 2-4 estados sTDA quase degenerados (comum neste
sistema), coeficientes quase empatados (ex. 0.55 vs 0.53) podem
corresponder a forças de oscilador completamente diferentes (ex.
f=0.06 vs f=0.31, 5x). **10/20 (50%) dos confôrmeros rotulados
originais do neutro tinham essa ambiguidade**, e em ~4 casos
(job_0011, job_0040, job_0041, job_0112) o algoritmo pegou
consistentemente a opção mais escura errada. Monoânion B e diânion
**não foram afetados** (0 mudanças) -- o resultado confirmado
Q²=0.412 (Δlog₁₀f do monoânion B) permanece válido. Monoânion A teve
mudanças (5/20) mas em valores já essencially zero (f~0.0001), sem
efeito prático.

**Correção**: escolher, entre os estados que batem a assinatura com
coeficiente ≥0.3 (não negligível), o de MAIOR FORÇA DE OSCILADOR --
não o de maior coeficiente. Alinha a escolha com o que B2 precisa
prever (brilho real), em vez de um desempate arbitrário que não
correlaciona com brilho. Consequência prática importante:
`job_0010`, `job_0038` e `job_0096` (do lote de extremos, antes
classificados como aparentes "mismatch tipo job_0023") **não eram
mismatches reais** -- a "verdade" usada na comparação estava errada
por causa do bug. Com a correção, a escolha por maior-f de produção
E o caráter correto COINCIDEM nesses 3 casos. **Só `job_0023`
continua sendo um mismatch de seleção de estado genuíno e confirmado**
-- taxa de mismatch real entre os 17 confôrmeros confirmados (9 do
lote de extremos + 8 do lote aleatório) = **1/17 = 5.9%** (IC 95%
Wilson [1.0%, 27.0%] -- amostra pequena, intervalo largo).

**Offset do neutro recalculado** (com os rótulos dos 20 originais
corrigidos):
```
delta_energy_eV(neutro) = 1.2113 ± 0.1723 eV   (era 1.3206 ± 0.2140)
delta_log10_f(neutro)   = 0.0218 ± 0.2906       (era 0.2066 ± 0.4305)
```
O viés médio de força de oscilador CAI de +0.207 pra quase zero
(+0.022) -- o valor antigo estava artificialmente inflado pelo bug.
Conclusão "sem sinal recuperável" (nenhuma feature geométrica passa de
Q²>0) se mantém idêntica nos dois alvos com os dados corrigidos --
não precisou redo a calibração por permutação.

**Testando dependência de magnitude** (`test_b2_osc_magnitude_recalibration.py`,
pool limpo N=28 = 20 originais corrigidos + 8 do lote aleatório,
ambos com sTDA casado por caráter, sem casos especiais nem
duplicatas):

```
log10(f_HSE06) = -0.203(±0.038) + 0.579(±0.059) x log10(f_sTDA)
```

Teste b=1: t=-7.14, dof=26, **p<0.0001** -- a inclinação difere
SIGNIFICATIVAMENTE de 1. Não é só um deslocamento de intercepto mal
calibrado -- é uma relação de COMPRESSÃO de magnitude real: o sTDA
exagera sistematicamente o quão mais brilhante um estado brilhante é
em relação a um mais escuro, comparado ao HSE06 real.

**Validação (protocolo padrão, sem pular etapa)**:
- Q²_LOO = **+0.5225** (pool de 28)
- Permutação (2000 embaralhamentos): percentil **100.0**, p<0.0005 --
  **CONFIRMADO**, sinal não explicável por acaso.
- **Checagem fora da amostra** (8 do lote de extremos, excluindo
  job_0023 -- mismatch de categoria diferente -- e job_0018 --
  duplicata): erro absoluto médio cai de **0.3555 (offset antigo) para
  0.1246 (novo modelo de inclinação)** -- redução de ~65%, mesmo
  nesses 8 pontos enriquecidos por viés de seleção de extremos
  (Seção 9.24, Tarefa 1). Modelo generaliza bem fora da amostra de
  ajuste.

**Os dois riscos, lado a lado, para o espectro final do neutro**:

| Risco | Fonte | Taxa/magnitude | Corrigível? |
|---|---|---|---|
| Erro de MAGNITUDE de f | offset aditivo constante não captura relação real (inclinação≠1) | erro médio absoluto ~0.36 sem correção, ~0.12 com o modelo de inclinação (redução 65%, validado fora da amostra) | **Sim** -- adotar `log10(f)=-0.203+0.579*log10(f_sTDA)` no lugar do offset plano |
| Erro de SELEÇÃO DE ESTADO | sTDA confiante (contraste alto) no estado canônico errado, HSE06 prefere um não-canônico (padrão job_0023) | 1/17 confirmado = 5.9% (IC 95% Wilson [1.0%, 27.0%]) | **Não** -- sem critério de detecção automática disponível (Seção 9.24, Tarefa 3: varredura de 112 confôrmeros deu 100% falso-positivo, sem poder discriminante) |

**Próximo passo natural (não decidido aqui)**: readotar o modelo de
inclinação para recalcular as previsões de força de oscilador dos 170
confôrmeros do neutro (Seção 9.21), o que deve reduzir bastante o
problema dos f>1.0 fisicamente implausíveis notado na Seção 9.22 --
mas a incerteza residual de ~5.9% de erro de seleção de estado
continua sem solução automática e deve ficar propagada explicitamente
no espectro final.

### 9.26. Reauditoria pós-bug: escopo do bug, 2 novos mismatches confirmados, reaplicação final (2026-09-11)

**Tarefa 1 — o bug NÃO afeta a classificação dos 170**: verificado por
inspeção do código -- `characterize_stda_selection_confidence.py`
(que gera contraste/E/f pra todo o ensemble) usa `max(window, key=f)`
puro, sem nenhum desempate por coeficiente entre estados candidatos; e
`analyze_stda_state_selection.py` (que gerou os casos de referência do
limiar de 1.26x) só olha o próprio estado de maior-f, também sem
desempate cruzado. **Nenhum dos dois usa a função com bug**
(`batch_stda_s1_character_matched.py`, usada só pra construir os
rótulos de treino do B2 a partir de HSE06 conhecido). Confirmado:
contraste/E/f dos 170 não mudam, não precisou reprocessar.

**Tarefa 2 — checagem rigorosa revela 2 mismatches novos, limiar de
contraste comprovadamente inútil pra esse risco**: em vez de confiar
na flag antiga e simplista de `analyze_stda_state_selection.py`
(que só testava "o próprio pick do sTDA é canônico?", sem comparar
contra o alvo real do HSE06), comparei diretamente produção
(maior-f-na-janela) contra a verdade corrigida (casamento por caráter
com o algoritmo já corrigido) nos 36 confôrmeros confirmados únicos.
Resultado: **3/36 mismatches genuínos** -- `job_0023` (já conhecido) +
**`job_0035` e `job_0036` (novos)** -- todos com a MESMA assinatura
(produção prevê brilhante canônico E≈2.74eV/f≈0.63, real é escuro
E≈2.32-2.33eV/f≈0.007-0.015) e todos com **contraste 2.49-3.57x**,
dentro do bucket "confiável" (>2.0x). Achado extra: os 2 confôrmeros
que definiam o limiar antigo (job_0011, job_0072) são FALSOS
POSITIVOS sob a checagem rigorosa -- a flag antiga os marcava porque
o pick do sTDA não era canônico, mas na verdade seus picks bateram
corretamente com o alvo real (que também não é canônico); a
comparação simplista não tinha como saber disso.

**Taxa recalculada**: 3/36 = **8.3%** (IC 95% Wilson [2.9%, 21.8%]),
substituindo a estimativa anterior 1/17=5.9%.

**Conclusão prática**: não existe limiar de contraste (binário ou
gradual) estatisticamente sustentável pra esse risco -- os 3 casos
confirmados sentam-se confortavelmente na faixa "confiável", igual à
varredura de 112 confôrmeros da Seção 9.24 (Tarefa 3, 100% flagados,
zero poder discriminante) já indicava. Tratamento adotado: risco de
mismatch de seleção de estado é uma taxa de fundo UNIFORME (8.3%,
IC [2.9%,21.8%]), aplicada a todo confôrmero previsto por surrogate
igualmente, independente do contraste -- não uma função do contraste.
O limiar de contraste 1.26x/2.0x continua válido pra seu propósito
ORIGINAL e diferente (quase-empates genuínos entre 2 estados
competindo, Seção 9.17/9.20), só não serve pra este risco.

**Tarefa 3 -- offset de energia recalculado**:
```
delta_energy_eV(neutro) = 1.2113 ± 0.1723 eV   (era 1.3206 ± 0.2140)
```
Mudança moderada (~8% na média, ~20% no desvio padrão) -- não
irrelevante, mas não muda a conclusão qualitativa (neutro sem modelo
recuperável em energia, tratamento offset-apenas mantido). Valor
atualizado adotado na reaplicação.

**Tarefa 4 -- reaplicação final aos 170** (`apply_b2_osc_surrogate.py`
v2: aceita múltiplas fontes de HSE06 real, usa o modelo de inclinação
log-log em vez do offset plano pra f, adiciona a nota de taxa de
mismatch de fundo):

| Categoria | N (v1, Seção 9.21) | N (v2, corrigido) |
|---|---|---|
| HSE06 real | 20 | **37** (+17 confirmados nos lotes) |
| Estado forçado (especial) | 2 | 2 |
| Surrogate confiável (contraste>2.0x) | 118 | 104 |
| Surrogate incerteza inflada (1.26-2.0x) | 22 | 19 |
| Excluído (≤1.26x) | 8 | 8 |

**f>1.0 fisicamente implausível, notado na Seção 9.22, desaparece
completamente** -- máximo agora 0.665, mínimo 0.015, zero casos acima
de 1.0 (eram vários antes). Confirma que o modelo de inclinação
resolve o problema de compressão de magnitude identificado na Seção
9.25.

**b2_osc_predictions_full_ensemble.csv (v1) está OBSOLETO** -- usar
`b2_osc_predictions_full_ensemble_v2.csv` a partir de agora pra
qualquer uso posterior (ranking espectral, seleção de lote, etc.).
Toda linha de surrogate carrega agora `state_selection_mismatch_rate`
e `state_selection_mismatch_ci` explícitos, pra propagação de
incerteza no espectro final.

### 9.27. Fechamento do B2: concentração de peso, propagação da taxa de mismatch, status final (2026-09-11)

**Tarefa 1 -- peso populacional dos 3 mismatches conhecidos**:

| job_dir | população | % do total | rank (de 170) | percentil |
|---|---|---|---|---|
| job_0023 | 0.0229 | 2.29% | 10º | top 5.9% |
| job_0036 | 0.0019 | 0.19% | 97º | top 57.1% |
| job_0035 | 0.0015 | 0.15% | 104º | top 61.2% |

Combinado: 2.63% da população total do neutro (vs. 1.76% esperado se
fossem 3/170 aleatórios -- levemente acima do proporcional, não
dramaticamente). **Não são um grupo concentrado em peso alto** --
`job_0035`/`job_0036` são medianos/baixos. Mas `job_0023` sozinho
**é** um caso de peso alto (10º maior população de 170) -- sua
correção importou de verdade pro espectro agregado.

Contexto de exposição: 50.4% da população total do neutro já é HSE06
real (37 confôrmeros confirmados), 40.9% ainda depende do surrogate
(exposta à taxa de fundo de 8.3%), 8.8% excluída.

**Tarefa 2 -- impacto contrafactual** (`compare_mismatch_counterfactual_impact.py`):
por confôrmero, o erro SERIA grande se não detectado (5.7x-31x na
força de oscilador individual), mas o impacto agregado no espectro
completo é pequeno -- **+2.84%** na intensidade integrada total, já
que mesmo o maior dos 3 (`job_0023`) é só 2.3% do peso total. Isso
dá uma cota superior REAL (não hipotética) pro tipo de distorção que
mismatches não detectados poderiam causar nos 40.9% de população
ainda dependente do surrogate -- embora não seja garantia de que
outro mismatch, se existir, teria peso/magnitude parecidos.

**Tratamento de incerteza adotado pro espectro final do neutro**:
- Espectro central: modelo de inclinação (v2), valores reais onde
  disponíveis (50.4% da população).
- Taxa de mismatch (8.3%, IC 95% Wilson [2.9%, 21.8%]) reportada como
  **limite de confiabilidade explícito**, não como barra de erro
  gaussiana -- é erro categórico (estado errado), não erro de
  magnitude, então não se soma em quadratura com a incerteza do
  offset/inclinação.
- Impacto máximo plausível já observado (não hipotético): +2.84% na
  intensidade agregada, calibrado pelos 3 casos conhecidos.

**Tarefa 3 -- status final consolidado do B2 (e B1, referência)**:

| Estado | N (ensemble) | B1 (IV, C=O) | B2 energia S1 | B2 força de osc. | Taxa de mismatch B2 | Limitação de escopo |
|---|---|---|---|---|---|---|
| Neutro | 170 | offset 28.222±0.676 cm⁻¹ | offset 1.2113±0.1723 eV | inclinação log-log: a=-0.203, b=0.579 (Q²=0.52, permutação confirmada) | 8.3% [2.9%,21.8%], uniforme, não-detectável por contraste | 37/170 (50.4% da população) já é HSE06 real; 8/170 (8.8% população) excluído por baixa confiança |
| Monoânion A | 143 | Bader(O carbonílico), Q²=0.277 (perm. p=0.017) | nenhum modelo sobrevive à calibração | nenhum modelo sobrevive à calibração | N/A -- critério de seleção de estado não se aplica (S1 real é escuro, ver Seção 9.17) | S1 escuro reduz a relevância prática da força de oscilador pro espectro visível |
| Monoânion B | 11 | distância O···H, Q²=0.594 | distância O···H, Q²=0.596 (perm. p=0.008, CONFIRMADO) | distância O···H, Q²=0.412 (perm. p=0.011, CONFIRMADO) | N/A -- 11/11 (100%) já é HSE06 real, nenhuma extrapolação necessária | faixa de treino da feature extremamente estreita (0.0008 Å) -- **não usar o modelo em confôrmeros novos** se a busca CREST for expandida |
| Diânion | 3 | fora de escopo (N insuficiente) | fora de escopo | fora de escopo | N/A | N=3 não permite nenhum protocolo estatístico deste projeto |

**B2 está fechado** para os 4 estados, nas condições documentadas
acima. Pendente (fora do escopo desta etapa): montagem final da
figura/tabela do espectro UV-Vis ponderado por população, e expansão
de aprendizado ativo pro restante do ensemble de 327 confôrmeros
(monoânion A e diânion ainda não passaram por nenhuma etapa de B2 além
da extração inicial de S1).

### 9.28. Decisão de parar o active learning de força de oscilador (2026-09-11)

**Decisão**: não perseguir mais rodadas de aprendizado ativo pra
reduzir o erro agregado de força de oscilador do neutro (+2.84% na
intensidade integrada, pior cenário contrafactual medido nos 3
mismatches confirmados, Seção 9.27).

**Justificativa**: benchmarks de força de oscilador em TD-DFT contra
referências de alto nível (CASPT2, CC2) reportam erro absoluto médio
(MAE) tipicamente entre 0.05 e 0.13 em f -- um erro relativo de
dezenas de por cento em transições comuns -- e mesmo comparações entre
métodos do estado da arte (DLPNO-STEOM-CCSD) e experimento tratam
desvio de intensidade relativa como detalhe secundário aceitável, não
como falha (ver "Benchmarking TD-DFT and Wave Function Methods for
Oscillator Strengths and Excited-State Dipole Moments", arXiv
2011.13233). O erro do surrogate deste projeto (+2.84% agregado) é uma
fração pequena do erro que já existe estruturalmente entre HSE06/TD-DFT
e o valor "verdadeiro" de referência wavefunction, e menor ainda frente
ao erro de HSE06 vs. experimento em posição de banda (0.1-0.3 eV, já
documentado nos benchmarks deste projeto, Seção 8).

**Declaração formal**: o erro residual de ~2.84% na intensidade
integrada é aceito como dentro do padrão da área para força de
oscilador; não justifica rodadas adicionais de active learning.

**Tabela final de status do B2** (repetida da Seção 9.27 pra
referência única):

| Estado | N (ensemble) | B1 (IV, C=O) | B2 energia S1 | B2 força de osc. | Taxa de mismatch B2 | Limitação de escopo |
|---|---|---|---|---|---|---|
| Neutro | 170 | offset 28.222±0.676 cm⁻¹ | offset 1.2113±0.1723 eV | inclinação log-log: a=-0.203, b=0.579 (Q²=0.52, permutação confirmada) | 8.3% [2.9%,21.8%], uniforme | 37/170 (50.4% população) já HSE06 real; 8/170 (8.8%) excluído |
| Monoânion A | 143 | Bader(O carbonílico), Q²=0.277 | nenhum modelo sobrevive | nenhum modelo sobrevive | N/A (S1 escuro) | força de osc. pouco relevante pro espectro visível |
| Monoânion B | 11 | distância O···H, Q²=0.594 | distância O···H, Q²=0.596 (CONFIRMADO) | distância O···H, Q²=0.412 (CONFIRMADO) | N/A (100% HSE06 real) | faixa de treino 0.0008 Å -- não extrapolar |
| Diânion | 3 | fora de escopo | fora de escopo | fora de escopo | N/A | N insuficiente |

### 9.29. Figura final do espectro UV-Vis ponderado e status para o manuscrito (2026-09-11)

**Bloqueios encontrados** (registrados aqui, não contornados por
suposição):
- `plano_especiacao_ph_alizarina.md` e `plot_speciation.py` não
  existem nesta máquina -- e não existe NENHUM dado de ΔG/pKa por
  espécie em lugar nenhum do projeto (a calibração de pKa nunca foi
  iniciada). **A figura de especiação vs. pH não pode ser construída**
  até essa etapa (separada, ainda não começada) ser feita.
- **Veusz não está instalado em nenhum dos 7 ambientes conda desta
  máquina** (checado em 2026-09-11). Os `.vsz` já existentes no
  projeto foram gerados numa máquina diferente (path embutido
  `/home/user/Documentos/.../Artepelin_C`). `plot_b2_final_spectrum_veusz.py`
  foi escrito seguindo o padrão de `plot_optical.py`, mas não pôde ser
  executado aqui -- rodar onde o Veusz estiver disponível.
  `plot_b2_final_spectrum_matplotlib.py` foi executado aqui como
  substituto imediato (`b2_final_spectrum.png`).

**Dataset consolidado** (`build_b2_consolidated_dataset.py`,
`b2_consolidated_dataset.csv`, 196 linhas): neutro (162 confôrmeros,
91.2% da população coberta), monoânion B (11, 100%), monoânion A (20,
26.6% -- só os rotulados, sem modelo validado pros outros 123),
diânion (3, 100%). Coluna `delta_g_eV` presente mas vazia em toda
linha -- pré-requisito de pKa não existe ainda.

**Figura final construída** (`build_b2_final_spectrum.py` +
`plot_b2_final_spectrum_matplotlib.py`, `b2_final_spectrum.png`):
envelope UV-Vis ponderado por população, neutro + monoânion B (escopo
explícito -- monoânion A/diânion não têm cobertura/modelo suficiente
pra entrar aqui). Neutro: pico em 3.97 eV; monoânion B: pico em 3.32
eV, mais intenso e mais estreito (população mais concentrada em poucos
confôrmeros de caráter parecido). Legenda inclui a taxa de mismatch do
neutro (8.3%, IC [2.9%,21.8%]) e a nota de faixa de treino do
monoânion B, como pedido.

**Tarefa 4 -- sumário de status pro manuscrito** (baseado na estrutura
já documentada neste arquivo, já que o `plano_especiacao_ph_alizarina.md`
referenciado não existe nesta máquina -- não posso confirmar contra
"Seções 6-9" literalmente):

**Já coberto por B1+B2:**
- Busca conformacional + filtro de população de Boltzmann (Camada A,
  Seções 5-6) -- completo pros 4 estados.
- Refino DFT barato (r2SCAN-3c) + reranking de população (Seção 7) --
  completo.
- Espectros de IV/UV-Vis de triagem, população-ponderados (PBE0/def2-SVP,
  Seção 8) -- completo, com a limitação documentada de modos
  imaginários sistemáticos (Seção 8.1-8.3).
- B1 (deslocamento do C=O por IV): fechado nos 4 estados (Seção 9.14).
- B2 (energia S1 + força de oscilador): fechado nos 4 estados, com
  correção de bug e taxa de mismatch residual documentada (Seções
  9.15-9.28).
- Espectro UV-Vis final ponderado por população (neutro + monoânion
  B): construído nesta seção.

**Ainda falta antes de escrever:**
- **Calibração de pKa / ΔG por espécie** -- prerequisito bloqueante
  pra qualquer figura de especiação vs. pH; protocolo de seleção de
  features já está congelado e pronto pra reuso
  (`feature_selection_protocol.md`), só não foi aplicado ainda.
- **Cobertura completa do monoânion A**: só 26.6% da população tem
  HSE06 real; sem modelo B2 validado pro resto -- decidir se vale a
  pena rotular mais confôrmeros (dado que S1 é escuro nesse estado,
  pode ser de baixa prioridade) ou aceitar cobertura parcial.
- **Seleção final (top-N) e reotimização de produção** nos parâmetros
  exatos do Paper 1 (Estágio 6) -- não iniciado.
- **Repetição pra ARS** (alizarina, diretório irmão já existe em
  `~/work_turbomole/alizarina/` com figuras próprias de IV/Bader/ordem
  de ligação, mas não passou pelo funil B1/B2 de conformação -- status
  exato não verificado nesta sessão).
- **Validação contra titulação experimental** (Camuri/Pazin, Paper 1)
  -- não iniciada, depende da calibração de pKa acima.
- **Descritores de reatividade de Fukui** -- mencionados no plano
  original, nenhum trabalho iniciado nesta sessão.
- Rodar `plot_b2_final_spectrum_veusz.py` numa máquina com Veusz, pra
  substituir a versão matplotlib pela figura de publicação.

### 9.30. Calibração de pKa -- levantamento (2026-09-11)

**Tarefa 1 -- cobertura do monoânion A**: dos 73.4% de população NÃO
cobertos, o maior confôrmero individual sozinho já é 5.08% (`job_0010`)
-- não é uma cauda desprezível. Dos 20 confôrmeros de maior população
geral, só 5 estão cobertos (incluindo os 3 maiores: `job_0093` 7.67%,
`job_0062` 7.21%, `job_0055` 5.53%). Ou seja: a semente de população
do desenho de aprendizado ativo pegou os líderes, mas a cauda logo
abaixo (vários confôrmeros de 2-5% cada) ficou de fora em peso
substancial. Não é "só cauda irrelevante".

Contexto útil pra decisão: o offset de energia do monoânion A (dos 20
cobertos) tem desvio padrão bem mais apertado que o do neutro (0.048
eV vs. 0.17 eV) -- um tratamento offset pros 73.4% restantes seria
proporcionalmente mais confiável ali do que foi pro neutro, mesmo sem
modelo por feature (Δlog₁₀f: 0.1256±0.3120, similar ao neutro).

**Tarefa 2 -- levantamento de pKa**:

*Literatura*: encontrado um valor citado de **pKa = 4.65** pra
Artepilin C em "pH-Dependence Cytotoxicity Evaluation of Artepillin C
against Tumor Cells" (Life, 2023), referenciando [17]/[18]. Rastreado
até a fonte primária (2026-09-14): **[17] Camuri, Costa, Ito, Pazin,
"Optical absorption and fluorescence spectroscopy studies of
Artepillin C..." Spectrochim. Acta A 2018, 198:71-77** e **[18]
Camuri, da Costa, Ito, Pazin, "pH and Charge Effects Behind the
Interaction of Artepillin C... " Photochem. Photobiol. 2019,
95:1345-1351** -- **o MESMO grupo (Camuri/Pazin) cujos espectros
experimentais já ancoram o Paper 1 inteiro**. Confirmado por busca
adicional: **pKa=4.65 é determinação experimental PRIMÁRIA da
CARBOXILA**, via titulação ácido-base acompanhada pela posição
espectral da banda de absorção UV próximo -- ou seja, medido pela
MESMA técnica (UV-Vis) que o B2 deste projeto simula, não uma citação
secundária de fonte não relacionada. Cobre o par neutro↔monoânion A
(o mais relevante -- "dominante em pH fisiológico", Seção 2).
**Nenhum pKa de fenol foi determinado nesses artigos** -- os pares
neutro↔monoânion B e monoânion↔diânion continuam sem nenhuma
referência experimental direta localizada.

*Vetor de features montado* (`build_pka_feature_vector.py`,
`pka_feature_vector.csv`, médias ponderadas por população dentro do
subconjunto rotulado de cada estado):

| Descritor | Disponível? | Nota |
|---|---|---|
| Gap HOMO-LUMO | **Sim** | extraído de graça do `ridft.out` já existente, sem cálculo novo |
| Carga Bader / ordem de ligação C=O carbonílico | Parcial | só neutro/monoânion A têm `population_feature.csv`; monoânion B/diânion nunca tiveram essa análise rodada |
| Carga/ordem de ligação no O-H IONIZÁVEL especificamente | **Não** | `population_feature.csv` rastreou o carbonilo (alvo do B1), não o O-H que de fato ioniza -- são átomos diferentes na maioria dos confôrmeros |
| ΔG via `freeh` | **Não** | `freeh` nunca foi rodado neste projeto |
| Momentos COSMO | **Não** | nenhum cálculo deste projeto usou solvatação implícita -- tudo em fase gasosa |

### 9.31. Extensão da cobertura B1 do monoânion A aos 143 confôrmeros (2026-09-11)

Em vez de offset/exclusão pra Tarefa 1 da Seção 9.30, o Roberto propôs
estender a ETAPA DE POPULAÇÃO (Bader) -- não a rotulagem HSE06 completa
-- aos 123 confôrmeros restantes, já que o modelo do B1
(`4.7223*bader_charge_carbonyl_o - 73.9249`, Q²=0.277, Seção 9.12) só
precisa da carga Bader, não do S1/TD-DFT do B2.

**Correção de premissa encontrada ao dimensionar**: `setup_population_job.py`
reaproveita o `mos` já convergido do diretório de ORIGEM -- que pros 20
rotulados originais era o `mos` em nível **HSE06** (não r2SCAN-3c). Os
123 restantes só tinham `mos` r2SCAN-3c -- precisou de um `ridft` novo
(SCF único, sem otimização geométrica, ~0.27 núcleo-horas/confôrmero)
antes da análise de população propriamente dita, não é totalmente
livre de HSE06 como a proposta original sugeria. Ainda assim, muito
mais barato que rotulagem completa (sem `escf`, sem `aoforce`).

**Execução** (aprovada, rodada): 123 `ridft`-only (job array 5133606,
123/123 sucesso) -> `setup_population_job.py` por confôrmero ->
população/Bader (job array 5134304, `dscf -proper` + `ridft -proper` +
`bader`, 123/123 sucesso, throttle 30 simultâneas) ->
`compute_b1_population_features.py` (123/123 extraído) ->
`apply_b1_model.py` com o modelo já calibrado e faixa de validade já
travada.

**Resultado**: **119/123 (96.7%) dentro da faixa de treino**
(interpolação, mesma confiança do Q²=0.277 medido); só **4/123** fora
da faixa, flagados `EXTRAPOLATION`, e **irrelevantes em peso**
(população combinada dos 4: ~1.4% -- job_0015 0.41%, job_0038 0.26%,
job_0049 0.60%, job_0100 0.12%). **Cobertura do B1 do monoânion A
agora é 143/143 (100% da população)**, era 26.6%.

**Escopo**: isso fecha a cobertura do **B1** (deslocamento IV do C=O)
do monoânion A, não do B2 -- o B2 desse estado continua sem modelo por
feature (Seção 9.18, S1 escuro), tratamento inalterado.

### 9.32. Calibração de pKa -- levantamento (continuação da Seção 9.30)

**Conclusão desta etapa (antes da Seção 9.33 abaixo) -- não calculei
pKa aqui ainda.** Cobertura experimental é PARCIAL (1 valor citado,
provavelmente só carboxila, fonte secundária não verificada) e os
descritores disponíveis são um subconjunto pequeno dos 5 originalmente
listados (só gap HOMO-LUMO cobre os 3 pares; carga/ordem de ligação
cobre só 1 dos 3). Estratégia decidida em seguida: usar o pKa de 4.65
como âncora só pro par carboxila (Seção 9.33), deixando os outros 2
pares (fenol, segunda desprotonação) sem calibração análoga por
enquanto.

### 9.33. Ciclo termodinâmico de pKa completo -- par carboxila, fechado (2026-09-15)

Metodologia rigorosa aprovada e executada por completo (Tarefas 1-4 da
solicitação): `freeh` (correção térmica, sem DFT novo, reaproveitando
Hessianos do B1) + COSMO(água, ε=78.4, single-point sobre a geometria
já otimizada, sem reotimizar) + ciclo termodinâmico completo com
ΔG_solv(H⁺)=-265.9 kcal/mol (Tissandier et al. 1998, *J. Phys. Chem.
A* 102, 7787 -- aproximação "cluster-pair", valor de consenso da
literatura) + G(H⁺,gás)=-6.28 kcal/mol (Sackur-Tetrode, valor padrão
usado na literatura de cálculo de pKa).

**Pegadinhas de infraestrutura encontradas e documentadas**:
- `freeh` é interativo (como `define`), com 3 prompts nesta ordem
  exata: número de simetria (RETURN aceita), fator de escala de
  frequência (RETURN aceita o padrão BP/SVP ~0.9914 -- não usado antes
  neste projeto, sinalizado como simplificação conhecida, não escolha
  deliberada), faixa de T/P (RETURN aceita 298.15K/0.1MPa, "q" na
  segunda vez que o prompt reaparece pra sair limpo). Rodar SEMPRE numa
  cópia, nunca no diretório de produção -- confirmado que roda no
  local mutando o `control` (`$last step freeh`, inofensivo mas
  desnecessário).
- `cosmoprep` tem uma sequência de ~14 prompts (documentada
  integralmente em `setup_cosmo_job.py`); bug encontrado: aceitar o
  padrão do prompt `amat` com "*" escreve `amat file=*` no `$cosmo`,
  que o `ridft` rejeita ("unknown keyword") -- removido por
  pós-processamento.
- Falha esporádica (1/40, `job_0132`) em que `cosmoprep` reportou
  "ended normally" mas NÃO escreveu o grupo `$cosmo` -- reproduzida
  2x pelo script de produção, mas funcionou numa execução manual mais
  cuidadosa (timeout maior). Provável problema de timing/carga no
  cluster durante execução concorrente, não uma incompatibilidade real
  com aquele confôrmero -- corrigido manualmente.
- Composto de referência (ácido acético, Tarefa 4) construído do zero
  (sem RDKit/OpenBabel disponíveis -- geometria inicial escrita à mão,
  refinada por otimização completa HSE06/def2-TZVP com os thresholds
  apertados de produção do Paper 1). Ácido acético neutro + acetato
  (base conjugada, geometria construída removendo o próton ácido da
  estrutura já otimizada) -- ambos otimizados, sem modos imaginários,
  `freeh` + COSMO completos. Custo total do par de referência: <0.5
  núcleo-horas (sistema de 7-8 átomos).

**Resultado do ciclo (par carboxila neutro↔monoânion A, N=20
rotulados de cada lado, ponderado por população dentro do
subconjunto)**:

| Sistema | pKa previsto | pKa experimental | Resíduo (pKa) | Resíduo (kcal/mol) |
|---|---|---|---|---|
| Ácido acético (referência) | 8.146 | 4.76 | -3.386 | -4.619 |
| **Artepilin C (carboxila)** | **8.689** | **4.65** | **-4.039** | **-5.510** |

**Os dois resíduos são próximos** (diferença de 0.65 unidades de pKa,
0.89 kcal/mol) -- evidência de erro sistemático TRANSFERÍVEL entre os
dois sistemas neste nível de teoria (HSE06-D3(BJ)/def2-TZVP + COSMO
single-point sem reotimização em solvente, sem correções além da OC).
Aplicando a correção calibrada no ácido acético (+3.386) à previsão
bruta da Artepilin C (8.689) → pKa corrigido ≈ **5.30**, a ~0.65
unidades do experimental (4.65) -- dentro da precisão típica esperada
de uma calibração por ponto único.

**Interpretação**: o erro residual (~3.4-4.0 unidades, direção
consistente: superestima o pKa em ambos) é compatível com limitações
já conhecidas e documentadas na literatura para ciclos de pKa via
COSMO single-point sem parametrização de cavidade específica pra
ânions nem correções explícitas de moléculas de água -- não é evidência
de erro no pipeline deste projeto. A concordância entre os dois
resíduos é o resultado mais importante desta etapa: sustenta usar o
offset calibrado no ácido acético como correção sistemática pro par
carboxila da Artepilin C, mas **não deve ser extrapolado sem novo teste
pros pares fenol (neutro↔monoânion B) ou segunda desprotonação
(monoânion↔diânion)**, que envolvem uma classe química diferente
(fenol, não carboxila) sem nenhuma calibração análoga feita ainda.

### 9.25. Lote de controle rodado — taxa de erro consolidada (17 pontos)

Lote de 8 rodado (`escf`-only, lista de jobs no filesystem
compartilhado desde o início desta vez). 8/8 exit 0, sem falha de
texto. Caráter orbital 8/8: todos batem HOMO->LUMO real com HOMO->LUMO
do sTDA -- nenhum mismatch tipo job_0023 nesta amostra. Magnitude de
f, 8/8: razão f_previsto/f_real entre 1.62x-2.63x em todos, sempre
superestimando -- viés sistemático, não erro aleatório.

Taxa consolidada (17 = 9 do lote anterior sem a duplicata job_0018 +
8 deste): mismatch de caráter tipo job_0023 = 1/17 = 5.9% (IC 95%
Wilson [1.0%, 27.0%]); razão f_previsto/f_real > 1.5x em 16/17 = 94%;
mediana da razão (excluindo job_0023) = 1.90x, média 2.01x
(min 1.22x, max 3.28x).

O achado mais importante não é o mismatch raro -- é o viés
sistemático. A taxa de mismatch (5.9%) fica abaixo de 10-15% na
estimativa pontual, mas o IC (limite superior 27%) não descarta uma
taxa real maior com essa amostra pequena. Mais relevante: 94% dos 17
confôrmeros confirmados mostram intensidade superestimada por ~2x
mesmo com seleção de estado perfeitamente correta -- não é problema
de seleção, é o offset em log-espaço (Δlog10f=+0.207, N=20) não
generalizando bem, consistente com os f previstos >1.0 já notados na
Seção 9.22. Decisão a preparar com o Roberto: tratar a intensidade do
neutro inteiro (não só os já flagados por contraste) com banda de
incerteza larga e viés conhecido, ou recalibrar o offset dependente
da magnitude de f do sTDA usando os 17 pontos novos.

### 9.34. Par fenol + especiação vs pH -- em andamento (2026-09-23)

Plano do usuário: (a) calibrar o par FENOL com um composto fenólico de
referência e (b) construir a figura de especiação vs pH em paralelo;
(c) Estágio 6 (top-N + reotimização, 800 passos se houver otimização)
adiado por depender de (a).

**Composto de referência fenólico**: fenol (pKa aquoso 9.99 a 25 °C;
também 9.90±0.06 a 23 °C em outra medição -- valor padrão da
literatura, ver Gross/Kass *JACS* 2002, 124, 2264 e PMC10376140) e
fenolato, HSE06-D3(BJ)/def2-TZVP, otimização completa (`jobex -c 800`,
thresholds apertados do Paper 1), Hessiano, `freeh`, COSMO(água).
Geometrias iniciais escritas à mão (hexágono planar; `jobex` refina).
Ácido acético/acetato agora persistidos em
`artepillin_C_studies/pka_reference/acetic_acid/` (antes só em `/tmp`).
**Ressalva de classe**: fenol simples NÃO tem os dois prenilas orto nem
o substituinte vinil-carboxilato para da artepilina C -- o offset
calibrado nele só é validado como transferível se um segundo composto
da mesma classe concordar (candidato: 2,6-dimetilfenol, pKa 10.6,
imita o impedimento orto; não feito).

**Regra do login node**: os testes COSMO (~8 min de parede, 4 threads) e o
pipeline do ácido acético (Seção 9.33) rodaram no login node -- violação
da regra do usuário (reincidência). Daqui pra frente todo cálculo vai por
`sbatch`: `job_reference_compound_pipeline.sh` (fenol/fenolato:
jobex→aoforce→freeh→COSMO→ridft), `job_cosmo_prep.sh` (preparo COSMO em
array, 3 tentativas -- ver bug esporádico do `cosmoprep`),
`job_cosmo_ridft.sh` (SCF COSMO dos confôrmeros). Jobs: **5205519**
(fenol/fenolato), **5205520** (preparo COSMO de 11 mon.B + 3 diânions,
14/14 ok na 1ª tentativa), **5205521** (SCF COSMO desses 14,
dependência afterok). Custo estimado ≈ 8–10 núcleo-horas no total
(pipeline fenol ~0.5 + COSMO ~7–9).

**Verificações feitas**: (i) 54/54 conformeros rotulados têm Hessiano
completo e **zero modos imaginários** (< −10 cm⁻¹) -- a preocupação de
que a geometria r2SCAN-3c não fosse mínimo HSE06 não se aplica ao
`freeh`; (ii) **entropia conformacional** (G_ens = ⟨G⟩ − T·S_conf,
S_conf = −R Σ p ln(p/g), com a degenerescência g dos rotâmeros do CREST
e as populações r2SCAN-3c do ensemble COMPLETO): N_eff ≈ 1052 (neutro),
1067 (mon.A), 630 (mon.B), 398 (diânion); efeito sobre o ΔG_aq de cada
par: N→A −0.008 kcal/mol (**−0.006 pKa**), N→B +0.303 kcal/mol
(**+0.22 pKa**), A→D +0.584 kcal/mol (**+0.43 pKa**). Pequeno frente à
incerteza de calibração (~0.65), NÃO usado no resultado primário
(mantida a média ponderada por população, consistente com o carboxila);
reportado como sensibilidade (`calibrate_pka_pair.py
--target-sconf-shift-kcal`).

**Scripts novos** (todos testados): `setup_hse06_opt_template.py`
(HSE06 + $statpt do Paper 1; bug evitado: filtro que apagaria a linha
`$closed shells`), `calibrate_pka_pair.py` (referência→resíduo→pKa
calibrado do alvo; regressão reproduz 8.689→5.302 do carboxila),
`compute_pka_thermodynamic_cycle.py` (fórmula refatorada em
`compute_cycle`, regressão reproduz 8.146/−3.386 do ác. acético),
`plot_speciation_vs_ph.py`.

**Especiação vs pH -- versão PARCIAL pronta** (`speciation_partial_NA.png`
/ `.dat`, script `plot_speciation_vs_ph.py`): só N⇌A, com pKa(N→A)
calibrado = 5.30 e marcador do experimental (4.65); banda ±0.65 =
espalhamento do offset (NÃO é IC). A região pH > 8 está sombreada como
NÃO VÁLIDA (fenol ainda não calibrado; o 8.0 é heurística, não valor
calculado). O modelo de 4 microestados
(w_A = pH−pKa_NA, w_B = pH−pKa_NB, w_D = 2pH−pKa_NA−pKa_AD; média das
duas rotas + erro de fechamento do ciclo quando ambas existem) está
implementado e será rodado quando os pKa fenólicos existirem.

**Resultado do par fenol (2026-09-23, jobs 5205519/20/21 concluídos, 0
falhas; fenol 19 ciclos, fenolato 14, 0 modos imaginários; 14/14 COSMO
sem falha de texto):**

| Referência | pKa bruto | pKa exp. | Resíduo (pKa) | Resíduo (kcal/mol) |
|---|---|---|---|---|
| ácido acético (0→−1) | 8.146 | 4.76 | −3.386 | −4.619 |
| **fenol (0→−1)** | 13.762 | 9.99 | **−3.772** | **−5.146** |
| artepilina, carboxila (0→−1) | 8.689 | 4.65 | −4.039 | −5.510 |

Os três resíduos ficam entre −3.4 e −4.0 pKa **apesar de serem duas
classes químicas** (ácido carboxílico vs. fenol) -- sugere que o erro
sistemático é dominado pela solvatação do ânion/termo do próton neste
esquema COSMO, não pelo grupo funcional.

**pKa microscópicos da artepilina C, calibrados** (offset carboxila =
ácido acético −3.386; offset fenol = fenol −3.772; entre parênteses =
sensibilidade com o termo −T·S_conf da Seção 9.34 acima):

| Par | pKa bruto | pKa calibrado | Confiança |
|---|---|---|---|
| N→A (carboxila) | 8.689 | **5.302** (5.297) | ancorado (exp. 4.65; +0.65) |
| N→B (fenol) | 12.040 | **8.268** (8.490) | razoável: mesma mudança de carga 0→−1 da referência |
| B→D (carboxila do fenolato) | 12.395 | 9.009 (9.208) | **BAIXA**: −1→−2, offset vem de 0→−1 |
| A→D (fenol do carboxilato) | 15.747 | 11.975 (12.403) | **BAIXA**: −1→−2, offset vem de 0→−1 |

Fechamento do ciclo (rota via A vs. via B) = 0.000 -- é identidade
por construção (os G são funções de estado, e cada offset é aplicado a
um par de cada classe), NÃO é validação independente.

**Tautômero B (fenolato)**: [B]/[A] = 10^(pKa_NA − pKa_NB) = **1.1×10⁻³**
(~0.1% do monoânion, independente do pH) -- confirma a expectativa
original da Seção 2 ("população desprezível") e diz que o UV-Vis
do monoânion B (Seção 9.27) só pesa pouco no espectro do monoânion
real.

**Alerta sobre o diânion**: a diferença pKa(A→D) − pKa(N→B) = 3.7
unidades é implausivelmente grande como interação eletrostática entre um
carboxilato e um fenol/fenolato a ~7 Å (estimativa qualitativa: em
diácidos pequenos costuma ficar abaixo de ~1 unidade; não verificado com
dado de literatura nesta sessão). Leitura mais provável: o COSMO
single-point subestima a solvatação do diânion (resíduo cresce com a
carga), então o offset de 0→−1 NÃO transfere pro passo −1→−2 e os
pKa envolvendo o diânion estão superestimados. **Teste barato e
direto**: uma referência com passo −1→−2 conhecido (ex. 4-hidroxibenzoato:
carboxilato + fenol, pKa₂ do fenol ≈ 9.3 -- valor a verificar na
literatura antes de usar; ~16 átomos, <1 núcleo-h) -- NÃO feito,
depende de decisão do usuário.

**Especiação vs pH -- versão completa gerada** (`speciation_full.png`
/ `.dat`, `plot_speciation_vs_ph.py`): N domina abaixo de pH ~5, **A
domina de pH ~6 a ~11 (0.99 em pH 7.4)**, D só passa a dominar acima de
pH ~12 (ponto de cruzamento A/D ≈ 11.98) -- mas essa região depende do
pKa(A→D) de baixa confiança. Banda ±0.65 = espalhamento do offset (NÃO
é IC). Legenda da figura já carrega o aviso do diânion.

### 9.35. Mais precisão: 2ª referência fenólica + referência −1→−2 (2026-09-24, em andamento)

Pedido do usuário: calcular a referência do passo −1→−2 e um segundo fenol
para melhorar a precisão. **Referências novas** (HSE06-D3(BJ)/def2-TZVP, mesma
receita: `jobex -c 800` → aoforce → `freeh` → COSMO(água)):
- **4-hidroxibenzoato** monoânion→dianião (fenol do carboxilato, −1→−2;
  ≈ o passo A→D da artepilina): pKa₂ exp. **9.40**, mas a literatura
  dispersa entre **9.3 e 9.49** (fontes diferentes; 4.54–4.58 p/ a carboxila)
  -- uso o ponto médio e propago o intervalo.
- **2,6-dimetilfenol** (imita os prenilas orto): pKa **10.63** a 25 °C
  (tabela de fenóis alquilados em *Phenolic Resins: A Century of Progress*,
  Tab. 2.1, achada por busca -- **fonte secundária ÚNICA, não conferida em
  fonte primária**; a mesma tabela dá 10.00 p/ fenol vs 9.99 acima).
  Valor "10.6" que eu tinha em mente antes era de memória e a busca só
  confirmou o do isômero 2,4 -- não usado.

**Aviso de infraestrutura**: o repo `scripts/` foi reorganizado por outra
sessão em 2026-09-24 (scripts rastreados → `scripts/scripts/`, ver cofre
`[[scripts-git-remote]]`); `setup_hse06_tzvp_template.py`/`xyz_to_coord.py`
mudaram de lugar -- `setup_hse06_opt_template.py` agora procura nos dois
lugares. Meus scripts do Paper 2 (não rastreados) seguem na raiz.

**Resultado parcial** (jobs 5209330: 4-HBA mono/di e 2,6-DMP fenolato ok;
2,6-DMP neutro estourou 3 h em 85 ciclos -> retomado, job 5209821):

| Referência | classe | pKa bruto | pKa exp. | Resíduo (pKa) |
|---|---|---|---|---|
| ácido acético | 0→−1 | 8.146 | 4.76 | −3.386 |
| fenol | 0→−1 | 13.762 | 9.99 | −3.772 |
| 4-hidroxibenzoato (fenol) | −1→−2 | 13.902 | 9.40 (9.30–9.49) | **−4.502** (−4.60…−4.41) |

Offset por classe (`summarize_pka_calibration.py`, média ± desvio amostral):
**0→−1: −3.579 ± 0.273 (N=2)**; **−1→−2: −4.502 (N=1, sem dispersão
estimável)**. O resíduo **cresce ~0.9 pKa com a carga** -- confirma a
suspeita da Seção 9.34 de que o offset de 0→−1 NÃO transfere pro
diânion (o COSMO single-point subestima mais a solvatação de espécies
mais carregadas). Com os offsets por classe: N→A **5.110** (exp. 4.65,
+0.46), N→B **8.461**, A→D **11.245** (era 11.975), B→D **7.894**
(era 9.009); incerteza da classe 0→−1 ±0.27, da classe −1→−2 ainda
desconhecida (N=1).

**Ressalva que continua aberta (e que liga ao Estágio 6)**: mesmo
corrigida, a lacuna pKa(A→D) − pKa(N→B) = 2.8 unidades é maior do que se
espera de uma interação eletrostática carboxilato…fenol (estimativa
qualitativa <~1; não verificada). Uma causa plausível NÃO capturada pelos
offsets: as referências foram **otimizadas em HSE06**, mas os confôrmeros
da artepilina C estão **congelados na geometria r2SCAN-3c** -- o erro de
não-relaxação (energia de relaxação r2SCAN-3c→HSE06, possivelmente
diferente entre estados de carga) entra só no lado da artepilina. Isso
significa que a reotimização HSE06 do Estágio 6 (adiado) teria efeito
direto na precisão do pKa; **magnitude não medida**.

**Resultado final da Seção 9.35 (2026-09-25; job 5209821 concluído, 2,6-DMP
neutro convergiu com +15 ciclos, ~100 no total; sem modos imaginários; modo
mais baixo 159 cm⁻¹; freeh e COSMO ok):**

| Referência | classe | pKa bruto | pKa exp. | Resíduo (pKa) | (kcal/mol) |
|---|---|---|---|---|---|
| ácido acético | 0→−1 | 8.146 | 4.76 | −3.386 | −4.619 |
| fenol | 0→−1 | 13.762 | 9.99 | −3.772 | −5.146 |
| 2,6-dimetilfenol | 0→−1 | 14.808 | 10.63* | −4.178 | −5.699 |
| 4-hidroxibenzoato (fenol) | −1→−2 | 13.902 | 9.40 (9.30–9.49) | −4.502 | −6.141 |

\* fonte secundária única, não conferida em fonte primária.

Offset por classe: **0→−1: −3.779 ± 0.396 (N=3; erro padrão da média
0.23)**; **−1→−2: −4.502 (N=1, dispersão não estimável; intervalo do pKa
experimental do 4-HBA propaga −4.60…−4.41)**. Com 3 referências a
dispersão da classe SUBIU (0.27→0.40): o resíduo parece crescer com o
tamanho/teor de alquilas (8 átomos −3.39; 13 átomos −3.77; 19 átomos −4.18;
artepilina, carboxila, 46 átomos, alvo fora da média: −4.04) -- hipótese
com N=3, não testada.

**pKa microscópicos calibrados da artepilina C** (entre () = com o termo
−T·S_conf da Seção 9.34, `--apply-sconf`):

| Par | bruto | calibrado | incerteza |
|---|---|---|---|
| N→A (carboxila) | 8.689 | **4.910** (4.904) | ±0.40 (classe); **exp. 4.65 → +0.26, dentro de 1σ, checagem FORA da amostra** (o N→A não entra no offset) |
| N→B (fenol) | 12.040 | **8.262** (8.484) | ±0.40 (classe) |
| A→D | 15.747 | **11.245** (11.673) | dispersão da classe −1→−2 NÃO medida (N=1) |
| B→D | 12.395 | **7.894** (8.093) | idem |

Fechamento do ciclo = 0.001 (identidade por construção). [B]/[A] =
10^(4.910−8.262) = **4.5×10⁻⁴**. Especiação final (`speciation_v3_4refs.png`
/ `.dat`, banda ±0.4 = dispersão da classe 0→−1, **assumida** também p/
−1→−2, NÃO é IC): N domina até pH ~4.9; **A domina de pH ~6 a ~10 (0.996
em pH 7.4)**; D cruza A em pH ≈ 11.2 e domina acima de ~12. Incerteza
prática: ±0.4 pKa nos pares 0→−1; desconhecida (provavelmente ≥ isso) nos
pares que envolvem o diânion.

**Limites que permanecem**: (i) classe −1→−2 com N=1; (ii) pKa exp. do
2,6-DMP de fonte única e do 4-HBA com dispersão de 0.19; (iii) lacuna
pK(A→D) − pK(N→B) = 3.0 unidades ainda parece grande (estimativa
qualitativa, não verificada) -- hipótese da não-relaxação r2SCAN-3c→HSE06
só do lado da artepilina, NÃO medida; (iv) mesmo nível de teoria (COSMO
single-point) sem moléculas de água explícitas. **Alavanca proposta (não
executada, requer aprovação de custo)**: estimar a energia de relaxação
r2SCAN-3c→HSE06 por passo de Newton, ΔE ≈ −½ gᵀH⁻¹g, usando o Hessiano
JÁ existente (aoforce do B1) + um único `rdgrad` por confôrmero no
nível HSE06 (custo estimado ~1–3 núcleo-h por confôrmero; o `gradient`
não está nos jobs rotulados) -- mediria diretamente se o erro difere
entre estados de carga, sem pagar reotimizações completas (~40–100
núcleo-h por confôrmero, estimativa).

### 9.36. Teste da energia de relaxação r2SCAN-3c→HSE06 -- hipótese REFUTADA (2026-09-25)

Hipótese (Seção 9.35): as referências de pKa foram otimizadas em HSE06,
mas os confôrmeros da artepilina C estão congelados na geometria r2SCAN-3c;
se a energia de relaxação diferisse entre estados de carga, entraria como
erro só do lado da artepilina, fora do alcance dos offsets de referência.

**Método** (`estimate_relaxation_energy.py`; `job_rdgrad_relax_test.sh`):
1 `rdgrad` HSE06-D3(BJ)/def2-TZVP por confôrmero, EM CÓPIA, na geometria
congelada; energia de relaxação harmônica por passo de Newton usando o
Hessiano projetado JÁ existente (aoforce do B1): dE = −½ Σ (vₖ·g)²/λₖ
sobre os modos não-TR (λ > corte). Amostra: os 5 confôrmeros de maior peso
por estado dentro do subconjunto rotulado (neutro, mon.A, mon.B; 77–94% da
população rotulada) + os 3 do diânion (100%) = 18. Job 5217058, 18/18 ok,
~22 min de CPU cada (**~6.6 núcleo-h no total**, abaixo dos ~20 estimados).
Gradiente residual (rms ~0.0024–0.0027 Ha/bohr, max 0.010–0.018) bem
acima do critério do Paper 1 (5e-4) -- i.e., de fato longe do mínimo HSE06.

**Resultado**: dE médio ponderado por estado (corte 1e-3 Ha/bohr²):
neutro −0.946 · mon.A −0.930 · mon.B −1.024 · diânion −1.140 kcal/mol --
**quase igual entre estados** (faixa por confôrmero −0.80…−1.18).
Deslocamento do pKa bruto (pKa_relaxado = pKa_congelado + deslocamento):

| Par | δ(ΔG) kcal/mol | δpKa (corte 1e-3) | faixa nos cortes 1e-4…5e-3 |
|---|---|---|---|
| N→A | +0.016 | +0.012 | +0.012…+0.034 |
| N→B | −0.078 | −0.057 | −0.028…−0.057 |
| A→D | −0.211 | −0.154 | −0.094…−0.154 |
| B→D | −0.116 | −0.085 | −0.033…−0.085 |

**Conclusão: a não-relaxação NÃO explica a lacuna pK(A→D) − pK(N→B)**
(o efeito é ≤0.15 pKa, dentro do ruído de ±0.4 da calibração) e **não
justifica o Estágio 6 (reotimização HSE06 completa) por motivo de
precisão do pKa**. Resultados primários NÃO alterados. Limites: (i)
estimativa harmônica -- passos de Newton de até 0.25–0.55 bohr nos modos
suaves (anarmônicos); o resultado só varia 0.03–0.15 pKa entre cortes;
(ii) só gradiente em fase gasosa -- relaxação na presença do campo de
reação COSMO (que poderia pesar mais nas espécies carregadas) NÃO testada;
(iii) amostra = confôrmeros de maior peso, não o ensemble inteiro.

**Suspeitos que restam** para a lacuna e para a dispersão da classe: (a)
erro do COSMO single-point que cresce com a carga (classe −1→−2 com N=1,
sem dispersão); (b) ausência de águas explícitas; (c) diânion com apenas 3
confôrmeros; (d) a própria lacuna pode ser parcialmente real -- a
expectativa de <~1 unidade foi uma estimativa qualitativa NÃO verificada
em literatura.

### 9.37. Impacto na precisão das lacunas de modelo dos monoânions + UV-Vis vs pH (2026-09-25)

Pedido do usuário: quantificar, por tipo, o impacto (i) da faixa de treino
estreita do monoânion B e (ii) da ausência de modelo para o monoânion A.
Script: `estimate_model_gap_impact.py` (saída `artepillin_C_studies/model_gap_impact.dat`).

**Fato que muda a pergunta:** o monoânion B tem os 11/11 confôrmeros em DFT
HSE06 real (aoforce + escf). O modelo de faixa estreita **nunca é aplicado a
nenhum confôrmero B** -> erro herdado do modelo = 0. O que resta é só
"quão bem o modelo generalizaria", o que é irrelevante para os resultados.
Monoânion A: 20/143 confôrmeros (26.6% da população) em DFT real; o modelo B1
foi estendido aos 143.

| Tipo | Monoânion B (11/11 real) | Monoânion A (20/143 real) |
|---|---|---|
| IV ν(C=O) | média 1739.25, sd 0.77, span 2.0 cm⁻¹; sem extrapolação | Δν rótulos −80.49 (sd 0.54); LOO offset-só 0.95 vs modelo Bader 0.77 (modelo compra 0.16 cm⁻¹); banda ponderada completa −80.12 (a extensão mudou a banda em +0.37 cm⁻¹ vs FWHM 10) -> **baixo** |
| UV-Vis | S1 3.317 eV (sd 0.007, span 0.031); f 0.804 (sd 0.045); LOO RMSE 0.011 eV; feature O···H (span 0.0008 Å) só divide o conjunto em 2 grupos (job_0008/9/10, f 0.727 vs 0.819) -> contraste de grupos, sem relação contínua, não interpretável; [B]/[A]=4.5e-4 -> erro diluído ~3e-6 eV | S1 2.574 eV é escuro (f≈2e-4); "primeiro estado brilhante" depende do limiar (3.70/3.91/4.38 eV para f≥0.02/0.05/0.1) -> **não usar**. Observáveis robustos: transição mais forte <6 eV 4.514 eV (sd 0.123; incerteza de cobertura 0.014 aleat./0.091 pior caso); centroide f-pond. 3–5.5 eV 4.592 eV (sd 0.062; 0.007/0.045) -> **baixo/limitado** |
| pKa/termoquímica (amostragem de G_aq) | totalmente coberto: 0 | neutro 20/170, monoA 20/143: sd de G_aq ≈0.53–0.55 kcal/mol; incerteza aleatória 0.05–0.06, pior caso correlacionado 0.39–0.40 kcal/mol. Par N→A: 0.061 pKa aleat., 0.580 pior caso; teste de representatividade r2SCAN-3c −0.062. N→B 0.040/0.284/+0.035; A→D 0.046/0.296/+0.097; B→D 0 |

**Conclusão:** impacto BAIXO para monoB (DFT real, peso desprezível) e
BAIXO/limitado para monoA (IV robusto; UV só via observáveis robustos; pKa
0.06 aleatório, ≤0.58 pior caso, comparável à dispersão de classe ±0.4).
Pela regra do usuário (impacto baixo -> executar a validação), segue a
Tarefa 2 (validação vs titulação Camuri/Pazin). Não é preciso rodar mais DFT
para o monoB.

**Bug encontrado e corrigido (Camada A UV-Vis):** `build_population_weighted_spectrum.py`
(`scripts/scripts/`, versionado; 2 linhas editadas) lia o grupo 2 do regex
(Hartree) como eV e o grupo 6 (f de velocidade) em vez de f de comprimento.
Corrigido para grupos 3 (eV) e 7 (f_len); os 4 espectros de
`layer_a/uvvis_population_weighted.dat` foram regenerados; os antigos foram
renomeados `.BUGGY_hartree_axis`. Limite: o escf PBE0/def2-SVP da Camada A
calculou só 10 estados (`$soes a 10`) -> para monoA (banda brilhante alta)
o espectro da Camada A é **não confiável**; usar os escf HSE06 (60 estados)
da Camada B. Máximos corrigidos Camada A: neutro 4.183 eV/296 nm, monoB
3.465 eV/358 nm, diânion 3.951 eV/314 nm, monoA quase escuro na janela.

**Espectro previsto vs pH** (`predict_uvvis_vs_ph.py`; sticks HSE06 reais,
σ=0.4 eV, pesos r2SCAN-3c renormalizados no conjunto com escf real;
especiação `speciation_v3_4refs.dat`; saída `uvvis_vs_ph.png/.dat`):

| Espécie | conf. reais | cobertura pop. | máx. (nm) | centroide 3–5.5 eV (nm) |
|---|---|---|---|---|
| N | 37/170 | 0.504 | 307 | 283 |
| A | 20/143 | 0.266 | 274 | 270 |
| B | 11/11 | 1.0 | 370 | 312 |
| D | 3/3 | 1.0 | 307 | 289 |

| pH | máx. (nm) | centroide (nm) |
|---|---|---|
| 2 | 307 | 283 |
| 4.65 | 297 | 279 |
| 7.4 | 274 | 270 |
| 10 | 275 | 271 |
| 12 | 304 | 287 |

Previsão qualitativa: **desvio para o azul (~307 -> 274 nm) ao desprotonar a
carboxila (pKa calc. 4.91)** e retorno ao vermelho acima de pH ~11 (diânion,
pKa_AD 11.2). É a assinatura a ser comparada com a titulação (pKa exp. 4.65
por posição de banda no UV próximo). **Comparação NÃO feita:** não há
espectros experimentais nesta máquina (os resumos de Camuri/Pazin só dão o
pKa 4.65). Necessário: máximos de banda (ou espectros digitalizados) vs pH.
Limites: σ=0.4 eV é largo; N usa só 50% da população e A 27% (incerteza de
cobertura acima); TD-DFT HSE06 tem erro sistemático de ~0.2–0.3 eV que
desloca todas as bandas, então comparar DESLOCAMENTOS, não posições absolutas.

**Tarefa 3 (2ª referência −1→−2):** job 5217096 (p-cumarato mono/di) rodando
(jobex ciclo 2 em 2026-09-25); ao terminar: adicionar linha
`4-coumarate_phenol,-1to-2,...,9.45` (fonte única Nordstrom & Lindberg 1965)
e refazer `summarize_pka_calibration.py` (classe −1→−2 passa de N=1 a N=2).
pKa verificados em fonte primária (dataset IUPAC Serjeant–Dempsey/Perrin):
2,6-dimetilfenol 10.59–10.63; 4-hidroxibenzoato (fenol) revisado para 9.36
[9.32, 9.49] (resíduo −4.462). Aviso de CSV: nomes sem vírgula.

### 9.38. Comparação com espectro experimental (1 pH) + cobertura A/dianion (2026-09-25)

Arquivo do usuário `artepillin_C_studies/uvvis_conf1_vs_exp.dat` (250–400 nm;
colunas calc_conf1_norm, calc_ensemble_norm, exp_norm; procedência/pH/estado
NÃO documentados no arquivo -- assumido neutro). Máx. calc (ensemble) 304.5 nm,
conf1 304.0 nm; máx. exp 310.0 nm -> **−5.5 nm (~−0.07 eV)**. Melhor ajuste
(alargamento gaussiano + deslocamento rígido, RMSE 250–330 nm): σ≈0.45 eV,
+10 nm, RMSE 0.05 (0.49 sem alargar). Discrepâncias de forma: FWHM exp ~65 nm
vs calc ~10 nm (alargamento do cálculo, não erro eletrônico); ombro exp 260–300
nm (0.4–0.9) e cauda vermelha 330–400 nm (0.7→0.07) ausentes no cálculo (causas
não testadas: vibrônica/solvente; espécies A ~274 nm / B ~370 nm se pH não for
ácido). Não testa o deslocamento com o pH (1 pH só).

**Cobertura:** os 4 ensembles foram cortados em 95% da população GFN2 acumulada
(neutro 170/337, A 143/322, B 11/37, diânion 3/21; 95.1–95.8%). Ou seja, "3/3"
do diânion = 95.3% da população GFN2 dos 21 confôrmeros CREST, NÃO cobertura
DFT. Monoânion A: top-N por população r2SCAN-3c cumulativa = 63% (N=20; +15 novos
sobre os 20 reais), 75% (30; +25), 82% (40; +35), 91% (60; +55); só 5 dos top-20
já são reais (o rotulado foi escolhido por aprendizado ativo, não por população).
Caveat aberto: CREST/GFN2 foi em vácuo -- representatividade do ensemble dos
ânions (sobretudo o diânion) em água NÃO testada.

### 9.39. Teste de ranking em água (GFN2/ALPB) + extensão top-30 do monoA (2026-09-25)

**Extensão do monoânion A (por POPULAÇÃO, não diversidade):** top-30 por
população de Boltzmann r2SCAN-3c (73.3%); 25 novos (`monoA_top30_new_jobs.txt`),
diretório `layer_b/top30_ext/{labeled_jobs,cosmo_jobs}` (cópias dos
`ridft_only_jobs` -- os dados do modelo B1 NÃO foram alterados). Plano SEM aoforce
(~32 núcleo-h/conf): ridft + escf 60 estados (~8) + COSMO (~1); IV via modelo
B1 (LOO 0.77 cm⁻¹); G_térmico transferido dos 20 reais (sd 0.15 kcal/mol,
faixa 0.56). Jobs: prep COSMO 5217396 -> escf 5217397 (dependência afterok) ->
COSMO ridft (submetido automaticamente após o prep).

**Reotimização GFN2/ALPB(água) de todos os confôrmeros CREST** (vácuo → água;
`job_xtb_alpb_rerank.sh`, `analyze_alpb_rerank.py`; arrays 5217422/3/4):

| Estado | pop. em água retida pelos confôrmeros do corte de 95% em vácuo | confôrmeros p/ 95% em água | entrariam | dist. de variação total vácuo↔água | Rg pond. vácuo→água (Å) |
|---|---|---|---|---|---|
| Diânion | 0.775 | 6 | 7, 4, 5 | 0.229 | 4.143→4.214 |
| Monoânion B | 0.875 | 12 | 12, 24, 26 | 0.192 | 4.212→4.281 |
| Monoânion A | 0.902 | 155 | 40 novos | 0.129 | 4.012→4.061 |

Leitura: a água comprime a dispersão de energias (grupos a +1.9 kcal/mol em
vácuo ficam a ~+0.9), então o corte em vácuo perde 10–22% da população em
água. **Hipótese do dobramento em vácuo NÃO sustentada** (Rg muda ≤1.7%).
Impacto estimado no pKa (ΔG_ens ≈ −RT ln(pop. retida)): diânion −0.15, B −0.08,
A −0.06 kcal/mol → δpKa ≤0.1 (A→D ≈ −0.07, B→D ≈ −0.05) -- **não explica a
lacuna pK(A→D)−pK(N→B)≈3.0**. Limites: nível GFN2/ALPB (erros de ~1 kcal/mol
plausíveis no ranking), degenerescência CREST mantida após reotimização, energias
DFT dos confôrmeros novos assumidas iguais às dos conformeros retidos, neutro
não testado.

### 9.40. Ampliação do ensemble guiada pela água (aprovada pelo usuário, 2026-09-25; em andamento)

Decisão do usuário: ampliar os ensembles com os confôrmeros necessários em
água e levar a comparação vácuo × água ao paper (mostra que a hipótese foi
considerada). Escopo: diânion (CREST 7, 4, 5) e monoânion B (CREST 12, 24, 26);
monoânion A segue só com o top-30 por população (os 40 confôrmeros extras em água
exigiriam otimização r2SCAN-3c própria; NÃO feito, opcional).
Consistência: mesma geometria inicial (CREST em vácuo), mesmo template e mesmo
teto de 10 ciclos do Layer A (`setup_water_extension.py`,
`layer_a/water_ext/{conformer_jobs,map.csv}`), r2SCAN-3c → HSE06 (ridft +
aoforce + escf) → COSMO. Pipeline encadeado em `chain_water_ext.sh`
(log `chain_water_ext.log`; opt 5217554). Depois: freeh, e pesos do ensemble
recalculados por G_aq (E_HSE06 + G_térmico + G_solv COSMO) sobre TODOS os
confôrmeros (mais fisico que pesos r2SCAN-3c em gás), com comparação vácuo×água
do pKa. Figura/tabela para o paper: `vacuum_vs_water_ensembles.png/.dat`
(`plot_vacuum_vs_water_ensembles.py`).

**Adendo 9.40 (usuário pediu, 2026-09-25):** fechar a ponta solta do monoânion A --
os 40 confôrmeros CREST fora do conjunto original que entram no corte de 95% em
água (índices em `monoA_water_ext_indices.txt`) foram montados por
`setup_water_extension.py` (`layer_a/water_ext/`, mesmo template/10 ciclos das
demais) e otimizados em r2SCAN-3c (fase gasosa, como toda a Camada A -- interpretado
assim; NÃO otimização com solvente). Job 5217609. Próximo passo: adicionar as
energias r2SCAN-3c à tabela de populações e refazer o teste de representatividade
(como nas §9.37) para o monoA; HSE06 nesses 40 só se o teste mostrar necessidade.
