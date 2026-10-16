# Protocolo de seleção de features/modelo para delta-learning (padrão do projeto)

## O que é isto

Um protocolo de 4 etapas, validado na prática durante o fechamento dos
modelos B1 de monoânion A/B (ver
[protonation_and_conformer_investigation.md](protonation_and_conformer_investigation.md),
Seções 9.7–9.13), pra decidir qual modelo usar quando se tem um alvo
caro (Δ = valor DFT − valor barato) e um punhado de features candidatas
sobre uma amostra pequena (N=10–30 pontos rotulados, típico de uma
etapa inicial de aprendizado ativo). **Reutilizar este protocolo
integralmente** ao fechar o B2 e ao calibrar pKa — não reinventar a
decisão a cada vez.

O motivo de existir: um GPR com várias features e poucos pontos rotulados
super-parametriza facilmente (visto na prática: o monoânion B tinha um
sinal linear real de Q²=0.594 com 1 feature, mas o GPR de 2 features
dava Q²=−0.061 nos mesmos dados — o modelo mais simples venceu por uma
margem grande). Este protocolo força passar pelas alternativas mais
simples antes de qualquer coisa mais flexível, e calibra o resultado
final contra o risco de comparação múltipla antes de aceitar qualquer
"sinal" como real.

## As 4 etapas (nesta ordem, não pular nenhuma)

### 1. Correlação individual de cada feature candidata vs. alvo

Correlação de Pearson simples (`scipy.stats.pearsonr`) de CADA feature
isoladamente contra o Δ-alvo. Reporta r e p-valor pra toda feature, não
só as que parecem promissoras — uma relação puramente não-linear pode
escapar do Pearson, mas isso já reduz a expectativa de que um modelo
linear vá aproveitá-la bem.

**Script:** `audit_b1_features.py`

### 2. OLS de 1 feature sob validação leave-one-out (LOO)

Pra cada feature candidata: regressão linear simples (2 parâmetros —
intercepto + inclinação), avaliada por LOO usando os MESMOS folds em
todo o projeto. Reporta Q² = 1 − SS_res/SS_tot (contra a média como
baseline) — **não MAE/RMSE sozinho**, e certamente não confiar só no ρ
de Spearman (ver a pegadinha do artefato de reversão à média, Seção
9.8: um preditor "média dos outros" sob LOO dá ρ=−1 exato sempre,
mesmo sem sinal nenhum).

Q² > 0 é a barra pra "essa feature sozinha tem sinal real aprendível".

**Scripts:** `test_b1_linear_baseline.py` (LOO-OLS + diagnóstico de
leverage/Cook's distance), `diagnose_b1_surrogate.py` (mesma lógica,
formato mais amplo incluindo o controle de reversão à média).

Para velocidade em muitas repetições (ex. teste de permutação da etapa
4), existe a versão fechada (exata, não precisa refit em loop) do LOO
de regressão de 1 feature via matriz hat — ver
`permutation_test_b1_features.py::loo_q2_closed_form`.

### 3. Combinação de features — só se justificada

Só tenta combinar 2+ features se **nenhuma isolada** já resolver bem
(Q²>0 razoável). Se duas ou mais isoladas passarem de Q²>0, tenta a
combinação linear (sem termo de interação) e **compara o Q² da
combinação contra a melhor isolada**:
- combinação melhora → usa a combinação;
- combinação piora ou empata → descarta a combinação, usa só a melhor
  isolada (foi exatamente o que aconteceu no monoânion B: a combinação
  das 2 features deu Q²=0.489, pior que a distância O···H sozinha,
  Q²=0.594 — a segunda feature diluiu o sinal, não somou).

**Script:** `select_b1_feature_model.py` (automatiza as etapas 1–3
inteiras, decide e reporta a recomendação).

### 4. Calibração por permutação — obrigatória antes de aceitar qualquer sinal

Embaralha o alvo Δ (mantém as features intactas), repete o
procedimento completo das etapas 1–3 (pega o melhor Q² entre as N
features candidatas) em cada repetição, **pelo menos 500 vezes**
(idealmente 1000). A distribuição resultante é a nula de "melhor Q² em
N tentativas, só por acaso, com esse tamanho de amostra".

Reporta, não só o Q² observado:
- percentil do Q² observado dentro da distribuição nula;
- Q² correspondente ao percentil 95 da nula (o "limiar honesto");
- p-valor empírico (fração de permutações com Q² ≥ observado).

**Regra de decisão:**
- Q² observado > percentil 95 da nula → sinal real, não explicado por
  acaso mesmo testando várias features (ex. monoânion A: Q²=0.277 no
  percentil 98.3%, p=0.017 — confirmado).
- Q² observado (mesmo negativo, "nenhuma feature resolveu") cai dentro
  do corpo comum da nula (aprox. percentil 10–90) → resultado
  consistente com ausência real de sinal, não com features ruins ou
  processo malfeito (ex. neutro: percentil 72.4%, p=0.276 —
  confirmado, fecha como "sem modelo recuperável", não como "precisa
  tentar mais features").
- Q² observado no percentil <10 da nula (pior que quase todo ruído
  puro) → anomalia genuína, investigar os dados antes de aceitar
  qualquer conclusão.
- Zona cinzenta (percentil 90–95): sinalizar explicitamente como
  fronteiriço, não tratar como seguro sem evidência adicional (mais
  pontos rotulados, ou razão mecanicista pra esperar aquela feature
  específica).

**Script:** `permutation_test_b1_features.py`

## Depois de fechar um modelo: faixa de validade

Um modelo linear fechado por este protocolo só é válido dentro da
faixa de valores da feature vista no treino. Ao aplicar em confôrmeros
novos (fora do conjunto rotulado), **sempre verificar se o valor da
feature cai dentro da faixa de treino** — fora dela é extrapolação, sem
garantia nenhuma mesmo que o Q² dentro da faixa tenha sido bom.

**Script:** `apply_b1_model.py` — aplica um modelo linear fechado
(`slope`, `intercept`, faixa de treino) a um novo conjunto de
confôrmeros, sinalizando cada previsão fora da faixa como extrapolação
explicitamente, em vez de devolver um número com a mesma confiança
implícita de uma interpolação.

## Quando (não) usar GPR

GPR com múltiplas features **nunca é a primeira tentativa**. Só
justificado se:
- (a) nenhum modelo linear (isolado ou combinado) atinge Q²>0 **e**
  há razão concreta pra suspeitar de não-linearidade (julgamento
  humano, este protocolo não decide isso sozinho); **ou**
- (b) N cresce o suficiente pra sustentar os hiperparâmetros extras
  sem repetir o problema do monoânion B original — regra prática:
  N ≥ 5× o número de hiperparâmetros do kernel (ex. GPR ARD de 2
  features tem ~4-6 hiperparâmetros → N≥20-30 mínimo razoável).

## Resumo de onde este protocolo já foi aplicado

| Uso | Resultado |
|---|---|
| B1 monoânion A | Bader(O carbonílico), Q²=0.277, calibrado (percentil 98.3%, p=0.017) |
| B1 monoânion B | distância O···H, Q²=0.594, validado por Cook's D + leave-one-out do r |
| B1 neutro | nenhuma feature passou — calibração confirma ausência real de sinal (percentil 72.4%, p=0.276) |
| B2 | ainda não aplicado |
| Calibração de pKa | ainda não aplicado |
