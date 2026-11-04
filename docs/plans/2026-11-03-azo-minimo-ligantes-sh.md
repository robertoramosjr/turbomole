# Azobenzeno: mínimo do trans → ligantes → SH com ligantes (2026-11-03)

Ordem combinada com o usuário em 2026-11-03. Complementa o plano v2 (`azobenzeno/docs/plans/2026-10-29-azobenzeno-plano-v2.md`, na access): B0, B3, Bv e §5 itens 2 e 7. Atualizar o estado aqui, no lugar; decisão nova vira uma linha em "Decisões".

**Nível de teoria (inalterado):** CAM-B3LYP-D3(BJ)/def2-TZVP, `$senex`, TM 8.0.1. Dinâmica: CAM-B3LYP/def2-SVP-TDA.
**Máquina:** harpia (`module load modmat/turbomole/8.0.1`), partição `par16`, 16 núcleos/job. Limite de RAM novo: 900 MB por CPU, logo `--mem=14400M` (16 núcleos). Consumo real do TM nunca medido aqui: conferir `MaxRSS` no primeiro job. Os `job.sh` antigos são da access (`short`, `$HOME/software/TURBOMOLE801`) e precisam ser adaptados.
**Regra:** custo em núcleo-h reportado ao usuário antes de cada `sbatch`; validar por texto (`GEO_OPT_CONVERGED`, `all done`, sem `ended abnormally`).

## Passo 1 — mínimo do trans (B0)
Contexto: cis é mínimo verdadeiro (+20,5 cm⁻¹). Trans: modo imaginário −33,7 → −32,9 → −29,2 → −19,5 (t1) → −23,2 cm⁻¹ (t2, grade 5). t1/t2 deslocaram 0,4 bohr e reotimizaram com `$statpt` do Paper 1 (thrmaxgrad 5e-4, thrrmsgrad 2e-4, thrmaxdispl 5e-4, thrrmsdispl 2e-4, threchange 1e-7). t3 foi montado sem `$senex` (energia 0,12 Eh acima, geometria congelada): descartado.
**Hipótese principal:** com ~20 cm⁻¹ a curvatura é ~1e-4 Eh/bohr² (estimativa grosseira, massa reduzida de alguns amu); 0,4 bohr dá gradiente ~4e-5, ~10× abaixo de `thrmaxgrad`. O `jobex` não enxerga o vale.
- [~] 1a. (RODADO 2026-11-03, `azo_trans_camb3lyp_tzvp_t4_tightopt`; job 344 travou em 209 ciclos com thrmaxgrad 1e-5: energia idêntica nos 10 dígitos, |dE/dxyz| 3,9e-5, MAX grad 1,3e-5 = piso de ruído da grade 5. Reprocessado com thrmaxgrad 2e-5 (job 347: opt em 1 ciclo + ridft + aoforce); aoforce do t2 levou 41 min/16 núcleos ≈ 11 núcleo-h, não 5). Aplicar thrmaxgrad ≈ 2e-5 (não 1e-5) nos ligantes. Original: copiar o t2 (`azo_trans_camb3lyp_tzvp_opt_disp_tight_t2`, partir de `coord`/`mos` finais) para pasta nova; `$statpt`: thrmaxgrad ≈ 1e-5, thrrmsgrad ≈ 5e-6, thrmaxdispl e thrrmsdispl ≈ 1e-4, threchange 1e-9; manter grade 5, `$scfconv 10`, `$denconv 1d-10`, `$senex`. `jobex -c 800`, depois `aoforce`.
- [ ] 1b. Varredura 1D (`displace_along_mode.py`, modo 1 do t2): ±1, ±2, ±3 bohr + ponto central, `ridft` com SCF apertado. Energia cai nos dois lados = sela real; poço plano dentro de ~1e-6 Eh = ruído.
- [ ] 1c. Se persistir: hessiana exata do aoforce no `$statpt` (`hssfreq`, `hssupdate`, manual 8.0 l.14439–14460), depois C2h/C2 (o modo quebra a planaridade?), depois ruído do `$senex` (grade maior, energia pontual).
- Saída: trans sem modo imaginário, ou modo declarado como artefato com a varredura como evidência. Fecha o plano v2 §5 item 2.
- Custo (estimativa, não medida): 1a ≈ 10–60 núcleo-h (critérios apertados pedem mais ciclos), 1b ≈ 3 núcleo-h, aoforce trans ≈ 5 núcleo-h.

## Passo 2 — os 4 ligantes
`azobenzeno/anchored/{cis,trans}_{sh,nh2}` (4-monossubstituídos, só `coord` idealizado de `add_azo_anchors.py`, nada calculado). Nunca os dois grupos na mesma molécula.
- [ ] 2a. `define` via template CAM-B3LYP-D3(BJ)/def2-TZVP + `$senex` (conferir o `control` à mão); critérios do passo 1a desde o início (modos rasos de S–H e N–H2).
- [ ] 2b. `jobex -c 800` → `ridft` → `aoforce` (sem imaginários; gate de estabilidade da regra 7 não é exigido para singleto de gap grande, mas conferir HOMO–LUMO).
- [ ] 2c. Novo Wigner a partir da hessiana de cada ligante (cis 500, trans 1000, como no azo puro; base da hessiana: decisão do usuário, §5 item 5 do plano v2).
- Custo: estimativa anterior ~160–320 núcleo-h para 4 moléculas (feita para HSE06; CAM e tight podem subir ×1,3–2). Aprovar antes do `sbatch`.

## Passo 3 — SH com ligantes
- [ ] 3a. Trazer/validar na harpia os scripts do SH (`seed_sh_log_groups.py`, `prepare_tm_s0_continuation_v3.py`, `job_tm_native_sh_s0_v3.sh`, `summarize_sh_traj.py`, `aggregate_sh.py`); `tm_native_sh/` (PROBLEMAS.md, JOBS.md, lotes) só existe na access.
- [ ] 3b. cis→trans 420 nm com cada ligante (protocolo congelado: gap 0,2, `iter300`, dt 0,5 fs, NSH trans 600, NS0 300, corte |CNNC| < 90° = cis, átomos CNNC a rever: o índice muda com os ligantes).
- [ ] 3c. trans→cis 425 nm em blocos de < 255 ciclos: orquestrador (P00a etapa 3) ainda não existe; `$nsteps` restante + poda de `$nac`/`$excitation`; só vale se os blocos chegarem a 200 quadros. Piloto trans 5×400 antes do lote.
- Gates de custo do plano 2026-10-19 (piloto ~1000 + lote ~8600 núcleo-h) continuam valendo e precisam do usuário.

## Decisões
- 2026-11-03 (usuário): seguir nesta ordem; trans (não cis) é a sela; ligantes já entram no SH porque mudam a transferência de carga.

## Estado
Passo 1a em andamento (job 347). CP2K roda na access (decisão do usuário 2026-11-03); a harpia fica só com Turbomole. Repo de scripts clonado em `~/work_turbomole/scripts` (harpia); cofre: `02_Memory/turbomole/azobenzeno.md` (Pendente).

## Registro 2026-11-04 — passo 1c lido (harpia 390–392) e teste `410`
- NumForce (390): sem imaginária (24,45; 65,57; 85,17; 88,98 cm⁻¹); aoforce g7 (392) e statpt+aoforce (391): −23,3. O modo imaginário do aoforce é o modo 4 do NumForce (+89,0), overlap 0,97. Hipótese: artefato da 2ª derivada analítica (CAM-B3LYP); NumForce = referência para ZPE/G.
- Termodinâmica NumForce (298,15 K, σ=1, RRHO próprio): ZPE 0,1928060 Eh, G−E_el 96,87 kcal/mol, G −572,5193464 Eh.
- Teste: job `410` (`azo_trans_t4_curvA`, array 18) = ΔE e gradiente ao longo do modo A (±1,2,3 amu^½·bohr). Próximo: ler, ajustar ω², registrar; se confirmado, hessiana por NumForce nos passos 2 (ligantes) e Wigner.

## Estado 2026-11-04 (noite)
- Passo 1 FECHADO: trans t4 é mínimo (curvatura +89 cm⁻¹, job 410); `aoforce` −23 = artefato analítico; NumForce = referência (freeh 429: G−E_el 96,86 kcal/mol). Falta só a causa raiz: job `430` (B3LYP aoforce, mesma geometria/grade).
- Próximo: passo 2 (ligantes) com NumForce no lugar de aoforce para a hessiana/ZPE/Wigner; custo ao usuário antes do `sbatch`.
