# hessian_check — aoforce × NumForce (azobenzeno trans, 2026-11-04)

Contexto: o `aoforce` dá −23 cm⁻¹ para o trans-azobenzeno CAM-B3LYP/def2-TZVP, o `NumForce` não.
Estes scripts documentam por que o NumForce é a referência (cofre `_log/azobenzeno.md` 2026-11-04).

1. `compare_hessian_modes.py --a azo_trans_t4_grid7_aoforce --b azo_trans_t4_numforce` — o modo imaginário é
   o mesmo modo (overlap 0,97) que o NumForce dá com +89 cm⁻¹.
2. `curvature_along_mode.py gen/analyze` — desloca só ao longo do modo; `actual -r; ridft; rdgrad` por caso
   (job array, 8 núcl., ~3 min/ponto); ΔE sobe nos dois sentidos e ω→+89 cm⁻¹ (k por gradientes).
3. `rrho_thermo.py` — conferência do `freeh` (ZPE bate com o do Turbomole em 1e-7 Eh).

Dependem de `tmhess.py` (mesma pasta). Python 3 + numpy. Hessiana: usar o bloco `$hessian (projected)`.
