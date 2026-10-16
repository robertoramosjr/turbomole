#!/usr/bin/env python3
"""
apply_b2_osc_surrogate.py

Applies B2's confirmed oscillator-strength/energy delta-learning
results to a full protonation-state ensemble, producing one row per
conformer with the predicted S1 energy + oscillator strength, which
model produced it, and every caveat needed before it enters a final
ranking.

v2 (2026-09-11, Sec. 9.25-9.26 of
protonation_and_conformer_investigation.md): supersedes the v1 flat
log10(f) offset with the validated log-log SLOPE+INTERCEPT model
(log10(f_HSE06) = --f-model-intercept + --f-model-slope *
log10(f_sTDA), slope significantly != 1, LOO Q2=0.52, permutation-
confirmed, out-of-sample error reduced ~65%% vs. the flat offset --
see test_b2_osc_magnitude_recalibration.py). Also accepts MULTIPLE
--labeled-jobs / --hse06-csv sources, since real HSE06 confirmations
now come from several batches (original labeled subset + extreme-value
confirmation batch + random-control batch), not just one.

Adds an explicit, UNIFORM background risk that is NOT contrast-
dependent: rigorously checking every confirmed conformer (production
max-f-in-window pick vs. the true character-matched state) found 3/36
genuine state-selection mismatches (job_0023, job_0035, job_0036), ALL
at f-contrast 2.5-3.6x -- solidly inside the "confident" (>2.0x)
bucket. This proves f-contrast cannot detect this failure mode at any
threshold (a screen across the other confident conformers for the same
orbital signature also came back 100% positive, i.e. no discriminating
power -- see Sec. 9.24 Task 3). So --state-selection-mismatch-rate (and
its Wilson CI bounds) is applied as a flat, contrast-independent
per-conformer risk note to every surrogate-predicted row, NOT as a
function of contrast. The separate --exclude-contrast /
--inflate-contrast flags remain for a genuinely different, still-valid
phenomenon (near-tied competing states -- Sec. 9.17/9.20) and are not
removed.

Per-conformer resolution order:
  1. Real HSE06 label found in any --labeled-jobs dir -> use it
     directly (method=hse06_ground_truth), no surrogate, no mismatch
     risk note (the true value is known, not predicted).
  2. Listed in --override-csv (manually confirmed special cases) ->
     the manually-forced (energy, osc) pair, then through the same
     energy offset / f slope-intercept model as any other surrogate
     conformer (method=surrogate_forced_state).
  3. f-contrast at or below --exclude-contrast -> EXCLUDED pending
     manual confirmation (method=excluded_low_confidence).
  4. Otherwise -> standard surrogate: confidence CSV's max-f-in-window
     value, energy offset + f slope-intercept model applied.
     uncertain_inflated=True if contrast <= --inflate-contrast.

Usage:
    python ~/work_turbomole/scripts/apply_b2_osc_surrogate.py \
        --state-label neutral \
        --confidence-csv layer_b/stda_confidence_full_ensemble.csv \
        --hse06-csv layer_b/s1_hse06.csv layer_b/s1_hse06_confirmation_batch.csv \
                    layer_b/s1_hse06_random_control_batch.csv \
        --labeled-jobs layer_b/labeled_jobs layer_b/confirmation_batch_jobs \
                       layer_b/random_control_batch_jobs \
        --override-csv neutral_special_overrides.csv \
        --offset-mean-energy 1.2113 --offset-std-energy 0.1723 \
        --f-model-intercept -0.2027 --f-model-slope 0.5794 \
        --exclude-contrast 1.26 --inflate-contrast 2.0 \
        --state-selection-mismatch-rate 0.0833 \
        --state-selection-mismatch-ci-low 0.0287 --state-selection-mismatch-ci-high 0.2183 \
        --osc-floor 1e-4 \
        --output layer_b/b2_osc_predictions_full_ensemble_v2.csv
"""

import argparse
import csv
import math
import os


def main():
    parser = argparse.ArgumentParser(
        description="Apply B2's oscillator-strength/energy delta-learning "
                    "result (v2: log-log slope+intercept for f, not a flat "
                    "offset) to a full protonation-state ensemble.")
    parser.add_argument("--state-label", required=True)
    parser.add_argument("--confidence-csv", required=True,
                         help="characterize_stda_selection_confidence.py output")
    parser.add_argument("--hse06-csv", nargs="+", required=True,
                         help="One or more batch_extract_s1_hse06.py outputs "
                              "(all confirmed HSE06 sources -- merged)")
    parser.add_argument("--labeled-jobs", nargs="+", required=True,
                         help="One or more folders of job_NNNN/ with real "
                              "HSE06 labels -- any job_dir found in any of "
                              "them bypasses the surrogate")
    parser.add_argument("--override-csv", default=None)
    parser.add_argument("--offset-mean-energy", type=float, required=True)
    parser.add_argument("--offset-std-energy", type=float, required=True)
    parser.add_argument("--f-model-intercept", type=float, required=True,
                         help="a in log10(f_HSE06) = a + b*log10(f_sTDA)")
    parser.add_argument("--f-model-slope", type=float, required=True,
                         help="b in log10(f_HSE06) = a + b*log10(f_sTDA)")
    parser.add_argument("--exclude-contrast", type=float, default=1.26)
    parser.add_argument("--inflate-contrast", type=float, default=2.0)
    parser.add_argument("--state-selection-mismatch-rate", type=float, required=True,
                         help="Flat, contrast-independent background rate of "
                              "genuine state-selection mismatch (job_0023-"
                              "type), applied to every surrogate row as a "
                              "note, NOT used to change the predicted value")
    parser.add_argument("--state-selection-mismatch-ci-low", type=float, required=True)
    parser.add_argument("--state-selection-mismatch-ci-high", type=float, required=True)
    parser.add_argument("--osc-floor", type=float, default=1e-4)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    labeled_dirs = set()
    for ld in args.labeled_jobs:
        if os.path.isdir(ld):
            labeled_dirs |= {d for d in os.listdir(ld)
                              if os.path.isdir(os.path.join(ld, d))}

    hse06 = {}
    for path in args.hse06_csv:
        with open(path) as f:
            for row in csv.DictReader(f):
                hse06[row["job_dir"]] = (float(row["s1_energy_eV"]), float(row["s1_osc_len"]))

    overrides = {}
    if args.override_csv:
        with open(args.override_csv) as f:
            for row in csv.DictReader(f):
                overrides[row["job_dir"]] = (
                    float(row["forced_energy_eV"]), float(row["forced_osc"]),
                    row["reason"],
                )

    with open(args.confidence_csv) as f:
        confidence = {row["job_dir"]: row for row in csv.DictReader(f)}

    all_job_dirs = sorted(set(confidence) | labeled_dirs)
    rows = []
    for jd in all_job_dirs:
        if jd in labeled_dirs and jd in hse06:
            e, osc = hse06[jd]
            rows.append({
                "job_dir": jd, "predicted_energy_eV": e, "predicted_osc": osc,
                "energy_std": 0.0, "log10_f_std": "",
                "method": "hse06_ground_truth", "special_treatment": "",
                "f_contrast_ratio": "", "uncertain_inflated": False,
                "excluded": False, "state_selection_mismatch_rate": "",
                "state_selection_mismatch_ci": "",
            })
            continue

        if jd not in confidence:
            print(f"WARNING: {jd} has no confidence-csv entry and no HSE06 "
                  f"label, skipping")
            continue
        contrast_str = confidence[jd]["f_contrast_ratio"]
        contrast = float(contrast_str) if contrast_str != "inf" else float("inf")

        if jd in overrides:
            e_stda, osc_stda, reason = overrides[jd]
            method = "surrogate_forced_state"
        elif contrast <= args.exclude_contrast:
            rows.append({
                "job_dir": jd, "predicted_energy_eV": "", "predicted_osc": "",
                "energy_std": "", "log10_f_std": "",
                "method": "excluded_low_confidence", "special_treatment": "",
                "f_contrast_ratio": f"{contrast:.3f}", "uncertain_inflated": False,
                "excluded": True, "state_selection_mismatch_rate": "",
                "state_selection_mismatch_ci": "",
            })
            continue
        else:
            e_stda = float(confidence[jd]["s1_energy_eV_stda"])
            osc_stda = float(confidence[jd]["s1_osc_stda"])
            reason = ""
            method = "surrogate_max_f"

        predicted_energy = e_stda + args.offset_mean_energy
        log10_f_stda = math.log10(max(osc_stda, args.osc_floor))
        predicted_log10_f = args.f_model_intercept + args.f_model_slope * log10_f_stda
        predicted_osc = 10 ** predicted_log10_f
        uncertain_inflated = contrast <= args.inflate_contrast

        rows.append({
            "job_dir": jd, "predicted_energy_eV": f"{predicted_energy:.4f}",
            "predicted_osc": f"{predicted_osc:.5f}",
            "energy_std": f"{args.offset_std_energy:.4f}",
            "log10_f_std": "model_residual_see_Sec_9.25",
            "method": method, "special_treatment": reason,
            "f_contrast_ratio": f"{contrast:.3f}" if contrast != float("inf") else "inf",
            "uncertain_inflated": uncertain_inflated, "excluded": False,
            "state_selection_mismatch_rate": args.state_selection_mismatch_rate,
            "state_selection_mismatch_ci": (
                f"[{args.state_selection_mismatch_ci_low:.4f},"
                f"{args.state_selection_mismatch_ci_high:.4f}]"),
        })

    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    n_gt = sum(1 for r in rows if r["method"] == "hse06_ground_truth")
    n_forced = sum(1 for r in rows if r["method"] == "surrogate_forced_state")
    n_excluded = sum(1 for r in rows if r["excluded"])
    n_std_confident = sum(1 for r in rows if r["method"] == "surrogate_max_f"
                            and not r["uncertain_inflated"])
    n_std_inflated = sum(1 for r in rows if r["method"] == "surrogate_max_f"
                           and r["uncertain_inflated"])
    print(f"=== {args.state_label} (N={len(rows)}) ===")
    print(f"  hse06_ground_truth:                        {n_gt}")
    print(f"  surrogate_forced_state:                     {n_forced}")
    print(f"  surrogate_max_f, confident (contrast>{args.inflate_contrast}):        {n_std_confident}")
    print(f"  surrogate_max_f, uncertain_inflated "
          f"({args.exclude_contrast}<contrast<={args.inflate_contrast}): {n_std_inflated}")
    print(f"  excluded_low_confidence:                    {n_excluded}")
    n_surrogate = n_forced + n_std_confident + n_std_inflated
    print(f"\n  NOTE: all {n_surrogate} surrogate-predicted rows (forced + "
          f"confident + inflated) carry a flat "
          f"state_selection_mismatch_rate={args.state_selection_mismatch_rate:.4f} "
          f"(95% CI [{args.state_selection_mismatch_ci_low:.4f}, "
          f"{args.state_selection_mismatch_ci_high:.4f}]) -- this is NOT "
          f"contrast-dependent (proven undetectable by contrast, Sec. 9.25) "
          f"and is NOT reflected in the point prediction, only as a caveat "
          f"for downstream uncertainty propagation.")
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
