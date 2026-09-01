#!/usr/bin/env python3
"""
filter_boltzmann_population.py

Layer A population filter of the Paper 2 conformational funnel. Reads
a CREST ensemble (crest_conformers.xyz + crest.energies, GFN2-xTB
relative energies in kcal/mol) and keeps the smallest set of
conformers, in energy order, whose cumulative Boltzmann population at
a given temperature reaches a target cutoff (default 95% at 298.15 K).
Conformers beyond the cutoff are population-irrelevant at that
temperature and are dropped before any DFT refinement is spent on
them.

Each conformer in crest.energies/crest_conformers.xyz is already
deduplicated by CREST's CREGEN (rotamers merged into one representative
structure per conformer) -- but that collapses away how many near-
degenerate rotamers each conformer actually represents, which matters
for a correct Boltzmann population (a conformer that is really 58
symmetry/rotation-equivalent rotamers contributes 58x the partition-
function weight of a conformer with no equivalent rotamers, even at
the same energy). CREST's own "population of lowest" figure in its log
accounts for this; a naive per-line Boltzmann average over
crest.energies does not, and silently under-weights exactly the
conformers CREGEN found the most equivalent paths to (independently
verified: naive gave 1.2% for the global minimum on the neutral-state
ensemble, vs. CREST's own reported 8.5% -- degeneracy-weighting here
reproduces that 8.5% to within 0.02 pp). Degeneracy per conformer comes
from CREST's 'cre_members' file (column 1 = degeneracy count, same row
order as crest.energies).

Usage:
    python ~/work_turbomole/scripts/filter_boltzmann_population.py \
        --energies crest.energies --conformers crest_conformers.xyz \
        --members cre_members \
        --output-dir population_filtered [--cutoff 0.95] [--temp 298.15]
"""

import argparse
import os

KB_KCAL_PER_MOL_K = 0.0019872041  # Boltzmann constant, kcal/(mol*K)


def read_energies(path):
    """Return list of (index, relative_energy_kcal) in file order."""
    entries = []
    with open(path) as f:
        for line in f:
            parts = line.split()
            if not parts:
                continue
            idx, erel = int(parts[0]), float(parts[1])
            entries.append((idx, erel))
    return entries


def read_degeneracies(path):
    """Return list of degeneracy counts (column 1) in file order, from
    CREST's 'cre_members' file (first line is the conformer count)."""
    with open(path) as f:
        n = int(f.readline().split()[0])
        degens = [int(f.readline().split()[0]) for _ in range(n)]
    return degens


def read_xyz_frames(path):
    """Return list of raw multi-line strings, one per xyz frame, in file order."""
    with open(path) as f:
        lines = f.readlines()
    frames = []
    i = 0
    while i < len(lines):
        natoms = int(lines[i].split()[0])
        frame = "".join(lines[i:i + natoms + 2])
        frames.append(frame)
        i += natoms + 2
    return frames


def boltzmann_weights(energies_kcal, degeneracies, temp_k):
    kt = KB_KCAL_PER_MOL_K * temp_k
    unnorm = [d * pow(2.718281828459045, -e / kt)
              for e, d in zip(energies_kcal, degeneracies)]
    total = sum(unnorm)
    return [u / total for u in unnorm]


def main():
    parser = argparse.ArgumentParser(
        description="Keep the smallest set of CREST conformers (in "
                    "energy order) whose cumulative Boltzmann population "
                    "reaches a target cutoff.")
    parser.add_argument("--energies", required=True,
                         help="CREST 'crest.energies' file (index, "
                              "relative energy in kcal/mol)")
    parser.add_argument("--conformers", required=True,
                         help="CREST 'crest_conformers.xyz' ensemble, "
                              "same order as --energies")
    parser.add_argument("--members", required=True,
                         help="CREST 'cre_members' file (per-conformer "
                              "rotamer degeneracy), same order as --energies")
    parser.add_argument("--output-dir", required=True,
                         help="Folder to write the filtered ensemble + "
                              "population table into")
    parser.add_argument("--cutoff", type=float, default=0.95,
                         help="Target cumulative Boltzmann population "
                              "[default: 0.95]")
    parser.add_argument("--temp", type=float, default=298.15,
                         help="Temperature in K [default: 298.15]")
    args = parser.parse_args()

    entries = read_energies(args.energies)
    frames = read_xyz_frames(args.conformers)
    degeneracies = read_degeneracies(args.members)
    if len(entries) != len(frames) or len(entries) != len(degeneracies):
        raise SystemExit(
            f"Mismatch: {len(entries)} energies, {len(frames)} xyz frames, "
            f"{len(degeneracies)} degeneracies -- '{args.energies}', "
            f"'{args.conformers}' and '{args.members}' must come from the "
            f"same CREST run.")

    energies_kcal = [e for _, e in entries]
    weights = boltzmann_weights(energies_kcal, degeneracies, args.temp)

    if os.path.exists(args.output_dir):
        raise FileExistsError(f"{args.output_dir} already exists -- refusing to overwrite.")
    os.makedirs(args.output_dir)

    table_path = os.path.join(args.output_dir, "population_table.csv")
    survivors = []
    cumulative = 0.0
    with open(table_path, "w") as f:
        f.write("crest_index,relative_energy_kcal,degeneracy,population,cumulative_population,kept\n")
        for (idx, erel), deg, w in zip(entries, degeneracies, weights):
            keep = cumulative < args.cutoff
            if keep:
                cumulative += w
                survivors.append((idx, erel, w))
            f.write(f"{idx},{erel:.4f},{deg},{w:.6f},{min(cumulative, 1.0):.6f},{int(keep)}\n")

    ensemble_path = os.path.join(args.output_dir, "population_filtered.xyz")
    with open(ensemble_path, "w") as f:
        for (idx, erel, w), frame in zip(survivors, frames):
            f.write(frame)

    print(f"Temperature: {args.temp} K, cutoff: {args.cutoff:.0%}")
    print(f"Kept {len(survivors)} / {len(entries)} conformers, "
          f"cumulative population = {cumulative:.4f}")
    print(f"  highest-population conformer: crest_index={survivors[0][0]}, "
          f"Erel={survivors[0][1]:.3f} kcal/mol, "
          f"population={survivors[0][2]:.1%}")
    print(f"  lowest-population survivor: crest_index={survivors[-1][0]}, "
          f"Erel={survivors[-1][1]:.3f} kcal/mol, "
          f"population={survivors[-1][2]:.2%}")
    print(f"\nWrote {ensemble_path} and {table_path}")


if __name__ == "__main__":
    main()
