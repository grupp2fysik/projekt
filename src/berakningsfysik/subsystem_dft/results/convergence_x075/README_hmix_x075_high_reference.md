# Delta H_mix high-reference test for TiAlN x = 0.75

## Purpose

This test was added after feedback from Adrian.

The goal was to compare the selected cheap production setup against a more accurate reference setup where both cutoff and k-point density are increased at the same time.

## Selected production setup

The selected cheap setup is:

- ecutwfc = 60 Ry
- ecutrho = 240 Ry
- TiAlN k-points: 2 2 2
- TiN/AlN k-points: 4 4 4

For x = 0.75 this gave:

- Delta H_mix = 481.014 meV/formula unit

## High-reference setup

The high-reference setup is:

- ecutwfc = 100 Ry
- ecutrho = 400 Ry
- TiAlN k-points: 4 4 4
- TiN/AlN k-points: 8 8 8

This was run only for x = 0.75 as an additional convergence check.

## Energies used

TiN:

- output: dft-sigma/TiN/output/scf/TiN_scf_ecut100_k888.out
- total energy: -848.68393036 Ry
- formula units: 4

AlN:

- output: dft-sigma/AlN/output/scf/AlN_scf_ecut100_k888.out
- total energy: -271.61329603 Ry
- formula units: 4

TiAlN x = 0.75:

- output: dft-sigma/TiAlN/output/scf/TiAlN_x075_scf_ecut100_k444.out
- total energy: -3325.92022291 Ry
- formula units: 32

## Result

High-reference result:

- Delta H_mix = 0.479352 eV/formula unit
- Delta H_mix = 479.352 meV/formula unit

Selected production setup result:

- Delta H_mix = 481.014 meV/formula unit

Difference:

- 1.662 meV/formula unit
- about 0.347 %

## Interpretation

The high-reference test confirms the earlier trend:

- the cutoff is already well converged at 60 Ry
- the k-point density has a larger effect than the cutoff
- the selected cheaper setup differs from the high-reference result by about 0.35 %

Because the project has limited CPU time, the selected production setup is still reasonable for calculating more x-values.

## Files

Input data for calculation:

- hmix_high_reference_input.csv

Python script:

- calculate_hmix_high_reference.py

Saved result:

- hmix_results_x075_high_reference.txt

Energy source summary:

- hmix_energy_sources_x075_high_reference.txt
