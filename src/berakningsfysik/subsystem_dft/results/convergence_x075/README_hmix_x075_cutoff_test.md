# Delta H_mix cutoff-test for TiAlN x = 0.75

## Purpose

This test was done to study how the mixing enthalpy, Delta H_mix, depends on the plane-wave cutoff energy ecutwfc.

The goal was to check whether a lower cutoff can be used for later TiAlN calculations without significantly changing Delta H_mix.

## Structures used

TiAlN:
- x = 0.75
- 64-atom SQS structure
- fixed volume
- relaxed structure from input/TiAlN_x075relaxedresult.in

Reference systems:
- cubic TiN, conventional cell, 8 atoms
- cubic AlN, conventional cell, 8 atoms

TiN and AlN were first run with fixed-cell relax. The final relaxed structures were then used for SCF calculations.

## Parameters

Cutoff values tested:
- ecutwfc = 60 Ry, ecutrho = 240 Ry
- ecutwfc = 80 Ry, ecutrho = 320 Ry
- ecutwfc = 100 Ry, ecutrho = 400 Ry

K-points:
- TiAlN: 2 2 2
- TiN: 4 4 4
- AlN: 4 4 4

This follows the agreed mapping where the larger TiAlN supercell uses a smaller k-point grid than the conventional TiN/AlN cells.

## Output locations

QE SCF output files:
- dft-sigma/TiN/output/scf/
- dft-sigma/AlN/output/scf/
- dft-sigma/TiAlN/output/scf/

Energy source summary:
- hmix_energy_sources_x075_cutoff_test.txt

CSV input for H_mix calculation:
- hmix_input.csv

Python script:
- calculate_hmix.py

Saved H_mix result:
- hmix_results_x075_cutoff_test.txt

## Result

Delta H_mix values:

- ecutwfc = 60 Ry:  481.014 meV/formula unit
- ecutwfc = 80 Ry:  480.884 meV/formula unit
- ecutwfc = 100 Ry: 480.922 meV/formula unit

Difference between 60 Ry and 100 Ry:

- 0.092 meV/formula unit
- about 0.019 %

## Conclusion

Delta H_mix is very stable between ecutwfc = 60 Ry and 100 Ry.

For saving compute time, ecutwfc = 60 Ry and ecutrho = 240 Ry appear sufficient for further calculations, based on this x = 0.75 cutoff test.

Next step: k-point convergence test.
