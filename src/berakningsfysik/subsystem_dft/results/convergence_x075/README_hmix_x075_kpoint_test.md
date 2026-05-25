# Delta H_mix k-point test for TiAlN x = 0.75

## Purpose

This test was done to study how the mixing enthalpy, Delta H_mix, depends on the k-point grid.

The cutoff was fixed to the value chosen from the cutoff test:

- ecutwfc = 60 Ry
- ecutrho = 240 Ry

## Structures used

TiAlN:
- x = 0.75
- 64-atom SQS structure
- fixed volume
- relaxed structure from input/TiAlN_x075relaxedresult.in

Reference systems:
- cubic TiN, conventional cell, 8 atoms
- cubic AlN, conventional cell, 8 atoms

## K-point grids tested

Low k-point level:

- TiAlN: 2 2 2
- TiN: 4 4 4
- AlN: 4 4 4

High k-point level:

- TiAlN: 4 4 4
- TiN: 8 8 8
- AlN: 8 8 8

This follows the agreed mapping where TiN and AlN use twice as dense k-point grids as TiAlN because their conventional cells are smaller.

## Output locations

QE SCF output files:

- dft-sigma/TiN/output/scf/
- dft-sigma/AlN/output/scf/
- dft-sigma/TiAlN/output/scf/

Energy source summary:

- hmix_energy_sources_x075_kpoint_test.txt

CSV input for H_mix calculation:

- hmix_kpoint_input.csv

Python script:

- calculate_hmix_kpoints.py

Saved H_mix result:

- hmix_results_x075_kpoint_test.txt

## Energy sources

Low k-point level:

- TiN: dft-sigma/TiN/output/scf/TiN_scf_ecut60_k444.out
- AlN: dft-sigma/AlN/output/scf/AlN_scf_ecut60_k444.out
- TiAlN: dft-sigma/TiAlN/output/scf/TiAlN_x075_scf_ecut60_k222.out

High k-point level:

- TiN: dft-sigma/TiN/output/scf/TiN_scf_ecut60_k888.out
- AlN: dft-sigma/AlN/output/scf/AlN_scf_ecut60_k888.out
- TiAlN: dft-sigma/TiAlN/output/scf/TiAlN_x075_scf_ecut60_k444.out

## Result

Delta H_mix values:

- low k-point level:  481.014 meV/formula unit
- high k-point level: 479.389 meV/formula unit

Difference:

- 1.625 meV/formula unit
- about 0.339 %

## Conclusion

The k-point grid has a larger effect on Delta H_mix than the cutoff test did, but the difference is still relatively small.

For x = 0.75:

- cutoff change from 60 Ry to 100 Ry changed Delta H_mix by about 0.019 %
- k-point change from low to high changed Delta H_mix by about 0.339 %

The high k-point result is more accurate, but the low k-point setup is much cheaper. This should be discussed before choosing final production parameters for all x-values.
