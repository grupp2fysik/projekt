# Selected DFT parameters for TiAlN production calculations

Based on the cutoff and k-point convergence tests for TiAlN x = 0.75, the following parameters were selected for further production calculations.

## Selected parameters

- ecutwfc = 60 Ry
- ecutrho = 240 Ry
- TiAlN k-points: 2 2 2
- TiN/AlN k-points: 4 4 4

## Reasoning

The cutoff test showed that Delta H_mix was stable between 60 Ry and 100 Ry.

- ecutwfc = 60 Ry: 481.014 meV/formula unit
- ecutwfc = 80 Ry: 480.884 meV/formula unit
- ecutwfc = 100 Ry: 480.922 meV/formula unit

The difference between 60 Ry and 100 Ry was about 0.019 %.

The k-point test showed a larger but still moderate difference:

- low k-point level: 481.014 meV/formula unit
- high k-point level: 479.389 meV/formula unit

The difference was about 0.34 %.

Because the project has limited compute time, the cheaper low-k setup was selected so that more x-values can be calculated.

## Current completed production point

The calculation for x = 0.75 is already completed using the selected parameters:

- TiAlN x = 0.75
- ecutwfc = 60 Ry
- ecutrho = 240 Ry
- k-points = 2 2 2

Output:

- dft-sigma/TiAlN/output/scf/TiAlN_x075_scf_ecut60_k222.out
