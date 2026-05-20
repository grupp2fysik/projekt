TiZrN extra DFT-beräkning.

System:
Ti1-xZrxN

Syfte:
Extra beräkning för kravet att motsvarande DFT-beräkningar ska göras för en annan binär legering.

Punkter:
x = 0.25, 0.50, 0.75

Metod:
- 64-atoms SQS-strukturer
- fixed-cell relaxation
- SCF på relaxerade strukturer
- ecutwfc = 60 Ry
- ecutrho = 240 Ry
- TiZrN k-points = 2 2 2
- ZrN endpoint k-points = 4 4 4

Mappar:
- sqs_strukturer: startstrukturer
- relax_input: QE-input för relaxation
- relaxerade_strukturer: extraherade slutstrukturer
- scf_input: QE-input för SCF
- zrn_endpoint: QE-input för ZrN
- scf_energier: slutenergier från outputfiler
- resultat: beräknad blandningsentalpi
