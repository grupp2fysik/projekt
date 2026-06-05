# Projekt i beräkningsfysik Grupp 2
===================================

Detta repo innehåller grunden för projektets Pythonkod, tester och datahantering.

För detaljer om hur man använder produkten, se användarhandledningen.

## Fil-Struktur

- src/ - Pythonkod och DFT filer
- tests/ - tester med `pytest`
- docs/ - projektdokumentation
- data/ - rådata, mellanformat och bearbetad data
- results/ - genererade resultat

## Pythonmiljö

Projektet använder `uv` för att hantera beroenden och virtuell miljö.

Installera allt med

> uv sync

## DFT-underlag
===============

Denna mapp innehåller DFT-underlaget för projektets TiAlN-beräkningar samt kompletterande resultat för bandstruktur, DOS och TiZrN.

### Fil-Struktur
subsystem_dft/
    input/
    results/

### input/

`input/` innehåller strukturer, Quantum ESPRESSO-inputfiler och jobbfiler för de slutliga TiAlN-beräkningarna.

Viktiga delar:

`SQSICET.py`
  Script som användes för att skapa SQS-strukturer.

`TiAlN_x*_sqs_64atoms.cif`
  SQS-strukturer för TiAlN i CIF-format.

`TiAlN_x*_sqs_64atoms_qe_structure.in`
  Samma SQS-strukturer i Quantum ESPRESSO-format.

`final_relax_sqs_inputs/`
  Inputfiler för fixed-cell relaxation.

`final_relax_sqs_jobs/`
  Jobbfiler för relaxation på Sigma.

`final_relaxed_sqs_structures/`
  Relaxerade strukturer efter fixed-cell relaxation.

`final_scf_sqs_inputs/`
  Inputfiler för slutliga SCF-beräkningar.

`final_scf_sqs_jobs/`
  Jobbfiler för SCF-beräkningar på Sigma.

`TiAlN_x075relaxedresult.in`
  Relaxerad struktur som användes i konvergenstesterna för x = 0.75.

### results/

`results/` innehåller outputfiler, sammanställda energier, blandningsentalpier, figurer och kompletterande beräkningar.

Viktiga delar:

`final_relax_sqs_outputs/`
  Quantum ESPRESSO-output från relaxationerna.

`final_scf_sqs_outputs/`
  Quantum ESPRESSO-output från SCF-beräkningarna.

`final_hmix_sqs/`
  Slutliga SCF-energier, beräknad blandningsentalpi och figurer för TiAlN.

`convergence_x075/`
  Konvergenstester för cutoff och k-punktsnät vid x = 0.75.

`optimalvolume/`
  Resultat kopplade till optimal volym.

`bandstructure/`
  Underlag och resultat för bandstrukturberäkningar.

`TiAlNdos/`
  Underlag och resultat för DOS-/NSCF-beräkningar.

`tizrn_hmix/`
  Kompletterande DFT-beräkning för TiZrN.

## Köra tester

Kör tester med

> uv run pytest

## Kodkvalitet och testtäckning

> uv run ruff check .

Kör `ruff` för att hitta stilproblem och enkla fel i koden.

> uv run pytest --cov=src

Kör tester med `pytest` och visar code coverage för koden i `src/`.