# Projekt i beräkningsfysik Grupp 2

Detta repo innehåller grunden för projektets Pythonkod, tester och datahantering.

## Struktur

- `src/` — Pythonkod
- `tests/` — tester med `pytest`
- `docs/` — projektdokumentation
- `data/` — rådata, mellanformat och bearbetad data
- `results/` — genererade resultat

## Pythonmiljö

Projektet använder `uv` för att hantera beroenden och virtuell miljö.

Installera allt med:

> uv sync

Kör tester med

> uv run pytest