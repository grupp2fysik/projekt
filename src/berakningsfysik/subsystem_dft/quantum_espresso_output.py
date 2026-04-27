from pathlib import Path
import re

import pandas as pd


TOTAL_ENERGY_PATTERN = re.compile(
    r"^\!\s+total energy\s+=\s+([-+]?\d+\.\d+)\s+Ry",
    re.MULTILINE,
)


def parse_total_energy_ry(file_path: str | Path) -> float:
    text = Path(file_path).read_text(encoding="utf-8")
    matches = TOTAL_ENERGY_PATTERN.findall(text)

    if not matches:
        raise ValueError(f"Could not find total energy in file: {file_path}")

    return float(matches[-1])


def parse_quantum_espresso_output(file_path: str | Path) -> dict[str, float | str]:
    return {
        "source_file": str(file_path),
        "total_energy_ry": parse_total_energy_ry(file_path),
    }


def results_to_dataframe(results: list[dict[str, float | str]]) -> pd.DataFrame:
    return pd.DataFrame(results)