from pathlib import Path

import pandas as pd
import pytest

from berakningsfysik.subsystem2.quantum_espresso_output import (
    parse_quantum_espresso_output,
    parse_total_energy_ry,
    results_to_dataframe,
)


QE_SCF_SAMPLE = Path("tests/data/subsystem2/qe_scf_sample.out")


def test_parse_total_energy_ry():
    energy = parse_total_energy_ry(QE_SCF_SAMPLE)
    assert energy == pytest.approx(-213.08859134)


def test_parse_quantum_espresso_output():
    result = parse_quantum_espresso_output(QE_SCF_SAMPLE)

    assert result["source_file"].endswith("qe_scf_sample.out")
    assert result["total_energy_ry"] == pytest.approx(-213.08859134)


def test_results_to_dataframe():
    results = [parse_quantum_espresso_output(QE_SCF_SAMPLE)]
    dataframe = results_to_dataframe(results)

    assert isinstance(dataframe, pd.DataFrame)
    assert list(dataframe.columns) == ["source_file", "total_energy_ry"]
    assert dataframe.loc[0, "total_energy_ry"] == pytest.approx(-213.08859134)