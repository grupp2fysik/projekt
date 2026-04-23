from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

RY_TO_EV = 13.60593
BOHR_TO_ANG = 0.529177210903
K_B_EV_PER_K = 8.6173e-5

TOTAL_ENERGY_RE = re.compile(
    r"^\s*!\s+total energy\s*=\s*([-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?)\s+Ry",
    re.MULTILINE,
)

NAT_RE = re.compile(
    r"number of atoms/cell\s*=\s*(\d+)",
    re.IGNORECASE,
)

CELL_BLOCK_RE = re.compile(
    r"CELL_PARAMETERS\s*\(([^)]*)\)\s*\n"
    r"\s*([^\n]+)\n"
    r"\s*([^\n]+)\n"
    r"\s*([^\n]+)",
    re.IGNORECASE,
)

ALAT_IN_HEADER_RE = re.compile(
    r"alat\s*=\s*([-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?)",
    re.IGNORECASE,
)

LATTICE_PARAMETER_RE = re.compile(
    r"lattice parameter \(alat\)\s*=\s*([-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?)\s*a\.u\.",
    re.IGNORECASE,
)

X_FROM_FILENAME_RE = re.compile(
    r"x\s*=\s*([-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?)",
    re.IGNORECASE,
)


@dataclass
class ParsedQEOutput:
    """Hämtad data från en Quantum ESPRESSO .out-fil."""
    path: Path
    x: float
    total_energy_ry: float
    natoms: int
    energy_ev_per_atom: float
    lattice_a_ang: float | float("nan")
    lattice_spread_ang: float | float("nan")


def _parse_vec(line: str) -> np.ndarray:
    """Omvandlar avlästa värden för en cellvektor från .out-fil till en numpy-vektor."""
    vals = [float(v) for v in line.split()]
    if len(vals) != 3:
        raise ValueError(f"Kunde inte tolka cellvektor från rad: {line!r}")
    return np.array(vals, dtype=float)


def _parse_x_from_filename(path: Path) -> float:
    """Läser ut sammansättningen x från filnamnet."""
    match = X_FROM_FILENAME_RE.search(path.name)
    if not match:
        raise ValueError(
            f"Kunde inte hitta sammansättningen x i filnamnet {path.name!r}. "
            f"Förväntar mig t.ex. 'TiAlN_x=0.375.out'."
        )
    x = float(match.group(1))
    if not (0.0 <= x <= 1.0):
        raise ValueError(f"x måste ligga i [0, 1], fick {x} från {path.name!r}")
    return x


def _extract_last_total_energy_ry(text: str) -> float:
    """Hittar sista förekomsten av '!total energy = ... Ry' och returnerar värdet."""
    matches = TOTAL_ENERGY_RE.findall(text)
    if not matches:
        raise ValueError("Ingen konvergerad totalenergi ('! total energy = ... Ry') hittades.")
    return float(matches[-1])


def _extract_natoms(text: str) -> int:
    """Hämtar antalet atomer per cell från .out-filen."""
    match = NAT_RE.search(text)
    if not match:
        raise ValueError("Kunde inte hitta 'number of atoms/cell = ...' i .out-filen.")
    return int(match.group(1))


def _extract_last_lattice_parameter_angstrom(text: str) -> tuple[float, float]:
    """Returnerar genomsnittlig gitterparameter i Å, från den sista CELL_PARAMETERS-blocket."""
    matches = list(CELL_BLOCK_RE.finditer(text))
    if not matches:
        return np.nan, np.nan

    m = matches[-1]
    unit_info = m.group(1).strip().lower()
    a1 = _parse_vec(m.group(2))
    a2 = _parse_vec(m.group(3))
    a3 = _parse_vec(m.group(4))
    cell = np.vstack([a1, a2, a3])

    if "angstrom" in unit_info:
        cell_ang = cell
    elif "bohr" in unit_info:
        cell_ang = cell * BOHR_TO_ANG
    elif "alat" in unit_info:
        alat_match = ALAT_IN_HEADER_RE.search(unit_info)
        if alat_match:
            alat_bohr = float(alat_match.group(1))
        else:
            global_alat = LATTICE_PARAMETER_RE.search(text)
            if not global_alat:
                raise ValueError(
                    "CELL_PARAMETERS är i alat-enheter men inget alat-värde kunde hittas."
                )
            alat_bohr = float(global_alat.group(1))
        cell_ang = cell * alat_bohr * BOHR_TO_ANG
    else:
        raise ValueError(
            f"Okänd enhet i CELL_PARAMETERS: {unit_info!r}. "
            "Stöd finns för angstrom, bohr och alat."
        )

    lengths = np.linalg.norm(cell_ang, axis=1)
    lattice_a_ang = float(np.mean(lengths))
    lattice_spread_ang = float(np.max(lengths) - np.min(lengths))
    return lattice_a_ang, lattice_spread_ang


def parse_qe_out(path: str | Path) -> ParsedQEOutput:
    """Läser en .out_fil och returnerar en ParsedQEOutput med data."""
    path = Path(path)
    text = path.read_text(encoding="utf-8", errors="ignore")

    x = _parse_x_from_filename(path)
    total_energy_ry = _extract_last_total_energy_ry(text)
    natoms = _extract_natoms(text)
    total_energy_ev = total_energy_ry * RY_TO_EV
    energy_ev_per_atom = total_energy_ev / natoms
    lattice_a_ang, lattice_spread_ang = _extract_last_lattice_parameter_angstrom(text)

    return ParsedQEOutput(
        path=path,
        x=x,
        total_energy_ry=total_energy_ry,
        natoms=natoms,
        energy_ev_per_atom=energy_ev_per_atom,
        lattice_a_ang=lattice_a_ang,
        lattice_spread_ang=lattice_spread_ang,
    )


def build_enthalpy_dataframe(
    directory: str | Path,
    glob_pattern: str = "*x=*.out",
) -> pd.DataFrame:
    """
    Skapar en dataframe mha parse_qe() för varje .out-fil i qe_outputs_katalogen.
    Beräknar även blandningsentalpin (eV/atom) med TiN (x=0) och Al(x=1) som referens.
    """
    directory = Path(directory)
    paths = sorted(directory.glob(glob_pattern))
    if not paths:
        raise FileNotFoundError(
            f"Hittade inga filer i {directory} som matchar mönstret {glob_pattern!r}"
        )

    records = [parse_qe_out(p) for p in paths]
    df = pd.DataFrame(
        {
            "file": [r.path.name for r in records],
            "x": [r.x for r in records],
            "total_energy_Ry": [r.total_energy_ry for r in records],
            "natoms": [r.natoms for r in records],
            "energy_eV_per_atom": [r.energy_ev_per_atom for r in records],
            "lattice_A": [r.lattice_a_ang for r in records],
            "lattice_spread_A": [r.lattice_spread_ang for r in records],
        }
    ).sort_values("x", kind="mergesort").reset_index(drop=True)

    tin_row = df.loc[np.isclose(df["x"].to_numpy(), 0.0)]
    aln_row = df.loc[np.isclose(df["x"].to_numpy(), 1.0)]

    if tin_row.empty or aln_row.empty:
        raise ValueError(
            "För att beräkna ΔH_mix behövs referensfiler för både x=0 (TiN) och x=1 (AlN)."
        )

    e_tin = float(tin_row["energy_eV_per_atom"].iloc[0])
    e_aln = float(aln_row["energy_eV_per_atom"].iloc[0])

    x = df["x"].to_numpy()
    e = df["energy_eV_per_atom"].to_numpy()
    hmix = e - ((1.0 - x) * e_tin + x * e_aln)

    df["H_mix_eV_per_atom"] = hmix
    return df


def rk_basis(x: np.ndarray, order: int) -> np.ndarray:
    """
    Bygger designmatrisen A för Redlich-Kister-basen ϕ_i(x) = x(1-x)*(2x-1)^i.
    Returnerar en array med dimension (len(x), order+1).
    """
    x = np.asarray(x, dtype=float)
    z = 2.0 * x - 1.0
    g = x * (1.0 - x)

    A = np.empty((x.size, order + 1), dtype=float)
    for i in range(order + 1):
        A[:, i] = g * z**i
    return A


def rk_basis_d1(x: np.ndarray, order: int) -> np.ndarray:
    """Beräknar blandningsentalpins första koncentrationsderivata."""
    x = np.asarray(x, dtype=float)
    z = 2.0 * x - 1.0
    g = x * (1.0 - x)
    gp = 1.0 - 2.0 * x

    A = np.empty((x.size, order + 1), dtype=float)
    for i in range(order + 1):
        h = z**i
        hp = np.zeros_like(x) if i == 0 else 2.0 * i * z ** (i - 1)
        A[:, i] = gp * h + g * hp
    return A


def rk_basis_d2(x: np.ndarray, order: int) -> np.ndarray:
    """Beräknar blandningsentalpins andra koncentrationsderivata."""
    x = np.asarray(x, dtype=float)
    z = 2.0 * x - 1.0
    g = x * (1.0 - x)
    gp = 1.0 - 2.0 * x
    gpp = -2.0

    # Varje rad i A blir den deriverade serieutvecklingen för rk-utvecklingen
    A = np.empty((x.size, order + 1), dtype=float)
    for i in range(order + 1):
        h = z**i
        hp = np.zeros_like(x) if i == 0 else 2.0 * i * z ** (i - 1)
        hpp = np.zeros_like(x) if i < 2 else 4.0 * i * (i - 1) * z ** (i - 2)
        A[:, i] = gpp * h + 2.0 * gp * hp + g * hpp
    return A


@dataclass
class RedlichKisterModel:
    """
    Redlich-Kister-modell för blandningsentalpin.
    ΔH_mix(x) = x(1-x) * Σ L_i * (2x-1)^i
    """

    coeffs: np.ndarray  # L_i-koefficienter [eV/atom]
    rmse: float         # Root mean square error för anpassningen

    @property
    def order(self) -> int:
        """Polynomets ordning, dvs högsta i."""
        return len(self.coeffs) - 1

    @classmethod
    def fit(cls, x: Iterable[float], hmix: Iterable[float], order: int = 3) -> "RedlichKisterModel":
        """Anpassar modellen till givna (x, ΔH_mix)-data med minsta kvadratmetoden."""
        x = np.asarray(list(x), dtype=float)
        y = np.asarray(list(hmix), dtype=float)

        mask = (x > 0.0) & (x < 1.0)
        x_fit = x[mask]
        y_fit = y[mask]

        if x_fit.size < order + 1:
            raise ValueError(
                f"För få inre datapunkter för ordning {order}. "
                f"Behöver minst {order + 1}, fick {x_fit.size}."
            )

        A = rk_basis(x_fit, order)
        coeffs, *_ = np.linalg.lstsq(A, y_fit, rcond=None)
        residuals = y_fit - A @ coeffs
        rmse = float(np.sqrt(np.mean(residuals**2)))

        return cls(coeffs=coeffs, rmse=rmse)

    def hmix(self, x: Iterable[float]) -> np.ndarray:
        """Beräknar ΔH_mix(x) från modellen."""
        x = np.asarray(list(np.atleast_1d(x)), dtype=float)
        return rk_basis(x, self.order) @ self.coeffs

    def d1(self, x: Iterable[float]) -> np.ndarray:
        """Beräknar första koncentrationsderivatan från modellen av ΔH_mix"""
        x = np.asarray(list(np.atleast_1d(x)), dtype=float)
        return rk_basis_d1(x, self.order) @ self.coeffs

    def d2(self, x: Iterable[float]) -> np.ndarray:
        """Andra koncentrationsderivatan från modellen av ΔH_mix"""
        x = np.asarray(list(np.atleast_1d(x)), dtype=float)
        return rk_basis_d2(x, self.order) @ self.coeffs

def main() -> None:
    """Gränssnitt för terminal. Läser QE-filer, anpassar RK_modell och sparar resultat."""
    parser = argparse.ArgumentParser(description="Interpolera blandningsentalpi från QE .out-filer.")
    parser.add_argument(
        "data_dir",
        nargs="?",
        default="qe_outputs",
        help="Katalog med .out-filer. Standard: qe_outputs",
    )
    parser.add_argument(
        "--glob",
        default="*x=*.out",
        help="Glob-mönster för filnamn. Standard: *x=*.out",
    )
    parser.add_argument(
        "--order",
        type=int,
        default=3,
        help="Ordning på Redlich-Kister-polynomet. Standard: 3",
    )
    parser.add_argument(
        "--save_model",
        action="store_true",
        help="Spara modellens koefficienter till en .npy-fil för senare användning.",
    )
    args = parser.parse_args()

    df = build_enthalpy_dataframe(args.data_dir, glob_pattern=args.glob)
    out_csv = Path(args.data_dir) / "enthalpy_dataset.csv"
    df.to_csv(out_csv, index=False)

    model = RedlichKisterModel.fit(
        x=df["x"].to_numpy(),
        hmix=df["H_mix_eV_per_atom"].to_numpy(),
        order=args.order,
    )

    print(f"Läste {len(df)} filer från: {Path(args.data_dir).resolve()}")
    print(f"Sparade dataset till: {out_csv}")
    print("Redlich-Kister-koefficienter L_i [eV/atom]:")
    for i, Li in enumerate(model.coeffs):
        print(f"L{i} = {Li:.8f}")
    print(f"RMSE = {model.rmse:.6e} eV/atom")

    if args.save_model:
        model_path = Path(args.data_dir) / "rk_coeffs.npy"
        np.save(model_path, model.coeffs)
        print(f"Sparade koefficienter till {model_path}")
    else:
        print("\n* För att spara modellen skriv: [interpolation.py --save_model] *\n")

if __name__ == "__main__":
    main()