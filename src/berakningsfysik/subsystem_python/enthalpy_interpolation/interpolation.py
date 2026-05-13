#!/usr/bin/env python3
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

DEFAULT_SYSTEM = "TiAlN"
DEFAULT_MODEL_DIRNAME = "rk_model"
DEFAULT_MODEL_NPZ_NAME = "rk_model.npz"
DEFAULT_COEFFS_NPY_NAME = "rk_coeffs.npy"
DEFAULT_MODEL_SUMMARY_NAME = "rk_model_summary.csv"


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
    """
    Hämtad data från en Quantum ESPRESSO .out-fil.
    """

    path: Path
    x: float
    total_energy_ry: float
    natoms: int
    energy_ev_per_atom: float
    lattice_a_ang: float
    lattice_spread_ang: float


def _parse_vec(line: str) -> np.ndarray:
    """
    Omvandlar avlästa värden för en cellvektor från .out-fil till en numpy-vektor.
    """

    vals = [float(v) for v in line.split()]
    if len(vals) != 3:
        raise ValueError(f"Kunde inte tolka cellvektor från rad: {line!r}")
    return np.array(vals, dtype=float)


def _parse_x_from_filename(path: Path) -> float:
    """
    Läser ut sammansättningen x från filnamnet.
    """

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
    """
    Hittar sista förekomsten av '! total energy = ... Ry' och returnerar värdet.
    """

    matches = TOTAL_ENERGY_RE.findall(text)
    if not matches:
        raise ValueError("Ingen konvergerad totalenergi ('! total energy = ... Ry') hittades.")
    return float(matches[-1])


def _extract_natoms(text: str) -> int:
    """
    Hämtar antalet atomer per cell från .out-filen.
    """

    match = NAT_RE.search(text)
    if not match:
        raise ValueError("Kunde inte hitta 'number of atoms/cell = ...' i .out-filen.")
    return int(match.group(1))


def _extract_last_lattice_parameter_angstrom(text: str) -> tuple[float, float]:
    """
    Returnerar genomsnittlig gitterparameter i Å, från sista CELL_PARAMETERS-blocket.
    """

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
    """
    Läser en .out-fil och returnerar ParsedQEOutput med avläst data.
    """

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


def _normalise_column_name(name: str) -> str:
    """
    Gör kolumnnamn robusta mot mellanslag, bindestreck och versaler.
    """

    return name.strip().lower().replace("-", "_").replace(" ", "_")


def _column_lookup(df: pd.DataFrame) -> dict[str, str]:
    """
    Returnerar mapping från normaliserat kolumnnamn till verkligt kolumnnamn.
    """

    return {_normalise_column_name(c): c for c in df.columns}


def _require_column(df: pd.DataFrame, column: str) -> str:
    """
    Hämtar verkligt kolumnnamn, case-insensitive och robust mot enkla namnvariationer.
    """

    lookup = _column_lookup(df)
    key = _normalise_column_name(column)
    if key not in lookup:
        raise ValueError(
            f"Saknar kolumnen {column!r}. Tillgängliga kolumner är: {list(df.columns)}"
        )
    return lookup[key]


def _hmix_to_ev_per_atom(
    df: pd.DataFrame,
    hmix_column: str,
    atoms_per_formula_unit: float = 2.0,
) -> np.ndarray:
    """
    Läser vald H_mix-kolumn och konverterar den till eV/atom.

    För Ti_{1-x}Al_xN är atoms_per_formula_unit = 2:
    en metallatom + en kväveatom per formelenhet.
    """

    if atoms_per_formula_unit <= 0:
        raise ValueError("atoms_per_formula_unit måste vara positivt.")

    actual_column = _require_column(df, hmix_column)
    key = _normalise_column_name(actual_column)
    values = pd.to_numeric(df[actual_column], errors="raise").to_numpy(dtype=float)

    if key in {"hmix_ev_per_atom", "h_mix_ev_per_atom", "deltah_ev_per_atom"}:
        return values

    if key in {"hmix_mev_per_atom", "h_mix_mev_per_atom", "deltah_mev_per_atom"}:
        return values * 1.0e-3

    if key in {"hmix_ry_per_atom", "h_mix_ry_per_atom", "deltah_ry_per_atom"}:
        return values * RY_TO_EV

    if key in {
        "hmix_ev_per_formula_unit",
        "h_mix_ev_per_formula_unit",
        "deltah_ev_per_formula_unit",
        "hmix_ev_per_fu",
        "h_mix_ev_per_fu",
        "deltah_ev_per_fu",
    }:
        return values / atoms_per_formula_unit

    raise ValueError(
        f"Vet inte hur kolumnen {actual_column!r} ska konverteras. "
        "Stödda H_mix-kolumner är t.ex. "
        "'hmix_ev_per_atom', 'hmix_mev_per_atom', "
        "'hmix_ry_per_atom' och 'hmix_ev_per_formula_unit'."
    )


def _validate_enthalpy_dataframe(df: pd.DataFrame) -> None:
    """
    Grundläggande sanity checks för x och H_mix.
    """

    x = df["x"].to_numpy(dtype=float)
    h = df["H_mix_eV_per_atom"].to_numpy(dtype=float)

    if not np.all(np.isfinite(x)):
        raise ValueError("Alla x-värden måste vara ändliga tal.")

    if not np.all((0.0 <= x) & (x <= 1.0)):
        raise ValueError("Alla x-värden måste ligga i intervallet [0, 1].")

    if not np.all(np.isfinite(h)):
        raise ValueError("Alla H_mix-värden måste vara ändliga tal.")

    duplicate_mask = df["x"].duplicated(keep=False)
    if duplicate_mask.any():
        bad_x = df.loc[duplicate_mask, "x"].to_list()
        raise ValueError(
            f"Hittade dubbla x-värden i dataset: {bad_x}. "
        )


def build_enthalpy_dataframe_from_table(
    table_path: str | Path,
    hmix_column: str = "hmix_ev_per_atom",
    atoms_per_formula_unit: float = 2.0,
) -> pd.DataFrame:
    """
    Läser en CSV-fil med färdigberäknad blandningsentalpi.

    Tabellen måste minst ha kolumnen 'x' och en H_mix-kolumn.
    Funktionen returnerar alltid H_mix både som eV/atom och eV/formelenhet.
    """

    table_path = Path(table_path)
    if not table_path.exists():
        raise FileNotFoundError(f"Hittar inte tabellfilen: {table_path}")

    df = pd.read_csv(table_path, comment="#")
    df.columns = [str(c).strip() for c in df.columns]

    x_col = _require_column(df, "x")

    out = df.copy()
    out["x"] = pd.to_numeric(out[x_col], errors="raise")

    hmix_ev_per_atom = _hmix_to_ev_per_atom(
        out,
        hmix_column=hmix_column,
        atoms_per_formula_unit=atoms_per_formula_unit,
    )

    out["H_mix_eV_per_atom"] = hmix_ev_per_atom
    out["H_mix_eV_per_formula_unit"] = hmix_ev_per_atom * atoms_per_formula_unit

    out = out.sort_values("x", kind="mergesort").reset_index(drop=True)
    _validate_enthalpy_dataframe(out)
    return out


def build_enthalpy_dataframe_from_qe(
    directory: str | Path,
    glob_pattern: str = "*x=*.out",
    atoms_per_formula_unit: float = 2.0,
) -> pd.DataFrame:
    """
    Skapar en dataframe från QE .out-filer.

    Denna gamla väg finns kvar för bakåtkompatibilitet:
    .out-filer -> totalenergi -> H_mix_eV_per_atom.
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
            "För att beräkna ΔH_mix från QE-totalenergier behövs referensfiler "
            "för både x=0 (TiN) och x=1 (AlN)."
        )

    e_tin = float(tin_row["energy_eV_per_atom"].iloc[0])
    e_aln = float(aln_row["energy_eV_per_atom"].iloc[0])

    x = df["x"].to_numpy(dtype=float)
    e = df["energy_eV_per_atom"].to_numpy(dtype=float)
    hmix = e - ((1.0 - x) * e_tin + x * e_aln)

    df["H_mix_eV_per_atom"] = hmix
    df["H_mix_eV_per_formula_unit"] = hmix * atoms_per_formula_unit

    _validate_enthalpy_dataframe(df)
    return df


def build_enthalpy_dataframe(
    data_source: str | Path,
    glob_pattern: str = "*x=*.out",
    hmix_column: str = "hmix_ev_per_atom",
    atoms_per_formula_unit: float = 2.0,
) -> pd.DataFrame:
    """
    Bygger ett entalpidataset från antingen:
      1. en CSV-fil med färdig H_mix, eller
      2. en katalog med QE .out-filer.

    Denna funktion behåller det gamla anropet:
        build_enthalpy_dataframe(directory, glob_pattern="*x=*.out")

    men tillåter nu också:
        build_enthalpy_dataframe("enthalpies.csv", hmix_column="hmix_ev_per_atom")
    """

    data_source = Path(data_source)

    if data_source.is_file():
        return build_enthalpy_dataframe_from_table(
            data_source,
            hmix_column=hmix_column,
            atoms_per_formula_unit=atoms_per_formula_unit,
        )

    if data_source.is_dir():
        return build_enthalpy_dataframe_from_qe(
            data_source,
            glob_pattern=glob_pattern,
            atoms_per_formula_unit=atoms_per_formula_unit,
        )

    raise FileNotFoundError(f"{data_source} är varken en fil eller en katalog.")


def rk_basis(x: np.ndarray, order: int) -> np.ndarray:
    """
    Bygger en matris A med serieutveckling av rk-basen på varje rad för Redlich-Kister-basen.

    Antalet rader = antalet blandningsentalpier (en för varje x).

    Antalet kolonner = 1 + i, där i är av ordningen av Redlich-Kister-basen.

    phi_i(x) = x(1-x)*(2x-1)^i.
    """

    x = np.asarray(x, dtype=float)
    z = 2.0 * x - 1.0
    g = x * (1.0 - x)

    A = np.empty((x.size, order + 1), dtype=float)
    for i in range(order + 1):
        A[:, i] = g * z**i

    return A


def rk_basis_d1(x: np.ndarray, order: int) -> np.ndarray:
    """
    Beräknar första koncentrationsderivatan av Redlich-Kister-basen.
    """

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
    """
    Beräknar andra koncentrationsderivatan av Redlich-Kister-basen.
    """

    x = np.asarray(x, dtype=float)
    z = 2.0 * x - 1.0
    g = x * (1.0 - x)
    gp = 1.0 - 2.0 * x
    gpp = -2.0

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

    Delta H_mix(x) = x(1-x) * sum_i L_i * (2x-1)^i
    """

    coeffs: np.ndarray
    rmse: float

    @property
    def order(self) -> int:
        """
        Polynomets ordning, dvs högsta i.
        """
        
        return len(self.coeffs) - 1

    @classmethod
    def fit(
        cls,
        x: Iterable[float],
        hmix: Iterable[float],
        order: int = 3,
    ) -> "RedlichKisterModel":
        """
        Anpassar modellen till givna (x, Delta H_mix)-data med minsta kvadratmetoden.
        """

        x = np.asarray(list(x), dtype=float)
        y = np.asarray(list(hmix), dtype=float)

        if x.size != y.size:
            raise ValueError(f"x och hmix måste ha samma längd, fick {x.size} och {y.size}.")

        if order < 0:
            raise ValueError("order måste vara >= 0.")

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
        """
        Beräknar Delta H_mix(x) från modellen.
        """

        x = np.asarray(np.atleast_1d(x), dtype=float)
        return rk_basis(x, self.order) @ self.coeffs

    def d1(self, x: Iterable[float]) -> np.ndarray:
        """
        Beräknar första koncentrationsderivatan av Delta H_mix.
        """

        x = np.asarray(np.atleast_1d(x), dtype=float)
        return rk_basis_d1(x, self.order) @ self.coeffs

    def d2(self, x: Iterable[float]) -> np.ndarray:
        """
        Beräknar andra koncentrationsderivatan av Delta H_mix.
        """

        x = np.asarray(np.atleast_1d(x), dtype=float)
        return rk_basis_d2(x, self.order) @ self.coeffs


def package_root() -> Path:
    """
    Returnerar enthalpy_interpolation-mappen.
    """

    return Path(__file__).resolve().parent


def subsystem_root() -> Path:
    """
    Returnerar subsystem_python-mappen.
    """

    return package_root().parent


def default_results_root() -> Path:
    """
    Returnerar subsystem_python/results.
    """

    return subsystem_root() / "results"


def default_system_dir(system: str) -> Path:
    """
    Returnerar t.ex. subsystem_python/results/TiAlN.
    """

    return default_results_root() / system


def default_rk_model_dir(system: str) -> Path:
    """
    Returnerar t.ex. subsystem_python/results/TiAlN/rk_model.
    """

    return default_system_dir(system) / "rk_model"


def save_rk_model(
    model: RedlichKisterModel,
    model_dir: str | Path,
    data_source: str | Path,
    energy_basis: str,
    hmix_fit_col: str,
    hmix_column: str,
    atoms_per_formula_unit: float,
) -> tuple[Path, Path, Path]:
    """
    Sparar RK-modellen strukturerat.

    Skapar:
      - rk_model.npz         : koefficienter + metadata
      - rk_coeffs.npy        : endast koefficienter, för bakåtkompatibilitet
      - rk_model_summary.csv : lättläst sammanfattning
    """

    model_dir = Path(model_dir)
    model_dir.mkdir(parents=True, exist_ok=True)

    npz_path = model_dir / DEFAULT_MODEL_NPZ_NAME
    npy_path = model_dir / DEFAULT_COEFFS_NPY_NAME
    summary_path = model_dir / DEFAULT_MODEL_SUMMARY_NAME

    np.savez(
        npz_path,
        coeffs=model.coeffs,
        rmse=model.rmse,
        order=model.order,
        energy_basis=energy_basis,
        hmix_fit_col=hmix_fit_col,
        hmix_column=hmix_column,
        atoms_per_formula_unit=atoms_per_formula_unit,
        data_source=str(Path(data_source).resolve()),
    )

    # Fortsätt spara koefficientvektor på rk_coeffs.npy för äldre kod som bara förväntar sig koefficienterna.
    np.save(npy_path, model.coeffs)

    summary = {
        "data_source": str(Path(data_source).resolve()),
        "energy_basis": energy_basis,
        "hmix_fit_col": hmix_fit_col,
        "hmix_column": hmix_column,
        "atoms_per_formula_unit": atoms_per_formula_unit,
        "order": model.order,
        "rmse": model.rmse,
    }

    for i, Li in enumerate(model.coeffs):
        summary[f"L{i}"] = Li

    pd.DataFrame([summary]).to_csv(summary_path, index=False)

    return npz_path, npy_path, summary_path


def main() -> None:
    """
    Terminalgränssnitt. Läser CSV eller QE-filer, anpassar RK-modell och sparar dataset.
    """

    parser = argparse.ArgumentParser(
        description="Interpolera blandningsentalpi från CSV-fil eller Quantum ESPRESSO .out-filer."
    )

    parser.add_argument(
        "data_source",
        nargs="?",
        default="qe_outputs/qe_outputs1",
        help="CSV-fil med x och H_mix, eller katalog med QE .out-filer.",
    )

    parser.add_argument(
        "--glob",
        default="*x=*.out",
        help="Glob-mönster om data_source är en QE-katalog. Standard: *x=*.out",
    )

    parser.add_argument(
        "--hmix-column",
        default="hmix_ev_per_atom",
        help=(
            "H_mix-kolumn i CSV-filen. Exempel: hmix_ev_per_atom, "
            "hmix_mev_per_atom, hmix_ry_per_atom eller hmix_ev_per_formula_unit."
        ),
    )

    parser.add_argument(
        "--atoms-per-formula-unit",
        type=float,
        default=2.0,
        help="Antal atomer per formelenhet. För Ti_{1-x}Al_xN är detta 2. Standard: 2.",
    )

    parser.add_argument(
        "--energy-basis",
        choices=["atom", "formula_unit"],
        default="atom",
        help=(
            "Vilken energibas RK-modellen ska anpassas mot. "
            "'atom' ger eV/atom. 'formula_unit' ger eV/formelenhet. Standard: atom."
        ),
    )

    parser.add_argument(
        "--order",
        type=int,
        default=3,
        help="Ordning på Redlich-Kister-polynomet. Standard: 3.",
    )

    parser.add_argument(
        "--save_model",
        action="store_true",
        help="Spara modellens koefficienter till en .npy-fil för senare användning.",
    )

    parser.add_argument(
        "--system",
        default=DEFAULT_SYSTEM,
        help="Materialsystem. Styr standardmapp för modellfiler. Standard: TiAlN.",
    )

    parser.add_argument(
        "--model-dir",
        default=None,
        help=(
            "Mapp där RK-modellen ska sparas. "
            "Standard: <system>/rk_model, t.ex. TiAlN/rk_model."
        ),
    )

    args = parser.parse_args()

    df = build_enthalpy_dataframe(
        args.data_source,
        glob_pattern=args.glob,
        hmix_column=args.hmix_column,
        atoms_per_formula_unit=args.atoms_per_formula_unit,
    )

    hmix_fit_col = {
        "atom": "H_mix_eV_per_atom",
        "formula_unit": "H_mix_eV_per_formula_unit",
    }[args.energy_basis]

    unit_label = {
        "atom": "eV/atom",
        "formula_unit": "eV/formelenhet",
    }[args.energy_basis]

    data_source = Path(args.data_source)
    out_dir = data_source if data_source.is_dir() else data_source.parent

    out_csv = out_dir / "enthalpy_dataset.csv"
    df.to_csv(out_csv, index=False)

    model = RedlichKisterModel.fit(
        x=df["x"].to_numpy(),
        hmix=df[hmix_fit_col].to_numpy(),
        order=args.order,
    )

    print(f"Läste {len(df)} datapunkter från: {data_source.resolve()}")
    print(f"Sparade dataset till: {out_csv}")
    print(f"Anpassningen gjordes mot kolumnen: {hmix_fit_col}")
    print(f"Redlich-Kister-koefficienter L_i [{unit_label}]:")

    for i, Li in enumerate(model.coeffs):
        print(f"L{i} = {Li:.8f}")

    print(f"RMSE = {model.rmse:.6e} {unit_label}")

    model_dir = Path(args.model_dir) if args.model_dir else default_rk_model_dir(args.system)

    if args.save_model:
        npz_path, npy_path, summary_path = save_rk_model(
            model=model,
            model_dir=model_dir,
            data_source=args.data_source,
            energy_basis=args.energy_basis,
            hmix_fit_col=hmix_fit_col,
            hmix_column=args.hmix_column,
            atoms_per_formula_unit=args.atoms_per_formula_unit,
        )

        print(f"Sparade RK-modell till: {npz_path}")
        print(f"Sparade koefficienter till: {npy_path}")
        print(f"Sparade modellsammanfattning till: {summary_path}")
    else:
        print(
            "\n* För att spara modellen skriv: "
            "[interpolation.py --save_model] *\n"
        )


if __name__ == "__main__":
    main()
