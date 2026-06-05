#!/bin/bash

set -euo pipefail

if [ "$#" -ne 4 ]; then
    echo "Usage:"
    echo "  $0 <case_name> <data_csv> <rk_order> <hmix_column>"
    echo ""
    echo "Example:"
    echo "  $0 TiAlN_ours_order5 validation/data/TiAlN_ours.csv 5 hmix_ev_per_formula_unit"
    exit 1
fi

case_name="$1"
data_csv="$2"
rk_order="$3"
hmix_column="$4"

alloy_name="TiAlN"
energy_basis="formula_unit"

case_dir="validation/results/${case_name}"
rk_dir="${case_dir}/rk_model"
thermo_dir="${case_dir}/thermodynamics"
phase_dir="${case_dir}/phase_curves"
gibbs_dir="${phase_dir}/Gibbs_plots"
enthalpy_plot_dir="${case_dir}/interpolation_plot"
phase_diagram_dir="${case_dir}/phase_diagram"
combined_plot_dir="${case_dir}/combined_gibbs_plot"

mkdir -p "$rk_dir"
mkdir -p "$thermo_dir"
mkdir -p "$phase_dir"
mkdir -p "$gibbs_dir"
mkdir -p "$enthalpy_plot_dir"
mkdir -p "$phase_diagram_dir"
mkdir -p "$combined_plot_dir"

echo "=============================================="
echo "Validation case: $case_name"
echo "Data:            $data_csv"
echo "RK order:        $rk_order"
echo "Hmix column:     $hmix_column"
echo "Energy basis:    $energy_basis"
echo "Output:          $case_dir"
echo "=============================================="

python3 enthalpy_interpolation/interpolation.py "$data_csv" \
    --hmix-column "$hmix_column" \
    --order "$rk_order" \
    --energy-basis "$energy_basis" \
    --system "$case_name" \
    --model-dir "$rk_dir" \
    --save_model

python3 build_dataframe.py "$alloy_name" \
    --model-path "$rk_dir" \
    --output "$thermo_dir/dataframe.csv" \
    --energy-basis "$energy_basis"

python3 find_phase_curves.py "$alloy_name" \
    --input "$thermo_dir/dataframe.csv" \
    --output "$phase_dir/curves.csv" \
    --gibbs-plot-dir "$gibbs_dir"

echo "Plotting phase diagram..."

if ! python3 plot_phase_diagram.py "$case_name" \
    --input "$phase_dir/curves.csv" \
    --output "$phase_diagram_dir/phase_diagram.png"; then
    echo "Warning: plot_phase_diagram.py failed for $case_name, but core validation data was generated."
fi

echo "Plotting Redlich-Kister interpolation..."

if ! python3 enthalpy_interpolation/plot_enthalpy.py "$data_csv" \
    --model-path "$rk_dir/rk_model.npz" \
    --hmix-column "$hmix_column" \
    --energy-basis "$energy_basis" \
    --system "$case_name" \
    --plot-dir "$enthalpy_plot_dir" \
    --output "$enthalpy_plot_dir/rk_interpolation_order${rk_order}.png"; then
    echo "Warning: plot_enthalpy.py failed for $case_name."
fi

echo "Plotting combined Hmix and Gibbs curves..."

if ! python3 validation/scripts/plot_hmix_gibbs_fu.py "$data_csv" \
    --model-path "$rk_dir/rk_model.npz" \
    --hmix-column "$hmix_column" \
    --output "$combined_plot_dir/hmix_gibbs_order${rk_order}.png"; then
    echo "Warning: combined Hmix/Gibbs plot failed for $case_name."
fi

echo "Done: $case_name"