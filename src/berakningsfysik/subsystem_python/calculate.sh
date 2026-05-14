#!/bin/bash

set -e

alloy_name=$1

uv run enthalpy_interpolation/interpolation.py "$alloy_name" \
    --hmix-column hmix_ev_per_atom \
    --save_model

uv run build_dataframe.py "$alloy_name"

uv run find_phase_curves_0.py "$alloy_name"

uv run plot_phase_diagram_0.py "$alloy_name"