#!/bin/bash

set -e 

alloy_name=$1

#uv run enthalpy_interpolation/interpolation.py --save_model

uv run build_dataframe.py "$alloy_name"

uv run find_phase_curves.py "$alloy_name"

uv run plot_phase_diagram.py "$alloy_name"

uv run plot_entropy.py "$alloy_name"

uv run plot_enthalpy.py "$alloy_name"