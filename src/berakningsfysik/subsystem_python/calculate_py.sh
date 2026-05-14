#!/bin/bash

set -e

alloy_name=$1

python3 enthalpy_interpolation/interpolation.py "$alloy_name" \
    --hmix-column hmix_ev_per_atom \
    --save_model

python3 build_dataframe.py "$alloy_name"

python3 find_phase_curves_0.py "$alloy_name"

python3 plot_phase_diagram_0.py "$alloy_name"

python3 plot_entropy.py "$alloy_name"