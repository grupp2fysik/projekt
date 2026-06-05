#!/bin/bash

set -euo pipefail

./validation/scripts/run_validation_case.sh \
    TiAlN_ours_order2 \
    validation/data/TiAlN_ours.csv \
    2 \
    hmix_ev_per_formula_unit

./validation/scripts/run_validation_case.sh \
    TiAlN_ours_order3 \
    validation/data/TiAlN_ours.csv \
    3 \
    hmix_ev_per_formula_unit

./validation/scripts/run_validation_case.sh \
    TiAlN_ours_order5 \
    validation/data/TiAlN_ours.csv \
    5 \
    hmix_ev_per_formula_unit

./validation/scripts/run_validation_case.sh \
    alling_order2 \
    validation/data/alling_random_alloy_fu.csv \
    2 \
    hmix_ev_per_formula_unit

./validation/scripts/run_validation_case.sh \
    alling_order3 \
    validation/data/alling_random_alloy_fu.csv \
    3 \
    hmix_ev_per_formula_unit

./validation/scripts/run_validation_case.sh \
    alling_order5 \
    validation/data/alling_random_alloy_fu.csv \
    5 \
    hmix_ev_per_formula_unit