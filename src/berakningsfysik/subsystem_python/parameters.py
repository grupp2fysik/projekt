import sys
from help_functions import *

spec_comps = [0.1, 0.3, 0.5, 0.7, 0.9]
alloy_name = sys.argv[1]
initial_columns = ["x", "deltaH", "deltaS"]
num_of_inter_points = 500

# just nu finns dessa i build_dataframe, find_phase_curves_0.py och plot_phase_diagram_0.py

DEFAULT_SYSTEM = "TiAlN"
DEFAULT_RESULTS_DIRNAME = "results"
DEFAULT_PHASE_CURVES_DIRNAME = "phase_curves"
DEFAULT_PHASE_DIAGRAM_DIRNAME = "phase_diagram"
DEFAULT_CURVES_NAME = "curves.csv"
DEFAULT_PHASE_DIAGRAM_NAME = "phase_diagram.png"


plots_dirname = f"plots/{alloy_name}"
n, temps, alloy_name, qe_dir = find_parameters(alloy_name)
