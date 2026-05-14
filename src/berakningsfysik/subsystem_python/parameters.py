import sys
from help_functions import *

spec_comps = [0.1, 0.3, 0.5, 0.7, 0.9]
alloy_name = sys.argv[1]
initial_columns = ["x", "deltaH", "deltaS"]
num_of_inter_points = 500

n, temps, alloy_name, qe_dir = find_parameters(alloy_name)
