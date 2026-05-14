import pandas as pd
from thermodynamics import *


def find_parameters(alloy_name):
    """Läser en csv-fil med materialspecifika 
    parametrar och returnerar dessa"""

    df = pd.read_csv(f"alloy_parameters/{alloy_name}.csv", header=None)
    parameters = df[1]
    n = parameters[0]
    T_start = parameters[1]
    T_end = parameters[2]

    T_start, T_end = find_temp_limits(T_start, T_end)

    qe_dir = parameters[3].strip()

    #hittar temperaturer som ska analyseras
    spec_temps = parameters[4]
    spec_temps_list = spec_temps.split()

    temps = [i for i in range(T_start, T_end, math.ceil((T_end - T_start)/30))]

    if T_end not in temps:
        temps.append(T_end)

    print("spec:", spec_temps)
    for index, temp in enumerate(spec_temps_list):

        check_if_valid_T(float(temp))
        if float(temp) not in temps:
            temps.append(float(temp))
    temps.append(20000)
    temps = [x for x in temps if x <= T_end]

    temps.sort()

    check_if_valid_n(n)
    n = int(n)

    print("temps:", temps)
    return n, temps, alloy_name, qe_dir