
import argparse
from pathlib import Path

import matplotlib.pyplot as plt

import numpy as np
import pandas as pd


from interpolation import RedlichKisterModel

def main():
    coeffs = np.load("qe_outputs/rk_coeffs.npy")
    print(coeffs)

if __name__ == "__main__":
    main()
