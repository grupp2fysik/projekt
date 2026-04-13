import pandas as pd
import numpy as np

from enthalpy_interpolation.interpolation import RedlichKisterModel

def main():
    coeffs = np.load("rk_coeffs.npy")
    print(coeffs)

if __name__ == "__main__":
    main()
