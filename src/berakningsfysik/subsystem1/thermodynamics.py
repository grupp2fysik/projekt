"""Enkla termodynamiska hjälpfunktioner."""

import math


def check_if_valid_x(x):
    """Kastar undantag om x är större än 1,
    eller mindre än 0."""
    if x < 0 or x > 1:
        raise Exception("Invalid composition x: unpermitted x if x < 0 or x > 1")

def check_if_valid_n(n):
    """n är antalet atomer per metallplats.
    Funktionen kastar undantag om n <= 0."""
    if n <= 0:
        raise Exception("Invalid number of atoms per metal space: unpermitted n if n <= 0")

def entropy_at_endpoint(x: float) -> float:
    """Returnerar 0 i ändpunkterna x=0 och x=1."""
    if x == 0 or x == 1:
        return 0.0
    return 1.0

def entropy_per_atom(x, n):
    """Returnerar entropin per metallplats (eV/(atom*K))
    x = komposition
    n = antal atomer per metallplats (2 för TiAlN)"""

    check_if_valid_n(n)
    check_if_valid_x(x)

    KB = 8.6173*10**(-5) # Boltzmanns konstant (eV/K)

    if entropy_at_endpoint(x) == 1:
        return (-KB*(x*math.log(x) + (1-x)*math.log(1-x))/n)
    else: 
        return entropy_at_endpoint(x)


