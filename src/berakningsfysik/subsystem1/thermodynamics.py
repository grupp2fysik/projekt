"""Enkla termodynamiska hjälpfunktioner."""

import math

def entropy_at_endpoint(x: float) -> float:
    """Returnerar 0 i ändpunkterna x=0 och x=1."""
    if x == 0 or x == 1:
        return 0.0
    return 1.0

def entropy_per_atom(x, n):
    """Returnerar entropin per metallplats (eV/(atom*K))
    x = komposition
    n = antal metallplatser (2 för TiAlN)"""

    if n < 0 or x < 0 or x > 1:
        raise Exception("Invalid input. Only valid arguments are 0 <= x <= 1 and n > 0.")

    KB = 8.6173*10**(-5) # Boltzmanns konstant (eV/K)

    if entropy_at_endpoint(x) == 1:
        return (-KB*(x*math.log(x) + (1-x)*math.log(1-x))/n)
    else: 
        return entropy_at_endpoint(x)
