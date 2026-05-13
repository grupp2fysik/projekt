"""
Enkla termodynamiska hjälpfunktioner.
"""

import math


def T_string(T):
    string = str(T)
    string = string.replace(".", ",")

    return string


def check_if_valid_T(T):

    try:
        T = float(T)
    except:
        raise ValueError(f"Temperaturer anges som tal. '{T}' är inte en giltig temperatur.")

    if T < 0:
        raise Exception(f"Fel i parameterfil: Alla temperaturer måste vara större än 0.")

def find_temp_limits(T_start: str, T_end:str):

    try:
        T_start = int(T_start)
    except ValueError:
        raise Exception(f"Fel i parameterfil: Gränser på analyserade temperaturområdet måste vara heltal.\
\n{T_start} är inte ett heltal.")
        

    try:
        T_end = int(T_end)
    except ValueError:
        raise Exception(f"Fel i parameterfil: Gränser på analyserade temperaturområdet måste vara heltal.\
\n{T_end} är inte ett heltal.")

        

    check_if_valid_T(T_start)
    check_if_valid_T(T_end)

    if T_end <= T_start:
        raise Exception(f"Fel i parameterfil: Starttemperatur måste vara mindre än sluttemperatur.")
    return T_start, T_end


def check_if_T_in_temps(T:int, temps:list):
    """Kastar undantag om T inte är i temps."""
    print("T", T)
    print("temps", temps)
    if T not in temps:
        raise Exception(f"Den specifika temperaturen {T}K finns inte definierad i din parameterfil.\
\nKan inte göra analys av temperaturer som inte finns definierade där.")


def check_if_valid_x(x):
    """
    Kastar undantag om x är större än 1, eller mindre än 0.
    """

    if x < 0 or x > 1:
        raise Exception("Ogiltig komposition x: otillåten x om x < 0 eller x > 1.")


def check_if_valid_n(n):
<<<<<<< HEAD
    """n är antalet atomer per metallplats.
    Funktionen kastar undantag om n <= 0."""
    try:
        n = int(n)
    except ValueError:
        raise Exception("Fel i parameterfilen: n måste vara heltal.")
=======
    """
    n är antalet atomer per metallplats.
    Funktionen kastar undantag om n <= 0.
    """
>>>>>>> a643f5610ac7ae9ea07bb7511f71a9eadcef5a30

    if n <= 0:
        raise Exception("Ogiltigt antal atomer per metallutrymme: otillåtet n om n <= 0.")


def entropy_at_endpoint(x: float) -> float:
    """
    Returnerar 0 i ändpunkterna x=0 och x=1.
    """

    if x == 0 or x == 1:
        return 0.0
    return 1.0


def entropy_per_atom(x, n):
    """
    Returnerar entropin per metallplats (eV/(atom*K)).
    x = komposition.
    n = antal atomer per metallplats (2 för TiAlN).
    """

    check_if_valid_n(n)
    check_if_valid_x(x)

    KB = 8.6173*10**(-5) # Boltzmanns konstant (eV/K)

    if entropy_at_endpoint(x) == 1:
        return (-KB*(x*math.log(x) + (1-x)*math.log(1-x))/n)
    else: 
        return entropy_at_endpoint(x)