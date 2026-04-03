import pytest
from berakningsfysik.subsystem1.thermodynamics import check_if_valid_x
from berakningsfysik.subsystem1.spinodal_functions import entropy_second_derivative

def test_entropy_second_derivative_valid_input():
    assert entropy_second_derivative(0.5) == pytest.approx(-4.0*8.6173*10**(-5))
    assert entropy_second_derivative(1/3) == pytest.approx(-4.5*8.6173*10**(-5))

def test_entropy_second_derivative_invalid_input():
    """Bör kasta undantag om x är ogiltigt"""
    with pytest.raises(Exception):
        entropy_second_derivative(1.2)
    with pytest.raises(Exception):
        entropy_second_derivative(-0.2)
    with pytest.raises(Exception):
        entropy_second_derivative(0)
    with pytest.raises(Exception):
        entropy_second_derivative(1)