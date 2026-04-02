import pytest
from berakningsfysik.subsystem1.spinodal_functions import entropy_second_derivative

def test_entropy_second_derivative_valid_input():
    assert entropy_second_derivative(0.5) == pytest.approx(-4.0*8.6173*10**(-5))

