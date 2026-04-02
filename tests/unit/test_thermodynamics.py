import pytest
from berakningsfysik.subsystem1.thermodynamics import entropy_at_endpoint
from berakningsfysik.subsystem1.thermodynamics import entropy_per_metalspace

def test_entropy_is_zero_at_x_zero():
    assert entropy_at_endpoint(0) == 0.0


def test_entropy_is_zero_at_x_one():
    assert entropy_at_endpoint(1) == 0.0


def test_entropy_per_metalspace_is_zero_at_x_one():
    assert entropy_per_metalspace(1, 2) == 0.0


def test_entropy_per_metalspace_is_zero_at_x_zero():
    assert entropy_per_metalspace(0, 2) == 0.0


def test_entropy_per_metalspace_invalid_input():
    """Ogiltiga argument bör lyfta undantag"""
    with pytest.raises(Exception):
        entropy_per_metalspace(-1, 2)
    with pytest.raises(Exception):
        entropy_per_metalspace(5, 9)
    with pytest.raises(Exception):
        entropy_per_metalspace(0.5, 0)
    with pytest.raises(Exception):
        entropy_per_metalspace(0.5, -1)


def test_entropy_per_metalspace_valid_input():
    assert entropy_per_metalspace(0.5, 2) == pytest.approx(2.9865286*10**(-5))
    assert entropy_per_metalspace(1/3, 1) == pytest.approx(5.485033542*10**(-5))


