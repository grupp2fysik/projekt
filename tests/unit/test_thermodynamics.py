from berakningsfysik.subsystem1.thermodynamics import entropy_at_endpoint


def test_entropy_is_zero_at_x_zero():
    assert entropy_at_endpoint(0) == 0.0


def test_entropy_is_zero_at_x_one():
    assert entropy_at_endpoint(1) == 0.0
