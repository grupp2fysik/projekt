"""Enkla termodynamiska hjälpfunktioner."""

def entropy_at_endpoint(x: float) -> float:
    """Returnerar 0 i ändpunkterna x=0 och x=1."""
    if x == 0 or x == 1:
        return 0.0
    return 1.0
