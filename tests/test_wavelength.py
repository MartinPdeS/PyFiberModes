import pytest
from scipy.constants import c, pi
from PyFiberModes.wavelength import Wavelength  # Replace `your_module_name` with the actual module name


def test_wavelength_initialization():
    """Test initialization with a direct wavelength value."""
    wl = Wavelength(500e-9)  # 500 nm
    assert wl == 500e-9
    assert wl.wavelength == 500e-9


def test_wavelength_initialization_k0():
    """Test initialization with wave number (k0)."""
    k0 = 2 * pi / 500e-9  # 500 nm
    wl = Wavelength(k0=k0)
    assert pytest.approx(wl, rel=1e-9) == 500e-9


def test_wavelength_initialization_omega():
    """Test initialization with angular frequency (omega)."""
    omega = c * 2 * pi / 500e-9  # 500 nm
    wl = Wavelength(omega=omega)
    assert pytest.approx(wl, rel=1e-9) == 500e-9


def test_wavelength_initialization_frequency():
    """Test initialization with frequency."""
    frequency = c / 500e-9  # 500 nm
    wl = Wavelength(frequency=frequency)
    assert pytest.approx(wl, rel=1e-9) == 500e-9


def test_wavelength_properties():
    """Test wavelength properties: k0, omega, frequency."""
    wl = Wavelength(500e-9)  # 500 nm

    assert pytest.approx(wl.k0, rel=1e-9) == 2 * pi / 500e-9
    assert pytest.approx(wl.omega, rel=1e-9) == c * 2 * pi / 500e-9
    assert pytest.approx(wl.frequency, rel=1e-9) == c / 500e-9


def test_string_representation():
    """Test string representation of wavelength."""
    wl = Wavelength(500e-9)  # 500 nm
    assert str(wl) == "500.00 nm"


def test_debug_representation():
    """Test debug representation (__repr__)."""
    wl = Wavelength(500e-9)  # 500 nm
    assert repr(wl) == "Wavelength(5.000000e-07 meters)"


def test_zero_wavelength():
    """Test properties with zero wavelength."""
    wl = Wavelength(0.0)

    assert wl.k0 == float('inf')
    assert wl.omega == float('inf')
    assert wl.frequency == float('inf')


def test_invalid_initialization():
    """Test invalid initializations."""
    with pytest.raises(TypeError):
        Wavelength(500e-9, k0=2 * pi / 500e-9)  # Multiple arguments

    with pytest.raises(TypeError):
        Wavelength()  # No arguments

    with pytest.raises(TypeError):
        Wavelength(unknown=1.0)  # Invalid keyword


def test_alias_properties():
    """Test alias properties for frequency and angular frequency."""
    wl = Wavelength(500e-9)  # 500 nm
    assert wl.f == wl.frequency
    assert wl.v == wl.frequency
    assert wl.w == wl.omega


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
