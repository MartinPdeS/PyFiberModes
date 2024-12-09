import pytest
import numpy as np
from unittest.mock import MagicMock
from PyFiberModes.field import Field
from PyFiberModes.coordinates import CartesianCoordinates, CylindricalCoordinates
from PyFiberModes.mode import Mode


@pytest.fixture
def mock_fiber():
    """
    Fixture to create a mock fiber object with required methods and attributes.
    """
    fiber = MagicMock()
    fiber.get_radial_field = MagicMock(return_value=(MagicMock(), MagicMock()))
    fiber.wavelength = 1.55e-6
    fiber.get_effective_index = MagicMock(return_value=1.45)
    return fiber


@pytest.fixture
def mock_mode():
    """
    Fixture to create a mock mode object with required attributes.
    """
    mode = MagicMock(spec=Mode)
    mode.nu = 1
    mode.family = "LP"
    return mode


@pytest.fixture
def field_instance(mock_fiber, mock_mode):
    """
    Fixture to create an instance of the Field class.
    """
    return Field(fiber=mock_fiber, mode=mock_mode, limit=10, n_point=101)


def test_field_initialization(field_instance):
    """
    Test initialization of the Field class.
    """
    assert field_instance.limit == 10
    assert field_instance.n_point == 101
    assert isinstance(field_instance.cartesian_coordinates, CartesianCoordinates)
    assert isinstance(field_instance.cylindrical_coordinates, CylindricalCoordinates)


def test_get_azimuthal_dependency(field_instance, mock_mode):
    """
    Test computation of azimuthal dependency.
    """
    dependency_f = field_instance.get_azimuthal_dependency(phi=0, dependency_type='f')
    assert isinstance(dependency_f, np.ndarray)
    assert dependency_f.shape == field_instance.cylindrical_coordinates.phi.shape

    dependency_g = field_instance.get_azimuthal_dependency(phi=0, dependency_type='g')
    assert isinstance(dependency_g, np.ndarray)
    assert dependency_g.shape == field_instance.cylindrical_coordinates.phi.shape

    with pytest.raises(ValueError):
        field_instance.get_azimuthal_dependency(phi=0, dependency_type='invalid')


def test_ex(field_instance):
    """
    Test computation of the electric field's x-component (Ex).
    """
    ex = field_instance.Ex(phi=0, theta=0)
    assert isinstance(ex, np.ndarray)
    assert ex.shape == field_instance.cartesian_coordinates.x.shape


def test_ey(field_instance):
    """
    Test computation of the electric field's y-component (Ey).
    """
    ey = field_instance.Ey(phi=0, theta=0)
    assert isinstance(ey, np.ndarray)
    assert ey.shape == field_instance.cartesian_coordinates.x.shape


def test_ez(field_instance):
    """
    Test computation of the electric field's z-component (Ez).
    """
    ez = field_instance.Ez(phi=0)
    assert isinstance(ez, np.ndarray)
    assert ez.shape == field_instance.cartesian_coordinates.x.shape


def test_get_intensity(field_instance):
    """
    Test the computation of the mode intensity.
    """
    intensity = field_instance.get_intensity()
    assert isinstance(intensity, float)
    assert intensity > 0


def test_get_effective_area(field_instance):
    """
    Test the computation of the effective area of the mode.
    """
    effective_area = field_instance.get_effective_area()
    assert isinstance(effective_area, float)
    assert effective_area > 0


def test_plot(field_instance):
    """
    Test the plotting functionality of the Field class.
    """
    fig = field_instance.plot(plot_type=['Ex', 'Ey'], show=False)
    assert fig is not None


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
