import pytest
import numpy as np
from types import SimpleNamespace
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


def test_radial_field_solution_is_cached(field_instance, mock_fiber):
    field_instance.Ex()
    first_call_count = mock_fiber.get_radial_field.call_count
    field_instance.Ex()
    assert first_call_count == max(257, 2 * field_instance.n_point)
    assert mock_fiber.get_radial_field.call_count == first_call_count


def test_all_radial_components_share_one_solver_pass(field_instance, mock_fiber):
    electric = SimpleNamespace(rho=1, phi=2, z=3)
    magnetic = SimpleNamespace(rho=4, phi=5, z=6)
    mock_fiber.get_radial_field.return_value = electric, magnetic

    fields = field_instance.get_components()
    expected_calls = max(257, 2 * field_instance.n_point)

    assert set(fields) == {"Ex", "Ey", "Ez", "Hx", "Hy", "Hz"}
    assert mock_fiber.get_radial_field.call_count == expected_calls
    assert field_instance.get_components()["Ex"] is fields["Ex"]
    assert mock_fiber.get_radial_field.call_count == expected_calls


def test_field_cache_is_invalidated_by_wavelength_change(field_instance, mock_fiber):
    electric = SimpleNamespace(rho=1, phi=2, z=3)
    magnetic = SimpleNamespace(rho=4, phi=5, z=6)
    mock_fiber.get_radial_field.return_value = electric, magnetic
    expected_calls = max(257, 2 * field_instance.n_point)

    field_instance.get_components(("Ex",))
    mock_fiber.wavelength = 1.31e-6
    field_instance.get_components(("Ex",))

    assert mock_fiber.get_radial_field.call_count == 2 * expected_calls


def test_get_components_rejects_unknown_names(field_instance):
    with pytest.raises(ValueError, match="Unknown field component"):
        field_instance.get_components(("not_a_component",))


def test_poynting_power_and_confinement(field_instance, monkeypatch):
    shape = field_instance.cartesian_coordinates.x.shape
    zero = np.zeros(shape)
    one = np.ones(shape)
    monkeypatch.setattr(field_instance, "Ex", lambda *args: one)
    monkeypatch.setattr(field_instance, "Ey", lambda *args: zero)
    monkeypatch.setattr(field_instance, "Ez", lambda *args: zero)
    monkeypatch.setattr(field_instance, "Hx", lambda *args: zero)
    monkeypatch.setattr(field_instance, "Hy", lambda *args: one)
    monkeypatch.setattr(field_instance, "Hz", lambda *args: zero)
    monkeypatch.setattr(field_instance, "Emod", lambda *args: one)
    sx, sy, sz = field_instance.get_poynting_vector()
    assert np.allclose(sx, 0)
    assert np.allclose(sy, 0)
    assert np.allclose(sz, 0.5)
    assert field_instance.get_power() > 0
    assert 0 < field_instance.get_confinement_factor(field_instance.limit / 2) < 1


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

    single_component_figure = field_instance.plot(plot_type=["Ex"], show=False)
    assert single_component_figure is not None


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
