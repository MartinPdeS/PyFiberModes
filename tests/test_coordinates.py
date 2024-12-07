import pytest
import numpy as np
from PyFiberModes.coordinates import CartesianCoordinates, CylindricalCoordinates  # Replace `coordinates` with your module name


def test_cartesian_to_cylindrical():
    """Test conversion from cartesian to cylindrical coordinates."""
    # Input data
    x = np.array([1, 0, -1])
    y = np.array([0, 1, 0])
    z = np.array([0, 0, 1])
    cartesian = CartesianCoordinates(x=x, y=y, z=z)

    # Conversion
    cylindrical = cartesian.to_cylindrical()

    # Expected results
    rho_expected = np.array([1, 1, 1])
    phi_expected = np.array([0, np.pi / 2, np.pi])

    # Assertions
    np.testing.assert_array_almost_equal(cylindrical.rho, rho_expected)
    np.testing.assert_array_almost_equal(cylindrical.phi, phi_expected)
    np.testing.assert_array_almost_equal(cylindrical.z, z)


def test_cylindrical_to_cartesian():
    """Test conversion from cylindrical to cartesian coordinates."""
    # Input data
    rho = np.array([1, 1, 1])
    phi = np.array([0, np.pi / 2, np.pi])
    z = np.array([0, 0, 1])
    cylindrical = CylindricalCoordinates(rho=rho, phi=phi, z=z)

    # Conversion
    cartesian = cylindrical.to_cartesian()

    # Expected results
    x_expected = np.array([1, 0, -1])
    y_expected = np.array([0, 1, 0])

    # Assertions
    np.testing.assert_array_almost_equal(cartesian.x, x_expected)
    np.testing.assert_array_almost_equal(cartesian.y, y_expected)
    np.testing.assert_array_almost_equal(cartesian.z, z)


def test_generate_from_square():
    """Test generation of structured 2D square coordinates."""
    length = 2.0
    center = (1.0, 1.0)
    n_points = 10

    # Generate square coordinates
    coords = CartesianCoordinates.generate_from_square(length, center, n_points)

    # Assertions
    assert coords.is_structured
    assert not coords.is_3D
    assert coords.x.shape == (n_points, n_points)
    assert coords.y.shape == (n_points, n_points)
    assert np.isclose(coords.x.mean(), center[0])
    assert np.isclose(coords.y.mean(), center[1])


def test_generate_from_boundaries():
    """Test generation of cartesian coordinates from boundaries."""
    x_limits = [-1, 1]
    y_limits = [-1, 1]
    z_limits = [-1, 1]
    x_points = 5
    y_points = 5
    z_points = 5

    # Generate coordinates
    coords = CartesianCoordinates.generate_from_boundaries(
        x_limits=x_limits,
        y_limits=y_limits,
        z_limits=z_limits,
        x_points=x_points,
        y_points=y_points,
        z_points=z_points
    )

    # Assertions
    assert coords.is_structured
    assert coords.is_3D
    assert coords.x.size == x_points
    assert coords.y.size == y_points
    assert coords.z.size == z_points


def test_shift_coordinates():
    """Test shifting cartesian coordinates."""
    x = np.array([1, 2, 3])
    y = np.array([4, 5, 6])
    z = np.array([7, 8, 9])
    shift = (1, -1, 2)

    cartesian = CartesianCoordinates(x=x, y=y, z=z)
    cartesian.shift_coordinates(shift)

    # Assertions
    np.testing.assert_array_equal(cartesian.x, x + shift[0])
    np.testing.assert_array_equal(cartesian.y, y + shift[1])
    np.testing.assert_array_equal(cartesian.z, z + shift[2])


def test_scale_coordinates():
    """Test scaling cartesian coordinates."""
    x = np.array([1, 2, 3])
    y = np.array([4, 5, 6])
    z = np.array([7, 8, 9])
    scale = (2, 0.5, 1.5)

    cartesian = CartesianCoordinates(x=x, y=y, z=z)
    cartesian.scale_coordinates(scale)

    # Assertions
    np.testing.assert_array_equal(cartesian.x, x * scale[0])
    np.testing.assert_array_equal(cartesian.y, y * scale[1])
    np.testing.assert_array_equal(cartesian.z, z * scale[2])


def test_centering():
    """Test centering cartesian coordinates."""
    x = np.array([1, 3, 5])
    y = np.array([2, 4, 6])
    z = np.array([3, 5, 7])

    cartesian = CartesianCoordinates(x=x, y=y, z=z)
    cartesian.centering()

    # Assertions
    np.testing.assert_array_almost_equal(cartesian.x.mean(), 0)
    np.testing.assert_array_almost_equal(cartesian.y.mean(), 0)
    np.testing.assert_array_almost_equal(cartesian.z.mean(), 0)


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
