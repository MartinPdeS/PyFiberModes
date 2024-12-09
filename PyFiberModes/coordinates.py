import numpy as np
from dataclasses import dataclass
from typing import Tuple


@dataclass
class CylindricalCoordinates:
    """
    Represents a set of points in cylindrical coordinates.

    Attributes
    ----------
    rho : np.ndarray
        Radial distance from the z-axis.
    phi : np.ndarray
        Azimuthal angle in radians.
    z : np.ndarray
        Height along the z-axis.
    """
    rho: np.ndarray
    phi: np.ndarray
    z: np.ndarray

    def to_cartesian(self) -> 'CartesianCoordinates':
        """
        Converts cylindrical coordinates to cartesian coordinates.

        Returns
        -------
        CartesianCoordinates
            Corresponding cartesian coordinates.
        """
        x = self.rho * np.cos(self.phi)
        y = self.rho * np.sin(self.phi)
        return CartesianCoordinates(x=x, y=y, z=self.z)

    def to_cylindrical(self):
        """
        Returns self in cylindrical format (no conversion needed).

        Returns
        -------
        Self
            The current instance in cylindrical format.
        """
        return self


@dataclass
class CartesianCoordinates:
    """
    Represents a set of points in cartesian coordinates.

    Attributes
    ----------
    x : np.ndarray
        x-coordinates (must be a 1D array).
    y : np.ndarray
        y-coordinates (must be a 1D array).
    z : np.ndarray
        z-coordinates (must be a 1D array).
    is_structured : bool
        Indicates if the mesh is structured. Default is False.
    is_3D : bool
        Indicates if the data represents 3D coordinates. Default is False.
    """
    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    is_structured: bool = False
    is_3D: bool = False

    def __post_init__(self) -> None:
        self.x = self.x.astype(float)
        self.y = self.y.astype(float)
        self.z = self.z.astype(float)

    @property
    def x_boundaries(self) -> Tuple[float, float]:
        """
        Calculates the minimum and maximum x-values.

        Returns
        -------
        Tuple[float, float]
            Minimum and maximum x-values.
        """
        return np.min(self.x), np.max(self.x)

    @property
    def y_boundaries(self) -> Tuple[float, float]:
        """
        Calculates the minimum and maximum y-values.

        Returns
        -------
        Tuple[float, float]
            Minimum and maximum y-values.
        """
        return np.min(self.y), np.max(self.y)

    @property
    def dx(self) -> float:
        """
        Calculates the grid spacing in the x-direction.

        Returns
        -------
        float
            Grid spacing in the x-direction.

        Raises
        ------
        ValueError
            If the mesh is unstructured.
        """
        if self.is_structured:
            x_min, x_max = self.x_boundaries
            return (x_max - x_min) / len(self.x)
        raise ValueError("dx value cannot be inferred from an unstructured mesh.")

    @property
    def dy(self) -> float:
        """
        Calculates the grid spacing in the x-direction.

        Returns
        -------
        float
            Grid spacing in the x-direction.

        Raises
        ------
        ValueError
            If the mesh is unstructured.
        """
        if self.is_structured:
            y_min, y_max = self.y_boundaries
            return (y_max - y_min) / len(self.y)
        raise ValueError("dx value cannot be inferred from an unstructured mesh.")

    @classmethod
    def generate_from_boundaries(cls, x_limits, y_limits, z_limits, x_points, y_points, z_points) -> 'CartesianCoordinates':
        """
        Generates a structured cartesian coordinate system from boundaries.

        Parameters
        ----------
        x_limits : list
            Limits of x-axis [min, max].
        y_limits : list
            Limits of y-axis [min, max].
        z_limits : list
            Limits of z-axis [min, max].
        x_points : int
            Number of points along x-axis.
        y_points : int
            Number of points along y-axis.
        z_points : int
            Number of points along z-axis.

        Returns
        -------
        CartesianCoordinates
            Structured cartesian coordinate system.
        """
        x = np.linspace(*x_limits, x_points)
        y = np.linspace(*y_limits, y_points)
        z = np.linspace(*z_limits, z_points)
        return cls(x=x, y=y, z=z, is_structured=True, is_3D=True)

    @classmethod
    def generate_from_square(
        cls,
        length: float,
        center: Tuple[float, float] = (0.0, 0.0),
        n_points: int = 100
    ) -> 'CartesianCoordinates':
        """
        Generates a structured 2D square coordinate system.

        Parameters
        ----------
        length : float
            The length of the square's side.
        center : Tuple[float, float], optional
            The coordinates of the square's center (x0, y0). Default is (0.0, 0.0).
        n_points : int, optional
            The number of points along each dimension (x, y). Default is 100.

        Returns
        -------
        CartesianCoordinates
            Cartesian coordinate system with x, y mesh grids and z set to 0.

        """
        x0, y0 = center

        # Generate x and y ranges centered at `center` with the specified length
        x = np.linspace(-length / 2, length / 2, n_points) + x0
        y = np.linspace(-length / 2, length / 2, n_points) + y0

        # Generate mesh grids
        x_mesh, y_mesh = np.meshgrid(x, y)

        # Create and return a CartesianCoordinates instance
        instance = cls(
            x=x_mesh.T,
            y=y_mesh.T,
            z=np.zeros_like(x_mesh.T),
            is_structured=True,
            is_3D=False
        )

        return instance

    def centering(self) -> None:
        """
        Recenters the coordinates so that their mean is zero along each axis.

        Returns
        -------
        None
            Modifies the CartesianCoordinates instance in place.
        """
        self.x -= self.x.mean()
        self.y -= self.y.mean()
        self.z -= self.z.mean()

    def to_cylindrical(self) -> CylindricalCoordinates:
        """
        Converts cartesian coordinates to cylindrical coordinates.

        Returns
        -------
        CylindricalCoordinates
            Corresponding cylindrical coordinates.
        """
        rho = np.sqrt(self.x**2 + self.y**2)
        phi = np.arctan2(self.y, self.x)
        return CylindricalCoordinates(rho=rho, phi=phi, z=self.z)

    def shift_coordinates(self, shift: Tuple[float, float, float]):
        """
        Shifts the coordinates by specified amounts.

        Parameters
        ----------
        shift : Tuple[float, float, float]
            Amounts to shift along x, y, z.

        Returns
        -------
        Self
            Updated cartesian coordinates.
        """
        self.x += shift[0]
        self.y += shift[1]
        self.z += shift[2]

    def scale_coordinates(self, scale: Tuple[float, float, float]) -> None:
        """
        Scales the coordinates by specified factors.

        Parameters
        ----------
        scale : Tuple[float, float, float]
            Scaling factors for x, y, and z coordinates.

        Returns
        -------
        None
            Modifies the CartesianCoordinates instance in place.
        """
        self.x *= scale[0]
        self.y *= scale[1]
        self.z *= scale[2]
