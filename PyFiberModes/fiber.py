#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
import logging
from itertools import count
from scipy.optimize import fixed_point
from dataclasses import dataclass, field
from scipy import constants

from PyFiberModes import fiber_geometry as geometry
from PyFiberModes.fiber_geometry.stepindex import StepIndex
from PyFiberModes import solver
from PyFiberModes import Wavelength, Mode, ModeFamily
from PyFiberModes.functions import get_derivative
from PyFiberModes.field import Field

from MPSTools.fiber_catalogue import loader


@dataclass
class Fiber(object):
    layer_names: list = field(default_factory=[])
    """ Name of each layers """
    layer_radius: list = field(default_factory=[])
    """ Radius of each layers """
    layer_types: list = field(default_factory=[])
    """ Type of each layers """
    material_types: list = field(default_factory=[])
    """ Material of each layers """
    material_parameters: list = field(default_factory=[])
    """ Refractive index of each layers """
    cutoff_solver: float = None
    """ Cutoff frequency class """
    neff_solver: float = None
    """ Neff class """

    logger = logging.getLogger(__name__)

    @property
    def n_layer(self) -> int:
        return len(self.layer_names)

    def __hash__(self):
        return hash(tuple(self.layer_radius))

    def __getitem__(self, index: int) -> object:
        return self.layers[index]

    def add_layer(self, name: str, radius: float, index: float, material_type: str, layer_type: str) -> None:
        self.layer_names.append(name)

        if name != 'cladding':
            self.layer_radius.append(radius)

        self.material_parameters.append(index)
        self.material_types.append(material_type)
        self.layer_types.append(layer_type)

    def initialize_layers(self) -> None:
        """
        Initializes the layers.

        :returns:   No returns
        :rtype:     None
        """
        self.layers = []

        enumerator = zip(
            self.layer_types,
            self.material_types,
            self.material_parameters
        )

        for idx, (layer_type, material_type, material_parameter) in enumerate(enumerator):

            if idx:
                radius_in = self.layer_radius[idx - 1]
            else:
                radius_in = 0

            if idx < len(self.layer_radius):
                radius_out = self.layer_radius[idx]
            else:
                radius_out = numpy.inf

            layer = StepIndex(
                radius_in=radius_in,
                radius_out=radius_out,
                material_type=material_type,
                material_parameter=[material_parameter],
            )

            self.layers.append(layer)

        self.set_solvers(
            cutoff_class=self.cutoff_solver,
            neff_class=self.neff_solver
        )

    def get_layer_name(self, layer_index: int) -> str:
        """
        Gets the layer name.

        :param      layer_index:  The layer index
        :type       layer_index:  int

        :returns:   The layer name.
        :rtype:     str
        """
        return self.layer_names[layer_index]

    def get_layer_at_radius(self, radius: float):
        """
        Gets the layer that is associated to a given radius.

        :param      radius:  The radius
        :type       radius:  float
        """
        radius = abs(radius)
        for idx, layer_radius in enumerate(self.layer_radius):
            if radius < layer_radius:
                return self.layers[idx]

        return self.layers[-1]

    def get_inner_radius(self, layer_idx: int) -> float:
        """
        Gets the radius of the inner most layer

        :param      layer_idx:  The layer index
        :type       layer_idx:  int

        :returns:   The inner radius.
        :rtype:     float
        """
        if layer_idx < 0:
            layer_idx = len(self.layer_radius) + layer_idx + 1

        if layer_idx != 0:
            return self.layer_radius[layer_idx - 1]

        return 0

    def get_outer_radius(self, layer_idx: int) -> float:
        """
        Gets the radius of the outer most layer

        :param      layer_idx:  The layer index
        :type       layer_idx:  int

        :returns:   The inner radius.
        :rtype:     float
        """
        if layer_idx < len(self.layer_radius):
            return self.layer_radius[layer_idx]
        return float("inf")

    def get_thickness(self, layer_idx: int) -> float:
        """
        Gets the thickness of a specific layer.

        :param      layer_idx:  The layer index
        :type       layer_idx:  int

        :returns:   The thickness.
        :rtype:     float
        """
        outer_radius = self.get_outer_radius(layer_idx=layer_idx)
        inner_radius = self.get_inner_radius(layer_idx=layer_idx)
        return outer_radius - inner_radius

    def get_fiber_radius(self) -> float:
        """
        Gets the fiber total radius taking account for all layers.

        :returns:   The fiber radius.
        :rtype:     float
        """
        layer_radius = [
            layer.radius_out for layer in self.layers
        ]

        largest_radius = numpy.max(layer_radius)

        return largest_radius

    def get_index_at_radius(self, radius: float, wavelength: float) -> float:
        """
        Gets the refractive index at a given radius.

        :param      radius:      The radius
        :type       radius:      float
        :param      wavelength:  The wavelength
        :type       wavelength:  float

        :returns:   The refractive index at given radius.
        :rtype:     float
        """
        layer = self.get_layer_at_radius(radius)

        return layer.index(radius, wavelength)

    def get_minimum_index(self, layer_idx: int, wavelength: float) -> float:
        """
        Gets the minimum refractive index of the layers.

        :param      layer_idx:   The layer index
        :type       layer_idx:   int
        :param      wavelength:  The wavelength
        :type       wavelength:  float

        :returns:   The minimum index.
        :rtype:     float
        """
        layer = self.layers[layer_idx]

        return layer.get_minimum_index(wavelength)

    def get_maximum_index(self, layer_idx: int, wavelength: float) -> float:
        """
        Gets the maximum refractive index of the layers.

        :param      layer_idx:   The layer index
        :type       layer_idx:   int
        :param      wavelength:  The wavelength
        :type       wavelength:  float

        :returns:   The minimum index.
        :rtype:     float
        """
        layer = self.layers[layer_idx]

        return layer.get_maximum_index(wavelength=wavelength)

    def find_cutoff_solver(self) -> solver.solver.FiberSolver:
        n_layers = len(self.layers)

        match n_layers:
            case 2:  # Standard Step-Index Fiber [SSIF|
                return solver.ssif.CutoffSolver
            case 3:  # Three-Layer Step-Index Fiber [TLSIF]
                return solver.tlsif.CutoffSolver
            case _:  # Multi-Layer Step-Index Fiber [MLSIF]
                return solver.solver.FiberSolver

    def find_neff_solver(self) -> solver.solver.FiberSolver:
        neff_class = solver.mlsif.NeffSolver

        number_of_layers = len(self.layers)

        if all(isinstance(layer, geometry.StepIndex) for layer in self.layers):
            if number_of_layers == 2:  # SSIF
                neff_class = solver.ssif.NeffSolver

        return neff_class

    def set_solvers(self, cutoff_class=None, neff_class=None) -> None:
        assert cutoff_class is None or issubclass(cutoff_class, solver.FiberSolver)

        assert neff_class is None or issubclass(neff_class, solver.FiberSolver)

        if cutoff_class is None:
            cutoff_class = self.find_cutoff_solver()

        if neff_class is None:
            neff_class = self.find_neff_solver()

        self.neff_solver = neff_class(self)
        self.cutoff_solver = cutoff_class(self)

    def get_NA(self, wavelength: float) -> float:
        """
        Gets the numerical aperture NA.

        :param      wavelength:  The wavelength
        :type       wavelength:  float

        :returns:   The numerical aperture.
        :rtype:     float
        """
        layers_maximum_index = [
            layer.get_maximum_index(wavelength=wavelength) for layer in self.layers
        ]

        n1 = max(layers_maximum_index)

        n2 = self.get_minimum_index(layer_idx=-1, wavelength=wavelength)

        return numpy.sqrt(n1**2 - n2**2)

    def get_V0(self, wavelength: float) -> float:
        """
        Gets the V0 parameter.

        :param      wavelength:  The wavelength
        :type       wavelength:  float

        :returns:   The parameter V0.
        :rtype:     float
        """
        wavelength = Wavelength(wavelength)

        NA = self.get_NA(wavelength=wavelength)

        inner_radius = self.get_inner_radius(layer_idx=-1)

        V0 = wavelength.k0 * inner_radius * NA

        return V0

    def V0_to_wavelength(self, V0: float, maxiter: int = 500, tol: float = 1e-15) -> float:
        """
        Convert V0 number to wavelength.
        An iterative method is used, since the index can be wavelength dependant.

        :param      V0:       The V0
        :type       V0:       float
        :param      maxiter:  The maxiter
        :type       maxiter:  int
        :param      tol:      The tolerance
        :type       tol:      float

        :returns:   The associated wavelength
        :rtype:     float
        """
        if V0 == 0:
            return float("inf")
        if numpy.isinf(V0):
            return 0

        def model(wl):
            return 2 * numpy.pi / V0 * b * self.get_NA(wavelength=wl)

        b = self.get_inner_radius(-1)

        wavelength = model(wl=1.55e-6)

        if abs(wavelength - model(wl=wavelength)) > tol:
            for w in (1.55e-6, 5e-6, 10e-6):
                try:
                    wavelength = fixed_point(model, w, xtol=tol, maxiter=maxiter)
                except RuntimeError:
                    # FIXME: What should we do if it does not converge?
                    self.logger.info(
                        f"V0_to_wavelength: did not converged from {w * 1e6}µm for {V0=} ({wavelength=})"
                    )
                if wavelength > 0:
                    break

        if wavelength == 0:
            self.logger.error(
                f"V0_to_wavelength: did not converged for {V0=} {wavelength=})"
            )

        return Wavelength(wavelength)

    def get_cutoff(self, mode: Mode) -> float:
        """
        Gets the cutoff wavelength of the fiber.

        :param      mode:  The mode
        :type       mode:  Mode

        :returns:   The cutoff wavelength.
        :rtype:     float
        """
        if mode in [Mode("HE", 1, 1), Mode("LP", 0, 1)]:
            return 0

        return self.cutoff_solver(mode=mode)

    def get_cutoff_wavelength(self, mode: Mode) -> float:
        """
        Gets the cutoff wavelength.

        :param      mode:  The mode
        :type       mode:  Mode

        :returns:   The cutoff wavelength.
        :rtype:     float
        """
        cutoff = self.get_cutoff(mode=mode)
        wavelength = self.V0_to_wavelength(V0=cutoff)
        return wavelength

    def get_effective_index(self,
            mode: Mode,
            wavelength: float,
            delta: float = 1e-6,
            lowbound: float = None) -> float:
        """
        Gets the effective index.

        :param      mode:        The mode
        :type       mode:        Mode
        :param      wavelength:  The wavelength
        :type       wavelength:  float
        :param      delta:       The delta
        :type       delta:       float
        :param      lowbound:    The lowbound
        :type       lowbound:    float

        :returns:   The effective index.
        :rtype:     float
        """
        wavelength = Wavelength(wavelength)

        neff = self.neff_solver(
            wavelength=wavelength,
            mode=mode,
            delta=delta,
            lowbound=lowbound
        )

        return neff

    def get_derivative_beta_vs_omega(self,
            omega: float,
            mode: Mode,
            order: float = 0,
            delta: float = 1e-6,
            lowbound: float = None) -> float:
        """
        Gets the derivative of beta vs omega.

        :param      omega:     The pulsation omega
        :type       omega:     float
        :param      mode:      The mode
        :type       mode:      Mode
        :param      order:     The order of the derivative
        :type       order:     float
        :param      delta:     The delta
        :type       delta:     float
        :param      lowbound:  The lowbound
        :type       lowbound:  float

        :returns:   The derivative beta vs omega.
        :rtype:     float
        """

        wavelength = Wavelength(omega=omega)

        if order == 0:
            neff = self.get_effective_index(
                mode=mode,
                wavelength=wavelength,
                delta=delta,
                lowbound=lowbound
            )

            return neff * wavelength.k0

        n_point = 5
        central_point = (n_point - 1) // 2
        delta = 1e12  # This value is critical for accurate computation
        lb = lowbound

        for i in range(n_point - 1, -1, -1):  # Precompute neff using previous wavelength
            o = omega + (i - central_point) * delta

            wavelength = Wavelength(omega=o)

            lb = self.get_effective_index(
                mode=mode,
                wavelength=wavelength,
                delta=delta,
                lowbound=lb
            )

            lb += delta * 1.1

        derivative = get_derivative(
            function=self.get_derivative_beta_vs_omega,
            x=omega,
            order=order,
            n_point=n_point,
            central_point=central_point,
            delta=delta,
            function_args=(mode, 0, delta, lowbound)
        )

        return derivative

    def get_normalized_beta(self,
            mode: Mode,
            wavelength: float,
            delta: float = 1e-6,
            lowbound: float = None) -> float:
        """
        Gets the normalized propagation constant [beta].

        :param      mode:        The mode
        :type       mode:        Mode
        :param      wavelength:  The wavelength
        :type       wavelength:  float
        :param      delta:       The delta
        :type       delta:       float
        :param      lowbound:    The lowbound
        :type       lowbound:    float

        :returns:   The normalized propagation constant.
        :rtype:     float
        """
        neff = self.get_effective_index(
            mode=mode,
            wavelength=wavelength,
            delta=delta,
            lowbound=lowbound
        )

        layer_maximum_index = [
            layer.get_maximum_index(wavelength) for layer in self.layers
        ]

        overall_maximum_index = max(layer_maximum_index)

        minimum_index = self.get_minimum_index(-1, wavelength)

        numerator = neff**2 - minimum_index**2

        denominator = overall_maximum_index**2 - minimum_index**2

        return numerator / denominator

    def get_phase_velocity(self,
            mode: Mode,
            wavelength: float,
            delta: float = 1e-6,
            lowbound=None) -> float:
        """
        Gets the phase velocity.

        :param      mode:        The mode
        :type       mode:        Mode
        :param      wavelength:  The wavelength
        :type       wavelength:  float
        :param      delta:       The delta
        :type       delta:       float
        :param      lowbound:    The lowbound
        :type       lowbound:    float

        :returns:   The phase velocity.
        :rtype:     float
        """
        n_eff = self.get_effective_index(
            mode=mode,
            wavelength=wavelength,
            delta=delta,
            lowbound=lowbound
        )

        return constants.c / n_eff

    def get_group_index(self,
            mode: Mode,
            wavelength: float,
            delta: float = 1e-6,
            lowbound: float = None) -> float:
        """
        Gets the group index.

        :param      mode:        The mode
        :type       mode:        Mode
        :param      wavelength:  The wavelength
        :type       wavelength:  float
        :param      delta:       The delta
        :type       delta:       float
        :param      lowbound:    The lowbound
        :type       lowbound:    float

        :returns:   The group index.
        :rtype:     float
        """
        wavelength = Wavelength(wavelength)

        beta = self.get_derivative_beta_vs_omega(
            omega=wavelength.omega,
            mode=mode,
            order=1,
            delta=delta,
            lowbound=lowbound
        )

        return beta * constants.c

    def get_groupe_velocity(self,
            mode: Mode,
            wavelength: float,
            delta: float = 1e-6,
            lowbound: float = None) -> float:
        """
        Gets the groupe velocity.

        :param      mode:        The mode
        :type       mode:        Mode
        :param      wavelength:  The wavelength
        :type       wavelength:  float
        :param      delta:       The delta
        :type       delta:       float
        :param      lowbound:    The lowbound
        :type       lowbound:    float

        :returns:   The groupe velocity.
        :rtype:     float
        """
        wavelength = Wavelength(wavelength)
        beta = self.get_derivative_beta_vs_omega(
            omega=wavelength.omega,
            mode=mode,
            order=1,
            delta=delta,
            lowbound=lowbound
        )

        return 1 / beta

    def get_dispersion(self,
            mode: Mode,
            wavelength: float,
            delta: float = 1e-6,
            lowbound: float = None) -> float:
        """
        Gets the modal dispersion.

        :param      mode:        The mode
        :type       mode:        Mode
        :param      wavelength:  The wavelength
        :type       wavelength:  float
        :param      delta:       The delta
        :type       delta:       float
        :param      lowbound:    The lowbound
        :type       lowbound:    float

        :returns:   The modal dispersion.
        :rtype:     float
        """
        wavelength = Wavelength(wavelength)
        beta = self.get_derivative_beta_vs_omega(
            omega=wavelength.omega,
            mode=mode,
            order=2,
            delta=delta,
            lowbound=lowbound
        )

        return -beta * 2 * numpy.pi * constants.c * 1e6 / wavelength**2

    def get_S_parameter(self,
            mode: Mode,
            wavelength: float,
            delta: float = 1e-6,
            lowbound: float = None) -> float:
        """
        Gets the s parameter.

        :param      mode:        The mode
        :type       mode:        Mode
        :param      wavelength:  The wavelength
        :type       wavelength:  float
        :param      delta:       The delta
        :type       delta:       float
        :param      lowbound:    The lowbound
        :type       lowbound:    float

        :returns:   The s parameter.
        :rtype:     float
        """
        wavelength = Wavelength(wavelength)

        beta = self.get_derivative_beta_vs_omega(
            omega=wavelength.omega,
            mode=mode,
            order=3,
            delta=delta,
            lowbound=lowbound
        )

        return 1e-3 * beta * (2 * numpy.pi * constants.c / wavelength**2)**2

    def get_vectorial_modes(self,
            wavelength: float,
            nu_max=None,
            m_max=None,
            delta: float = 1e-6) -> set:
        """
        Gets the family of vectorial modes.

        :param      wavelength:  The wavelength
        :type       wavelength:  float
        :param      nu_max:      The numax
        :type       nu_max:      int
        :param      m_max:       The mmax
        :type       m_max:       int
        :param      delta:       The delta
        :type       delta:       float

        :returns:   The vectorial modes.
        :rtype:     set
        """
        families = (
            ModeFamily.HE,
            ModeFamily.EH,
            ModeFamily.TE,
            ModeFamily.TM
        )

        modes = self.get_modes_from_familly(
            families=families,
            wavelength=wavelength,
            nu_max=nu_max,
            m_max=m_max,
            delta=delta
        )

        return modes

    def get_LP_modes(self,
            wavelength: float,
            ellmax: int = None,
            m_max: int = None,
            delta: float = 1e-6) -> set:
        """
        Gets the family of LP modes.

        :param      wavelength:  The wavelength
        :type       wavelength:  float
        :param      ellmax:      The ellmax
        :type       ellmax:      int
        :param      mmax:        The mmax
        :type       mmax:        int
        :param      delta:       The delta
        :type       delta:       float

        :returns:   The lp modes.
        :rtype:     set
        """
        families = (ModeFamily.LP,)

        modes = self.get_modes_from_familly(
            families=families,
            wavelength=wavelength,
            nu_max=ellmax,
            m_max=m_max,
            delta=delta
        )

        return modes

    def get_modes_from_familly(self,
            families,
            wavelength: float,
            nu_max: int = numpy.inf,
            m_max: int = numpy.inf,
            delta: float = 1e-6) -> set:
        """
        Find all modes of given families, within given constraints

        :param      families:    The families
        :type       families:    object
        :param      wavelength:  The wavelength
        :type       wavelength:  float
        :param      numax:       The radial number nu maximum to reach
        :type       numax:       int
        :param      mmax:        The azimuthal number m maximum to reach
        :type       mmax:        int
        :param      delta:       The delta
        :type       delta:       float

        :returns:   The modes from familly.
        :rtype:     set
        """
        modes = set()
        v0 = self.get_V0(wavelength=wavelength)

        for family in families:
            for nu in count(0):

                try:
                    _mmax = m_max[nu]
                except IndexError:
                    _mmax = m_max[-1]
                except TypeError:
                    _mmax = m_max

                if family in [ModeFamily.TE, ModeFamily.TM] and nu > 0:
                    break

                if family in [ModeFamily.HE, ModeFamily.EH] and nu == 0:
                    continue

                if nu > nu_max:
                    break

                for m in count(1):
                    if m > _mmax:
                        break

                    mode = Mode(family, nu, m)

                    try:
                        if self.get_cutoff(mode=mode) > v0:
                            break

                    except (NotImplementedError, ValueError):
                        neff = self.get_effective_index(
                            mode=mode,
                            wavelength=wavelength,
                            delta=delta
                        )

                        if numpy.isnan(neff):
                            break

                    modes.add(mode)

                if m == 1:
                    break
        return modes

    def get_field(self,
            mode: Mode,
            wavelength: float,
            limit: float = None,
            n_point: int = 101) -> Field:
        """
        Get field class

        :param      mode:        The mode
        :type       mode:        Mode
        :param      wavelength:  The wavelength
        :type       wavelength:  float
        :param      limit:       The limit boundary
        :type       limit:       float
        :param      n_point:     The number of point for axis discreditization
        :type       n_point:     int

        :returns:   The field instance of the mode.
        :rtype:     Field
        """
        if limit is None:
            limit = self.get_fiber_radius()

        field = Field(
            fiber=self,
            mode=mode,
            wavelength=wavelength,
            limit=limit,
            n_point=n_point
        )

        return field

    def get_radial_field(self,
            mode: Mode,
            wavelength: float,
            radius: float) -> float:
        """
        Gets the mode field without the azimuthal component.

        :param      mode:        The mode
        :type       mode:        Mode
        :param      wavelength:  The wavelength
        :type       wavelength:  float
        :param      radius:      The radius
        :type       radius:      float

        :returns:   The radial field.
        :rtype:     float
        """
        neff = self.get_effective_index(
            mode=mode,
            wavelength=wavelength
        )

        kwargs = dict(
            wavelength=wavelength,
            nu=mode.nu,
            neff=neff,
            radius=radius
        )

        match mode.family:
            case ModeFamily.LP:
                return self.neff_solver.get_LP_field(**kwargs)
            case ModeFamily.TE:
                return self.neff_solver.get_TE_field(**kwargs)
            case ModeFamily.TM:
                return self.neff_solver.get_TM_field(**kwargs)
            case ModeFamily.EH:
                return self.neff_solver.get_EH_field(**kwargs)
            case ModeFamily.HE:
                return self.neff_solver.get_HE_field(**kwargs)


def load_fiber(fiber_name: str, wavelength: float = None):
    """
    Loads a fiber as type that suit PyFiberModes.

    :param      fiber_name:  The fiber name
    :type       fiber_name:  str
    :param      wavelength:  The wavelength
    :type       wavelength:  float

    :returns:   { description_of_the_return_value }
    :rtype:     { return_type_description }
    """

    fiber_dict = loader.load_fiber_as_dict(
        fiber_name=fiber_name,
        wavelength=wavelength,
        order='out-to-in'
    )

    index_list = []
    name_list = []
    radius_list = []
    material_types_list = []
    layer_type_list = []

    for _, layer in fiber_dict['layers'].items():
        if layer.get('name') in ['air']:
            continue
        index_list.append([layer.get('index')])
        name_list.append(layer.get('name'))
        radius_list.append(layer.get('radius'))
        material_types_list.append('Fixed')
        layer_type_list.append('StepIndex')

    fiber = Fiber(
        layer_radius=radius_list,
        layer_types=layer_type_list,
        material_types=material_types_list,
        material_parameters=index_list,
        layer_names=name_list,
    )

    return fiber

# -
