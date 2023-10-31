from PyFiberModes import fiber_geometry as geometry
from PyFiberModes import solver
from math import sqrt, isnan, isinf, pi
from PyFiberModes import Wavelength, Mode, ModeFamily
from PyFiberModes import constants
from PyFiberModes.functions import derivative
from PyFiberModes.field import Field
from itertools import count
import logging
from scipy.optimize import fixed_point
from functools import lru_cache
from dataclasses import dataclass


# @dataclass
class Fiber(object):
    # layer_radius: list
    # layer_type: list
    # fp: object
    # material_parameter: list
    # layer_index: list
    # names: list
    # Cutoff = None
    # Neff = None

    logger = logging.getLogger(__name__)

    def __init__(self,
            layer_radius: list,
            layer_types: list,
            fp,
            material_parameters: list,
            layer_indexes,
            layer_names: list,
            Cutoff=None,
            Neff=None):

        self.layer_radius = layer_radius
        self.layer_names = layer_names

        self.layers = []

        enumerator = zip(layer_types, fp, material_parameters, layer_indexes)

        for idx, (layer_type, fp_, m_, mp_) in enumerate(enumerator):
            ri = self.layer_radius[idx - 1] if idx else 0

            ro = self.layer_radius[idx] if idx < len(layer_radius) else float("inf")

            layer = geometry.__dict__[layer_type](
                ri,
                ro,
                *fp_,
                m=m_,
                mp=mp_,
                cm=material_parameters[-1],
                cmp=layer_indexes[-1]
            )

            self.layers.append(layer)

        self.cutoff_cache = {
            Mode("HE", 1, 1): 0,
            Mode("LP", 0, 1): 0
        }

        self.ne_cache = {}

        self.set_solvers(Cutoff, Neff)

    def __len__(self):
        return len(self.layers)

    def __str__(self):
        s = "Fiber {\n"
        for i, layer in enumerate(self.layers):
            geom = str(layer)
            radius = self.get_outer_radius(i)
            radius = '' if isinf(radius) else f' {radius * 1e6} µm'
            name = self.get_layer_name(i)
            name = f' "{name}"'
            s += f"    {geom}{radius}{name}\n"
        s += "}"
        return s

    def fixedMatFiber(self, wavelength: float):
        layer_type, fp, material_parameter, layer_index = [], [], [], []

        for layer in self.layers:
            layer_type.append(layer.__class__.__name__)
            fp.append(layer._fp)
            material_parameter.append("Fixed")
            layer_index.append([layer._m.n(wavelength, *layer._mp)])

        fiber = Fiber(
            layer_radius=self.layer_radius,
            layer_type=layer_type,
            fp=fp,
            material_parameter=material_parameter,
            layer_index=layer_index,
            layer_names=self.layer_names,
            Cutoff=self._cutoff.__class__,
            Neff=self._neff.__class__
        )

        return fiber

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
        radius = abs(radius)
        for idx, r_ in enumerate(self.layer_radius):
            if radius < r_:
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
        return self.get_outer_radius(layer_idx) - self.get_inner_radius(layer_idx)

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

    def find_cutoff_solver(self):
        cutoff = solver.solver.FiberSolver
        if all(isinstance(layer, geometry.StepIndex)
               for layer in self.layers):
            nlayers = len(self)
            if nlayers == 2:  # SSIF
                cutoff = solver.ssif.Cutoff
            elif nlayers == 3:
                cutoff = solver.tlsif.Cutoff
        return cutoff

    def find_n_eff_solver(self) -> float:
        neff = solver.mlsif.Neff
        if all(isinstance(layer, geometry.StepIndex)
               for layer in self.layers):
            nlayers = len(self)
            if nlayers == 2:  # SSIF
                neff = solver.ssif.Neff
        return neff

    def set_solvers(self, Cutoff=None, Neff=None) -> None:
        assert Cutoff is None or issubclass(Cutoff, solver.FiberSolver)
        assert Neff is None or issubclass(Neff, solver.FiberSolver)
        if Cutoff is None:
            Cutoff = self.find_cutoff_solver()
        self._cutoff = Cutoff(self)
        if Neff is None:
            Neff = self.find_n_eff_solver()
        self._neff = Neff(self)

    def set_ne_cache(self, wavelength: float, mode, neff) -> None:
        try:
            self.ne_cache[wavelength][mode] = neff
        except KeyError:
            self.ne_cache[wavelength] = {mode: neff}

    def get_NA(self, wavelength: float) -> float:
        n1 = max(layer.get_maximum_index(wavelength) for layer in self.layers)
        n2 = self.get_minimum_index(-1, wavelength)
        return sqrt(n1 * n1 - n2 * n2)

    def get_V0(self, wavelength: float):
        wavelength = Wavelength(wavelength)
        return wavelength.k0 * self.get_inner_radius(-1) * self.get_NA(wavelength)

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
        if isinf(V0):
            return 0

        def model(x):
            return 2 * pi / V0 * b * self.get_NA(x)

        b = self.get_inner_radius(-1)

        wavelength = model(1.55e-6)
        if abs(wavelength - model(wavelength)) > tol:
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

    def cutoff(self, mode: Mode):
        try:
            return self.cutoff_cache[mode]
        except KeyError:
            cutoff = self._cutoff(mode=mode)
            self.cutoff_cache[mode] = cutoff
            return cutoff

    def get_cutoff_wavelength(self, mode: Mode) -> float:
        cutoff = self.cutoff(mode=mode)
        return self.V0_to_wavelength(cutoff)

    def neff(self, mode: Mode, wavelength: float, delta: float = 1e-6, lowbound=None) -> float:
        try:
            return self.ne_cache[wavelength][mode]
        except KeyError:
            neff = self._neff(Wavelength(wavelength), mode, delta, lowbound)
            self.set_ne_cache(wavelength, mode, neff)
            return neff

    def beta(self, omega: float, mode: Mode, p: float = 0, delta: float = 1e-6, lowbound=None):
        wl = Wavelength(omega=omega)
        if p == 0:
            neff = self.neff(mode, wl, delta, lowbound)
            return neff * wl.k0

        m = 5
        j = (m - 1) // 2
        h = 1e12  # This value is critical for accurate computation
        lb = lowbound
        for i in range(m - 1, -1, -1):
            # Precompute neff using previous wavelength
            o = omega + (i - j) * h
            wavelength = Wavelength(omega=o)
            lb = self.neff(mode, wavelength, delta, lb) + delta * 1.1

        return derivative(
            self.beta, omega, p, m, j, h, mode, 0, delta, lowbound)

    def get_normalized_beta(self, mode: Mode, wavelength: float, delta: float = 1e-6, lowbound=None) -> float:
        """
        Normalized propagation constant
        """
        neff = self.neff(
            mode=mode,
            wavelength=wavelength,
            delta=delta,
            lowbound=lowbound
        )

        overall_maximum_index = max(layer.get_maximum_index(wavelength) for layer in self.layers)

        minimum_index = self.get_minimum_index(-1, wavelength)

        numerator = neff**2 - minimum_index**2

        denominator = overall_maximum_index**2 - minimum_index**2

        return numerator / denominator

    def get_phase_velocity(self, mode: Mode, wavelength: float, delta: float = 1e-6, lowbound=None) -> float:
        n_eff = self.neff(
            mode=mode,
            wavelength=wavelength,
            delta=delta,
            lowbound=lowbound
        )

        return constants.c / n_eff

    def get_group_index(self, mode: Mode, wavelength: float, delta: float = 1e-6, lowbound=None) -> float:
        wavelength = Wavelength(wavelength)
        beta = self.beta(
            omega=wavelength.omega,
            mode=mode,
            p=1,
            delta=delta,
            lowbound=lowbound
        )

        return beta * constants.c

    def get_groupe_velocity(self, mode: Mode, wavelength: float, delta: float = 1e-6, lowbound=None) -> float:
        wavelength = Wavelength(wavelength)
        beta = self.beta(
            omega=wavelength.omega,
            mode=mode,
            p=1,
            delta=delta,
            lowbound=lowbound
        )

        return 1 / beta

    def D(self, mode: Mode, wavelength: float, delta: float = 1e-6, lowbound=None) -> float:
        wavelength = Wavelength(wavelength)
        beta = self.beta(
            omega=wavelength.omega,
            mode=mode,
            p=2,
            delta=delta,
            lowbound=lowbound
        )

        return -beta * 2 * pi * constants.c * 1e6 / wavelength**2

    def S(self, mode: Mode, wavelength: float, delta: float = 1e-6, lowbound=None) -> float:
        wavelength = Wavelength(wavelength)
        beta = self.beta(
            omega=wavelength.omega,
            mode=mode,
            p=3,
            delta=delta,
            lowbound=lowbound
        )

        return 1e-3 * beta * (2 * pi * constants.c / wavelength**2)**2

    def get_vectorial_modes(self, wavelength: float, numax=None, mmax=None, delta=1e-6):
        families = (ModeFamily.HE, ModeFamily.EH, ModeFamily.TE, ModeFamily.TM)

        modes = self.get_modes_from_familly(
            families=families,
            wavelength=wavelength,
            numax=numax,
            mmax=mmax,
            delta=delta
        )

        return modes

    def get_LP_modes(self, wavelength: float, ellmax=None, mmax=None, delta=1e-6):
        families = (ModeFamily.LP,)

        modes = self.get_modes_from_familly(
            families=families,
            wavelength=wavelength,
            ellmax=ellmax,
            mmax=mmax,
            delta=delta
        )

        return modes

    def get_modes_from_familly(self, families, wavelength: float, numax=None, mmax=None, delta: float = 1e-6):
        """
        Find all modes of given families, within given constraints

        """
        modes = set()
        v0 = self.get_V0(wavelength=wavelength)
        for fam in families:
            for nu in count(0):
                try:
                    _mmax = mmax[nu]
                except IndexError:
                    _mmax = mmax[-1]
                except TypeError:
                    _mmax = mmax

                if (fam is ModeFamily.TE or fam is ModeFamily.TM) and nu > 0:
                    break
                if (fam is ModeFamily.HE or fam is ModeFamily.EH) and nu == 0:
                    continue
                if numax is not None and nu > numax:
                    break
                for m in count(1):
                    if _mmax is not None and m > _mmax:
                        break
                    mode = Mode(fam, nu, m)
                    try:
                        co = self.cutoff(mode)
                        if co > v0:
                            break
                    except (NotImplementedError, ValueError):
                        neff = self.neff(mode, wavelength, delta)
                        if isnan(neff):
                            break
                    modes.add(mode)
                if m == 1:
                    break
        return modes

    def field(self, mode: Mode, wavelength: float, r, n_point: int = 101):
        """
        Return electro-magnetic field.

        """
        return Field(self, mode, wavelength, r, n_point)

    @lru_cache(maxsize=None)
    def _rfield(self, mode: Mode, wavelength: float, r):
        neff = self.neff(
            mode=mode,
            wavelength=wavelength
        )

        fct = {ModeFamily.LP: self._neff._lpfield,
               ModeFamily.TE: self._neff._tefield,
               ModeFamily.TM: self._neff._tmfield,
               ModeFamily.EH: self._neff._ehfield,
               ModeFamily.HE: self._neff._hefield}

        return fct[mode.family](wavelength, mode.nu, neff, r)
