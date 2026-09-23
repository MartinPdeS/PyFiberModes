"""Public API for circular optical-fiber mode analysis."""

from PyFiberModes.mode import Family, Mode
from PyFiberModes.mode_instances import (
    EH11,
    HE11,
    HE12,
    HE21,
    HE22,
    HE31,
    HE32,
    LP01,
    LP02,
    LP03,
    LP11,
    LP12,
    LP21,
    LP22,
    LP31,
    LP32,
    LP41,
    TE01,
    TM01,
)
from PyFiberModes.factory import FiberFactory
from PyFiberModes.fiber import Fiber, load_fiber
from PyFiberModes.field import Field
from PyFiberModes.coordinates import CartesianCoordinates, CylindricalCoordinates
from PyFiberModes.exceptions import (
    ConvergenceError,
    SolverError,
    UnsupportedGeometryError,
    ValidationError,
)
from PyFiberModes.models import LayerSpec, SolverSettings
from PyFiberModes.solver.results import SolverResult
from PyFiberModes.analysis import ModeSweepResult, find_modes, sweep_modes
from PyFiberModes.optimization import DesignParameter, DesignResult, optimize_fiber
from PyFiberModes.propagation import CoupledModeSystem, PropagationResult, coupling_matrix, overlap
from PyFiberModes.materials import (
    CallableIndex,
    ConstantIndex,
    FusedSilica,
    MaterialRegistry,
    RefractiveIndexModel,
)


try:
    from ._version import version as __version__  # noqa: F401

except ImportError:
    __version__ = "0.0.0"


__all__ = [
    'Mode',
    'Family',
    'Fiber',
    'load_fiber',
    'FiberFactory',
    'Field',
    'CartesianCoordinates', 'CylindricalCoordinates',
    'LayerSpec', 'SolverSettings', 'SolverResult',
    'ValidationError', 'SolverError', 'ConvergenceError', 'UnsupportedGeometryError',
    'ModeSweepResult', 'find_modes', 'sweep_modes',
    'DesignParameter', 'DesignResult', 'optimize_fiber',
    'CoupledModeSystem', 'PropagationResult', 'coupling_matrix', 'overlap',
    'RefractiveIndexModel', 'ConstantIndex', 'CallableIndex', 'FusedSilica',
    'MaterialRegistry',
    'HE11',
    'HE12',
    'HE21',
    'HE22',
    'HE31',
    'HE32',
    'EH11',
    'TE01',
    'TM01',
    'LP01',
    'LP11',
    'LP21',
    'LP02',
    'LP31',
    'LP12',
    'LP22',
    'LP03',
    'LP32',
    'LP41',
]
