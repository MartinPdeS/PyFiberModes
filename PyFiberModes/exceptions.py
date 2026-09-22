"""Domain-specific exceptions raised by PyFiberModes."""


class PyFiberModesError(Exception):
    """Base class for package-specific errors."""


class ValidationError(PyFiberModesError, ValueError):
    """Indicate invalid physical geometry or solver configuration."""


class SolverError(PyFiberModesError, RuntimeError):
    """Indicate that a modal solver could not complete normally."""


class ConvergenceError(SolverError):
    """Indicate that a numerical solver did not converge to a valid root."""
