"""Public unit contract and numerical parity across the SI solver boundary."""

import numpy as np
import pytest
import pint
from TypedUnit import ureg as typed_registry

from PyFiberModes import (
    Fiber,
    FiberFactory,
    LP01,
    Field,
    FusedSilica,
    ConstantIndex,
    CartesianCoordinates,
    CoupledModeSystem,
    DesignParameter,
    optimize_fiber,
    coupling_matrix,
    overlap,
    load_fiber,
    sweep_modes,
)
from PyFiberModes.source import Source
from PyFiberModes.units import Length, RefractiveIndex, MODAL_UNITS, ureg


@pytest.fixture
def fiber():
    """Load a reference fiber with an explicitly nanometer-valued wavelength."""
    return load_fiber(fiber_name="SMF28", wavelength=1550 * ureg.nanometer)


def test_shared_registry_and_typed_units(fiber):
    """Use exactly the same registry and type conventions as PyMieSim."""
    assert ureg is typed_registry
    assert isinstance(fiber.wavelength, Length)
    assert isinstance(fiber.get_effective_index(LP01), RefractiveIndex)
    assert fiber.wavelength.to("nanometer").magnitude == pytest.approx(1550)


@pytest.mark.parametrize("value", [1550e-9, None, [1550e-9]])
def test_wavelength_requires_quantity(value):
    """Never silently interpret a bare numerical wavelength as meters."""
    with pytest.raises(TypeError, match="wavelength"):
        Fiber(wavelength=value)
    with pytest.raises(TypeError, match="wavelength"):
        load_fiber(fiber_name="SMF28", wavelength=value)


@pytest.mark.parametrize("value", [0, -1, np.inf, np.nan])
def test_wavelength_must_be_positive_finite_scalar(value, fiber):
    """Preserve physical validation after unit conversion."""
    with pytest.raises(ValueError):
        Fiber(wavelength=value * ureg.meter)
    before = fiber.wavelength
    with pytest.raises(ValueError):
        fiber.update_wavelength(value * ureg.meter)
    assert fiber.wavelength == before


def test_wrong_dimensions_and_array_scalar_rejected(fiber):
    """Reject incompatible physical dimensions and arrays for scalar inputs."""
    with pytest.raises(pint.DimensionalityError):
        Fiber(wavelength=1 * ureg.second)
    with pytest.raises(ValueError, match="scalar"):
        Fiber(wavelength=np.array([1, 2]) * ureg.micrometer)
    with pytest.raises(pint.DimensionalityError):
        fiber.get_index_profile(np.array([1, 2]) * ureg.second)
    with pytest.raises(TypeError):
        fiber.wavelength = 1e-6


def test_explicit_layer_units_and_geometry_mutations():
    """Normalize geometry/index quantities and validate later mutations."""
    fiber = Fiber(wavelength=1.55 * ureg.micrometer)
    with pytest.raises(TypeError, match="radius"):
        fiber.add_layer(name="core", radius=4e-6, index=1.46 * ureg.dimensionless)
    with pytest.raises(TypeError, match="index"):
        fiber.add_layer(name="core", radius=4 * ureg.micrometer, index=1.46)
    with pytest.raises(pint.DimensionalityError):
        fiber.add_layer(name="core", radius=4 * ureg.second, index=1.46 * ureg.dimensionless)
    fiber.add_layer(name="core", radius=4 * ureg.micrometer, index=1.46 * ureg.dimensionless)
    fiber.add_layer(name="cladding", radius=np.inf * ureg.meter, index=1.44 * ureg.dimensionless)
    fiber.initialize_layers()
    assert fiber.radius.to("micrometer").magnitude == pytest.approx(4)
    before = fiber.get_effective_index(LP01)
    fiber[0].refractive_index = 1.47 * ureg.dimensionless
    assert fiber.get_effective_index(LP01) > before
    with pytest.raises(TypeError):
        fiber[0].radius_out = 5e-6
    assert fiber.scale(factor=0.5).radius == fiber.radius / 2


@pytest.mark.parametrize("metric,unit", MODAL_UNITS.items())
def test_modal_results_have_units_and_preserve_numerical_values(fiber, metric, unit):
    """Protect SI results, including dispersion's historical reporting scales."""
    method = f"get_{metric}"
    actual = getattr(fiber, method)(LP01)
    reference = getattr(fiber._si, method)(LP01)
    assert isinstance(actual, pint.Quantity)
    assert actual.units == ureg.Unit(unit)
    assert actual.to(unit).magnitude == pytest.approx(reference, nan_ok=True)


def test_units_convert_and_failures_keep_diagnostics(fiber):
    """Expose convertible velocities, delays, dispersion and missing modes."""
    assert fiber.get_propagation_constant(LP01).to("1 / kilometer").magnitude > 0
    assert fiber.get_phase_velocity(LP01).to("kilometer / second").magnitude > 0
    assert fiber.get_group_velocity_dispersion(LP01).check("[time] ** 2 / [length]")
    assert np.isinf(fiber.get_mode_cutoff_wavelength(LP01).magnitude)
    solved = fiber.solve_effective_index(LP01)
    assert solved.unwrap() == fiber.get_effective_index(LP01)
    assert fiber.analysis.effective_index_result(LP01).value == solved.value


def test_source_conversions():
    """Attach physical units to source frequency, angular frequency and k0."""
    source = Source(wavelength=1 * ureg.micrometer)
    assert source.frequency.to("terahertz").magnitude == pytest.approx(299.792458)
    assert source.omega.to("radian / second").magnitude == pytest.approx(2 * np.pi * 299792458e6)
    assert source.wavenumber.to("radian / meter").magnitude == pytest.approx(2 * np.pi * 1e6)
    with pytest.raises(TypeError):
        source.wavelength = 1e-6


def test_factory_preserves_quantity_arrays():
    """Generate fibers from radii expressed in micrometers."""
    factory = FiberFactory(wavelength=1550 * ureg.nanometer)
    factory.add_layer(
        name="core", radius=np.array([4, 5]) * ureg.micrometer, index=1.46 * ureg.dimensionless
    )
    factory.add_layer(name="cladding", index=1.44 * ureg.dimensionless)
    fibers = list(factory)
    assert len(fibers) == 2
    assert fibers[1].radius.to("micrometer").magnitude == pytest.approx(5)


def test_sweep_and_records_retain_units(fiber):
    """Return quantity-valued arrays and records without mutating the template."""
    result = sweep_modes(fiber, np.array([1.31, 1.55]) * ureg.micrometer, modes=[LP01])
    assert result.parameters.units == ureg.micrometer
    assert result["propagation_constant"].units == 1 / ureg.meter
    record = result.as_records()[0]
    assert isinstance(record["wavelength"], Length)
    assert record["propagation_constant"].check("1 / [length]")
    assert fiber.wavelength.to("nanometer").magnitude == pytest.approx(1550)
    with pytest.raises(TypeError):
        fiber.sweep([1e-6, 2e-6], modes=[LP01])
    with pytest.raises(pint.DimensionalityError):
        fiber.sweep(np.array([1, 2]) * ureg.second, modes=[LP01])


def test_custom_parameter_sweep_has_consistent_units(fiber):
    """Pass quantity-valued custom sweep coordinates to user setters."""
    result = fiber.sweep(
        [4 * ureg.micrometer, 5000 * ureg.nanometer],
        modes=[LP01],
        parameter_name="radius",
        setter=lambda current, value: (
            setattr(current[0], "radius_out", value),
            setattr(current[1], "radius_in", value),
        ),
    )
    np.testing.assert_allclose(result.parameters.to("micrometer").magnitude, [4, 5])


def test_fields_and_coordinates_have_physical_units(fiber):
    """Return field amplitudes, positions, areas and integrated power with units."""
    field = fiber.get_mode_field(mode=LP01, limit=15 * ureg.micrometer, n_point=15)
    assert isinstance(field, Field)
    assert field.cartesian_coordinates.x.check("[length]")
    assert field.cartesian_coordinates.dx.check("[length]")
    assert field.cylindrical_coordinates.phi.units == ureg.radian
    assert field.Ex(phi=90 * ureg.degree).check("volt / meter")
    assert field.Hy().check("ampere / meter")
    assert field.get_effective_area().check("[length] ** 2")
    assert field.get_power().check("watt")
    assert field.get_poynting_vector()[2].check("watt / meter ** 2")
    assert field.get_integrale_square(field.Ex()).check("volt ** 2")
    assert field.get_confinement_factor(4 * ureg.micrometer).dimensionless
    assert overlap(field, field).magnitude == pytest.approx(1)
    assert coupling_matrix([field, field]).dimensionless
    assert coupling_matrix([field, field], perturbation=2 / ureg.meter).check("1 / [length]")
    with pytest.raises(TypeError):
        field.Ex(phi=0)
    electric, magnetic = fiber.get_radial_fields(LP01, np.array([0, 2]) * ureg.micrometer)
    assert electric.rho.check("volt / meter")
    assert magnetic.z.check("ampere / meter")


def test_public_coordinate_input_contract():
    """Reject bare coordinate grids and preserve physical spacing."""
    with pytest.raises(TypeError):
        CartesianCoordinates(x=[0], y=[0], z=[0])
    coordinates = CartesianCoordinates.generate_from_square(
        length=20 * ureg.micrometer, n_points=11
    )
    assert coordinates.dx.to("micrometer").magnitude == pytest.approx(2)


def test_material_models_require_lengths_and_return_indices():
    """Material models use the same public unit convention as fibers."""
    assert FusedSilica().refractive_index(1550 * ureg.nanometer).magnitude == pytest.approx(
        1.4440236217
    )
    assert (
        ConstantIndex(value=1.45 * ureg.dimensionless)
        .refractive_index(1 * ureg.micrometer)
        .magnitude
        == 1.45
    )
    with pytest.raises(TypeError):
        FusedSilica().refractive_index(1550e-9)
    with pytest.raises(TypeError):
        ConstantIndex(value=1.45)


def test_propagation_checks_inverse_lengths_and_returns_positions():
    """Only convert physical propagation values to SI at the ODE boundary."""
    system = CoupledModeSystem(
        betas=np.array([0, 0]) / ureg.meter, coupling=np.array([[0, 0.1], [0.1, 0]]) / ureg.meter
    )
    result = system.propagate(initial=[1, 0], z=np.linspace(0, 10, 11) * ureg.centimeter)
    assert result.z.units == ureg.meter
    assert result.z[-1].magnitude == pytest.approx(0.1)
    assert result.amplitudes.dimensionless
    np.testing.assert_allclose(result.powers.sum(axis=1).magnitude, 1, atol=1e-7)
    with pytest.raises(TypeError):
        CoupledModeSystem(betas=[0], coupling=[[0]])
    with pytest.raises(pint.DimensionalityError):
        system.propagate(initial=[1, 0], z=np.array([0, 1]) * ureg.second)


def test_optimization_preserves_parameter_units(fiber):
    """Normalize mixed-unit bounds for SciPy and restore quantities for setters."""
    parameter = DesignParameter.layer_radius(0, (3 * ureg.micrometer, 5000 * ureg.nanometer))
    result = optimize_fiber(
        fiber,
        [parameter],
        objective=lambda current: (current[0].radius_out.to("micrometer").magnitude - 4.5) ** 2,
        method="L-BFGS-B",
    )
    assert result.parameters[parameter.name].to("micrometer").magnitude == pytest.approx(
        4.5, abs=1e-4
    )
    assert fiber[0].radius_out.to("micrometer").magnitude == pytest.approx(4.1)
