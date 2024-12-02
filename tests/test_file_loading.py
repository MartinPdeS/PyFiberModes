import pytest
import yaml
import numpy as np
from PyOptik import MaterialBank
from PyFiberModes.loader import (  # Replace `your_module_name` with the actual module name
    get_fiber_file_path,
    load_yaml_configuration,
    calculate_layer_index,
    process_layers,
    cleanup_layers,
    reorder_layers,
    load_fiber_as_dict,
)

WAVELENGTH = 1550e-9
FUSEDSILICA_RI = MaterialBank.fused_silica.compute_refractive_index(WAVELENGTH)
CROWN_RI = MaterialBank.crown.compute_refractive_index(WAVELENGTH)


@pytest.fixture
def mock_fiber_file(tmp_path):
    """Fixture to create a mock fiber YAML file."""
    fiber_content = {
        'layers': {
            '1': {'material': 'fused_silica'},
            '2': {'NA': 0.1, 'index': 1.6},
            '3': {'index': 1.7},
        }
    }
    mock_file = tmp_path / "fiber_files" / "mock_fiber.yaml"
    mock_file.parent.mkdir(parents=True, exist_ok=True)
    with mock_file.open('w') as f:
        yaml.dump(fiber_content, f)
    return mock_file


def test_get_fiber_file_path():
    """Test file path resolution for a fiber."""
    fiber_name = "example_fiber"
    assert get_fiber_file_path(fiber_name), "Error while loading fiber file."


def test_load_yaml_configuration(mock_fiber_file):
    """Test loading of YAML configuration."""
    config = load_yaml_configuration(mock_fiber_file)
    assert "layers" in config
    assert config["layers"]["1"]["material"] == "fused_silica"


def test_calculate_layer_index():
    """Test calculation of layer index."""
    layer = {'material': 'fused_silica'}
    assert calculate_layer_index(layer, wavelength=WAVELENGTH) == FUSEDSILICA_RI

    layer_with_na = {'NA': 0.1}
    outer_layer = {'index': 1.5}
    assert np.isclose(calculate_layer_index(layer_with_na, wavelength=None, outer_layer=outer_layer), 1.503329637837291)


def test_process_layers():
    """Test processing of layers."""
    layers = {
        '1': {'material': 'fused_silica'},
        '2': {'NA': 0.1},
        '3': {'index': 1.7},
    }

    processed = process_layers(layers, wavelength=WAVELENGTH)

    assert '1' in processed
    assert processed['1']['index'] == FUSEDSILICA_RI
    assert np.isclose(processed['2']['index'], 1.447482027535058)
    assert processed['3']['index'] == 1.7


def test_cleanup_layers():
    """Test cleanup of layers."""
    layers = {
        '1': {'index': 1.55, 'material': 'Material1'},
        '2': {'index': 1.503, 'NA': 0.1},
    }
    cleaned = cleanup_layers(layers)

    for layer in cleaned.values():
        assert 'material' not in layer
        assert 'NA' not in layer


def test_reorder_layers():
    """Test reordering of layers."""
    layers = {
        '1': {'index': 1.5},
        '2': {'index': 1.6},
        '3': {'index': 1.7},
    }

    reordered = reorder_layers(layers, order='out-to-in')
    expected_keys = list(reversed(layers.keys()))
    assert list(reordered.keys()) == expected_keys

    reordered = reorder_layers(layers, order='in-to-out')
    assert list(reordered.keys()) == list(layers.keys())


def test_load_fiber_as_dict(mock_fiber_file):
    """Test loading and processing of fiber configuration."""
    result = load_fiber_as_dict('SMF28', wavelength=WAVELENGTH, order='out-to-in')

    assert 'layers' in result
    layers = result['layers']
    assert len(layers) == 2
    assert layers['layer_1']['index'] == 1.444023621703261
    assert 'material' not in layers['layer_1']
    assert 'NA' not in layers['layer_2']


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
