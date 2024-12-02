#!/usr/bin/env python
# -*- coding: utf-8 -*-

import yaml
import numpy as np
from pathlib import Path
from PyOptik import MaterialBank


def get_fiber_file_path(fiber_name: str) -> Path:
    """
    Get the file path for the specified fiber configuration file.

    Parameters
    ----------
    fiber_name : str
        The name of the fiber.

    Returns
    -------
    Path
        Path object pointing to the fiber's YAML file.
    """
    return Path(__file__).parent / 'fiber_files' / f'{fiber_name}.yaml'


def load_yaml_configuration(file_path: Path) -> dict:
    """
    Load and parse a YAML configuration file.

    Parameters
    ----------
    file_path : Path
        Path object pointing to the YAML file.

    Returns
    -------
    dict
        Parsed content of the YAML file.
    """
    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")

    with file_path.open('r') as file:
        return yaml.safe_load(file)


def calculate_layer_index(layer: dict, wavelength: float, outer_layer: dict = None) -> float:
    """
    Calculate the refractive index for a layer.

    Parameters
    ----------
    layer : dict
        Dictionary containing layer properties.
    wavelength : float
        Wavelength for material index calculation (optional).
    outer_layer : dict, optional
        The outer layer dictionary, if any, for NA-based calculations.

    Returns
    -------
    float
        Computed refractive index for the layer.
    """
    if 'material' in layer and wavelength:
        return getattr(MaterialBank, layer['material']).compute_refractive_index(wavelength)
    elif 'NA' in layer and outer_layer:
        return np.sqrt(layer['NA']**2 + outer_layer['index']**2)
    return layer.get('index')


def process_layers(layers: dict, wavelength: float = None) -> dict:
    """
    Process and calculate refractive indices for all layers.

    Parameters
    ----------
    layers : dict
        Dictionary of layers from the configuration.
    wavelength : float, optional
        Wavelength for material index calculation.

    Returns
    -------
    dict
        Processed layers with computed indices.
    """
    processed_layers = {}
    outer_layer = None

    for idx, layer in layers.items():
        layer_index = calculate_layer_index(layer, wavelength, outer_layer)
        processed_layers[idx] = {**layer, 'index': layer_index}
        outer_layer = processed_layers[idx]

    return processed_layers


def cleanup_layers(layers: dict) -> dict:
    """
    Remove unnecessary keys from layer dictionaries.

    Parameters
    ----------
    layers : dict
        Dictionary of processed layers.

    Returns
    -------
    dict
        Cleaned-up dictionary of layers.
    """
    for layer in layers.values():
        layer.pop('NA', None)
        layer.pop('material', None)
    return layers


def reorder_layers(layers: dict, order: str) -> dict:
    """
    Reorder layers based on the specified order.

    Parameters
    ----------
    layers : dict
        Dictionary of layers to reorder.
    order : str
        Order to apply ('in-to-out' or 'out-to-in').

    Returns
    -------
    dict
        Reordered dictionary of layers.
    """
    if order == 'out-to-in':
        return dict(reversed(list(layers.items())))
    return layers


def load_fiber_as_dict(fiber_name: str, wavelength: float = None, order: str = 'in-to-out') -> dict:
    """
    Load and process a fiber configuration file.

    Parameters
    ----------
    fiber_name : str
        Name of the fiber file (without extension).
    wavelength : float, optional
        Wavelength for refractive index calculations.
    order : str
        Order of the layers ('in-to-out' or 'out-to-in').

    Returns
    -------
    dict
        Dictionary with processed fiber layers.
    """
    file_path = get_fiber_file_path(fiber_name)

    if not file_path.exists():
        available_files = [f.stem for f in file_path.parent.glob('*.yaml')]
        available_list = '\n'.join(available_files)
        raise FileNotFoundError(
            f"Fiber file '{fiber_name}.yaml' not found. Available fibers:\n{available_list}"
        )

    config = load_yaml_configuration(file_path)
    layers = config.get('layers', {})

    processed_layers = process_layers(layers, wavelength)
    cleaned_layers = cleanup_layers(processed_layers)
    reordered_layers = reorder_layers(cleaned_layers, order)

    return {'layers': reordered_layers}
