import os

import yaml

from .._emodlib_py.malaria import (
    DebugAntigens,
    Infection,
    IntrahostComponent,
    MalariaConfig,
    Susceptibility,
)
from ..params import Params, deep_update


def _config_yaml_property(self):
    """Return YAML string representation of the config's source params."""
    return yaml.dump(dict(self._params), default_flow_style=False, sort_keys=False)


def _config_update(self, params):
    """
    Update this config with new parameters (nested merge on top of current).

    This is the instance-based equivalent of the old update_params() method.
    Modifies the config in-place and returns self for chaining.

    Args:
        params: Dict of parameters to merge on top of current config.
                Can be nested dict matching config.yml structure.

    Returns:
        self (for method chaining)

    Example:
        config = create_config({'Run_Number': 1})
        config.update({'infection_params': {'Antigen_Switch_Rate': 1e-8}})
    """
    merged = deep_update(self._params, params)
    self.configure(merged)
    self._params = merged
    return self


# Add yaml property and update method to MalariaConfig
MalariaConfig.yaml = property(_config_yaml_property)
MalariaConfig.update = _config_update


def params_from_default_file():
    """Load default parameters from the bundled config.yml file."""
    path = os.path.join(os.path.realpath(os.path.dirname(__file__)), "config.yml")
    return Params.from_yaml(path)


# Store default params at class level for reference
IntrahostComponent.default_params = params_from_default_file()


def create_config(params=None):
    """
    Create a new MalariaConfig instance.

    Args:
        params: Optional dict of parameters to override defaults.
                Can be a flat dict or nested dict matching config.yml structure.

    Returns:
        A configured MalariaConfig instance.

    Example:
        config = create_config({'Run_Number': 42})
        ic = IntrahostComponent.create(config)
    """
    if params is None:
        merged = Params(IntrahostComponent.default_params)
    else:
        merged = deep_update(IntrahostComponent.default_params, params)
    cfg = MalariaConfig.from_params(merged)
    cfg._params = merged  # Store for yaml property
    return cfg


__all__ = ["IntrahostComponent", "Susceptibility", "Infection", "MalariaConfig", "DebugAntigens", "create_config"]
