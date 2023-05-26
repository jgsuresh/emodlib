import os

from .._emodlib_py.malaria import Infection, IntrahostComponent, Susceptibility
from ..params import Params, set_params, update_params


def params_from_default_file():
    path = os.path.join(os.path.realpath(os.path.dirname(__file__)), "config.yml")
    return Params.from_yaml(path)


IntrahostComponent.default_params = params_from_default_file()

# monkey-patch set_params + update_params
# nested on top of defaults + current values, respectively
IntrahostComponent.set_params = set_params
IntrahostComponent.update_params = update_params

# initialize default parameters
IntrahostComponent.set_params()


__all__ = ["IntrahostComponent", "Susceptibility", "Infection"]
