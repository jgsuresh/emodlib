from typing import Any, Dict, TypeVar

import yaml


class Params(dict):
    @classmethod
    def from_yaml(cls, path):
        with open(path) as cfg:
            params = yaml.load(cfg, Loader=yaml.FullLoader)
        return cls(params)

    @property
    def yaml(self):
        return yaml.dump(self)


yaml.add_representer(
    Params,
    lambda dumper, data: dumper.represent_mapping(
        "tag:yaml.org,2002:map", data.items()
    ),
)


KeyType = TypeVar("KeyType")


def deep_update(
    mapping: Dict[KeyType, Any], *updating_mappings: Dict[KeyType, Any]
) -> Dict[KeyType, Any]:
    updated_mapping = mapping.copy()
    for updating_mapping in updating_mappings:
        for k, v in updating_mapping.items():
            if (
                k in updated_mapping
                and isinstance(updated_mapping[k], dict)
                and isinstance(v, dict)
            ):
                updated_mapping[k] = deep_update(updated_mapping[k], v)
            else:
                updated_mapping[k] = v
    return Params(updated_mapping)


@classmethod
def update_params(cls, params={}):
    """Nested update on top of current parameters"""
    cfg = deep_update(cls.params, params)
    cls.params = cfg
    cls._configure_from_params(cfg)  # setting static variables in bound C++ classes


@classmethod
def set_params(cls, params={}):
    """Nested update on top of default parameters"""
    cfg = deep_update(cls.default_params, params)
    cls.params = cfg
    cls._configure_from_params(cfg)  # setting static variables in bound C++ classes
