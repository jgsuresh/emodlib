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
