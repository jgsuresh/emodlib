# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2025-12-16

### Changed
- **BREAKING**: Refactored configuration from static class variables to instance-based `MalariaConfig` for thread-safe parallel simulations
- `IntrahostComponent.create()` now requires a config argument
- Removed `IntrahostComponent.set_params()` and `update_params()` static methods
- Added `create_config(params)` function to create config instances
- Added `config.update(params)` method for incremental parameter updates
- Added `config.yaml` property for YAML representation

### Added
- `MalariaConfig` class that owns all configuration state including RNG and SUID generator
- `create_config()` helper function with defaults merging
- Thread-safety tests for parallel simulation isolation
- [Migration guide](MIGRATION.md) for upgrading from v0.0.4

## [0.0.4] - 2025-11-12

### Added
- Python 3.12 and 3.13 support and compatibility

### Changed
- Updated minimum Python version from 3.7 to 3.9 (Python 3.7 and 3.8 are EOL)
- Updated CI/CD workflows to build and test with Python 3.9, 3.11, and 3.12
- Updated package classifiers to support Python 3.9-3.13
- Updated artifact upload/download actions from v3 to v4
- Updated cibuildwheel to v2.21.3 for improved Windows build support
- Updated pypa/gh-action-pypi-publish from v1.6.4 to v1.13.0

## [0.0.3] - Previous release
