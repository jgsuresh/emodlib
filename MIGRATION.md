# Migration Guide: v0.0.4 to v0.1.0

## Breaking Changes

Version 0.1.0 introduces a thread-safe, instance-based configuration system. The previous static configuration pattern has been removed.

### Configuration API Changes

**v0.0.4 (Old)**
```python
from emodlib.malaria import IntrahostComponent

# Set params globally (affected all instances)
IntrahostComponent.set_params({'Run_Number': 42})
IntrahostComponent.update_params({'infection_params': {'Antigen_Switch_Rate': 1e-9}})

# Create instance (used global params)
ic = IntrahostComponent.create()
```

**v0.1.0 (New)**
```python
from emodlib.malaria import IntrahostComponent, create_config

# Create config instance (thread-safe, independent)
config = create_config({'Run_Number': 42})

# Or with nested params
config = create_config({
    'Run_Number': 42,
    'infection_params': {'Antigen_Switch_Rate': 1e-9}
})

# Create instance with config
ic = IntrahostComponent.create(config)
```

### Removed Functions

| Removed | Replacement |
|---------|-------------|
| `IntrahostComponent.set_params(params)` | `config = create_config(params)` |
| `IntrahostComponent.update_params(params)` | `config.update(params)` |
| `IntrahostComponent.params` | `config.yaml` or config properties |

### Incremental Updates with config.update()

The `update()` method provides the same incremental merge behavior as the old `update_params()`:

```python
# Create config with initial params
config = create_config({'Run_Number': 1, 'Max_Individual_Infections': 3})

# Later, update just one param (keeps Run_Number=1, Max_Individual_Infections=3)
config.update({'infection_params': {'Antigen_Switch_Rate': 1e-8}})

# Method chaining also works
config.update({'Run_Number': 2}).update({'infection_params': {'Antigen_Switch_Rate': 1e-9}})
```

### Creating Susceptibility and Infection Objects

**v0.0.4**
```python
from emodlib.malaria import Susceptibility, Infection

s = Susceptibility.create()
inf = Infection.create(susceptibility=s, hepatocytes=1)
```

**v0.1.0**
```python
from emodlib.malaria import Susceptibility, Infection, create_config

config = create_config()
s = Susceptibility.create(config)
inf = Infection.create(susceptibility=s, config=config, hepatocytes=1)
```

### Viewing Configuration

**v0.0.4**
```python
print(IntrahostComponent.params.yaml)
```

**v0.1.0**
```python
config = create_config()
print(config.yaml)
```

### Multi-Simulation Workflows

**v0.0.4** - Simulations with different parameters required careful sequencing:
```python
IntrahostComponent.set_params({'Run_Number': 1})
ic1 = IntrahostComponent.create()

IntrahostComponent.set_params({'Run_Number': 2})  # Changed global state!
ic2 = IntrahostComponent.create()
```

**v0.1.0** - Each config is independent (thread-safe):
```python
config1 = create_config({'Run_Number': 1})
config2 = create_config({'Run_Number': 2})

ic1 = IntrahostComponent.create(config1)
ic2 = IntrahostComponent.create(config2)

# Can run in parallel without interference
```

## Quick Migration Checklist

1. Add `create_config` to imports:
   ```python
   from emodlib.malaria import IntrahostComponent, create_config
   ```

2. Replace `set_params()` with `create_config()`:
   ```python
   # Before
   IntrahostComponent.set_params({'Run_Number': 42})

   # After
   config = create_config({'Run_Number': 42})
   ```

3. Replace `update_params()` with `config.update()`:
   ```python
   # Before
   IntrahostComponent.update_params({'infection_params': {'X': 1}})

   # After
   config.update({'infection_params': {'X': 1}})
   ```

4. Pass config to `create()`:
   ```python
   # Before
   ic = IntrahostComponent.create()

   # After
   ic = IntrahostComponent.create(config)
   ```

5. For `Susceptibility.create()` and `Infection.create()`, add config parameter.
