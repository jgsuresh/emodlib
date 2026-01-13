# Vector Population Module Guide

The `emodlib.vector` module provides mosquito population dynamics with cohort-based lifecycle tracking, temperature-dependent development, and malaria transmission.

## Overview

This module implements:
- Cohort-based vector population tracking (not individual mosquitoes)
- Temperature-dependent development via Arrhenius equations
- Age-dependent mortality via Styer et al. (2007) formula
- Human-vector-human malaria transmission cycle

Default parameters are calibrated for *Anopheles gambiae* in tropical conditions.

## Quick Start

```python
from emodlib.vector import VectorConfig, VectorPopulation

# Create configuration
config = VectorConfig()
config.carrying_capacity = 5000

# Initialize population
pop = VectorPopulation()
pop.initialize(config, seed=42)

# Run simulation
for day in range(100):
    pop.update(dt=1.0, temperature=25.0, human_infectiousness=[0.1, 0.2])

print(f"Adults: {pop.adult_population:.0f}")
print(f"Infectious: {pop.infectious_count:.0f}")
```

## Classes

### VectorConfig

Configuration parameters for vector population dynamics.

**Arrhenius Parameters** (temperature-dependent development):
- `aquatic_arrhenius1`, `aquatic_arrhenius2`: Larval development rate
- `infected_arrhenius1`, `infected_arrhenius2`: EIP (extrinsic incubation period)

**Mortality**:
- `adult_life_expectancy`: Mean adult lifespan (default: 21 days)
- `aquatic_mortality_rate`: Daily aquatic stage mortality (default: 0.2)
- `enable_age_dependent_mortality`: Use Styer senescence model (default: True)

**Feeding**:
- `anthropophily`: Fraction of feeds on humans (default: 0.9)
- `days_between_feeds`: Gonotrophic cycle length (default: 3 days)

**Transmission**:
- `transmission_rate`: V->H probability per infectious bite (default: 0.02)
- `acquire_modifier`: H->V probability modifier (default: 0.5)
- `sporozoites_per_bite`: Average sporozoites delivered (default: 11)

**Helper Methods**:
```python
config.get_larval_progress(temp_celsius, dt)   # Larval development per timestep
config.get_infected_progress(temp_celsius, dt) # EIP progress per timestep
config.get_survival_probability(age, base_mort) # Daily survival probability
VectorConfig.get_age_dependent_mortality(age)   # Styer mortality rate
```

### VectorCohort

Represents a group of mosquitoes sharing the same state.

**Attributes**:
- `population`: Number of mosquitoes (float for fractional mortality)
- `age`: Age in days
- `progress`: Development progress [0, 1+]
- `state`: VectorState enum
- `days_since_infection`: Days since infected (-1 if susceptible)

**State Queries**:
```python
cohort.is_susceptible()  # Can be infected
cohort.is_exposed()      # Infected, in EIP
cohort.is_infectious()   # Can transmit
cohort.is_aquatic()      # Egg/larva/immature
cohort.is_adult()        # Any adult state
```

### VectorPopulation

Manages cohort queues through lifecycle stages.

**Lifecycle Stages**:
```
EGG -> LARVA -> IMMATURE -> ADULT -> INFECTED -> INFECTIOUS
```

**Key Methods**:
```python
pop.initialize(config, seed)
pop.update(dt, temperature, human_infectiousness)
pop.get_sporozoite_challenges(n_humans)
```

**Population Queries**:
```python
pop.adult_population      # All adults (susceptible + exposed + infectious)
pop.susceptible_count     # Can be infected
pop.exposed_count         # In EIP
pop.infectious_count      # Can transmit
pop.infectious_fraction   # infectious / adults
```

### TransmissionSimulation

Coordinates human intrahost models with vector population.

```python
from emodlib.malaria import MalariaConfig
from emodlib.vector import VectorConfig, TransmissionSimulation

malaria_config = MalariaConfig()
vector_config = VectorConfig()
vector_config.carrying_capacity = 2000

sim = TransmissionSimulation()
sim.initialize(malaria_config, vector_config, n_humans=100, seed=42)

# Seed initial infection
sim.seed_infections(n_infections=5)

# Run simulation
for day in range(365):
    sim.update(dt=1.0, temperature=25.0)

    if day % 30 == 0:
        print(f"Day {day}: prevalence={sim.prevalence:.1%}, "
              f"EIR={sim.vector_infectious_fraction:.3f}")
```

## Temperature-Dependent Development

Development rates follow Arrhenius equations:

```
progress = A * exp(-B / (T + 273.15)) * dt
```

**EIP (Extrinsic Incubation Period)**:
| Temperature | EIP Duration |
|-------------|--------------|
| 20C         | ~25 days     |
| 25C         | ~19 days     |
| 30C         | ~16 days     |
| 35C         | ~13 days     |

## Age-Dependent Mortality

Styer et al. (2007) senescence model:

```
e = exp(0.2 * age)
mortality = 0.006 * e / (1 + 0.045 * (e - 1))
```

Combined with base mortality (1/life_expectancy), daily survival probability:

```python
survival = exp(-(base_mortality + styer_mortality))
```

## Transmission Cycle

1. **H->V**: Susceptible mosquitoes feeding on infectious humans become infected
   - `prob = biting_rate * anthropophily * infectiousness * acquire_modifier`

2. **EIP**: Infected mosquitoes progress through extrinsic incubation
   - Temperature-dependent via Arrhenius equation
   - Become infectious when progress >= 1.0

3. **V->H**: Infectious mosquitoes challenge humans with sporozoites
   - Poisson-distributed bites per human
   - Each bite has `transmission_rate` probability of success

## Lifecycle Mode

By default, adults emerge directly (`enable_lifecycle=False`). For full lifecycle:

```python
config.enable_lifecycle = True
```

This enables:
- Adults lay eggs based on feeding cycle
- Eggs hatch to larvae after `egg_hatch_duration` days
- Larvae develop (temperature-dependent) to immature
- Immature emerge as adults

## Example: Comparing EIP at Different Temperatures

```python
from emodlib.vector import VectorConfig, VectorPopulation

config = VectorConfig()
config.carrying_capacity = 1000
config.initial_infected_fraction = 0.1

for temp in [20, 25, 30]:
    pop = VectorPopulation()
    pop.initialize(config, seed=42)

    # Track when infectious mosquitoes appear
    for day in range(40):
        pop.update(1.0, temp, [])
        if pop.infectious_count > 10:
            print(f"T={temp}C: First infectious at day {day}")
            break
```

## References

- Styer LM et al. (2007) "Mosquito survival and population growth"
- Detinova TS (1962) "Age-grouping methods in Diptera"
- EMOD documentation: https://docs.idmod.org/projects/emod-malaria/
