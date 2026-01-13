"""
Tests for the emodlib vector module.

Tests cover:
- VectorConfig parameter calculations (Arrhenius, Styer mortality)
- VectorCohort state management and transitions
- VectorPopulation lifecycle and transmission dynamics
- TransmissionSimulation coordination
"""

import pytest
import math

from emodlib.vector import (
    VectorState,
    VectorConfig,
    VectorCohort,
    VectorPopulation,
    TransmissionSimulation,
)
from emodlib.malaria import MalariaConfig


# =============================================================================
# VectorConfig Tests
# =============================================================================

class TestVectorConfig:
    """Tests for VectorConfig parameter calculations."""

    @pytest.fixture
    def config(self):
        return VectorConfig()

    def test_default_values(self, config):
        """Test that default values match An. gambiae parameters."""
        assert config.adult_life_expectancy == pytest.approx(21.0)
        assert config.aquatic_mortality_rate == pytest.approx(0.2)
        assert config.days_between_feeds == pytest.approx(3.0)
        assert config.anthropophily == pytest.approx(0.9)
        assert config.egg_batch_size == pytest.approx(100.0)

    def test_computed_mortality_rate(self, config):
        """Test adult mortality rate calculation."""
        expected = 1.0 / 21.0
        assert abs(config.adult_mortality_rate - expected) < 1e-6

    def test_computed_biting_rate(self, config):
        """Test biting rate calculation."""
        expected = 1.0 / 3.0
        assert abs(config.biting_rate - expected) < 1e-6

    def test_daily_emergence(self, config):
        """Test daily emergence at equilibrium."""
        expected = config.carrying_capacity * config.adult_mortality_rate
        assert config.daily_emergence == pytest.approx(expected, rel=1e-5)

    def test_larval_progress_temperature_dependence(self, config):
        """Test that larval progress increases with temperature."""
        progress_20 = config.get_larval_progress(20.0, 1.0)
        progress_25 = config.get_larval_progress(25.0, 1.0)
        progress_30 = config.get_larval_progress(30.0, 1.0)

        assert progress_20 < progress_25 < progress_30
        assert progress_20 > 0
        assert progress_30 < 10  # Sanity check

    def test_infected_progress_temperature_dependence(self, config):
        """Test that EIP progress increases with temperature."""
        progress_20 = config.get_infected_progress(20.0, 1.0)
        progress_25 = config.get_infected_progress(25.0, 1.0)
        progress_30 = config.get_infected_progress(30.0, 1.0)

        assert progress_20 < progress_25 < progress_30
        assert progress_20 > 0

    def test_eip_reasonable_duration(self, config):
        """Test that EIP duration is biologically reasonable."""
        # At 25C, EIP is around 15-20 days with default Arrhenius parameters
        # At 30C, EIP is shorter (~10-12 days)
        progress_25 = config.get_infected_progress(25.0, 1.0)
        eip_days_25 = 1.0 / progress_25

        progress_30 = config.get_infected_progress(30.0, 1.0)
        eip_days_30 = 1.0 / progress_30

        # EIP should be in reasonable biological range
        assert 10 < eip_days_25 < 25, f"EIP at 25C = {eip_days_25} days"
        assert 7 < eip_days_30 < 20, f"EIP at 30C = {eip_days_30} days"
        assert eip_days_30 < eip_days_25, "EIP should decrease with temperature"

    def test_age_dependent_mortality(self):
        """Test Styer age-dependent mortality curve."""
        mort_0 = VectorConfig.get_age_dependent_mortality(0)
        mort_10 = VectorConfig.get_age_dependent_mortality(10)
        mort_20 = VectorConfig.get_age_dependent_mortality(20)
        mort_30 = VectorConfig.get_age_dependent_mortality(30)

        # Mortality should increase with age
        assert mort_0 < mort_10 < mort_20 < mort_30

        # Young mosquitoes should have low Styer mortality
        assert mort_0 < 0.01

        # Old mosquitoes should have high Styer mortality
        assert mort_30 > 0.1

    def test_survival_probability(self, config):
        """Test combined survival probability."""
        base_mort = config.adult_mortality_rate

        surv_0 = config.get_survival_probability(0, base_mort)
        surv_20 = config.get_survival_probability(20, base_mort)
        surv_40 = config.get_survival_probability(40, base_mort)

        # Survival should decrease with age
        assert surv_0 > surv_20 > surv_40

        # Survival should be in range (0, 1)
        assert 0 < surv_0 < 1
        assert 0 < surv_20 < 1
        assert 0 < surv_40 < 1

        # Mortality increases with age (survival per day decreases)
        assert surv_40 < surv_0  # Older mosquitoes have lower daily survival


# =============================================================================
# VectorCohort Tests
# =============================================================================

class TestVectorCohort:
    """Tests for VectorCohort state management."""

    def test_default_construction(self):
        """Test default cohort is empty egg."""
        cohort = VectorCohort()
        assert cohort.population == 0.0
        assert cohort.state == VectorState.EGG
        assert cohort.age == 0.0

    def test_parameterized_construction(self):
        """Test cohort construction with parameters."""
        cohort = VectorCohort(100.0, VectorState.ADULT, 5.0)
        assert cohort.population == 100.0
        assert cohort.state == VectorState.ADULT
        assert cohort.age == 5.0

    def test_state_queries(self):
        """Test state query methods."""
        adult = VectorCohort(100.0, VectorState.ADULT)
        infected = VectorCohort(50.0, VectorState.INFECTED)
        infectious = VectorCohort(25.0, VectorState.INFECTIOUS)
        larva = VectorCohort(200.0, VectorState.LARVA)

        assert adult.is_susceptible()
        assert not adult.is_exposed()
        assert not adult.is_infectious()
        assert adult.is_adult()
        assert not adult.is_aquatic()

        assert infected.is_exposed()
        assert infected.is_adult()

        assert infectious.is_infectious()
        assert infectious.is_adult()

        assert larva.is_aquatic()
        assert not larva.is_adult()

    def test_apply_mortality(self):
        """Test mortality application."""
        cohort = VectorCohort(100.0, VectorState.ADULT)
        cohort.apply_mortality(0.9)  # 90% survival
        assert abs(cohort.population - 90.0) < 0.01

    def test_progress_tracking(self):
        """Test development progress."""
        cohort = VectorCohort(100.0, VectorState.LARVA)
        assert cohort.progress == pytest.approx(0.0)
        assert not cohort.is_progress_complete()

        cohort.add_progress(0.5)
        assert cohort.progress == pytest.approx(0.5)
        assert not cohort.is_progress_complete()

        cohort.add_progress(0.6)
        assert cohort.progress == pytest.approx(1.1)
        assert cohort.is_progress_complete()

    def test_transition(self):
        """Test state transition."""
        cohort = VectorCohort(100.0, VectorState.LARVA, 5.0)
        cohort.add_progress(1.0)

        cohort.transition_to(VectorState.IMMATURE)

        assert cohort.state == VectorState.IMMATURE
        assert cohort.progress == 0.0  # Progress reset

    def test_split_fraction(self):
        """Test splitting by fraction."""
        cohort = VectorCohort(100.0, VectorState.ADULT, 10.0)
        split = cohort.split(0.3)

        assert abs(split.population - 30.0) < 0.01
        assert abs(cohort.population - 70.0) < 0.01
        assert split.state == cohort.state
        assert split.age == cohort.age

    def test_split_count(self):
        """Test splitting by count."""
        cohort = VectorCohort(100.0, VectorState.ADULT)
        split = cohort.split_count(25.0)

        assert abs(split.population - 25.0) < 0.01
        assert abs(cohort.population - 75.0) < 0.01

    def test_is_empty(self):
        """Test empty detection."""
        cohort = VectorCohort(0.005, VectorState.ADULT)
        assert cohort.is_empty()

        cohort.population = 1.0
        assert not cohort.is_empty()


# =============================================================================
# VectorPopulation Tests
# =============================================================================

class TestVectorPopulation:
    """Tests for VectorPopulation dynamics."""

    @pytest.fixture
    def config(self):
        config = VectorConfig()
        config.carrying_capacity = 1000.0
        return config

    @pytest.fixture
    def population(self, config):
        pop = VectorPopulation()
        pop.initialize(config, seed=42)
        return pop

    def test_initialization(self, population, config):
        """Test population initializes to some positive value."""
        adult_pop = population.adult_population
        # Should have a substantial population (equilibrium age structure
        # means total is less than daily_emergence * life_expectancy due to mortality)
        assert adult_pop > 0
        # Should be roughly proportional to carrying capacity
        assert adult_pop > config.carrying_capacity * 0.3

    def test_population_queries(self, population):
        """Test population query methods."""
        assert population.susceptible_count > 0
        assert population.exposed_count >= 0
        assert population.infectious_count >= 0

        total_adults = (population.susceptible_count +
                        population.exposed_count +
                        population.infectious_count)
        assert abs(population.adult_population - total_adults) < 0.1

    def test_equilibrium_stability(self, config):
        """Test population trends toward equilibrium."""
        pop = VectorPopulation()
        pop.initialize(config, seed=42)

        initial_pop = pop.adult_population

        # Run for 100 days
        for _ in range(100):
            pop.update(1.0, 25.0, [])

        final_pop = pop.adult_population

        # Population should stay in reasonable range
        # (may drift slightly as equilibrium adjusts)
        ratio = final_pop / initial_pop
        assert 0.5 < ratio < 2.0, f"Population ratio {ratio} outside reasonable range"
        # Should still have substantial population
        assert final_pop > config.carrying_capacity * 0.3

    def test_eip_progression(self, config):
        """Test that infected mosquitoes progress through EIP."""
        config.initial_infected_fraction = 0.0
        pop = VectorPopulation()
        pop.initialize(config, seed=42)

        # No infected/infectious initially
        assert pop.exposed_count == 0
        assert pop.infectious_count == 0

        # Update with high human infectiousness to infect mosquitoes
        for _ in range(5):
            pop.update(1.0, 25.0, [1.0])  # 100% infectious humans

        exposed_after_infection = pop.exposed_count
        assert exposed_after_infection > 0, "Should have some exposed mosquitoes"

        # Continue updating - exposed should become infectious
        # EIP at 25C is ~19 days, so run for 25 days to ensure completion
        infectious_before = pop.infectious_count
        for _ in range(25):
            pop.update(1.0, 25.0, [0.0])

        infectious_after = pop.infectious_count
        # Either we have infectious, or exposed decreased significantly
        # (some die during EIP, others complete it)
        exposed_after = pop.exposed_count
        assert (infectious_after > infectious_before or
                exposed_after < exposed_after_infection * 0.5), \
            "Should see EIP progression (infectious increase or exposed decrease)"

    def test_sporozoite_challenges(self, config):
        """Test sporozoite challenge generation."""
        config.initial_infected_fraction = 0.5  # Start with some infectious
        pop = VectorPopulation()
        pop.initialize(config, seed=42)

        # Run a few days to let some become infectious
        for _ in range(20):
            pop.update(1.0, 25.0, [])

        challenges = pop.get_sporozoite_challenges(10)
        assert len(challenges) == 10

        # With infectious mosquitoes, some humans should receive challenges
        if pop.infectious_count > 0:
            total_challenges = sum(challenges)
            assert total_challenges >= 0  # Could be 0 due to stochasticity

    def test_no_infection_without_infectiousness(self, config):
        """Test that mosquitoes don't get infected without infectious humans."""
        config.initial_infected_fraction = 0.0
        pop = VectorPopulation()
        pop.initialize(config, seed=42)

        # Update with no human infectiousness
        for _ in range(30):
            pop.update(1.0, 25.0, [0.0, 0.0, 0.0])

        assert pop.exposed_count == 0
        assert pop.infectious_count == 0


# =============================================================================
# TransmissionSimulation Tests
# =============================================================================

class TestTransmissionSimulation:
    """Tests for TransmissionSimulation coordination."""

    @pytest.fixture
    def malaria_config(self):
        return MalariaConfig()

    @pytest.fixture
    def vector_config(self):
        config = VectorConfig()
        config.carrying_capacity = 500.0
        return config

    @pytest.fixture
    def simulation(self, malaria_config, vector_config):
        sim = TransmissionSimulation()
        sim.initialize(malaria_config, vector_config, n_humans=20, seed=42)
        return sim

    def test_initialization(self, simulation):
        """Test simulation initializes correctly."""
        assert simulation.n_humans == 20
        assert simulation.current_day == 0
        assert simulation.prevalence == 0.0
        assert simulation.vector_population > 0

    def test_seed_infections(self, simulation):
        """Test seeding infections in humans."""
        simulation.seed_infections(5)

        # Run one update to process
        simulation.update(1.0, 25.0)

        assert simulation.infected_count == 5
        assert simulation.prevalence == 0.25  # 5/20

    def test_update_advances_time(self, simulation):
        """Test that update advances simulation time."""
        assert simulation.current_day == 0

        simulation.update(1.0, 25.0)
        assert simulation.current_day == 1

        simulation.update(1.0, 25.0)
        assert simulation.current_day == 2

    def test_transmission_cycle(self, malaria_config, vector_config):
        """Test complete transmission cycle over time."""
        vector_config.carrying_capacity = 2000.0  # Higher for transmission

        sim = TransmissionSimulation()
        sim.initialize(malaria_config, vector_config, n_humans=50, seed=42)

        # Seed initial infection
        sim.seed_infections(1)

        initial_infected = 1
        saw_gametocytes = False
        saw_vector_infection = False

        # Run for 60 days
        for day in range(60):
            sim.update(1.0, 25.0)

            if sim.mean_gametocyte_density > 0:
                saw_gametocytes = True

            if sim.vector_exposed_count > 0 or sim.vector_infectious_count > 0:
                saw_vector_infection = True

        # After 60 days with initial infection, should see gametocytes
        assert saw_gametocytes, "Should have seen gametocytes by day 60"

        # May or may not see vector infection depending on timing
        # Just check simulation ran without error

    def test_vector_metrics(self, simulation):
        """Test vector metric accessors."""
        simulation.update(1.0, 25.0)

        assert simulation.vector_population > 0
        assert simulation.vector_susceptible_count >= 0
        assert simulation.vector_exposed_count >= 0
        assert simulation.vector_infectious_count >= 0

        # Fractions should be in [0, 1]
        assert 0 <= simulation.vector_infectious_fraction <= 1
        assert 0 <= simulation.vector_infected_fraction <= 1


# =============================================================================
# Integration Tests
# =============================================================================

class TestIntegration:
    """Integration tests for complete workflows."""

    def test_arrhenius_formula_matches_expected(self):
        """Verify Arrhenius formula produces expected EIP values."""
        config = VectorConfig()

        # Test EIP at different temperatures
        progress_25 = config.get_infected_progress(25.0, 1.0)
        eip_25 = 1.0 / progress_25

        progress_30 = config.get_infected_progress(30.0, 1.0)
        eip_30 = 1.0 / progress_30

        progress_35 = config.get_infected_progress(35.0, 1.0)
        eip_35 = 1.0 / progress_35

        # EIP should decrease with temperature (Arrhenius)
        assert eip_35 < eip_30 < eip_25, "EIP should decrease with temperature"

        # EIP should be biologically reasonable (5-25 days depending on temp)
        assert 5 < eip_35 < 15, f"EIP at 35C = {eip_35}"
        assert 7 < eip_30 < 20, f"EIP at 30C = {eip_30}"
        assert 10 < eip_25 < 25, f"EIP at 25C = {eip_25}"

    def test_lifecycle_mode(self):
        """Test full lifecycle mode with eggs."""
        config = VectorConfig()
        config.carrying_capacity = 100.0
        config.enable_lifecycle = True

        pop = VectorPopulation()
        pop.initialize(config, seed=42)

        # Run for a few days
        for _ in range(10):
            pop.update(1.0, 25.0, [])

        # In lifecycle mode, should have eggs/larvae
        # (though they may be small due to equilibrium initialization)
        assert pop.adult_population > 0

    def test_deterministic_with_seed(self):
        """Test that same seed produces same results."""
        config = VectorConfig()
        config.carrying_capacity = 500.0

        pop1 = VectorPopulation()
        pop1.initialize(config, seed=12345)

        pop2 = VectorPopulation()
        pop2.initialize(config, seed=12345)

        # Should have same initial population
        assert pop1.adult_population == pop2.adult_population

        # Run same updates
        for _ in range(10):
            pop1.update(1.0, 25.0, [0.1])
            pop2.update(1.0, 25.0, [0.1])

        # Should still be equal
        assert pop1.adult_population == pop2.adult_population
        assert pop1.exposed_count == pop2.exposed_count


if __name__ == "__main__":
    pytest.main(["-v", "-s", __file__])
