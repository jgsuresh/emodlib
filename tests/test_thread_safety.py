"""
Thread-safety tests for the malaria intrahost model.

These tests verify that multiple simulations can run in parallel with independent
configurations without interfering with each other.
"""

import concurrent.futures

import pytest

from emodlib.malaria import IntrahostComponent, create_config


def run_simulation(seed, n_steps=100):
    """Run a simulation with the given seed and return results."""
    config = create_config({'Run_Number': seed})
    ic = IntrahostComponent.create(config)
    ic.challenge()

    results = []
    for _ in range(n_steps):
        ic.update(dt=1)
        results.append(ic.parasite_density)

    return results


def test_parallel_simulations_are_deterministic():
    """
    Run same simulation twice sequentially - should get identical results
    since they use the same seed.
    """
    results1 = run_simulation(seed=42)
    results2 = run_simulation(seed=42)

    assert results1 == results2, "Sequential runs with same seed should be deterministic"


def test_parallel_simulations_are_independent():
    """
    Run same simulation twice in parallel with same seed.
    With instance-based config, both should get identical results.
    """
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
        future1 = executor.submit(run_simulation, seed=42)
        future2 = executor.submit(run_simulation, seed=42)
        results1 = future1.result()
        results2 = future2.result()

    assert results1 == results2, "Parallel runs with same seed should produce identical results"


def test_different_seeds_produce_different_results():
    """Different seeds should produce different results."""
    results1 = run_simulation(seed=1)
    results2 = run_simulation(seed=2)

    assert results1 != results2, "Different seeds should produce different results"


def test_different_configs_are_isolated():
    """
    Different configs should not affect each other.
    This is the key thread-safety test.
    """
    config1 = create_config({'Run_Number': 1, 'Max_Individual_Infections': 3})
    config2 = create_config({'Run_Number': 2, 'Max_Individual_Infections': 10})

    ic1 = IntrahostComponent.create(config1)
    ic2 = IntrahostComponent.create(config2)

    # Challenge both many times
    for _ in range(20):
        ic1.challenge()
        ic2.challenge()

    # Verify each uses its own max_ind_inf config
    assert ic1.n_infections <= 3, f"ic1 should have at most 3 infections, got {ic1.n_infections}"
    assert ic2.n_infections <= 10, f"ic2 should have at most 10 infections, got {ic2.n_infections}"

    # Since ic2 has a higher limit and same number of challenges, it should have more infections
    # (probabilistically, though not guaranteed - let's just check the bounds)
    print(f"ic1 infections: {ic1.n_infections}, ic2 infections: {ic2.n_infections}")


def test_multiple_parallel_simulations():
    """Run many simulations in parallel."""
    num_simulations = 10
    seeds = list(range(num_simulations))

    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
        futures = [executor.submit(run_simulation, seed=s, n_steps=50) for s in seeds]
        results = [f.result() for f in futures]

    # All results should be different (different seeds)
    for i, r1 in enumerate(results):
        for j, r2 in enumerate(results):
            if i != j:
                assert r1 != r2, f"Results {i} and {j} should be different"

    print(f"Successfully ran {num_simulations} parallel simulations")


if __name__ == "__main__":
    pytest.main(["-vv", "-s", __file__])
