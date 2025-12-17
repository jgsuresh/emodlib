import os

import pytest
import yaml

from emodlib.malaria import IntrahostComponent, create_config


def describe(c, t=None):
    s = "t=%d: " % t if t is not None else ""
    s += "(asexual, mature gametocyte, fever) = ({:0.2f}, {:0.3f}, {:0.1f})".format(
        c.parasite_density, c.gametocyte_density, c.fever_temperature
    )
    print(s)


def params_from_default_file():
    with open(
        os.path.join(
            os.path.realpath(os.path.dirname(__file__)),
            "..",
            "src",
            "emodlib",
            "malaria",
            "config.yml",
        )
    ) as cfg:
        params = yaml.load(cfg, Loader=yaml.FullLoader)

    return params


def test_config():
    print("Load default model parameters...\n")
    params = params_from_default_file()

    # Test that we can create a config and modify it
    config = create_config()
    assert config.random_seed == params["Run_Number"]

    run_number = 123456
    config2 = create_config({'Run_Number': run_number})
    assert config2.random_seed == run_number

    # Test that modifying merozoites_per_schizont works
    config3 = create_config({
        'Run_Number': run_number,
        'infection_params': {'Merozoites_Per_Schizont': 10}
    })
    assert config3.random_seed == run_number
    assert config3.merozoites_per_schizont == 10


def test_bindings():
    print("Create config...")
    config = create_config()

    print("Create...")
    ic = IntrahostComponent.create(config)

    print("Challenge...")
    ic.challenge()

    print("Update...")
    for t in range(35):
        ic.update(dt=1)
        describe(ic, t)

        msp = ic.infections[0].msp_antibody
        print(
            "\tanti-MSP (capacity, concentration): %0.2f, %0.2f"
            % (msp.antibody_capacity, msp.antibody_concentration)
        )

    assert ic.parasite_density > 0
    assert ic.gametocyte_density > 0
    assert ic.infectiousness > 0

    assert ic.infections[0].msp_antibody.antibody_capacity > 0
    assert ic.infections[0].msp_antibody.antibody_concentration > 0

    print("Treat...")
    ic.treat()
    describe(ic)

    assert ic.parasite_density == 0
    assert ic.gametocyte_density == 0


def test_infection_clearance():
    print("Create config...")
    config = create_config()

    print("Create...")
    ic = IntrahostComponent.create(config)

    print("Update...")
    n_infections_cache = 0
    n_cleared = 0
    for t in range(365 * 5):
        if t % 180 == 0:
            ic.challenge()
            print("%d infections after t=%d challenge" % (ic.n_infections, t))
        ic.update(dt=1)
        n_infections = ic.n_infections
        if n_infections_cache != n_infections:
            if n_infections < n_infections_cache:
                n_cleared += 1
            n_infections_cache = n_infections
            print("%d infections at t=%d" % (n_infections, t))

    print("%d total cleared infections" % n_cleared)
    assert n_cleared > 0


def test_max_infections():
    print("Load default model parameters...\n")
    params = params_from_default_file()

    print("Create config...")
    config = create_config()

    print("Create...")
    ic = IntrahostComponent.create(config)

    print("Challenge...")
    for i in range(10):
        ic.challenge()
        print("%d infections after %d challenges" % (ic.n_infections, i))
        assert ic.n_infections <= params["Max_Individual_Infections"]


if __name__ == "__main__":
    pytest.main(["-vv", "-s", __file__])
