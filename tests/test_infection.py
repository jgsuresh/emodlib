import pytest

from emodlib.malaria import Infection, IntrahostComponent, Susceptibility, create_config


@pytest.fixture
def config():
    return create_config()


@pytest.fixture
def default_susceptibility(config):
    yield Susceptibility.create(config)


@pytest.fixture
def default_infection(default_susceptibility, config):
    yield Infection.create(susceptibility=default_susceptibility, config=config, hepatocytes=1)


def test_infection(default_infection, config):
    infs = [default_infection for _ in range(10)]

    msp_types = [i.msp_type for i in infs]
    print(msp_types)
    assert all([t < config.falciparum_MSP_variants for t in msp_types])

    major_types = infs[0].pfemp1_major_types
    print(major_types)
    assert len(major_types) == 50
    assert all([t < config.falciparum_PfEMP1_variants for t in major_types])


def test_antibody(default_infection):
    inf = default_infection
    for _ in range(10):
        inf.update(dt=1)
        print("MSP antigen count = %0.2f" % inf.msp_antibody.antigen_count)
    assert inf.msp_antibody.antigen_count > 0


if __name__ == "__main__":
    pytest.main(["-vv", "-s", __file__])
