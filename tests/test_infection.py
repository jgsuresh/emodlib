import pytest
from test_host import params_from_test_file

from emodlib.malaria import Infection, IntrahostComponent, Susceptibility


@pytest.fixture
def params():
    p = params_from_test_file()
    IntrahostComponent.configure(p)  # required for Infection::Initialize params?
    return p


@pytest.fixture
def default_susceptibility(params):
    yield Susceptibility.create()


@pytest.fixture
def default_infection(default_susceptibility):
    yield Infection.create(susceptibility=default_susceptibility, hepatocytes=1)


def test_infection(default_infection, params):
    infs = [default_infection for _ in range(10)]

    msp_types = [i.msp_type for i in infs]
    print(msp_types)
    assert all([t < params["Falciparum_MSP_Variants"] for t in msp_types])

    major_types = infs[0].pfemp1_major_types
    print(major_types)
    assert len(major_types) == 50
    assert all([t < params["Falciparum_PfEMP1_Variants"] for t in major_types])


def test_antibody(default_infection):
    inf = default_infection
    for _ in range(10):
        inf.update(dt=1)
        print("MSP antigen count = %0.2f" % inf.msp_antibody.antigen_count)
    assert inf.msp_antibody.antigen_count > 0


if __name__ == "__main__":
    pytest.main(["-vv", "-s", __file__])
