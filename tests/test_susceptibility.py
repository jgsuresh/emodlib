import pytest

from emodlib.malaria import IntrahostComponent, Susceptibility


def test_aging():
    s = Susceptibility.create()

    s.age = 30
    print("immune system age = %d days" % s.age)

    s.update(dt=10)
    print("immune system age = %d days" % s.age)

    assert s.age == 40


def test_maternal_antibodies():
    s = Susceptibility.create()

    s.maternal_antibody_strength = 0.8

    for t in range(5):
        s.update(dt=7)
        print(
            "maternal antibody strength remaining = %0.2f"
            % s.maternal_antibody_strength
        )

    assert s.maternal_antibody_strength < 0.8


def test_immune_init():
    print("Set default parameters...")
    IntrahostComponent.set_params()

    ic = IntrahostComponent.create()
    s = ic.susceptibility

    print("immune system age = %d days" % s.age)
    assert s.age == 20 * 365

    print("maternal antibodies strength = %0.2f" % s.maternal_antibody_strength)
    assert s.maternal_antibody_strength == 0

    print("pyrogenic threshold = %d" % s.pyrogenic_threshold)
    assert s.pyrogenic_threshold == 15000

    print("fever_kill_rate = %0.2f" % s.fever_kill_rate)
    assert pytest.approx(s.fever_kill_rate, 1e-6) == 1.4

    # ------

    s.age = 0  # reset immune age
    s.pyrogenic_threshold = 30000

    ic.update(dt=10)  # update the individual

    assert ic.susceptibility.age == 10
    assert ic.susceptibility.pyrogenic_threshold == 30000


if __name__ == "__main__":
    pytest.main(["-vv", "-s", __file__])
