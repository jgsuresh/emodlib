import yaml

from emodlib.malaria import IntrahostComponent

from test_bindings import params_from_test_file


def test_infection_clearance():

    print('Model parameters...\n')
    params = params_from_test_file()
    print(yaml.dump(params))

    print('Configure...')
    IntrahostComponent.configure(params)

    print('Create...')
    ic = IntrahostComponent.create()

    print('Update...')
    n_infections_cache = 0
    n_cleared = 0
    for t in range(365*5):
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

    print('%d total cleared infections' % n_cleared)
    assert n_cleared > 0


if __name__ == '__main__':

    test_infection_clearance()