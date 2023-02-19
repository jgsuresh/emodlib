import yaml

from emodlib.malaria import IntrahostComponent

from test_bindings import params_from_test_file


def test_max_infections():

    print('Model parameters...\n')
    params = params_from_test_file()
    print(yaml.dump(params))

    print('Configure...')
    IntrahostComponent.configure(params)

    print('Create...')
    ic = IntrahostComponent.create()

    print('Challenge...')
    for i in range(10):
        ic.challenge()
        print('%d infections after %d challenges' % (ic.n_infections, i))
        assert ic.n_infections <= params['Max_Individual_Infections']


if __name__ == '__main__':
    test_max_infections()
