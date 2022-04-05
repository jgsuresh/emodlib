import yaml

from emodlib.malaria import IntrahostComponent

from test_bindings import params

if __name__ == '__main__':

    print('Model parameters...\n')
    print(yaml.dump(params))

    print('Configure...')
    IntrahostComponent.configure(params)

    print('Create...')
    ic = IntrahostComponent.create()

    print('Challenge...')
    for i in range(10):
        ic.challenge()
        print('%d infections after %d challenges' % (ic.n_infections, i))