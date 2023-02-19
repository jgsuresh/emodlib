import os

import yaml

from emodlib.malaria import IntrahostComponent


def describe(c, t=None):
    s = 't=%d: ' % t if t is not None else ''
    s += '(asexual, mature gametocyte, fever) = (%0.2f, %0.3f, %0.1f)' % \
        (c.parasite_density, c.gametocyte_density, c.fever_temperature)
    print(s)


def params_from_test_file():

    with open(os.path.join(os.path.realpath(os.path.dirname(__file__)), 'config.yaml')) as cfg:
        params = yaml.load(cfg, Loader=yaml.FullLoader)

    return params


def test_bindings():
    
    print('Model parameters...\n')
    params = params_from_test_file()
    print(yaml.dump(params))

    print('Configure...')
    IntrahostComponent.configure(params)

    print('Create...')
    ic = IntrahostComponent.create()

    print('Challenge...')
    ic.challenge()
    
    print('Update...')
    for t in range(30):
        ic.update(dt=1)
        describe(ic, t)

    assert ic.parasite_density > 0
    assert ic.gametocyte_density > 0

    print('Treat...')
    ic.treat()
    describe(ic)

    assert ic.parasite_density == 0
    assert ic.gametocyte_density == 0


if __name__ == '__main__':

    test_bindings()
