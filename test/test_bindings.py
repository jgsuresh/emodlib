import yaml

from emodlib.malaria import *


print('Configuring...\n')

params = dict(
    infection_params=dict(
        increment_parasite=12.3,
        increment_gametocyte=5.6
    ),
    susceptibility_params=dict(
        increment_fever=0.12
    )
)

print(yaml.dump(params))

IntrahostComponent.configure(params)


def describe(c):
    print('(p,g,f) = (%0.2f, %0.2f, %0.2f)' %
        (c.parasite_density, c.gametocyte_density, c.fever_temperature))


print('Create...')
ic = IntrahostComponent.create()
describe(ic)

print('Update...')
ic.update()
describe(ic)

print('Challenge + update...')
ic.challenge()
ic.update()
describe(ic)

print('Treat...')
ic.treat()
describe(ic)
