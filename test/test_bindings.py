import yaml

from emodlib.malaria import *


def describe(c, t=None):
    s = 't=%d: ' % t if t is not None else ''
    s += '(asexual, mature gametocyte, fever) = (%0.2f, %0.3f, %0.1f)' % \
        (c.parasite_density, c.gametocyte_density, c.fever_temperature)
    print(s)


params = dict(
    Run_Number=12345,
    
    Max_Individual_Infections=5,

    Falciparum_MSP_Variants=32,
    Falciparum_Nonspecific_Types=76,
    Falciparum_PfEMP1_Variants=1070,
    
    infection_params=dict(
        Base_Incubation_Period=7,
        Antibody_IRBC_Kill_Rate=1.596,
        Nonspecific_Antigenicity_Factor=0.415111634,
        MSP1_Merozoite_Kill_Fraction=0.511735322,
        Gametocyte_Stage_Survival_Rate=0.588569307,
        Base_Gametocyte_Fraction_Male=0.2,
        Base_Gametocyte_Production_Rate=0.06150582,
        Antigen_Switch_Rate=pow(10, -9.116590124),
        Merozoites_Per_Hepatocyte=15000,
        Merozoites_Per_Schizont=16,
        RBC_Destruction_Multiplier=3.29,
        Number_Of_Asexual_Cycles_Without_Gametocytes=1
    ),
    
    susceptibility_params=dict(
        Antibody_Memory_Level=0.34,
        Max_MSP1_Antibody_Growthrate=0.045,
        Antibody_Stimulation_C50=30,
        Antibody_Capacity_Growth_Rate=0.09,
        Min_Adapted_Response=0.05,
        Nonspecific_Antibody_Growth_Rate_Factor=0.5,
        Antibody_CSP_Decay_Days=90,
        Maternal_Antibody_Decay_Rate=0.01,
        Pyrogenic_Threshold=1.5e4,
        Fever_IRBC_Kill_Rate=1.4,
        Erythropoiesis_Anemia_Effect=3.5
    )
)


if __name__ == '__main__':

    print('Model parameters...\n')
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

    print('Treat...')
    ic.treat()
    describe(ic)
