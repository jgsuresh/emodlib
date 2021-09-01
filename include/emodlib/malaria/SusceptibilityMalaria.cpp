/**
 * @file SusceptibilityMalaria.cpp
 *
 * @brief Malaria susceptibility implementation
 */

#include "SusceptibilityMalaria.h"


namespace emodlib
{

    namespace malaria
    {
    
        float  Susceptibility::params::memory_level                      = 0.0f;
        float  Susceptibility::params::hyperimmune_decay_rate            = 0.0f;
        float  Susceptibility::params::MSP1_antibody_growthrate          = 0.0f;
        float  Susceptibility::params::antibody_stimulation_c50          = 0.0f;
        float  Susceptibility::params::antibody_capacity_growthrate      = 0.0f;
        float  Susceptibility::params::minimum_adapted_response          = 0.0f;
        float  Susceptibility::params::non_specific_growth               = 0.0f;
        float  Susceptibility::params::antibody_csp_decay_days           = 0.0f;
    
        bool   Susceptibility::params::enable_maternal_antibodies_transmission  = false;
        MaternalAntibodiesType::Enum Susceptibility::params::maternal_antibodies_type = MaternalAntibodiesType::OFF;
        float  Susceptibility::params::maternal_antibody_protection      = 0.0f;
        float  Susceptibility::params::maternal_antibody_decay_rate      = 0.0f;
        InnateImmuneVariationType::Enum Susceptibility::params::innate_immune_variation_type = InnateImmuneVariationType::NONE;
        float  Susceptibility::params::base_gametocyte_mosquito_survival = 1.0f;
        float  Susceptibility::params::cytokine_gametocyte_inactivation  = 1.0f;
        float  Susceptibility::params::erythropoiesis_anemia_effect      = 0.0f;
        float  Susceptibility::params::pyrogenic_threshold               = 0.0f;
        float  Susceptibility::params::fever_IRBC_killrate               = 0.0f;
    
        void Susceptibility::params::Configure(const ParamSet& pset)
        {

        }
    
        Susceptibility::Susceptibility()
            : fever_temperature(37.0f)
        {

        }

        void Susceptibility::Update()
        {

        }

        float Susceptibility::GetFeverTemperature() const
        {
            return fever_temperature;
        }

    }

}
