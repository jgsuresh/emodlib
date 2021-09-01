/**
 * @file InfectionMalaria.cpp
 *
 * @brief Malaria infection implementation
 */

#include "InfectionMalaria.h"


namespace emodlib
{

    namespace malaria
    {

        ParasiteSwitchType::Enum Infection::params::parasite_switch_type = ParasiteSwitchType::RATE_PER_PARASITE_7VARS;
        MalariaStrains::Enum     Infection::params::malaria_strains = MalariaStrains::FALCIPARUM_RANDOM_STRAIN;

        float Infection::params::antibody_IRBC_killrate = 0.0f;
        float Infection::params::MSP1_merozoite_kill = 0.0f;
        float Infection::params::gametocyte_stage_survival = 0.0f;
        float Infection::params::base_gametocyte_sexratio = 0.0f;
        float Infection::params::base_gametocyte_production = 0.0f;
        float Infection::params::antigen_switch_rate = 0.0f;
        float Infection::params::merozoites_per_hepatocyte = 0.0f;
        float Infection::params::merozoites_per_schizont = 0.0f;
        float Infection::params::non_specific_antigenicity = 0.0f;
        float Infection::params::RBC_destruction_multiplier = 0.0f;
        int   Infection::params::n_asexual_cycles_wo_gametocytes = 0;
    
        void Infection::params::Configure(const ParamSet& pset)
        {

        }
    
        Infection::Infection()
            : parasite_density(0)
            , gametocyte_density(0)
        {

        }

        void Infection::Update()
        {

        }

        float Infection::GetParasiteDensity() const
        {
            return parasite_density;
        }

        float Infection::GetGametocyteDensity() const
        {
            return gametocyte_density;
        }

    }

}
