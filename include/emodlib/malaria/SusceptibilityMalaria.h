/**
 * @file SusceptibilityMalaria.h
 *
 * @brief Malaria susceptibility interface
 */

#pragma once

#include "emodlib/ParamSet.h"

#include "MalariaEnums.h"
//#include "MalariaContexts.h"
//#include "IMalariaAntibody.h"

namespace emodlib
{

    namespace malaria
    {

        class Susceptibility
        {

        public:
            struct params
            {
                // Used by MalariaAntibody for Decay and Update functions
                static float memory_level;
                static float hyperimmune_decay_rate;
                static float MSP1_antibody_growthrate;
                static float antibody_stimulation_c50;
                static float antibody_capacity_growthrate;
                static float minimum_adapted_response;
                static float non_specific_growth;
                static float antibody_csp_decay_days;
                
                // Used only by Susceptibility
                // TODO: expose demographic initialization variables to Python layer and remove from C++
                static bool enable_maternal_antibodies_transmission;
                static MaternalAntibodiesType::Enum maternal_antibodies_type; // <--
                static float maternal_antibody_protection;
                static float maternal_antibody_decay_rate;
                static InnateImmuneVariationType::Enum innate_immune_variation_type; // <--
                static float base_gametocyte_mosquito_survival;
                static float cytokine_gametocyte_inactivation;
                static float erythropoiesis_anemia_effect;
                static float pyrogenic_threshold;
                static float fever_IRBC_killrate;
                
                static void Configure(const ParamSet& pset);
            };
            
            Susceptibility();

            void Update();

            float GetFeverTemperature() const;

        private:

            float fever_temperature;
        };

    }

}
