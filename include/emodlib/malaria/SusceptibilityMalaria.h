/**
 * @file SusceptibilityMalaria.h
 *
 * @brief Malaria susceptibility interface
 */

#pragma once

#include "emodlib/ParamSet.h"

#include "MalariaEnums.h"
#include "IMalariaAntibody.h"


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
            
            
            static Susceptibility *Create();
            IMalariaAntibody* RegisterAntibody(MalariaAntibodyType::Enum type, int variant, float capacity=0.0f);
            
            void Update();

            float GetFeverTemperature() const;

            
        private:

            float age;  // TODO: access from DemographicComponent * -- was in base class as parallel variable
            
            // containers for antibody objects
            int32_t m_antigenic_flag;
            float m_maternal_antibody_strength;
            IMalariaAntibody* m_CSP_antibody;
            std::vector<IMalariaAntibody*> m_active_MSP_antibodies;
            std::vector<IMalariaAntibody*> m_active_PfEMP1_minor_antibodies;
            std::vector<IMalariaAntibody*> m_active_PfEMP1_major_antibodies;

            // RBC information
            int64_t m_RBC;
            int64_t m_RBCcapacity;
            int64_t m_RBCproduction;   // how many RBC's a person should have /120 (AVERAGE_RBC_LIFESPAN)
            float   m_inv_microliters_blood; // ==/(age dependent estimate of blood volume)

            // symptomatic variables
            float m_cytokines;
            float m_ind_pyrogenic_threshold;
            float m_ind_fever_kill_rate;
            float m_cytokine_stimulation;
            float m_parasite_density;
            
            float fever_temperature;
            
            
            Susceptibility();
            
            // TODO: expose heterogeneity initialization + interaction with DemographicComponent via Python layer
            void Initialize();  // <-- original base class: Initialize(age, immmod, riskmod)
            void  recalculateBloodCapacity( float _age );
            
        };

    }

}
