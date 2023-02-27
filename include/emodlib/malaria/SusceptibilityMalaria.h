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
                // Used in MalariaAntibody for boost-decay functions
                static float memory_level;
                static float hyperimmune_decay_rate;
                static float MSP1_antibody_growthrate;
                static float antibody_stimulation_c50;
                static float antibody_capacity_growthrate;
                static float minimum_adapted_response;
                static float non_specific_growth;
                static float antibody_csp_decay_days;

                // Used in Susceptibility for:
                // ...maternal protection
                // static bool enable_maternal_antibodies_transmission;
                // static MaternalAntibodiesType::Enum maternal_antibodies_type; // TODO: emodlib#9 (innate init)
                // static float maternal_antibody_protection;
                static float maternal_antibody_decay_rate;

                // ...innate immunity effects
                // static InnateImmuneVariationType::Enum innate_immune_variation_type; // TODO: emodlib#9 (innate init)
                static float pyrogenic_threshold;
                static float fever_IRBC_killrate;

                // ... infectiousness calculations
                // static float base_gametocyte_mosquito_survival;  // TODO: emodlib#7 (infectiousness calculations)
                // static float cytokine_gametocyte_inactivation;

                // ... red blood cell effects
                static float erythropoiesis_anemia_effect;

                static void Configure(const ParamSet& pset);
            };


            static Susceptibility *Create();
            IMalariaAntibody* RegisterAntibody(MalariaAntibodyType::Enum type, int variant, float capacity=0.0f);
            void UpdateActiveAntibody( pfemp1_antibody_t &pfemp1_variant, int minor_variant, int major_variant );
            void remove_RBCs(int64_t infectedAsexual, int64_t infectedGametocytes, double RBC_destruction_multiplier);

            void Update(float dt);
            void SetAntigenPresent();

            long long get_RBC_count() const;
            float get_inv_microliters_blood() const;
            double get_RBC_availability() const;
            float get_fever() const;
            float get_fever_celsius() const;
            float get_fever_killing_rate() const;
            float get_parasite_density() const;
            float get_maternal_antibodies() const;

            float get_age() const;
            void set_age(float _age);

            float get_maternal_antibody_strength() const;
            void set_maternal_antibody_strength(float _matAb);

            float get_pyrogenic_threshold() const;
            void set_pyrogenic_threshold(float _threshold);

            float get_fever_kill_rate() const;
            void set_fever_kill_rate(float _rate);

        private:

            float age;  // TODO: emodlib#10 (demographic components)

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


            Susceptibility();
            void Initialize();  // TODO: emodlib#9 (innate init) + emodlib#10 (demographic/transmission components)

            void recalculateBloodCapacity( float _age );
            void updateImmunityCSP( float dt );
            void updateImmunityMSP( float dt, float& temp_cytokine_stimulation );
            void updateImmunityPfEMP1Minor( float dt );
            void updateImmunityPfEMP1Major( float dt );
            void decayAllAntibodies( float dt );

        };

    }

}
