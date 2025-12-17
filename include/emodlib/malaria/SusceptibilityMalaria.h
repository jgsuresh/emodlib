/**
 * @file SusceptibilityMalaria.h
 *
 * @brief Malaria susceptibility interface
 */

#pragma once

#include <vector>

#include "MalariaEnums.h"
#include "IMalariaAntibody.h"


namespace emodlib
{

    namespace malaria
    {

        struct MalariaConfig;  // forward declaration

        class Susceptibility
        {

        public:

            static Susceptibility *Create(MalariaConfig* config);
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
            float get_cytokines() const;
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

            MalariaConfig* GetConfig() const;

        private:

            MalariaConfig* m_config;  // non-owning pointer to configuration

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


            Susceptibility(MalariaConfig* config);
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
