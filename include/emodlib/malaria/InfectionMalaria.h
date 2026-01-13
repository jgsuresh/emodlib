/**
 * @file InfectionMalaria.h
 *
 * @brief Malaria infection interface
 */

#pragma once

#include <vector>

#include "emodlib/utils/suids.hpp"

#include "Malaria.h"
#include "MalariaEnums.h"
#include "IMalariaAntibody.h"


namespace emodlib
{

    namespace malaria
    {

        struct MalariaConfig;  // forward declaration
        class Susceptibility;

        // Debug structure for forcing specific antigens (for Python comparison)
        struct DebugAntigens {
            int32_t msp_type;
            int32_t nonspec_type;
            std::vector<int32_t> irbc_types;        // size = CLONAL_PfEMP1_VARIANTS (50)
            std::vector<int32_t> minor_epitope_types; // size = CLONAL_PfEMP1_VARIANTS (50)
        };

        class Infection
        {

        public:

            static Infection *Create(Susceptibility* _susceptibility, MalariaConfig* config, int initial_hepatocytes=1);
            // Debug version that accepts explicit antigens
            static Infection *CreateWithAntigens(Susceptibility* _susceptibility, MalariaConfig* config,
                                                  const DebugAntigens& antigens, int initial_hepatocytes=1);

            void Update(float dt);

            suids::suid GetSuid() const;
            int64_t get_MaleGametocytes(int stage) const;
            int64_t get_FemaleGametocytes(int stage) const;
            float get_asexual_density() const;
            float get_mature_gametocyte_density() const;
            bool IsCleared() const;

            int32_t get_msp_type() const;
            std::vector<int32_t> get_pfemp1_major_types() const;

            IMalariaAntibody* get_msp_antibody() const;

            // Additional getters for detailed state exposure (for validation)
            std::vector<int64_t> get_irbc_counts() const;
            int32_t get_hepatocyte_count() const;
            float get_liver_stage_timer() const;
            int32_t get_asexual_phase() const;
            double get_asexual_cycle_timer() const;
            int32_t get_asexual_cycle_count() const;
            std::vector<int64_t> get_all_male_gametocytes() const;
            std::vector<int64_t> get_all_female_gametocytes() const;
            std::vector<int32_t> get_minor_epitope_types() const;

        private:

            suids::suid suid; // unique id of this infection within the system

            float m_liver_stage_timer;
            double m_IRBCtimer;
            int32_t m_hepatocytes;
            AsexualCycleStatus::Enum m_asexual_phase;
            int32_t m_asexual_cycle_count;

            int32_t m_MSPtype;        // allow variation in MSP from clone to clone
            int32_t m_nonspectype;    // what is the set of minor_epitope_types
            int32_t m_minor_epitope_type[CLONAL_PfEMP1_VARIANTS];
            int32_t m_IRBCtype[CLONAL_PfEMP1_VARIANTS];

            IMalariaAntibody* m_MSP_antibody;
            std::vector< pfemp1_antibody_t > m_PfEMP1_antibodies;

            std::vector<int64_t> m_IRBC_count;
            int64_t m_malegametocytes[GametocyteStages::Count];
            int64_t m_femalegametocytes[GametocyteStages::Count];

            // placeholders for infection-level variation in merozoite-to-gametocyte dynamics
            double m_gametorate;
            double m_gametosexratio;

            Susceptibility* immunity;
            MalariaConfig* m_config;  // non-owning pointer to configuration

            // EMOD 2.22: Time tracking for antibody decay
            float m_current_time;

            Infection(MalariaConfig* config);
            void Initialize(Susceptibility* _susceptibility, int initial_hepatocytes);
            void InitializeWithAntigens(Susceptibility* _susceptibility, const DebugAntigens& antigens, int initial_hepatocytes);

            void malariaProcessHepatocytes(float dt);
            void processEndOfAsexualCycle(float dt);
            void malariaIRBCAntigenSwitch(double merozoitesurvival = 1.0);
            void malariaCycleGametocytes(double merozoitesurvival = 1.0);
            void malariaImmuneStimulation(float dt);
            void malariaImmunityIRBCKill(float dt);
            void malariaImmunityGametocyteKill(float dt);
            void malariaCheckInfectionStatus(float dt);  // TODO: emodlib#3 (InfectionStateChange::Cleared)
            void apply_MatureGametocyteKillProbability(float pkill);

        };

    }

}
