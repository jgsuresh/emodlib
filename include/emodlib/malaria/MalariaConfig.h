/**
 * @file MalariaConfig.h
 *
 * @brief Malaria configuration - instance-based parameters for thread-safe operation
 */

#pragma once

#include <memory>

#include "emodlib/ParamSet.h"
#include "emodlib/utils/RANDOM.h"
#include "emodlib/utils/suids.hpp"

#include "Malaria.h"


namespace emodlib
{

    namespace malaria
    {

        /**
         * @brief Configuration object that owns all malaria model parameters.
         *
         * This replaces the static params structs in IntrahostComponent, Susceptibility,
         * and Infection classes. Each IntrahostComponent instance holds a shared_ptr
         * to a MalariaConfig, enabling thread-safe operation with independent configs.
         */
        struct MalariaConfig
        {
            // =====================================================
            // IntrahostComponent-level parameters
            // =====================================================

            int randomSeed = 0;
            int max_ind_inf = 1;

            int falciparumMSPVars = DEFAULT_MSP_VARIANTS;
            int falciparumNonSpecTypes = DEFAULT_NONSPECIFIC_TYPES;
            int falciparumPfEMP1Vars = DEFAULT_PFEMP1_VARIANTS;

            float base_gametocyte_mosquito_survival = DEFAULT_BASE_GAMETOCYTE_MOSQUITO_SURVIVAL;
            float cytokine_gametocyte_inactivation = DEFAULT_CYTOKINE_GAMETOCYTE_INACTIVATION;

            // =====================================================
            // Susceptibility parameters
            // =====================================================

            float memory_level = 0.2f;
            float hyperimmune_decay_rate = 0.0f;  // derived from memory_level
            float MSP1_antibody_growthrate = 0.02f;
            float antibody_stimulation_c50 = 10.0f;
            float antibody_capacity_growthrate = 0.1f;
            float minimum_adapted_response = 0.02f;
            float non_specific_growth = 0.5f;
            float antibody_csp_decay_days = DEFAULT_ANTIBODY_CSP_DECAY_DAYS;
            float antibody_days_to_long_term_decay = 730.0f;   // 2 years before long-term decay starts
            float antibody_long_term_decay_days = 7300.0f;     // 20 year decay time constant

            float maternal_antibody_decay_rate = 0.01f;

            float PfHRP2_boost_rate = 7.0e-14f;    // Picograms/iRBC/day
            float PfHRP2_decay_rate = 0.172f;      // Fraction per day (3.67 day half-life)

            float pyrogenic_threshold = 1000.0f;
            float fever_IRBC_killrate = DEFAULT_FEVER_IRBC_KILL_RATE;

            float erythropoiesis_anemia_effect = 3.5f;

            // =====================================================
            // Infection parameters
            // =====================================================

            float incubation_period = 7.0f;
            float antibody_IRBC_killrate = DEFAULT_ANTIBODY_IRBC_KILLRATE;
            float non_specific_antigenicity = DEFAULT_NON_SPECIFIC_ANTIGENICITY;
            float MSP1_merozoite_kill = DEFAULT_MSP1_MEROZOITE_KILL;
            float gametocyte_stage_survival = DEFAULT_GAMETOCYTE_STAGE_SURVIVAL;
            float base_gametocyte_sexratio = DEFAULT_BASE_GAMETOCYTE_SEX_RATIO;
            float base_gametocyte_production = DEFAULT_BASE_GAMETOCYTE_PRODUCTION;
            float antigen_switch_rate = DEFAULT_ANTIGEN_SWITCH_RATE;
            float merozoites_per_hepatocyte = DEFAULT_MEROZOITES_PER_HEPATOCYTE;
            float merozoites_per_schizont = DEFAULT_MEROZOITES_PER_SCHIZONT;
            float RBC_destruction_multiplier = DEFAULT_RBC_DESTRUCTION_MULTIPLIER;
            int n_asexual_cycles_wo_gametocytes = DEFAULT_ASEXUAL_CYCLES_WITHOUT_GAMETOCYTES;

            // =====================================================
            // Instance state (owned by config)
            // =====================================================

            std::shared_ptr<RANDOMBASE> rng;
            suids::distributed_generator suidGenerator;

            // =====================================================
            // Methods
            // =====================================================

            MalariaConfig();

            /**
             * @brief Configure from a ParamSet (Python dict/YAML).
             *
             * This replaces the static Configure() methods in the params structs.
             */
            void Configure(const ParamSet& pset);

            /**
             * @brief Create a configured MalariaConfig from a ParamSet.
             */
            static std::shared_ptr<MalariaConfig> FromParamSet(const ParamSet& pset);
        };

    }

}
