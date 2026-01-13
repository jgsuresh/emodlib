/**
 * @file VectorConfig.h
 *
 * @brief Vector population configuration - species parameters and model settings.
 *
 * Default values are for Anopheles gambiae based on EMOD parameters.
 * Arrhenius equations control temperature-dependent development rates.
 */

#pragma once

#include <cmath>

namespace emodlib
{
    namespace vector
    {
        /**
         * @brief Configuration for vector population dynamics.
         *
         * Contains species-specific biological parameters and model settings.
         * Default values represent Anopheles gambiae in tropical conditions.
         */
        struct VectorConfig
        {
            // =====================================================
            // Arrhenius parameters for temperature-dependent rates
            // =====================================================

            // Larval/aquatic development: progress = A1 * exp(-A2 / (T + 273.15)) * dt
            float aquatic_arrhenius1 = 1.5e6f;   // A coefficient
            float aquatic_arrhenius2 = 5000.0f;  // B coefficient (activation energy)

            // Infected (EIP) development: progress = A1 * exp(-A2 / (T + 273.15)) * dt
            float infected_arrhenius1 = 6500.0f;
            float infected_arrhenius2 = 3500.0f;

            // Feeding cycle duration (optional temperature dependence)
            float cycle_arrhenius1 = 4.09e10f;
            float cycle_arrhenius2 = 7740.0f;

            // =====================================================
            // Life expectancy and mortality
            // =====================================================

            float adult_life_expectancy = 21.0f;   // Days (An. gambiae field estimate)
            float male_life_expectancy = 10.0f;    // Days
            float aquatic_mortality_rate = 0.2f;   // Per day (20%)
            float immature_duration = 2.0f;        // Days in immature stage

            // Age-dependent mortality (Styer et al. 2007)
            // mortality_from_age = 0.006 * e / (1 + 0.045 * (e - 1)) where e = exp(0.2 * age)
            bool enable_age_dependent_mortality = true;

            // =====================================================
            // Feeding and host-seeking behavior
            // =====================================================

            float anthropophily = 0.9f;            // Fraction of feeds on humans (vs animals)
            float indoor_feeding_fraction = 1.0f;  // Fraction feeding indoors
            float days_between_feeds = 3.0f;       // Gonotrophic cycle length

            // =====================================================
            // Reproduction
            // =====================================================

            float egg_batch_size = 100.0f;         // Eggs per oviposition
            float infected_egg_batch_factor = 0.8f; // Reduction when infected
            float egg_survival_rate = 0.99f;       // Daily egg survival

            // =====================================================
            // Transmission parameters
            // =====================================================

            float transmission_rate = 0.02f;       // V->H per infectious bite
            float acquire_modifier = 0.5f;         // H->V probability modifier
            float sporozoites_per_bite = 11.0f;    // Average sporozoites delivered

            // =====================================================
            // Population parameters
            // =====================================================

            float carrying_capacity = 10000.0f;    // Equilibrium adult population
            float initial_infected_fraction = 0.0f; // For endemic initialization

            // =====================================================
            // Computed values (derived from above)
            // =====================================================

            /**
             * @brief Get daily adult mortality rate.
             */
            float GetAdultMortalityRate() const
            {
                return 1.0f / adult_life_expectancy;
            }

            /**
             * @brief Get daily male mortality rate.
             */
            float GetMaleMortalityRate() const
            {
                return 1.0f / male_life_expectancy;
            }

            /**
             * @brief Get biting rate (bites per mosquito per day).
             */
            float GetBitingRate() const
            {
                return 1.0f / days_between_feeds;
            }

            /**
             * @brief Get daily emergence rate to maintain carrying capacity.
             *
             * At equilibrium: emergence = carrying_capacity * mortality_rate
             */
            float GetDailyEmergence() const
            {
                return carrying_capacity * GetAdultMortalityRate();
            }

            // =====================================================
            // Arrhenius helper functions
            // =====================================================

            /**
             * @brief Calculate larval development progress for one timestep.
             *
             * @param temp_celsius Temperature in Celsius
             * @param dt Timestep in days
             * @return Progress [0,1] toward next stage
             */
            float GetLarvalProgress(float temp_celsius, float dt) const
            {
                float temp_kelvin = temp_celsius + 273.15f;
                return aquatic_arrhenius1 * std::exp(-aquatic_arrhenius2 / temp_kelvin) * dt;
            }

            /**
             * @brief Calculate infected (EIP) progress for one timestep.
             *
             * @param temp_celsius Temperature in Celsius
             * @param dt Timestep in days
             * @return Progress [0,1] toward becoming infectious
             */
            float GetInfectedProgress(float temp_celsius, float dt) const
            {
                float temp_kelvin = temp_celsius + 273.15f;
                return infected_arrhenius1 * std::exp(-infected_arrhenius2 / temp_kelvin) * dt;
            }

            /**
             * @brief Calculate age-dependent mortality rate (Styer et al. 2007).
             *
             * Based on Aedes aegypti lab data, commonly applied to Anopheles.
             * Represents senescence-related mortality increase.
             *
             * @param age_days Age in days
             * @return Additional mortality rate (per day)
             */
            static float GetAgeDependentMortality(float age_days)
            {
                float e = std::exp(0.2f * age_days);
                return 0.006f * e / (1.0f + 0.045f * (e - 1.0f));
            }

            /**
             * @brief Get combined survival probability for one timestep.
             *
             * Combines baseline mortality with age-dependent component.
             *
             * @param age_days Age in days
             * @param base_mortality_rate Baseline daily mortality rate
             * @return Probability of surviving one day
             */
            float GetSurvivalProbability(float age_days, float base_mortality_rate) const
            {
                float total_mortality = base_mortality_rate;
                if (enable_age_dependent_mortality)
                {
                    total_mortality += GetAgeDependentMortality(age_days);
                }
                return std::exp(-total_mortality);
            }
        };

    }  // namespace vector

}  // namespace emodlib
