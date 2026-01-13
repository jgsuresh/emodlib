/**
 * @file VectorCohort.h
 *
 * @brief Cohort data structure for tracking groups of mosquitoes.
 *
 * A cohort represents a group of mosquitoes that share the same state
 * (age, infection status, development progress). This is more efficient
 * than tracking individual mosquitoes while preserving key dynamics.
 */

#pragma once

#include <cstdint>

namespace emodlib
{
    namespace vector
    {
        /**
         * @brief Lifecycle state of a vector cohort.
         */
        enum class VectorState : uint8_t
        {
            EGG,        ///< Egg stage
            LARVA,      ///< Larval stage (aquatic)
            IMMATURE,   ///< Immature/pupal stage
            ADULT,      ///< Adult, susceptible (can be infected)
            INFECTED,   ///< Infected, in EIP (not yet infectious)
            INFECTIOUS, ///< Infectious (can transmit to humans)
            MALE        ///< Adult male (for completeness, simplified)
        };

        /**
         * @brief A cohort of mosquitoes sharing the same state.
         *
         * Cohorts track population count, age, development progress,
         * and infection status. The population is tracked as a float
         * to handle fractional mortality/splitting without rounding errors.
         */
        struct VectorCohort
        {
            // =====================================================
            // Core state
            // =====================================================

            float population;           ///< Number of mosquitoes in cohort
            float age;                  ///< Age in days
            float progress;             ///< Development progress [0, 1+]
            VectorState state;          ///< Current lifecycle state
            int days_since_infection;   ///< Days since infected (-1 if susceptible)

            // =====================================================
            // Constructors
            // =====================================================

            /**
             * @brief Default constructor - creates empty cohort.
             */
            VectorCohort()
                : population(0.0f)
                , age(0.0f)
                , progress(0.0f)
                , state(VectorState::EGG)
                , days_since_infection(-1)
            {
            }

            /**
             * @brief Construct a cohort with specified parameters.
             */
            VectorCohort(float pop, VectorState s, float age_days = 0.0f)
                : population(pop)
                , age(age_days)
                , progress(0.0f)
                , state(s)
                , days_since_infection(s == VectorState::INFECTED ? 0 : -1)
            {
            }

            // =====================================================
            // State queries
            // =====================================================

            /**
             * @brief Check if cohort is susceptible (can be infected).
             */
            bool IsSusceptible() const
            {
                return state == VectorState::ADULT && days_since_infection < 0;
            }

            /**
             * @brief Check if cohort is exposed (infected but not yet infectious).
             */
            bool IsExposed() const
            {
                return state == VectorState::INFECTED;
            }

            /**
             * @brief Check if cohort is infectious (can transmit).
             */
            bool IsInfectious() const
            {
                return state == VectorState::INFECTIOUS;
            }

            /**
             * @brief Check if cohort is in aquatic stage (egg, larva, immature).
             */
            bool IsAquatic() const
            {
                return state == VectorState::EGG ||
                       state == VectorState::LARVA ||
                       state == VectorState::IMMATURE;
            }

            /**
             * @brief Check if cohort is an adult (any adult state).
             */
            bool IsAdult() const
            {
                return state == VectorState::ADULT ||
                       state == VectorState::INFECTED ||
                       state == VectorState::INFECTIOUS ||
                       state == VectorState::MALE;
            }

            /**
             * @brief Check if cohort is empty (negligible population).
             */
            bool IsEmpty() const
            {
                return population < 0.01f;
            }

            // =====================================================
            // State transitions
            // =====================================================

            /**
             * @brief Apply mortality - reduce population by survival probability.
             *
             * @param survival_prob Probability of surviving (0 to 1)
             */
            void ApplyMortality(float survival_prob)
            {
                population *= survival_prob;
            }

            /**
             * @brief Increment age by timestep.
             *
             * @param dt Timestep in days
             */
            void IncrementAge(float dt)
            {
                age += dt;
            }

            /**
             * @brief Add development progress.
             *
             * @param delta_progress Progress to add
             */
            void AddProgress(float delta_progress)
            {
                progress += delta_progress;
            }

            /**
             * @brief Check if development is complete (ready to transition).
             */
            bool IsProgressComplete() const
            {
                return progress >= 1.0f;
            }

            /**
             * @brief Reset progress after stage transition.
             */
            void ResetProgress()
            {
                progress = 0.0f;
            }

            /**
             * @brief Transition to next lifecycle state.
             *
             * @param new_state New state after transition
             */
            void TransitionTo(VectorState new_state)
            {
                state = new_state;
                progress = 0.0f;

                // Reset age for certain transitions
                if (new_state == VectorState::ADULT ||
                    new_state == VectorState::MALE)
                {
                    age = 0.0f;
                }

                // Set infection tracking for infected state
                if (new_state == VectorState::INFECTED)
                {
                    days_since_infection = 0;
                }
            }

            /**
             * @brief Increment infection days (for EIP tracking).
             *
             * @param dt Timestep in days (typically 1)
             */
            void IncrementInfectionDays(int dt = 1)
            {
                if (days_since_infection >= 0)
                {
                    days_since_infection += dt;
                }
            }

            // =====================================================
            // Splitting and merging
            // =====================================================

            /**
             * @brief Split off a fraction of this cohort.
             *
             * @param fraction Fraction to split off (0 to 1)
             * @return New cohort with the split population
             */
            VectorCohort Split(float fraction)
            {
                VectorCohort split_cohort = *this;
                float split_pop = population * fraction;
                split_cohort.population = split_pop;
                population -= split_pop;
                return split_cohort;
            }

            /**
             * @brief Split off a specific number from this cohort.
             *
             * @param count Number to split off
             * @return New cohort with the specified count
             */
            VectorCohort SplitCount(float count)
            {
                VectorCohort split_cohort = *this;
                float actual_count = (count < population) ? count : population;
                split_cohort.population = actual_count;
                population -= actual_count;
                return split_cohort;
            }

            /**
             * @brief Merge another cohort into this one.
             *
             * Only merges if states match. Ages are population-weighted averaged.
             *
             * @param other Cohort to merge
             * @return true if merged, false if incompatible
             */
            bool Merge(const VectorCohort& other)
            {
                if (state != other.state ||
                    days_since_infection != other.days_since_infection)
                {
                    return false;
                }

                // Population-weighted average of age and progress
                float total_pop = population + other.population;
                if (total_pop > 0.0f)
                {
                    age = (age * population + other.age * other.population) / total_pop;
                    progress = (progress * population + other.progress * other.population) / total_pop;
                }
                population = total_pop;
                return true;
            }
        };

    }  // namespace vector

}  // namespace emodlib
