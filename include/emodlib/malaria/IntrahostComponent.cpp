/**
 * @file IntrahostComponent.cpp
 *
 * @brief Malaria intrahost component implementation
 */

#include "IntrahostComponent.h"

#include "emodlib/utils/Common.h"
#include "emodlib/utils/Sigmoid.h"


namespace emodlib
{

    namespace malaria
    {

        IntrahostComponent::IntrahostComponent(std::shared_ptr<MalariaConfig> config)
            : m_config(config)
            , susceptibility(nullptr)
            , infections()
        {

        }

        IntrahostComponent* IntrahostComponent::Create(std::shared_ptr<MalariaConfig> config)
        {
            IntrahostComponent* ic = new IntrahostComponent(config);
            ic->susceptibility = Susceptibility::Create(config.get());
            return ic;
        }

        // TODO: emodlib#7 (infectiousness calculations)

        void IntrahostComponent::Update(float dt)
        {
            // TODO: emodlib#5 (mature gametocyte decay) + emodlib#4 (mature gametocyte drug killing)

            // Sub-step to match EMOD's Infection_Updates_Per_Timestep = 8.
            // Without sub-stepping, large parasite bursts can deplete all RBCs
            // in a single step before erythropoiesis feedback kicks in.
            int n_substeps = 8;
            float sub_dt = dt / n_substeps;

            for (int step = 0; step < n_substeps; step++) {
                // Match EMOD order (Individual.cpp): infections first, then susceptibility
                for (auto it = infections.begin(); it != infections.end();) {
                    (*it)->Update(sub_dt);

                    // TODO: emodlib#3 (InfectionStateChange::Cleared)

                    if ((*it)->IsCleared()) {
                        delete *it;
                        it = infections.erase(it);
                        continue;
                    }

                    ++it;
                }

                susceptibility->Update(sub_dt);
            }
        }

        void IntrahostComponent::Challenge()
        {
            if (infections.size() < m_config->max_ind_inf) {
                Infection* inf = Infection::Create(susceptibility, m_config.get());
                infections.push_back(inf);  // TODO: emodlib#2 (Max_Individual_Infections)
            }
        }

        void IntrahostComponent::ChallengeWithAntigens(const DebugAntigens& antigens)
        {
            if (infections.size() < m_config->max_ind_inf) {
                Infection* inf = Infection::CreateWithAntigens(susceptibility, m_config.get(), antigens);
                infections.push_back(inf);
            }
        }

        bool IntrahostComponent::ChallengeWithSporozoites(int n_sporozoites)
        {
            if (infections.size() >= m_config->max_ind_inf) return false;

            float survival_prob = m_config->base_sporozoite_survival_fraction;

            // CSP antibody reduces sporozoite survival
            IMalariaAntibody* csp = susceptibility->get_csp_antibody();
            float csp_conc = csp->GetAntibodyConcentration();
            if (csp_conc > 0) {
                survival_prob *= (1.0f - Sigmoid::variableWidthSigmoid(
                    log10f(csp_conc),
                    log10f(m_config->antibody_csp_killing_threshold),
                    m_config->antibody_csp_killing_invwidth
                ));
            }

            // Auto-boost CSP antibody from exposure (matches EMOD behavior)
            // Use a slow growth rate (~3 year time constant: 0.001/day)
            csp->UpdateAntibodyCapacityByRate(1.0f, 0.001f);
            csp->SetAntigenicPresence(true);

            int hepatocytes = m_config->rng->Poisson(n_sporozoites * survival_prob);

            if (hepatocytes > 0) {
                Infection* inf = Infection::Create(susceptibility, m_config.get(), hepatocytes);
                infections.push_back(inf);
                return true;
            }
            return false;
        }

        bool IntrahostComponent::ChallengeWithBites(int n_bites)
        {
            int n_sporozoites = static_cast<int>(n_bites * m_config->mean_sporozoites_per_bite);
            return ChallengeWithSporozoites(n_sporozoites);
        }

        void IntrahostComponent::Treat()
        {
            infections.clear();  // TODO: emodlib#4 (asexual drug killing) + emodlib#3 (InfectionStateChange::Cleared)
        }

        int IntrahostComponent::GetNumInfections() const
        {
            return infections.size();
        }

        float IntrahostComponent::GetParasiteDensity() const
        {
            float total = 0.0f;
            for (auto* inf: infections) {
                total += inf->get_asexual_density();
            }
            return total;
        }

        float IntrahostComponent::GetGametocyteDensity() const
        {
            float total = 0.0f;
            for (auto* inf: infections) {
                total += inf->get_mature_gametocyte_density();  // TODO: emodlib#5 (mature gametocyte decay)
            }
            return total;
        }

        float IntrahostComponent::GetFeverTemperature() const
        {
            return susceptibility->get_fever_celsius();
        }

        float IntrahostComponent::GetInfectiousness() const
        {
            float cytokines = susceptibility->get_cytokines();
            double fever_effect = Sigmoid::basic_sigmoid(m_config->cytokine_gametocyte_inactivation, cytokines);
            return EXPCDF(-GetGametocyteDensity() * MICROLITERS_PER_BLOODMEAL * m_config->base_gametocyte_mosquito_survival * (1.0 - fever_effect));
        }

        Susceptibility* IntrahostComponent::GetSusceptibility() const
        {
            return susceptibility;
        }

        std::list<Infection*> IntrahostComponent::GetInfections() const
        {
            return infections;
        }

        std::shared_ptr<MalariaConfig> IntrahostComponent::GetConfig() const
        {
            return m_config;
        }

        // Aggregate infection state getters (for validation)

        std::vector<int64_t> IntrahostComponent::get_total_irbc_counts() const
        {
            std::vector<int64_t> totals(CLONAL_PfEMP1_VARIANTS, 0);
            for (auto* inf : infections)
            {
                auto counts = inf->get_irbc_counts();
                for (size_t i = 0; i < counts.size(); i++)
                {
                    totals[i] += counts[i];
                }
            }
            return totals;
        }

        int64_t IntrahostComponent::get_total_hepatocytes() const
        {
            int64_t total = 0;
            for (auto* inf : infections)
            {
                total += inf->get_hepatocyte_count();
            }
            return total;
        }

        std::vector<int64_t> IntrahostComponent::get_total_male_gametocytes() const
        {
            std::vector<int64_t> totals(GametocyteStages::Count, 0);
            for (auto* inf : infections)
            {
                auto counts = inf->get_all_male_gametocytes();
                for (size_t i = 0; i < counts.size(); i++)
                {
                    totals[i] += counts[i];
                }
            }
            return totals;
        }

        std::vector<int64_t> IntrahostComponent::get_total_female_gametocytes() const
        {
            std::vector<int64_t> totals(GametocyteStages::Count, 0);
            for (auto* inf : infections)
            {
                auto counts = inf->get_all_female_gametocytes();
                for (size_t i = 0; i < counts.size(); i++)
                {
                    totals[i] += counts[i];
                }
            }
            return totals;
        }

        int32_t IntrahostComponent::get_total_asexual_cycles() const
        {
            int32_t total = 0;
            for (auto* inf : infections)
            {
                total += inf->get_asexual_cycle_count();
            }
            return total;
        }

    }

}
