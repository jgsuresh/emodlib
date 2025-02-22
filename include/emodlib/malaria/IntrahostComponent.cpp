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

        int IntrahostComponent::params::randomSeed = 0;

        int IntrahostComponent::params::max_ind_inf = 1;


        int IntrahostComponent::params::falciparumMSPVars = DEFAULT_MSP_VARIANTS;
        int IntrahostComponent::params::falciparumNonSpecTypes = DEFAULT_NONSPECIFIC_TYPES;
        int IntrahostComponent::params::falciparumPfEMP1Vars = DEFAULT_PFEMP1_VARIANTS;

        // TODO: emodlib#7 (infectiousness calculations)
        float  IntrahostComponent::params::base_gametocyte_mosquito_survival = DEFAULT_BASE_GAMETOCYTE_MOSQUITO_SURVIVAL;
        float  IntrahostComponent::params::cytokine_gametocyte_inactivation  = DEFAULT_CYTOKINE_GAMETOCYTE_INACTIVATION;

        std::shared_ptr<RANDOMBASE> IntrahostComponent::p_rng = nullptr;


        void IntrahostComponent::params::Configure(const ParamSet& pset)
        {
            randomSeed = pset["Run_Number"].cast<int>();
            IntrahostComponent::p_rng = std::shared_ptr<RANDOMBASE>(new PSEUDO_DES(randomSeed, 256));

            max_ind_inf = pset["Max_Individual_Infections"].cast<int>();

            falciparumMSPVars = pset["Falciparum_MSP_Variants"].cast<int>();
            falciparumNonSpecTypes = pset["Falciparum_Nonspecific_Types"].cast<int>();
            falciparumPfEMP1Vars = pset["Falciparum_PfEMP1_Variants"].cast<int>();

            // TODO: emodlib#7 (infectiousness calculations)
            base_gametocyte_mosquito_survival = pset["Base_Gametocyte_Mosquito_Survival_Rate"].cast<float>();
            cytokine_gametocyte_inactivation = pset["Cytokine_Gametocyte_Inactivation"].cast<float>();

            Infection::params::Configure(pset["infection_params"]);
            Susceptibility::params::Configure(pset["susceptibility_params"]);
        }


        IntrahostComponent::IntrahostComponent()
            : susceptibility(nullptr)
            , infections()
        {

        }

        IntrahostComponent* IntrahostComponent::Create()
        {
            IntrahostComponent* ic = new IntrahostComponent();
            ic->susceptibility = Susceptibility::Create();
            return ic;
        }

        // TODO: emodlib#7 (infectiousness calculations)

        void IntrahostComponent::Update(float dt)
        {
            // TODO: emodlib#5 (mature gametocyte decay) + emodlib#4 (mature gametocyte drug killing)

            susceptibility->Update(dt);

            for (auto it = infections.begin(); it != infections.end();) {

                (*it)->Update(dt);

                // TODO: emodlib#3 (InfectionStateChange::Cleared)

                if ((*it)->IsCleared()) {
                    delete *it;
                    it = infections.erase(it);
                    continue;
                }

                ++it;
            }
        }

        void IntrahostComponent::Challenge()
        {
            if (infections.size() < params::max_ind_inf) {
                Infection* inf = Infection::Create(susceptibility);
                infections.push_back(inf);  // TODO: emodlib#2 (Max_Individual_Infections)
            }
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
            double fever_effect = Sigmoid::basic_sigmoid(params::cytokine_gametocyte_inactivation, cytokines);
            return EXPCDF(-GetGametocyteDensity() * MICROLITERS_PER_BLOODMEAL * params::base_gametocyte_mosquito_survival * (1.0 - fever_effect));
        }

        Susceptibility* IntrahostComponent::GetSusceptibility() const
        {
            return susceptibility;
        }

        std::list<Infection*> IntrahostComponent::GetInfections() const
        {
            return infections;
        }

    }

}
