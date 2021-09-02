/**
 * @file IntrahostComponent.cpp
 *
 * @brief Malaria intrahost component implementation
 */

#include "IntrahostComponent.h"


namespace emodlib
{

    namespace malaria
    {

        int IntrahostComponent::params::randomSeed = 0;
    
        int IntrahostComponent::params::falciparumMSPVars = DEFAULT_MSP_VARIANTS;
        int IntrahostComponent::params::falciparumNonSpecTypes = DEFAULT_NONSPECIFIC_TYPES;
        int IntrahostComponent::params::falciparumPfEMP1Vars = DEFAULT_PFEMP1_VARIANTS;
    
    
        std::shared_ptr<RANDOMBASE> IntrahostComponent::p_rng = nullptr;

    
        void IntrahostComponent::params::Configure(const ParamSet& pset)
        {
            randomSeed = pset["Run_Number"].cast<int>();
            IntrahostComponent::p_rng = std::shared_ptr<RANDOMBASE>(new PSEUDO_DES(randomSeed, 256));
            
            falciparumMSPVars = pset["Falciparum_MSP_Variants"].cast<int>();
            falciparumNonSpecTypes = pset["Falciparum_Nonspecific_Types"].cast<int>();
            falciparumPfEMP1Vars = pset["Falciparum_PfEMP1_Variants"].cast<int>();

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

        void IntrahostComponent::Update(float dt)
        {
            susceptibility->Update(dt);

            for (auto* inf : infections) {
                inf->Update(dt);
            }
        }

        void IntrahostComponent::Challenge()
        {
            Infection* inf = Infection::Create(susceptibility);
            infections.push_back(inf);
        }

        void IntrahostComponent::Treat()
        {
            infections.clear();
        }

        float IntrahostComponent::GetParasiteDensity() const
        {
            float total = 0.0f;
            for (auto* inf: infections) {
                total += inf->GetParasiteDensity();
            }
            return total;
        }

        float IntrahostComponent::GetGametocyteDensity() const
        {
            float total = 0.0f;
            for (auto* inf: infections) {
                total += inf->GetGametocyteDensity();
            }
            return total;
        }

        float IntrahostComponent::GetFeverTemperature() const
        {
            return susceptibility->get_fever_celsius();
        }

    }

}
