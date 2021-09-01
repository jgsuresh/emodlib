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

        void IntrahostComponent::params::Configure(const ParamSet& pset)
        {
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
            ic->susceptibility = new Susceptibility();
            return ic;
        }

        void IntrahostComponent::Update()
        {
            susceptibility->Update();

            for (auto* inf : infections) {
                inf->Update();
            }
        }

        void IntrahostComponent::Challenge()
        {
            infections.push_back(new Infection());
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
            return susceptibility->GetFeverTemperature();
        }

    }

}
