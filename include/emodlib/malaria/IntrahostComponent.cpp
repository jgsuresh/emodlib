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

        float IntrahostComponent::increment_parasite = 1.0f;
        float IntrahostComponent::increment_gametocyte = 1.0f;
        float IntrahostComponent::increment_fever = 1.0f;
    
        void IntrahostComponent::Configure(const ParamSet& pset)
        {
            increment_parasite = pset["increment_parasite"].cast<float>();
            increment_gametocyte = pset["increment_gametocyte"].cast<float>();
            increment_fever = pset["increment_fever"].cast<float>();
        }
    
        IntrahostComponent::IntrahostComponent()
            : susceptibility(nullptr)
            , infections()
        {

        }

        IntrahostComponent* IntrahostComponent::Create()
        {
            IntrahostComponent* ic = new IntrahostComponent();
            ic->susceptibility = new SusceptibilityMalaria();
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
            infections.push_back(new InfectionMalaria());
        }

        void IntrahostComponent::Treat()
        {
            infections.clear();
        }

        float IntrahostComponent::GetParasiteDensity()
        {
            float total = 0.0f;
            for (auto* inf: infections) {
                total += inf->GetParasiteDensity();
            }
            return total;
        }

        float IntrahostComponent::GetGametocyteDensity()
        {
            float total = 0.0f;
            for (auto* inf: infections) {
                total += inf->GetGametocyteDensity();
            }
            return total;
        }

        float IntrahostComponent::GetFeverTemperature()
        {
            return susceptibility->GetFeverTemperature();
        }

    }

}
