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
            int total = 0.0f;
            for (auto* inf: infections) {
                total += inf->GetParasiteDensity();
            }
            return total;
        }

        float IntrahostComponent::GetGametocyteDensity()
        {
            int total = 0.0f;
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
