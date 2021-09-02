#pragma once

#include <math.h>

namespace emodlib
{

    struct Sigmoid
    {
        inline static double basic_sigmoid ( double threshold = 100.0, double variable = 0.0 )
        {
            return (variable > 0) ? (variable / (threshold + variable)) : 0.0;
        }


        inline static float sigmoid(float x)
        {
            return 0.5f + 0.5f * tanh(x / 2.0f); // instead of 1/(1+exp(x)), prevents exp() from overflowing
        }
    };

}
