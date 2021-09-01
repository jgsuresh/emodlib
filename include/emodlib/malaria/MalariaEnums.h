/**
 * @file MalariaEnums.h
 *
 * @brief Malaria enumerations
 */

#pragma once

#include "emodlib/utils/EnumSupport.h"

namespace emodlib
{

    namespace malaria
    {

        // ENUM defs for MALARIA_STRAINS
        // FALCIPARUM_RANDOM_STRAIN is default
        ENUM_DEFINE(MalariaStrains,
            ENUM_VALUE_SPEC(FALCIPARUM_NONRANDOM_STRAIN                         , 11)
            ENUM_VALUE_SPEC(FALCIPARUM_RANDOM50_STRAIN                          , 12)
            ENUM_VALUE_SPEC(FALCIPARUM_RANDOM_STRAIN                            , 13)
            ENUM_VALUE_SPEC(FALCIPARUM_STRAIN_GENERATOR                         , 20))
        
        // ENUM defs for PARASITE_SWITCH_TYPE
        // RATE_PER_PARASITE_7VARS is default
        ENUM_DEFINE(ParasiteSwitchType,
            ENUM_VALUE_SPEC(CONSTANT_SWITCH_RATE_2VARS                          , 0)
            ENUM_VALUE_SPEC(RATE_PER_PARASITE_7VARS                             , 1)
            ENUM_VALUE_SPEC(RATE_PER_PARASITE_5VARS_DECAYING                    , 2))
    
        // ENUM defs for MATERNAL_ANTIBODIES_TYPE
        // SIMPLE_WANING draws a PfEMP1 antibody strength is initialized
        //               at a fraction of the mother's level and wanes exponentially
        // CONSTANT_INITIAL_IMMUNITY ignores the mother's level (or mean adult population level) of immunity
        //               and gives some constant initial immunity level to all children, which wanes exponentially.
        ENUM_DEFINE(MaternalAntibodiesType,
            ENUM_VALUE_SPEC(OFF,                       0)
            ENUM_VALUE_SPEC(SIMPLE_WANING,             1)
            ENUM_VALUE_SPEC(CONSTANT_INITIAL_IMMUNITY, 2))
    
        ENUM_DEFINE(InnateImmuneVariationType,
            ENUM_VALUE_SPEC(NONE, 0)
            ENUM_VALUE_SPEC(PYROGENIC_THRESHOLD, 1)
            ENUM_VALUE_SPEC(CYTOKINE_KILLING, 2)
            ENUM_VALUE_SPEC(PYROGENIC_THRESHOLD_VS_AGE, 3))
    
    }

}
