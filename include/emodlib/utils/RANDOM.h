/******************************Module*Header*******************************\
* Module Name: random.hxx                                                  *
*                                                                          *
* A random number class.  It depends on a 32 bit ULONG type and 16 bit     *
* USHORT type.                                                             *
*                                                                          *
* Created: 25-Mar-1994 12:34:57                                            *
* Author: Charles Whitmer [chuckwh]                                        *
*                                                                          *
* Copyright (c) 1994 Charles Whitmer                                       *
\**************************************************************************/

#pragma once

#include <stdint.h>
#include <vector>
#include <set>


namespace emodlib
{

    // ------------------------------------------------------------------------
    // --- RANDOMBASE
    // ------------------------------------------------------------------------
    class RANDOMBASE
    {
        
    public:

        RANDOMBASE( size_t nCache );
        virtual ~RANDOMBASE();

        uint32_t ul();  // Returns a random 32 bit number.
        float e();      // Returns a randon float between 0 and 1.

        // Finds an uniformally distributed number less than N
        // The 16-bit version requires one less random number than
        // the 32-bit version, i.e. less overhead.
        uint16_t uniformZeroToN16( uint16_t N );
        uint32_t uniformZeroToN32( uint32_t N );

        double ee();
        double eGauss();    // Returns a normal deviate.

        // Added by Philip Eckhoff, Poisson takes in a rate, and returns the number of events in unit time
        // Or equivalently, takes in rate*time and returns number of events in that time
        // Poisson uses a Gaussian approximation for large lambda, while Poisson_true is the fully accurate Poisson
        // expdist takes in a rate and returns the sample from an exponential distribution
        uint64_t Poisson(double=1.0);
        uint32_t Poisson_true(double=1.0);

    protected:

        virtual void fill_bits();
        void bits_to_float();

        size_t    cache_count;
        size_t    index;
        uint32_t* random_bits;
        float*    random_floats;
            
        bool   bGauss;
        double eGauss_;

    };


    // ------------------------------------------------------------------------
    // --- PSEUDO_DES
    // ------------------------------------------------------------------------
    // Numerical Recipes in C, 2nd ed. Press, William H. et. al, 1992.

    class PSEUDO_DES : public RANDOMBASE
    {

    public:
        PSEUDO_DES( uint64_t iSequence = 0, size_t nCache = 0 );
        ~PSEUDO_DES();

    protected:
        virtual void fill_bits() override;

        uint32_t iSeq;
        uint32_t iNum;
    };

}
