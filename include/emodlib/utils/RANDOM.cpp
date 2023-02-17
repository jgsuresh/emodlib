#include "RANDOM.h"

#include <assert.h>
#include <math.h>

#ifndef WIN32
#include <memory.h>    // memset
#include <climits>     // UINT_MAX
#if defined(__APPLE__) && defined(__arm64__)  // intended for Apple M1 hardware
#include "sse2neon.h"
#else
#include <emmintrin.h> // __m128i
#include <wmmintrin.h> // AES
#include <smmintrin.h> // _mm_insert_epi64
#include <tmmintrin.h> // _mm_shuffle_epi8
#endif
#endif

#define PRNG_COUNT  (1<<20) // Let's start with ~1 million


namespace emodlib
{

    // ----------------------------------------------------------------------------
    // --- RANDOMBASE
    // ----------------------------------------------------------------------------

    RANDOMBASE::RANDOMBASE( size_t nCache )
        : cache_count( nCache )
        , index( UINT_MAX )   // Make sure fill_bits() is called...
        , random_bits( nullptr )
        , random_floats( nullptr )
        , bGauss( false )
        , eGauss_( 0.0f )
        {
            if( cache_count == 0 )
            {
                cache_count = PRNG_COUNT;
            }
            random_bits = reinterpret_cast<uint32_t*>(malloc( cache_count * sizeof( uint32_t ) ));
            random_floats = reinterpret_cast<float*>(malloc( cache_count * sizeof( float ) ));
        }

    RANDOMBASE::~RANDOMBASE()
    {
    #ifdef _DEBUG
        free(random_bits);
        free(random_floats);
    #endif
    }

    uint32_t RANDOMBASE::ul()
    {
        if (index >= cache_count)
        {
            fill_bits();
            bits_to_float();
            index = 0;
        }

        return random_bits[index++];
    }

    float RANDOMBASE::e()
    {
        if (index >= cache_count)
        {
            fill_bits();
            bits_to_float();
            index = 0;
        }

        return random_floats[index++];
    }

    // Finds an uniformally distributed number between 0 (inclusive) and N (exclusive)
    uint16_t RANDOMBASE::uniformZeroToN16( uint16_t N )
    {
        uint32_t ulA = ul();
        uint32_t ll = (ulA & 0xFFFFL) * N;
        ll >>= 16;
        ll += (ulA >> 16) * N;
        return (uint16_t)(ll >> 16);
    }

    // Finds an uniformally distributed number between 0 (inclusive) and N (exclusive)
    uint32_t RANDOMBASE::uniformZeroToN32( uint32_t N )
    {
        uint64_t ulA = uint64_t( ul() );
        uint64_t ulB = uint64_t( ul() );
        ulB <<= 32;
        ulA += ulB;
        uint64_t ll = (ulA & 0xFFFFFFFFL) * N;
        ll >>= 32;
        return uint32_t( ll );
    }

    void RANDOMBASE::fill_bits()
    {
        assert(false);
    }

    #define FLOAT_EXP   8
    #define DOUBLE_EXP 11

    void RANDOMBASE::bits_to_float()
    {
        __m128i m = _mm_set1_epi32(0x007FFFFF);
        __m128i o = _mm_set1_epi32(0x00000001);
        __m128 f = _mm_set1_ps(1.0f);
        __m128i fi = _mm_castps_si128(f);
        for (size_t i = 0; i < cache_count; i += 4)
        {
            __m128i x = _mm_load_si128(reinterpret_cast<__m128i const*>(random_bits+i));    // x = bits
    //        x = _mm_and_si128(x, m);                                    // x &= 0x007FFFFF
            x = _mm_srli_epi32(x, (FLOAT_EXP+1));                       // x = x >> 9 (we just want the 23 mantissa bits)
            x = _mm_or_si128(x, o);                                     // x |= 0x00000001
            __m128i y = _mm_or_si128(fi, x);                            // y = fi | x
            __m128 z = _mm_castsi128_ps(y);                             // z = y interpreted as floating point
            z = _mm_sub_ps(z, f);                                       // z -= 1.0f
            _mm_store_ps(random_floats + i, z);
        }
    }

    /******************************Public*Routine******************************\
    * eGauss()                                                                 *
    *                                                                          *
    * Returns a normal deviate for each call.  It generates two deviates at a  *
    * time and caches one for the next call.                                   *
    *                                                                          *
    *  Sat 26-Mar-1994 14:10:12 -by- Charles Whitmer [chuckwh]                 *
    * I actually wrote this back in 1982 for some physics stuff!               *
    \**************************************************************************/

    double RANDOMBASE::eGauss()
    {
        if (bGauss)
        {
            bGauss = false;
            return eGauss_;
        }

        double rad, norm;
        double s, r1, r2;

        rad = -2.0 * log(ee());
        do
        {
            r1 = ee() - 0.5;
            r2 = ee() - 0.5;
            s = r1 * r1 + r2 * r2;
        }
        while (s > 0.25);
        norm = sqrt(rad / s);
        eGauss_ = r1 * norm;
        bGauss = true;
        return r2 * norm;
    }

    double RANDOMBASE::ee()
    {
        union
        {
            double ee;
            struct
            {
                uint32_t Low;
                uint32_t High;
            };
        } ee_ul;

        uint32_t ll = ul();    // Choose a random 32 bits.

        ee_ul.ee = 1.0;
        ee_ul.High += (ll >> (DOUBLE_EXP + 1));
        ee_ul.Low = (ll << (31 - DOUBLE_EXP)) + (1 << (30 - DOUBLE_EXP));

        return ee_ul.ee - 1.0;
    }

    // Poisson() added by Philip Eckhoff, uses Gaussian approximation for ratetime>10
    uint64_t RANDOMBASE::Poisson(double ratetime)
    {
        if (ratetime <= 0)
        {
            return 0;
        }
        uint64_t events = 0;
        double Time = 0;
        double tempval;
        if (ratetime < 10)
        {
            while (Time < 1)
            {
                Time += -log(e()) / ratetime;
                if (Time < 1)
                {
                    events++;
                }
            }
        }
        else
        {
            tempval = (eGauss() * sqrt(ratetime) + ratetime + .5);
            if (tempval < 0)
            {
                events = 0;
            }
            else
            {
                events = uint64_t(tempval);
            }
        }
        return events;
    }

    // Poisson_true added by Philip Eckhoff, actual Poisson, without approximation
    uint32_t RANDOMBASE::Poisson_true(double ratetime)
    {
        if (ratetime <= 0)
        {
            return 0;
        }
        uint32_t events = 0;
        double Time = 0;
        while (Time < 1)
        {
            Time += -log(e()) / ratetime;
            if (Time < 1)
            {
                events++;
            }
        }
        return events;
    }

    // ----------------------------------------------------------------------------
    // --- PSEUDO_DES
    // ----------------------------------------------------------------------------

    PSEUDO_DES::PSEUDO_DES( uint64_t iSequence, size_t nCache )
        : RANDOMBASE( nCache )
        , iSeq( uint32_t( iSequence & 0xFFFFFFFF ) ) // lower 32-bits
        , iNum( uint32_t( iSequence >> 32        ) ) // upper 32-bits
    {
    }

    PSEUDO_DES::~PSEUDO_DES()
    {
    }

    const uint32_t c1[4] = {0xBAA96887L, 0x1E17D32CL, 0x03BCDC3CL, 0x0F33D1B2L};
    const uint32_t c2[4] = {0x4B0F3B58L, 0xE874F0C3L, 0x6955C5A6L, 0x55A7CA46L};

    #define HI(x) ((uint32_t) ((uint16_t*) &x)[1])
    #define LO(x) ((uint32_t) ((uint16_t*) &x)[0])
    #define XCHG(x) ((LO(x) << 16) | HI(x))

    void PSEUDO_DES::fill_bits()
    {
        uint32_t kk[3];
        uint32_t iA;
        uint32_t iB;
    #ifdef _DEBUG
        uint32_t ul;
    #endif

        for (size_t i = 0; i < cache_count; ++i)
        {
            iA = iNum ^ c1[0];
            iB = LO(iA) * LO(iA) + ~(HI(iA) * HI(iA));
            kk[0] = iSeq ^ ((XCHG(iB) ^ c2[0]) + LO(iA) * HI(iA));

            iA = kk[0] ^ c1[1];
            iB = LO(iA) * LO(iA) + ~(HI(iA) * HI(iA));
            kk[1] = iNum ^ ((XCHG(iB) ^ c2[1]) + LO(iA) * HI(iA));

            iNum++;
            if (iNum == 0)
                iSeq++;

            iA = kk[1] ^ c1[2];
            iB = LO(iA) * LO(iA) + ~(HI(iA) * HI(iA));
            kk[2] = kk[0] ^ ((XCHG(iB) ^ c2[2]) + LO(iA) * HI(iA));

            iA = kk[2] ^ c1[3];
            iB = LO(iA) * LO(iA) + ~(HI(iA) * HI(iA));

            random_bits[i] =
    #ifdef _DEBUG
                ul =
    #endif
                kk[1] ^ ((XCHG(iB) ^ c2[3]) + LO(iA) * HI(iA));
        }
    }

}
