set(EMODLIB_MALARIA_TARGET emodlib-malaria)

add_library(${EMODLIB_MALARIA_TARGET} SHARED
    IntrahostComponent.cpp
    InfectionMalaria.cpp
    SusceptibilityMalaria.cpp)

# include top-level directory to allow full include paths
target_include_directories(${EMODLIB_MALARIA_TARGET} PUBLIC ${CMAKE_SOURCE_DIR})

target_compile_features(${EMODLIB_MALARIA_TARGET} PUBLIC cxx_std_14)
