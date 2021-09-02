set(EMODLIB_MALARIA_TARGET emodlib-malaria)

find_package(pybind11 CONFIG REQUIRED)
message(STATUS "Found pybind11 v${pybind11_VERSION}: ${pybind11_INCLUDE_DIRS}")

add_library(${EMODLIB_MALARIA_TARGET} SHARED
    IntrahostComponent.cpp
    InfectionMalaria.cpp
    SusceptibilityMalaria.cpp
    MalariaAntibody.cpp)

# include top-level directory to allow full include paths
target_include_directories(${EMODLIB_MALARIA_TARGET} PUBLIC ${CMAKE_SOURCE_DIR}/include)

target_link_libraries(${EMODLIB_MALARIA_TARGET} PRIVATE pybind11::embed emodlib-utils)

target_compile_features(${EMODLIB_MALARIA_TARGET} PUBLIC cxx_std_14)
