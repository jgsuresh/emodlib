cmake_minimum_required(VERSION 3.1)

set(CMAKE_MODULE_PATH
    ${CMAKE_SOURCE_DIR}/cmake
    ${CMAKE_MODULE_PATH})
    
# Dependencies
# ============

# -- Build Python bindings independently of emodlib
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    project(emodlib-python)

    find_package(emodlib-malaria REQUIRED CONFIG)
endif()

# -- Other deps

find_package(pybind11 CONFIG REQUIRED)
message(STATUS "Found pybind11 v${pybind11_VERSION}: ${pybind11_INCLUDE_DIRS}")

# Installation directory
# ======================

# -- Set the Python site-package installation directory
find_package(PythonInterp)
execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_python_lib; import os;print(get_python_lib())"
                OUTPUT_VARIABLE PY_SP_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)

set(SP_DIR "${PY_SP_DIR}" CACHE PATH "Site-package directory")
message("-- CMake install prefix: " ${CMAKE_INSTALL_PREFIX})
message("  -> site-packages: " ${SP_DIR})

# Build
# =====

set(PY_EMODLIB _emodlib_py)

# -- emodlib
set(EMODLIB_OBJECTS
  ${CMAKE_SOURCE_DIR}/include/emodlib/malaria/InfectionMalaria.cpp
  ${CMAKE_SOURCE_DIR}/include/emodlib/malaria/SusceptibilityMalaria.cpp
  ${CMAKE_SOURCE_DIR}/include/emodlib/malaria/IntrahostComponent.cpp
  ${CMAKE_SOURCE_DIR}/include/emodlib/malaria/MalariaAntibody.cpp
  ${CMAKE_SOURCE_DIR}/include/emodlib/utils/RANDOM.cpp
  ${CMAKE_SOURCE_DIR}/include/emodlib/utils/suids
)

pybind11_add_module(${PY_EMODLIB} src/main.cpp ${EMODLIB_OBJECTS})

# include top-level directory to allow full include paths
target_include_directories(${PY_EMODLIB} PUBLIC ${CMAKE_SOURCE_DIR}/include)

#target_link_libraries(${PY_EMODLIB} PUBLIC emodlib-malaria)

target_compile_features(${PY_EMODLIB} PUBLIC cxx_std_14)

# Installation
# ============

install(TARGETS ${PY_FASTSCAPELIB}
        DESTINATION ${SP_DIR})
