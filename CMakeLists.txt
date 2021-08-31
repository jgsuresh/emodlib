cmake_minimum_required(VERSION 3.1)

project(emodlib)

set(EMODLIB_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

OPTION(BUILD_TESTS "build emodlib test suite" OFF)
OPTION(BUILD_BENCHMARK "build emodlib benchmark suite" OFF)
OPTION(BUILD_PYTHON_MODULE "build emodlib Python bindings" ON)

set(CMAKE_MODULE_PATH
    ${CMAKE_SOURCE_DIR}/cmake
    ${CMAKE_MODULE_PATH})

# Versioning
# ==========

include(VersionUtils)
set_version_str("${EMODLIB_INCLUDE_DIR}")
message(STATUS "Building emodlib v${VERSION_STR}")

# Configure/Build
# ===============

# -- emodlib

add_subdirectory(include/emodlib)

# -- optional subdirectories

if(BUILD_TESTS)
  enable_testing()
  add_subdirectory(test)
endif()

if(BUILD_BENCHMARK)
  add_subdirectory(benchmark)
endif()

if(BUILD_PYTHON_MODULE)
  add_subdirectory(python)
endif()
