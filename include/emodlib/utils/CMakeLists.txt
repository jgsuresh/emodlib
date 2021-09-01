set(EMODLIB_UTILS_TARGET emodlib-utils)

add_library(${EMODLIB_UTILS_TARGET} STATIC suids.cpp RANDOM.cpp)

target_compile_features(${EMODLIB_UTILS_TARGET} PUBLIC cxx_std_14)
