cmake_minimum_required(VERSION 3.1)

function(set_version_str INCLUDEDIR)
    file(STRINGS "${INCLUDEDIR}/emodlib/version.hpp" emod_version_defines
        REGEX "#define EMODLIB_VERSION_(MAJOR|MINOR|PATCH)")
    foreach(ver ${emod_version_defines})
        if(ver MATCHES "#define EMODLIB_VERSION_(MAJOR|MINOR|PATCH) +([^ ]+)$")
            set(EMODLIB_VERSION_${CMAKE_MATCH_1} "${CMAKE_MATCH_2}" CACHE INTERNAL "")
        endif()
    endforeach()
    set(VERSION_STR
        ${EMODLIB_VERSION_MAJOR}.${EMODLIB_VERSION_MINOR}.${EMODLIB_VERSION_PATCH} PARENT_SCOPE)
endfunction()
