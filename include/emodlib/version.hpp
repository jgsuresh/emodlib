#pragma once

#define EMODLIB_VERSION_MAJOR 0
#define EMODLIB_VERSION_MINOR 1
#define EMODLIB_VERSION_PATCH 0

#include <string>

namespace emodlib
{
    namespace version
    {
        constexpr int version_major = EMODLIB_VERSION_MAJOR;
        constexpr int version_minor = EMODLIB_VERSION_MINOR;
        constexpr int version_patch = EMODLIB_VERSION_PATCH;
        static const std::string version_str = std::to_string(version_major) + "." +
                                               std::to_string(version_minor) + "." +
                                               std::to_string(version_patch);     
    }
}
