//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//

#if !defined(KRATOS_FILESYSTEM)
#define  KRATOS_FILESYSTEM

// System / External includes
#if defined(__cplusplus) && __cplusplus >= 201703L && defined(__has_include) && __has_include(<filesystem>)
    #include <filesystem>
    #define KRATOS_FILESYSTEM_AVAILABLE
#else
    //---------------------------------------------------------------------------------------
    // fs_fwd.hpp - The forwarding header for the header/implementation seperated usage of
    //              ghc::filesystem.
    // This file can be include at any place, where ghc::filesystem api is needed while
    // not bleeding implementation details (e.g. system includes) into the global namespace,
    // as long as one cpp includes fs_impl.hpp to deliver the matching implementations.
    //---------------------------------------------------------------------------------------
    #include "ghc/fs_fwd.hpp"
#endif

namespace Kratos {

// alias to make it available as Kratos::filesystem
#ifdef KRATOS_FILESYSTEM_AVAILABLE
    namespace filesystem = std::filesystem;
#else
    namespace filesystem = ghc::filesystem;
#endif
}

#endif // KRATOS_FILESYSTEM defined