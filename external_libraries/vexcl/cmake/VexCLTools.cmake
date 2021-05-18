# Distributed under the MIT License.
# Copyright (c) 2012-2018 Denis Demidov <dennis.demidov@gmail.com>

function(vexcl_add_executables targetname)
    add_library(${targetname} INTERFACE)

    if(TARGET VexCL::OpenCL)
        add_executable(${targetname}_cl ${ARGN})
        target_link_libraries(${targetname}_cl PUBLIC VexCL::OpenCL ${targetname})
        set_target_properties(${targetname}_cl PROPERTIES CXX_EXTENSIONS OFF)
    endif()
    if(TARGET VexCL::Compute)
        add_executable(${targetname}_comp ${ARGN})
        target_link_libraries(${targetname}_comp PUBLIC VexCL::Compute ${targetname})
        set_target_properties(${targetname}_comp PROPERTIES CXX_EXTENSIONS OFF)
    endif()
    if(TARGET VexCL::CUDA)
        add_executable(${targetname}_cuda ${ARGN})
        target_link_libraries(${targetname}_cuda PUBLIC VexCL::CUDA ${targetname})
        set_target_properties(${targetname}_cuda PROPERTIES CXX_EXTENSIONS OFF)
    endif()
    if(TARGET VexCL::JIT)
        add_executable(${targetname}_jit ${ARGN})
        target_link_libraries(${targetname}_jit PUBLIC VexCL::JIT ${targetname})
        set_target_properties(${targetname}_jit PROPERTIES CXX_EXTENSIONS OFF)
    endif()
endfunction()

function(vexcl_add_libraries targetname typeoflib)
    add_library(${targetname} INTERFACE)
    add_library(${targetname}::Common ALIAS ${targetname})

    if(TARGET VexCL::OpenCL)
        add_library(${targetname}_cl ${TYPEOFLIB} ${ARGN})
        target_link_libraries(${targetname}_cl PUBLIC VexCL::OpenCL ${targetname})
        set_target_properties(${targetname}_cl PROPERTIES CXX_EXTENSIONS OFF)
        add_library(${targetname}::OpenCl ALIAS ${targetname}_cl)
    endif()
    if(TARGET VexCL::Compute)
        add_library(${targetname}_comp ${TYPEOFLIB} ${ARGN})
        target_link_libraries(${targetname}_comp PUBLIC VexCL::Compute ${targetname})
        set_target_properties(${targetname}_comp PROPERTIES CXX_EXTENSIONS OFF)
        add_library(${targetname}::Compute ALIAS ${targetname}_comp)
    endif()
    if(TARGET VexCL::CUDA)
        add_library(${targetname}_cuda ${TYPEOFLIB} ${ARGN})
        target_link_libraries(${targetname}_cuda PUBLIC VexCL::CUDA ${targetname})
        set_target_properties(${targetname}_cuda PROPERTIES CXX_EXTENSIONS OFF)
        add_library(${targetname}::CUDA ALIAS ${targetname}_cuda)
    endif()
    if(TARGET VexCL::JIT)
        add_library(${targetname}_jit ${TYPEOFLIB} ${ARGN})
        target_link_libraries(${targetname}_jit PUBLIC VexCL::JIT ${targetname})
        set_target_properties(${targetname}_jit PROPERTIES CXX_EXTENSIONS OFF)
        add_library(${targetname}::JIT ALIAS ${targetname}_jit)
    endif()
endfunction()

