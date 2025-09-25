#-------------------------------------------------------------------------------
# This file is part of the RBVMS application.
#
# This is the superbuild file.
#
# RBVMS is free software; you can redistribute it and/or modify it under the
# terms of the BSD-3 license. We welcome feedback and contributions, see file
# CONTRIBUTING.md for details.
#-------------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.8...4.99)
project("Rigid-Body Dynamics" LANGUAGES CXX)

#-------------------------------------------------------------------------------
# Prohibit in-source builds
#-------------------------------------------------------------------------------
if ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message(FATAL_ERROR "In-source builds are prohibited.")
endif ()

#-------------------------------------------------------------------------------
# Defaults
#-------------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#-------------------------------------------------------------------------------
# Add Precice
#-------------------------------------------------------------------------------
find_package(precice REQUIRED
             HINTS "${PRECICE_DIR}" "${PRECICE_DIR}/lib/cmake/precice")
list(APPEND ${LIBRARIES} "precice::precice")
message(STATUS "Found precice config in: ${precice_DIR} (version ${precice_VERSION})")

get_target_property(PRECICE_INC_PATH  "precice::precice" INTERFACE_INCLUDE_DIRECTORIES)
include_directories(${PRECICE_INC_PATH})
set(LIBRARIES "precice::precice")

#-------------------------------------------------------------------------------
# Add Json
#-------------------------------------------------------------------------------
find_package(jsoncpp REQUIRED)
get_target_property(JSON_INC_PATH jsoncpp_lib INTERFACE_INCLUDE_DIRECTORIES)
include_directories(${JSON_INC_PATH})
list(APPEND LIBRARIES "jsoncpp_lib")

#-------------------------------------------------------------------------------
# Build executable
#-------------------------------------------------------------------------------
add_executable(rbd rbd.cpp linalg.hpp  rigid-body.hpp)
target_link_libraries(rbd PRIVATE  ${LIBRARIES})
install(TARGETS rbd DESTINATION bin)

#-------------------------------------------------------------------------------
# Add test
#-------------------------------------------------------------------------------
configure_file(${CMAKE_SOURCE_DIR}/rb.json ${CMAKE_BINARY_DIR}/rb.json COPYONLY)
configure_file(${CMAKE_SOURCE_DIR}/rb_matrix.json ${CMAKE_BINARY_DIR}/rb_matrix.json COPYONLY)
configure_file(${CMAKE_SOURCE_DIR}/rb_angles.json ${CMAKE_BINARY_DIR}/rb_angles.json COPYONLY)
configure_file(${CMAKE_SOURCE_DIR}/rb_radians.json ${CMAKE_BINARY_DIR}/rb_radians.json COPYONLY)
enable_testing()
add_test(NAME rbd         COMMAND rbd rb.json)
add_test(NAME rbd_matrix  COMMAND rbd rb_matrix.json)
add_test(NAME rbd_angles  COMMAND rbd rb_angles.json)
add_test(NAME rbd_radians COMMAND rbd rb_radians.json)


