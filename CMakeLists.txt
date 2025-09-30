#-------------------------------------------------------------------------------
# This file is part of the RBD application
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
add_subdirectory(src)

add_executable(precice_rbd src/precice_rbd.cpp ${SOURCES} ${HEADERS})
target_link_libraries(precice_rbd PRIVATE rbd ${LIBRARIES})
install(TARGETS precice_rbd DESTINATION bin)

add_executable(test_rbd src/test_rbd.cpp ${SOURCES} ${HEADERS})
target_link_libraries(test_rbd PRIVATE rbd ${LIBRARIES})
install(TARGETS test_rbd DESTINATION bin)

#-------------------------------------------------------------------------------
# Add examples
#-------------------------------------------------------------------------------
enable_testing()
add_subdirectory(examples)


