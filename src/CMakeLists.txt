#-------------------------------------------------------------------------------
# This file is part of the RBD application
#-------------------------------------------------------------------------------
set(SOURCES
   vector.cpp
   matrix.cpp
   state.cpp
   rigid-body.cpp)

set(HEADERS
   vector.hpp
   matrix.hpp
   state.hpp
   rigid-body.hpp)

add_library(rbd ${SOURCES} ${HEADERS})
target_link_libraries(rbd PRIVATE ${LIBRARIES})
