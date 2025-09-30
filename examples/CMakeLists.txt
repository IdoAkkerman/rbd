#-------------------------------------------------------------------------------
# This file is part of the RBD application
#-------------------------------------------------------------------------------
add_custom_target(examples ALL)
add_custom_command(TARGET examples
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/examples
                   COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/examples/* ${CMAKE_BINARY_DIR}/examples)

add_test(NAME rbd_2d      COMMAND ../rbd rb_2d.json)
add_test(NAME rbd_matrix  COMMAND ../rbd rb_matrix.json)
add_test(NAME rbd_angles  COMMAND ../rbd rb_angles.json)
add_test(NAME rbd_radians COMMAND ../rbd rb_radians.json)
