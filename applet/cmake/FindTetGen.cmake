set(TETGEN_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/externals/tetgen")
include_directories(${TETGEN_INCLUDE_DIR})
add_subdirectory(${TETGEN_INCLUDE_DIR})
target_link_libraries(tetgen ${LIBRARIES})