get_filename_component(EXERCISENAME ${CMAKE_CURRENT_LIST_DIR} NAME)

file(GLOB_RECURSE SOURCES "*.cpp")
file(GLOB_RECURSE HEADERS "*.h")
file(GLOB_RECURSE SHADERS "*.glsl")

add_executable(${EXERCISENAME} ${SOURCES} ${HEADERS} ${SHADERS} Edge.h PDConstraint.h Quaternion.cpp Quaternion.h)
target_link_libraries(${EXERCISENAME} surface_mesh)
target_link_libraries(${EXERCISENAME} nanogui ${NANOGUI_EXTRA_LIBS})

if(NOT MSVC)
    find_package(OpenMP REQUIRED)
    target_link_libraries(${EXERCISENAME} ${OpenMP_CXX_LIBRARIES})
endif()