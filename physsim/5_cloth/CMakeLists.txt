set(EXECUTABLE_NAME physsim_5_cloth)

# Add your project files
file(GLOB SRCFILES *.cpp)
file(GLOB HFILES *.hpp)

add_executable(${EXECUTABLE_NAME} ${SRCFILES} ${HFILES})
target_link_libraries(${EXECUTABLE_NAME} PRIVATE physsim_common)
target_include_directories(${EXECUTABLE_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
set_target_properties(${EXECUTABLE_NAME} PROPERTIES FOLDER "physsim")

include(../../cmake/generate_init.cmake)
generate_init_function(${EXECUTABLE_NAME})
