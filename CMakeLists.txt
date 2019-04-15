#top dir cmake project for hackrf_hunt

cmake_minimum_required(VERSION 2.8)
project (HackRF C)

set(CMAKE_C_FLAGS "$ENV{CFLAGS}" CACHE STRING "C Flags")

add_subdirectory(hackrf-hunt)

########################################################################
# Create uninstall target
########################################################################

configure_file(
    ${PROJECT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake
@ONLY)


add_custom_target(uninstall
    ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake
)
