cmake_minimum_required(VERSION 2.8)

project(elwave-3d-cg CXX)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build: Debug | Release | RelWithDebInfo" FORCE)
  set(DEBUG_WAVE OFF)
endif()
string(TOUPPER ${CMAKE_BUILD_TYPE} BUILD_TYPE)

if(BUILD_TYPE MATCHES DEBUG)
  set(DEBUG_WAVE ON)
endif()

set(MFEM_DIR "/u/artemyev/projects/mfem/")
include_directories("${MFEM_DIR}")
link_directories("${MFEM_DIR}/build")

file(GLOB SRC_LIST "${PROJECT_SOURCE_DIR}/*.cpp")
file(GLOB HDR_LIST "${PROJECT_SOURCE_DIR}/*.hpp")

configure_file(
  "${PROJECT_SOURCE_DIR}/config.hpp.in"
  "${PROJECT_SOURCE_DIR}/config.hpp")

add_executable(${PROJECT_NAME} ${SRC_LIST} ${HDR_LIST})
target_link_libraries(${PROJECT_NAME} mfem rt)

