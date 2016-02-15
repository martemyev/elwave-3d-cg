cmake_minimum_required(VERSION 2.8)


project(elwave-3d-cg CXX)


if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build: Debug | Release | RelWithDebInfo | PDebug | PRelease" FORCE)
endif()
string(TOUPPER ${CMAKE_BUILD_TYPE} BUILD_TYPE)


set(AUX_CXX_FLAGS "-std=c++0x -Wall") # use this with MFEM
if(NOT MFEM_DIR)
  set(MFEM_DIR "$ENV{HOME}/projects/mfem" CACHE PATH "Path to the MFEM root directory" FORCE)
endif()
message(STATUS "MFEM_DIR ${MFEM_DIR}")
if(NOT MFEM_BUILD)
  if(BUILD_TYPE STREQUAL "DEBUG")
    set(MFEM_BUILD "serialdebug" CACHE STRING
        "Name of MFEM build (it's a subdirectory with config file and library")
  elseif(BUILD_TYPE STREQUAL "RELWITHDEBINFO")
    set(MFEM_BUILD "serial")
  elseif(BUILD_TYPE STREQUAL "RELEASE")
    set(MFEM_BUILD "serial")
  elseif(BUILD_TYPE STREQUAL "PDEBUG")
    set(MFEM_BUILD "paralleldebug")
  elseif(BUILD_TYPE STREQUAL "PRELEASE")
    set(MFEM_BUILD "parallel")
  else()
    message(FATAL_ERROR "Unknown BUILD_TYPE")
  endif()
endif()
message(STATUS "MFEM_BUILD ${MFEM_BUILD}")
include_directories("${MFEM_DIR}") # for "mfem.hpp"
include_directories("${MFEM_DIR}/${MFEM_BUILD}") # for "config.hpp"
link_directories("${MFEM_DIR}/${MFEM_BUILD}")


if(BUILD_TYPE STREQUAL "PDEBUG" OR BUILD_TYPE STREQUAL "PRELEASE")
  find_package(MPI REQUIRED)
  include_directories(${MPI_CXX_INCLUDE_PATH})
  #set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
  # Hypre
  find_path(HYPRE_DIR include/HYPRE.h HINTS ENV HYPRE_DIR PATHS $ENV{HOME}/hypre DOC "Hypre Directory")
  if(EXISTS ${HYPRE_DIR}/include/HYPRE.h)
    set(HYPRE_FOUND YES)
    message(STATUS "Hypre found: ${HYPRE_DIR}")
    set(HYPRE_INCLUDES ${HYPRE_DIR})
    find_path(HYPRE_INCLUDE_DIR HYPRE.h HINTS "${HYPRE_DIR}" PATH_SUFFIXES include NO_DEFAULT_PATH)
    list(APPEND HYPRE_INCLUDES ${HYPRE_INCLUDE_DIR})
    include_directories(${HYPRE_INCLUDES})
    file(GLOB HYPRE_LIBRARIES RELATIVE "${HYPRE_DIR}/lib" "${HYPRE_DIR}/lib/libHYPRE*.a")
    message(STATUS "Hypre lib: ${HYPRE_LIBRARIES}")
  else()
    message(FATAL_ERROR "Hypre is NOT found")
    set(HYPRE_FOUND NO)
  endif()
  # Metis
  find_path(METIS_DIR include/metis.h HINTS ENV METIS_DIR PATHS $ENV{HOME}/metis DOC "Metis Directory")
  if(EXISTS ${METIS_DIR}/include/metis.h)
    set(METIS_FOUND YES)
    message(STATUS "Metis found: ${METIS_DIR}")
    set(METIS_INCLUDES ${METIS_DIR})
    find_path(METIS_INCLUDE_DIR metis.h HINTS "${METIS_DIR}" PATH_SUFFIXES include NO_DEFAULT_PATH)
    list(APPEND METIS_INCLUDES ${METIS_INCLUDE_DIR})
    include_directories(${METIS_INCLUDES})
    file(GLOB METIS_LIBRARIES RELATIVE "${METIS_DIR}/lib" "${METIS_DIR}/lib/libmetis*.a")
    message(STATUS "Metis lib: ${METIS_LIBRARIES}")
    set(MFEM_USE_METIS_5 ON)
  else()
    message(FATAL_ERROR "Metis is NOT found")
    set(METIS_FOUND NO)
  endif()
endif()


if(NOT(BUILD_TYPE STREQUAL "RELEASE") AND
   NOT(BUILD_TYPE STREQUAL "PRELEASE"))
  if(NOT(${CMAKE_CXX_COMPILER_ID} STREQUAL "MSVC"))
    set(AUX_CXX_FLAGS "${AUX_CXX_FLAGS} -fno-inline")
  endif()
endif()

if(BUILD_TYPE STREQUAL "DEBUG" OR BUILD_TYPE STREQUAL "PDEBUG")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${AUX_CXX_FLAGS}")
elseif(BUILD_TYPE STREQUAL "RELWITHDEBINFO")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ${AUX_CXX_FLAGS}")
elseif(BUILD_TYPE STREQUAL "RELEASE" OR BUILD_TYPE STREQUAL "PRELEASE")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${AUX_CXX_FLAGS}")
endif()


file(GLOB SRC_LIST "${PROJECT_SOURCE_DIR}/src/*.cpp")
file(GLOB HDR_LIST "${PROJECT_SOURCE_DIR}/src/*.hpp")
include_directories("${PROJECT_SOURCE_DIR}/src") # for "config.hpp"


configure_file(
  "${PROJECT_SOURCE_DIR}/config_elwave.hpp.in"
  "${PROJECT_SOURCE_DIR}/src/config_elwave.hpp")


message("\n--------------- configuration information ----------------")
  message("build type     : " ${BUILD_TYPE})
  message("compiler       : " ${CMAKE_CXX_COMPILER})
  message("hypre dir      : " ${HYPRE_DIR})
  message("hypre lib      : " ${HYPRE_LIBRARIES})
  message("metis dir      : " ${METIS_DIR})
  message("metis lib      : " ${METIS_LIBRARIES})
if(BUILD_TYPE STREQUAL "DEBUG")
  message("compiler flags : " ${CMAKE_CXX_FLAGS_DEBUG})
elseif(BUILD_TYPE STREQUAL "RELWITHDEBINFO")
  message("compiler flags : " ${CMAKE_CXX_FLAGS_RELWITHDEBINFO})
elseif(BUILD_TYPE STREQUAL "RELEASE")
  message("compiler flags : " ${CMAKE_CXX_FLAGS_RELEASE})
endif()


add_executable(${PROJECT_NAME} ${SRC_LIST} ${HDR_LIST})
target_link_libraries(${PROJECT_NAME} mfem rt)
if(BUILD_TYPE STREQUAL "PDEBUG" OR BUILD_TYPE STREQUAL "PRELEASE")
  target_link_libraries(${PROJECT_NAME} ${MPI_CXX_LIBRARIES} HYPRE metis)
endif()

