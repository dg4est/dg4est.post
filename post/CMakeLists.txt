CMAKE_MINIMUM_REQUIRED(VERSION 3.1)
PROJECT(Post.Flow)

# --------------------------------------------------------------------------- #
enable_language(Fortran)
include_directories(${MPI_INCLUDE_PATH})
# --------------------------------------------------------------------------- #
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
# --------------------------------------------------------------------------- #
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# --------------------------------------------------------------------------- #
# set C++17 standard for OCCA compatibility
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# --------------------------------------------------------------------------- #
# force GNU linker to include all shared libraries
if (NOT APPLE)
  set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-as-needed")
else ()
  set(CMAKE_MACOSX_RPATH ON)
endif ()

set(occa_dir  "" CACHE PATH "occa  install directory location")
# --------------------------------------------------------------------------- #
# set OCCA directory
set(ENV{OCCA_DIR} "${occa_dir}")
message("OCCA_DIR: " $ENV{OCCA_DIR})
# --------------------------------------------------------------------------- #
option(ENABLE_OpenMP OFF)
option(ENABLE_CUDA OFF)
option(ENABLE_HIP OFF)

if (ENABLE_OpenMP)
  find_package(OpenMP REQUIRED)
endif()

if (OPENMP_FOUND)
  message(">>> Enabling OpenMP!")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS}; -fopenmp)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
# --------------------------------------------------------------------------- #
# ====================== #
# CUDA/HIP Specific Code #
# ====================== #
set(cuda_flag 0)
if(ENABLE_CUDA)
  set(cuda_flag 1)
endif()

set(hip_flag 0)
if(ENABLE_HIP)
  set(hip_flag 1)
endif()

if((${cuda_flag} EQUAL ${hip_flag}) AND (${cuda_flag} EQUAL 1))
  message(FATAL_ERROR "[ERROR] Cannot have CUDA and HIP specifices enabled at the same time! Set one to off.")
endif()

if(ENABLE_CUDA)
  find_package(CUDA REQUIRED)
  include_directories(${CUDA_INCLUDE_DIRS})
  add_definitions(-D__DG_CUDA__)
  message(">>> Enabling CUDA specifics...")
endif()

if(ENABLE_HIP)
  find_package(HIP REQUIRED)
  include_directories(${HIP_INCLUDE_DIRS})
  add_definitions(-D__DG_HIP__)

  execute_process (
    COMMAND bash -c "hipconfig --platform"
    OUTPUT_VARIABLE hip_platform
  )
  message(">>> Enabling HIP specifics... " ${hip_platform})
  if(${hip_platform} MATCHES "hcc")
    add_definitions(-D__HIP_PLATFORM_HCC__)
  elseif(${hip_platform} MATCHES "nvcc")
    add_definitions(-D__HIP_PLATFORM_NVCC__)
  endif()
endif()
# --------------------------------------------------------------------------- #
# ======================================= #
if (CMAKE_C_COMPILER_ID MATCHES "Intel")
  # Intel optimization flag
  add_definitions( -xCore-AVX2 )
endif ()

# Library header files to be included
include_directories($ENV{OCCA_DIR}/include)

# Libraries to be linked
find_library(occa_lb NAMES occa PATHS $ENV{OCCA_DIR}/lib)
# --------------------------------------------------------------------------- #
# ======================================= #
# Analysis Source
add_subdirectory(ExtractAnalysis)
# ======================================= #
