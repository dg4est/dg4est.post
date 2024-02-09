# set C++17 standard for smart pointers
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")

# --------------------------------------------------------------------------- #
# load BLAS package
find_package(BLAS REQUIRED)
message(STATUS "BLAS found: ${BLAS_LIBRARIES}")

option(ENABLE_OpenMP ON)
if (ENABLE_OpenMP)
  find_package(OpenMP)
endif()

if (OPENMP_FOUND)
  message(">>> Enabling OpenMP!")
  message(STATUS "OpenMP found: ${OPENMP_LIBRARIES}")
 #set(CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS}; -fopenmp)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# --------------------------------------------------------------------------- #
# force GNU linker to include all shared libraries
if (NOT APPLE)
  set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-as-needed")
else ()
  set(CMAKE_MACOSX_RPATH ON)
endif ()

find_library(occa_lb   NAMES occa   PATHS $ENV{OCCA_DIR}/lib)
# --------------------------------------------------------------------------- #

# header file includes
include_directories(include)
include_directories(src)
include_directories($ENV{OCCA_DIR}/include)

set(SRC
    src/analysis_read_inputs.cxx
    src/analysis_binary2ascii.cxx
    src/analysis_temporal_mean.cxx
    src/analysis_reynolds_stresses.cxx
    src/analysis_les_subgridscales.cxx
    src/analysis_pod_output.cxx
    src/analysis_pod_traditional.cxx
    src/analysis_pod_snapshot.cxx
    main.cxx
)

# Fortran POD library
add_subdirectory(src/F90Kernels)

add_executable(PostFlow.exe ${SRC})
target_link_libraries(PostFlow.exe math-kernels ${occa_lb} m)

install(TARGETS PostFlow.exe
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib/
        ARCHIVE DESTINATION lib/
)
