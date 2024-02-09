# source files
set(F90_KERNELS
    my_kinddefs.F90
    timer_module.F90
    vars_analysis_module.F90
    POD_trad.F90
    POD_snap.F90
)

# build shared library
add_library(math-kernels SHARED ${F90_KERNELS})

if(BLAS_FOUND)
    include_directories(${BLAS_INCLUDE_DIRS})
    target_link_libraries(math-kernels ${BLAS_LIBRARIES})
endif(BLAS_FOUND)

# install math_kernels
install(TARGETS math-kernels
        RUNTIME DESTINATION bin/
        LIBRARY DESTINATION lib/
        ARCHIVE DESTINATION lib/
)