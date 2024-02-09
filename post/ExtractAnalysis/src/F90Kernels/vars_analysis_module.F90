module vars_analysis_module
    use my_kinddefs
    implicit none

    real(dp),pointer :: work(:)
    integer(i4)      :: work_size = 1
    logical          :: work_alloc = .false.
end module