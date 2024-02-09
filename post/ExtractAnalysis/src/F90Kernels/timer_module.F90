!
! File:   timer_module.F90
! Author: akirby
!
! Created on January 24, 2024, 7:50 PM
!
module timer_module
    implicit none
    contains

    subroutine wtime(t)
        use my_kinddefs
        implicit none

        real(rp),intent(out) :: t

        integer(i4) :: clock_max
        integer(i4) :: clock_rate
        integer(i4) :: clock_reading

        call system_clock(clock_reading,clock_rate,clock_max)

        t = real(clock_reading,rp) &
          / real(clock_rate,rp)
    end subroutine
end module