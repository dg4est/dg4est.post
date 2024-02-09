!***************************************************
!  Module my_kinddefs:
!     Specifies machine independent precision of
!     double and single precision real variables
!     and 4 types of integer variables.  The reals
!     are made accessible through variables sp and dp
!     while the integers are i1, i2 i4, i8, each of
!     which has that many digits behind the decimal place
!**********************************************************

module my_kinddefs
    implicit none

    public rp, sp, dp, cp, i1, i2, i4, i8, i16
    public onethird,twothird,fourthird

    private
    !---> All these a based on the default IEEE 754 format

    !---> Single Precision Real: IEEE 754
    integer, parameter :: sp=selected_real_kind(6,37)

    !---> Double Precision Real: IEEE 754
    integer, parameter :: dp=selected_real_kind(15,307)

    !---> Double Precision Complex: IEEE 754
    integer, parameter :: cp=selected_real_kind(15,307)

    !---> Quad Precision Real: IEEE
    integer, parameter :: qp=selected_real_kind(34,4931)

    !---> Single decimal place integer (1 byte)
    integer, parameter :: i1=selected_int_kind(2)

    !---> 2 decimal place integer (2 bytes):
    integer, parameter :: i2=selected_int_kind(3)

    !---> 4 decimal place integer ( 4 bytes):
    integer, parameter :: i4=selected_int_kind(5)

    !---> 8 decimal place integer (8 bytes):
    integer, parameter :: i8=selected_int_kind(10)

    !---> 12 decimal place integer:
    integer,parameter :: i16=selected_int_kind(16)

#ifdef SINGLE_PRECISION
    !---> Single Precision Real: IEEE 754
    integer, parameter :: rp=selected_real_kind(6,37)
#else
    !---> Double Precision Real: IEEE 754
    integer, parameter :: rp=selected_real_kind(15,307)
#endif

    real(rp),parameter :: onethird  = 1.0_rp / 3.0_rp
    real(rp),parameter :: twothird  = 2.0_rp / 3.0_rp
    real(rp),parameter :: fourthird = 4.0_rp / 3.0_rp
end module
