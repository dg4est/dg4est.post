!
! File:   POD_snap.F90
! Author: akirby
!
! Created on January 23, 2024, 2:10 PM
!
subroutine POD_snap_covariance(pod_var,nfield,npt,dU,Y_data) bind(C,name="POD_snap_covariance_")
    use my_kinddefs
    use iso_c_binding
    implicit none

    integer(i4),intent(in) :: pod_var
    integer(i4),intent(in) :: nfield
    integer(i4),intent(in) :: npt
    real(dp),   intent(in) :: dU(nfield,npt)
    real(dp),  intent(out) :: Y_data(npt)

    integer(i4) :: i

    !==============!
    ! Snapshot POD !
    !==============!
!$OMP PARALLEL DO private(i)
    do i = 1,npt
        Y_data(i) = dU(pod_var,i)
    end do
!$OMP END PARALLEL DO
end subroutine

subroutine POD_snap_solve(npt,nstep,nstepMone,oneOnstepMone,KT,Y,R_tot,eigval,energy_total) bind(c,name="POD_snap_solve_")
    use my_kinddefs
    use timer_module
    use iso_c_binding
    use vars_analysis_module
    implicit none

    !!!!!!!!!!!!!!!!!!!!!! NOTE !!!!!!!!!!!!!!!!!!!!!
    ! This function should be only called PER field !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    integer(i4),intent(in) :: npt
    integer(i4),intent(in) :: nstep
    real(dp),   intent(in) :: nstepMone
    real(dp),   intent(in) :: oneOnstepMone
    real(dp),intent(inout) :: KT(nstep,nstep)
    real(dp),intent(inout) :: Y(npt,nstep)
    real(dp),intent(inout) :: R_tot(npt,nstep)
    real(dp),  intent(out) :: eigval(nstep)
    real(dp),  intent(out) :: energy_total

    character(1):: transA,transB
    integer(i4) :: lwork,info,i,j
    integer(i4) :: M,N,K,LDA,LDB,LDC
    real(dp)    :: alpha,beta
    real(dp)    :: w_t(4)

    real(dp) :: check1,check2,check3,check4,check5
    real(dp) :: check1d,check2d,check3d,check4d,check5d
    real(dp) :: t1,t2

    !========================================== NOTE ===========================================!
    ! This function computes the eigenvalues and eigenvectors of the Covariance Matrix.         !
    ! We use the BLAS (Basic Linear Algebra Subroutines) package to perform the computation.    !
    !                                                                                           !
    ! DSYEV computes all eigenvalues and, optionally, eigenvectors of a real symmetric matrix A.!
    ! POD Ref: https://sebastianraschka.com/Articles/2015_pca_in_3_steps.html#covariance-matrix !
    !===========================================================================================!
    ! DSYEV Ref: https://netlib.org/lapack/explore-html-3.6.1/d2/d8a/group__double_s_yeigen_ga442c43fca5493590f8f26cf42fed4044.html

    !==============!
    ! Snapshot POD !
    !==============!
    !------------------!
    ! Form K_transpose !
    !------------------!
    transA = 'T'
    transB = 'N'
    alpha = oneOnstepMone !non-biasing for discrete variance
    beta = 0.0_dp
    M = nstep      ! rows    of op(A) => rows    of Y^T -- Y^T(nstep,npt)
    N = nstep      ! columns of op(B) => columns of Y   -- Y  (npt,nstep)
    K = npt        ! columns of op(A) => columns of Y^T -- Y^T(nstep,npt)
    LDA = npt      ! Lead dimension of A = Y
    LDB = npt      ! Lead dimension of B = Y
    LDC = nstep    ! K_transpose (nstep,nstep)

    !form KT = oneOntimesMone * Y^T * Y = C(nstep,nstep)
    call wtime(t1)
    call dgemm(transA, & ! TRANSA: op(A) ['T'] transpose A, ['N'] no transpose of A.
               transB, & ! TRANSB: op(B) ['T'] transpose B, ['N'] no transpose of B.
               M,      & ! M: number  of rows  of the  matrix op(A) and of the matrix C.
               N,      & ! N: number  of columns of the matrix op(B) and the number of columns of the matrix C.
               K,      & ! K: number of columns of the matrix op(A) and the number of rows of the matrix op(B)
               alpha,  & ! ALPHA specifies the scalar alpha: C := alpha*op( A )*op( B ) + beta*C,
               Y,      & ! A: Matrix array.
               LDA,    & ! LDA: specifies the first dimension of A.
               Y,      & ! B: Matrix array.
               LDB,    & ! LDB: specifies the first dimension of B.
               beta,   & ! BETA: specifies the scalar  beta: C := alpha*op( A )*op( B ) + beta*C.
               KT,     & ! C: Output Matrix array.
               LDC)      ! LDC: specifies the first dimension of C.
    call wtime(t2)
    print*,"[Wake Analysis]  Matrix Multiply KT - Time (seconds): ",t2-t1

    if(.not.work_alloc)then
        allocate(work(1))
        work_alloc = .true.
    end if

    !-----------------------------!
    ! Find Optimal Workspace Size !
    !-----------------------------!
    call wtime(t1)
    lwork = -1 ! workspace query is assumed;
    call dsyev('V',   & ! JOBZ: ['N'] Compute eigenvalues only; ['V'] Compute eigenvalues and eigenvectors.
               'L',   & ! UPLO: ['U'] Upper triangle of A is stored; ['L'] Lower triangle of A is stored.
               nstep, & ! The order of the matrix A.  N >= 0.
               KT,    & ! The Matrix A, dimension (LDA, npt).
               nstep, & ! LDA is INTEGER. The leading dimension of the array A.  LDA >= max(1,N).
               eigval,& ! W is DOUBLE PRECISION array, dimension (N). If INFO = 0, the eigenvalues in ascending order.
               w_t,   & ! WORK dbl prec array: dimension (MAX(1,LWORK)). On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
               lwork, & ! LWORK = -1: then workspace query is assumed; the routine only calculates the optimal size of the WORK array.
               INFO)    ! INFO: [0] successful exit; [< 0] i-th argument illegal value; [> 0] algorithm failed to converge.

    lwork = int(w_t(1)) ! optimal size stored in first entry of WORK array
    if(work_size < lwork)then
        write(*,*) "New lwork: ", lwork,"MB:",lwork*sizeof(real(dp))/1024./1024.
        deallocate(work)
        allocate(work(lwork))
        work_size = lwork
    end if

    !--------------------------!
    ! Solve Eigenvalue Problem !
    !--------------------------!
    write(*,"(XA)") "[POD] Solving Eigenvalues/vectors..."
    call dsyev('V',   & ! JOBZ: ['N'] Compute eigenvalues only; ['V'] Compute eigenvalues and eigenvectors.
               'L',   & ! UPLO: ['U'] Upper triangle of A is stored; ['L'] Lower triangle of A is stored.
               nstep, & ! The order of the matrix A.  N >= 0.
               KT,    & ! The Matrix A, dimension (LDA, npt).
               nstep, & ! LDA is INTEGER. The leading dimension of the array A.  LDA >= max(1,N).
               eigval,& ! W is DOUBLE PRECISION array, dimension (N). If INFO = 0, the eigenvalues in ascending order.
               work,  & ! Work array
               lwork, & ! The length of the array WORK.
               info)    ! INFO: [0] successful exit; [< 0] i-th argument illegal value; [> 0] algorithm failed to converge.

    ! check for errors
    if (info /= 0 ) then
        write(*,*) "error in DSYEV"
        stop
    end if
    call wtime(t2)
    print*,"[Wake Analysis] DGYEV eig/vec - Time (seconds): ",t2-t1

    !------------!
    ! Form R_tot !
    !------------!
    transA = 'N'
    transB = 'N'
    alpha = 1.0_dp/sqrt(nstepMone) !non-biasing for discrete variance
    beta = 0.0_dp
    M = npt       ! rows    of op(A) => rows    of Y
    N = nstep     ! columns of op(B) => columns of K_transpose
    K = nstep     ! columns of op(A) => columns of Y
    LDA = npt     ! Lead dimension of A = Y
    LDB = nstep   ! Lead dimension of B = K_transpose
    LDC = npt     ! R_tot (npt,nstep)

    !form R_tot = oneOntimesMone * Y * KT(eigenvectors) = R_tot(nstep,nstep)
    call wtime(t1)
    call dgemm(transA,transB,M,N,K,alpha,Y,LDA,KT,LDB,beta,R_tot,LDC)
    call wtime(t2)
    print*,"[Wake Analysis] Matrix Multiply R_TOT- Time (seconds): ",t2-t1

    !-------------------------------------!
    ! Rescale to get correct eigenvectors !
    !-------------------------------------!
!$OMP PARALLEL DO private(i)
    do i = 1,nstep
        R_tot(:,i) = R_tot(:,i)/sqrt(abs(eigval(i)))
    end do
!$OMP END PARALLEL DO

    energy_total = sum(eigval)
    eigval = eigval/energy_total

    write(*,"(XA,F8.3)") "[POD] Energy Total: ",energy_total
    do i = nstep,nstep-9,-1
        write(*,"(7XA,I2,A,F6.3,A)") "Mode [",nstep-i+1,"]: ",eigval(i)*100.0,"%"
    end do

    !-----------------!
    ! Verify Solution !
    !-----------------!
#ifdef DEBUG
        !check if eigenvector dotted with it self is 1.0
        Check1  = dot_product(R_tot(:,nstep),R_tot(:,nstep))
        Check2  = dot_product(R_tot(:,nstep-2),R_tot(:,nstep-2))
        Check3  = dot_product(R_tot(:,nstep-3),R_tot(:,nstep-3))
        Check4  = dot_product(R_tot(:,nstep-4),R_tot(:,nstep-4))
        Check5  = dot_product(R_tot(:,nstep-5),R_tot(:,nstep-5))

        !check if eigenvector dotted diff eigenvector is orthogonal
        Check1d = dot_product(R_tot(:,nstep-1),R_tot(:,nstep-2))
        Check2d = dot_product(R_tot(:,nstep-2),R_tot(:,nstep-3))
        Check3d = dot_product(R_tot(:,nstep-3),R_tot(:,nstep-4))
        Check4d = dot_product(R_tot(:,nstep-4),R_tot(:,nstep-5))
        Check5d = dot_product(R_tot(:,nstep-5),R_tot(:,nstep-6))

        write(*,*) "POD check"

        write(*,*) "check1: ", check1
        write(*,*) "check2: ", check2
        write(*,*) "check3: ", check3
        write(*,*) "check4: ", check4
        write(*,*) "check5: ", check5

        write(*,*) "check1d: ", check1d
        write(*,*) "check2d: ", check2d
        write(*,*) "check3d: ", check3d
        write(*,*) "check4d: ", check4d
        write(*,*) "check5d: ", check5d
        write(*,*)
#endif

    ! clean up memory
    !if(work_alloc)then
    !    deallocate(work)
    !    work_alloc = .false.
    !endif

    !do i = 0,nstep-1
    !    write(*,*) maxval(abs(R_tot(:,npt-i)) - abs(R_tot(:,nstep-i)))
    !end do
end subroutine

subroutine POD_snap_reconstruct_solution(iextract,istep,pod_var,nfield,plotdim,npt,nx,ny,nz,nstep,nmode,POD_plot_nmodes,     &
                                         energy_total,eig,coordinates,U,Ubar,dU,R_tot,R_a_p,TVC,mode,UROM,plot_TVC,plot_ROM) &
    bind(c,name="POD_snap_reconstruct_solution_")
    use my_kinddefs
    implicit none

    integer(i4),intent(in) :: iextract
    integer(i4),intent(in) :: istep
    integer(i4),intent(in) :: pod_var
    integer(i4),intent(in) :: nfield
    integer(i4),intent(in) :: plotdim
    integer(i4),intent(in) :: npt
    integer(i4),intent(in) :: nx
    integer(i4),intent(in) :: ny
    integer(i4),intent(in) :: nz
    integer(i4),intent(in) :: nstep
    integer(i4),intent(in) :: nmode
    integer(i4),intent(in) :: POD_plot_nmodes
    real(dp),   intent(in) :: energy_total
    real(dp),   intent(in) :: eig(npt)
    real(dp),   intent(in) :: coordinates(3,npt)
    real(dp),   intent(in) :: U(nfield,npt)
    real(dp),   intent(in) :: Ubar(nfield,npt)
    real(dp),   intent(in) :: dU(nfield,npt)
    real(dp),   intent(in) :: R_tot(npt,nstep)
    real(dp),intent(inout) :: R_a_p(npt)
    real(dp),  intent(out) :: TVC(nmode)
    real(dp),  intent(out) :: mode(npt)
    real(dp),  intent(out) :: UROM(npt)
    integer(i4),intent(in) :: plot_TVC
    integer(i4),intent(in) :: plot_ROM

    !integer(i4),intent(in) :: POD_plot_timeVcoeffs

    integer(i4) :: i,j,n,iMode
    integer(i4) :: plot_count
    integer(i4) :: modeStop
    integer(i4) :: idiff(1)
    real(dp)    :: fdiff
    real(dp)    :: t1,t2
    real(dp)    :: perc_modes

    character(80) :: filename
    character(80) :: variables
    character(40) :: file_label
    character(3)  :: extract_num_str
    character(5)  :: file_num_str
    character(4)  :: mode_num_str
    character(5)  :: field_str
    character(20) :: temp

    Character(20) :: FMT

    ! calculate mode end mode index
    modeStop = nstep - nmode + 1

    !field
    write(extract_num_str,'(I3.3)') iextract
    write(file_num_str,'(I5.5)') istep

    select case(pod_var)
        case(1)
            field_str = "RHO"
            variables = 'VARIABLES = "x", "y", "z", "RHO"'
        case(2)
            field_str = "U"
            variables = adjustr(trim('VARIABLES = "x", "y", "z", "U"'))
        case(3)
            field_str = "V"
            variables = adjustr(trim('VARIABLES = "x", "y", "z", "V"'))
        case(4)
            field_str = "W"
            variables = adjustr(trim('VARIABLES = "x", "y", "z", "W"'))
        case(5)
            field_str = "P"
            variables = adjustr(trim('VARIABLES = "x", "y", "z", "P"'))
        case default
            field_str = "?"
            variables = adjustr(trim('VARIABLES = "x", "y", "z", "?"'))
    end select

    !reconstruct all energies at this time state
    !time varying modal coefficients: a_p(nmode)
    !NOTE: These are the time varying coefficients because
    !      we are projecting onto the i-th (in time) instantaneous
    !      flow solution. Plot the first m-coefficients in time by
    !      appending the time varying coefficients per mode to a
    !      file for each instantaneous flow solution.

    !==============!
    ! Snapshot POD !
    !==============!
    R_a_p(:) = 0.0_dp
    do j = 1,nstep
    do i = 1,npt
        R_a_p(j) = R_a_p(j) + R_tot(i,j)*dU(pod_var,i)
    end do
    end do

    !---------------------------!
    ! time varying coefficients !
    !---------------------------!
    if(plot_TVC==1)then
        iMode = 1
        do j = nstep,modeStop,-1 !NOTE: reverse order do to sorting of eigenvectors from BLAS
            TVC(iMode) = R_a_p(j) !MERGE((R_a_p(j)/sqrt((eig(j)*energy_total))),0.0_dp, abs(eig(j)) > 1E-12)
            iMode=iMode+1
        end do
    end if

    !---------------------!
    ! flow reconstruction !
    !---------------------!
    if(plot_ROM==1)then
        plot_count=0
        UROM(:) = 0.0_dp
        do j = nstep,modeStop,-1
            mode(:) = 0.0_dp
            do i = 1,npt
                mode(i) = mode(i) + R_tot(i,j)*R_a_p(j)
                UROM(i) = UROM(i) + mode(i)
            end do

            ! POD output(mode) : plot current mode
            if(plot_count < POD_plot_nmodes)then
                write(mode_num_str,'(I0)') nstep-j+1

                filename = adjustr(trim("analysis/tec/pod/modes/extract"))  // &
                           adjustr(trim(extract_num_str))                   // &
                           adjustr(trim("_"))                               // &
                           adjustr(trim(file_num_str))                      // &
                           adjustr(trim(".POD.mode"))                       // &
                           adjustr(trim(mode_num_str))                      // &
                           adjustr(trim("."))                               // &
                           adjustr(trim(field_str))                         // &
                           adjustr(trim(".tec"))

                call POD_snap_output(filename, &
                                     variables,&
                                     plotdim,npt,nx,ny,nz,    &
                                     mode,coordinates)
                plot_count=plot_count+1
            endif
        end do

        ! add mean-flow for final reconstruction
!$OMP PARALLEL DO private(i)
        do i = 1,npt
            UROM(i) = UROM(i) + Ubar(pod_var,i)
        end do
!$OMP END PARALLEL DO

        ! compute percent-error formula
        fdiff = maxval(abs((UROM - U(pod_var,:))))

        perc_modes = real(nmode,dp)/(real(nstep,dp))*100.0
        WRITE(temp,'(F6.2)') perc_modes
        WRITE(*,"(A,ES12.6,A,A,A,I0,A,I0,A)") ", Max Error: ",fdiff,&
                                               "  (",trim(adjustl(temp)),"% of modes [",nmode,"/",nstep,"])"
    end if
end subroutine

subroutine POD_snap_output(filename,variables,plotdim,npt,nx,ny,nz,U,coordinates)
    use my_kinddefs
    implicit none

    character(80),intent(in) :: filename
    character(80),intent(in) :: variables
    integer(i4),  intent(in) :: plotdim
    integer(i4),  intent(in) :: npt
    integer(i4),  intent(in) :: nx
    integer(i4),  intent(in) :: ny
    integer(i4),  intent(in) :: nz
    real(dp),     intent(in) :: U(npt)
    real(dp),     intent(in) :: coordinates(3,npt)

    integer(i4)     :: i
    integer(i4)     :: iunit = 1111

   !print*,">>> FILE: ",adjustr(trim(filename))," ",adjustr(trim(variables))," ",npt,nx,ny,nz
    open(unit=iunit,file=adjustr(trim(filename)),status='replace',form='formatted')
    write(iunit,"(A)") 'TITLE = "Wake Analysis POD"'
    write(iunit,"(A)") adjustr(trim(variables))

    if(plotdim==1) write(iunit,"(A22,I0,A4)") 'ZONE T="Only Zone", I=',nx,', F=POINT'
    if(plotdim==2) write(iunit,"(A22,I0,A4,I0,A9)") 'ZONE T="Only Zone", I=',nx,', J=',ny,', F=POINT'
    if(plotdim==3) write(iunit,"(A22,I0,A4,I0,A9,I0,A9)") 'ZONE T="Only Zone", I=',nx,', J=',ny,', K=',nz,', F=POINT'

    do i = 1,npt
        write(iunit,"(4(ES15.8,X))") coordinates(1:3,i),U(i)
    end do
    close(iunit)
end subroutine