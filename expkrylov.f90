module expkrylov
! ************************************************************************************************
! This is the source code of the Exponential Integrator with Schur-Krylov Approximation.
!
! Developed by 
!     Zaigang LIU             
!         @ Institute of Engineering Thermophysics, Chinese Academy of Sciences, Beijing 100190, China
!     Jean-Louis CONSALVI     
!         @ IUSTI UMR 7343, CNRS, Aix-Marseille Université, Marseille 13453, France
!     Wenjun KONG             
!         @ Institute of Engineering Thermophysics, Chinese Academy of Sciences, Beijing 100190, China
!
! Please cite this reference:
!     Zaigang Liu, Jean-L. Consalvi, Wenjun Kong,
!     An Exponential Integrator with Schur–Krylov Approximation to accelerate combustion chemistry 
!     computation,
!     Combustion and Flame,
!     Volume 203,
!     2019,
!     Pages 180-189,
!     ISSN 0010-2180,
!     https://doi.org/10.1016/j.combustflame.2019.01.031.
!     (http://www.sciencedirect.com/science/article/pii/S0010218019300495)
!     Abstract: The Exponential Integrator with Schur–Krylov Approximation (EISKA) algorithm was 
!     developed for combustion applications. This algorithm combines the advantages of the explicit 
!     large step advancement of the exponential schemes and the dimension reduction effect of the 
!     Krylov subspace approximation, and was improved by introducing the Schur decomposition to control 
!     the rounding error. The EISKA based on the SpeedCHEM (SC) package was implemented to simulate a 
!     methane partially stirred reactor (PaSR) with pair-wise mixing model by considering the 
!     mechanisms of Li et al., GRI-Mech 3.0 and USC Mech II. Accuracy and computational efficiency 
!     of EISKA are systematically compared with those of DVODE. In the case of the Li mechanism which 
!     is a priori sufficiently small to be handled directly in combustion simulations, the 
!     computations were accelerated by a factor of 1.99 without losing accuracy. In the cases of 
!     GRI-Mech 3.0 and USC Mech II which are significantly larger than the Li mechanism, chemical 
!     reduction methods, namely the Correlated Dynamic Adaptive Chemistry (CoDAC) and the 
!     Multi-timescale (MTS) method were coupled with either DVODE or EISKA. The results show that 
!     the EISKA is faster than DVODE either with or without chemical reduction methods. Model results 
!     show that the best strategy is to use EISKA without any reduction method which leads to the 
!     same accuracy as compared to DVODE and acceleration factors of 2.61 and 2.19 for GRI-Mech 3.0 
!     and USC Mech II, respectively.
!     Keywords: Chemistry computation acceleration; Exponential integrator; Krylov subspace; Mechanism reduction
!
! Please follow updates or send comments through the web page <https://github.com/codepublzg/expintschurkrylov.git.>
! Published online @2019.02.23 by Zaigang Liu.
! ************************************************************************************************
!
! Usage:
!
!     Step 1. Invoke this module:
!
!         use expkrylov
!
!     Step 2. Create a option data type, e.g.:
!
!         type(expkry_opts) :: expkry_opt
!
!     Step 3. Set options at a proper place in your code, e.g.:
!
!         expkry_opt = expkry_setopts ( atol = 1d-8, rtol = 1d-4, ...)
!
!     Step 4. define external subroutines:
!
!         * Perform Matrix-Vector Multiplication: y = A * x
!         subroutine matvec (n, x, y)
!             implicit none
!             integer :: n
!             real(WP),dimension(n),intent(in)  :: x
!             real(WP),dimension(n),intent(out) :: y
!         end subroutine matvec
!
!         * Update reaction rate vector and, if necessary, please
!             1)  check the mass fraction is realistic, and 
!             2)  enforce conservation of the mass fraction.
!         subroutine updateyw (n, phiw, y, w, ierr)
!             ! update 
!             implicit none
!             integer :: n, ierr
!             real(WP),dimension(n),intent(in)    :: phiw
!             real(WP),dimension(n),intent(inout) :: y
!             real(WP),dimension(n),intent(out)   :: w
!         end subroutine updateyw
!
!         * Compute Jacobian matrix
!         subroutine obtainJac (method, n, y, L)
!             implicit none
!             integer :: method, n
!             real(WP),dimension(n),intent(in)    :: y
!             real(WP),dimension(n,n),intent(out),optional :: L
!         end subroutine obtainJac
!
!     Step 5. Voilà:
!
!         call expkry_solve (n, h, y, lwork, rwrk, expkry_opt, matvec, updateyw, obtainJac, status, info)
!
!     For more details please read the code and its comments, or contact the authors.
! ************************************************************************************************
    use expmA
    implicit none
    private

    integer, parameter :: WP = kind(1.0D0)
    real(WP), parameter :: eps = epsilon(1D0)

    type expkry_opts
        integer :: krylovDimThreshold
        integer :: MaxKrylovSize
        integer :: MinKrylovSize
        integer :: Maxsubsteps
        integer :: MaxcheckEigmax
        integer :: MaxReject
        integer :: JACinterv
        integer :: ioutfile
        real(WP) :: atol
        real(WP) :: rtol
        real(WP) :: uround
        real(WP) :: break_tol
        real(WP) :: gamma
        real(WP) :: delta
        logical :: SchurDecomp
        logical :: ComputSchurQuotient
        logical :: InfoOutput
        logical :: ErrGiveUp
        logical :: reorth
        logical :: expkrylov_init = .false.
        
    end type expkry_opts
    
    public :: expkry_opts, expkry_setopts, expkry_lwork, expkry_solve

    contains

    !#########################################
    ! usage: option = setopts(...)
    function expkry_setopts(   &
        krylovDimThreshold,  &
        MaxKrylovSize,       &
        MinKrylovSize,       &
        Maxsubsteps,         &
        MaxcheckEigmax,      &
        MaxReject,           &
        JACinterv,           &
        ioutfile,            &
        atol, rtol,          &
        uround, break_tol,   &
        gamma, delta,        &
        SchurDecomp,         &
        ComputSchurQuotient, &
        ErrGiveUp,           &
        reorth,              &
        InfoOutput           &
        ) result(option)
    implicit none

    type(expkry_opts) :: option
    integer, optional :: krylovDimThreshold, MaxKrylovSize, MinKrylovSize, MaxReject
    integer, optional :: Maxsubsteps, MaxcheckEigmax, JACinterv, ioutfile
    real(WP), optional :: atol, rtol, break_tol, gamma, delta, uround
    logical, optional :: SchurDecomp, ComputSchurQuotient, ErrGiveUp, reorth, InfoOutput

    !if (present(krylovDimThreshold)) then
    !    option%krylovDimThreshold = krylovDimThreshold
    !else
        option%krylovDimThreshold = 1
    !end if
    
    if (present(MaxKrylovSize)) then
        option%MaxKrylovSize = MaxKrylovSize
    else
        option%MaxKrylovSize = 200
    end if
    
    if (present(MinKrylovSize)) then
        option%MinKrylovSize = MinKrylovSize
    else
        option%MinKrylovSize = 10
    end if
    
    if (present(Maxsubsteps)) then
        option%Maxsubsteps = Maxsubsteps
    else
        option%Maxsubsteps = 1000
    end if
    
    if (present(MaxcheckEigmax)) then
        option%MaxcheckEigmax = MaxcheckEigmax
    else
        option%MaxcheckEigmax = 5
    end if
    
    if (present(MaxReject)) then
        option%MaxReject = MaxReject
    else
        option%MaxReject = 100
    end if
    
    if (present(JACinterv)) then
        option%JACinterv = JACinterv
    else
        option%JACinterv = 50
    end if       
    
    if(present(uround)) then
        option%uround = uround
    else
        option%uround = 1d-15
    end if           
    
    if(present(atol)) then
        option%atol = atol
    else
        option%atol = 1d-8
    end if
    
    option%atol = max(option%uround,option%atol)
    
    if(present(rtol)) then
        option%rtol = rtol
    else
        option%rtol = 1d-4
    end if 
    
    option%rtol = max(option%rtol,option%atol)

    if(present(break_tol)) then
        option%break_tol = break_tol
    else
        option%break_tol = 1d-8 !option%rtol
    end if
    
    if(present(gamma)) then
        option%gamma = gamma
    else
        option%gamma = 0.8
    end if 
    
    if(present(delta)) then
        option%delta = delta
    else
        option%delta = 1.2
    end if 
    
    if(present(SchurDecomp)) then
        option%SchurDecomp = SchurDecomp
    else
        option%SchurDecomp = .true.
    end if  
    
    if(present(ComputSchurQuotient)) then
        option%ComputSchurQuotient = ComputSchurQuotient &
                                    .and. option%SchurDecomp
    else
        option%ComputSchurQuotient = .true. .and. option%SchurDecomp
    end if  
    
    
    if(present(ErrGiveUp)) then
        option%ErrGiveUp = ErrGiveUp
    else
        option%ErrGiveUp = .false.
    end if  

    if(present(reorth)) then
        option%reorth = reorth
    else
        option%reorth = .false.
    end if  
    
    if(present(InfoOutput)) then
        option%InfoOutput = InfoOutput
    else
        option%InfoOutput = .false.
    end if  

    if (present(ioutfile)) then
        option%ioutfile = ioutfile
    else
        option%ioutfile = 6
    end if
    
    option%expkrylov_init = .true.
    
    end function expkry_setopts
    !#########################################

    !#########################################
    function expkry_lwork(n, option)
    implicit none
    integer :: expkry_lwork, n, p, ldm
    type(expkry_opts) :: option
    
    if (n .lt. option%krylovDimThreshold) then
        if (option%SchurDecomp) then
            expkry_lwork = 12*n*n + 98*n + 4167
        else
            expkry_lwork = 9*n*n + 20*n + 7
        end if
    else
        p = 2
        ldm = min(option%MaxKrylovSize,n-1) + p
        if (option%SchurDecomp) then
            expkry_lwork = 83*ldm + 3*n + ldm*n + ldm*p + 12*ldm*ldm + 4160
        else
            expkry_lwork = 4*ldm + 3*n + ldm*n + ldm*p + 9*ldm*ldm
        end if        
    end if

    end function
    !#########################################        
        
    !#########################################
    function MV_blas (n, M, vec)
    implicit none
    integer :: n
    real(WP), dimension(n,n) :: M
    real(WP), dimension(n)   :: vec, MV_blas
    
    call DGEMV ('N', n, n, 1d0, M, n, vec, 1, 0d0, MV_blas, 1)
    
    end function MV_blas
    !#########################################

    !#########################################
    function MTV_blas (n, M, vec)
    implicit none
    integer :: n
    real(WP), dimension(n,n) :: M
    real(WP), dimension(n)   :: vec, MTV_blas
    
    call DGEMV ('T', n, n, 1d0, M, n, vec, 1, 0d0, MTV_blas, 1)
    
    end function MTV_blas
    !#########################################
    
    !#########################################
    subroutine rcond(N, A, lda, condest, anorm, lwork, work)
    ! workspace : (N+4)*N
    ! flops : (2/3)*n**3 + 2*n*n
    implicit none
    integer :: lwork, N, lda
    real(WP) :: condest, anorm, A(lda,N)
    real(WP) :: work(lwork)

    real(WP),dimension(:,:),pointer :: B
    integer :: IPIV(N), IWORK(N), INFO, iworkNow
    
    ! check input
    if (lwork .lt. (N+4)*N) then
        write(*,"('subroutine rcond: not enough workspace:', 2I8)") lwork,  (N+4)*N
        stop
    end if
    
    iworkNow = 1
    B => assignWorkspace(iworkNow, lwork, work, N, N) ! + N*N
    
    !B = A
    call DLACPY( 'A', n, n, A, lda, B, n) !
    anorm = Mat1norm(n, n, B, n)
    call DGETRF( N, N, B, N, IPIV, INFO )
    call DGECON( '1', N, B, N, anorm, condest, WORK(iworkNow), IWORK, INFO)
    
    if (INFO .ne. 0) stop 'Error in subroutine rcond.'
    
    B => nullWorkspace(iworkNow, N, N, 'B')
    
    end subroutine rcond
    !#########################################
    
    !#########################################
    function tolyra (n, y, m, option)
    implicit none
    integer :: n,i, m
    real(WP) :: y(n), tolyra(n)
    type(expkry_opts),intent(in) :: option
    
    !uround = dble(m**3)*eps
    do i = 1, n
        tolyra(i) = abs(y(i))*option%rtol + option%atol !, option%uround)
    end do    
    
    end function tolyra
    !######################################### 
    
    !#########################################
    function setTimestepLimitNorm (n, Norm, tmax, tol, uround) result(tlimit)
    implicit none
    integer :: n
    real(WP):: Norm, tmax, tol, tlimit, uround
    
    !uround = dble(n**3)*eps
    tlimit = min(tol/uround/Norm,tmax)
    
    end function setTimestepLimitNorm
    !#########################################
    
    !#########################################
    function setTimestepLimitEigmax (n, eigmax, tmax, tol, uround) result(tlimit)
    implicit none
    integer :: n
    real(WP):: eigmax, tmax, tol, tlimit, uround
    
    if (eigmax .gt. 0) then
        !uround = dble(n**3)*eps
        tlimit = min(log(tol/uround)/eigmax, tmax)
    else
        tlimit = tmax
    end if
    
    end function setTimestepLimitEigmax
    !#########################################
    
    !#########################################
    subroutine buildMbar (n, t, M, ldm, vec, Mbar, p, lhess)
    ! flops: n^2
    implicit none
    integer :: n, p, ldm
    real(WP) :: t
    real(WP), dimension(ldm,n),intent(in) :: M
    real(WP), dimension(n),intent(in) :: vec
    real(WP), dimension(n+p,n+p),intent(out) :: Mbar  
    logical :: lhess
    integer :: info
    
    integer :: i,j
    
    call DLACPY( 'A', n, n, M, ldm, Mbar, n+p) ! Mbar(1:n,1:n) = M
    call DCOPY(n,vec,1,Mbar(1,n+1),1) ! Mbar(1:n,n+1) = vec
    do j = n+2, n+p
        Mbar(1:n,j) = 0.0_WP
    end do
    do i = n+1, n+p 
        Mbar(i, 1:n+p) = 0.0_WP
    end do
    do i = 2, p
        Mbar(n+i-1, n+i) = 1.0_WP
    end do

    ! Mbar = t*Mbar
    if (lhess) then
        call DLASCL( 'H', 0, 0, 1.0_WP, t, n+p, n+p, Mbar, n+p, INFO )
    else
        call DLASCL( 'G', 0, 0, 1.0_WP, t, n+p, n+p, Mbar, n+p, INFO )
    end if

    end subroutine buildMbar
    !#########################################
    
    !#########################################
    subroutine solvePhiw (n, tstep, L, ldl, w, p, phiwn, iscale, lphiwrk, phiwrk, option)
    !
    ! lphiwrk = 7*(n+p)*(n+p)
    ! flops: 4*(n+p) + 2*(n+p)^3*s + 27*(n+p)^2 + (44*(n+p)^3)/3 = (n+p)^2 + 4*(n+p) + 2*(n+p)^3*s + 26*(n+p)^2 + (44*(n+p)^3)/3
    implicit none
    type(expkry_opts) :: option
    integer :: n, lphiwrk, p, ldl
    real(WP) :: tstep
    real(WP), dimension(ldl,n), intent(in) :: L
    real(WP), dimension(n), intent(in) :: w
    real(WP), dimension(n,p), intent(out)  :: phiwn
    real(WP), dimension(lphiwrk) :: phiwrk
    
    integer :: iposition, iscale
    real(WP), dimension(:,:),pointer :: Lbar, expmLbar
    
    integer :: nlexpm, iflag
    real(WP), dimension(:),pointer :: expmwrk
    type(expmA_opts) :: expmOption
    logical :: lhess
    
    character(len=*),parameter ::   &
        fmt_expm  = "('       Expm Success, dimension : ', I8)",   & ! n+p
        fmt_lhess = "('       - Hessenberg in expm ? ', L8)",   & !
        fmt_iscal = "('       - Scaling in expm    : ', I8)"      !
    
    ! check input
    if (lphiwrk .lt. 7*(n+p)*(n+p)) stop 'solvePhiw : not enough workspcae.'
    
    iposition = 1
    if (option%ComputSchurQuotient) then
        expmOption = expmSetOptions(SchurQuotient = .true.)
        lhess = .true.
    else
        expmOption =expmSetOptions(performSchur = .false.)
        lhess = .false.
    end if
    
    Lbar => assignWorkspace(iposition, lphiwrk, phiwrk, n+p, n+p)      ! lwork: + (n+p)*(n+p)
    expmLbar =>  assignWorkspace(iposition, lphiwrk, phiwrk, n+p, n+p) ! lwork: + (n+p)*(n+p)
    
    nlexpm = 5*(n+p)*(n+p)
    expmwrk =>  assignWorkspace(iposition, lphiwrk, phiwrk, nlexpm) ! lwork: + nlexpm
    
    call buildMbar (n, tstep, L, ldl, w, Lbar, p, lhess) ! flops: (n+p)^2
    
    !call print_dense_to_file(L, ldl, n, 'L')
    !call print_dense_to_file(Lbar, n+p, n+p, 'Lbar')

    ! flops: 4*(n+p) + 2*(n+p)^3*s + 26*(n+p)^2 + (44*(n+p)^3)/3
    call dexpm_pade (n+p, Lbar, expmLbar, iscale, nlexpm, expmwrk, expmOption, iflag)
    if (iflag .ne. 0) stop 'solvePhiw : Error in dexpm_pade.'
    
    if (option%InfoOutput) then
        write(option%ioutfile,fmt_expm) n+p
        write(option%ioutfile,fmt_lhess) lhess
        write(option%ioutfile,fmt_iscal) iscale
    end if
    
    !call print_dense_to_file(expmLbar, n+p, n+p, 'expmLbar')
    
    call DLACPY( 'A', n, p, expmLbar(1,n+1), n+p, phiwn, n) ! phiwn(1:n,1:p) = expmLbar(1:n, n+1:n+p)
    
    expmwrk  => nullWorkspace(iposition, nlexpm, 'expmwrk')  
    expmLbar => nullWorkspace(iposition, n+p, n+p, 'expmLbar')  
    Lbar     => nullWorkspace(iposition, n+p, n+p, 'Lbar')  

    end subroutine solvePhiw
    !#########################################
    
    
    !#########################################
    subroutine expkry_solve (n, h, y, lwork, rwrk, option, &
        matvec, updateyw, obtainJac, status, info)

    ! solve the problem
    !   y = y0 + phi1(h*L)*(h*w)
    !   y(n), w(n), L(n,n)
    !
    ! Please provide subroutines:
    !
    !    subroutine matvec (n, x, y)
    !        implicit none
    !        integer :: n
    !        real(WP),dimension(n),intent(in)  :: x
    !        real(WP),dimension(n),intent(out) :: y
    !    end subroutine matvec
    !
    !    subroutine updateyw (n, phiw, y, w, ierr)
    !        implicit none
    !        integer :: n, ierr
    !        real(WP),dimension(n),intent(in)    :: phiw
    !        real(WP),dimension(n),intent(inout) :: y
    !        real(WP),dimension(n),intent(out)   :: w
    !    end subroutine updateyw
    !
    !    subroutine obtainJac (method, n, y, L)
    !        implicit none
    !        integer :: method, n
    !        real(WP),dimension(n),intent(in)    :: y
    !        real(WP),dimension(n,n),intent(out),optional :: L
    !    end subroutine obtainJac
    !
    implicit none
    
    integer :: n, info
    real(WP) :: h
    real(WP),dimension(n) :: y, ybk, w, phiw
    real(WP),dimension(1,1) :: L
    integer :: lwork, ierr
    real(WP),dimension(lwork) :: rwrk
    type(expkry_opts) :: option
    integer,dimension(3) :: status
    external :: matvec, updateyw, obtainJac
    logical :: FileIsOpened
    
    if (.not. option%expkrylov_init) stop 'Please call expkry_setopts first.'
    
    info = 0
    status = 0
    ybk = y
    
    if(option%InfoOutput) then
        if (option%ioutfile .ne. 6) then
            inquire(option%ioutfile, opened=FileIsOpened)
            if (.not.FileIsOpened) then
                open(option%ioutfile, file='screenOutput')
            end if
        end if
    end if
    
    if (n .eq. 1) then
        phiw = 0.0_WP
        call updateyw (n, phiw, y, w, ierr)
        call obtainJac (1, n, y, L)
        y(1) = y(1) + ((exp(h*L(1,1))-1.0_WP)/L(1,1))*w(1)  
    ! elseif (n .lt. option%krylovDimThreshold) then
    ! *** Deprecated for low efficiency ***
    !     if(option%InfoOutput) then
    !         write(option%ioutfile,"('############ Using directSolve...  ############')")
    !         write(option%ioutfile,"('Dimension: ', I10)") n
    !         write(option%ioutfile,"('Time step: ', ES10.3)") h
    !     end if
    !     call directSolve(n, h, y, lwork, rwrk, option, updateyw, obtainJac, status, info)
    !     if(info .ne. 0 ) then
    !         if (option%InfoOutput) write(option%ioutfile,"('directSolve failed. Change to KrylovSolve ... ')")
    !         y = ybk
    !     end if
    !     if(option%InfoOutput) write(option%ioutfile,"('++++++++++++ directSolve finished. ++++++++++++')")
    end if
    
    if ((n .gt. option%krylovDimThreshold) .or. info .ne. 0) then
        if(option%InfoOutput) then
            write(option%ioutfile,"('############ Using KrylovSolve...  ############')")
            write(option%ioutfile,"('Dimension: ', I10)") n
            write(option%ioutfile,"('Time step: ', ES10.3)") h
        end if
        call KrylovSolve(n, h, y, lwork, rwrk, option, matvec, updateyw, obtainJac, status, info)
        if(option%InfoOutput) write(option%ioutfile,"('++++++++++++ KrylovSolve finished. ++++++++++++')")
    end if
        
    if (info .ne. 0) then
        if(option%InfoOutput) then
            write(option%ioutfile,"('+++++++++!!! expkry_solve Fail. !!!++++++++')")
        end if
        y = 1.0_WP/eps
    end if
    
    if(option%InfoOutput) then
        if (.not. FileIsOpened) then
            close(option%ioutfile)
        end if
    end if
    
    end subroutine expkry_solve
    !#########################################
        
    !#########################################  
    ! *** Deprecated for low efficiency ***  
    ! subroutine directSolve (n, h, y, lwork, rwrk, option, updateyw, obtainJac, status, info)  
    ! ! info = 1 -- try Krylov
    ! ! info = 2 -- failure
    ! implicit none
    ! integer :: n, info
    ! integer,dimension(3) :: status
    ! real(WP) :: h
    ! real(WP),dimension(n) :: y
    ! integer :: lwork
    ! real(WP),dimension(lwork) :: rwrk
    ! type(expkry_opts) :: option        
    ! external :: updateyw, obtainJac
    
   
    ! end subroutine directSolve
    !#########################################   

    !#########################################   
    subroutine Arnoldi(n, b, mStart, mEnd, matvec, option, V, H, ldm, beta, hm1m)
    ! flops: 3*n + n*(mEnd - mStart)*(2*mEnd - 2*mStart + 2*n + 3) = 2*n + n + n*(mEnd - mStart)*(2*mEnd - 2*mStart + 2*n + 3)
    implicit none
    integer :: n, mStart, mEnd, m, ldm
    type(expkry_opts) :: option  
    real(WP),dimension(n),  intent(in) :: b
    real(WP),dimension(n,ldm+1) :: V
    real(WP),dimension(ldm,ldm) :: H
    real(WP), intent(out) :: beta, hm1m
    external :: matvec
    real(WP),external :: DDOT, DNRM2
    
    integer :: i, j
    
    ! check input
    if(mStart .gt. mEnd) then
        write(option%ioutfile,"('Arnoldi: mStart is greater than mEnd: ', I8, I8)") mStart, mEnd
        stop
    end if
    if (mEnd .gt. n-1) then
        write(option%ioutfile,"('Arnoldi: mEnd is greater than n-1: ', I8, I8)") mEnd, n-1
        stop
    end if
    
    ! initialize
    beta = DNRM2(n,b,1)   ! flops: 2*n
    if(beta .lt. option%uround) then
        mEnd = 0
        V = 0
        H = 0
        hm1m = 0
        return
    end if
    call DCOPY(n,b,1,V(1:n,1),1)
    call DSCAL(n,1.0d0/beta,V(1:n,1),1)  ! v1 = b / beta ! flops: n
    
    m = mEnd
    ! start looping
    ! flops: n*(mEnd - mStart)*(2*mEnd - 2*mStart + 2*n + 3) = (mEnd - mStart)*(2*n^2+ (mEnd - mStart)*2*n + 2*n + n)
    do j = mStart, m
        call matvec(n, V(1:n,j), V(1:n,j+1)) ! p = A*vj ! flops: 2*n^2
        ! flops(j = mStart:mEnd): j*4*n = (mEnd - mStart)*2*n
        do i = 1, j
            H(i,j) = DDOT(n,V(1:n,i),1,V(1:n,j+1),1)       ! hij = vT*p ! flops: 2*n
            call DAXPY(n,-H(i,j),V(1:n,i),1,V(1:n,j+1),1)   ! p = -hij*v + p ! flops: 2*n
        end do
        ! new hm1m
        hm1m = DNRM2(n,V(1:n,j+1),1) ! flops: 2*n
        call DSCAL(n,1.0d0/hm1m,V(1:n,j+1),1)  ! vj+1 = p / hj+1,j ! flops: n
        if(j .ne. m)  H(j+1,j) = hm1m
        ! if 'happy breakdown'?
        if (abs(hm1m) .le. option%break_tol) then
            mEnd = j
            exit
        end if
    end do
    
    end subroutine Arnoldi
    !#########################################

    
    !#########################################
    subroutine dcpArnoldi (n, v, mStart, mEnd, matvec, option, Q, H, ldm, beta, hm1m)
    
    ! The Arnoldi algorithm with complete reorthogonalization and partial pivoting
    ! Reference: 
    !    A.S. Hodel, P. Misra
    !    Partial pivoting in the computation of Krylov subspaces of large sparse systems
    !    42nd IEEE Conference on Decision and Control, IEEE, Maui, HI, USA (2003), 10.1109/CDC.2003.1273062
    !    paper 2878-2883. Print ISBN: 0-7803-7924-1; Print ISSN: 0191-2216

    ! n - dimension of square matrix A
    ! v - vector, dimension n, normalized, i.e., norm(v)=1
    integer :: n, mStart, mEnd, m, ldm
    type(expkry_opts) :: option  
    real(WP),dimension(n),  intent(in) :: v
    real(WP),dimension(n,ldm+1) :: Q
    real(WP),dimension(ldm,ldm) :: H
    real(WP) :: brktol, beta, hm1m
        
    external :: matvec

    ! operating variables
    integer :: i,j,k
    integer :: maxid(1)
    real(WP), dimension(n,n) :: U
    real(WP), dimension(n) :: alpha, w, wpp, tmpw, Up, signHistory
    real(WP) :: signTmp1, signTmp2
    real(WP) :: rtmp, zero, normRest
    integer :: iter, itmp, nRemain        
    integer,dimension(n) :: pivot_vec
    real(WP),external :: DNRM2


    if (n .le. 0) return

    !if (lrwk .lt. n*(n+6) ) then
    !write(*,*) 'Not enough space for dcpArnoldi:',lrwk,n*(n+6)
    !end if

    ! zero
    zero = option%uround
    
    beta = DNRM2(n,v,1)   ! flops: 2*n
    if(beta .lt. zero) then
        mEnd = 0
        Q = 0
        H = 0
        hm1m = 0
        return
    end if
    
    m = mEnd
    if (mStart .eq. 1) then
        ! 1. Set up vector of pivot points.
        ! pivot_vec(j): real position of j.th element in w
        do i = 1, n
            pivot_vec(i) = i
            w(i) = v(i)/beta
        end do

        Q(:,1) = w !-w

        ! Locate max magnitude element in w.
        tmpw = abs(w)
        maxid = maxloc(tmpw)

        ! See if need to change the pivot list.
        if (maxid(1) .ne. 1 ) then
            ! swap pivot_vec(maxid) and pivot_vec(1)
            pivot_vec(maxid(1)) = 1
            pivot_vec(1) = maxid(1)
        end if

        ! Householder reflection.
        wpp = w(pivot_vec)
        call housh(n, wpp, 1, Up, alpha(1), signHistory(1), zero)
        U(pivot_vec(1:n),1) = Up(:)
    else
        write(*,*) 'mStart must be 1 for reorthogonalization'
        stop
    end if

    ! loop
    ! start from 1, not mStart?
    do j = mStart, m
            
        nRemain = n-j
            
        call matvec(n, Q(1:n,j), w)

        ! Project off of previous vectors.
        signTmp1 = 1d0
        do i = 1, j
            Up(:) = U(:,i)
            rtmp = alpha(i)*dot_product(Up,w)
            w = rtmp*Up(:)-w !
            w = signHistory(i)*w
            signTmp1 = signTmp1*signHistory(i)
        end do

        signTmp2 = 1d0
        do i = 1, j
            signTmp2 = signTmp2*signHistory(i)
            H(i,j) = (-1)**(i+j)*signTmp1*signTmp2*w(pivot_vec(i))
        end do

        ! short index
        wpp(1:nRemain) = w(pivot_vec(j+1:n))

        ! norm of short

        normRest = 0d0
        do i = 1, nRemain
            normRest = normRest + wpp(i)*wpp(i)
        end do
        normRest = sqrt(normRest)

        ! find max
        tmpw(1:j) = 0d0
        tmpw(j+1:n) = abs(wpp(1:nRemain))
        maxid = maxloc(tmpw) ! in pivot_vec order

        ! swap
        if (maxid(1) .ne. j+1) then
            itmp = pivot_vec(j+1)
            pivot_vec(j+1) = pivot_vec(maxid(1))
            pivot_vec(maxid(1)) = itmp
        end if

        ! pivot_vec changed
        wpp(1:nRemain) = w(pivot_vec(j+1:n))

        !
        call housh(nRemain, wpp, 1, Up, alpha(j+1),signHistory(j+1),zero)
        U(pivot_vec(1:j),  j+1) = 0d0
        U(pivot_vec(j+1:n),j+1) = Up(1:nRemain)

        !
        H(j+1,j) = normRest !
        hm1m = normRest
      
        ! Construct next Q and multiply.
        Q(:,j+1) = 0d0
        Q(pivot_vec(j+1),j+1) = 1d0
        do i = j+1, 1, -1
            wpp(1:n-i+1) = Q(pivot_vec(i:n),j+1)
            Up(1:n-i+1) = U(pivot_vec(i:n),i)
            rtmp = alpha(i)*dot_product(Up(1:n-i+1), wpp(1:n-i+1))
            Q(pivot_vec(i:n),j+1) = rtmp*Up(1:n-i+1) - wpp(1:n-i+1) !
            Q(:,j+1) = signHistory(i)*Q(:,j+1)
        end do

        if ( abs(H(j+1,j)) .lt. option%break_tol) then
            ! break down
            mEnd = j
            exit
        end if

    end do 
    
    end subroutine dcpArnoldi
    !#########################################
    
    
    !#########################################
    subroutine housh (nx, x, nj, housv, tau, tauSign, zero)
    ! [housv, tau] = housholder(x, nj, zero)
    ! (I-tau*housv*housv')x = -sigma*ej
    integer :: nx, nj
    real(8), dimension(nx) :: x, housv
    real(8) :: zero, tau, tauSign

    integer :: nv
    real(8) :: ita, sigma
    real(WP),external :: DNRM2


    housv = x
    ita = maxval(abs(housv))
    if (ita .ne. 0d0) then
        housv = housv / ita
        !sigma = 0d0
        !do nv = 1, nx
        !    sigma = sigma + housv(nv)*housv(nv)
        !end do
        !sigma = sqrt(sigma)
        sigma = DNRM2( nx, housv, 1)
        
        if (sigma .gt. zero) then
            tau = 1.0d0 / (sigma * (sigma + abs(housv(nj))))
            tauSign = sign(1d0,housv(nj))
            housv(nj) = housv(nj) + tauSign*sigma
            ! sigma = sigma*ita
        else
            tau = 0d0
        end if
    else
        tau = 0d0
    end if    
    
    end subroutine housh
    !#########################################
    
    !#########################################
    subroutine AdaptNextStep(omega, omega_old, tstep, tstep_old, h, tnow, tlimit,  &
        m, m_old, mLowlimit, mUplimit, n, s, fixmFlag, mfixed, option, imode, factort, factorm)
    implicit none

    real(WP) :: omega, omega_old, tstep, tstep_old, h, tnow, tlimit, factort
    integer :: imode, m, m_old, n, s, mLowlimit, mUplimit, factorm, mfixed
    type(expkry_opts) :: option 
    logical :: fixmFlag

    real(WP) :: order, kest    
    
    real(WP) :: omegaFac, timeFac, tmpFac1, tmpFac2
    integer  :: mDiff, mOpt, mOptmid
    real(WP) :: timeOpt, timeOptmid
    
    real(WP) :: cost1, cost2, cost3, const
    integer :: p
    
    character(len=*),parameter :: &
        fmt_opt   = "('       - Optimized m and dt: ', I8, ES15.6)",               & ! mOpt, timeOpt
        fmt_mid   = "('       - Mid       m and dt: ', I8, ES15.6)",               & ! mOpt, timeOpt
        fmt_mfix  = "('       - m is fixed. time step will be changed.')",         & !
        fmt_modem = "('       - Maximum m reached, consider to increase mmax.')",  & !
        fmt_mode1 = "('       - Change time next step, costs are: ', 3ES15.6)",    & !
        fmt_mode2 = "('       - Change m    next step, costs are: ', 3ES15.6)",    & !
        fmt_mode3 = "('       - Change both next step, costs are: ', 3ES15.6)"       !
    
    !--- happy breakdown
    if (imode .eq. 3) return
    
    ! factors
    omegaFac = omega / omega_old
    timeFac  = tstep / tstep_old
    mDiff    = m_old - m
    
    !--- cases
    if (imode .eq. 0) then
        ! initial guess
        order = max(1.0_WP, 0.25_WP*m) + 1.0_WP
        kest  = 2.0_WP       
    elseif (imode .eq. 1) then
        ! change time       
        order = 1.0_WP + max(1.0_WP, dlog(omegaFac)/dlog(timeFac))
        kest  = 2.0_WP 
    elseif (imode .eq. 2) then
        ! change m
        if (mDiff .eq. 0) then
            write(option%ioutfile,"('AdaptNextStep : mDiff == 0')")
            stop
        end if  
        order = max(1.0_WP, 0.25_WP*m) + 1.0_WP
        kest  = max(1.1_WP, (omegaFac)**(1.0_WP/(mDiff)))  
    end if
    
    !--- optimized time and m
    
    tmpFac1 = omega/option%gamma
    
    tmpFac2 = tmpFac1**(-1.0_WP/order) ! minus
    tmpFac2 = max(0.2_WP, min(tmpFac2, 10.0_WP))
    timeOpt = tstep*tmpFac2
    timeOpt = min(timeOpt, h-tnow, tlimit)
    timeOpt = max(eps, timeOpt)
    timeOptmid = sqrt(tstep*timeOpt)
      
    mOpt = ceiling(dble(m)+log(tmpFac1)/log(kest))
    mOpt = max(mOpt, floor(3d0*dble(m)/4d0), mLowlimit)
    mOpt = min(mOpt, mUplimit, ceiling(4d0*dble(m)/3d0))
    mOptmid = (m+mOpt)/2
    
    if(option%InfoOutput) then
        write(option%ioutfile,fmt_opt) mOpt, timeOpt
        write(option%ioutfile,fmt_mid) mOptmid, timeOptmid
    end if
    
    !--- computational cost throughout the rest steps (flops)
    ! flops: 
    !                  |---------Arnoldi--------|  |------------schur---------|  |----- rcond  ----| |---------------------expm--------------------------| |----err&update------|
    !   schur         : 3*n + n*m*(2*m + 2*n + 3) + 10*m^3 + 2*m^3 + 2*m^2 + m + (2/3)*m^3 + 2*m*m   + 4*(m+p) + 2*(m+p)^3*s + 27*(m+p)^2 + (44*(m+p)^3)/3 + 5*n + 2*n^2 + update
    !                 = 5*m + 8*n + 4*p + 2*s*(m + p)^3 + 27*(m + p)^2 + (44*(m + p)^3)/3 + 4*m^2 + (38*m^3)/3 + 2*n^2 + m*n*(2*m + 2*n + 3) + update
    !   without schur : 3*n + n*m*(2*m + 2*n + 3)                              + (2/3)*m^3 + 2*m*m   + 4*(m+p) + 2*(m+p)^3*s + 27*(m+p)^2 + (44*(m+p)^3)/3 + 5*n + 2*n^2 + update
    !                 = 4*m + 8*n + 4*p + 2*s*(m + p)^3 + 27*(m + p)^2 + (44*(m + p)^3)/3 + 2*m^2 + (2*m^3)/3 + 2*n^2 + m*n*(2*m + 2*n + 3) + update
    p = 2
    const = dble(n)*dble(m)**3 ! *** A manually defined tuning constant based on experience. This maybe not the best choice.
    ! change time step 
    if (option%SchurDecomp) then
        cost1 = dble(5*m + 8*n + 4*p + 2*s*(m + p)**3 + 27*(m + p)**2 + &
            (44*(m + p)**3)/3 + 4*m**2 + (38*m**3)/3 + 2*n**2 + &
            m*n*(2*m + 2*n + 3)) + const
    else
        cost1 = dble(4*m + 8*n + 4*p + 2*s*(m + p)**3 + 27*(m + p)**2 + &
            (44*(m + p)**3)/3 + 2*m**2 + (2*m**3)/3 + 2*n**2 + &
            m*n*(2*m + 2*n + 3)) + const
    end if
    !cost1 = cost1 * dble(ceiling((h-tnow)/timeOpt))
    cost1 = cost1 * ((h-tnow)/timeOpt)
    if (cost1 .lt. 0.0_WP) then
        write(option%ioutfile,*) cost1, n, m, h-tnow, tstep
        stop 'cost1 < 0'
    end if
    
    ! change m
    if (option%SchurDecomp) then
        cost2 = dble(5*mOpt + 8*n + 4*p + 2*s*(mOpt + p)**3 + 27*(mOpt + p)**2 + &
            (44*(mOpt + p)**3)/3 + 4*mOpt**2 + (38*mOpt**3)/3 + 2*n**2 + &
            mOpt*n*(2*mOpt + 2*n + 3)) + const
    else
        cost2 = dble(4*mOpt + 8*n + 4*p + 2*s*(mOpt + p)**3 + 27*(mOpt + p)**2 + &
            (44*(mOpt + p)**3)/3 + 2*mOpt**2 + (2*mOpt**3)/3 + 2*n**2 + &
            mOpt*n*(2*mOpt + 2*n + 3)) + const
    end if
    !cost2 = cost2 * dble(ceiling((h-tnow)/tstep))
    cost2 = cost2 * ((h-tnow)/tstep)
    if (cost2 .lt. 0.0_WP) then
        write(option%ioutfile,*) cost2, h-tnow, tstep
        stop 'cost2 < 0'
    end if
    
    ! change both *** testing feature, the efficiency is to be evaluated
    if (option%SchurDecomp) then
        cost3 = dble(5*mOptmid + 8*n + 4*p + 2*s*(mOptmid + p)**3 + 27*(mOptmid + p)**2 + &
            (44*(mOptmid + p)**3)/3 + 4*mOptmid**2 + (38*mOptmid**3)/3 + 2*n**2 + &
            mOptmid*n*(2*mOptmid + 2*n + 3)) + const
    else
        cost3 = dble(4*mOptmid + 8*n + 4*p + 2*s*(mOptmid + p)**3 + 27*(mOptmid + p)**2 + &
            (44*(mOptmid + p)**3)/3 + 2*mOptmid**2 + (2*mOptmid**3)/3 + 2*n**2 + &
            mOptmid*n*(2*mOptmid + 2*n + 3)) + const
    end if
    !cost3 = cost3 * dble(ceiling((h-tnow)/timeOptmid))  
    cost3 = cost3 * ((h-tnow)/timeOptmid)
    if (cost3 .lt. 0.0_WP) then
        write(option%ioutfile,*) cost3, h-tnow, tstep
        stop 'cost3 < 0'
    end if    
    
    if (fixmFlag) then
        imode = 1
        factort = timeOpt/tstep
        factorm = mfixed-m
        if(option%InfoOutput) write(option%ioutfile, fmt_mfix) 
    elseif (cost3 .lt. cost2 .and. cost3 .lt. cost1) then
        ! change both
        imode = 3
        factort = timeOptmid/tstep
        factorm = mOptmid-m
        if(option%InfoOutput) write(option%ioutfile, fmt_mode3) cost1, cost2, cost3
    elseif (mOpt .ge. m .and. m .eq. mUplimit) then
        ! change time
        imode = 1
        factort = timeOpt/tstep
        factorm = 0
        if(option%InfoOutput) write(option%ioutfile, fmt_modem)   
    elseif (cost1 .le. cost2) then
        ! change time
        imode = 1
        factort = timeOpt/tstep
        factorm = 0
        if(option%InfoOutput) write(option%ioutfile, fmt_mode1) cost1, cost2, cost3
    else
        ! change m
        imode = 2
        factort = 1.0_WP
        factorm = mOpt-m
        if(option%InfoOutput) write(option%ioutfile, fmt_mode2) cost1, cost2, cost3
    end if    
    
    end subroutine AdaptNextStep
    !#########################################
    
    !#########################################
    subroutine KrylovSolve(n, h, y, lwork, rwrk, option, matvec, updateyw, obtainJac, status, info)
    implicit none
    integer :: n, info
    integer,dimension(3) :: status
    real(WP) :: h
    real(WP),dimension(n) :: y
    integer :: lwork
    real(WP),dimension(lwork) :: rwrk
    type(expkry_opts) :: option        
    external :: matvec, updateyw, obtainJac 
    
    ! general control
    integer :: iPosKry, lworkRequired, lphiwrk
    integer :: i
    logical :: ArnoldiFlag
    
    ! Arnoldi
    integer ::  mStart, m, mEnd, ldm, m_old, p, mUplimit, mLowlimit, mEig
    real(WP),dimension(:),  pointer :: w, phiw, phiwrk, vec
    real(WP),dimension(:,:),pointer :: phiwm
    real(WP),dimension(:,:),pointer :: V, Hmax, Hm
    real(WP) :: beta, hm1m
    
    ! schur
    real(WP),dimension(:,:),pointer :: UHmax, UHm, THmax, THm
    real(WP),dimension(:),pointer :: eigrH, eigiH, TauH, rwrkSch, workrcond
    integer :: lworkSch, lworkrcond
    real(WP) :: eigabsmin, spectradius, eigmax, condNum, norm  
    
    ! Sub timestep
    real(WP) :: tnow, tlimit, tstep, tstep_old, toly(n), tolymin
    integer :: ierr, iscale, istep, icheckEigmax, mfixed, istepExpectRemain
    logical :: fixmFlag, stopCheckEigmin, stopCheckcond, stopCheckEigmax
    
    ! Error control
    real(WP) :: betahphi2, errloc, omega, omega_old, factort
    integer :: imode, ntry, factorm

    ! info output
    character(len=*),parameter ::   &

        fmt_strt   = "('   --- Start of an attemption, No.: ', I8)",    & ! nreject
        fmt_Arnds  = "('       Arnoldi from ',I8, ' to ', I8)",         & ! mStart mEnd
        fmt_brk    = "('   *** Krylov subspace breakdown.',A1)",        &
        fmt_Arndf  = "('       Arnoldi finished, dimension: ', I8)",    & ! m
        fmt_schur  = "('       Schur decompostion finished.')",         & ! condNum, eigmax
        fmt_eigmin = "('       Zero eigenvalue found. m decrease: ', I8)", &
        fmt_fixm   = "('       Use fixed m: ', I8)",                    &
        fmt_norm   = "('       1-Norm Number         : ', ES15.6)",     &
        fmt_eigmax = "('       Maxmimum eigenvalue   : ', ES15.6)",     &
        fmt_tlimit = "('       Time step limit       : ', ES15.6 )",    & ! tlimit
        fmt_stepExp= "('       Steps expected        : ', I8 )",        & ! istepExpectRemain
        fmt_stopeig= "('        Stop checking max eigenvalues, ', I8)",  &
        fmt_m1     = "('       1-Norm is too large. change m and try agin: ', I8)", & ! m     
        fmt_m2     = "('       Eigmax is too large. change m and try agin: ', I8)", & ! m  
        fmt_tstep  = "('       Time step in this try : ', ES15.6)",      & ! tstep
        fmt_m      = "('       m         in this try : ', I8)",          & ! m
        fmt_errs   = "('       Errors: ', 3ES15.6)",                    & ! errloc omega
        fmt_omega  = "('       Local error and Omega: ', 2ES15.6)",     & ! errloc omega
        fmt_rej    = "('   *** Time step rejected.',A1)",               &
        fmt_mode1  = "('       Change time step to: ', ES15.6)",        & ! tstep
        fmt_mode2  = "('       Change m    step to: ', I8)",            &     ! m 
        fmt_subn   = "('   Time advanced. Substep ',I8)",               & ! istep
        fmt_subnow = "('   Time now      : ', ES15.6)",                 & ! tnow
        fmt_substp = "('   Sub time step : ', ES15.6 )" ,               & ! tstep
        fmt_subm   = "('   m of this step : ', I8)",                    &   
        fmt_stperr = "('   Error estimate this step:', ES15.6)",        &  !   
        fmt_stpend = "('   -------------------------------------------')"
        
    real(WP),external :: DDOT
    
    ! check input
    p = 2
    ldm = min(option%MaxKrylovSize, n-1) + p
    lphiwrk = 7*ldm*ldm
    lworkSch = ldm*ldm + 75*ldm + 4160
    
    if (option%SchurDecomp) then
        ! 2*n + n*(ldm+1) + ldm*ldm + ldm + ldm*p + 7*ldm*ldm + 2*ldm*ldm + 3*ldm + ldm*ldm + 75*ldm + 4160 + (ldm+4)*ldm
        lworkRequired = 83*ldm + 3*n + ldm*n + ldm*p + 12*ldm*ldm + 4160
    else
        ! 2*n + p*ldm + n*(ldm+1) + ldm*ldm + 7*ldm*ldm + (ldm+4)*ldm
        lworkRequired = 4*ldm + 3*n + ldm*n + ldm*p + 9*ldm*ldm
    end if
    if (lwork .lt. lworkRequired) then
        write(option%ioutfile,"('KrylovSolve : Not enough workspace (Provided/Required)', I8, I8)")lwork, lworkRequired
        stop
    end if
    
    info = 0
    status = 0
    status(1) = 2
    
    ! pointers
    iPosKry = 1
    w    => assignWorkspace(iPosKry, lwork, rwrk, n)         ! lwork: + n
    phiw => assignWorkspace(iPosKry, lwork, rwrk, n)         ! lwork: + n
    V    => assignWorkspace(iPosKry, lwork, rwrk, n, ldm+1)    ! lwork: + n*(ldm+1)
    Hmax => assignWorkspace(iPosKry, lwork, rwrk, ldm, ldm)  ! ldm*ldm
    
    if (option%SchurDecomp) then
        UHmax => assignWorkspace(iPosKry, lwork, rwrk, ldm, ldm)   ! lwork: + ldm*ldm
        THmax => assignWorkspace(iPosKry, lwork, rwrk, ldm, ldm)   ! lwork: + ldm*ldm
        eigrH  => assignWorkspace(iPosKry, lwork, rwrk, ldm)    ! lwork: + ldm
        eigiH  => assignWorkspace(iPosKry, lwork, rwrk, ldm)    ! lwork: + ldm
        TauH   => assignWorkspace(iPosKry, lwork, rwrk, ldm)    ! lwork: + ldm
    end if
    
    ! initialization
    tnow  = 0.0_WP
    tstep  = h
    tlimit = h
    phiw = 0.0_WP
    V = 0.0_WP
    Hmax = 0.0_WP
    mUplimit = min(option%MaxKrylovSize, n-1)
    mLowlimit = min(option%MinKrylovSize, n-1)
    m = mLowlimit
    fixmFlag = .false.
    mfixed = 0
    mStart = 1
    mEnd   = 0
    omega = 0.0_WP
    omega_old = 0.0_WP
    ArnoldiFlag = .true.
    imode = 0
    istep = 1
    ntry = 1
    icheckEigmax = 0
    stopCheckEigmin = .false.
    stopCheckcond = .false.
    stopCheckEigmax = .false.
    m_old = 0
    tstep_old = 0.0_WP
    istepExpectRemain = 1

    call updateyw(n, phiw, y, w, ierr) ! get the initial w 
    if (ierr.ne. 0) then
        info = 2
        if(option%InfoOutput) write(option%ioutfile,"('KrylovSolve : Error in getting the initial w.')")
        return
    end if
    !
    ! main loop:
    do while (tnow+eps .lt. h)
        
        if (option%InfoOutput) write(option%ioutfile,fmt_strt) ntry
        ! if too many rejects
        if (ntry .gt. option%MaxReject) then
            info = 2
            if(option%InfoOutput) write(option%ioutfile,"('Too manny rejects in this step: ', I8 )") ntry
            return
        end if
        
        if (ArnoldiFlag) then
            
            if(option%InfoOutput) write(option%ioutfile,fmt_Arnds) mStart, m
            if (mEnd .lt. m) then
                !if (mStart .ne. 1) then
                !    call print_dense_to_file(w, n, 1, 'wvec')
                !end if
                if (option%reorth .or. m .gt. 50) then
                    mEnd = m
                    mStart = 1
                    call dcpArnoldi(n, w, mStart, mEnd, matvec, option, V, Hmax, ldm, beta, hm1m)
                    !if (mStart .ne. 1) then
                    !    call print_dense_to_file(Hmax, ldm, ldm, 'Hmaxdcp')
                    !    call print_dense_to_file(V, n, ldm+1, 'Vdcp')
                    !end if
                else
                    ! flops: 3*n + n*m*(2*m + 2*n + 3)
                    mEnd = m
                    call Arnoldi   (n, w, mStart, mEnd, matvec, option, V, Hmax, ldm, beta, hm1m)
                    !call print_dense_to_file(Hmax,  ldm, m, 'Hmax')
                end if
                
                if (m .ne. mEnd) then
                    if(option%InfoOutput) write(option%ioutfile,fmt_brk) ' '
                    m = mEnd  
                    if (m .lt. mLowlimit) then
                        mLowlimit = m
                    end if
                end if
                ! mEnd is the maximum m that is computed.
            else
                ! do not need to compute, retrieve from storage
            end if

            if(option%InfoOutput) write(option%ioutfile,fmt_Arndf) m
            if (m .eq. 0) then
                if(option%InfoOutput) write(option%ioutfile,"('m = 0, return without doing anything.')")
                return
            end if
            
            ! Schur
            ! flops: 10*m^3
            if (option%SchurDecomp) then 
                lworkSch = m*m + 75*m + 4160
                rwrkSch  => assignWorkspace(iPosKry, lwork, rwrk, lworkSch)  ! lwork: + lworkSch  
                call realSchur (m, Hmax, ldm, UHmax, ldm, THmax, ldm, eigrH, eigiH, TauH, lworkSch, rwrkSch, .true., info) 
                rwrkSch  => nullWorkspace(iPosKry, lworkSch, 'rwrkSch')
                if (info .ne. 0) then
                    info = 2
                    if(option%InfoOutput) then
                        write(option%ioutfile,"('KrylovSolve: Error in realSchur', 3I8)") info, mStart, mEnd
                        call print_dense_to_file(Hmax,  ldm, m,   'ErrorCheck_Hmax')
                        call print_dense_to_file(V,     n,   m+1, 'ErrorCheck_V')
                        call print_dense_to_file(THmax, ldm, m,   'ErrorCheck_THmax')
                        call print_dense_to_file(UHmax, ldm, m,   'ErrorCheck_UHmax')
                        call print_dense_to_file(w, n, 1, 'ErrorCheck_w')
                    end if
                    return
                end if
                !condNum = eigabsmax/max(eigabsmin,eps)
                if(option%InfoOutput) write(option%ioutfile,fmt_schur)
                
                eigmax = maxval(eigrH(1:m)) 
                spectradius =  maxval(abs(eigrH(1:m))) 
                if (option%ComputSchurQuotient) then
                    norm = Mat1norm(m, m, THmax, ldm)
                else
                    norm = Mat1norm(m, m, Hmax, ldm)
                end if
                
                ! 1.  check zero eigenvalue : flops: 8*m
                if (.not.stopCheckEigmin) then
                    mEig = m
                    do i = 1, m
                        if (sqrt(eigrH(i)*eigrH(i)+eigiH(i)*eigiH(i)) .lt. norm*option%uround) then
                            mEig = mEig-1
                        end if
                    end do
                    if (mEig .lt. m) then
                        m = mEig
                        mUplimit = mEig
                        if (option%InfoOutput) write(option%ioutfile,fmt_eigmin) mEig
                        if (mEig .lt. mLowlimit) then
                            mLowlimit = mEig
                            fixmFlag = .true.
                            mfixed = mEig
                            stopCheckEigmin = .true.
                            if (option%InfoOutput) write(option%ioutfile,fmt_fixm) mfixed
                        end if
                        mStart = mEnd
                        ArnoldiFlag = .true.
                        ntry = ntry + 1 
                        cycle
                    end if
                end if
                
                toly = tolyra (n, y, m, option)
                tolymin = minval(toly)  
                if(option%InfoOutput) write(option%ioutfile,"('       Min tol: ', ES15.6)") tolymin
                
                ! 2. check maximum positive real eigenvalue  
                tlimit = min(tlimit, h-tnow)
                tlimit = setTimestepLimitEigmax (m, eigmax, tlimit, tolymin,  option%uround)
                istepExpectRemain = ceiling((h-tnow)/(tlimit+option%uround))
                if (istepExpectRemain .lt. 0) then
                    write(*,*) "Consider to increase atol."
                    stop 
                end if
                    
                if (option%InfoOutput) then
                    write(option%ioutfile,fmt_eigmax) eigmax 
                    write(option%ioutfile,fmt_tlimit) tlimit
                    write(option%ioutfile,fmt_stepExp) istepExpectRemain 
                end if
                if (.not.stopCheckEigmax) then
                    if (istepExpectRemain .gt. option%Maxsubsteps-istep) then
                        ! too many steps
                        if (icheckEigmax .ge. option%MaxcheckEigmax) then    
                            ! eigmax may come from the original matrix
                            ! stop checking eigenvalue
                            if (option%InfoOutput) write(option%ioutfile,fmt_stopeig) icheckEigmax
                            stopCheckEigmax = .true.
                        else
                            icheckEigmax = icheckEigmax + 1
                            if (m .eq.mUplimit) then
                                mUplimit = 0.8 * m
                                m = mUplimit
                                if (m .lt. mLowlimit) then
                                    mLowlimit = m
                                    fixmFlag = .true.
                                    mfixed = m
                                    stopCheckEigmax = .true.
                                    if (option%InfoOutput) write(option%ioutfile,fmt_fixm) mfixed
                                end if
                            else
                                m = m + 1
                            end if
                            if (option%InfoOutput) write(option%ioutfile,fmt_m2) m
                            mStart = mEnd
                            ArnoldiFlag = .true.
                            ntry = ntry + 1 
                            cycle
                        end if
                    end if
                end if
                
                
                ! 3. check Matrix 1-norm
                tlimit = setTimestepLimitNorm (m, norm, tlimit, tolymin,  option%uround)
                istepExpectRemain = ceiling((h-tnow)/(tlimit+option%uround))
                if (option%InfoOutput) then
                    write(option%ioutfile,fmt_norm)   norm 
                    write(option%ioutfile,fmt_tlimit) tlimit 
                    write(option%ioutfile,fmt_stepExp) istepExpectRemain 
                end if
                if (.not. stopCheckcond) then
                    if ( istepExpectRemain .gt. option%Maxsubsteps-istep ) then
                        ! too many steps
                        mUplimit = 0.8 * m
                        m = mUplimit
                        if (option%InfoOutput) write(option%ioutfile,fmt_m1) m
                        if (m .lt. mLowlimit) then
                            mLowlimit = m
                            fixmFlag = .true.
                            mfixed = m
                            stopCheckcond = .true.
                            if (option%InfoOutput) write(option%ioutfile,fmt_fixm) mfixed
                        end if
                        mStart = mEnd
                        ArnoldiFlag = .true.
                        ntry = ntry + 1 
                        cycle
                    end if
                end if                
                
            else
                !uround = dble(m**3)*eps
                norm   = option%uround!eps
                eigmax = option%uround!eps
            end if
        end if
        ! m for the current try is set.
    
        if (tlimit .lt. 10.0_WP*eps) then
            info = 2
            write(option%ioutfile,"('Cannot reach the required tolerance:',ES15.6)") tlimit
            return
        end if       
        if (m .eq. m_old .and. imode .eq. 2) then
            imode = 0
        end if
            
        ! limit time step
        !tstep = min(tstep, tlimit) 
        tstep = min(tstep, (h-tnow)/dble(istepExpectRemain))

        ! update step expected
        istepExpectRemain = ceiling((h-tnow)/(tstep+option%uround))
            
        if(option%InfoOutput) then
            write(option%ioutfile,fmt_tstep) tstep
            write(option%ioutfile,fmt_m)     m             
            write(option%ioutfile,"(3I8, ES15.6)") istep, ntry, m, tstep
            write(option%ioutfile,fmt_stepExp) istepExpectRemain 
        end if
        
        vec  => assignWorkspace(iPosKry, lwork, rwrk, m)    
        phiwm=> assignWorkspace(iPosKry, lwork, rwrk, m, p)          
        lphiwrk = 7*(m+p)*(m+p)
        phiwrk => assignWorkspace(iPosKry, lwork, rwrk, lphiwrk)        
        if (option%ComputSchurQuotient) then
            ! phiwm = UH*phi(UH)*UH^T*w
            ! vec = UH^T * vec, vec = e1
            call DCOPY(m,UHmax,ldm,vec,1) !  vec(1:m) = UH(1,1:m)
            call solvePhiw(m, tstep, THmax, ldm, vec, p, phiwm, iscale, lphiwrk, phiwrk, option)
            ! phiwm = UH * phiwm
            vec = phiwm(1:m,1)
            call DGEMV ('N',m,m,1d0,UHmax,ldm,vec,1,0d0,phiwm,1) ! flops: 2*m^2 
            !phiwm(1:m,1) = MV_blas(m, UHm, phiwm(1:m,1)) ! flops: 2*m^2 
            phiwm(m,2) = DDOT(m,UHmax(m,1),ldm,phiwm(1,2),1) ! UH(m,1:m)*phiwm(1:m,2) ! flops: m
        else
            vec = 0.0_WP; vec(1) = 1
            call solvePhiw(m, tstep, Hmax,  ldm, vec, p, phiwm, iscale, lphiwrk, phiwrk, option)
        end if
        phiwrk => nullWorkspace(iPosKry, lphiwrk, 'phiwrk')  
       
        !call print_dense_to_file(UHmax, ldm, m, 'UHmax')
        ! local error for this time step
        betahphi2 = beta * hm1m * phiwm(m,2)
        errloc = max(abs(betahphi2), option%uround*norm*tstep, option%uround*exp(eigmax*tstep))
        
        if(option%InfoOutput) write(option%ioutfile, fmt_errs) abs(betahphi2), option%uround*norm*tstep, option%uround*exp(eigmax*tstep)

        toly = tolyra (n, y, mUplimit, option)
        tolymin = minval(toly) 
        ! omega = error/tolerance
        ! omega(dt) = errloc(dt)/(tolerance*dt/h)
        omega_old = omega
        omega = 0.0_WP
        ! flops : 5*n
        do i = 1, n
            omega = max(omega, abs(V(i,m+1)/toly(i)))
        end do        
        omega = errloc*omega !*(h/tstep)
        
        if(option%InfoOutput) write(option%ioutfile, fmt_omega) errloc, omega
    
        call AdaptNextStep(omega, omega_old, tstep, tstep_old, h, tnow, tlimit,  &
            m, m_old, mLowlimit, mUplimit, n, iscale, fixmFlag, mfixed, option, imode, factort, factorm)
        
        ! error check
        if (omega .gt. option%delta) then
            ntry = ntry + 1 
            if(option%InfoOutput) write(option%ioutfile, fmt_rej) ' '

            phiwm=> nullWorkspace(iPosKry, m, p, 'phiwm')
            vec  => nullWorkspace(iPosKry, m, 'vec') 
                        
            if (imode .eq. 1) then
                ! change time
                ArnoldiFlag = .false.
                tstep_old = tstep
                tstep = tstep * factort
                if(option%InfoOutput) write(option%ioutfile, fmt_mode1) tstep
            else
                ! change m, imode = 2 or 3
                if (imode .eq. 3) then
                   ! change both
                    tstep_old = tstep
                    tstep = tstep * factort
                    if(option%InfoOutput) write(option%ioutfile, fmt_mode1) tstep     
                end if
                ArnoldiFlag = .true.
                mStart = mEnd
                m_old = m
                m = m + factorm
                if(option%InfoOutput) write(option%ioutfile, fmt_mode2) m
            end if
            cycle ! try again
        end if
        
        ! update results 

        ! compute w_int(w) = beta*Vm*[t_step*phi1(t_step*Hbar)*e1]
        phiw = V(1:n,m+1)
        call DGEMV('n',n,m,beta,V,n,phiwm(1:m,1),1, betahphi2,phiw,1) ! flops: 2*n^2
        
        call updateyw(n, phiw, y, w, ierr)
        
        if (ierr .ne. 0) then
            ! y is not updated
            if (.not. option%ErrGiveUp) then
                !tlimit = 0.1_WP*tstep
                istepExpectRemain = istepExpectRemain * 5
                tlimit = (h-tnow)/dble(istepExpectRemain) 
                if(istepExpectRemain .lt. option%Maxsubsteps-istep) then
                    if(option%InfoOutput) then
                        write(option%ioutfile,"('KrylovSolve : Error in updating y and w.')") ! or change and cycle? 
                        write(option%ioutfile,"('Further reduce tstep to ', ES15.6)") tlimit
                        write(option%ioutfile,"('Further increase   m to ', I8)") min(m*2, mUplimit)
                        write(option%ioutfile,fmt_stepExp) istepExpectRemain 
                    end if
                    phiwm=> nullWorkspace(iPosKry, m, p, 'phiwm')
                    vec  => nullWorkspace(iPosKry, m, 'vec') 
                    tstep = tlimit
                    m = min(m*2, mUplimit)
                    mStart = mEnd
                    ArnoldiFlag = .true.
                    imode = 0
                    ntry = ntry + 1 
                    cycle
                else
                    if(option%InfoOutput) then
                        write(option%ioutfile,"('Too many steps will be resulted', ES15.6)") istepExpectRemain
                    end if
                end if
            end if
            
            ! give up
            if(option%InfoOutput) then
                write(option%ioutfile,"('KrylovSolve : Error in updating y and w. Giving up ...')") 
            end if
            info = 3
            phiwm=> nullWorkspace(iPosKry, m, p, 'phiwm')
            vec  => nullWorkspace(iPosKry, m, 'vec')    
            if (option%SchurDecomp) then
                ! schur decomposition
                TauH   => nullWorkspace(iPosKry, ldm, 'TauH')  
                eigiH  => nullWorkspace(iPosKry, ldm, 'eigiH')   
                eigrH  => nullWorkspace(iPosKry, ldm, 'eigrH')  
                THmax  => nullWorkspace(iPosKry, ldm, ldm, 'THmax')   
                UHmax  => nullWorkspace(iPosKry, ldm, ldm, 'UHmax')  
            end if
            Hmax => nullWorkspace(iPosKry, ldm, ldm, 'Hmax') 
            V    => nullWorkspace(iPosKry, n, ldm+1, 'V') 
            phiw => nullWorkspace(iPosKry, n, 'phiw')         ! lwork: - n
            w    => nullWorkspace(iPosKry, n, 'w')         ! lwork: - n                
            return                
        end if
        
        ! update this step 
        tnow = tnow + tstep
        
        if(option%InfoOutput)then
            write(option%ioutfile,fmt_subn)   istep
            write(option%ioutfile,fmt_subnow) tnow
            write(option%ioutfile,fmt_substp) tstep
            write(option%ioutfile,fmt_subm)   m
            write(option%ioutfile,fmt_stperr) errloc
            write(option%ioutfile,fmt_stpend) 
        end if
    
        status(2) = max(status(2), m)
        status(3) = status(3) + 1
        
        istep = istep + 1
        
        if ( mod(istep, option%JACinterv) .eq. 0) then
            ! update Jacobian
            call obtainJac (2, n, y)
        end if
        
        ! reset error controlling parameters
        tlimit = h-tnow
        istepExpectRemain = 1
        ntry = 1
        tstep_old = 0.0_WP
        m_old = 0
        if (.not.fixmFlag) then
            mUplimit = min(option%MaxKrylovSize, n-1)
            mLowlimit = min(option%MinKrylovSize, n-1)
        end if
        omega = 0.0_WP
        omega_old = 0.0_WP
        mStart = 1
        mEnd   = 0
        ArnoldiFlag = .true.
        V = 0.0_WP
        Hmax = 0.0_WP
        icheckEigmax = 0
        !fixmFlag
        !mfixed
        !stopCheckEigmin = .false.
        !stopCheckcond = .false.
        !stopCheckEigmax = .false.
        phiwm=> nullWorkspace(iPosKry, m, p, 'phiwm')
        vec  => nullWorkspace(iPosKry, m, 'vec') 


        ! next step
        if (istep .ge. option%Maxsubsteps)then
            info = 4
            return
        end if
        
        tstep = min(tstep * factort, h-tnow) ! tnow has changed
        m = m + factorm
        imode = 0

    end do

    ! finalize 
    if (option%SchurDecomp) then
        ! schur decomposition
        TauH   => nullWorkspace(iPosKry, ldm, 'TauH')  
        eigiH  => nullWorkspace(iPosKry, ldm, 'eigiH')   
        eigrH  => nullWorkspace(iPosKry, ldm, 'eigrH')  
        THmax  => nullWorkspace(iPosKry, ldm, ldm, 'THmax')   
        UHmax  => nullWorkspace(iPosKry, ldm, ldm, 'UHmax')  
    end if
    Hmax => nullWorkspace(iPosKry, ldm, ldm, 'Hmax') 
    V    => nullWorkspace(iPosKry, n, ldm+1, 'V') 
    phiw => nullWorkspace(iPosKry, n, 'phiw')         ! lwork: - n
    w    => nullWorkspace(iPosKry, n, 'w')         ! lwork: - n
    

    end subroutine KrylovSolve
    !#########################################

end module
