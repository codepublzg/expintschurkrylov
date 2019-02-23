module expmA
  ! ************************************************************************************************
  ! This is a part of the source code of the Exponential Integrator with Schur-Krylov Approximation.
  ! Objective:  To calculate the exponential of matrix, expm(A).
  ! Method:     The Padé approximation with scaling and squaring.
  ! Reference:  N.J. Higham, The scaling and squaring method for the matrix exponential revisited. 
  !             SIAM J. Matrix Anal. A., 26 (2005), pp. 1179-1193. https://doi.org/10.1137/04061101X
  ! Please follow updates or send comments through the web page 
  !     <https://github.com/codepublzg/expintschurkrylov.git.>
  ! Developed and Published online @2019.02.23 by Zaigang Liu.
  ! ************************************************************************************************
  ! usage: 
  !
  !   Step 1. Invoke this module:
  !
  !     use expmA
  !
  !   Step 2. Create a expmA_opts data type, e.g.:
  !
  !     type(expmA_opts) :: expmOption
  !
  !   Step 3. Set options at a proper place in your code, e.g.:
  !
  !     expmOption = expmSetOptions(performSchur  = .true., Hessenberg = .false.)
  !
  !   Step 4. Voilà:
  !
  !     call dexpm_pade (N, A, F, iscale, nrLength, rwrksp, expmOption, iflag)
  !
  ! ************************************************************************************************
  implicit none
  private

  ! double precision 
  integer, parameter :: WP = kind(1d0)
    
  type expmA_opts
    logical :: performSchur
    logical :: Hessenberg
    logical :: SchurQuotient
    logical :: opt_init = .false.
  end type expmA_opts
  
  interface assignWorkspace
    module procedure assignWorkspace1D
    module procedure assignWorkspace2D
  end interface assignWorkspace
  
  interface nullWorkspace
    module procedure nullWorkspace1D
    module procedure nullWorkspace2D
  end interface nullWorkspace

  public :: expmA_opts, expmSetOptions, realSchur, dexpm_pade
  public :: assignWorkspace, nullWorkspace, Matmul_blas, MatmulT_blas
  public :: Mat1norm

  contains

  !#########################################
  ! performSchur -- Logical, TRUE:Perform Schur decomposiiton before Padé approximation.
  ! Hessenberg   -- Logical, TRUE:The input matrix is an upper Hessenberg matrix.
  ! SchurQuotient-- Logical, TRUE:The input matrix is in Schur form.
  function expmSetOptions  (  &
    performSchur,   &
    Hessenberg,     &
    SchurQuotient   &
    ) result(option)
    implicit none
    
    type(expmA_opts) :: option
    logical,  optional, intent(in) :: &
        performSchur, Hessenberg, SchurQuotient

    if (present(performSchur)) then
        option%performSchur = performSchur
    else
        option%performSchur = .true.
    end if
    if (present(Hessenberg)) then
        option%Hessenberg = Hessenberg
    else
        option%Hessenberg = .false.
    end if
    if (present(SchurQuotient)) then
        option%SchurQuotient = SchurQuotient
    else
        option%SchurQuotient = .true.
    end if
    
    if (option%SchurQuotient) option%performSchur  = .false.

    option%opt_init = .true.
    
  end function expmSetOptions
  !#########################################

  !#########################################
  function assignWorkspace1D(iPosNow, lwrk, wrk, nele) result(p)
    implicit none
    integer, parameter :: WP = kind(1d0)
    integer :: iPosNow, lwrk, nele
    real(WP), dimension(lwrk), target  :: wrk
    real(WP), dimension(:), pointer :: p

    p(1:nele) => wrk(iPosNow:iPosNow+nele-1)
    iPosNow = iPosNow + nele
    if (iPosNow-1 .gt. lwrk) then
        write (*,"('assignWorkspace1D : No more space, iPosNow ',I8, ' lwrk', I8)") iPosNow, lwrk
        stop 
    end if
  end function assignWorkspace1D
  !#########################################

  !#########################################
  function assignWorkspace2D(iPosNow, lwrk, wrk, mele, nele) result(p)
    implicit none
    integer, parameter :: WP = kind(1d0)
    integer :: iPosNow, lwrk, mele, nele
    real(WP), dimension(lwrk), target :: wrk
    real(WP), dimension(:,:), pointer :: p

    p(1:mele,1:nele) => wrk(iPosNow:iPosNow+mele*nele-1)
    iPosNow = iPosNow + mele*nele
    if (iPosNow-1 .gt. lwrk) then
        write (*,"('assignWorkspace2D : No more space, iPosNow ',I8, ' lwrk', I8)") iPosNow, lwrk
        stop 
    end if
    
  end function assignWorkspace2D
  !#########################################

  !#########################################
  function nullWorkspace1D(iPosNow, nele, label) result(p)
    implicit none
    integer, parameter :: WP = kind(1d0)
    integer :: iPosNow, lwrk, nele
    character(len=*) :: label
    real(WP), dimension(:), pointer :: p

    nullify(p)
    iPosNow = iPosNow - nele
    if (iPosNow .lt. 1) then
        write(*,"('nullWorkspace1D : Multi nullified, iPosNow ',I8)") iPosNow
        write(*,"('Pointer label: ', A)") label
        stop
    end if

  end function nullWorkspace1D
  !#########################################

  !#########################################
  function nullWorkspace2D(iPosNow, mele,nele, label) result(p)
    implicit none
    integer, parameter :: WP = kind(1d0)
    integer :: iPosNow, lwrk, mele, nele
    character(len=*) :: label
    real(WP), dimension(:,:), pointer :: p

    nullify(p)
    iPosNow = iPosNow - mele*nele
    if (iPosNow .lt. 1) then
        write(*,"('nullWorkspace2D : Multi nullified, iPosNow ',I8)") iPosNow
        write(*,"('Pointer label: ', A)") label
        stop
    end if

  end function nullWorkspace2D
  !#########################################

  !#########################################
  subroutine log2fr (x, f, e)
    implicit none
  
    real(WP) :: x, f
    integer :: e
    
    e = ceiling(log(x)/log(2d0))
    f = x/dble(2**e)

  end subroutine log2fr
  !#########################################

  !#########################################
  ! C = A*B
  function Matmul_blas (N, A, B) result(C)
    implicit none
    integer :: N
    real(WP),dimension(N,N) :: A, B, C
  
    call DGEMM('N','N',N, N, N, 1d0, A, N, B, N, 0d0, C, N)

  end function Matmul_blas
  !#########################################

  !#########################################
  ! C = A*B^T
  function MatmulT_blas (N, A, B) result(C)
    implicit none
    integer :: N
    real(WP),dimension(N,N) :: A, B, C
  
    call DGEMM('N','T',N, N, N, 1d0, A, N, B, N, 0d0, C, N)

  end function MatmulT_blas
  !#########################################
  
    !#########################################
    function MatInfNorm(m, n, A)
      implicit none
  
      real(WP) :: sumn, MatInfNorm, A(m,n)
      integer :: m, n, i, j
  
      n = size(A,2)
      MatInfNorm = 0

      do i = 1, m
        sumn = 0d0
        do j = 1, n
          sumn = sumn + abs(A(i,j))
        end do
        if (isnan(sumn)) then
          write(*,*) 'Nan in MatInfNorm'
          stop
        else if (sumn .gt. MatInfNorm) then
          MatInfNorm = sumn
        end if
      end do      

    end function MatInfNorm
    !#########################################

    !#########################################
    function Mat1norm(m, n, A, lda)
      implicit none
  
      real(WP) :: summ, Mat1norm, A(lda,n)
      integer :: lda, m, n, i, j
  
      Mat1norm = 0
  
      do j = 1, n
        summ = 0d0
        do i = 1, m
          summ = summ + abs(A(i,j))
        end do
        if (isnan(summ)) then
          write(*,*), 'error in Mat1norm'
          stop
        else if (summ .gt. Mat1norm) then
          Mat1norm = summ
        end if
      end do      

      return
    end function Mat1norm
    !#########################################
     
  !#########################################
  subroutine realSchur (N, A, lda, U, ldu, T, ldt, eigr, eigi, Tau, lworkin, work, lhess, info)
    ! T = U^T*A*T
    ! lworkin at least:
    !   not Hessenberg : (65+N)*64+N*N+11*N = N*N + 75*N + 4160
    !   Hessenberg     : 11*N if lhess = .false.
    ! flops: 
    !   not Hessenberg : (10/3)*n^3 + 10*n^3 + 2*n^3
    !   Hessenberg     : 10*n^3
    !
    implicit none
    integer :: N, lworkin, lwork, lR, lworkH, lworkS
    integer :: lda, ldu, ldt
    real(WP),dimension(lda,N) :: A
    real(WP),dimension(ldu,N) :: U
    real(WP),dimension(ldt,N) :: T
    real(WP),dimension(N) :: eigr, eigi, Tau
    real(WP),dimension(lworkin),target :: work
    logical :: lhess

    real(WP),dimension(:,:),pointer :: R

    integer :: info
    
    info = 0
    
    !---check workspace length
    lworkS = 11*N
    if (.not.lhess) then
      lR = N*N
      lworkH = (65+N)*64
      lwork = lR+lworkH+lworkS
    else
      lR = 0
      lworkH = 0
      lwork = lworkS
    end if
    if (lworkin .lt. lwork) stop 'Not enough workspace for Schur : realSchur'

    !---calling Lapack subroutine to perform Schur decomp.
    ! T = A
    call DLACPY( 'A', N, N, A, lda, T, ldt)

    if (.not.lhess) then
      ! H = Q1^T*A*Q1     
      R(1:N,1:N) => work(1:N*N) ! store the reflectors
      call DGEHRD (N, 1, N, T, ldt, Tau, work(lR+1), lworkH, info) ! flops : (10/3)*n^3
      if (info .ne. 0) then
        info = 1
        nullify(R)
        return
      end if
      R = T
    end if

    ! T =Q2^T*H*Q2 
    call DHSEQR ('S', 'I', N, 1, N, T, ldt, eigr, eigi, U, ldu,  &
        work(lR+lworkH+1), lworkS, info)                        ! flops: 10*n^3
    if (info .ne. 0) then
        info = 2
        nullify(R)
        return
    end if

    if (.not.lhess) then
      ! U = Q1*Q2
      CALL DORMHR('L', 'N',N ,N, 1, N, R, N, Tau, U, ldu,     &
        work(lR+1), lworkH, info)                               ! flops: 2*n^3
      if (info .ne. 0) then
          info = 3
          nullify(R)
          return
      end if
      nullify(R)
    end if

  end subroutine realSchur
  !#########################################

  !#########################################
  subroutine dexpm_pade (N, A, F, iscale, nrLength, rwrksp, option, iflag)
    implicit none
    ! comput matrix exponential using matlab expm method. lzg
    ! 
    ! N -- (in)dimension of the square matrix
    ! A -- input matrix, double
    ! F -- output matrix, double
    ! iscale   -- scale of matrix A, integer, see squaring and scaling method
    ! nrLength -- integer, length of the rwrksp: 7*N*N + 78*N + 4160 (Schur), 4*N*N (no Schur)
    ! rwrksp   -- double workspace at length of at least nrLength
    ! iflag    -- error flag, integer
    ! 
    ! flops: 
    !   s = ceiling(log(Norm/5.37)/log(2))
    !   SchurQuotient :      4*n + 2*n^3*s + 26*n^2 + (44*n^3)/3 = 2*n + n^2 + (44*n^3)/3 + 24*n^2 + 2*n + s*2*n^3 + n^2
    !   not SchurQuotient : 
    !     with schur:
    !       Hessenberg :     4*n + 2*n^3*s + 26*n^2 + (74*n^3)/3 = 10*n^3 + 2*n + n^2 + (44*n^3)/3 + 24*n^2 + 2*n + s*2*n^3 + n^2
    !       not Hessenberg : 4*n + 2*n^3*s + 26*n^2 + 30*n^3     = (10/3)*n^3 + 10*n^3 + 2*n^3 + 2*n + n^2 + (44*n^3)/3 + 24*n^2 + 2*n + s*2*n^3 + n^2
    !     without schur:     4*n + 2*n^3*s + 26*n^2 + (44*n^3)/3 = 2*n + n^2 + (44*n^3)/3 + 24*n^2 + 2*n + s*2*n^3 + n^2
    ! 
    ! inputs and outputs
    integer,  intent(in)  :: N
    real(WP), intent(in)  :: A(N,N)
    real(WP), intent(out) :: F(N,N)
    integer :: iscale, iflag, info
    real(WP), intent(inout),  target :: rwrksp(nrLength)
    type(expmA_opts) :: option
    
    ! General control
    integer :: iworkNow, nrLength, nrLrequired

    ! Pade algorithm
    integer,parameter :: m_vals(5) = (/3, 5, 7, 9, 13/)
    real(WP) :: theta(5) = (/ &
      !0d0, &!3.650024139523051e-008
      !0d0, &!5.317232856892575e-004
      1.495585217958292d-002,  &! m_vals = 3
      !0d0, &!8.536352760102745e-002
      2.539398330063230d-001,  &! m_vals = 5
      !0d0, &!5.414660951208968e-001
      9.504178996162932d-001,  &! m_vals = 7
      !0d0, &!1.473163964234804e+000
      2.097847961257068d+000,  &! m_vals = 9
      !0d0, &!2.811644121620263e+000
      !0d0, &!3.602330066265032e+000
      !0d0, &!4.458935413036850e+000
      5.371920351148152d+000  &! m_vals = 13
      /)
    real(WP) :: normT, t
    integer  :: i,j,s

    
    ! Lapack arrays
    integer :: lwork, info_lpk
    real(WP), dimension(:), pointer   :: Tau, work, eigr, eigi
    real(WP), dimension(:,:), pointer :: UM, TM, X2, X3, X4, U
    
    ! :::::::::: start :::::::::::

    !------- input check 
    iflag = 0
    iscale = 0
    if (.not. option%opt_init) stop 'Please call expmSetOptions first.'
    
    ! length of the double workspace
    if (option%performSchur) then
      nrLrequired = 7*N*N + 78*N + 4160 ! 5*N*N + (N*N+3*N) + (N*N + 75*N + 4160)
    else
      nrLrequired = 5*N*N
    end if
    if (nrLength .lt. nrLrequired) then
      write(*,"('Workspace is not large enough (now/required): ',I8,I8)") nrLength, nrLrequired
      stop
    end if

    !------- workspace pointers 
    iworkNow = 1
    X2  => assignWorkspace(iworkNow, nrLength, rwrksp, N, N) ! nrLrequired:  N*N
    X3  => assignWorkspace(iworkNow, nrLength, rwrksp, N, N) ! nrLrequired:  N*N
    X4  => assignWorkspace(iworkNow, nrLength, rwrksp, N, N) ! nrLrequired:  N*N
    U   => assignWorkspace(iworkNow, nrLength, rwrksp, N, N) ! nrLrequired:  N*N
    TM => assignWorkspace(iworkNow, nrLength, rwrksp, N, N)  ! nrLrequired:  N*N
    
    !------- perform schur decomposiiton
    if (option%performSchur) then
      UM => assignWorkspace(iworkNow, nrLength, rwrksp, N, N)   ! nrLrequired:  N*N
      eigr => assignWorkspace(iworkNow, nrLength, rwrksp, N)    ! nrLrequired:  N
      eigi => assignWorkspace(iworkNow, nrLength, rwrksp, N)    ! nrLrequired:  N
      Tau  => assignWorkspace(iworkNow, nrLength, rwrksp, N)    ! nrLrequired:  N
      lwork = N*N + 75*N + 4160
      work => assignWorkspace(iworkNow, nrLength, rwrksp, lwork) ! nrLrequired:  lwork = N*N + 75*N + 4160
      call realSchur (N, A, N, UM, N, TM, N, eigr, eigi, Tau, lwork, work, option%Hessenberg, iflag)
      work => nullWorkspace(iworkNow, lwork, 'work')           ! nrLrequired: -lwork = N*N + 75*N + 4160
      Tau  => nullWorkspace(iworkNow, N, 'Tau')                ! nrLrequired: -N
      if (iflag .ne. 0) then
        iflag = 1
        if (option%performSchur) then
            ! F = UM*F*UM^T  
            !call DGEMM('N','N',N, N, N, 1d0, UM, N, F, N, 0d0, U, N) ! F = UM*F ! dont "call dgemm(...UM,...,F,...,F,...)"
            !call DGEMM('N','T',N, N, N, 1d0, U, N, UM, N, 0d0, F, N) ! F = F*UM^T
            F = Matmul_blas (N, UM, F)
            F = MatmulT_blas(N, F, UM)
            eigi => nullWorkspace(iworkNow, N, 'eigi')    ! nrLrequired:  -N
            eigr => nullWorkspace(iworkNow, N, 'eigr')    ! nrLrequired:  -N
            UM => nullWorkspace(iworkNow, N, N, 'UM')     ! nrLrequired:  -N*N
        end if
        ! clear pointers
        TM  => nullWorkspace(iworkNow, N, N, 'TM')  ! nrLrequired:  -N*N
        U   => nullWorkspace(iworkNow, N, N, 'U')   ! nrLrequired:  -N*N
        X4  => nullWorkspace(iworkNow, N, N, 'X4')  ! nrLrequired:  -N*N
        X3  => nullWorkspace(iworkNow, N, N, 'X3')  ! nrLrequired:  -N*N
        X2  => nullWorkspace(iworkNow, N, N, 'X2')  ! nrLrequired:  -N*N
        return
      end if
    else
      ! TM = A
      call DLACPY( 'A', N, N, A, N, TM, N )
    end if

    !------- norm or spectral radius
    if (option%SchurQuotient .or. option%performSchur) then
        !flops: 2*n
        normT = abs(TM(1,1))
        do i = 2, N
            if (TM(i,i) .eq. TM(i-1,i-1)) then
                normT = max(sqrt(TM(i,i)*TM(i,i) +abs(TM(i-1,i)*TM(i,i-1))), normT)
            else
                normT = max(abs(TM(i,i)),normT)
            end if
        end do             
    else
        normT = Mat1norm(N, N, TM, N) !flops: 2*n
    end if

    !------- Pade approximation

    if (normT .le. theta(5)) then

      ! no scaling and squaring is required.
      do i = 1, 5
        if(normT .lt. theta(i))then
          !    PadeDegree (X, m, Y, X2, X3, X4, U)
          call PadeDegree(TM, m_vals(i), F, X2, X3, X4, U)
          exit
        end if
      end do

    else

      ! scaling
      call log2fr(normT/theta(5), t, s)
      if (t .eq. 0.5d0) then
        s = s - 1
      end if
      iscale = s   
      ! rtmp = 1d0/dble(2**s)
      
      ! TM = TM / dble(2**s)
      ! flops n^2
      if (option%SchurQuotient .or. option%performSchur) then
          call DLASCL( 'H', 0, 0, dble(2**s), 1.0_WP, N, N, TM, N, info )
      else
          call DLASCL( 'G', 0, 0, dble(2**s), 1.0_WP, N, N, TM, N, info )
      end if
      if (info .ne. 0) then
          write(6,"('Dimension and number of scale ', 2I8)") N, s
      end if
      
      call PadeDegree(TM, m_vals(5), F, X2, X3, X4, U)

      ! squaring
      ! flops: s*2*n^3
      do j = 1, s
        F = Matmul_blas(N, F, F)
      end do 

    end if

    if (option%performSchur) then
      ! F = UM*F*UM^T  
      !call DGEMM('N','N',N, N, N, 1d0, UM, N, F, N, 0d0, U, N) ! F = UM*F ! dont "call dgemm(...UM,...,F,...,F,...)"
      !call DGEMM('N','T',N, N, N, 1d0, U, N, UM, N, 0d0, F, N) ! F = F*UM^T
      F = Matmul_blas (N, UM, F)
      F = MatmulT_blas(N, F, UM)
      eigi => nullWorkspace(iworkNow, N, 'eigi')    ! nrLrequired:  -N
      eigr => nullWorkspace(iworkNow, N, 'eigr')    ! nrLrequired:  -N
      UM => nullWorkspace(iworkNow, N, N, 'UM')     ! nrLrequired:  -N*N

    end if

    ! clear pointers
    TM  => nullWorkspace(iworkNow, N, N, 'TM')  ! nrLrequired:  -N*N
    U   => nullWorkspace(iworkNow, N, N, 'U')   ! nrLrequired:  -N*N
    X4  => nullWorkspace(iworkNow, N, N, 'X4')  ! nrLrequired:  -N*N
    X3  => nullWorkspace(iworkNow, N, N, 'X3')  ! nrLrequired:  -N*N
    X2  => nullWorkspace(iworkNow, N, N, 'X2')  ! nrLrequired:  -N*N
    
    do i = 1, N
        do j = 1, N
            if (isnan(F(i,j))) iflag = 9
        end do
    end do
    
 
    contains
    

    subroutine PadeDegree (X, m, Y, X2, X3, X4, U)
      ! Pade approximant to exponential.
      ! flops if m==13 : (44*n^3)/3 + 24*n^2 + 2*n = 6*2*n^3+12*n^2+12*n^2+2*n+(2/3)*n^3+2*n^3
      implicit double precision (C)
      
      real(WP),dimension(N,N) :: X, Y, X2, X3, X4, U

      integer,dimension(N) :: ipiv      
      real(WP) :: d,uk
      integer :: i,j,k,m,info,rank
      
      ! Constant coefficients
      C31 = 120D0;
      C32 = 60D0;
      C33 = 12D0;
      C51 = 30240D0;
      C52 = 15120D0;
      C53 = 3360D0;
      C54 = 420D0;
      C55 = 30D0;
      C71 = 17297280D0;
      C72 = 8648640D0;
      C73 = 1995840D0;
      C74 = 277200D0;
      C75 = 25200D0;
      C76 = 1512D0;
      C77 = 56D0;
      C91 = 17643225600D0;
      C92 = 8821612800D0;
      C93 = 2075673600D0;
      C94 = 302702400D0;
      C95 = 30270240D0;
      C96 = 2162160D0;
      C97 = 110880D0;
      C98 = 3960D0;
      C99 = 90D0;
      C131 = 64764752532480000D0;
      C132 = 32382376266240000D0;
      C133 = 7771770303897600D0;
      C134 = 1187353796428800D0;
      C135 = 129060195264000D0;
      C136 = 10559470521600D0;
      C137 = 670442572800D0;
      C138 = 33522128640D0;
      C139 = 1323241920D0;
      C1310 = 40840800D0;
      C1311 = 960960D0;
      C1312 = 16380D0;
      C1313 = 182D0;
      C1314 = 1D0;

      ! Compute the approximation.
      ! A switch statement on m would be clearer, but this way we can share some
      ! code between the cases to do matrix powers.
      
      ! M*M = 6
      ! M+M = 12
      ! t*M = 12
      ! M+D = 2
      ! flops: 6*2*n^3+12*n^2+12*n^2+2*n
      
      X2 = Matmul_blas(N,X,X) 
      if (m .eq. 3) then
        U = X2
        do k = 1, N
          U(k,k) = U(k,k) + C32
        end do
        U = Matmul_blas(N, X, U)
        Y = C33*X2
        d = C31
      else
        X3 = Matmul_blas(N, X2, X2)
        if (m .eq. 5)then
          U = X3 + C54*X2
          do k = 1, N
            U(k,k) = U(k,k) + C52
          end do
          U = Matmul_blas(N, X, U)
          Y = C55*X3 + C53*X2
          d = C51
        else
          X4 = Matmul_blas(N, X3,X2)
          if (m .eq. 7)then
            U = X4 + C76*X3 + C74*X2
            do k = 1, N
              U(k,k) = U(k,k) +C72
            end do
            U = Matmul_blas(N, X, U)
            Y = C77*X4 + C75*X3 + C73*X2
            d = C71
          else if (m .eq. 9) then  
            Y = Matmul_blas(N, X4,X2)
            U = Y + C98*X4 + C96*X3 + C94*X2
            do k = 1, N
                U(k,k) = U(k,k) + C92
            end do
            U = Matmul_blas(N, X, U)
            Y = C99*Y + C97*X4 + C95*X3 + C93*X2
            d = C91           
          else ! m .eq. 13  
            U = C138*X4 + C136*X3 + C134*X2
            do k = 1, N
                U(k,k) = U(k,k) + C132
            end do
            Y = C1314*X4 + (C1312*X3) + (C1310*X2)
            Y = Matmul_blas(N, X4, Y)
            U = U + Y
            U = Matmul_blas(N, X, U)
            Y = C1313*X4 + C1311*X3 + C139*X2
            Y = Matmul_blas(N, X4, Y)
            Y = Y + C137*X4 + C135*X3 + C133*X2
            d = C131           
          end if
        end if
      end if        
      
      do k = 1, N
        Y(k,k) = Y(k,k) + d
      end do
      do j = 1, N
        do i = 1, N
          uk = U(i,j)
          U(i,j) = Y(i,j) - uk
          Y(i,j) = Y(i,j) + uk  
        end do
      end do

      ! Y = U\Y
      ! solve U*Y = Y   
      
      ! SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
      ! SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
      !  DGETRS solves a system of linear equations
      !     A * X = B  or  A' * X = B
      !  with a general N-by-N matrix A using the LU factorization computed
      !  by DGETRF.      
      call DGETRF(N, N, U, N, ipiv, info)                     ! flops: (2/3)*n^3
      if (info .ne. 0) stop 'Error in DGETRF : PadeDegree'
      call DGETRS('N',N, N, U, N, ipiv, Y, N, info)           ! flops: 2*n^3
      if (info .ne. 0) stop 'Error in DGETRS : PadeDegree'

    end subroutine PadeDegree

    ! contains end

  end subroutine dexpm_pade
  !#########################################

end module expmA
