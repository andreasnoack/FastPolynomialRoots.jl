 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Aurentz² Mach³ Vandebril³ Watkins²
!
! ²Dept. Mathematics, Washington State University
! ³Dept. Computer Science, KU Leuven
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Last modified 22 August 2014
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Test of interesting polynomias
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! arguments
!
!   mi     start index (no. polynom)
!          default  1  
!
!   mj     stop index (no. polynom)
!          default  48  
!
!   mk     no. of runs
!          default  1  
!
!   nobalnceor   0  balance some
!                1  balance all
!                2  balance none
!          default  0  
!
! 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! No.\ & Description   &  Deg.\ & Roots   % number in rootrace
! 
! 1  & Wilkinson polynomial & 10 & 1,...,10
! 2  & Wilkinson polynomial & 15 & 1,...,15
! 3  & Wilkinson polynomial & 20 & 1,...,20
! 4  & scaled and shifted Wilkinson poly. & 20 & -2.1,-1.9,...,1.7
! 5  & reverse Wilkinson polynomial & 10 & 1,1/2,...,1/10
! 6  & reverse Wilkinson polynomial & 15 & 1,1/2,...,1/15
! 7  & reverse Wilkinson polynomial & 20 & 1,1/2,...,1/20
! 8  & prescribed roots of varying scale & 20 & 2^{-10},2^{-9},...,2^{9}
! 9  & prescribed roots of varying scale -3 & 20 & 2^{-10}-3,...,2^{9}-3
! 10 & Chebyshev polynomial & 20 & \cos(\frac{2j-1}{40}\pi)
! 11 & z^{20} + z^{19}+ ... + z + 1& 20 & \cos(\frac{2j}{21}\pi)
! 
! MPSolve
! 12 & trv\_m, C.\ Traverso& 24 & known
! 13 & mand31 Mandelbrot example (k=5)&31&known
! 14 & mand63 Mandelbrot example (k=6)&63&known
! 
! 15 & polynomial from V.\ Noferini & 12 & almost random, z_{1} scaled with 1e+12, z_{2} with 1e+9
! 16 & polynomial from V.\ Noferini & 35 & almost random, z_{1} scaled with 1e+12, z_{2} with 1e+9
! 
! Jenkins and Traub 1970
! 17 & p_{1}(z) with a=1e-8 & 3 & 1e-8, -1e-8, 1
! 18 & p_{1}(z) with a=1e-15 & 3 & 1e-15, -1e-15, 1
! 19 & p_{1}(z) with a=1e+8 & 3 & 1e+8, -1e+8, 1
! 20 & p_{1}(z) with a=1e+15 & 3 & 1e+15, -1e+15, 1
! 21 & p_{3}(z) underflow test & 10 & 1e-1,...,1e-10  
! 22 & p_{3}(z) underflow test & 20 & 1e-1,...,1e-20  
! 23 & p_{10}(z) deflation test a=10e+3 & 3 & 1, 1e+3, 1e-3  
! 24 & p_{10}(z) deflation test a=10e+6 & 3 & 1, 1e+6, 1e-6  
! 25 & p_{10}(z) deflation test a=10e+9 & 3 & 1, 1e+9, 1e-9 
! 26 & p_{11}(z) deflation test m=15 & 60 & \exp(\frac{ik\pi}{2m}), 0.9\exp(\frac{ik\pi}{2m})
! 
! 27 & Bernoulli polynomial (k=20) & 20 & ---4
! 28 & truncated exponential (k=20)& 20 & ---3
! 
! used in Bevilacqua, Del Corso, and Gemignani 2014
! 29--33 & p_{1}(z) with m = 10, 20, 30, 256, 512 & 2m & ---   
! 34--38 & p_{2}(z) with m = 10, 20, 30, 256, 512 & 2m & ---   
! 39--43 & p_{3}(z) with m+1 = 20, 40, 60, 512, 1024, \lambda = 0.9 & m+1 & ---   
! 44--48 & p_{3}(z) with m+1 = 20, 40, 60, 512, 1024, \lambda = 0.999 & m+1 & ---   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Remark: In the paper we include a comparison
!         with BEGG (and BBEGG). Since we cannot
!         redistribute their code, these 
!         comparisons have been removed.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program main 
  
  use mpmodule
  implicit none 
  
  complex(kind(1d0)), allocatable :: poly(:),roots(:),allroots(:,:),npoly(:),poly2(:)
  complex(kind(1d0)), allocatable :: B(:,:),Y(:,:),cwork(:)  
  integer, allocatable :: iterations(:),its(:), seed(:)
  double precision, allocatable :: res(:,:), err(:), relerr(:)
  double precision :: time,begg_time,lapack_time,pdble,error(9),alpha, hd
  double precision, allocatable :: Q(:),D(:),C(:),B2(:)
  double precision, allocatable :: rcoeffs(:),icoeffs(:),rwork(:)
  double precision, allocatable :: reigs(:),ieigs(:),rnew(:),inew(:)
  complex(kind(1d0)), allocatable :: eigs(:)
  double precision, allocatable ::residuals(:,:)
  character (len=100) :: Desc(200)

  complex(kind(1d0)) :: coeffs

  integer :: ii,jj,kk,ll,mm,N,zero,flag,info=0,nnew, kj, kl,mj,mi,mk,ry, rsize
  integer :: klm
  integer :: clock_start,clock_end,clock_rate, newtnum
  integer :: Deg(200),num_trials(200), lapack_trials
  character (len=*), parameter :: path = "./"

  integer :: eigsknown = 0, nobalance = 0, nobalanceor = 0
  complex(kind(1d0)), allocatable :: exacteigs(:)
  integer, allocatable :: eigtaken(:)
  double precision :: pi = 3.14159265358979323d0, hA, lambda, normofp
  double precision :: relforwarderror(5), backwarderror(5), back2(5), back3(5)
  double precision, allocatable :: xr_db(:),xi_db(:),c_db(:)

  character(len=32) :: arg

  ! BLAS
  integer :: idamax
  double precision :: dznrm2


  call mpinit

  if (iargc()>0) then
     call getarg(1, arg)
     read (arg,'(I10)') mi
     if (iargc()>1) then
        call getarg(2, arg)
        read (arg,'(I10)') mj
        if (iargc()>2) then
           call getarg(3, arg)
           read (arg,'(I10)') mk
           if (iargc()>3) then
              call getarg(4, arg)
              read (arg,'(I10)') nobalanceor
              if (iargc()>4) then
                 call getarg(5, arg)
                 read (arg,'(I10)') ry
              else
                 ry = 0
              end if
           else
              nobalanceor = 0
              ry = 0
           end if
        else
           ry = 0
           nobalanceor = 0
           mk = 1
        end if
     else
        ry = 0
        nobalanceor = 0
        mk = 1              
        mj = mi
     end if
  else
     ry = 0
     nobalanceor = 0
     mk = 1
     mj = 27
     mi = 1
  end if
  
  open (unit=7, file="sp_degrees.txt", status='unknown', position="append")
  open (unit=31, file="sp_table1.txt", status='unknown', position="append")
  open (unit=32, file="sp_table2.txt", status='unknown', position="append")
  
  Deg = 0
  Desc = ""

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Wilkinson Polynomial degree 10
  ! roots 1,2,...,10
  num_trials(1) = 2**(0) 
  Deg(1) = 10
  Desc(1) = "Wilkinson polynomial"
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Wilkinson Polynomial degree 15
  ! roots 1,2,...,15
  num_trials(2) = 2**(0) 
  Deg(2) = 15
  Desc(2) = "Wilkinson polynomial"
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Wilkinson Polynomial
  ! roots 1,2,...,20
  num_trials(3) = 2**(0) 
  Deg(3) = 20
  Desc(3) = "Wilkinson Polynomial"
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! monic polynomial roots [-2.1:0.2:1.7]
  num_trials(4) = 2**(0)
  Deg(4) = 20
  Desc(4) = "roots at [-2.1:0.2:1.7]"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  num_trials(5) = 2**(0)
  Deg(5) = 10
  Desc(5) = "reverse Wilkinson degree 10"
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  num_trials(6) = 2**(0)
  Deg(6) = 15
  Desc(6) = "reverse Wilkinson degree 15"
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  num_trials(7) = 2**(0)
  Deg(7) = 20
  Desc(7) = "reverse Wilkinson degree 20"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! univariate polynomial with zeros 2^-10, 2^-9, ..., 2^8, 2^9
  num_trials(8) = 2**(0)
  Deg(8) = 20
  Desc(8) = "univariate polynomial with zeros 2^-10, 2^-9, ..., 2^8, 2^9"

  ! univariate polynomial with zeros 2^-10, 2^-9, ..., 2^8, 2^9 shifted
  num_trials(9) = 2**(0)
  Deg(9) = 20
  Desc(9) = "univariate polynomial with zeros 2^-10, 2^-9, ..., 2^8, 2^9 shifted with 3"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Chebyshev polynomial of degree 20 
  num_trials(10) = 2**(0)
  Deg(10) = 20
  Desc(10) = "Chebyshev polynomial"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! z^20 + z^19 + ... + z + 1
  num_trials(11) = 2**(0)
  Deg(11) = 20
  Desc(11) = "z^20 + z^19 + ... + z + 1"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! C. Traverso 24 MPSolve
  num_trials(12) = 2**(0) 
  Deg(12) = 24
  Desc(12) = "trv_m 24 example from MPSolve"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Mandelbrot 31 MPSolve
  num_trials(13) = 2**(0) 
  Deg(13) = 31
  Desc(13) = "Mandelbrot example from MPSolve"
  ! Mandelbrot 63 MPSolve
  num_trials(14) = 2**(0) 
  Deg(14) = 63
  Desc(14) = "Mandelbrot example from MPSolve"
  ! we are not good at mandelbrot 63, others are useless

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  num_trials(15) = 2**(0)
  Deg(15) = 12
  Desc(15) = "Vanni Noferini's example degree 12"
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  num_trials(16) = 2**(0)
  Deg(16) = 35
  Desc(16) = "Vanni Noferini's example degree 35"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  num_trials(17) = 2**(0)
  Deg(17) = 3
  Desc(17) = "cubic polynomial small roots"
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  num_trials(18) = 2**(0)
  Deg(18) = 3
  Desc(18) = "cubic polynomial very small roots"
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  num_trials(19) = 2**(0)
  Deg(19) = 3
  Desc(19) = "cubic polynomial large roots"
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  num_trials(20) = 2**(0)
  Deg(20) = 3
  Desc(20) = "cubic polynomial very large roots" 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  num_trials(21) = 2**(0)
  Deg(21) = 10
  Desc(21) = "underflow test degree 10"
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  num_trials(22) = 2**(0)
  Deg(22) = 20
  Desc(22) = "underflow test degree 20"
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  num_trials(23) = 2**(0)
  Deg(23) = 3
  Desc(23) = "deflation test A = 10d3"
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  num_trials(24) = 2**(0)
  Deg(24) = 3
  Desc(24) = "deflation test A = 10d6"
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  num_trials(25) = 2**(0)
  Deg(25) = 3
  Desc(25) = "deflation test A = 10d9"
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  num_trials(26) = 2**(0)
  Deg(26) = 60
  Desc(26) = "deflation test M = 15"
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! ROOTS UNKNOWN
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Bernoulli polynomial of degree 20  
  num_trials(27) = 2**(0)
  Deg(27) = 20
  Desc(27) = "Bernoulli polynomial"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! p(z) = (20!) sum_{k=0}^{20} z^k/k! ! truncated exponential degree 20
  num_trials(28) = 2**(0)
  Deg(28) = 20
  Desc(28) = "truncated exp"



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  num_trials(29) = 2**(0)
  Deg(29) = 20
  Desc(29) = "BDCG P1 N = 20"
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  num_trials(30) = 2**(0)
  Deg(30) = 40
  Desc(30) = "BDCG P1 N = 40"
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  num_trials(31) = 2**(0)
  Deg(31) = 60
  Desc(31) = "BDCG P1 N = 60"
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  num_trials(32) = 2**(0)
  Deg(32) = 512
  Desc(32) = "BDCG P1 N = 512"
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  num_trials(33) = 2**(0)
  Deg(33) = 1024
  Desc(33) = "BDCG P1 N = 1024"
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  num_trials(34) = 2**(0)
  Deg(34) = 20
  Desc(34) = "BDCG P2 N = 20"
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  num_trials(35) = 2**(0)
  Deg(35) = 40
  Desc(35) = "BDCG P2 N = 40"
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  num_trials(36) = 2**(0)
  Deg(36) = 60
  Desc(36) = "BDCG P2 N = 60"
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  num_trials(37) = 2**(0)
  Deg(37) = 512
  Desc(37) = "BDCG P2 N = 512"
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  num_trials(38) = 2**(0)
  Deg(38) = 1024
  Desc(38) = "BDCG P2 N = 1024"
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  num_trials(39) = 2**(0)
  Deg(39) = 20
  Desc(39) = "BDCG P3 N = 20"
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  num_trials(40) = 2**(0)
  Deg(40) = 40
  Desc(40) = "BDCG P3 N = 40"
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  num_trials(41) = 2**(0)
  Deg(41) = 60
  Desc(41) = "BDCG P3 N = 60"
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  num_trials(42) = 2**(0)
  Deg(42) = 512
  Desc(42) = "BDCG P3 N = 512"
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  num_trials(43) = 2**(0)
  Deg(43) = 1024
  Desc(43) = "BDCG P3 N = 1024"
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  num_trials(44) = 2**(0)
  Deg(44) = 20
  Desc(44) = "BDCG P3 N = 20"
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  num_trials(45) = 2**(0)
  Deg(45) = 40
  Desc(45) = "BDCG P3 N = 40"
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  num_trials(46) = 2**(0)
  Deg(46) = 60
  Desc(46) = "BDCG P3 N = 60"
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  num_trials(47) = 2**(0)
  Deg(47) = 512
  Desc(47) = "BDCG P3 N = 512"
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  num_trials(48) = 2**(0)
  Deg(48) = 1024
  Desc(48) = "BDCG P3 N = 1024"
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  

  call RANDOM_SEED(size = rsize)
  allocate(seed(rsize))
  call RANDOM_SEED(GET = seed)
  ii = 1
  seed(1) = ii
  seed(2) = ii+1000000
  seed(3) = ii+2100000
  seed(4) = ii+3210000
  seed(5) = ii+43210000
  seed(6) = ii+5432100
  seed(7) = ii+6543210
  seed(8) = ii+7654321
  seed(9) = ii+8765432
  seed(10) = ii+9876543
  seed(11) = ii+10987654
  seed(12) = ii+11109876
  call RANDOM_SEED(PUT = seed)

  do kk=mi,mj
     !num_trials(kk) = 2**(10-kk)
     !if (num_trials(kk)<1) then
     if (mk>0) then
        num_trials(kk) = mk
     end if
     !end if
     !Deg(kk) = 2**(kk+1)+1
     write (7,*) Deg(kk), num_trials(kk), Desc(kk)
  end do
  close(7)

  do ll=mi,mj
     
     if (Deg(ll)>0) then
        
        time = 0d0
        
        N = Deg(ll)
        
        write(*,*) "Current N =",N
        
        ! allocate memory
        allocate(Q(3*n),D(2*(n+1)),C(3*n),B2(3*n),rcoeffs(n))
        allocate(icoeffs(n),rnew(n),inew(n),npoly(n))
        allocate(its(n),reigs(n),ieigs(n))
        allocate(residuals(n,3*(newtnum+1)))
     
     allocate(poly(N),roots(N),allroots(N,2),iterations(N),res(N,6))
     allocate(B(N,N),Y(N,N),cwork(5*N),rwork(2*N))

     allocate(exacteigs(N),eigtaken(N))
     allocate(poly2(N),err(N),relerr(N))
     allocate(xr_db(N),xi_db(N),c_db(N))

     nobalance = 0
     rcoeffs = 0d0
     icoeffs = 0d0

     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     select case (ll)


        case (1) ! Wilkinson Polynomial !check
           print*, "Wilkinson Polynomial degree 10"
           eigsknown = 1
           exacteigs(1)  = complex(1d0,0d0)
           exacteigs(2)  = complex(2d0,0d0)
           exacteigs(3)  = complex(3d0,0d0)
           exacteigs(4)  = complex(4d0,0d0)
           exacteigs(5)  = complex(5d0,0d0)
           exacteigs(6)  = complex(6d0,0d0)
           exacteigs(7)  = complex(7d0,0d0)
           exacteigs(8)  = complex(8d0,0d0)
           exacteigs(9)  = complex(9d0,0d0)
           exacteigs(10) = complex(10d0,0d0)
           do ii=1,n
              xr_db(ii) = dble(exacteigs(ii))
              xi_db(ii) = 0d0
           end do
           call rootstocoeffs(n,xr_db,xi_db,c_db)
           RCOEFFS = c_db

        case (2) ! Wilkinson Polynomial !check
           print*, "Wilkinson Polynomial degree 15"
           eigsknown = 1
           exacteigs(1)  = complex(1d0,0d0)
           exacteigs(2)  = complex(2d0,0d0)
           exacteigs(3)  = complex(3d0,0d0)
           exacteigs(4)  = complex(4d0,0d0)
           exacteigs(5)  = complex(5d0,0d0)
           exacteigs(6)  = complex(6d0,0d0)
           exacteigs(7)  = complex(7d0,0d0)
           exacteigs(8)  = complex(8d0,0d0)
           exacteigs(9)  = complex(9d0,0d0)
           exacteigs(10) = complex(10d0,0d0)
           exacteigs(11) = complex(11d0,0d0)
           exacteigs(12) = complex(12d0,0d0)
           exacteigs(13) = complex(13d0,0d0)
           exacteigs(14) = complex(14d0,0d0)
           exacteigs(15) = complex(15d0,0d0)
           do ii=1,n
              xr_db(ii) = dble(exacteigs(ii))
              xi_db(ii) = imag(exacteigs(ii))
           end do
           call rootstocoeffs(n,xr_db,xi_db,c_db)
           RCOEFFS = c_db

        case (3) ! Wilkinson Polynomial
           print*, "Wilkinson Polynomial"
           eigsknown = 1
           exacteigs(1)  = complex(1d0,0d0)
           exacteigs(2)  = complex(2d0,0d0)
           exacteigs(3)  = complex(3d0,0d0)
           exacteigs(4)  = complex(4d0,0d0)
           exacteigs(5)  = complex(5d0,0d0)
           exacteigs(6)  = complex(6d0,0d0)
           exacteigs(7)  = complex(7d0,0d0)
           exacteigs(8)  = complex(8d0,0d0)
           exacteigs(9)  = complex(9d0,0d0)
           exacteigs(10) = complex(10d0,0d0)
           exacteigs(11) = complex(11d0,0d0)
           exacteigs(12) = complex(12d0,0d0)
           exacteigs(13) = complex(13d0,0d0)
           exacteigs(14) = complex(14d0,0d0)
           exacteigs(15) = complex(15d0,0d0)
           exacteigs(16) = complex(16d0,0d0)
           exacteigs(17) = complex(17d0,0d0)
           exacteigs(18) = complex(18d0,0d0)
           exacteigs(19) = complex(19d0,0d0)
           exacteigs(20) = complex(20d0,0d0)
           do ii=1,n
              xr_db(ii) = dble(exacteigs(ii))
              xi_db(ii) = 0d0
           end do
           call rootstocoeffs(n,xr_db,xi_db,c_db)
           RCOEFFS = c_db

        case (4)  ! monic polynomial roots [-2.1:0.2:1.7]
           print*, "monic polynomial roots [-2.1:0.2:1.7]"
           eigsknown = 1           
           exacteigs(1)  = complex(-2.1d0,0d0)
           exacteigs(2)  = complex(-1.9d0,0d0)
           exacteigs(3)  = complex(-1.7d0,0d0)
           exacteigs(4)  = complex(-1.5d0,0d0)
           exacteigs(5)  = complex(-1.3d0,0d0)
           exacteigs(6)  = complex(-1.1d0,0d0)
           exacteigs(7)  = complex(-0.9d0,0d0)
           exacteigs(8)  = complex(-0.7d0,0d0)
           exacteigs(9)  = complex(-0.5d0,0d0)
           exacteigs(10) = complex(-0.3d0,0d0)
           exacteigs(11) = complex(-0.1d0,0d0)
           exacteigs(12) = complex(0.1d0,0d0)
           exacteigs(13) = complex(0.3d0,0d0)
           exacteigs(14) = complex(0.5d0,0d0)
           exacteigs(15) = complex(0.7d0,0d0)
           exacteigs(16) = complex(0.9d0,0d0)
           exacteigs(17) = complex(1.1d0,0d0)
           exacteigs(18) = complex(1.3d0,0d0)
           exacteigs(19) = complex(1.5d0,0d0)
           exacteigs(20) = complex(1.7d0,0d0)
           do ii=1,n
              xr_db(ii) = dble(exacteigs(ii))
              xi_db(ii) = 0d0
           end do
           call rootstocoeffs(n,xr_db,xi_db,c_db)
           RCOEFFS = c_db

        case (5,6,7)  ! inverse Wilkinson polynomial degree 20
           print*, "inverse Wilkinson polynomial degree", Deg(ll)
           do ii=1,n
              exacteigs(ii) = complex(1d0/ii,0d0)
           end do
           eigsknown = 1
           nobalance = 1
           do ii=1,n
              xr_db(ii) = dble(exacteigs(ii))
              xi_db(ii) = 0d0
           end do
           call rootstocoeffs(n,xr_db,xi_db,c_db)
           RCOEFFS = c_db

        case (8)  ! univariate polynomial with zeros 2^-10, 2^-9, ..., 2^8, 2^9
           print*, "univariate polynomial with zeros 2^-10, 2^-9, ..., 2^8, 2^9"         
           eigsknown = 1
           nobalance = 1
           exacteigs(1)  = complex(1d0/1024d0,0d0)
           exacteigs(2)  = complex(1d0/512d0,0d0)
           exacteigs(3)  = complex(1d0/256d0,0d0)
           exacteigs(4)  = complex(1d0/128d0,0d0)
           exacteigs(5)  = complex(1d0/64d0,0d0)
           exacteigs(6)  = complex(1d0/32d0,0d0)
           exacteigs(7)  = complex(0.0625d0,0d0)
           exacteigs(8)  = complex(0.125d0,0d0)
           exacteigs(9)  = complex(0.25d0,0d0)
           exacteigs(10) = complex(0.5d0,0d0)
           exacteigs(11) = complex(1d0,0d0)
           exacteigs(12) = complex(2d0,0d0)
           exacteigs(13) = complex(4d0,0d0)
           exacteigs(14) = complex(8d0,0d0)
           exacteigs(15) = complex(16d0,0d0)
           exacteigs(16) = complex(32d0,0d0)
           exacteigs(17) = complex(64d0,0d0)
           exacteigs(18) = complex(128d0,0d0)
           exacteigs(19) = complex(256d0,0d0)
           exacteigs(20) = complex(512d0,0d0)
           do ii=1,n
              xr_db(ii) = dble(exacteigs(ii))
              xi_db(ii) = imag(exacteigs(ii))
           end do
           call rootstocoeffs(n,xr_db,xi_db,c_db)
           print*, c_db
           RCOEFFS = c_db

        case (9)  ! univariate polynomial with zeros 2^-10, 2^-9, ..., 2^8, 2^9
           print*, "univariate polynomial with zeros 2^-10, 2^-9, ..., 2^8, 2^9"
           eigsknown = 1
           nobalance = 1
           exacteigs(1)  = complex(1d0/1024d0-3d0,0d0)
           exacteigs(2)  = complex(1d0/512d0-3d0,0d0)
           exacteigs(3)  = complex(1d0/256d0-3d0,0d0)
           exacteigs(4)  = complex(1d0/128d0-3d0,0d0)
           exacteigs(5)  = complex(1d0/64d0-3d0,0d0)
           exacteigs(6)  = complex(1d0/32d0-3d0,0d0)
           exacteigs(7)  = complex(0.0625d0-3d0,0d0)
           exacteigs(8)  = complex(0.125d0-3d0,0d0)
           exacteigs(9)  = complex(0.25d0-3d0,0d0)
           exacteigs(10) = complex(0.5d0-3d0,0d0)
           exacteigs(11) = complex(1d0-3d0,0d0)
           exacteigs(12) = complex(2d0-3d0,0d0)
           exacteigs(13) = complex(4d0-3d0,0d0)
           exacteigs(14) = complex(8d0-3d0,0d0)
           exacteigs(15) = complex(16d0-3d0,0d0)
           exacteigs(16) = complex(32d0-3d0,0d0)
           exacteigs(17) = complex(64d0-3d0,0d0)
           exacteigs(18) = complex(128d0-3d0,0d0)
           exacteigs(19) = complex(256d0-3d0,0d0)
           exacteigs(20) = complex(512d0-3d0,0d0)
           do ii=1,n
              xr_db(ii) = dble(exacteigs(ii))
              xi_db(ii) = imag(exacteigs(ii))
           end do
           call rootstocoeffs(n,xr_db,xi_db,c_db)
           print*, c_db
           RCOEFFS = c_db

        case (10) ! Chebyshev polynomial of degree 20 
           print*, "Chebyshev polynomial of degree 20"
           RCOEFFS(1) = 0d0
           RCOEFFS(2) =-2621440d0/524288d0
           RCOEFFS(3) = 0d0
           RCOEFFS(4) =+5570560d0/524288d0
           RCOEFFS(5) = 0d0
           RCOEFFS(6) =-6553600d0/524288d0
           RCOEFFS(7) = 0d0
           RCOEFFS(8) =+4659200d0/524288d0
           RCOEFFS(9) = 0d0
           RCOEFFS(10)=-2050048d0/524288d0
           RCOEFFS(11)= 0d0
           RCOEFFS(12)=+549120d0/524288d0
           RCOEFFS(13)= 0d0
           RCOEFFS(14)=-84480d0/524288d0
           RCOEFFS(15)= 0d0
           RCOEFFS(16)=+6600d0/524288d0
           RCOEFFS(17)= 0d0
           RCOEFFS(18)=-200d0/524288d0
           RCOEFFS(19)= 0d0
           RCOEFFS(20)= 1d0/524288d0   !1
           eigsknown = 1
           exacteigs(1)  = complex(cos(1d0/40d0*pi),0d0)
           exacteigs(2)  = complex(cos(3d0/40d0*pi),0d0)
           exacteigs(3)  = complex(cos(5d0/40d0*pi),0d0)
           exacteigs(4)  = complex(cos(7d0/40d0*pi),0d0)
           exacteigs(5)  = complex(cos(9d0/40d0*pi),0d0)
           exacteigs(6)  = complex(cos(11d0/40d0*pi),0d0)
           exacteigs(7)  = complex(cos(13d0/40d0*pi),0d0)
           exacteigs(8)  = complex(cos(15d0/40d0*pi),0d0)
           exacteigs(9)  = complex(cos(17d0/40d0*pi),0d0)
           exacteigs(10)  = complex(cos(19d0/40d0*pi),0d0)
           exacteigs(11)  = complex(cos(21d0/40d0*pi),0d0)
           exacteigs(12)  = complex(cos(23d0/40d0*pi),0d0)
           exacteigs(13)  = complex(cos(25d0/40d0*pi),0d0)
           exacteigs(14)  = complex(cos(27d0/40d0*pi),0d0)
           exacteigs(15)  = complex(cos(29d0/40d0*pi),0d0)
           exacteigs(16)  = complex(cos(31d0/40d0*pi),0d0)
           exacteigs(17)  = complex(cos(33d0/40d0*pi),0d0)
           exacteigs(18)  = complex(cos(35d0/40d0*pi),0d0)
           exacteigs(19)  = complex(cos(37d0/40d0*pi),0d0)
           exacteigs(20)  = complex(cos(39d0/40d0*pi),0d0)         

        case (11)   ! z^20 + z^19 + ... + z + 1
           print*, "z^20 + z^19 + ... + z + 1"
           RCOEFFS(20)=1d0
           RCOEFFS(19)=1d0
           RCOEFFS(18)=1d0
           RCOEFFS(17)=1d0
           RCOEFFS(16)=1d0
           RCOEFFS(15)=1d0
           RCOEFFS(14)=1d0
           RCOEFFS(13)=1d0
           RCOEFFS(12)=1d0
           RCOEFFS(11)=1d0
           RCOEFFS(10)=1d0
           RCOEFFS(9)=1d0
           RCOEFFS(8)=1d0
           RCOEFFS(7)=1d0
           RCOEFFS(6)=1d0
           RCOEFFS(5)=1d0
           RCOEFFS(4)=1d0
           RCOEFFS(3)=1d0
           RCOEFFS(2)=1d0
           RCOEFFS(1)=1d0           
           eigsknown = 1
           ! roots of unity 21 without 1
           exacteigs(1)  = complex(cos(2d0/21d0*pi),sin(2d0/21d0*pi))
           exacteigs(2)  = complex(cos(4d0/21d0*pi),sin(4d0/21d0*pi))
           exacteigs(3)  = complex(cos(6d0/21d0*pi),sin(6d0/21d0*pi))
           exacteigs(4)  = complex(cos(8d0/21d0*pi),sin(8d0/21d0*pi))
           exacteigs(5)  = complex(cos(10d0/21d0*pi),sin(10d0/21d0*pi))
           exacteigs(6)  = complex(cos(12d0/21d0*pi),sin(12d0/21d0*pi))
           exacteigs(7)  = complex(cos(14d0/21d0*pi),sin(14d0/21d0*pi))
           exacteigs(8)  = complex(cos(16d0/21d0*pi),sin(16d0/21d0*pi))
           exacteigs(9)  = complex(cos(18d0/21d0*pi),sin(18d0/21d0*pi))
           exacteigs(10)  = complex(cos(20d0/21d0*pi),sin(20d0/21d0*pi))
           exacteigs(11)  = complex(cos(22d0/21d0*pi),sin(22d0/21d0*pi))
           exacteigs(12)  = complex(cos(24d0/21d0*pi),sin(24d0/21d0*pi))
           exacteigs(13)  = complex(cos(26d0/21d0*pi),sin(26d0/21d0*pi))
           exacteigs(14)  = complex(cos(28d0/21d0*pi),sin(28d0/21d0*pi))
           exacteigs(15)  = complex(cos(30d0/21d0*pi),sin(30d0/21d0*pi))
           exacteigs(16)  = complex(cos(32d0/21d0*pi),sin(32d0/21d0*pi)) 
           exacteigs(17)  = complex(cos(34d0/21d0*pi),sin(34d0/21d0*pi))
           exacteigs(18)  = complex(cos(36d0/21d0*pi),sin(36d0/21d0*pi))
           exacteigs(19)  = complex(cos(38d0/21d0*pi),sin(38d0/21d0*pi))
           exacteigs(20)  = complex(cos(40d0/21d0*pi),sin(40d0/21d0*pi))

        case (12)
           print*, "C. Traverso 24 MPSolve"
           RCOEFFS(24)=-54765291428198020791747503747742749163073958404455022926495744d0
           RCOEFFS(23)=-4052135566767965847649766745769409681058667331648450681896960d0
           RCOEFFS(22)=-31969984081155943263834965670035075493639295858977076674560d0
           RCOEFFS(21)=575060225471570237690073740639182419333523437771848417280d0
           RCOEFFS(20)=7337981286595499156409929740830030318565357725459415040d0
           RCOEFFS(19)=6611223380089859336490797585290455483968982077145088d0
           RCOEFFS(18)=-195514288747757987122118583800597358656801082441728d0
           RCOEFFS(17)=-726907419403715013562762609680450059293446635520d0
           RCOEFFS(16)=197178719520196724204974332265013056299335680d0
           RCOEFFS(15)=5968852409133617129605588058090797893943296d0
           RCOEFFS(14)=16576506891508825500182005531742679597056d0
           RCOEFFS(13)=23375026506968330494765978581548924928d0
           RCOEFFS(12)=2206941937668751746514177591607296d0
           RCOEFFS(11)=-75617855277818001758431020580864d0
           RCOEFFS(10)=-204797687173976372829472423936d0
           RCOEFFS(9)= -143150263927579584306872320d0
           RCOEFFS(8)=  20214880144364480233472d0
           RCOEFFS(7)=  453786251090072698880d0
           RCOEFFS(6)=  1265052493274939392d0
           RCOEFFS(5)= -968887355572224d0
           RCOEFFS(4)=  1015406084096d0
           RCOEFFS(3)= -3949133824d0
           RCOEFFS(2)=  3284992d0
           RCOEFFS(1)= -1728d0
           eigsknown = 1         
           exacteigs( 1 )=complex(-3.52d2, 0d0)
           exacteigs( 2 )=complex(-3.52d2, 0d0)
           exacteigs( 3 )=complex(-2.8371450777d2, -2.9920517772d2)
           exacteigs( 4 )=complex(-2.8371450777d2,  2.9920517772d2)
           exacteigs( 5 )=complex(-2.7867414048d2,  6.1005469197d2)
           exacteigs( 6 )=complex(-2.7867414048d2, -6.1005469197d2)
           exacteigs( 7 )=complex(-2.74892372d2, 0d0)
           exacteigs( 8 )=complex(-2.014171531d2, 0d0)
           exacteigs( 9 )=complex(-1.255366582d2, 0d0)
           exacteigs( 10 )=complex(-9.599999999d1, 0d0)
           exacteigs( 11 )=complex(-8.8692435121d1,  5.5009607430d2)
           exacteigs( 12 )=complex(-8.869243512d1, -5.5009607430d2)
           exacteigs( 13 )=complex(-1.6000000000d1, 0d0)
           exacteigs( 14 )=complex( 8.23178509855d1, 0d0)
           exacteigs( 15 )=complex( 8.8692435121d1, -5.50096074303d2)
           exacteigs( 16 )=complex( 8.8692435121d1,  5.5009607430d2)
           exacteigs( 17 )=complex( 1.9293739373d2,  1.60865921259d3)
           exacteigs( 18 )=complex( 1.929373937d2, -1.6086592125d3)
           exacteigs( 19 )=complex( 2.0141715312d2, 0d0)
           exacteigs( 20 )=complex( 2.7489237213d2, 0d0)
           exacteigs( 21 )=complex( 7.52d2, 0d0)
           exacteigs( 22 )=complex( 7.52d2, 0d0)
           exacteigs( 23 )=complex( 9.1106065d2,  1.5722d0)
           exacteigs( 24 )=complex( 9.1106065d2, -1.5722d0)

       case (13) ! Mandelbort 31 MPSolve
           print*, "Mandelbrot 31 MPSolve"
           RCOEFFS(31) = 1d0
           RCOEFFS(30) = 1d0
           RCOEFFS(29) = 2d0
           RCOEFFS(28) = 5d0
           RCOEFFS(27) = 14d0
           RCOEFFS(26) = 42d0
           RCOEFFS(25) = 100d0
           RCOEFFS(24) = 221d0
           RCOEFFS(23) = 470d0
           RCOEFFS(22) = 958d0
           RCOEFFS(21) = 1860d0
           RCOEFFS(20) = 3434d0
           RCOEFFS(19) = 6036d0
           RCOEFFS(18) = 10068d0
           RCOEFFS(17) = 15864d0
           RCOEFFS(16) = 23461d0
           RCOEFFS(15) = 32398d0
           RCOEFFS(14) = 41658d0
           RCOEFFS(13) = 49700d0
           RCOEFFS(12) = 54746d0
           RCOEFFS(11) = 55308d0
           RCOEFFS(10) = 50788d0
           RCOEFFS(9) = 41944d0
           RCOEFFS(8) = 30782d0
           RCOEFFS(7) = 19788d0
           RCOEFFS(6) = 10948d0
           RCOEFFS(5) = 5096d0
           RCOEFFS(4) = 1932d0
           RCOEFFS(3) = 568d0
           RCOEFFS(2) = 120d0
           RCOEFFS(1) = 16d0
           eigsknown = 1          
           exacteigs( 1 ) = complex(-1.996376137,0d0)
           exacteigs( 2 ) = complex(-1.966773216,0d0)
           exacteigs( 3 ) = complex(-1.907280091,0d0)
           exacteigs( 4 ) = complex(-1.772892903,0d0)
           exacteigs( 5 ) = complex(-1.754877666,0d0)
           exacteigs( 6 ) = complex(-1.47601464272,0d0)
           exacteigs( 7 ) = complex(-1.284084925525, 4.272688960406d-1)
           exacteigs( 8 ) = complex(-1.284084925525,-4.272688960406d-1)
           exacteigs( 9 ) = complex(-1.138000666650,-2.403324012620d-1)
           exacteigs( 10 )= complex(-1.138000666650, 2.403324012620d-1)
           exacteigs( 11 )= complex(-1d0,0d0)
           exacteigs( 12 )= complex(-5.968916446451269d-1, 6.629807445770295d-1)
           exacteigs( 13 )= complex(-5.968916446451269d-1,-6.629807445770295d-1)
           exacteigs( 14 )= complex(-2.17526747030511d-1,-1.11445426587329)
           exacteigs( 15 )= complex(-2.17526747030511d-1, 1.11445426587329)
           exacteigs( 16 )= complex(-1.6359826155202d-1, 1.09778064288827)
           exacteigs( 17 )= complex(-1.6359826155202d-1,-1.09778064288827)
           exacteigs( 18 )= complex(-1.225611668766536d-1,-7.4486176661974423d-1)
           exacteigs( 19 )= complex(-1.225611668766536d-1, 7.4486176661974423d-1)
           exacteigs( 20 )= complex(-1.13418655949436d-1,-8.605694725015730d-1)
           exacteigs( 21 )= complex(-1.13418655949436d-1,8.605694725015730d-1)
           exacteigs( 22 )= complex(-1.5570386020902d-2, 1.020497366498289d0)
           exacteigs( 23 )= complex(-1.5570386020902d-2,-1.020497366498289d0)
           exacteigs( 24 )= complex(3.59892739012579001d-1, 6.84762020211812856d-1)
           exacteigs( 25 )= complex(3.59892739012579001d-1,-6.84762020211812856d-1)
           exacteigs( 26 )= complex(3.8900684056977123543d-1,-2.1585065087081910777d-1)
           exacteigs( 27 )= complex(3.8900684056977123543d-1, 2.1585065087081910777d-1)
           exacteigs( 28 )= complex(3.96534570032415023d-1, 6.04181810488988837d-1)
           exacteigs( 29 )= complex(3.96534570032415023d-1,-6.04181810488988837d-1)
           exacteigs( 30 )= complex(4.433256333996235387d-1, 3.729624166628465083d-1)
           exacteigs( 31 )= complex(4.433256333996235387d-1,-3.729624166628465083d-1)

      case (14) ! Mandelbort 63 MPSolve
           print*, "Mandelbrot 63 MPSolve"
           RCOEFFS(63) = 1d0
           RCOEFFS(62) = 1d0
           RCOEFFS(61) = 2d0
           RCOEFFS(60) = 5d0
           RCOEFFS(59) = 14d0
           RCOEFFS(58) = 42d0
           RCOEFFS(57) = 132d0
           RCOEFFS(56) = 365d0
           RCOEFFS(55) = 950d0
           RCOEFFS(54) = 2398d0
           RCOEFFS(53) = 5916d0
           RCOEFFS(52) = 14290d0
           RCOEFFS(51) = 33708d0
           RCOEFFS(50) = 77684d0
           RCOEFFS(49) = 175048d0
           RCOEFFS(48) = 385741d0
           RCOEFFS(47) = 831014d0
           RCOEFFS(46) = 1749654d0
           RCOEFFS(45) = 3598964d0
           RCOEFFS(44) = 7228014d0
           RCOEFFS(43) = 14162220d0
           RCOEFFS(42) = 27049196d0
           RCOEFFS(41) = 50323496d0
           RCOEFFS(40) = 91143114d0
           RCOEFFS(39) = 160617860d0
           RCOEFFS(38) = 275276716d0
           RCOEFFS(37) = 458591432d0
           RCOEFFS(36) = 742179284d0
           RCOEFFS(35) = 1166067016d0
           RCOEFFS(34) = 1777171560d0
           RCOEFFS(33) = 2625062128d0
           RCOEFFS(32) = 3754272037d0
           RCOEFFS(31) = 5193067630d0
           RCOEFFS(30) = 6939692682d0
           RCOEFFS(29) = 8948546308d0
           RCOEFFS(28) = 11120136162d0
           RCOEFFS(27) = 13299362332d0
           RCOEFFS(26) = 15286065700d0
           RCOEFFS(25) = 16859410792d0
           RCOEFFS(24) = 17813777994d0
           RCOEFFS(23) = 17999433372d0
           RCOEFFS(22) = 17357937708d0
           RCOEFFS(21) = 15941684776d0
           RCOEFFS(20) = 13910043524d0
           RCOEFFS(19) = 11500901864d0
           RCOEFFS(18) = 8984070856d0
           RCOEFFS(17) = 6609143792d0
           RCOEFFS(16) = 4562339774d0
           RCOEFFS(15) = 2943492972d0
           RCOEFFS(14) = 1766948340d0
           RCOEFFS(13) = 981900168d0
           RCOEFFS(12) = 502196500d0
           RCOEFFS(11) = 234813592d0
           RCOEFFS(10) = 99582920d0
           RCOEFFS(9) = 37945904d0
           RCOEFFS(8) = 12843980d0
           RCOEFFS(7) = 3807704d0
           RCOEFFS(6) = 971272d0
           RCOEFFS(5) = 208336d0
           RCOEFFS(4) = 36440d0
           RCOEFFS(3) = 4976d0
           RCOEFFS(2) = 496d0
           RCOEFFS(1) = 32d0
           eigsknown = 1
           exacteigs( 1 )=complex(-1.999095682327018473210d0,0d0)
           exacteigs( 2 )=complex(-1.9918141725491222157325609498622881d0,0d0)
           exacteigs( 3 )=complex(-1.977179587006257387346088520662828616836d0,0d0)
           exacteigs( 4 )=complex(-1.953705894284396245427622199013653238901d0,0d0)
           exacteigs( 5 )=complex(-1.927147709363950262460068188946594278007d0,0d0)
           exacteigs( 6 )=complex(-1.8848035715866817923294780929158396496359d0,0d0)
           exacteigs( 7 )=complex(-1.8323152027512291920848975260425181432293d0,0d0)
           exacteigs( 8 )=complex(-1.76926167027683114607548022863625740038777d0, 5.6919500395600315304900187298015859319654d-2)
           exacteigs( 9 )=complex(-1.76926167027683114607548022863625740038777d0,-5.6919500395600315304900187298015859319654d-2)
           exacteigs( 10 )=complex(-1.674066091474787971565296029172325596206403d0,0d0)
           exacteigs( 11 )=complex(-1.5748891397523009698199655524959742837719482d0,0d0)
           exacteigs( 12 )=complex(-1.408446485740072654917577008805998851928020904d0, &
                -1.36171997304659915684707793608163610038822995d-1)
           exacteigs( 13 )=complex(-1.408446485740072654917577008805998851928020904d0, &
                1.36171997304659915684707793608163610038822995d-1)
           exacteigs( 14 )=complex(-1.29255806103352208716418470636149411998013630326d0, &
                4.3819881608663183712973712432734844004535476504d-1)
           exacteigs( 15 )=complex(-1.29255806103352208716418470636149411998013630326d0, &
                -4.3819881608663183712973712432734844004535476504d-1)
           exacteigs( 16 )=complex(-1.26228728143847254301011194120806575232050489502d0, &
                4.0810432411269038329016065742601506306041169168d-1)
           exacteigs( 17 )=complex(-1.26228728143847254301011194120806575232050489502d0, &
                -4.0810432411269038329016065742601506306041169168d-1)
           exacteigs( 18 )=complex(-1.25273588401203794629581100256433997387062287256d0, &
                -3.4247064788975089386187578687092843396383393805d-1)
           exacteigs( 19 )=complex(-1.25273588401203794629581100256433997387062287256d0, &
                3.4247064788975089386187578687092843396383393805d-1)
           exacteigs( 20 )=complex(-1.02819385245481759930249745596731843328070508279421d0, &
                -3.61376517118561592479460832997830315786692639704085d-1)
           exacteigs( 21 )=complex(-1.02819385245481759930249745596731843328070508279421d0, &
                3.61376517118561592479460832997830315786692639704085d-1)
           exacteigs( 22 )=complex(-6.23532485956252757990016587001026776428072703359878868d-1, &
                6.81064414225239608090835812686561539088332735217609127d-1)
           exacteigs( 23 )=complex(-6.23532485956252757990016587001026776428072703359878868d-1, &
                -6.81064414225239608090835812686561539088332735217609127d-1)
           exacteigs( 24 )=complex(-6.2243629504129358796016350694723840189750985673649588591d-1, &
                4.2487843647562918431157443880525338683545992964599689876d-1)
           exacteigs( 25 )=complex(-6.2243629504129358796016350694723840189750985673649588591d-1, &
                -4.2487843647562918431157443880525338683545992964599689876d-1)
           exacteigs( 26 )=complex(-5.308278048599427289214772971196026578135170949646890946d-1, &
                6.682887255592057714440924655647011851367651843270734380d-1)
           exacteigs( 27 )=complex(-5.308278048599427289214772971196026578135170949646890946d-1, &
                -6.682887255592057714440924655647011851367651843270734380d-1)
           exacteigs( 28 )=complex(-2.72102461488938894219383324518026874585894699621947085d-1, &
                -8.42364690294128145503155708242929569550778268698265965d-1)
           exacteigs( 29 )=complex(-2.72102461488938894219383324518026874585894699621947085d-1, &
                8.42364690294128145503155708242929569550778268698265965d-1)
           exacteigs( 30 )=complex(-2.24915951286740054685326255204118310792682454680693d-1, &
                1.11626015745499183500126825424467009109873946082435d0)
           exacteigs( 31 )=complex(-2.24915951286740054685326255204118310792682454680693d-1, &
                -1.11626015745499183500126825424467009109873946082435d0)
           exacteigs( 32 )=complex(-2.0728383545566641282413385018667121332401155604017d-1, &
                1.11748077249496291137377567312207879579746389236127d0)
           exacteigs( 33 )=complex(-2.0728383545566641282413385018667121332401155604017d-1, &
                -1.11748077249496291137377567312207879579746389236127d0)
           exacteigs( 34 )=complex(-1.7457822113571696945156643266162905020167505710204d-1, &
                1.07142767145403118922964631021955987671322451961088d0)
           exacteigs( 35 )=complex(-1.7457822113571696945156643266162905020167505710204d-1, &
                -1.07142767145403118922964631021955987671322451961088d0)
           exacteigs( 36 )=complex(-1.57516053475965356164335109644674141293297577896685d-1, &
                -1.10900651411360717797175198615475582901468585712356d0) 
           exacteigs( 37 )=complex(-1.57516053475965356164335109644674141293297577896685d-1, &
                1.10900651411360717797175198615475582901468585712356d0)
           exacteigs( 38 )=complex(-1.274999735463630001995395653459879637298616757217284d-1, &
                9.874609094894567922074076807929788675642068522522938d-1)
           exacteigs( 39 )=complex(-1.274999735463630001995395653459879637298616757217284d-1, &
                -9.874609094894567922074076807929788675642068522522938d-1)
           exacteigs( 40 )=complex(-1.42334819203540667677618453136202688358025283954839d-2, &
                -1.0329147752136441093950134026551104360994260360822540d0)
           exacteigs( 41 )=complex(-1.42334819203540667677618453136202688358025283954839d-2, &
                1.0329147752136441093950134026551104360994260360822540d0)
           exacteigs( 42 )=complex(-6.98356849626139181796649107548406610452886379651341d-3, &
                -1.0036038622882895485307049669513531297649273745391915d0)
           exacteigs( 43 )=complex(-6.98356849626139181796649107548406610452886379651341d-3, &
                1.0036038622882895485307049669513531297649273745391915d0)
           exacteigs( 44 )=complex( 1.4895466603687646529815779208794106185666477731693128d-2, &
                -8.481487619084165277193311117832376290806619901265058603d-1)
           exacteigs( 45 )=complex( 1.4895466603687646529815779208794106185666477731693128d-2, &
                8.481487619084165277193311117832376290806619901265058603d-1)
           exacteigs( 46 )=complex( 1.211927861059064863147044434105037593859287800520963579338d-1, &
                6.1061169221075421167538724415035774824319702690063863369691d-1)
           exacteigs( 47 )=complex( 1.211927861059064863147044434105037593859287800520963579338d-1, &
                -6.1061169221075421167538724415035774824319702690063863369691d-1)
           exacteigs( 48 )=complex( 3.52482539722363278193253964052161589243593334212239870706d-1, &
                -6.98337239583330331258141954760484537633150485928512286760d-1)
           exacteigs( 49 )=complex( 3.52482539722363278193253964052161589243593334212239870706d-1, &
                6.98337239583330331258141954760484537633150485928512286760d-1)
           exacteigs( 50 )=complex( 3.7600868184676755970480431772902286888800357334637481632029182d-1, &
                -1.4474937132163286474711018201298830556966056842762643026975894d-1)
           exacteigs( 51 )=complex( 3.7600868184676755970480431772902286888800357334637481632029182d-1, &
                1.4474937132163286474711018201298830556966056842762643026975894d-1)
           exacteigs( 52 )=complex( 3.76893240379311323690004017968512473363482317941533875341d-1, &
                6.78568693190448141957540792996773280196881194582788907016d-1)
           exacteigs( 53 )=complex( 3.76893240379311323690004017968512473363482317941533875341d-1, &
                -6.78568693190448141957540792996773280196881194582788907016d-1)
           exacteigs( 54 )=complex( 3.865391765961580265082930869043677799284877313516569138807d-1, &
                5.693247113031029032137923571351905081619323911951388853856d-1)
           exacteigs( 55 )=complex( 3.865391765961580265082930869043677799284877313516569138807d-1, &
                -5.693247113031029032137923571351905081619323911951388853856d-1)
           exacteigs( 56 )=complex( 4.12916024722700479197334566382612257174765142865547121703d-1, &
                6.148067601433856949545497204007997358291659758563137777616d-1)
           exacteigs( 57 )=complex( 4.12916024722700479197334566382612257174765142865547121703d-1, &
                -6.148067601433856949545497204007997358291659758563137777616d-1)
           exacteigs( 58 )=complex( 4.3237619264199450782466964808692137388785063987699403620424125d-1, &
                2.267599044353486186978765599716989721202321914603899690444951d-1)
           exacteigs( 59 )=complex( 4.3237619264199450782466964808692137388785063987699403620424125d-1, &
                -2.267599044353486186978765599716989721202321914603899690444951d-1)
           exacteigs( 60 )=complex( 4.52774498724915493508803077732546131473562000961307327749350d-1, &
                -3.96170128033165002412596877271155937712569079351815707744770d-1)
           exacteigs( 61 )=complex( 4.52774498724915493508803077732546131473562000961307327749350d-1, &
                3.96170128033165002412596877271155937712569079351815707744770d-1)
           exacteigs( 62 )=complex( 4.56823285823316651283953236253270107801699459631329688710054d-1, &
                3.47758700883481983632188723200264206004781117755664551397643d-1)
           exacteigs( 63 )=complex( 4.56823285823316651283953236253270107801699459631329688710054d-1, &
                -3.47758700883481983632188723200264206004781117755664551397643d-1)

        case (15,16)  ! Vanni Noferini's example
           print*, "Vanni Noferini's example degree 12 or 35"
           do ii=1,n
              call cnormalpoly(1,normofp)
              exacteigs(ii) = complex(dble(normofp),0d0)
           end do
           exacteigs(1) = exacteigs(1)*1d9
           exacteigs(2) = exacteigs(2)*1d12
           eigsknown = 1
           nobalance = 1
           do ii=1,n
              xr_db(ii) = dble(exacteigs(ii))
              xi_db(ii) = 0d0
           end do
           call rootstocoeffs(n,xr_db,xi_db,c_db)
           RCOEFFS = c_db


        case (17)  
           print*, "cubic rcoeffsnomial small roots"
           exacteigs(1) = complex(1d0,0d0)
           exacteigs(2) = complex(1d-8,0d0)
           exacteigs(3) = complex(-1d-8,0d0)
           eigsknown = 1
           nobalance = 0
           do ii=1,n
              xr_db(ii) = dble(exacteigs(ii))
              xi_db(ii) = imag(exacteigs(ii))
           end do
           call rootstocoeffs(n,xr_db,xi_db,c_db)
           RCOEFFS = c_db

        case (18)  
           print*, "cubic rcoeffsnomial very small roots"
           exacteigs(1) = complex(1d0,0d0)
           exacteigs(2) = complex(1d-15,0d0)
           exacteigs(3) = complex(-1d-15,0d0)
           
           eigsknown = 1
           nobalance = 0
           do ii=1,n
              xr_db(ii) = dble(exacteigs(ii))
              xi_db(ii) = imag(exacteigs(ii))
           end do
           call rootstocoeffs(n,xr_db,xi_db,c_db)
           RCOEFFS = c_db

        case (19)  
           print*, "cubic rcoeffsnomial large roots"
           exacteigs(1) = complex(1d0,0d0)
           exacteigs(2) = complex(1d+8,0d0)
           exacteigs(3) = complex(-1d+8,0d0)
           eigsknown = 1
           nobalance = 0
           do ii=1,n
              xr_db(ii) = dble(exacteigs(ii))
              xi_db(ii) = imag(exacteigs(ii))
           end do
           call rootstocoeffs(n,xr_db,xi_db,c_db)
           RCOEFFS = c_db

        case (20)  
           print*, "cubic rcoeffsnomial very large roots"
           exacteigs(1) = complex(1d0,0d0)
           exacteigs(2) = complex(1d+15,0d0)
           exacteigs(3) = complex(-1d+15,0d0)
           eigsknown = 1
           nobalance = 0
           do ii=1,n
              xr_db(ii) = dble(exacteigs(ii))
              xi_db(ii) = imag(exacteigs(ii))
           end do
           call rootstocoeffs(n,xr_db,xi_db,c_db)
           RCOEFFS = c_db

        case (21,22)  
           print*, "underflow test"
           do ii=1,n
              exacteigs(ii) = complex(10d0**(-ii),0d0)
           end do
           eigsknown = 1
           nobalance = 0
           do ii=1,n
              xr_db(ii) = dble(exacteigs(ii))
              xi_db(ii) = imag(exacteigs(ii))
           end do
           call rootstocoeffs(n,xr_db,xi_db,c_db)
           RCOEFFS = c_db


        case (23,24,25)  
           print*, "deflation stability test"
           if (ll == 23) then
              hA = 10d3
           else if (ll == 24) then
              hA = 10d6
           else if (ll == 25) then
              hA = 10d9
           else
              print*, "ERROR"
           end if
           exacteigs(1) = complex(1,0d0)
           exacteigs(2) = complex(hA,0d0)
           exacteigs(3) = complex(1d0/hA,0d0)
           eigsknown = 1
           nobalance = 0
           do ii=1,n
              xr_db(ii) = dble(exacteigs(ii))
              xi_db(ii) = imag(exacteigs(ii))
           end do
           call rootstocoeffs(n,xr_db,xi_db,c_db)
           RCOEFFS = c_db

        case (26)  
           print*, "deflation stability test M = 15"
           do ii=-14,0
              exacteigs(15+ii) = 0.9d0*exp(complex(0d0,ii*pi/2d0/15d0))
           end do
           do ii=1,15
              exacteigs(15+ii) = conjg(exacteigs(16-ii))
           end do
           do ii=16,30
              exacteigs(15+ii) = exp(complex(0d0,ii*pi/2d0/15d0))
           end do
           do ii=31,45
              exacteigs(15+ii) = conjg(exacteigs(76-ii))
           end do
           eigsknown = 1
           nobalance = 0
           do ii=1,n
              xr_db(ii) = dble(exacteigs(ii))
              xi_db(ii) = imag(exacteigs(ii))
           end do
           call rootstocoeffs(n,xr_db,xi_db,c_db)
           RCOEFFS = c_db

        case (27)  ! Bernoulli rcoeffsnomial of degree 20  
           print*, "Bernoulli rcoeffsnomial of degree 20"
           RCOEFFS(20) = -174611d0/330d0    ! 1
           RCOEFFS(19) = 0d0                ! z^1  = 0d0
           RCOEFFS(18) = +219335d0/21d0     ! z^2
           RCOEFFS(17) = 0d0                ! z^3  = 0d0
           RCOEFFS(16) = -68723d0/2d0       ! z^4
           RCOEFFS(15) = 0d0                ! z^5  = 0d0
           RCOEFFS(14) = +45220d0           ! z^6
           RCOEFFS(13) = 0d0                ! z^7  = 0d0   
           RCOEFFS(12) = -223193d0/7d0      ! z^8
           RCOEFFS(11)= 0d0                 ! z^9
           RCOEFFS(10)= +41990d0/3d0        ! z^10
           RCOEFFS( 9)= 0d0                 ! z^11
           RCOEFFS( 8)= -4199d0             ! z^12
           RCOEFFS( 7)= 0d0                 ! z^13
           RCOEFFS( 6)= +6460d0/7d0         ! z^14
           RCOEFFS( 5)= 0d0                 ! z^15
           RCOEFFS( 4)= -323d0/2d0          ! z^16
           RCOEFFS( 3)= 0d0                 ! z^17
           RCOEFFS( 2)= +95d0/3d0           ! z^18
           RCOEFFS( 1)= -10d0               ! z^19
           eigsknown = 0

        case (28)   ! p(z) = (20!) sum_{k=0}^{20} z^k/k!
           print*, "p(z) = (20!) sum_{k=0}^{20} z^k/k!"
           RCOEFFS(1) = 0.2d+2
           RCOEFFS(2) = 0.38d+3
           RCOEFFS(3) = 0.684d+4
           RCOEFFS(4) = 116280d+0
           RCOEFFS(5) = 1860480d+0
           RCOEFFS(6) = 27907200d+0
           RCOEFFS(7) = 390700800d+0
           RCOEFFS(8) = 5.0791104d+9
           RCOEFFS(9) = 6.09493248d+10
           RCOEFFS(10)= 6.704425728000000d+11
           RCOEFFS(11)= 6.704425728000000d+12
           RCOEFFS(12)= 6.033983155200000d+13
           RCOEFFS(13)= 4.827186524160000d+14
           RCOEFFS(14)= 3.379030566912000e+15
           RCOEFFS(15)= 2.027418340147200e+16
           RCOEFFS(16)= 1.013709170073600e+17
           RCOEFFS(17)= 4.054836680294400e+17
           RCOEFFS(18)= 1.216451004088320e+18
           RCOEFFS(19)= 2.432902008176640e+18
           RCOEFFS(20)= 2.432902008176640e+18
           eigsknown = 0

        case (29,30,31,32,33)  
           print*, "Bevilacqua, Del Corso, Gemignani P1"
           eigsknown = 0
           nobalance = 0
           RCOEFFS = 0d0
           RCOEFFS(N) = 1 
           RCOEFFS(N/2) = (1d0*N)/(1d0*N+1d0) + (1d0*N+1d0)/(1d0*N)

        case (34,35,36,37,38)  
           print*, "Bevilacqua, Del Corso, Gemignani P2"
           eigsknown = 0
           nobalance = 0
           do ii = 1,n/2-1
              RCOEFFS(ii) = 1d0/n*(n+ii)
              RCOEFFS(n-ii) = 1d0/n*(n+ii)
           end do
           RCOEFFS(n)= 1d0
           RCOEFFS(n/2) = (n+1d0)/(1d0*n)

        case (39,40,41,42,43,44,45,46,47,48)
           print*, "Bevilacqua, Del Corso, Gemignani P3"
           eigsknown = 0
           nobalance = 0
           if (ll<=43) then
              lambda = 0.9
           else
              lambda = 0.999
           end if
           RCOEFFS = 0d0
           RCOEFFS(N) = -1d0
           RCOEFFS(N-1) = (lambda+1d0)/(1d0-lambda)
           RCOEFFS(1) = -(lambda+1d0)/(1d0-lambda)

        end select


        if (nobalanceor > 0) then
           if (nobalanceor == 1) then
              nobalance = 0
           else if (nobalanceor == 2) then
              nobalance = 1
           end if
           print*, "nobalance", nobalance
        end if
        if (nobalance == 0) then
           call balance(n,rcoeffs,icoeffs,nnew,rnew,inew,alpha)
        else
           alpha = 1d0
           nnew=n
           do ii=1,N
              rnew(ii)=rcoeffs(ii)
              inew(ii)=icoeffs(ii)
           end do
        end if
        print*, "alpha", alpha

        do ii=1,N
           poly(ii) = complex(rcoeffs(ii),icoeffs(ii))       
           npoly(ii) = complex(rnew(ii),inew(ii))       
        end do
     

        error = 0d0
        
    
        normofp = dsqrt(dznrm2(N,poly,1)**2 + 1d0)


     ! ZAMVW    
     write(*,*) "ZAMVW"
     CALL SYSTEM_CLOCK(COUNT_RATE=clock_rate) 
     CALL SYSTEM_CLOCK(COUNT=clock_start)
         
     do ii=1,num_trials(ll)
          
        call zamvw(n,rnew(1:n),inew(1:n),reigs,ieigs,its,ry)
        print*, ry
        do jj=1,n
           ieigs(jj) = alpha*ieigs(jj)
           reigs(jj) = alpha*reigs(jj)
           roots(jj) = complex(reigs(jj),ieigs(jj))
           if (num_trials(ll)==1) then
              print*, roots(jj)
           end if
	end do
 
        call RESCHECK(0,N,1,1,POLY(1:N),COEFFS,ROOTS,ALLROOTS,RES)
        
        do jj=1,6
           do kk=1,N
              if(res(kk,jj) > error(jj))then
                 error(jj) = res(kk,jj)
              end if
           end do
        end do
        
        error(7) = 0
        error(8) = 0
        
     end do
     CALL SYSTEM_CLOCK(COUNT=clock_end)  
     call backward(n,poly,reigs,ieigs,poly2,err,relerr)
     time = dble(clock_end - clock_start)/dble(clock_rate) 
     time = time/dble(num_trials(ll))

     if (eigsknown == 1) then
        ! sum of abs error, take each eigenvalue only ones
        do kk=1,N
           eigtaken(kk) = 0
        end do
        kl = 2
        error(7) = 0d0
        do kj=1,N
           jj = kl-1
           do kk = kl,N
              if ( (abs(roots(kj)-exacteigs(kk))<abs(roots(kj)-exacteigs(jj))) .AND. (eigtaken(kk)==0) ) then
                 jj = kk
              end if
           end do
           eigtaken(jj)=1
           if (num_trials(ll)==1) then
              print*, kj, jj, roots(kj), exacteigs(jj), abs(roots(kj)- exacteigs(jj))
           end if
           error(7)=error(7)+abs(roots(kj)-exacteigs(jj))
           if (jj==kl-1) then
              kl = kl + 1
           end if
        end do
           
        print*, "sum of abs err", error(7)
        error(8) = 0d0
        error(9) = 0d0
        do kj=1,N
           hd = abs(roots(kj)-exacteigs(1))
           klm = 1
           do jj = 2,N
              if (hd>abs(roots(kj)-exacteigs(jj))) then
                 hd = abs(roots(kj)-exacteigs(jj))
                 klm = jj
              end if
           end do
           if (error(8)<hd) then
              error(8) = hd             
           end if
           if (error(9)<hd/abs(exacteigs(klm))) then
              error(9) = hd/abs(exacteigs(klm))             
           end if

        end do
        print*, "max abs err   ", error(8)
        print*, "max rel err   ", error(9)
     end if
              
     print*, error(:)

     backwarderror(1)  = error(1)
     relforwarderror(1)= error(9)
     
     jj=idamax(N,relerr,1)
     back2(1) = relerr(jj)
     jj=idamax(N,err,1)
     back3(1) = err(jj)



     ! LAPACK
     write(*,*) "LAPACK"
     error = 0d0
     res = 0d0
     
     CALL SYSTEM_CLOCK(COUNT_RATE=clock_rate) 
     CALL SYSTEM_CLOCK(COUNT=clock_start)
     
     if (N<4200) then
        lapack_trials=num_trials(ll)
     else
        lapack_trials=1
        if (N>=2000) then
           lapack_trials=0
        end if
     end if
     do ii=1,lapack_trials
        
        B = complex(0.d0,0.d0)
        B(1,:) = -npoly(1:N)
        do kk=1,(N-1)
           B(kk+1,kk) = complex(1.d0,0.d0)
        end do
        
        call ZHSEQR('E','N',N,1,N,B,N,roots,Y,N,cwork,N,info)
        print*, "info ZHSEQR", info
        do jj=1,n
           roots(jj) = alpha*roots(jj)
           reigs(jj) = dble(roots(jj))
           ieigs(jj) = imag(roots(jj))
           if (num_trials(ll)==1) then
              print*, roots(jj)
           end if
	end do
        call RESCHECK(0,N,1,1,POLY(1:N),COEFFS,ROOTS,ALLROOTS,RES)
        
        do jj=1,6
           do kk=1,N
              if(res(kk,jj) > error(jj))then
                 error(jj) = res(kk,jj)
              end if
           end do
        end do
     end do

     CALL SYSTEM_CLOCK(COUNT=clock_end)  
     call backward(n,poly,reigs,ieigs,poly2,err,relerr)
     time = dble(clock_end - clock_start)/dble(clock_rate) 
     time = time/dble(lapack_trials)

     error(7) = 0
     error(8) = 0
     
     if (eigsknown == 1) then
        do kk=1,N
           eigtaken(kk) = 0
        end do
        kl = 2
        error(7)=0
        do kj=1,N
           jj = kl-1
           do kk = kl,N
              if ( (abs(roots(kj)-exacteigs(kk))<abs(roots(kj)-exacteigs(jj))) .AND. (eigtaken(kk)==0) ) then
                 jj = kk
              end if
           end do
           eigtaken(jj)=1 
           if (num_trials(ll)==1) then
              print*, kj, jj, roots(kj), exacteigs(jj), abs(roots(kj)- exacteigs(jj))
           end if
           error(7)=error(7)+abs(roots(kj)-exacteigs(jj))
           if (jj==kl-1) then
              kl = kl + 1
           end if
        end do
        
        print*, "sum of abs err", error(7)
        error(8) = 0d0
        error(9) = 0d0
        do kj=1,N
           hd = abs(roots(kj)-exacteigs(1))
           klm = 1
           do jj = 2,N
              if (hd>abs(roots(kj)-exacteigs(jj))) then
                 hd = abs(roots(kj)-exacteigs(jj))
                 klm = jj
              end if
           end do
           if (error(8)<hd) then
              error(8) = hd             
           end if
           if (error(9)<hd/abs(exacteigs(klm))) then
              error(9) = hd/abs(exacteigs(klm))             
           end if

        end do
        print*, "max abs err   ", error(8)
        print*, "max rel err   ", error(9)
     end if
     
     
     print*, error(:)

     backwarderror(2)  = error(1)
     relforwarderror(2)= error(9)
     jj=idamax(N,relerr,1)
     back2(2) = relerr(jj)
     jj=idamax(N,err,1)
     back3(2) = err(jj)


     ! LAPACK
     write(*,*) "LAPACK"
     error = 0d0
     res = 0d0
     
     CALL SYSTEM_CLOCK(COUNT_RATE=clock_rate) 
     CALL SYSTEM_CLOCK(COUNT=clock_start)
     
     if (N<4200) then
        lapack_trials=num_trials(ll)
     else
        lapack_trials=1
        if (N>=2000) then
           lapack_trials=0
        end if
     end if
     do ii=1,lapack_trials
        
        B = complex(0.d0,0.d0)
        B(1,:) = -npoly(1:N)
        do kk=1,(N-1)
           B(kk+1,kk) = complex(1.d0,0.d0)
        end do
        
        call ZGEEV('N','N',N,B,N,roots,Y,N,Y,N,cwork,5*N,rwork,info)

        print*, "info ZGEEV", info
        do jj=1,n
           roots(jj) = alpha*roots(jj)
           reigs(jj) = dble(roots(jj))
           ieigs(jj) = imag(roots(jj))
           if (num_trials(ll)==1) then
              print*, roots(jj)
           end if
	end do
        call RESCHECK(0,N,1,1,POLY(1:N),COEFFS,ROOTS,ALLROOTS,RES)
        
        do jj=1,6
           do kk=1,N
              if(res(kk,jj) > error(jj))then
                 error(jj) = res(kk,jj)
              end if
           end do
        end do
     end do

     CALL SYSTEM_CLOCK(COUNT=clock_end)  
     call backward(n,poly,reigs,ieigs,poly2,err,relerr)
     time = dble(clock_end - clock_start)/dble(clock_rate) 
     time = time/dble(lapack_trials)

     error(7) = 0
     error(8) = 0
     
     if (eigsknown == 1) then
        do kk=1,N
           eigtaken(kk) = 0
        end do
        kl = 2
        error(7)=0
        do kj=1,N
           jj = kl-1
           do kk = kl,N
              if ( (abs(roots(kj)-exacteigs(kk))<abs(roots(kj)-exacteigs(jj))) .AND. (eigtaken(kk)==0) ) then
                 jj = kk
              end if
           end do
           eigtaken(jj)=1 
           if (num_trials(ll)==1) then
              print*, kj, jj, roots(kj), exacteigs(jj), abs(roots(kj)- exacteigs(jj))
           end if
           error(7)=error(7)+abs(roots(kj)-exacteigs(jj))
           if (jj==kl-1) then
              kl = kl + 1
           end if
        end do
        
        print*, "sum of abs err", error(7)
        error(8) = 0d0
        error(9) = 0d0
        do kj=1,N
           hd = abs(roots(kj)-exacteigs(1))
           klm = 1
           do jj = 2,N
              if (hd>abs(roots(kj)-exacteigs(jj))) then
                 hd = abs(roots(kj)-exacteigs(jj))
                 klm = jj
              end if
           end do
           if (error(8)<hd) then
              error(8) = hd             
           end if
           if (error(9)<hd/abs(exacteigs(klm))) then
              error(9) = hd/abs(exacteigs(klm))             
           end if

        end do
        print*, "max abs err   ", error(8)
        print*, "max rel err   ", error(9)
     end if
     
     
     print*, error(:)

     backwarderror(5)  = error(1)
     relforwarderror(5)= error(9)
     jj=idamax(N,relerr,1)
     back2(5) = relerr(jj)
     jj=idamax(N,err,1)
     back3(5) = err(jj)


     print*, "Relative Forward Error"
     write(*,"(I2,1x,A,1x,ES10.4E2,1x,A,1x,ES10.4E2,1x,A,1x,&
          &ES10.4E2)"), ll, "&", &
          & relforwarderror(1), "&", relforwarderror(2),&
          & "&", relforwarderror(5) 

     print*, "Backward Error (Abs Error in Coefficients divided by norm of coefficients) "
     write(*,"(I2,1x,A,1x,ES10.4E2,1x,A,1x,ES10.4E2,1x,A,1x,&
          &ES10.4E2)"), ll, "&", &
          & back3(1)/normofp, "&", back3(2)/normofp,&
          & "&", back3(5)/normofp

     print*, ""
     print*, ""

     print*, desc(ll)

     write(*,"(I2,1x,A,1x,ES10.4E2,1x,A,1x,ES10.4E2,1x,A,1x,&
          &ES10.4E2,1x,A)"), ll, "&", &
          & relforwarderror(1), "&", relforwarderror(2), "&", relforwarderror(5),"\\%"
     write(31,"(I2,1x,A,1x,ES10.4E2,1x,A,1x,ES10.4E2,1x,A,1x&
          &ES10.4E2,1x,A)"), ll, "&", &
          & relforwarderror(1), "&", relforwarderror(2), "&", relforwarderror(5),"\\%"
     write(*,"(I2,1x,A,1x,ES10.4E2,1x,A,1x,ES10.4E2,1x,A,1x,&
          &ES10.4E2,1x,A,1x,ES10.4E2)"), ll, "&", &
          & back3(1)/normofp, "&", back3(2)/normofp, "&", back3(5)/normofp,"\\%",normofp
     write(32,"(I2,1x,A,1x,ES10.4E2,1x,A,1x,ES10.4E2,1x,A,1x&
          &ES10.4E2,1x,A,1x,ES10.4E2)"), ll, "&", &
          & back3(1)/normofp, "&", back3(2)/normofp, "&", back3(5)/normofp,"\\%",normofp

     close(31)
     close(32)
     open (unit=31, file="sp_table1.txt", status='unknown', position="append")
     open (unit=32, file="sp_table2.txt", status='unknown', position="append")

     ! free memory
     deallocate(Q,D,C,B2,rcoeffs,icoeffs)
     deallocate(its,reigs,ieigs,rnew,inew,npoly)
     deallocate(residuals)
     deallocate(poly,roots,allroots,iterations,res)
     deallocate(B,Y,cwork,rwork)
     deallocate(exacteigs, eigtaken)
     deallocate(poly2,err,relerr)
     deallocate(xr_db,xi_db,c_db)
  end if
     
  end do
  
  close(31)
  close(32)

end program main
