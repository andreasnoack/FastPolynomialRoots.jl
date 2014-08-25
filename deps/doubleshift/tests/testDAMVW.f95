!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Aurentz² Mach³ Vandebril³ Watkins²
!
! ²Dept. Mathematics, Washington State University
! ³Dept. Computer Science, KU Leuven
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Last modified 22 August 2014
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! Test program for random polynomials of given 
! size
!
! The roots are computed by Francis's implicitly
! shifted QR algorithm via the Companion matrix.
! The rank-structure in the iterates is used.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! input parameter
!
! 1) problem size, default 4096
!
! 2) seed for random number generator, 
!    default random seed based on CPU clock, 
!    set fixed seed for reproducibility
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program testDAMVW

  implicit none
  
  ! input variables
  integer :: N, FLAG, NEWTNUM, rsize
  complex(kind(1d0)), allocatable :: COEFFS(:), ALLROOTS(:,:), ROOTS(:), WPOLY(:)
  double precision, allocatable :: POLY(:), REIGS(:), IEIGS(:), RESIDUALS(:,:)
  integer, allocatable :: ITS(:),  seed(:)
  
  ! compute variables
  integer :: ii, noits, mri, mri1, mri2, mri3, kk
  integer :: clock_start, clock_end, clock_rate 
  real :: time
  double precision :: rpart, ipart, temp, mr, mr1, mr2, mr3
  character(len=32) :: arg
  
  FLAG = 1
  if (iargc()>0) then
     if (iargc()>2) then
        FLAG=101
     end if
     call RANDOM_SEED(size = rsize)
     allocate(seed(rsize))
     call RANDOM_SEED(GET = seed)

     !print*, iargc()
     !print*, arg
     if (iargc()>1) then
        call getarg(2, arg)
        print*, arg

        read (arg,'(I10)') ii
        !print*, seed
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
     else
        call init_random_seed()
     end if

     call getarg(1, arg)
     read (arg,'(I10)') kk
  else
     kk = 4096
  end if

  N = kk
  print*, N

  NEWTNUM = 1
  
  ! open (unit=7, file='poly.txt', status='unknown')
  ! read(7,*) N
  
  allocate(POLY(N),REIGS(N),IEIGS(N),ITS(N),COEFFS(1),ALLROOTS(N,NEWTNUM+1),RESIDUALS(N,3*(NEWTNUM+1)))
  allocate(WPOLY(N),ROOTS(N))
  
  
  call DNORMALPOLY(N,POLY)

  ! roots of x^n - 1 = 0
  ! POLY=0d0
  ! POLY(N)=-1d0
  
  if (iargc()>1) then
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
  end if
  
  !open (unit=7, file='poly.txt', status='unknown')
  !write(7,*) N
  !do ii=1,N
  !   write(7,*) poly(ii)
  !end do
  !close(7)
 
  ! start timer
  CALL SYSTEM_CLOCK(COUNT_RATE=clock_rate) 
  CALL SYSTEM_CLOCK(COUNT=clock_start)
  ! compute roots
  call DAMVW(N,POLY,REIGS,IEIGS,ITS,FLAG)
  
  do ii=1,N
     ROOTS(ii) = complex(REIGS(ii),IEIGS(ii))
     WPOLY(ii) = complex(poly(ii),0d0)
  end do
  ! check roots
  call RESCHECK(0,N,0,NEWTNUM,WPOLY,COEFFS,ROOTS,ALLROOTS,RESIDUALS)
  ! stop timer
  CALL SYSTEM_CLOCK(COUNT=clock_end)  
  time = real(clock_end - clock_start)/real(clock_rate)
  
  print*, "Residuals"
  temp = 0
  mr = 0.d0
  mr1 = 0.d0
  mr2 = 0.d0
  mr3 = 0.d0
  mri = 0
  mri1 = 0
  mri2 = 0
  mri3 = 0
  do ii=1,N
     temp = temp + abs(RESIDUALS(ii,1))**2
     if (RESIDUALS(ii,1)>mr) then
        mr3 = mr2
        mr2 = mr1
        mr1 = mr
        mr = RESIDUALS(ii,1)
        mri3 = mri2
        mri2 = mri1
        mri1 = mri
        mri = ii
     else if (RESIDUALS(ii,1)>mr1) then
        mr3 = mr2
        mr2 = mr1
        mr1 = RESIDUALS(ii,1)
        mri3 = mri2
        mri2 = mri1
        mri1 = ii
     else if (RESIDUALS(ii,1)>mr2) then
        mr3 = mr2
        mr2 = RESIDUALS(ii,1)
        mri3 = mri2
        mri2 = ii
     else if (RESIDUALS(ii,1)>mr3) then
        mr3 = RESIDUALS(ii,1)
        mri3 = ii
     end if
  end do
  print*, ""
  print*, "||Residual||_2 = ", dsqrt(temp)
  print*, "max(Residual)  = ", mr, mr1, mr2, mr3
  print*, "imax(Residual) = ", mri, mri1, mri2, mri3
  
  print*, "its",  its(n),its(n-1),its(n-2),its(n-3),its(n-4),its(n-5),its(n-6)
  print*, "its",  its(n-7),its(n-8),its(n-9),its(n-10),its(n-11),its(n-12),its(n-13)
  
  print*, "#IT   = ", noits
  print*, "#IT/N = ", real(noits)/real(N)
  
  print*, "N =",N
  print*,'Total time =', time, 'secs'
  

  if (N<=20) then
     do ii=1,N
        print*, REIGS(ii), IEIGS(ii)
     end do
  end if

  deallocate(POLY,REIGS,IEIGS,ITS,COEFFS,ALLROOTS,RESIDUALS,WPOLY,ROOTS);

end program
