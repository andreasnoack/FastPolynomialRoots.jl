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
! Test program for random polynomials with complex
! coefficients of given size
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
program blocktest

  implicit none
  
  ! compute variables
  integer :: n, nnew, ii, strt, newtnum
  integer, allocatable :: its(:), seed(:)
  double precision, allocatable :: Q(:),D(:),C(:),B(:)
  double precision, allocatable :: rcoeffs(:),icoeffs(:)
  double precision, allocatable :: rnew(:),inew(:)
  double precision, allocatable :: reigs(:),ieigs(:)
  complex(kind(1d0)), allocatable :: poly(:),eigs(:),allroots(:,:)
  double precision, allocatable ::residuals(:,:)
  complex(kind(1d0)) :: block(2,2),coeffs
  double precision :: time, error1, error2, alpha
  integer :: clock_start, clock_end, clock_rate, rsize, flag
  character(len=32) :: arg
  

  ! initial random seed
  ! call init_random_seed()
  
  if (iargc()>0) then
     call RANDOM_SEED(size = rsize)
     allocate(seed(rsize))
     call RANDOM_SEED(GET = seed)
     
     print*, iargc()
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
        print*, "random seed"
        call init_random_seed()
     end if
     
     call getarg(1, arg)
     read (arg,'(I10)') n
  else
     ! set degree
     n = 2**12
  end if
  
		
	! set newtnum
	newtnum = 1
	
	! allocate memory
	allocate(Q(3*n),D(2*n+2),C(3*n),B(3*n),rcoeffs(n),icoeffs(n))
	allocate(rnew(n),inew(n))
	allocate(its(n),reigs(n),ieigs(n))
	allocate(poly(n),eigs(n),allroots(n,newtnum+1),residuals(n,3*(newtnum+1)))
	
	! initialize arrays
	Q = 0d0
	D = 0d0
	C = 0d0
	B = 0d0
	rcoeffs = 0d0
	icoeffs = 0d0
	rnew = 0d0
	inew = 0d0
	its = 0
	reigs = 0d0
	ieigs = 0d0
	poly = complex(0d0,0d0)
	eigs = complex(0d0,0d0)
	allroots = complex(0d0,0d0)
	residuals = 0d0
	
	! build random poly
        !call init_random_seed()
	call normalpoly(n,rcoeffs,icoeffs)
!	rcoeffs = 0d0
!	icoeffs = 0d0
	!icoeffs(n) = -1d0!/sqrt(2d0)
!	rcoeffs(n) = -1d0!/sqrt(2d0)

	! print poly
!	print*,"poly"
!	do ii=1,n
        ii = 1
	print*,rcoeffs(ii),icoeffs(ii)
!	end do
!	print*,""

	! balance
	call balance(n,rcoeffs,icoeffs,nnew,rnew,inew, alpha)
	print*,"alpha =",alpha
	print*,""
	
	! print poly
!	print*,"poly"
!	do ii=1,nnew
!		print*,rnew(ii),inew(ii)
!	end do
!	print*,""
	
	! check new degree
	if(n /= nnew)then
		print*,"poly has trivial zeros"
	end if

	! store in complex array
	do ii=1,nnew
		poly(ii) = complex(rcoeffs(ii),icoeffs(ii))
	end do
	
	! print random poly
!	print*,"coeffs"
!	do ii=1,n
!		print*,rcoeffs(ii),icoeffs(ii)
!	end do
!	print*,""
	
	! start timer
	CALL SYSTEM_CLOCK(COUNT_RATE=clock_rate) 
	CALL SYSTEM_CLOCK(COUNT=clock_start)
	! factor companion matrix
	call factor(nnew,rnew,inew,Q,D,C,B)
	
!	print*,"Q"
!	do ii=1,n
!		strt = 3*(ii-1)
!		print*,Q(strt+1),Q(strt+2)
!		print*,Q(strt+3),0d0
!		print*,""
!	end do
!	print*,""
	
!	print*,"D"
!	do ii=1,n+1
!		strt = 2*(ii-1)
!		print*,D(strt+1),D(strt+2)
!		print*,""
!	end do
!	print*,""
	
!	print*,"C"
!	do ii=1,n
!		strt = 3*(ii-1)
!		print*,C(strt+1),C(strt+2)
!		print*,C(strt+3),0d0
!		print*,""
!	end do
!	print*,""
	
!	print*,"B"
!	do ii=1,n
!		strt = 3*(ii-1)
!		print*,B(strt+1),B(strt+2)
!		print*,B(strt+3),0d0
!		print*,""
!	end do
!	print*,""

  ! store in complex array
  do ii=1,nnew
     poly(ii) = complex(rcoeffs(ii),icoeffs(ii))
  end do
	
  ! start timer
  CALL SYSTEM_CLOCK(COUNT_RATE=clock_rate) 
  CALL SYSTEM_CLOCK(COUNT=clock_start)
  ! factor companion matrix
  call factor(nnew,rnew,inew,Q,D,C,B)
	
	! compute eigs
	call zamvw2(nnew,Q,D,C,B,reigs,ieigs,its,flag,nnew-1)
	
	! store eigs in complex array
	do ii=1,n
		eigs(ii) = alpha*complex(reigs(ii),ieigs(ii))
	end do
	
	! compute residuals
	call rescheck(0,nnew,0,1,poly,coeffs,eigs,allroots,residuals)
	
	! stop timer
	CALL SYSTEM_CLOCK(COUNT=clock_end)  
	time = dble(clock_end - clock_start)/dble(clock_rate)
	print*, "n =",n
	print*,'Total time =', time, 'secs'
	
	! compute worst error
	error1 = 0d0
	do ii=1,nnew
		if(residuals(ii,1) > error1)then
			error1 = residuals(ii,1)
		end if
	end do
	
	! compute worst error
	error2 = 0d0
	do ii=1,nnew
		if(residuals(ii,4) > error2)then
			error2 = residuals(ii,4)
		end if
	end do
	
	! print worst error
	print*,"worst error =",error1,error2
	print*,""
	
	! print residuals
!	print*,"residuals,reigs,ieigs"
!	do ii=1,n
!		print*,residuals(ii,1:3),reigs(ii),ieigs(ii)
!	end do
!	print*,""

	! print output
!	print*,"reigs,ieigs,its"
!	do ii=1,n
!		print*,reigs(ii),ieigs(ii),its(ii)
!	end do
!	print*,""

	! free memory
	deallocate(Q,D,C,B,rcoeffs,icoeffs)
	deallocate(its,reigs,ieigs)
	deallocate(poly,eigs,allroots,residuals)
	deallocate(rnew,inew)
	
end program
