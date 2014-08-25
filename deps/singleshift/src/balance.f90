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
! balance the coefficients of a complex polynomial
! simple scaling by p*(z) = p(alpha z)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                !  
!  !!   This is a simple balancing routine.   !! !
!  !!   Numerical experiments show that it    !! !
!  !!   is not always advantageous. Further   !! !
!       research is necessary.                   ! 
!  !!   Use with caution.                     !! !         
!                                                !  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!
! n       problem size
! 
! rcoeffs real parts coefficients
! icoeffs imag parts coefficients
!
! nnew    deflated problem size
!
! rnew    real parts balanced coefficients
! inew    imag parts balanced coefficients
!
! alpha   scaling parameter 
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine balance(n,rcoeffs,icoeffs,nnew,rnew,inew,alpha)

	implicit none
	
	! input variables
	integer, intent(in) :: n
	double precision, intent(in) :: rcoeffs(n), icoeffs(n)
	integer, intent(inout) :: nnew
	double precision, intent(inout) :: rnew(n), inew(n), alpha
	
	! compute variables
	integer :: ii
	double precision :: nrm, a, b
	
	! check size
	if(n < 3)then
		write(*,*) "n should be at least 3!"
		stop
	end if
	
	! check for zeros
	nnew = 0
	do ii=1,n
		nrm = abs(complex(rcoeffs(n+1-ii),icoeffs(n+1-ii)))
		if(nrm /= 0)then
			nnew = n+1-ii
			exit
		end if
	end do
	
	! check deflated size
	if(nnew == 0)then
		write(*,*) "enter a non-zero polynomial"
		return
	end if		

	! compute balance factor
	nrm = abs(complex(rcoeffs(nnew),icoeffs(nnew)))
	a = 1d0/dble(nnew)
	alpha = nrm**(a)
	a = 1d0/alpha
	b = a
	
	! compute new coefficients
	do ii=1,nnew
		rnew(ii) = b*rcoeffs(ii)
		inew(ii) = b*icoeffs(ii)
		b = a*b
	end do
        print*, "balancing alpha", alpha

end subroutine
