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
! Balancing of a real polynomial
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
! N        (in) problem size
!
! rcoeffs  (in) real coefficients of the polynomial
! 
! nnew     (out) deflated problem size
! 
! rnew     (out) balanced coefficients
! 
! alpha    (out) balancing parameter
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine balance(n,rcoeffs,nnew,rnew,alpha)

	implicit none
	
	! input variables
	integer, intent(in) :: n
	double precision, intent(in) :: rcoeffs(n)
	integer, intent(inout) :: nnew
	double precision, intent(inout) :: rnew(n), alpha
	
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
		nrm = abs(rcoeffs(n+1-ii))
		if(nrm /= 0)then
			nnew = n+1-ii
			exit
		end if
	end do
	
	! check deflated size
	if(nnew == 0)then
		write(*,*) "enter a non-zero polynomial"
		stop
	end if		

	! compute balance factor
	nrm = abs(rcoeffs(nnew))
	a = 1d0/dble(nnew)
	alpha = nrm**(a)
	a = 1d0/alpha
	b = a
	
	! compute new coefficients
	do ii=1,nnew
		rnew(ii) = b*rcoeffs(ii)
		b = a*b
	end do
        print*, "balancing alpha", alpha

end subroutine
