!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! computes n random complex coefficients
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine normalpoly(n,rcoeffs,icoeffs)

	implicit none
	
	! input variables
	integer, intent(in) :: n
	double precision, intent(inout) :: rcoeffs(n), icoeffs(n) 
	
	! compute variables
	double precision :: u,v,s,pi = 3.141592653589793239d0
	integer :: ii,jj

	do ii=1,n
		do jj=1,100
			
			call random_number(u)
			call random_number(v)
	
			s = u**2 + v**2
	
			if(s > 0 .and. s < 1)then				
				rcoeffs(ii) = dcos(2.d0*pi*v)*dsqrt(-2.d0*dlog(u))
				icoeffs(ii) = dsin(2.d0*pi*v)*dsqrt(-2.d0*dlog(u))
				exit
			end if
		end do

	end do

end subroutine

