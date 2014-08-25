! ****************************************************************************
!
! This subroutine generates a one dimensional complex array whose entries are 
! normally distributed with mean 0 and variance 1 in both the real and imaginary parts
!
! ****************************************************************************
subroutine cnormalpoly(degree,poly)

  implicit none
  
  integer, intent(in) :: degree
  complex(kind(1d0)), intent(inout) :: poly(degree) 
	
  double precision :: u,v,s,pi = 3.141592653589793239d0
  integer :: i,j
  
  do i=1,degree
     do j=1,20
        
        call random_number(u)
        call random_number(v)
	
        s = u**2 + v**2
	
        if(s > 0 .and. s < 1)then				
           poly(i) = complex(dcos(2.d0*pi*v)*dsqrt(-2.d0*dlog(u)),dsin(2.d0*pi*v)*dsqrt(-2.d0*dlog(u)))
           exit
        end if
     end do
     
  end do


end subroutine

