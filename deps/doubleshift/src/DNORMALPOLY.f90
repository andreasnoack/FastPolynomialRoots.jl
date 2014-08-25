!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Aurentz² Mach³ Vandebril³ Watkins²
!
! ²Dept. Mathematics, Washington State University
! ³Dept. Computer Science, KU Leuven
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Last modified 22 August 2014
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! D Normally Distributed Polynomial
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine generates double precision, normally distributed
! coefficients for a monic polynomial of degree N. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! N		degree of the polynomial
!
! POLY		array containing coefficients of P(x),
! 		POLY = [a_N-1, ... , a_0]
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine DNORMALPOLY(N,POLY)

  implicit none
  
  ! input variables
  integer, intent(in) :: N
  double precision, intent(inout) :: POLY(N) 
  
  ! compute variables
	double precision :: u,v,s,pi = 3.141592653589793239d0
 integer :: ii,jj
 
 ! generate normally distributed numbers using Box-Muller
 do ii=1,N
    do jj=1,200
       
       call random_number(u)
       call random_number(v)
       
       s = u**2 + v**2
       
       if(s > 0d0 .and. s < 1d0)then				
          POLY(ii) = dcos(2.d0*pi*v)*dsqrt(-2.d0*dlog(u))
          exit
       end if
    end do

 end do
 

end subroutine
