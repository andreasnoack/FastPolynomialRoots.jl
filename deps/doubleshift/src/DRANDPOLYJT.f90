!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Aurentz² Mach³ Vandebril³ Watkins²
!
! ²Dept. Mathematics, Washington State University
! ³Dept. Computer Science, KU Leuven
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Last modified 22 August 2014
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! D RANDom POLYnomial following Jenkins Traub (iv)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine generates random polynomial as 
! (iv) in [Jenkins, Traub 1970]:
!
! "(iv) polynomials whose coefficients are chosen
!  randomly by taking the mantissa and exponents
!  from separate uniform distributions. The 
!  resulting polynomials have widely varying zeros 
!  and hence yield a reasonable test that the 
!  program has wide applications."
!
! [Jenkins, Traub 1970] M. A. Jenkins and J. F. 
!    Traub, Principles for testing polynomial 
!    zerofinding programs, ACM Transactions on 
!    Mathematical Software, 1 (1975), pp. 26–34.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! N		degree of the polynomial
!
! POLY		array coefficients of P(x),
! 		POLY = [a_N-1, ... , a_0]
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine DRANDPOLYJT(N,POLY,R)

  implicit none
  
  ! input variables
  integer, intent(in) :: N
  double precision, intent(inout) :: POLY(N), R
  
  ! compute variables
  double precision :: u,v,s,pi = 3.141592653589793239d0
  integer :: ii,jj
  
  
  do ii=1,N
     
     call random_number(u) ! uniform distribution in [0,1)
     call random_number(v) ! uniform distribution in [0,1)
     
     POLY(ii) = (2d0*u-1d0) * 10**(2d0*R*v-R)
     
  end do
  
  
end subroutine DRANDPOLYJT
