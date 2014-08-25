!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Aurentz² Mach³ Vandebril³ Watkins²
!
! ²Dept. Mathematics, Washington State University
! ³Dept. Computer Science, KU Leuven
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Last modified 22 August 2014
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! D Aurentz Mach Vandebril Watkins
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute coefficients from roots and check backward error 
! only for polynomials with real coefficients!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! N		degree of the polynomial
!
! POLY		array containing coefficients of P(x),
! 		POLY = [a_N-1, ... , a_0]
!
! RROOTS	array for real part of eigenvalues
!
! IROOTS	array for imaginary part of eigenvalues
!
! POLY2		array with coefficients obtained from roots
!
! ERR           array, absolute error |POLY(i)-POLY2(i)|
!
! RELERR        array, relative error |POLY(i)-POLY2(i)|/|POLY(i)|
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine backward(n, poly, rroots, iroots, poly2, err, relerr)
  
  use mpmodule
  implicit none

  !input variables
  integer, intent(in):: n
  !complex(kind(1d0)), intent(in) :: poly(n)
  double precision, intent(in) :: poly(n), rroots(n), iroots(n)
  !complex(kind(1d0)), intent(out) :: poly2(n)
  double precision, intent(inout) :: poly2(n), err(n), relerr(n)

  ! compute variables
  integer :: ii
  
  call rootstocoeffs(n, rroots, iroots, poly2)
  do ii=1,n
     err(ii)=abs(poly(ii)-poly2(ii))
     relerr(ii)=err(ii)/abs(poly(ii))
  end do
  
  print*, "poly", poly
  print*, "poly2", poly2
  

  print*, "abs err", err
  print*, "rel err", relerr


end subroutine backward
