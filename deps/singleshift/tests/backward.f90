!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Aurentz² Mach³ Vandebril³ Watkins²
!
! ²Dept. Mathematics, Washington State University
! ³Dept. Computer Science, KU Leuven
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Last modified 22 Augst 2014
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! D Aurentz Mach Vandebril Watkins
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute coefficients from and check backward error 
! only for polynomials with real coefficients!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! N		degree of the polynomial
!
! POLY		complex array containing coefficients of P(x),
! 		POLY = [a_N-1, ... , a_0]
!
! RROOTS	real array for real part of eigenvalues
!
! IROOTS	real array for imaginary part of eigenvalues
!
! POLY2		complex array with coefficients obtained from roots
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
  complex(kind(1d0)), intent(in) :: poly(n)
  double precision, intent(in) :: rroots(n), iroots(n)
  complex(kind(1d0)), intent(out) :: poly2(n)
  double precision, intent(inout) :: err(n), relerr(n)

  ! compute variables
  integer :: ii
  double precision :: poly3r(n),poly3i(n)


!!$  do ii=1,n
!!$     print*, rroots(ii), iroots(ii)
!!$  end do

  call zrootstocoeffs(n, rroots, iroots, poly3r, poly3i)
  do ii=1,n
     poly2(ii)=complex(poly3r(ii),poly3i(ii))
     err(ii)=abs(poly(ii)-poly2(ii))
     relerr(ii)=err(ii)/abs(poly(ii))
  end do
  
  print*, "poly", poly
  print*, "poly2", poly2
  

  print*, "abs err", err
  print*, "rel err", relerr


end subroutine backward


subroutine zrootstocoeffs(n,rroots,iroots,rcoeffs,icoeffs)
  use mpmodule
  implicit none

  integer, intent(in):: n
  double precision, intent(in) :: rroots(n), iroots(n)
  double precision, intent(inout) :: rcoeffs(n),icoeffs(n)
  type (mp_real) :: rr(n), ir(n), rc(n), ic(n), ralpha, ialpha

  ! compute variables
  integer :: ii, jj

  ! convert to mp-type
  do ii=1,n
     rr(ii) = rroots(ii)
     ir(ii) = iroots(ii)
     rc(ii)  = 0.0d0
     ic(ii)  = 0.0d0
  end do
  
  ii=1
  
  do while (ii .LE. n)
     ralpha = -rr(ii)
     ialpha = -ir(ii)
     do jj=ii,1,-1
        if (jj==1) then
           rc(jj) = rc(jj) + ralpha*1d0
           ic(jj) = ic(jj) + ialpha*1d0
        else
           rc(jj) = rc(jj) + ralpha*rc(jj-1) - ialpha*ic(jj-1)
           ic(jj) = ic(jj) + ralpha*ic(jj-1) + ialpha*rc(jj-1)
        end if
     end do
     ii=ii+1
  enddo
  

  ! convert to double precision
  do ii=1,n
     rcoeffs(ii) = rc(ii)
     icoeffs(ii) = ic(ii)
  end do

end subroutine zrootstocoeffs
