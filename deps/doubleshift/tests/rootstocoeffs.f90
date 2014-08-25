!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Last modified 22 August 2014
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! Compute the real coefficients of a polynomial
! from its roots
! roots have to be real or in conj complex pairs
! 
! uses MPFUN
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! n        (in) problem size
!
! rroots   (in) real parts of the roots
! iroots   (in) imag parts of the roots
! 
! coeffs   (out) real coefficients of the polynomial
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rootstocoeffs(n,rroots,iroots,coeffs)
  use mpmodule
  implicit none

  integer, intent(in):: n
  double precision, intent(in) :: rroots(n), iroots(n)
  double precision, intent(inout) :: coeffs(n)
  type (mp_real) :: rr(n), ir(n), c(n), alpha, beta, zero

  ! compute variables
  integer :: ii, jj

  ! convert to mp-type
  zero = 0.0d0
  do ii=1,n
     rr(ii) = rroots(ii)
     ir(ii) = iroots(ii)
     c(ii)  = 0.0d0
  end do
  
  ii=1
  
  do while (ii .LE. n)
     if (ir(ii) .EQ. zero) then
        alpha = -rr(ii)
        do jj=ii,1,-1
           if (jj==1) then
              c(jj) = c(jj) + alpha*1d0
           else
              c(jj) = c(jj) + alpha*c(jj-1)
           end if
        end do
        ii=ii+1
     else
        alpha = - rr(ii)*2
        beta = rr(ii)**2 + ir(ii)**2
        do jj=ii+1,1,-1
           if (jj == 2) then
              c(jj) = c(jj) + alpha*c(jj-1) + beta
           elseif (jj == 1) then
              c(jj) = c(jj) + alpha 
           else 
              c(jj) = c(jj) + alpha*c(jj-1) + beta*c(jj-2)
           endif
        enddo
        ii=ii+2
     endif
  enddo
  

  ! convert to double precision
  do ii=1,n
     coeffs(ii) = c(ii)
  end do

end subroutine rootstocoeffs
