!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Aurentz² Mach³ Vandebril³ Watkins²
!
! ²Dept. Mathematics, Washington State University
! ³Dept. Computer Science, KU Leuven
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Last modified 22 August 2014
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! D Givens Rotation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine computes c and s such that,
!
! -s*a + c*b = 0, |c|^2 + |s|^2 = 1,
!
! and 
!
! r = sqrt(|a|^2 + |b|^2).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! a,b		input coefficients
!
! c,s		givens rotation generators
!
! r		output norm
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine DGR(a,b,c,s,r)

  implicit none
  
  ! input variables
  double precision, intent(in) :: a,b
  double precision, intent(inout) :: c,s,r

  ! tol to avoid unnecessary sqrt computations not implemented
  ! double precision :: tol
  ! set tol
  ! tol = epsilon(1d0)
  
  if (b == 0) then
     if (a<0) then 
        c = -1d0
        s = 0d0
        r = -a
     else
        c = 1d0
        s = 0d0
        r = a
     endif
  else if (dabs(a) >= dabs(b)) then
     s = b/a
     r = dsqrt(1.d0 + s**2)
     if (a<0) then
        c = -1.d0/r
        s =  s*c
        r = -a*r
     else
        c =  1.d0/r
        s =  s*c
        r =  a*r
     end if
  else
     c = a/b;
     r = dsqrt(1.d0 + c**2)
     if (b<0) then
        s = -1.d0/r
        c =  c*s
        r = -b*r
     else
        s =  1.d0/r
        c =  c*s
        r =  b*r
     end if
  end if
           
end subroutine
