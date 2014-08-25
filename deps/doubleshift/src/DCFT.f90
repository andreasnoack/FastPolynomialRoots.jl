!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Aurentz² Mach³ Vandebril³ Watkins²
!
! ²Dept. Mathematics, Washington State University
! ³Dept. Computer Science, KU Leuven
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Last modified 22 August 2014
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! D Compute First Transformation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine computes the first two givens rotation to initialize 
! the bulge chase. The shifts rho1 and rho2 are expected to be both 
! real or a complex conjugate pair.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! N		problem size
!
! str		index for first relevant column of A
!
! Q,C,B		generators for A
!
! rrho1, rrho2	real parts of shifts
!
! irho1, irho2	imaginary parts of shifts
!
! B1, B2	generators for the first two givens transforms
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine DCFT(N,str,QCB,rrho1,irho1,rrho2,irho2,B1,B2)
  
  implicit none
  
  ! input variables
  integer, intent(in) :: N, str
  double precision, intent(in) :: QCB(6*N)
  double precision, intent(inout) :: B1(2), B2(2)
  double precision, intent(in) :: rrho1, irho1, rrho2, irho2
  
  ! compute variables
  double precision :: scrap, COL(3), scrap2
  double precision :: TEMP(3,2), T(3,2)
  
  ! compure first two columns of A
  TEMP = 0d0
  call DCDB(N,str,T,QCB)
  TEMP(1:2,1:2) = T(1:2,1:2)
  call DCDB(N,str+1,T,QCB)
  TEMP(3,2) = T(2,1)

  ! compute (A-rho1)(A-rho2)e_1
  COL(1) = TEMP(1,1)*TEMP(1,1) + TEMP(1,2)*TEMP(2,1) + rrho1*rrho2 - irho1*irho2 - TEMP(1,1)*(rrho1+rrho2)
  COL(2) = TEMP(2,1)*(TEMP(1,1)+TEMP(2,2)-(rrho1+rrho2))
  COL(3) = TEMP(2,1)*TEMP(3,2)
  
  ! compute first two givens rotations
  call DGR(COL(2),COL(3),B1(1),B1(2),scrap)
  call DGR(COL(1),scrap,B2(1),B2(2),scrap2)

end subroutine

