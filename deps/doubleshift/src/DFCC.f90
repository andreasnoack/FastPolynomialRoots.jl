!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Aurentz² Mach³ Vandebril³ Watkins²
!
! ²Dept. Mathematics, Washington State University
! ³Dept. Computer Science, KU Leuven
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Last modified 22 August 2014
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! D Factor Column Companion
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine computes a factorization of the column companion 
! matrix for P(x),
!
! P(x) = x^N-1 + a_N-2 x^N-2 + ... + a_1 x + a_0.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! N-1		degree of the polynomial
!
! POLY		array containing coefficients of P(x),
! 		POLY = [a_N-2, ... , a_0]
!
! QCB		array of generators for A
!
! ALPHA         parameter for balancing, currently not implemented               
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine DFCC(N,POLY,QCB,ALPHA)

	implicit none
  
	! input variables
	integer, intent(in) :: N
	double precision, intent(inout) :: POLY(N)
	double precision, intent(inout) :: QCB(6*N)
	double precision, intent(out) :: ALPHA
  
	! compute variables
	integer :: ii,strt
	double precision :: temp1,temp2
  
	! balance parameter (currently not used)
	ALPHA = 1d0

	! build Q
	QCB = 0d0
	do ii=1,N-1
 		strt = 6*(ii-1)
		QCB(strt+2) = 1d0
	end do
	strt = 6*(N-1)
	QCB(strt+1) = 1d0
  
	! build C and B
	strt = 6*(N-1)+2
	temp1 = dble((-1)**(N-1))
	call DGR(dble((-1)**(N))*POLY(N),temp1,QCB(strt+1),QCB(strt+2),temp2)
	QCB(strt+3) = QCB(strt+2)*dble((-1)**(N))
	QCB(strt+4) = QCB(strt+1)*dble((-1)**(N))

	do ii=2,N
		strt = 6*(N-ii)+2
		temp1 = temp2
		call DGR(-POLY(ii-1),temp1,QCB(strt+1),QCB(strt+2),temp2)
		QCB(strt+3) = QCB(strt+1)
		QCB(strt+4) = -QCB(strt+2)
	end do
  
end subroutine
