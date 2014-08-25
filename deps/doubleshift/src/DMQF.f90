!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Aurentz² Mach³ Vandebril³ Watkins²
!
! ²Dept. Mathematics, Washington State University
! ³Dept. Computer Science, KU Leuven
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Last modified 22 August 2014
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! D Modified Quadratic Formula
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine computes the eigenvalues of a 2x2 matrix using the 
! modified quadratic formula.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! BLOCK		2x2 block matrix
!
! re1, re2	real parts of eig1 and eig2
!
! ie1, rie2	imaginary parts of eig1 and eig2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine DMQF(BLOCK,re1,ie1,re2,ie2)

	implicit none

	! input variables
	double precision, intent(in) :: BLOCK(2,2)
	double precision, intent(inout) :: re1, ie1, re2, ie2

	! compute variables
	double precision :: trace, detm, disc

	! compute intermediate values
	trace = BLOCK(1,1) + BLOCK(2,2)
	detm = BLOCK(1,1)*BLOCK(2,2) - BLOCK(2,1)*BLOCK(1,2)
	disc = trace*trace - 4d0*detm

	! compute e1 and e2
	! complex eigenvalues
	if(disc < 0)then
		re1 = trace/2d0
		ie1 = sqrt(-disc)/2d0
		re2 = re1
		ie2 = -ie1
	! real eigenvalues
	else if(abs(trace+sqrt(disc)) > abs(trace-sqrt(disc)))then
		if(abs(trace+sqrt(disc)) == 0)then
			re1 = 0d0
			ie1 = 0d0
			re2 = 0d0
			ie2 = 0d0
		else
			re1 = (trace+sqrt(disc))/2d0
			ie1 = 0d0
			re2 = detm/re1
			ie2 = 0d0
		end if
	else
		if(abs(trace-sqrt(disc)) == 0)then
			re1 = 0d0
			ie1 = 0d0
			re2 = 0d0
			ie2 = 0d0
		else
			re1 = (trace-sqrt(disc))/2d0
			ie1 = 0d0
			re2 = detm/re1
			ie2 = 0d0
		end if
	end if

end subroutine

