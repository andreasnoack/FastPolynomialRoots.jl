!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Aurentz² Mach³ Vandebril³ Watkins²
!
! ²Dept. Mathematics, Washington State University
! ³Dept. Computer Science, KU Leuven
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Last modified 22 August 2014
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! D Fuse Givens Rotations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine fuses two givens rotations Q1, Q2 and stores the 
! output in Q1 or Q2 depending on FTYPE
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! FTYPE		FTYPE = 0 if stored in Q1, FTYPE = 1 if stored in Q2
!
! Q1, Q2	generators for two givens rotations
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine DFGR(FTYPE,Q1,Q2)

	implicit none

	! input variables
	integer, intent(in) :: FTYPE
	double precision, intent(inout) :: Q1(2), Q2(2)

	! compute variables
	double precision :: temp

	! store in Q1
	if(FTYPE == 0)then
		! compute new generators
		temp = Q1(1)*Q2(1) - Q1(2)*Q2(2)
		Q1(2) = Q1(2)*Q2(1) + Q1(1)*Q2(2)
		Q1(1) = temp

		! enforce orthogonality
		!temp = dsqrt(Q1(1)**2 + Q1(2)**2)
		!Q1(1) = Q1(1)/temp
		!Q1(2) = Q1(2)/temp
	! store in Q2
	else if(FTYPE == 1)then
		! compute new generators
		temp = Q1(1)*Q2(1) - Q1(2)*Q2(2)
		Q2(2) = Q1(2)*Q2(1) + Q1(1)*Q2(2)
		Q2(1) = temp

		! enforce orthogonality
		!temp = dsqrt(Q2(1)**2 + Q2(2)**2)
		!Q2(1) = Q2(1)/temp
		!Q2(2) = Q2(2)/temp
	! bad input
	else
		print*, "Not a valid input for FTYPE!"		
		return
	end if

end subroutine
