!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Aurentz² Mach³ Vandebril³ Watkins²
!
! ²Dept. Mathematics, Washington State University
! ³Dept. Computer Science, KU Leuven
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Last modified 22 August 2014
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! build initial bulge from shifts
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! n         problem size
!
! strt      start index of current block
! 
! bulge(3)  (out) rotation, initial bulge
!
! shift     complex shift
!
! Q,D,C,B   generators of A
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine buildbulge(n,strt,bulge,shift,Q,D,C,B)

	implicit none

	! input variables
	integer, intent(in) :: n, strt
	double precision, intent(in) :: Q(3*n), D(2*n), C(3*n), B(3*n)
	complex(kind(1d0)), intent(in) :: shift
	double precision, intent(inout) :: bulge(3)

	! compute variables
	complex(kind(1d0)) :: block(2,2)

	! top block
	call diagblock(n,strt,block,Q,D,C,B)
	
	! shift
	block(1,1) = block(1,1) - shift
	
	! bulge
	call crgivens(block(1,1),block(2,1),bulge)

end subroutine
