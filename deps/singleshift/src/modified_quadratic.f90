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
! modified quadratic formula to compute eigenvalues
! of 2x2 matrix
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! block(2,2)   complex 2x2 matrix
!
! e1,e2        complex, eigenvalues of block
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine modified_quadratic(BLOCK,e1,e2)

	implicit none

	! input variables
	complex(kind(1d0)), intent(in) :: BLOCK(2,2)
	complex(kind(1d0)), intent(inout) :: e1, e2

	! compute variables
	complex(kind(1d0)) :: trace, detm, disc

	! compute intermediate values
	trace = BLOCK(1,1) + BLOCK(2,2)
	detm = BLOCK(1,1)*BLOCK(2,2) - BLOCK(2,1)*BLOCK(1,2)
	disc = zsqrt(trace*trace - 4d0*detm)

	! compute e1 and e2
	if(zabs(trace+disc) > zabs(trace-disc))then
		if(zabs(trace+disc) == 0)then
			e1 = complex(0d0,0d0)
			e2 = complex(0d0,0d0)
		else
			e1 = (trace+disc)/complex(2d0,0d0)
			e2 = detm/e1
		end if
	else
		if(zabs(trace-disc) == 0)then
			e1 = complex(0d0,0d0)
			e2 = complex(0d0,0d0)
		else
			e1 = (trace-disc)/complex(2d0,0d0)
			e2 = detm/e1
		end if
	end if

end subroutine
