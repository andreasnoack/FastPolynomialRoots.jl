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
! generates a rotation G zeroing b if applied to
! (a;b)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! a         complex
!
! b         complex
! 
! G(3)      rotation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine crgivens(a,b,G)

	implicit none

	! input variables
	complex(kind(1d0)), intent(inout) :: a,b
	double precision, intent(inout) :: G(3)
	
	! compute variables
	double precision :: t1r, t1i, t2r, t2i
	double precision :: nrm, T(2)
	complex(kind(1d0)) :: temp
	
	! BLAS DNRM2
	double precision :: dnrm2
	
	! check for 0
        nrm = abs(b)
	if(nrm == 0)then
		G(1) = 1d0
		G(2) = 0d0
		G(3) = 0d0
		return
		
	end if

        t2r = dble(b)
	t2i = dimag(b)
        nrm = abs(b)
        t2r = t2r/nrm
        t2i = t2i/nrm
 
	! store nrm
	G(3) = nrm
	
	! compute complex part
	G(1) = dble(a*complex(t2r,-t2i))
	G(2) = dimag(a*complex(t2r,-t2i))
	
	! normalize
	nrm = dnrm2(3,G,1)
	G = G/nrm
	
	! update a and b
	a = a*complex(G(1),-G(2)) + b*complex(G(3),0d0)
	b = complex(0d0,0d0)		

end subroutine
