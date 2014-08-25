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
! chase bulge thru A
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! n         problem size
!
! strt      start index of current block
! 
! stp       stop index of current block
!
! bulge(3)  rotation, initial bulge
!
! Q,D,C,B   generators of A
!
! tr        first tr rotations in B are equal to C*
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine chasebulge(n,strt,stp,bulge,Q,D,C,B,tr)

	implicit none

	! input variables
	integer, intent(in) :: n, strt, stp
	double precision, intent(inout) :: Q(3*n), D(2*n+2), C(3*n), B(3*n)
	double precision, intent(inout) :: bulge(3)
	integer, intent(inout) :: tr

	! compute variables
	integer :: ind1,ind2,ii,str,ll
	double precision :: binv(3)
	
	double precision :: T1(3), T2(3), T3(3), D1(4), error
	double precision :: Q2(6), D2(6), bulge2(3)
	
	! bulge inverse
	binv(1) = bulge(1)
	binv(2) = -bulge(2)
	binv(3) = -bulge(3)
	
	! fusion at top
	call fuse(n,strt,stp,Q,D,binv,1)
	
	! main chasing loop
	do ii=strt,(stp-1)
	
		! set indices for B and C
		ind1 = 3*(ii-1) + 1
		ind2 = ind1+2
	
  
                if (ii<tr) then
                   ! B(ii) equals C(ii)*
                   bulge2 = bulge
                   ! through B
                   call dto4(B(ind1:ind2),B((ind1+3):(ind2+3)),bulge2)
                   
                   C(ind1) = B(ind1)
                   C(ind1+1) = -B(ind1+1)
                   C(ind1+2) = -B(ind1+2)
                   C(ind1+3) = B(ind1+3)
                   C(ind1+4) = -B(ind1+4)
                   C(ind1+5) = -B(ind1+5)
                else
                   ! through B
                   call dto4(B(ind1:ind2),B((ind1+3):(ind2+3)),bulge)
                   
                   ! through C
                   call dto4(C((ind1+3):(ind2+3)),C(ind1:ind2),bulge)
		end if

		! through diag
		call throughdiag(n,ii,D,bulge)

		! through Q
		call dto4(Q(ind1:ind2),Q((ind1+3):(ind2+3)),bulge)
		
	end do

	! set indices for B and C
	ind1 = 3*(stp-1) + 1
	ind2 = ind1+2
	
	! through B
	call dto4(B(ind1:ind2),B((ind1+3):(ind2+3)),bulge)

	! through C
	call dto4(C((ind1+3):(ind2+3)),C(ind1:ind2),bulge)

	! through diag
	call throughdiag(n,stp,D,bulge)

	! fusion at bottom
	call fuse(n,strt,stp,Q,D,bulge,0)

end subroutine
