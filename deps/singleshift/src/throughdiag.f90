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
! pass rotation B through diagonal D
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! n         problem size
!
! k         index of bulge
! 
!
! D         nx2 real matrix, generator of A, 
!
! B(3)      bulge to pass thru diagonal
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine throughdiag(n,k,D,B)

	implicit none

	! input variables
	integer, intent(in) :: n,k
	double precision, intent(inout) :: D(2*n+2), B(3)
	
	! compute variables
	integer :: strt
	double precision :: c1r, c1i, s1
	double precision :: d1r, d1i
	double precision :: d2r, d2i
	double precision :: nrm
	double precision :: tol
	
	! set tol 
	tol = epsilon(1d0)
	
	! set inputs
	c1r = B(1)
	c1i = B(2)
	s1 = B(3)
	
	! retrieve D
	strt = 2*(k-1)
	d1r = D(strt+1)
	d1i = D(strt+2)
	d2r = D(strt+3)
	d2i = D(strt+4)
	
	! pass through diagonal
	nrm = (d1r*d2r + d1i*d2i)*c1r - (-d1r*d2i + d1i*d2r)*c1i
	c1i = (d1r*d2r + d1i*d2i)*c1i + (-d1r*d2i + d1i*d2r)*c1r
	c1r = nrm
		
	! renormalize
	nrm = c1r*c1r + c1i*c1i + s1*s1
	if(abs(nrm-1d0) > tol)then
		nrm = sqrt(nrm)
		c1r = c1r/nrm
		c1i = c1i/nrm
		s1 = s1/nrm
	end if
	
	! set B
	B(1) = c1r
	B(2) = c1i
	B(3) = s1
	
	! set D
	strt = 2*(k-1)
	D(strt+1) = d2r 
	D(strt+2) = d2i
	D(strt+3) = d1r
	D(strt+4) = d1i
	
end subroutine
