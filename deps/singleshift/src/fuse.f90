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
! fuse bulge B into sequence of rotators Q with
! diagonal D
!
! the diagonal D is necessary to keep the rotations
! in Q in the form complex cosine, real sine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! n         problem size
!
! strt      start index of current block
! 
! stp       stop index of current block
!
! Q,D       generators of A
!
! B(3)      bulge to fuse
!
! flag      0: bulge from the right
!           1: bulge from the left
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fuse(n,strt,stp,Q,D,B,flag)

	implicit none

	! input variables
	integer, intent(in) :: n,strt,stp
	double precision, intent(inout) :: Q(3*n),D(2*n+2)
	double precision, intent(in) :: B(3)
	integer, intent(in) :: flag
	
	! compute variables
	integer :: k,ii
	double precision :: c1r, c1i, s1
	double precision :: c2r, c2i, s2
	double precision :: c3r, c3i, s3r, s3i
	double precision :: d1r, d1i
	double precision :: d2r, d2i
	double precision :: phr, phi
	double precision :: nrm
	double precision :: tol
	
	! set tol 
	tol = epsilon(1d0)
	
	! set inputs	
	c2r = B(1)
	c2i = B(2)
	s2 = B(3)
	
	! fusion from the right
	if(flag == 0)then
		! retrieve Q	
		k = 3*(stp-1)
		c1r = Q(k+1)
		c1i = Q(k+2)
		s1 = Q(k+3)
	
		! retrieve D
		k = 2*(stp-1)
		d1r = D(k+1)
		d1i = D(k+2)
		d2r = D(k+3)
		d2i = D(k+4)
		
		! compute givens product
		c3r = c1r*c2r - c1i*c2i - s1*s2
		c3i = c1r*c2i + c1i*c2r
		s3r = s1*c2r + s2*c1r
		s3i = s1*c2i - s2*c1i
		
		! compute phase
		nrm = abs(complex(s3r,s3i))
		if(nrm /= 0)then
			phr = s3r/nrm
			phi = s3i/nrm
		else
			phr = 1d0
			phi = 0d0
		end if
		
		! update Q
		c2r = c3r*phr + c3i*phi
		c2i = -c3r*phi + c3i*phr
		s2 = nrm
		nrm = c2r*c2r + c2i*c2i + s2*s2
		if(abs(nrm-1d0) > tol)then
			nrm = sqrt(nrm)
			c2r = c2r/nrm
			c2i = c2i/nrm
			s2 = s2/nrm
		end if
		k = 3*(stp-1)
		Q(k+1) = c2r
		Q(k+2) = c2i
		Q(k+3) = s2
		
		! update D
		c1r = phr*d1r - phi*d1i
		c1i = phr*d1i + phi*d1r
		nrm = c1r*c1r + c1i*c1i
		if(abs(nrm-1d0) > tol)then
			nrm = sqrt(nrm)
			c1r = c1r/nrm
			c1i = c1i/nrm
		end if
		c2r = phr*d2r + phi*d2i
		c2i = phr*d2i - phi*d2r
		nrm = c2r*c2r + c2i*c2i
		if(abs(nrm-1d0) > tol)then
			nrm = sqrt(nrm)
			c2r = c2r/nrm
			c2i = c2i/nrm
		end if
		k = 2*(stp-1)
		D(k+1) = c1r 
		D(k+2) = c1i
		D(k+3) = c2r
		D(k+4) = c2i
	
	! fusion from the left
	else
		! retrieve Q	
		k = 3*(strt-1)
		c1r = Q(k+1)
		c1i = Q(k+2)
		s1 = Q(k+3)
	
		! retrieve D
		k = 2*(strt-1)
		d1r = D(k+1)
		d1i = D(k+2)
		k = 2*(stp)
		d2r = D(k+1)
		d2i = D(k+2)
	
		! compute givens product
		c3r = c1r*c2r - c1i*c2i - s1*s2
		c3i = c1r*c2i + c1i*c2r
		s3r = s1*c2r + s2*c1r
		s3i = -(s1*c2i - s2*c1i)
		
		! compute phase
		nrm = abs(complex(s3r,s3i))
		if(nrm /= 0)then
			phr = s3r/nrm
			phi = s3i/nrm
		else
			phr = 1d0
			phi = 0d0
		end if
		
		! update Q
		c2r = c3r*phr + c3i*phi
		c2i = -c3r*phi + c3i*phr
		s2 = nrm
		nrm = c2r*c2r + c2i*c2i + s2*s2
		if(abs(nrm-1d0) > tol)then
			nrm = sqrt(nrm)
			c2r = c2r/nrm
			c2i = c2i/nrm
			s2 = s2/nrm
		end if
		k = 3*(strt-1)
		Q(k+1) = c2r
		Q(k+2) = c2i
		Q(k+3) = s2
		
		do ii=(strt+1),stp
			k = 3*(ii-1)
			c1r = Q(k+1)
			c1i = Q(k+2)
			s1 = Q(k+3)
			nrm = c1r*phr + c1i*phi
			c1i = -c1r*phi + c1i*phr
			c1r = nrm
			nrm = c1r*c1r + c1i*c1i + s1*s1
			if(abs(nrm-1d0) > tol)then
				nrm = sqrt(nrm)
				c1r = c1r/nrm
				c1i = c1i/nrm
				s1 = s1/nrm
			end if
			Q(k+1) = c1r
			Q(k+2) = c1i
			Q(k+3) = s1
		end do
		
		! update D
		c1r = phr*d1r - phi*d1i
		c1i = phr*d1i + phi*d1r
		nrm = c1r*c1r + c1i*c1i
		if(abs(nrm-1d0) > tol)then
			nrm = sqrt(nrm)
			c1r = c1r/nrm
			c1i = c1i/nrm
		end if
		c2r = phr*d2r + phi*d2i
		c2i = phr*d2i - phi*d2r
		nrm = c2r*c2r + c2i*c2i
		if(abs(nrm-1d0) > tol)then
			nrm = sqrt(nrm)
			c2r = c2r/nrm
			c2i = c2i/nrm
		end if
		k = 2*(strt-1)
		D(k+1) = c1r 
		D(k+2) = c1i
		k = 2*(stp)
		D(k+1) = c2r
		D(k+2) = c2i
	
	end if

end subroutine
