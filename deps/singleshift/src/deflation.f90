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
! checks for deflations and deflates if possible
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! n         problem size
!
! strt      start index of current block
! 
! stp       stop index of current block
!
! zero      index of last zero above current block
!
! Q,D,C,B   generators of A
!
! its, itcnt iteration counter
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine deflation(n,strt,stp,zero,Q,D,C,B,its,itcnt)

	implicit none

	! input variables
	integer, intent(in) :: n
	integer, intent(inout) :: strt, stp, zero, itcnt
	double precision, intent(inout) :: Q(3*n), D(2*n+2), C(3*n), B(3*n)
	integer, intent(inout) :: its(n)

	! compute variables
	integer :: ii, ind1, ll, jj, k
	double precision :: tol, nrm
	double precision :: d1r, d1i, c1r, c1i, s
	
	! set tolerance
	tol = 1d0*epsilon(1d0)

	! check for deflation
	do ii=1,stp
		ind1 = 3*(stp-ii)

		nrm = abs(Q(ind1+3))
		
		if(nrm < tol)then
		
			! set sub-diagonal to 0
			Q(ind1+3) = 0d0
			
			! update first diagonal
			c1r = Q(ind1+1)
			c1i = Q(ind1+2)
			Q(ind1+1) = 1d0
			Q(ind1+2) = 0d0
			
			ind1 = 2*(stp-ii)	
			d1r = D(ind1+1)
			d1i = D(ind1+2)
			
			nrm = c1r*d1r - c1i*d1i
			d1i = c1r*d1i + c1i*d1r
			d1r = nrm
			nrm = d1r*d1r + d1i*d1i
			if(abs(nrm-1d0) > tol)then
				nrm = sqrt(nrm)
				d1r = d1r/nrm
				d1i = d1i/nrm
			end if
			
			D(ind1+1) = d1r
			D(ind1+2) = d1i
			
			! 1x1 deflation
			if(ii == 1)then
				! update second diagonal
				ind1 = 2*(stp-ii)	
				d1r = D(ind1+3)
				d1i = D(ind1+4)
			
				nrm = c1r*d1r + c1i*d1i
				d1i = c1r*d1i - c1i*d1r
				d1r = nrm
				nrm = d1r*d1r + d1i*d1i
				if(abs(nrm-1d0) > tol)then
					nrm = sqrt(nrm)
					d1r = d1r/nrm
					d1i = d1i/nrm
				end if
			
				D(ind1+3) = d1r
				D(ind1+4) = d1i			
			
			! 2x2 or bigger
			else
				! update Q
				do ll=(stp+1-ii),(stp-1)
					ind1 = 3*(ll)
					d1r = Q(ind1+1)
					d1i = Q(ind1+2)
					s = Q(ind1+3)
			
					nrm = c1r*d1r + c1i*d1i
					d1i = c1r*d1i - c1i*d1r
					d1r = nrm
					nrm = d1r*d1r + d1i*d1i + s*s
					if(abs(nrm-1d0) > tol)then
						nrm = sqrt(nrm)
						d1r = d1r/nrm
						d1i = d1i/nrm
						s = s/nrm
					end if
			
					Q(ind1+1) = d1r
					Q(ind1+2) = d1i
					Q(ind1+3) = s
				end do			
			
				! update second diagonal
				ind1 = 2*(stp)	
				d1r = D(ind1+1)
				d1i = D(ind1+2)
			
				nrm = c1r*d1r + c1i*d1i
				d1i = c1r*d1i - c1i*d1r
				d1r = nrm
				nrm = d1r*d1r + d1i*d1i
				if(abs(nrm-1d0) > tol)then
					nrm = sqrt(nrm)
					d1r = d1r/nrm
					d1i = d1i/nrm
				end if
			
				D(ind1+1) = d1r
				D(ind1+2) = d1i
			end if			
	
			! update indices
			zero = stp+1-ii
			strt = zero + 1

			! store it_count
			its(zero) = itcnt
			itcnt = 0
			

			exit
		end if
	end do
		
end subroutine
