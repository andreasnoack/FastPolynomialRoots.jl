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
! computes generators of Companion matrix of p
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! n         problem size
!
! rcoeffs, icoeffs  ... coefficients of p:
!
! a_j = rcoeffs(j) + i*icoeffs(j)
!
! p(z) = z^n + a_1 z^{n-1} + a_2 z^(n-2} + ... + a_{n-1} z + a_n
!
!
! Q,D,C,B   (out) generators of A
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine factor(n,rcoeffs,icoeffs,Q,D,C,B)

	implicit none
	
	! input variables
	integer, intent(in) :: n
	double precision, intent(in) :: rcoeffs(n), icoeffs(n)
	double precision, intent(inout) :: Q(3*n),D(2*(n+1)),C(3*n),B(3*n)
	
	! compute variables
	integer :: ii,strt,stp
	complex(kind(1d0)) :: t1, t2
	double precision :: phr, phi, nrm
	
	! check n
	if(n <= 2)then
		write(*,*) "In factor: n must be > 2!"
		stop
	end if
	
	! set Q
	do ii=1,(n-1)
		Q(3*(ii-1) + 1) = 0d0
		Q(3*(ii-1) + 2) = 0d0
		Q(3*(ii-1) + 3) = 1d0
	end do
	Q(3*(n-1) + 1) = 1d0
	Q(3*(n-1) + 2) = 0d0
	Q(3*(n-1) + 3) = 0d0
	
	! set D
	do ii=1,n+1
		D(2*(ii-1) + 1) = 1d0
		D(2*(ii-1) + 2) = 0d0
	end do
	
	! initialize B and C
	t1 = complex((-1d0)**(n),0d0)*complex(rcoeffs(n),icoeffs(n))
	t2 = complex((-1d0)**(n-1),0d0)
	strt = 3*(n-1) + 1
	stp = strt + 2
	call crgivens(t1,t2,C(strt:stp))
	B(strt) = C(strt)
	B(strt+1) = -C(strt+1)
	B(strt+2) = -C(strt+2)
	
	do ii=2,n
		t2 = t1
		t1 = -complex(rcoeffs(ii-1),icoeffs(ii-1))
		strt = 3*(n-ii) + 1
		stp = strt + 2
		call crgivens(t1,t2,C(strt:stp))
		B(strt) = C(strt)
		B(strt+1) = -C(strt+1)
		B(strt+2) = -C(strt+2)

	end do
	
	! B_n
	strt = 3*(n-1) + 1
	t1 = complex(B(strt),B(strt+1))
	nrm = abs(t1)
	if(nrm == 0)then
		phr = 1d0
		phi = 0d0
	else
		phr = B(strt)/nrm
		phi = B(strt+1)/nrm
	end if
	B(strt) = -B(strt+2)*phr
	B(strt+1) = -B(strt+2)*phi
	B(strt+2) = nrm
	
	! update D
	nrm = (-1d0)**(n)
	nrm = nrm/abs(nrm)
	D(2*(n) + 1) = nrm*phr
	D(2*(n) + 2) = nrm*phi
	D(2*(n-2) + 1) = nrm*phr
	D(2*(n-2) + 2) = -nrm*phi

end subroutine
