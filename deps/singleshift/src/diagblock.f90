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
! computes A(k:k+1,k:k+1) from the generators
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! n         problem size
!
! k         index
!
! block(2,2) (out) block
!
! Q,D,C,B   generators of A
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine diagblock(n,k,block,Q,D,C,B)

	implicit none

	! input variables
	integer, intent(in) :: n, k
	double precision, intent(in) :: Q(3*n), D(2*n), C(3*n), B(3*n)
	complex(kind(1d0)), intent(inout) :: block(2,2)

	! compute variables
	integer :: strt
	complex(kind(1d0)) :: temp, R(3,2)
	complex(kind(1d0)) :: H(2,2)
	
	! check k
	if(k > n-1)then
		write(*,*) 'k must be <= n-1 in diagblock'
   		stop
	end if
	
	R = complex(0d0,0d0)

	! case K = 1
	if(k == 1)then
		! R
		R(1,1) = -complex(B(3)/C(3),0d0)
		R(2,2) = -complex(B(6)/C(6),0d0)
		R(1,2) = (complex(-B(1),B(2))*complex(B(4),B(5)) + R(2,2)*complex(C(1),C(2))*complex(C(4),-C(5)))/complex(C(3),0d0)

		! diag
		R(1,:) = complex(D(1),D(2))*R(1,:)
		R(2,:) = complex(D(3),D(4))*R(2,:)
		
		! block
		R(2,2) = complex(Q(4),Q(5))*R(2,2)
			
		block(1,1) = complex(Q(1),Q(2))
		block(2,1) = complex(Q(3),0d0)
		block(1,2) = complex(-Q(3),0d0)
		block(2,2) = complex(Q(1),-Q(2))
		
		block = matmul(block,R(1:2,1:2))

	! other cases
	else
	
		! first column of R
		strt = 3*(k-1) + 1
		R(2,1) = complex(-B(strt+2)/C(strt+2),0d0)
		R(1,1) = (complex(-B(strt-3),B(strt-2))*complex(B(strt),B(strt+1)) &
			+ R(2,1)*complex(C(strt-3),C(strt-2))*complex(C(strt),-C(strt+1)))/complex(C(strt-1),0d0)
			
		! second column of R
		strt = 3*k + 1
		R(3,2) = complex(-B(strt+2)/C(strt+2),0d0)
		
		R(2,2) = (complex(-B(strt-3),B(strt-2))*complex(B(strt),B(strt+1)) &
			+ R(3,2)*complex(C(strt-3),C(strt-2))*complex(C(strt),-C(strt+1)))/complex(C(strt-1),0d0)
			
		R(1,2) = (complex(B(strt-6),-B(strt-5))*complex(B(strt-1),0d0)*complex(B(strt),B(strt+1)) - &
			complex(C(strt-6),C(strt-5))/complex(C(strt-1),0d0)* &
			(complex(C(strt-3),-C(strt-2))*complex(B(strt-3),-B(strt-2))*complex(B(strt),B(strt+1)) - &
			complex(C(strt),-C(strt+1))*R(3,2)))/complex(C(strt-4),0d0)

		! diag
		strt = 2*(k-1) + 1
		R(1,:) = complex(D(strt-2),D(strt-1))*R(1,:)
		R(2,:) = complex(D(strt),D(strt+1))*R(2,:)
		R(3,:) = complex(D(strt+2),D(strt+3))*R(3,:)

		! start index for Q
		strt = 3*(k-1) + 1

		! block	
		R(3,2) = complex(Q(strt+3),Q(strt+4))*R(3,2)

		block(1,1) = complex(Q(strt),Q(strt+1))
		block(2,1) = complex(Q(strt+2),0d0)
		block(1,2) = complex(-Q(strt+2),0d0)
		block(2,2) = complex(Q(strt),-Q(strt+1))
	
		R(2:3,1:2) = matmul(block,R(2:3,1:2))

		block(1,1) = complex(Q(strt-3),Q(strt-2))
		block(2,1) = complex(Q(strt-1),0d0)
		block(1,2) = complex(-Q(strt-1),0d0)
		block(2,2) = complex(Q(strt-3),-Q(strt-2))
			
		R(1:2,1:2) = matmul(block,R(1:2,1:2))
	
		block = R(2:3,1:2)

	end if

end subroutine
