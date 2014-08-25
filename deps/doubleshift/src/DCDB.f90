!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Aurentz² Mach³ Vandebril³ Watkins²
!
! ²Dept. Mathematics, Washington State University
! ³Dept. Computer Science, KU Leuven
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Last modified 22 August 2014
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! D Compute Diagonal Blocks
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine computes A(K:K+2,K:K+1)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! N		size of problem
!
! K		desired block index
!
! A		storage location for A(K:K+1,K:K+1)
!
! QCB		generators for A
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine DCDB(N,K,A,QCB)
  
  implicit none
  
  ! input variables
  integer, intent(in) :: N, K
  double precision, intent(in) :: QCB(6*N)
  double precision, intent(inout) :: A(3,2)
  
  ! compute variables
  integer :: ii, ind
  double precision :: temp, R(3,2)
  
	! initialize
	A = 0d0
	R = 0d0
  
	! case K = 1
	if(K == 1)then

		! first column of R
		R(1,1) = -QCB(6)/QCB(4)
		
		! second column of R
		R(2,2) = -QCB(12)/QCB(10)
		R(1,2) = -(QCB(5)*QCB(11) - R(2,2)*QCB(3)*QCB(9))/QCB(4)
     
		! A
		! index for first Q
		ind = 6*(k)
		R(2,2) = QCB(ind+1)*R(2,2)
		
		! index for second Q
		ind = 6*(k-1)
		A(1,1) = QCB(ind+1)
		A(2,1) = QCB(ind+2)
		A(1,2) = -QCB(ind+2)
		A(2,2) = QCB(ind+1)
		R(1:2,:) = matmul(A(1:2,1:2),R(1:2,:))
		

		
		A(1:2,:) = R(1:2,:)
     
	! other cases
	else

		! first column of R
		ind = 6*(k-1)
		R(2,1) = -QCB(ind+6)/QCB(ind+4)
		R(1,1) = -(QCB(ind-1)*QCB(ind+5) - R(2,1)*QCB(ind-3)*QCB(ind+3))/QCB(ind-2)
		
		! second column of R
		ind = 6*(k)
		R(3,2) = -QCB(ind+6)/QCB(ind+4)
		R(2,2) = -(QCB(ind-1)*QCB(ind+5) - R(3,2)*QCB(ind-3)*QCB(ind+3))/QCB(ind-2)
		R(1,2) = (QCB(ind-7)*QCB(ind)*QCB(ind+5) - QCB(ind-9)*(QCB(ind-3)*QCB(ind-1)*QCB(ind+5) &
			- QCB(ind+3)*R(3,2))/QCB(ind-2))/QCB(ind-8)
     
		! A = QR
		! index for first Q
		ind = 6*(k)
		R(3,2) = QCB(ind+1)*R(3,2)
		
		! index for second Q
		ind = 6*(k-1)
		A(1,1) = QCB(ind+1)
		A(2,1) = QCB(ind+2)
		A(1,2) = -QCB(ind+2)
		A(2,2) = QCB(ind+1)
		R(2:3,:) = matmul(A(1:2,1:2),R(2:3,:))
		
		! index for third Q
		ind = 6*(k-2)
		A(1,1) = QCB(ind+1)
		A(2,1) = QCB(ind+2)
		A(1,2) = -QCB(ind+2)
		A(2,2) = QCB(ind+1)
		R(1:2,:) = matmul(A(1:2,1:2),R(1:2,:))
     
		A(1:2,:) = R(2:3,:)
     
	end if
  
end subroutine DCDB
