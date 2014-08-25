!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Aurentz² Mach³ Vandebril³ Watkins²
!
! ²Dept. Mathematics, Washington State University
! ³Dept. Computer Science, KU Leuven
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Last modified 22 August 2014
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! D Check For Deflation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! N		size of problem
!
! strt		start index of current block
!
! stp           stop index of current block
!
! zero          last zero above the current block
!
! QCB		generators for A
!
! its           number of iterations 
!
! itcnt         iteration counter
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine DCFD(N,strt,stp,zero,QCB,its,itcnt)

	implicit none
	
	! input variables
	integer, intent(in) :: N
	integer, intent(inout) :: strt, stp , zero, its(N), itcnt
	double precision, intent(inout) :: QCB(6*N)
  
	! compute variables
	integer :: ii,ind
	double precision :: tol
	
	! set tol
	tol = epsilon(1d0)

	! loop for deflation
	do ii=1,stp
		ind = 6*(stp-ii)
		if(abs(QCB(ind+2)) < tol)then
			! set sub-diagonal to 0
			QCB(ind+2) = 0d0
			QCB(ind+1) = QCB(ind+1)/abs(QCB(ind+1))
           
			! update indices
			zero = stp+1-ii
			strt = zero + 1
           
			! store it_count
			ITS(zero) = itcnt
			itcnt = 0

			exit
			
		end if
		
	end do

end subroutine
