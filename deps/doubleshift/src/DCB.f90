!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Aurentz² Mach³ Vandebril³ Watkins²
!
! ²Dept. Mathematics, Washington State University
! ³Dept. Computer Science, KU Leuven
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Last modified 22 August 2014
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! D Chase Bulge
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine chases the bulge.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! N		problem size
!
! str		index for first relevant column of A
!
! stp		index for last relevant block of A
!
! QCB		generators for A
!
! WORK		workspace to compute blocks
!
! B2,B3		generators for the first givens rotations
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine DCB(N,str,stp,QCB,B1,B2,tr)

	implicit none
  
	! input variables
	integer, intent(in) :: N, str, stp
	double precision, intent(inout) :: QCB(6*N), B1(2), B2(2)
	integer, intent(in) :: tr
	  
	! compute variables
	integer :: ii, jj, ind
	double precision :: TEMP(2), B3(2), B1b(2), B2b(2)
  
 	! starting at the top 
	if(str == 1)then
		TEMP(1) = B2(1)
		TEMP(2) = -B2(2)
		B3(1) = B1(1)
		B3(2) = -B1(2)
		call DGTO2(TEMP,B3,QCB(1:2))
		call DFGR(1,B3,QCB(7:8))
		B3 = QCB(1:2)
		QCB(1:2) = TEMP
				
	! otherwise
	else
		ind = 6*(str-1)
		TEMP(1) = B2(1)
		TEMP(2) = -QCB(ind-5)*B2(2)
		B3(1) = B1(1)
		B3(2) = -B1(2)
		call DGTO2(TEMP,B3,QCB((ind+1):(ind+2)))
		call DFGR(1,B3,QCB((ind+7):(ind+8)))
		B3 = QCB((ind+1):(ind+2))
		QCB((ind+1):(ind+2)) = TEMP
	end if
  
	! main chasing loop
	do ii=str,(stp-2)
		
		! set ind
		ind = 6*(ii-1)
  
  
                if (ii<tr-1) then
                   ! B and C* are equal 
                   B1b = B1
                   B2b = B2
                   ! through B
                   call DGTO2(QCB((ind+11):(ind+12)),QCB((ind+17):(ind+18)),B1b)
                   call DGTO2(QCB((ind+5):(ind+6)),QCB((ind+11):(ind+12)),B2b)
                   QCB(ind+3)  =  QCB(ind+5)
                   QCB(ind+4)  = -QCB(ind+6) 
                   QCB(ind+9)  =  QCB(ind+11)
                   QCB(ind+10) = -QCB(ind+12)
                   QCB(ind+15) =  QCB(ind+17)
                   QCB(ind+16) = -QCB(ind+18)

                else                 
                   ! through B
                   call DGTO2(QCB((ind+11):(ind+12)),QCB((ind+17):(ind+18)),B1)
                   call DGTO2(QCB((ind+5):(ind+6)),QCB((ind+11):(ind+12)),B2)
                   
                   ! through C
                   call DGTO2(QCB((ind+15):(ind+16)),QCB((ind+9):(ind+10)),B1)
                   call DGTO2(QCB((ind+9):(ind+10)),QCB((ind+3):(ind+4)),B2)
                end if

		! through Q
		call DGTO2(QCB((ind+7):(ind+8)),QCB((ind+13):(ind+14)),B1)
		call DGTO2(QCB((ind+1):(ind+2)),QCB((ind+7):(ind+8)),B2)
		
		! push B3 down
		call DGTO2(B3,B1,B2)
		TEMP = B2
		B2 = B3
		B3 = B1
		B1 = TEMP
	end do
	
	! set ind
	ind = 6*(stp-2)
  
	! through B
	call DGTO2(QCB((ind+11):(ind+12)),QCB((ind+17):(ind+18)),B1)
	call DGTO2(QCB((ind+5):(ind+6)),QCB((ind+11):(ind+12)),B2)

	! through C
	call DGTO2(QCB((ind+15):(ind+16)),QCB((ind+9):(ind+10)),B1)
	call DGTO2(QCB((ind+9):(ind+10)),QCB((ind+3):(ind+4)),B2)

	! through Q
	B1(2) = QCB(6*stp+1)*B1(2)
	call DFGR(0,QCB((ind+7):(ind+8)),B1)
	call DGTO2(QCB((ind+1):(ind+2)),QCB((ind+7):(ind+8)),B2)
	call DFGR(0,B3,B2)
	
	! last bulge
	ind = 6*(stp-1)
	call DGTO2(QCB((ind+5):(ind+6)),QCB((ind+11):(ind+12)),B3)
	call DGTO2(QCB((ind+9):(ind+10)),QCB((ind+3):(ind+4)),B3)
	B3(2) = QCB(6*stp+1)*B3(2)
	call DFGR(0,QCB((ind+1):(ind+2)),B3)

end subroutine

