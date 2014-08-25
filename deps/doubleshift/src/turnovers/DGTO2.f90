!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Aurentz² Mach³ Vandebril³ Watkins²
!
! ²Dept. Mathematics, Washington State University
! ³Dept. Computer Science, KU Leuven
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Last modified 22 August 2014
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! D Givens Turn Over
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine passes the givens rotation Q3 from the left down 
! through Q1 and Q2.
!
! It overwrites Q1, Q2 and Q3.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Q1, Q2, Q3	generators for two givens rotations
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine DGTO2(Q1,Q2,Q3)

	implicit none
  
	! input variables
	double precision, intent(inout) :: Q1(2), Q2(2), Q3(2)

	! compute variables
	double precision :: tol, nrm, T(3), dnrm2
	double precision :: a, b 
	double precision :: c1, s1
	double precision :: c2, s2
	double precision :: c3, s3
	double precision :: c4, s4
	double precision :: c5, s5
	double precision :: c6, s6
	
	! set tol
	tol = epsilon(1d0)

	! set local variables
	c1 = Q1(1)
	s1 = Q1(2)
	c2 = Q2(1)
	s2 = Q2(2)
	c3 = Q3(1)
	s3 = Q3(2)
	
	! initialize c4 and s4
	a = s1*c3 + c1*c2*s3
	b = s2*s3
        !nrm = a*a + b*b
        !if (abs(nrm-1d0)<tol) then
        !   c4 = a
        !   s4 = b
        !   nrm = 1d0
        !else
        call rot1(a,b,c4,s4,nrm)
        !end if
	
	! initialize c5 and s5
	a = c1*c3 - s1*c2*s3
	b = nrm
        !nrm = a*a + b*b
        !if (abs(nrm-1d0)<tol) then
        !   c5 = a
        !   s5 = b
        !else
        call rot2(a,b,c5,s5)
        !end if
	
	! second column
        ! T(1) = -c1*s3 - s1*c2*c3
        ! T(2) = -s1*s3 + c1*c2*c3
        ! T(3) = s2*c3

	! update second column
        ! b = -T(2)*s4 + T(3)*c4
        ! T(2) = T(2)*c4 + T(3)*s4
        ! a = -T(1)*s5 + T(2)*c5

	! third column
	T(1) = s1*s2
	T(2) = -c1*s2
	T(3) = c2

	! update third column
	a = -T(2)*s4 + T(3)*c4
	T(2) = T(2)*c4 + T(3)*s4
	b = -(-T(1)*s5 + T(2)*c5)
	
	! initialize c6 and s6
        !nrm = a*a + b*b
        !if (abs(nrm-1d0)<tol) then
        !   c6 = a
        !   s6 = b
        !else
        call rot2(a,b,c6,s6)	
        !end if

	! set output
	Q1(1) = c5
	Q1(2) = s5
	Q2(1) = c6
	Q2(2) = s6
	Q3(1) = c4
	Q3(2) = s4
	
     
end subroutine DGTO2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! D Givens Rotation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine computes a givens rotation G1 zeroing b
!
!   [ c -s ] [ a ] = [ r ]
!   [ s  c ] [ b ] = [ 0 ]
! 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! in
! a    scalar  
! b    scalar
!
! out
! c    cosine of rotation
! s    sine of rotation
! r    norm of [ a; b ]
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rot1(a,b,c,s,r)

  implicit none
  
  ! input variables
  double precision, intent(in) :: a,b
  double precision, intent(inout) :: c,s,r
  
	if (b == 0 .AND. a < 0) then
		c = -1d0
        s = 0d0
        r = -a
        
	else if (b == 0) then
        c = 1d0
        s = 0d0
        r = a
        
	else if (abs(a) >= abs(b)) then

		s = b/a
		r = sqrt(1.d0 + s*s)

		if (a<0) then
			c= -1.d0/r
			s= s*c
			r=-a*r
		else
			c = 1.d0/r
			s = s*c
			r = a*r
		end if

	else

		c = a/b;
		r = sqrt(1.d0 + c*c)

		if (b<0) then
			s =-1.d0/r
			c = c*s
			r =-b*r
		else
			s = 1.d0/r
			c = c*s
			r = b*r
		end if

  end if

end subroutine rot1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! D Givens Rotation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine computes a givens rotation G1 zeroing b
!
!   [ c -s ] [ a ] = [ r ]
!   [ s  c ] [ b ] = [ 0 ]
! 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! in
! a    scalar  
! b    scalar
!
! out
! c    cosine of rotation
! s    sine of rotation
! 
! Remark: Faster than rot1 since r is not computed.
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rot2(a,b,c,s)

	implicit none
  
	! input variables
	double precision, intent(in) :: a,b
	double precision, intent(inout) :: c,s
  
	! compute variables
	double precision :: r
  
	if (b == 0 .AND. a < 0) then
		c = -1d0
        s = 0d0
        
	else if (b == 0) then
        c = 1d0
        s = 0d0
        
	else if (abs(a) >= abs(b)) then

		s = b/a
		r = sqrt(1.d0 + s*s)

		if (a<0) then
			c= -1.d0/r
			s= s*c
		else
			c = 1.d0/r
			s = s*c
		end if

	else

		c = a/b;
		r = sqrt(1.d0 + c*c)

		if (b<0) then
			s =-1.d0/r
			c = c*s
		else
			s = 1.d0/r
			c = c*s
		end if

  end if

end subroutine rot2
