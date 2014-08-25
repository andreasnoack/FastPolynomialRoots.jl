!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Aurentz² Mach³ Vandebril³ Watkins²
!
! ²Dept. Mathematics, Washington State University
! ³Dept. Computer Science, KU Leuven
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Last modified 22 August 2014
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Computes Givens rotation turnover
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! input/output
! G1, G2, B     rotations given as
!               cr + i*ci, s, with cr, ci, s \in R
!               
!
! G1    B         G1
! G1 G2 B  =>  B  G1 G2
!    G2        B     G2
!
! or
!    G2        B     G2
! G1 G2 B  =>  B  G1 G2
! G1    B         G1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dto4(G1,G2,B)

	implicit none

	! input variables
	double precision, intent(inout) :: G1(3)
	double precision, intent(inout) :: G2(3)
	double precision, intent(inout) :: B(3)

	! compute variables
	double precision :: nrm, tol, dnrm2, T(3)
	
	double precision :: c1r
	double precision :: c1i
	double precision :: s1
	double precision :: c2r
	double precision :: c2i
	double precision :: s2
	double precision :: c3r
	double precision :: c3i
	double precision :: s3
	
	double precision :: c4r
	double precision :: c4i
	double precision :: s4
	double precision :: c5r
	double precision :: c5i
	double precision :: s5
	double precision :: c6r
	double precision :: c6i
	double precision :: s6
	
	! set tol
	tol = epsilon(1d0)
        !tol = 3e-14
	
	! set local variables
	c1r = G1(1)
	c1i = G1(2)
	s1 = G1(3)
	c2r = G2(1)
	c2i = G2(2)
	s2 = G2(3)
	c3r = B(1)
	c3i = B(2)
	s3 = B(3)	

	! initialize c4r, c4i and s4
	T(1) = s1*c3r + (c1r*c2r + c1i*c2i)*s3
	T(2) = s1*c3i + (-c1i*c2r + c1r*c2i)*s3
	T(3) = s2*s3
        nrm = T(1)*T(1) + T(2)*T(2) + T(3)*T(3)
        if (dabs(nrm-1d0)<tol) then
           c4r = T(1)
           c4i = T(2)
           s4  = T(3)
           nrm = 1d0
        else
           call rot3(T(1),T(2),T(3),c4r,c4i,s4,nrm)
        end if
	! initialize c5r, c5i and s5
	T(1) = c1r*c3r - c1i*c3i - s1*c2r*s3
	T(2) = c1r*c3i + c1i*c3r - s1*c2i*s3
	T(3) = nrm
        nrm = T(1)*T(1) + T(2)*T(2) + T(3)*T(3)
        if (dabs(nrm-1d0)<tol) then
           c5r = T(1)
           c5i = T(2)
           s5  = T(3)
        else
           call rot4(T(1),T(2),T(3),c5r,c5i,s5)	
        end if
	
	! initialize c6r, c6i and s6
	T(1) = c2r*c4r + c2i*c4i + c1r*s2*s4
	T(2) = c2i*c4r - c2r*c4i + c1i*s2*s4
	T(3) = s1*s2*s5 - (c5r*c2r + c5i*c2i)*s4 + s2*(c1r*(c4r*c5r + c4i*c5i) + c1i*(c4r*c5i - c4i*c5r))
        nrm = T(1)*T(1) + T(2)*T(2) + T(3)*T(3)
        if (dabs(nrm-1d0)<tol) then
           c6r = T(1)
           c6i = T(2)
           s6  = T(3)
        else
           call rot4(T(1),T(2),T(3),c6r,c6i,s6)		
        end if

	! store output
	G1(1) = c5r
	G1(2) = c5i
	G1(3) = s5
	G2(1) = c6r
	G2(2) = c6i
	G2(3) = s6
	B(1) = c4r
	B(2) = c4i
	B(3) = s4
	
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute Givens rotation zeroing b
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! G1 [ ar + i*ai ] = [ nrm ]
! G1 [    b      ] = [     ]
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! all variables real
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rot3(ar,ai,b,cr,ci,s,nrm)

	implicit none
	
	! input variables
	double precision, intent(inout) :: ar, ai, b, cr, ci, s, nrm
	
	! compute variables
	double precision :: nar, nai, nb, tol
	
	! set tol
	tol = epsilon(1d0)

	! set local variables
	nar = abs(ar)
	nai = abs(ai)
	nb = abs(b)
	
	if(nar == 0 .AND. nai == 0 .AND. nb == 0)then
		cr = 1d0
		ci = 0d0
		s = 0d0
		nrm = 0d0
	else if(nb == 0 .AND. nar > nai)then
		ai = ai/ar
		nrm = sqrt(1d0 + ai*ai)
		if(ar < 0)then
			cr = -1d0/nrm
			ci = ai*cr
			s = 0d0
			nrm = -ar*nrm
		else
			cr = 1d0/nrm
			ci = ai*cr
			s = 0d0
			nrm = ar*nrm	
		end if
	else if(nb == 0)then
		ar = ar/ai
		nrm = sqrt(1d0 + ar*ar)
		if(ai < 0)then
			ci = -1d0/nrm
			cr = ar*ci
			s = 0d0
			nrm = -ai*nrm
		else
			ci = 1d0/nrm
			cr = ar*ci
			s = 0d0
			nrm = ai*nrm
		end if
	else if(nar >= nb .AND. nar >= nai)then
		b = b/ar
		ai = ai/ar
                nrm = sqrt(1d0 + b*b + ai*ai)
		if(ar < 0)then
			cr = -1d0/nrm
			ci = ai*cr
			s = b*cr
			nrm = -ar*nrm
		else
			cr = 1d0/nrm
			ci = ai*cr
			s = b*cr
			nrm = ar*nrm	
		end if
	else if(nai >= nb .AND. nai >= nar)then
		b = b/ai
		ar = ar/ai
		nrm = sqrt(1d0 + b*b + ar*ar)
		if(ai < 0)then
			ci = -1d0/nrm
			cr = ar*ci
			s = b*ci
			nrm = -ai*nrm
		else
			ci = 1d0/nrm
			cr = ar*ci
			s = b*ci
			nrm = ai*nrm	
		end if
	else
		ar = ar/b
		ai = ai/b
		nrm = sqrt(1d0 + ai*ai + ar*ar)
		if(b < 0)then
			s = -1d0/nrm
			cr = ar*s
			ci = ai*s
			nrm = -b*nrm
		else
			s = 1d0/nrm
			cr = ar*s
			ci = ai*s
			nrm = b*nrm
		end if	
	end if


end subroutine rot3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute Givens rotation zeroing b
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! identical with rot3, but nrm is not computed
! => faster
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rot4(ar,ai,b,cr,ci,s)

        implicit none
	
	! input variables
	double precision, intent(inout) :: ar, ai, b, cr, ci, s
 
        ! compute variables
        double precision :: nar, nai, nb, nrm, tol

	! set tol
	tol = epsilon(1d0)
	
	! set local variables
	nar = abs(ar)
        nai = abs(ai)
	nb = abs(b)
	
	if(nar == 0 .AND. nai == 0 .AND. nb == 0)then
		cr = 1d0
		ci = 0d0
		s = 0d0
		nrm = 0d0
	else if(nb == 0 .AND. nar > nai)then
		ai = ai/ar
		cr = 1d0/sqrt(1d0 + ai*ai)
		if(ar < 0)then
			cr = -cr
		end if
                ci = ai*cr
		s = 0d0
	else if(nb == 0)then
		ar = ar/ai
		ci = 1d0/sqrt(1d0 + ar*ar)
		if(ai < 0)then
			ci = -ci
		end if
                cr = ar*ci
		s = 0d0
	else if(nar >= nb .AND. nar >= nai)then
		b = b/ar
		ai = ai/ar
                cr = 1d0/sqrt(1d0 + b*b + ai*ai)
                if(ar < 0)then
                   cr = -cr
                !   cr = -1d0/sqrt(1d0 + b*b + ai*ai)
                !else 
                !   cr = 1d0/sqrt(1d0 + b*b + ai*ai)
		end if
                ci = ai*cr
                s = b*cr
	else if(nai >= nb .AND. nai >= nar)then
		b = b/ai
		ar = ar/ai
		ci = 1d0/sqrt(1d0 + b*b + ar*ar)
		if(ai < 0)then
			ci = -ci
		end if
 		cr = ar*ci
                s = b*ci
             else 
		ar = ar/b
		ai = ai/b
		s = 1.d0/sqrt(1d0 + ai*ai + ar*ar)
		if(b < 0)then
                   s = -s
                end if
                cr = ar*s
                ci = ai*s
	end if


end subroutine rot4
