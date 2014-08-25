!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Aurentz² Vandebril³ Watkins²
!
! ²Dept. Mathematics, Washington State University
! ³Dept. Computer Science, KU Leuven
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Residual check, includes newton steps
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine RESCHECK(BANDSWITCH,N,K,NEWTNUM,POLY,COEFFS,ROOTS,ALLROOTS,RESIDUALS)

	implicit none

	integer, intent(in) :: BANDSWITCH,N,K,newtnum
	complex(kind(1d0)), intent(in) :: poly(N),roots(N),coeffs(N,K)
	double precision, intent(inout) :: residuals(N,3*(newtnum+1))
	complex(kind(1d0)), intent(inout) :: allroots(N,newtnum+1)

	if(BANDSWITCH == 0)then
	
		call COMPRESCHECK(N,newtnum,poly,roots,allroots,residuals)

	else if(BANDSWITCH == 1)then
		
		call CONGRESCHECK(N,K,NEWTNUM,POLY,COEFFS,ROOTS,ALLROOTS,RESIDUALS)

	else

		write(*,*) "Not a valid argument for BANDSWITCH!"
		return
		
	end if


end subroutine

! ********************************************** 
!  August 17, 2012
! ********************************************** 
! 
! This subroutine computes three residuals for each computed root, lambda, 
! of a polynomial P(x). It is also capable of applying an arbitrary number
! of Newton iterations to each computed root. 
!
! The residuals are |P(lambda)/P'(lambda)|, |P(lambda)/P'(lambda)/lambda|, and 
! ||Cv-lambda v||/||C||/||v||, in the inifinity norm where C is the Companion Matrix
! and v is the eigenvectro associated with lambda.
!
! **************************************************************
!
! Input variables:
!
! POLY         complex array of length N, containing the nonleading 
!              coefficients of P (ordered with decreasing N). P must be
!	       monic
!
! N            N of polynomial 
!
!
! ROOTS        complex array of length N, containing the computed roots
!              of Poly
!
! NEWTNUM      a non-negative integer specifying the number of Newton iterations
!              zero is acceptable.
!
! Output variables:
!
! ALLROOTS     complex array of dimension (N,NEWTNUM+1) contains the original roots 
!	       as well as any newton corrections
!
! RESIDUALS    double precision array of dimension (N,3*(NEWTNUM+1)).
!	       each row of RESIDUALS corresponds to one root of POLY.
!	       columns are in sets of three, with the columns 1, 2 and 3 corresponding
!	       to the three residuals mentioned above respectively. Every set of three
!              columns corresponds to a Newton iteration except for the first one.
!
! ***************************************************************

subroutine COMPRESCHECK(N,newtnum,poly,roots,allroots,residuals)

	implicit none

	integer, intent(in) :: N,newtnum
	complex(kind(1d0)), intent(in) :: poly(N),roots(N)
	double precision, intent(inout) :: residuals(N,3*(newtnum+1))
	complex(kind(1d0)), intent(inout) :: allroots(N,newtnum+1)

	integer ii,jj,kk
	double precision :: Cnorm
	complex(kind(1d0)) :: f, fprime, lambda

	! initialize allroots
	allroots = complex(0d0,0d0)
	allroots(:,1) = roots

	! Matrix infinity norms
	Cnorm = 0d0
	do ii=1,N
		Cnorm = Cnorm + abs(poly(ii))
	end do

	Cnorm = dmax1(1d0,Cnorm)

	! Function evaluations and Newton Corrections
	do ii=1,N

		do jj=1,(newtnum+1)
			! Roots inside or on the unit circle
			if(abs(allroots(ii,jj)) <= 1d0)then
				! function evals
				lambda = allroots(ii,jj)
				f = lambda + poly(1)
				fprime = complex(dble(N),0d0)*f - poly(1)
				do kk=2,(N-1)
					f = lambda*f + poly(kk)
					fprime = lambda*fprime + complex(dble(N-kk),0d0)*poly(kk)
				end do
				f = f*lambda + poly(N)
	
				! Store residuals
				residuals(ii,3*(jj-1)+1) = abs(f/fprime)
				residuals(ii,3*(jj-1)+2) = abs(f/fprime/lambda)
				residuals(ii,3*(jj-1)+3) = abs(f)/Cnorm

				! Newton correction
				if((newtnum+1-jj) > 0)then
					lambda = lambda - f/fprime
					allroots(ii,jj+1) = lambda
				end if
			! Roots outside the unit circle
			else
				! function evals
				lambda = complex(1d0,0d0)/allroots(ii,jj)
				f = poly(N)*lambda + poly(N-1)
				fprime = complex(dble(N),0d0)*f - poly(N-1)
				do kk=2,(N-1)
					f = lambda*f + poly(N-kk)
					fprime = lambda*fprime + complex(dble(N-kk),0d0)*poly(N-kk)
				end do
				f = f*lambda + complex(1d0,0d0)
	
				! Store residuals
				residuals(ii,3*(jj-1)+1) = abs(f/fprime*lambda*lambda)
				residuals(ii,3*(jj-1)+2) = abs(f/fprime*lambda)
				residuals(ii,3*(jj-1)+3) = abs(f*lambda)/Cnorm

				! Newton correction
				if((newtnum+1-jj) > 0)then
					lambda = lambda - f/fprime
					allroots(ii,jj+1) = complex(1d0,0d0)/lambda
				end if
			end if
		end do
	end do


end subroutine


! ********************************************** 
!  August 17, 2012
! ********************************************** 
! 
! This subroutine computes three residuals for each computed root, lambda, 
! of a polynomial P(x). It is also capable of applying an arbitrary number
! of Newton iterations to each computed root. 
!
! The residuals are |P(lambda)/P'(lambda)|, |P(lambda)/P'(lambda)/lambda|, and 
! ||Cv-lambda v||/||C||/||v||, in the inifinity norm where C is the Companion Matrix
! and v is the eigenvectro associated with lambda.
!
! **************************************************************
!
! Input variables:
!
! N            degree of polynomial 
!
! K            Bandwidth of upper-triangular matrix 
!
! POLY         complex array of length N, containing the nonleading 
!              coefficients of P (ordered with decreasing N). P must be
!	       monic
!
! COEFFS       complex array of dimension (N,K) containing the recursion
!	       coefficients for the polynomial basis, COEFFS(1,1) must be 1
!
!
! ROOTS        complex array of length N, containing the computed roots
!              of Poly
!
! NEWTNUM      a non-negative integer specifying the number of Newton iterations
!              zero is acceptable.
!
! Output variables:
!
! ALLROOTS     complex array of dimension (N,NEWTNUM+1) contains the original roots 
!	       as well as any newton corrections
!
! RESIDUALS    double precision array of dimension (N,3*(NEWTNUM+1)).
!	       each row of RESIDUALS corresponds to one root of POLY.
!	       columns are in sets of three, with the columns 1, 2 and 3 corresponding
!	       to the three residuals mentioned above respectively. Every set of three
!              columns corresponds to a Newton iteration except for the first one.
!
! ***************************************************************

subroutine CONGRESCHECK(N,K,NEWTNUM,POLY,COEFFS,ROOTS,ALLROOTS,RESIDUALS)

	implicit none

	integer, intent(in) :: N,K,newtnum
	complex(kind(1d0)), intent(in) :: poly(N),roots(N),coeffs(N,K)
	double precision, intent(inout) :: residuals(N,3*(newtnum+1))
	complex(kind(1d0)), intent(inout) :: allroots(N,newtnum+1)

	integer ii,jj,kk,ll,length
	double precision :: Cnorm,temp,Pnorm
	complex(kind(1d0)) :: f, fprime, lambda
	complex(kind(1d0)), allocatable :: P(:),Pprime(:)

	! allocate memory
	allocate(P(N+1),Pprime(N+1))

	! initialize residuals
	residuals = 10d0

	! initialize allroots
	allroots = complex(0d0,0d0)
	allroots(:,1) = roots

	! Matrix infinity norm
	P = poly
	P(1:(K-1)) = P(1:(K-1))- Coeffs(1,2:K)

	Cnorm = 0d0
	do ii=1,N
		Cnorm = Cnorm + abs(P(ii))
	end do

	do ii=2,N
		temp = 0d0
		length = min(K,N+2-ii)
		do jj=1,length
			temp = temp + abs(Coeffs(ii,jj))
		end do

		if(temp > Cnorm)then
			Cnorm = temp
		end if
	end do

	! Function evaluations and Newton Corrections
	do ii=1,N

		do jj=1,(newtnum+1)
			! Roots inside or on the unit circle
			if(abs(allroots(ii,jj)) <= 1d0)then
				! function evals
				lambda = allroots(ii,jj)
			
				! P_k(lambda)
				P(1) = complex(1d0,0d0)

				if(K == 1)then
					P(2) = lambda/Coeffs(N,1)
				else
					P(2) = (lambda - Coeffs(N,2))/Coeffs(N,1)
				end if

				do kk=1,(N-1)
					if(K == 1)then
						P(kk+2) = lambda*P(kk+1)/Coeffs(N-kk,1)
					else
						P(kk+2) = (lambda - Coeffs(N-kk,2))*P(kk+1)

						length=min(kk,K-2)
						do ll=1,length
							P(kk+2) = P(kk+2) - Coeffs(N-kk,2+ll)*P(kk+1-ll)
						end do

						P(kk+2) = P(kk+2)/Coeffs(N-kk,1)
					end if
				end do

				! P'_k(lambda)
				Pprime(1) = complex(0d0,0d0)
	
				Pprime(2) = complex(1d0,0d0)/Coeffs(N,1)

				do kk=1,(N-1)
					if(K == 1)then
						Pprime(kk+2) = (lambda*Pprime(kk+1) + P(kk+1))/Coeffs(N-kk,1)
					else
						Pprime(kk+2) = (lambda - Coeffs(N-kk,2))*Pprime(kk+1) + P(kk+1)
						length=min(kk,K-2)
						do ll=1,length
							Pprime(kk+2) = Pprime(kk+2) - Coeffs(N-kk,2+ll)*Pprime(kk+1-ll)
						end do

						Pprime(kk+2) = Pprime(kk+2)/Coeffs(N-kk,1)

					end if
				end do

				! compute vector norm
				Pnorm = abs(P(1))
				do kk=2,N+1
					temp = abs(P(kk))
					if(temp > Pnorm)then
						Pnorm = temp
					end if
				end do

				!P(lambda) and P'(lambda)
				f = P(N+1)
				fprime = Pprime(N+1)
				do kk=1,N
					f = f + poly(kk)*P(N+1-kk)
					fprime = fprime + poly(kk)*Pprime(N+1-kk)
				end do
	

!write(*,*) "f' =",fprime
!write(*,*) "f =",f
!P(1) = (lambda - Coeffs(1,2))/Coeffs(1,1) + poly(1)
!P(2) = (lambda - Coeffs(2,2))*P(1)/Coeffs(2,1) - Coeffs(1,3)/Coeffs(1,1) + poly(2) 
!do kk = 3,N
!P(kk) = (lambda - Coeffs(kk,2))*P(kk-1)/Coeffs(kk,1) - Coeffs(kk-1,3)*P(kk-2)/Coeffs(kk-1,1) + poly(kk) 
!end do
!write(*,*) "f =",P(N)
!return
				! Store residuals
				residuals(ii,3*(jj-1)+1) = abs(f/fprime)
				residuals(ii,3*(jj-1)+2) = abs(f/fprime/lambda)
				residuals(ii,3*(jj-1)+3) = abs(f)/Cnorm/Pnorm

				! Newton correction
				if((newtnum+1-jj) > 0)then
					lambda = lambda - f/fprime
					allroots(ii,jj+1) = lambda
				end if
			! Roots outside the unit circle
			else
				! function evals
				lambda = allroots(ii,jj)
			
				! P_k(lambda)
				P(1) = complex(1d0,0d0)/lambda
				Pprime(1) = complex(0d0,0d0)

				if(K == 1)then
					P(2) = P(1)/Coeffs(N,1)
					Pprime(2) = complex(1d0,0d0)/Coeffs(N,1)/lambda/lambda
				else
					P(2) = (lambda - Coeffs(N,2))*P(1)/Coeffs(N,1)/lambda
					Pprime(2) = complex(1d0,0d0)/Coeffs(N,1)/lambda/lambda
				end if

				P(1) = P(1)/lambda
				
				do kk=1,(N-1)
					if(K == 1)then
						P(kk+2) = P(kk+1)/Coeffs(N-kk,1)
						Pprime(kk+2) = (lambda*Pprime(kk+1) + P(kk+1))/Coeffs(N-kk,1)
					else
						P(kk+2) = (lambda - Coeffs(N-kk,2))*P(kk+1)/lambda
						Pprime(kk+2) = (lambda - Coeffs(N-kk,2))*Pprime(kk+1) + P(kk+1)

						length=min(kk,K-2)
						do ll=1,length
							P(kk+2) = P(kk+2) - Coeffs(N-kk,2+ll)*P(kk+1-ll)/lambda
							Pprime(kk+2) = Pprime(kk+2) - Coeffs(N-kk,2+ll)*Pprime(kk+1-ll)
						end do

						P(kk+2) = P(kk+2)/Coeffs(N-kk,1)
						Pprime(kk+2) = Pprime(kk+2)/Coeffs(N-kk,1)
						
					end if

					P(1:kk+1) = P(1:kk+1)/lambda
					Pprime(2:kk+2) = Pprime(2:kk+2)/lambda

				end do

				! compute vector norm
				Pnorm = abs(P(1))
				do kk=2,N+1
					temp = abs(P(kk))
					if(temp > Pnorm)then
						Pnorm = temp
					end if
				end do

				!P(lambda) and P'(lambda)
				f = P(N+1)
				fprime = Pprime(N+1)

				do kk=1,N
					f = f + poly(kk)*P(N+1-kk)
					fprime = fprime + poly(kk)*Pprime(N+1-kk)
				end do

				! Store residuals
				residuals(ii,3*(jj-1)+1) = abs(f/fprime)
				residuals(ii,3*(jj-1)+2) = abs(f/fprime/lambda)
				residuals(ii,3*(jj-1)+3) = abs(f)/Cnorm/Pnorm

				! Newton correction
				if((newtnum+1-jj) > 0)then
					lambda = lambda - f/fprime
					allroots(ii,jj+1) = lambda
				end if
			end if
		end do
	end do

	! free memory
	deallocate(P,Pprime)


end subroutine
