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
! zamvw runs Francis's implicity shifted QR 
! algorithm on a factored form of the companion
! matrix, the rank structure in the upper 
! triangular is exploited
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! n         problem size
!
! Q,D,C,B   generators of A
!
! reigs     (out) real parts of eigenvalues
! ieigs     (out) imag parts of eigenvalues
!
! its       (out) array, number of iteration
!           between two subsequent deflations
!
! flag      error flag
!           0    no error, all eigenvalues found
!           k>0  QR algorithm did not converge,
!                k eigenvalues are found (first k 
!                entries of reigs,ieigs)
!
! tr        first tr rotations in B and C* are equal
!
! rayleigh  0 (default) Wilkinson shift
!           1 Rayleigh shift
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine zamvw2(n,Q,D,C,B,reigs,ieigs,its,flag,tr,rayleigh)

	implicit none

	! input variables
	integer, intent(in) :: n
	double precision, intent(inout) :: Q(3*n), D(2*n+2), C(3*n), B(3*n)
	double precision, intent(inout) :: reigs(n), ieigs(n)
	integer, intent(inout) :: its(n), flag
        integer, optional :: rayleigh
	integer, intent(inout) :: tr

	! compute variables
	integer :: ii, jj, kk, ind1, ind2, ll, strt, k, ray, bc
	integer :: start_index, stop_index, zero_index, it_max, it_count
	double precision :: tol, nrm
	complex(kind(1d0)) :: shift, block(2,2), e1, e2
	double precision :: bulge(3), s1, s2

        flag = 0
        bc = 0 

        ! rayleigh: 0 for Wilkinson shift, 1 for Rayleigh shift
        if (present(rayleigh)) then
           if (rayleigh .NE. 0) then
              ray = 1
           else
              ray = 0
           end if
        else
           ray = 0
        end if

	! set tol	
	tol = epsilon(1d0)

	! check to make sure it's worth the effort
	if(n <= 2)then
		print*, "n =", n
		print*, "n should be atleast 3 to use this algorithm!"
		stop
	end if   

	! initialize storage
	its = 0

	! initialize indices
	start_index = 1
	stop_index = n-1
	zero_index = 0
	it_max = 30*n
	it_count = 0

	! loop for bulgechasing
	do kk=1,it_max
		! check for completion
		if(stop_index <= 0)then
			!print*, "Algorithm is complete!"
			exit
		end if

		! check for deflation
		call deflation(n,start_index,stop_index,zero_index,Q,D,C,B,its,it_count)
		
		! if 1x1 block remove and check again 
		if(stop_index == zero_index)then
			! get 2x2 block
			call diagblock(n,stop_index,block,Q,D,C,B)
	
			! zero at top
			if(stop_index == 1)then
				! store the eigenvalue
				reigs(stop_index) = dble(block(1,1))
				ieigs(stop_index) = dimag(block(1,1))
				reigs(stop_index+1) = dble(block(2,2))
				ieigs(stop_index+1) = dimag(block(2,2))

				! update stop_index
				stop_index = 0

			! anywhere else
			else
				! store the eigenvalue
				reigs(stop_index+1) = dble(block(2,2))
				ieigs(stop_index+1) = dimag(block(2,2))

				! update indices
				stop_index = stop_index - 1
				zero_index = 0
				start_index = 1

			end if
			

		! if 2x2 block remove and check again
		else if(stop_index-1 == zero_index)then
			! get 2x2 block
			call diagblock(n,stop_index,block,Q,D,C,B) 
			
			! zero at top
			if(stop_index == 2)then
				! store the eigenvalues
				call modified_quadratic(block,e1,e2)
				reigs(stop_index) = dble(e1)
				ieigs(stop_index) = dimag(e1)
				reigs(stop_index+1) = dble(e2)
				ieigs(stop_index+1) = dimag(e2)
				call diagblock(n,1,block,Q,D,C,B)
				reigs(1) = dble(block(1,1))
				ieigs(1) = dimag(block(1,1))

				! update indices
				stop_index = 0
				!zero_index = 0
				!start_index = 0 

			! otherwise
			else
				! store the eigenvalues
				call modified_quadratic(block,e1,e2)
				reigs(stop_index) = dble(e1)
				ieigs(stop_index) = dimag(e1)
				reigs(stop_index+1) = dble(e2)
				ieigs(stop_index+1) = dimag(e2)

				! update indices
				stop_index = stop_index - 2
				zero_index = 0
				start_index = 1 

			end if

		! if greater than 2x2 chase a bulge and check again
		else

			! it_count
			it_count = it_count + 1	

			! compute first transformation
			if(kk == 1)then
				call normalpoly(1,s1,s2)
                                if (ray == 0) then                                                              
                                   shift = complex(s1,s2)
                                else
                                   shift = complex(s1,0d0)
                                end if
			else if(mod(it_count,15) == 0)then
				call normalpoly(1,s1,s2)
                                if (ray == 0) then                                                              
                                   shift = complex(s1,s2)
                                else
                                   shift = complex(s1,0d0)
                                end if
				print*, "Random shift!", shift
			else
				call diagblock(n,stop_index,block,Q,D,C,B)

                                if (ray == 0) then              
                                   call modified_quadratic(block,e1,e2)
                                   if(zabs(block(2,2)-e1) < zabs(block(2,2)-e2))then
                                      shift = e1
                                   else
                                      shift = e2
                                   end if
                                else
                                   shift = block(2,2)
                                end if
			end if

			! build bulge
			call buildbulge(n,start_index,bulge,shift,Q,D,C,B)

			! chase bulge
			call chasebulge(n,start_index,stop_index,bulge,Q,D,C,B,tr)
                        bc = bc + 1
                        tr = tr - 1

		end if
	end do
 
   if (kk>=it_max-1) then
     if (stop_index < N-1) then
        ! there some found eigenvalues, but not all have been found
        FLAG = N - 1 - stop_index
        print*, "QR algorithm did not converged within 30*N&
             & iterations, although FLAG = ", FLAG ,&
             & "eigenvalues have been found. This is a very rare case."
        print*, "Try to increase it_max &
             & or consider a bug-report to email:&
             & thomas.mach+zamvw.bugreport@gmail.com."
        do ii=1,FLAG
           reigs(ii) = reigs(stop_index+1+ii)
           ieigs(ii) = ieigs(stop_index+1+ii)
           reigs(stop_index+1+ii) = 0d0
           ieigs(stop_index+1+ii) = 0d0
        end do
     end if
  end if
  !print*, "Total number of bulgechases: ", bc 
end subroutine
