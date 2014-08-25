!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Aurentz² Mach³ Vandebril³ Watkins²
!
! ²Dept. Mathematics, Washington State University
! ³Dept. Computer Science, KU Leuven
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Last modified 22 August 2014
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! D Aurentz Mach Vandebril Watkins
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Real Doubleshift Code
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine computes the eigenvalues of the companion matrix for P(x),
!
! P(x) = x^N + a_N-1 x^N-1 + ... + a_1 x + a_0,
! 
! using a variant of Francis' real, doubleshift algorithm that 
! exploits the rank-structure in the upper triangular part.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! N		degree of the polynomial
!
! POLY		array containing coefficients of P(x),
! 		POLY = [a_N-1, ... , a_0]
!
! REIGS		array for real part of eigenvalues
!
! IEIGS		array for imaginary part of eigenvalues
!
! ITS		array for iteration counts
!
! FLAG		flag for errors
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine DAMVW(NP,POLY,REIGS,IEIGS,ITS,FLAG)
  
  implicit none
  
  ! input variables
  integer, intent(in) :: NP
  double precision, intent(in) :: POLY(NP)
  double precision, intent(inout) :: REIGS(NP), IEIGS(NP)
  integer, intent(inout) :: ITS(NP), FLAG
  
  ! compute variables
  integer :: ii, jj, kk, N, strt, nnew, tr
  integer :: start_index, stop_index, zero_index, it_max, it_count, chase_count
  double precision :: tol = 1d-16
  double precision :: scrap, num, temp(3,2), B1(2), B2(2), re1, ie1, re2, ie2
  double precision :: Q1(2), Q2(2), Q3(2)
  double precision, allocatable :: QCB(:)
  double precision :: ALPHA, nrm, trace, disc, detm
  
  FLAG = 0

  ! check to make sure it's worth the effort
  if(NP <= 0)then
     print*, "N =", NP
     print*, "N should be at least 3 to use this algorithm!"
     FLAG = -1
     return
  end if

  REIGS = 0d0
  IEIGS = 0d0

  ! check for deflatable zero coefficients 
  N = NP
  nnew = 0
  do ii=1,n
     nrm = abs(POLY(N+1-ii))
     if ( nrm /= 0) then
        nnew = n+1 - ii
        exit
     end if
  end do
  
  N = nnew
  print*, N
  !print*, POLY
  if (N == 0) then
     ! all coefficients 0 => all roots 0
     FLAG = 0 
     return
  end if
  if (N == 1) then
     ! it remains a polynomial of degree 1 
     reigs(1) = -POLY(1)
     FLAG = 0
     return
  end if
  if (N == 2) then
     ! it remains a polynomial of degree 2
     ! use modified quadratic formula
     ! compute intermediate values
     trace = -POLY(1)
     detm = POLY(2)
     disc = trace*trace - 4d0*detm
     
     ! compute e1 and e2
     ! complex eigenvalues
     if(disc < 0)then
        reigs(1) = trace/2d0
        ieigs(1) = sqrt(-disc)/2d0
        reigs(2) = reigs(1)
        ieigs(2) = -ieigs(1)
	! real eignevalues
     else if(abs(trace+sqrt(disc)) > abs(trace-sqrt(disc)))then
        if(abs(trace+sqrt(disc)) == 0)then
           reigs(1) = 0d0
           ieigs(1) = 0d0
           reigs(2) = 0d0
           ieigs(2) = 0d0
        else
           reigs(1) = (trace+sqrt(disc))/2d0
           ieigs(1) = 0d0
           reigs(2) = detm/reigs(1)
           ieigs(2) = 0d0
        end if
     else
        if(abs(trace-sqrt(disc)) == 0)then
           reigs(1) = 0d0
           ieigs(1) = 0d0
           reigs(2) = 0d0
           ieigs(2) = 0d0
        else
           reigs(1) = (trace-sqrt(disc))/2d0
           ieigs(1) = 0d0
           reigs(2) = detm/reigs(1)
           ieigs(2) = 0d0
        end if
     end if
     FLAG = 0
     return
  end if

  ! remaining polynomial has a degree larger than 2
  
  ! allocate memory
  allocate(QCB(6*N))
  
  ! factor column companion matrix
  call DFCC(N,POLY,QCB,ALPHA)
  tr = n-2
  
  ! initialize storage
  ITS = 0
  REIGS = 0d0
  IEIGS = 0d0   
  
  ! initialize indices
  start_index = 1
  stop_index = N-1
  zero_index = 0
  it_max = 30*N
  it_count = 0
  chase_count = 0

  ! loop for bulge chasing
  do kk=1,it_max

     ! check for completion
     if(stop_index <= 0)then
        !print*, "Algorithm is complete!"
        exit
     end if
       
     
     ! check for deflation
     call DCFD(N,start_index,stop_index,zero_index,QCB,its,it_count)
     
     ! if 1x1 block remove and check again 
     if(stop_index == zero_index)then
        ! get 2x2 block
        call DCDB(N,stop_index,TEMP,QCB) 
        
        ! zero at top
        if(stop_index == 1)then
           ! store the eigenvalue
           REIGS(stop_index) = TEMP(1,1)
           REIGS(stop_index+1) = TEMP(2,2)
           
           ! update stop_index
           stop_index = 0
           
        ! anywhere else
        else
           ! store the eigenvalue
           REIGS(stop_index+1) = TEMP(2,2)
           
           ! update indices
           stop_index = stop_index - 1
           zero_index = 0
           start_index = 1
           
        end if
        
		! if 2x2 block remove and check again
		else if(stop_index-1 == zero_index)then
			! get 2x2 block
			call DCDB(N,stop_index,TEMP,QCB) 
        
		    ! zero at top
		    if(stop_index == 2)then
		       ! store the eigenvalues
		       call DMQF(TEMP(1:2,:),REIGS(stop_index),IEIGS(stop_index),REIGS(stop_index+1),IEIGS(stop_index+1))
		       call DCDB(N,1,TEMP,QCB) 
		       REIGS(1) = TEMP(1,1)
		       
		       ! update indices
		       stop_index = stop_index - 2
		       zero_index = 0
		       start_index = 1 
		       
		       ! otherwise
		    else
		       ! store the eigenvalues
		       call DMQF(TEMP(1:2,:),REIGS(stop_index),IEIGS(stop_index),REIGS(stop_index+1),IEIGS(stop_index+1))
		       
		       ! update indices
		       stop_index = stop_index - 2
		       zero_index = 0
		       start_index = 1 
		       
		    end if
        
		! if greater than 2x2 chase a bulge and check again
		else
        
		    ! it_count
		    it_count = it_count + 1	
		    
		    ! compute shifts
		    if(kk == 1) then          
		       call DCDB(N,stop_index,TEMP,QCB) 
		       call DMQF(TEMP(1:2,:),re1,ie1,re2,ie2)
		    elseif (mod(it_count,15) == 0)then
				call dnormalpoly(1,re1)
				call dnormalpoly(1,ie1)
				re2 = re1
				ie2 = -ie1
				print*, "Random shift!"
		    else
		       call DCDB(N,stop_index,TEMP,QCB) 
		       call DMQF(TEMP(1:2,:),re1,ie1,re2,ie2)
		    end if

		    ! build bulge
		    call DCFT(N,start_index,QCB,re1,ie1,re2,ie2,B1,B2)
	       
			! chase bulge
		    chase_count = chase_count + 1
		    call DCB(N,start_index,stop_index,QCB,B1,B2,tr)
                    tr = tr - 2
		end if
	end do

  !print*, chase_count

  if (kk>=it_max-1) then
     if (stop_index < N-1) then
        ! some eigenvalues have been found, but not all of them
        ! this is a rare case
        FLAG = N - 1 - stop_index
        print*, "QR algorithm did not converged within 30*N&
             & iterations, although FLAG = ", FLAG ,&
             & "eigenvalues have been found.&
             & This is a very rare case."
        print*, "Try to increase it_max &
             & or consider a bug-report to email:&
             & thomas.mach+damvw.bugreport@gmail.com."
        do ii=1,FLAG
           reigs(ii) = reigs(stop_index+1+ii)
           ieigs(ii) = ieigs(stop_index+1+ii)
           reigs(stop_index+1+ii) = 0d0
           ieigs(stop_index+1+ii) = 0d0
        end do
     end if
     ! debugging

     !print*, kk
     !print*, it_max
     !print*, start_index, stop_index
     !print*, reigs
     !print*, ieigs

     !do ii=1,N-1
     !print*, ""
     !   print*, poly(ii)
     !end do
     !print*, ""
  end if
  
  ! free memory
  deallocate(QCB)
  
  
end subroutine
