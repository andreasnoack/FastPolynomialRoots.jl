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
! Rootrace 
!
! comparing DAMVW with LAPACK's DHSEQR
! polynomials with random coefficients
!
!
! a * exp(b)
! a random uniform in [-1,1]
! b random uniform in [-R,R]
!
! R given as argument, default R=5
!
!!!!!!!!!!!!!!!!!!!!
!
! Random polynomials chosen as (iv) in 
! Jenkins, Traub 1970
!
! "(iv) polynomials whose coefficients are chosen
!  randomly by taking the mantissa and exponents
!  from seprate uniform distributions. The resulting
!  polynomials have widely varying zeros and hence
!  yield a reasonable test that the program has wide
!  applications."
!
! [Jenkins, Traub 1970] M. A. Jenkins and J. F. 
!    Traub, Principles for testing polynomial 
!    zerofinding programs, ACM Transactions on 
!    Mathematical Software, 1 (1975), pp. 26–34.
!
! 
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! Remark: In the paper we include a comparison
!         with BBEGG and BEGG. Since we cannot
!         distribute their code, these 
!         comparisons have been removed.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program main 
   
  implicit none 
  
  complex(kind(1d0)), allocatable :: poly(:),roots(:),allroots(:,:)
  complex(kind(1d0)), allocatable :: B(:,:),Y(:,:),cwork(:)  
  integer, allocatable :: iterations(:),its(:)
  double precision, allocatable :: res(:,:)
  double precision :: time,begg_time,lapack_time,pdble,error(6)
  double precision, allocatable :: Q(:),D(:),C(:),B2(:)
  double precision, allocatable :: rcoeffs(:),icoeffs(:),rnew(:),inew(:),alpha(:)
  double precision, allocatable :: reigs(:),ieigs(:)
  !double precision, allocatable :: rn(:),in(:)
  complex(kind(1d0)), allocatable :: eigs(:)
  double precision, allocatable ::residuals(:,:)

  complex(kind(1d0)) :: coeffs

  integer :: ii,jj,kk,ll,mm,N,zero,flag,info=0,str,stp,nnew,it,R
  integer :: clock_start,clock_end,clock_rate, newtnum, c1,c2
  integer :: Deg(20),num_trials(20), lapack_trials, num_trials_error(22)
  character (len=*), parameter :: path = "./"
  character(len=32) :: arg
   
  open (unit=7, file=path//"degrees.txt", status='unknown', position="append")
  open (unit=8, file=path//"damvw_times.txt", status='unknown', position="append")
  open (unit=18, file=path//"damvw_errors.txt", status='unknown', position="append")
  open (unit=28, file=path//"fig_damvw_times.txt", status='unknown', position="append")
  open (unit=38, file=path//"fig_damvw_errors.txt", status='unknown', position="append")
  open (unit=11, file=path//"lapack_times.txt", status='unknown', position="append")
  open (unit=21, file=path//"lapack_errors.txt", status='unknown', position="append")
  open (unit=31, file=path//"fig_lapack_times.txt", status='unknown', position="append")
  open (unit=41, file=path//"fig_lapack_errors.txt", status='unknown', position="append")
  
  if (iargc()>0) then
     call getarg(1, arg)
     read (arg,'(I10)') R
  else
     R = 5d0
  end if
 
  num_trials(1) = 2**(14)
  num_trials(2) = 2**(14)
  num_trials(3) = 2**(14)
  num_trials(4) = 2**(13)
  num_trials(5) = 2**(12)
  num_trials(6) = 2**(11)
  num_trials(7) = 2**(10)
  num_trials(8) = 2**(09)
  num_trials(9) = 2**(08)
  num_trials(10) = 2**(07)
  num_trials(11) = 2**(06)
  num_trials(12) = 2**(05)
  num_trials(13) = 2**(04)
  num_trials(14) = 2**(03)
  num_trials(15) = 2**(02)
  num_trials(16) = 2**(01)
  num_trials(17) = 2**(00)
  num_trials_error = 10

  Deg(1)  = 6
  Deg(2)  = 7
  Deg(3)  = 8
  Deg(4)  = 10
  Deg(5)  = 12
  Deg(6)  = 14
  Deg(7)  = 16
  Deg(8)  = 32
  Deg(9)  = 64
  Deg(10) = 128
  Deg(11) = 256
  Deg(12) = 512
  Deg(13) = 1024
  Deg(14) = 2048
  Deg(15) = 4096
  Deg(16) = 8192
  Deg(17) = 16384
  
  do kk=1,17
     write (7,*) Deg(kk), num_trials(kk)
     print*, num_trials(kk)
  end do
  
  call init_random_seed()


  write (28,*) "% random polynomials, JT", R
  write (38,*) "% random polynomials, JT", R
  write (31,*) "% random polynomials, JT", R
  write (41,*) "% random polynomials, JT", R

  write (28,*) "\addplot coordinates{ % AMVW"
  write (38,*) "\addplot coordinates{ % AMVW"

  write (31,*) "\addplot coordinates{ % LAPACK ZHSEQR" 
  write (41,*) "\addplot coordinates{ % LAPACK ZHSEQR"
  
  ! set newtnum
  newtnum = 1

  do ll=1,3
     
     time = 0d0
     
     N = Deg(ll)



     if (num_trials(ll) > num_trials_error(ll)) then
        mm = num_trials(ll)
     else
        mm = num_trials_error(ll)
     end if

     ! allocate memory
     allocate(Q(3*n),D(2*(n+1)),C(3*n),B2(3*n),rcoeffs(n*mm))
     allocate(icoeffs(n*mm))
     allocate(rnew(n*mm))
     allocate(inew(n*mm))
     allocate(its(n),reigs(n),ieigs(n),alpha(mm))
     allocate(residuals(n,3*(newtnum+1)))
     
     write(*,*) "Current degree =",N, num_trials(ll), mm
     
     allocate(poly((N+1)*mm),roots(N),allroots(N,2),iterations(N),res(N,6))
     allocate(B(N,N),Y(N,N),cwork(N))
 
     rcoeffs=0d0
     icoeffs=0d0
    
     do ii=1,mm
        rcoeffs(N*ii)=-1d0
        
        do jj=1,N
           poly(jj+(ii-1)*n) = complex(rcoeffs(jj+(ii-1)*n),icoeffs(jj+(ii-1)*n))
        end do
     end do

     error = 0d0
     
     ! ZAMVW
     write(*,*) "ZAMVW"
     CALL SYSTEM_CLOCK(COUNT_RATE=clock_rate) 
     CALL SYSTEM_CLOCK(COUNT=clock_start)
     
     do ii=1,num_trials(ll)
 	call factor(n,rcoeffs((ii-1)*n+1:(ii)*n),icoeffs((ii-1)*n+1:(ii)*n),Q,D,C,B2)
        call zamvw2(n,Q,D,C,B2,reigs,ieigs,its,flag,n-1,0)       

     end do
     

     CALL SYSTEM_CLOCK(COUNT=clock_end)  
     time = dble(clock_end - clock_start)/dble(clock_rate) 
     time = time/dble(num_trials(ll))

     do ii=1,num_trials_error(ll)
	call factor(n,rcoeffs((ii-1)*n+1:(ii)*n),icoeffs((ii-1)*n+1:(ii)*n),Q,D,C,B2)
	call zamvw2(n,Q,D,C,B2,reigs,ieigs,its,flag,n-1,0)

        it = 0
        do jj=1,n
           it = it + its(jj)
           roots(jj) = complex(reigs(jj),ieigs(jj))
           !print*, roots(jj), abs(roots(jj))
	end do
        call RESCHECK(0,N,0,1,POLY((ii-1)*n+1:ii*n),COEFFS,ROOTS,ALLROOTS,RES)
        
        do jj=1,6
           do kk=1,N
              if(res(kk,jj) > error(jj))then
                 error(jj) = res(kk,jj)
              end if
           end do
        end do
        
     end do


     print*, time
     print*, error(:)

     write (8,*) time
     write (18,*) error(:)

     write(*,"(A,I7,1x,A,1x,ES10.4E2,1x,A)"), "(", deg(ll), ",", &
          & time, ")%"
     write(28,"(A,I7,1x,A,1x,ES10.4E2,1x,A)"), "(", deg(ll), ",", &
          & time, ")%"

     write(*,"(A,I7,1x,A,1x,ES10.4E2,1x,A)"), "(", deg(ll), ",", &
          & error(1), ")%"
     write(38,"(A,I7,1x,A,1x,ES10.4E2,1x,A)"), "(", deg(ll), ",", &
          & error(1), ")%"

     if (ll == 15) then 
        write(28,"(A,1x,F6.3,1x,A,1x,I7,1x,A,1x,ES10.4E2,1x,A)"),  &
             &"%\node[coordinate,pin=below:{AMVW:", time,&
             &"}] at (axis cs:",deg(ll),",",time,"){};"
        write(*,"(A,1x,F6.3,1x,A,1x,I7,1x,A,1x,ES10.4E2,1x,A)"),  &
             &"%\node[coordinate,pin=below:{AMVW:", time,&
             &"}] at (axis cs:",deg(ll),",",time,"){};"
     end if


     
     ! LAPACK  
     if (ll<=14) then
        write(*,*) "LAPACK"
        error = 0d0
        res = 0d0
        CALL SYSTEM_CLOCK(COUNT_RATE=clock_rate) 
        CALL SYSTEM_CLOCK(COUNT=clock_start)
        do ii=1,num_trials(ll)
           do jj=1,N
              poly(jj) = complex(rcoeffs(jj+(ii-1)*n),icoeffs(jj+(ii-1)*n))
           end do
        
           B = complex(0.d0,0.d0)
           B(1,:) = -poly(2:N+1)
           do kk=1,(N-1)
              B(kk+1,kk) = complex(1.d0,0.d0)
           end do
           
           call ZHSEQR('E','N',N,1,N,B,N,roots,Y,N,cwork,N,info)
        end do
        
        CALL SYSTEM_CLOCK(COUNT=clock_end)  
        time = dble(clock_end - clock_start)/dble(clock_rate) 
        time = time/dble(num_trials(ll))
        do ii=1,num_trials_error(ll)
           do jj=1,N
              poly(jj) = complex(rcoeffs(jj+(ii-1)*n),icoeffs(jj+(ii-1)*n))
           end do
        
           B = complex(0.d0,0.d0)
           B(1,:) = -poly(2:N+1)
           do kk=1,(N-1)
              B(kk+1,kk) = complex(1.d0,0.d0)
           end do
           
           call ZHSEQR('E','N',N,1,N,B,N,roots,Y,N,cwork,N,info)
           call RESCHECK(0,N,0,1,POLY(2:N+1),COEFFS,ROOTS,ALLROOTS,RES)
           
           do jj=1,6
              do kk=1,N
                 if(res(kk,jj) > error(jj))then
                    error(jj) = res(kk,jj)
                 end if
              end do
           end do
        end do
        
        print*, time
        print*, error(:)
        write (11,*) time
        write (21,*) error(:)
        write(*,"(A,I7,1x,A,1x,ES10.4E2,1x,A)"), "(", deg(ll), ",", &
             & time, ")%"
        write(31,"(A,I7,1x,A,1x,ES10.4E2,1x,A)"), "(", deg(ll), ",", &
             & time, ")%"
        
        write(*,"(A,I7,1x,A,1x,ES10.4E2,1x,A)"), "(", deg(ll), ",", &
             & error(1), ")%"
        write(41,"(A,I7,1x,A,1x,ES10.4E2,1x,A)"), "(", deg(ll), ",", &
             & error(1), ")%"
     end if
     
     
     ! free memory
     deallocate(Q,D,C,B2,rcoeffs,icoeffs)
     deallocate(its,reigs,ieigs,rnew,inew)
     deallocate(residuals)
     deallocate(poly,roots,allroots,iterations,res)
     deallocate(B,Y,cwork,alpha)

  end do

  write (28,*) "};"
  write (38,*) "};"
  write (31,*) "};"
  write (41,*) "};"
  
  write (28,*) ""
  write (38,*) ""
  write (31,*) ""
  write (41,*) ""
  

  close(8)
  close(18) 
  close(28)
  close(38) 
  close(11)
  close(21)
  close(31)
  close(41)


end program
