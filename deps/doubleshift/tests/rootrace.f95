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
! comparing DAMVW with LAPACK's DHSEQR and DGEEV
! polynomials with random coefficients
! 
! Remark: In the paper we include a comparison
!         with CGXZ (and BBEGG). Since we cannot
!         redistribute their code, these 
!         comparisons have been removed.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program main 
  
  implicit none 
   
  double precision, allocatable :: poly(:),C(:,:),rroots(:),iroots(:),Z(:,:),work(:)  
  double precision, allocatable :: res(:,:)
  complex(kind(1d0)), allocatable :: wpoly(:),wcoeffs(:,:),wroots(:),wallroots(:,:)  
  integer, allocatable :: iterations(:) 
  double precision :: time,error(6)
   
  integer :: ii,jj,kk,ll,mm,N,flag,info=0
  integer :: clock_start,clock_end,clock_rate 
  integer :: Deg(22),num_trials(22),num_trials_error(22)
  ! set path
  character (len=*), parameter :: path = "./"
  
  ! Files with fig_* contain the LaTeX lines of the Figures in the paper.
  open (unit=7, file=path//"degrees.txt", status='unknown', position="append")
  open (unit=8, file=path//"damvw_times.txt", status='unknown', position="append")
  open (unit=18, file=path//"damvw_errors.txt", status='unknown', position="append")
  open (unit=28, file=path//"fig_damvw_times.txt", status='unknown', position="append")
  open (unit=38, file=path//"fig_damvw_errors.txt", status='unknown', position="append")
  open (unit=11, file=path//"lapack_times.txt", status='unknown', position="append")
  open (unit=21, file=path//"lapack_errors.txt", status='unknown', position="append")
  open (unit=31, file=path//"fig_lapack_times.txt", status='unknown', position="append")
  open (unit=41, file=path//"fig_lapack_errors.txt", status='unknown', position="append")
  open (unit=33, file=path//"fig_dgeev_times.txt", status='unknown', position="append")
  open (unit=43, file=path//"fig_dgeev_errors.txt", status='unknown', position="append")

 
  ! the number of runs is for small polynomials higher
  ! => there are enough runs for small polynomials to 
  !    get timings
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
  ! => errors depend on the number of runs!
  num_trials_error = 10

  Deg(1)  = 6
  Deg(2)  = 7
  Deg(3)  = 8
  Deg(4)  = 10
  Deg(5)  = 12
  Deg(6)  = 14
  Deg(7)  = 16
  ! O(n²) for DAMVW and O(n³) for LAPACK 
  Deg(8)  = 32
  Deg(9)  = 64
  Deg(10) = 128
  Deg(11) = 256   
  Deg(12) = 512
  Deg(13) = 1024
  Deg(14) = 2048
  ! LAPACK is turnoff for larger polynomials, see below
  Deg(15) = 4096
  Deg(16) = 8192
  Deg(17) = 16384
  
  do kk=1,17
     write (7,*) Deg(kk), num_trials(kk)
  end do
  close(7)

  call init_random_seed()



  write (28,*) "% random polynomials, normal distributed coefficients"
  write (38,*) "% random polynomials, normal distributed coefficients"
  write (31,*) "% random polynomials, normal distributed coefficients"
  write (41,*) "% random polynomials, normal distributed coefficients"
  write (33,*) "% random polynomials, normal distributed coefficients"
  write (43,*) "% random polynomials, normal distributed coefficients"

  write (28,*) "\addplot coordinates{ % AMVW"
  write (38,*) "\addplot coordinates{ % AMVW"

  write (31,*) "\addplot coordinates{ % LAPACK DHSEQR"
  write (41,*) "\addplot coordinates{ % LAPACK DHSEQR"

  write (33,*) "\addplot coordinates{ % LAPACK DGEEV"
  write (43,*) "\addplot coordinates{ % LAPACK DGEEV"

  do ll=1,17
   
     time = 0d0
     
     N = Deg(ll)
     
     write(*,*) "Current N =",N, num_trials(ll)


     if (num_trials(ll) > num_trials_error(ll)) then
        mm = num_trials(ll)
     else
        mm = num_trials_error(ll)
     end if
     
     allocate(poly(N*mm),iterations(N),res(N,6))
     allocate(C(N,N),Z(N,N),work(5*N),rroots(N),iroots(N))
     allocate(wpoly(N),wCOEFFS(N,1),wROOTS(N),wALLROOTS(N,2))
     
     do ii=1,N
        wcoeffs(ii,1) = complex(1d0,0d0)
     end do
     
     error = 0d0
     res = 0d0
     ! choose random polynomials for all test runs
     call dnormalpoly(N*mm,poly(1:mm))
     
     
     ! DAMVW     
     write(*,*) "DAMVW"
     ! start timer
     CALL SYSTEM_CLOCK(COUNT_RATE=clock_rate) 
     CALL SYSTEM_CLOCK(COUNT=clock_start)
     do ii=1,num_trials(ll)
        
        call DAMVW(N,POLY((ii-1)*N+1:ii*N),RROOTS,IROOTS,ITERATIONS,FLAG)
        
     end do
     
     CALL SYSTEM_CLOCK(COUNT=clock_end)  
     time = dble(clock_end - clock_start)/dble(clock_rate) 
     time = time/dble(num_trials(ll))
         
     do ii=1,num_trials_error(ll)
          
        ! compute roots
        call DAMVW(N,POLY((ii-1)*N+1:ii*N),rroots,iroots,iterations,flag)
        
        do jj=1,N
           wpoly(jj) = complex(poly((ii-1)*N+jj),0d0)
           wroots(jj) = complex(rroots(jj),iroots(jj))
        end do

        ! check residual
        call RESCHECK(0,N,0,1,wpoly,wcoeffs,wroots,wallroots,res)
        
        do jj=1,6
           do kk=1,N
              if(res(kk,jj) > error(jj))then
                 error(jj) = res(kk,jj)
              end if
           end do
        end do
        
     end do
     
     
     write (8,*) time
     print*, time
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

     if (ll<=14) then
        ! LAPACK
        write(*,*) "LAPACK"
        error = 0d0
        res = 0d0
        
        CALL SYSTEM_CLOCK(COUNT_RATE=clock_rate) 
        CALL SYSTEM_CLOCK(COUNT=clock_start)
     
        do ii=1,num_trials(ll)
           
           C = 0d0
           C(1,:) = -POLY((ii-1)*N+1:ii*N)
           do kk=1,(N-1)
              C(kk+1,kk) = 1.d0
           end do
           
           call DHSEQR('E','N',N,1,N,C,N,rroots,iroots,Z,N,work,N,info)
           
        end do
        
        CALL SYSTEM_CLOCK(COUNT=clock_end)  
        time = dble(clock_end - clock_start)/dble(clock_rate) 
        time = time/dble(num_trials(ll))
        
        
        do ii=1,num_trials_error(ll)
           
           C = 0d0
           C(1,:) = -POLY((ii-1)*N+1:ii*N)
           do kk=1,(N-1)
              C(kk+1,kk) = 1.d0
           end do
           
           call DHSEQR('E','N',N,1,N,C,N,rroots,iroots,Z,N,work,N,info)
           
           do jj=1,N
              wpoly(jj) = complex(POLY((ii-1)*N+jj),0d0)
              wroots(jj) = complex(rroots(jj),iroots(jj))
           end do
           
           call RESCHECK(0,N,0,1,wPOLY,wCOEFFS,wROOTS,wALLROOTS,RES)
           
           do jj=1,6
              do kk=1,N
                 if(res(kk,jj) > error(jj))then
                    error(jj) = res(kk,jj)
                 end if
              end do
           end do
           
        end do
        write (11,*) time
        print*, time
        write (21,*) error(:)
        
        write(*,"(A,I7,1x,A,1x,ES10.4E2,1x,A)"), "(", deg(ll), ",", &
             & time, ")%"
        write(31,"(A,I7,1x,A,1x,ES10.4E2,1x,A)"), "(", deg(ll), ",", &
             & time, ")%"
        
        write(*,"(A,I7,1x,A,1x,ES10.4E2,1x,A)"), "(", deg(ll), ",", &
             & error(1), ")%"
        write(41,"(A,I7,1x,A,1x,ES10.4E2,1x,A)"), "(", deg(ll), ",", &
             & error(1), ")%"


        write(*,*) "LAPACK"
        error = 0d0
        res = 0d0
        
        CALL SYSTEM_CLOCK(COUNT_RATE=clock_rate) 
        CALL SYSTEM_CLOCK(COUNT=clock_start)

        do ii=1,num_trials(ll)
           
           C = 0d0
           C(1,:) = -POLY((ii-1)*N+1:ii*N)
           do kk=1,(N-1)
              C(kk+1,kk) = 1.d0
           end do
           
           call DGEEV('N','N',N,C,N,rroots,iroots,Z,N,Z,N,work,5*N,info)
           
        end do
        
        CALL SYSTEM_CLOCK(COUNT=clock_end)  
        time = dble(clock_end - clock_start)/dble(clock_rate) 
        time = time/dble(num_trials(ll))
        
        
        do ii=1,num_trials_error(ll)
           
           C = 0d0
           C(1,:) = -POLY((ii-1)*N+1:ii*N)
           do kk=1,(N-1)
              C(kk+1,kk) = 1.d0
           end do

           call DGEEV('N','N',N,C,N,rroots,iroots,Z,N,Z,N,work,5*N,info)
           
           do jj=1,N
              wpoly(jj) = complex(POLY((ii-1)*N+jj),0d0)
              wroots(jj) = complex(rroots(jj),iroots(jj))
           end do
           
           call RESCHECK(0,N,0,1,wpoly,wcoeffs,wroots,wallroots,res)
           
           do jj=1,6
              do kk=1,N
                 if(res(kk,jj) > error(jj))then
                    error(jj) = res(kk,jj)
                 end if
              end do
           end do
           
        end do
        write (11,*) time
        print*, time
        write (21,*) error(:)
        
        write(*,"(A,I7,1x,A,1x,ES10.4E2,1x,A)"), "(", deg(ll), ",", &
             & time, ")%"
        write(33,"(A,I7,1x,A,1x,ES10.4E2,1x,A)"), "(", deg(ll), ",", &
             & time, ")%"
        
        write(*,"(A,I7,1x,A,1x,ES10.4E2,1x,A)"), "(", deg(ll), ",", &
             & error(1), ")%"
        write(43,"(A,I7,1x,A,1x,ES10.4E2,1x,A)"), "(", deg(ll), ",", &
             & error(1), ")%"
     end if
     
     deallocate(poly,iterations,res)
     deallocate(C,Z,work,rroots,iroots)
     deallocate(wpoly,wcoeffs,wroots,wallroots)
     
  end do

  write (28,*) "};"
  write (38,*) "};"
  write (31,*) "};"
  write (41,*) "};"
  write (33,*) "};"
  write (43,*) "};"

  write (28,*) ""
  write (38,*) ""
  write (31,*) ""
  write (41,*) ""
  write (33,*) ""
  write (43,*) ""

  close(8)
  close(18) 
  close(28)
  close(38)
  close(11)
  close(21)
  close(31)
  close(41)
  close(13)
  close(23)
  close(33)
  close(43)


end program main
