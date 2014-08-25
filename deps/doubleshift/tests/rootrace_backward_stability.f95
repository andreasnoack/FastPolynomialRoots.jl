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
! special backward stability test Table 8.16 and
! Table 8.17 in the paper
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program main 
  
  implicit none 
   
  double precision, allocatable :: poly(:),C(:,:),rroots(:),iroots(:),Z(:,:),work(:)  
  double precision, allocatable :: coeffs(:,:),rallroots(:,:),iallroots(:,:),res(:,:)
  complex(kind(1d0)), allocatable :: wpoly(:),wCOEFFS(:,:),wROOTS(:),wALLROOTS(:,:)  
  integer, allocatable :: iterations(:) 
  double precision :: time,error(6),lerror(6), lzerror(6), alpha
   
  integer :: ii,jj,kk,ll,mm,N,zero,flag,info=0
  integer :: clock_start,clock_end,clock_rate 
  integer :: Deg(400), num_trials(400), lapack_trials
  real :: sss_time
  character (len=*), parameter :: path = "/home/thomasm/data/"
  double precision :: normx(400),azero(400), normofp

  ! BLAS functions
  double precision :: dznrm2, dnrm2
  
  open (unit=7, file="table_bs.txt", status='unknown')

  write (7,*) "Table 1"
  
  do kk=1,1!,20
     do ll=1,20
        num_trials(20*(kk-1)+ll) = 2**(00)
        Deg(20*(kk-1)+ll) = 32
        normx(ll) = 10d0**(ll-8)
        azero(ll) = 0d0        
     end do
  end do

  call init_random_seed()

  do ll=8,20
   
     time = 0d0
     
     N = Deg(ll)
     
     print*, "Current N =",N, num_trials(ll)
     
     allocate(poly(N*num_trials(ll)),coeffs(N,1),iterations(N),res(N,6))
     allocate(C(N,N),Z(N,N),work(5*N),rroots(N),iroots(N))
     allocate(wpoly(N),wCOEFFS(N,1),wROOTS(N),wALLROOTS(N,2))
     
     do ii=1,N
        coeffs(ii,1) = 1d0
        wcoeffs(ii,1) = complex(1d0,0d0)
     end do
    
     error = 0d0


     call dnormalpoly(N*num_trials(ll),poly(1:N*num_trials(ll)))
     do ii=1,num_trials(ll)
        if (normx(ll)/=0d0) then
           alpha = normx(ll)/dnrm2(N,poly((ii-1)*N+1:ii*N),1)
           call dscal(N,alpha,poly((ii-1)*N+1:ii*N),1)
        end if
        if (azero(ll)/=0d0) then
           POLY(ii*N)=POLY(ii*N)/abs(POLY(ii*N))*azero(ll)
        end if
        do jj=1,N
           wpoly(jj) = complex(POLY((ii-1)*N+jj),0d0)
        end do
     end do

     normofp = dsqrt(dznrm2(N,wpoly,1)**2 + 1d0)


     ! DAMVW
     
     write(*,*) "DAMVW"
         
     do ii=1,num_trials(ll)
          
        call DAMVW(N,POLY((ii-1)*N+1:ii*N),RROOTS,IROOTS,ITERATIONS,FLAG)
        
        do jj=1,N
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
     

     ! LAPACK
     write(*,*) "LAPACK"
     lerror = 0d0
     res = 0d0

     do ii=1,num_trials(ll)
        C = 0d0
        ! quadratic formula test from Edelman-Murakami-1995
        if (n==2) then
           C(1,2) = -POLY(2)
           C(2,2) = -POLY(1)
        else if (n==4) then
           C(1,4) = -POLY(4)
           C(2,4) = -POLY(3)
           C(3,4) = -POLY(2)
           C(4,4) = -POLY(1)
        else if (n==3) then
           C(1,3) = -POLY(3)
           C(2,3) = -POLY(2)
           C(3,3) = -POLY(1)
        else
           C(1,:) = -POLY((ii-1)*N+1:ii*N) 
        end if

        do kk=1,(N-1)
           C(kk+1,kk) = 1.d0
        end do
        
        call DHSEQR('E','N',N,1,N,C,N,rroots,iroots,Z,N,work,N,info)
        
        do jj=1,N
           !wpoly(jj) = complex(POLY((ii-1)*N+jj),0d0)
           wroots(jj) = complex(rroots(jj),iroots(jj))
        end do
        
        call RESCHECK(0,N,0,1,wPOLY,wCOEFFS,wROOTS,wALLROOTS,RES)
        print*, "lapack info", info
        do jj=1,6
           do kk=1,N
              if(res(kk,jj) > lerror(jj))then
                 lerror(jj) = res(kk,jj)
              end if
           end do
        end do        
     end do



     ! LAPACK
     write(*,*) "LAPACK"
     lzerror = 0d0
     res = 0d0

     do ii=1,num_trials(ll)
        C = 0d0
        ! quadratic formula test from Edelman-Murakami-1995
        if (n==2) then
           C(1,2) = -POLY(2)
           C(2,2) = -POLY(1)
        else if (n==4) then
           C(1,4) = -POLY(4)
           C(2,4) = -POLY(3)
           C(3,4) = -POLY(2)
           C(4,4) = -POLY(1)
        else if (n==3) then
           C(1,3) = -POLY(3)
           C(2,3) = -POLY(2)
           C(3,3) = -POLY(1)
        else
           C(1,:) = -POLY((ii-1)*N+1:ii*N) 
        end if

        do kk=1,(N-1)
           C(kk+1,kk) = 1.d0
        end do
        
        !call DHSEQR('E','N',N,1,N,C,N,rroots,iroots,Z,N,work,N,info)
        call DGEEV('N','N',N,C,N,rroots,iroots,Z,N,Z,N,work,5*N,info)
        
        do jj=1,N
           !wpoly(jj) = complex(POLY((ii-1)*N+jj),0d0)
           wroots(jj) = complex(rroots(jj),iroots(jj))
        end do
        
        call RESCHECK(0,N,0,1,wPOLY,wCOEFFS,wROOTS,wALLROOTS,RES)
        print*, "lapack info", info
        do jj=1,6
           do kk=1,N
              if(res(kk,jj) > lzerror(jj))then
                 lzerror(jj) = res(kk,jj)
              end if
           end do
        end do        
     end do
     
     


     write (7,"(ES10.0E2,1x,A,1x,ES10.4E2,1x,A,1x,&
          &ES10.4E2,1x,A,1x,ES10.4E2,1x,A,1x,ES10.4E2,1x,A,I3)"), normx(ll),&
          & "&", error(1), "&", lerror(1), "&", lzerror(1), "&", normofp, "\\%", ll
     write (*,"(ES10.0E2,1x,A,1x,ES10.4E2,1x,A,1x,&
          &ES10.4E2,1x,A,1x,ES10.4E2,1x,A,1x,ES10.4E2,1x,A,I3)"), normx(ll),&
          & "&", error(1), "&", lerror(1), "&", lzerror(1), "&", normofp, "\\%", ll

     deallocate(poly,coeffs,iterations,res)
     deallocate(C,Z,work,rroots,iroots)
     deallocate(wpoly,wCOEFFS,wROOTS,wALLROOTS)
     
  end do

  write (7,*) "Table 2"

  do kk=1,1!,20
     do ll=1,20
        num_trials(20*(kk-1)+ll) = 2**(00)
        Deg(20*(kk-1)+ll) = 512
        normx(ll) = 0d0
        azero(ll) = 10d0**(ll-8)        
     end do
  end do

  do ll=1,10
   
     time = 0d0
     
     N = Deg(ll)
     
     print*, "Current N =",N, num_trials(ll)
     
     allocate(poly(N*num_trials(ll)),coeffs(N,1),iterations(N),res(N,6))
     allocate(C(N,N),Z(N,N),work(5*N),rroots(N),iroots(N))
     allocate(wpoly(N),wCOEFFS(N,1),wROOTS(N),wALLROOTS(N,2))
     
     do ii=1,N
        coeffs(ii,1) = 1d0
        wcoeffs(ii,1) = complex(1d0,0d0)
     end do
    
     error = 0d0


     call dnormalpoly(N*num_trials(ll),poly(1:N*num_trials(ll)))
     do ii=1,num_trials(ll)
        if (normx(ll)/=0d0) then
           alpha = normx(ll)/dnrm2(N,poly((ii-1)*N+1:ii*N),1)
           call dscal(N,alpha,poly((ii-1)*N+1:ii*N),1)
        end if
        if (azero(ll)/=0d0) then
           POLY(ii*N)=POLY(ii*N)/abs(POLY(ii*N))*azero(ll)
        end if
        do jj=1,N
           wpoly(jj) = complex(POLY((ii-1)*N+jj),0d0)
        end do
     end do

     normofp = dsqrt(dznrm2(N,wpoly,1)**2 + 1d0)

     ! DAMVW    
     write(*,*) "DAMVW"
         
     do ii=1,num_trials(ll)
          
        call DAMVW(N,POLY((ii-1)*N+1:ii*N),RROOTS,IROOTS,ITERATIONS,FLAG)
        
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
     

     ! LAPACK
     write(*,*) "LAPACK"
     lerror = 0d0
     res = 0d0

     do ii=1,num_trials(ll)
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
        print*, "lapack info", info
        do jj=1,6
           do kk=1,N
              if(res(kk,jj) > lerror(jj))then
                 lerror(jj) = res(kk,jj)
              end if
           end do
        end do        
     end do
     
     
     ! LAPACK
     write(*,*) "LAPACK"
     lzerror = 0d0
     res = 0d0

     do ii=1,num_trials(ll)
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
        
        call RESCHECK(0,N,0,1,wPOLY,wCOEFFS,wROOTS,wALLROOTS,RES)
        print*, "lapack info", info
        do jj=1,6
           do kk=1,N
              if(res(kk,jj) > lzerror(jj))then
                 lzerror(jj) = res(kk,jj)
              end if
           end do
        end do        
     end do
     
     


     write (7,"(ES10.0E2,1x,A,1x,ES10.4E2,1x,A,1x,&
          &ES10.4E2,1x,A,1x,ES10.4E2,1x,A,1x,ES10.4E2,1x,A,I3)"), azero(ll),&
          & "&", error(1), "&", lerror(1), "&", lzerror(1), "&", normofp, "\\%", ll
     write (*,"(ES10.0E2,1x,A,1x,ES10.4E2,1x,A,1x,&
          &ES10.4E2,1x,A,1x,ES10.4E2,1x,A,1x,ES10.4E2,1x,A,I3)"), azero(ll),&
          & "&", error(1), "&", lerror(1), "&", lzerror(1), "&", normofp, "\\%", ll

     deallocate(poly,coeffs,iterations,res)
     deallocate(C,Z,work,rroots,iroots)
     deallocate(wpoly,wCOEFFS,wROOTS,wALLROOTS)
     
  end do
  close(7)


end program main
