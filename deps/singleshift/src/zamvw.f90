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
! zamvw computes roots of a complex polynomial
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! n         problem size
!
! rcoeffs   real parts polynomial coefficients
! icoeffs   imag parts polynomial coefficients
! 
! a_j = rcoeffs(j) + i*icoeffs(j)
! p(z) = z^n + a_1 z^{n-1} + a_2 z^(n-2} + ... + a_{n-1} z + a_n
!
!
! reigs     real parts roots
! ieigs     imag parts roots
!
! its       array, no of iteration between subsequent deflation
!
! flag      error flag
!           0    no error, all eigenvalues found
!           k>0  QR algorithm did not converge,
!                k eigenvalues are found (first k 
!                entries of reigs,ieigs)
!           
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine zamvw(n,rcoeffs,icoeffs,reigs,ieigs,its,flag)
  implicit none

  ! input variables
  integer, intent(in) :: n
  double precision, intent(in) :: rcoeffs(n), icoeffs(n)
  double precision, intent(inout) :: reigs(n), ieigs(n)
  integer, intent(inout) :: its(n), flag
  

  ! compute variables
  double precision :: Q(3*n),D(2*(n+1)),C(3*n),B(3*n)
  double precision :: nrm
  complex(kind(1d0)) :: trace, detm, disc, e1, e2
  integer :: ii, nnew, ry
  
!   print*, "here"
  ry = flag
  flag = 0
  
  !print*, rcoeffs, icoeffs

  reigs = 0d0
  ieigs = 0d0

  ! check dimension
  if (n<=0) then
     flag = -1
     return
  end if
  

  do ii=1,n
     nrm = abs(complex(rcoeffs(n+1-ii),icoeffs(n+1-ii)))
     if(nrm /= 0)then
        nnew = n+1-ii
        exit
     end if
  end do

  if (nnew <= 0) then
     ! polynomial of zero coefficients, return zeros as roots
     flag = 0
     return
  end if
  

  if (nnew == 1) then
     ! it remains a polynomial of degree 1 
     reigs(1) = -rcoeffs(1)
     ieigs(1) = -icoeffs(1)
     FLAG = 0
     return
  end if

  if (nnew == 2) then
     ! modified quadratic formula
     trace = -complex(rcoeffs(1),icoeffs(1))
     detm = complex(rcoeffs(2),icoeffs(2))
     disc = zsqrt(trace*trace - 4d0*detm)
     
     ! compute e1 and e2
     if(zabs(trace+disc) > zabs(trace-disc))then
        if(zabs(trace+disc) == 0)then
           reigs(1) = 0d0
           ieigs(1) = 0d0
           reigs(2) = 0d0
           ieigs(2) = 0d0
        else
           e1 = (trace+disc)/complex(2d0,0d0)
           e2 = detm/e1
           reigs(1) = dble(e1)
           ieigs(1) = imag(e1)
           reigs(2) = dble(e2)
           ieigs(2) = imag(e2)
        end if
     else
        if(zabs(trace-disc) == 0)then
           reigs(1) = 0d0
           ieigs(1) = 0d0
           reigs(2) = 0d0
           ieigs(2) = 0d0
        else
           e1 = (trace-disc)/complex(2d0,0d0)
           e2 = detm/e1
           reigs(1) = dble(e1)
           ieigs(1) = imag(e1)
           reigs(2) = dble(e2)
           ieigs(2) = imag(e2)
        end if
     end if
     flag = 0
     return
  end if

  call factor(n,rcoeffs,icoeffs,Q,D,C,B)
  call zamvw2(n,Q,D,C,B,reigs,ieigs,its,flag,n-1,ry)
  

end subroutine zamvw
