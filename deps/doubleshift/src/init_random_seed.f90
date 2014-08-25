!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize Random Seed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine initializes the random number generator using the 
! CPU clock.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine INIT_RANDOM_SEED()

	implicit none

	! compute variables
        integer :: ii, n, clock
        integer, allocatable :: seed(:)
          
        call RANDOM_SEED(size = n)
        allocate(seed(n))
          
        call SYSTEM_CLOCK(COUNT=clock)
          
        seed = clock + 37 * (/ (ii - 1, ii = 1, n) /)
        call RANDOM_SEED(PUT = seed)
          
	deallocate(seed)

end subroutine
