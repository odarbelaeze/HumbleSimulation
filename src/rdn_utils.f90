!
!   Copyright 2012 Oscar David Arbeláez, Juan Camilo Henao
!
!   This file is part of HumbleSimulation.
!
!   HumbleSimulation is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
! 
!   HumbleSimulation is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
! 
!   You should have received a copy of the GNU General Public License
!   along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
!

module rdn_utils

implicit none

contains

subroutine rdn_vec (rv, norm)

    real, dimension (3), intent (out) :: rv
    real, intent (in) :: norm

    call random_number (rv)
    rv = rv - 0.5
    rv = norm * rv / sqrt (sum (rv ** 2))

end subroutine rdn_vec

subroutine rand_changed_vec (rcv, vec, r)

    real, dimension (3), intent (out) :: rcv
    real, dimension (3), intent (in)  :: vec
    real, dimension (3) :: rtemp
    real, intent (in) :: r
    real :: norm

    norm = sqrt (sum (vec ** 2))

    call random_number (rtemp)
    rtemp = rtemp - 0.5
    rtemp = r * rtemp / sqrt (sum (rtemp ** 2))
    rcv = vec + rtemp

    rcv = norm * rcv / sqrt (sum (rcv ** 2))

end subroutine rand_changed_vec

subroutine cheap_rand_changed_vec (rcv, vec)

    real, dimension (3), intent (out) :: rcv
    real, dimension (3), intent (in)  :: vec
    real, dimension (3) :: rtemp
    real :: norm

    norm = sqrt (sum (vec ** 2))

    call random_number (rtemp)
    rtemp = 2.0 * (rtemp - 0.5)
    
    rcv = vec + rtemp
    rcv = norm * rcv / sqrt (sum (rcv ** 2))

end subroutine cheap_rand_changed_vec

!TODO: Change subroutine init_random_seed.

subroutine init_random_seed ()

    integer :: i, n, clock
    integer, allocatable, dimension (:) :: seed

    call random_seed (size = n)
    allocate (seed (n))

    call system_clock (count = clock)
    seed = clock + 37 * (/ (i - 1, i = 1, n) /)

    call random_seed(put = seed)
    deallocate (seed)

end subroutine init_random_seed

end module rdn_utils
