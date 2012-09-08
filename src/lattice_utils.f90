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

module lattice_utils

implicit none

contains

subroutine lt_sc (r, w, h, l, n)

    real, allocatable, dimension (:,:), intent (out) :: r
    integer, intent (in) :: w, h, l
    integer, optional, intent (out) :: n
    integer :: id = 1, i, j, k

    n = w * h * l
    allocate (r(n,3))

    do i = 0, w - 1
        do j = 0, h - 1
            do k = 0, l - 1
                r(id, :)  = (/ i, j, k /) + (/ 0.0, 0.0, 0.0 /); id = id + 1
            end do
        end do
    end do



end subroutine lt_sc

subroutine lt_bcc (r, w, h, l, n)

    real, allocatable, dimension (:,:), intent (out) :: r
    integer, intent (in) :: w, h, l
    integer, optional, intent (out) :: n
    integer :: id = 1, i, j, k
    n = 2 * w * h * l
    allocate (r(n,3))

    do i = 0, w - 1
        do j = 0, h - 1
            do k = 0, l - 1
                r(id, :)  = (/ i, j, k /) + (/ 0.0, 0.0, 0.0 /); id = id + 1
                r(id, :)  = (/ i, j, k /) + (/ 0.5, 0.5, 0.5 /); id = id + 1
            end do
        end do
    end do

end subroutine lt_bcc

subroutine lt_fcc (r, w, h, l, n)

    real, allocatable, dimension (:,:), intent (out) :: r
    integer, intent (in) :: w, h, l
    integer, optional, intent (out) :: n
    integer :: id = 1, i, j, k
    n = 4 * w * h * l
    allocate (r(n,3))

    do i = 0, w - 1
        do j = 0, h - 1
            do k = 0, l - 1
                r(id, :)  = (/ i, j, k /) + (/ 0.0, 0.0, 0.0 /); id = id + 1
                r(id, :)  = (/ i, j, k /) + (/ 0.5, 0.0, 0.0 /); id = id + 1
                r(id, :)  = (/ i, j, k /) + (/ 0.0, 0.5, 0.0 /); id = id + 1
                r(id, :)  = (/ i, j, k /) + (/ 0.0, 0.0, 0.5 /); id = id + 1
            end do
        end do
    end do

end subroutine lt_fcc

subroutine lt_hcp (r, w, h, l, n)

    real, allocatable, dimension (:,:), intent (out) :: r
    integer, intent (in) :: w, h, l
    integer, optional, intent (out) :: n
    integer :: id = 1, i, j, k
    n = 6 * w * h * l
    allocate (r(n,3))

    do i = 0, w - 1
        do j = 0, h - 1
            do k = 0, l - 1
                r(id, :)  = (/ i, j, k /) + (/ 0.0, 0.0, 0.0 /); id = id + 1
                r(id, :)  = (/ i, j, k /) + (/ 0.5, 0.5, 0.0 /); id = id + 1
                r(id, :)  = (/ i, j, k /) + (/ 0.5 * sqrt(2.0), 0.5 * sqrt(2.0), 0.5 /); id = id + 1
                r(id, :)  = (/ i, j, k /) + (/ 0.5 * sqrt(2.0) + 0.5, 0.5 * sqrt(2.0), 0.5 /); id = id + 1
                r(id, :)  = (/ i, j, k /) + (/ 0.5 * sqrt(2.0), 0.5 * sqrt(2.0) + 0.5, 0.5 /); id = id + 1
                r(id, :)  = (/ i, j, k /) + (/ 0.5 * sqrt(2.0), 0.5 * sqrt(2.0), 0.0 /) + 0.5; id = id + 1
            end do
        end do
    end do

end subroutine lt_hcp

end module lattice_utils
