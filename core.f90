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

program core

use rdn_utils
use lattice_utils

implicit none

! Counters

integer :: i, j, k, id

! Lattice stuff

real, dimension (:), allocatable :: nnb
real, dimension (:,:), allocatable :: nbh, d, r
integer :: w, h, l, n
integer, dimension(3) :: pbc
character :: lt_type

! Physical stuff

real, dimension (:,:), allocatable :: s
real, dimension (3) :: B, mag, magt
real :: E1, E2, E, Et
real :: T_max, T_min, T, J_ex, dT
real :: B_max, B_min, norm_B, dnorm_B, norm_s

! Monte Carlo stuff

integer :: mcs_max, mcs_c

! Namelists

namelist /Lattice/ r, w, h, l, n, pbc
namelist /Physical/ norm_s, T_max, T_min, B_max, B_min, J_ex
namelist /Monte/ mcs_c, mcs_max

! Setting variables

w = 10; h = 10; l = 4; lt_type = 'b'
B = (/ 0,0,1 /)
T_max = 500; T_min = 10; J_ex = 5.5; dT = 10
B_max = 6; B_min = - B_max; dnorm_B = 0.1; norm_s = 1.0
mcs_max = 50000; mcs_c = 25000

! Creating lattice

write (*,*) "Creating latice...", lt_type

select case (lt_type)
    case ('s') ! Stands for sc
        call lt_sc (r, w, h, l, n)
        write (*,*) "Lattice created succesfully: ", lt_type
    case ('b') ! Stands for bcc
        call lt_bcc (r, w, h, l, n)
        write (*,*) "Lattice created succesfully: ", lt_type
    case ('f') ! Stands for fcc
        call lt_fcc (r, w, h, l, n)
        write (*,*) "Lattice created succesfully: ", lt_type
    case ('h') ! Stands for hcp
        call lt_hcp (r, w, h, l, n)
        write (*,*) "Lattice created succesfully: ", lt_type
    case default
        write (*,*) "Lattice not known: ", lt_type
end select

end program core
