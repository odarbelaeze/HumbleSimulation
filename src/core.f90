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

use lattice_utils
use core_utils

implicit none

! Counters

integer :: i, j, k, id

! Lattice stuff

integer, dimension (:), allocatable :: nnb
integer, dimension (:,:), allocatable :: nbh
real, dimension (:,:), allocatable :: d, r
integer :: w, h, l, n
integer, dimension(3) :: pbc
character :: lt_type
real :: rc

! Physical stuff

real, dimension (:,:), allocatable :: s
real, dimension (:),   allocatable :: smp
real, dimension (3) :: dir_B, mag, magt, ek
real :: T_max, T_min, T, J_ex, dT
real :: norm_B_max, norm_B, dnorm_B, norm_s, norm_k, smp_val

! Monte Carlo stuff

integer :: mcs_max, mcs_c, mcs

! Boltzmann constant

real, parameter :: k_B = 8.617E-2

! File names

integer :: argc
character (len = 50) :: f_in, work, f_out
logical :: ex

! Measure instruments

real, dimension(:,:), allocatable :: magnetization
real, dimension(:),   allocatable :: energy

! Namelists

namelist /Lattice/ w, h, l, pbc, lt_type, rc
namelist /Physical/ T_max, T_min, dT, T,                    &
                    norm_B_max, dir_B, norm_B, dnorm_B,     &
                    J_ex, norm_s, smp_val, norm_k
namelist /Monte/ mcs_c, mcs_max

! Setting default variables

w = 10; h = 10; l = 4; pbc = 1.0; lt_type = 's'; rc = 1.1
T_max = 1000; T_min = 0; dT = 10; T = 200
norm_B_max = 20; dir_B = (/ 0,0,1 /); norm_B = 10.0; dnorm_B = 0.5
J_ex = 5.5; norm_s = 1.5; smp_val = 3.0; norm_k = 4.0
mcs_c = 5000; mcs_max = 10000

! Checking comand line arguments

argc = iargc()

if (argc .lt. 2) then
    write(*,*) "# You must especify at least a file input and a work for output..."
    call exit(1)
end if

call getarg(1, f_in)
call getarg(2, work)

if (access(trim(f_in), 'r') .ne. 0) then
    write(*,*) "# You don't have the permissions to read: ", trim(f_in)
    call exit(2)
end if

! if (access(trim(path) // "/.", ' ') .ne. 0) then
!     write(*,*) "# The folder: ", trim(path), " doesn't exist and will be created."
!     call system("mkdir " // trim(path))
! end if

! Reading namelists from file

open(1, file = trim(f_in), status = 'old')
read(1,nml = Lattice)
read(1,nml = Physical)
read(1,nml = Monte)
close(1)

!TODO: Validate all the input variables (retvalue 4).

! Creating lattice

select case (lt_type)
    case ('s') ! Stands for sc
        call lt_sc (r, w, h, l, n)
        write (*,*) "# Lattice created succesfully: ", lt_type
    case ('b') ! Stands for bcc
        call lt_bcc (r, w, h, l, n)
        write (*,*) "# Lattice created succesfully: ", lt_type
    case ('f') ! Stands for fcc
        call lt_fcc (r, w, h, l, n)
        write (*,*) "# Lattice created succesfully: ", lt_type
    case ('h') ! Stands for hcp
        call lt_hcp (r, w, h, l, n)
        write (*,*) "# Lattice created succesfully: ", lt_type
    case default
        write (*,*) "# Lattice not known: ", lt_type
        call exit(4)
end select

! Plotting lattice file

! open (1, file = trim(path) // "/structure.out", status = "unknown")
! do i = 1, n
!     write (1,*) r(i,:)
! end do
! close(1)

! Calculating boundary conditions

allocate (nnb(n), d(n,n))
nnb = 0

do i = 1, n
    do j = 1, n
        d(i,j) = sqrt (sum (min (abs(r(j,:) - r(i,:)), abs((pbc * (/w,h,l/)) - (r(j,:) - r(i,:)))) ** 2))
        if (d(i,j) .le. rc .and. i .ne. j) nnb(i) = nnb(i) + 1
    end do
end do

allocate (nbh (n, maxval(nnb)))

do i = 1, n
    id = 1
    do j = 1, n
        if (d(i,j) .le. rc .and. i .ne. j) then
            nbh(i, id) = j
            id = id + 1
        end if
    end do
end do

call init_random_seed()
    
allocate (s(n,3), smp(n))

do i = 1, n
    call rdn_vec (s(i,:), norm_s)
end do

smp = smp_val

! Deallocating unnecesary variables

deallocate (d, r)

! magnetization 

ek  = (/ 0, 0, 1 /)

if (trim(work) == "hloop") then

    write(*,*) "# Writing hysteresys loop..."
    write(*,*) "# With velocity: ", mcs_c

    do norm_B = 0, norm_B_max, dnorm_B

        call measure(                       &
            T, norm_B, dir_B, mcs_c,        &
            s, smp,                         &
    	    ek, norm_k,                     &
            n, nnb, nbh,                    &
            energy, magnetization           &
        )
        
        write (*,*) norm_B, sum(energy) / mcs_c, sum(magnetization(:,3)) / mcs_c

    end do

    do norm_B = norm_B_max, - norm_B_max, - dnorm_B

        call measure(                       &
            T, norm_B, dir_B, mcs_c,        &
            s, smp,                         &
    	    ek, norm_k,                     &
            n, nnb, nbh,                    &
            energy, magnetization           &
        )
        
        write (*,*) norm_B, sum(energy) / mcs_c, sum(magnetization(:,3)) / mcs_c

    end do

    do norm_B = - norm_B_max, norm_B_max, dnorm_B

        call measure(                       &
            T, norm_B, dir_B, mcs_c,        &
            s, smp,                         &
    	    ek, norm_k,                     &
            n, nnb, nbh,                    &
            energy, magnetization           &
        )
        
        write (*,*) norm_B, sum(energy) / mcs_c, sum(magnetization(:,3)) / mcs_c

    end do

else if (trim(work) == "magcurv") then

    write(*,*) "# Writing magnetization curves..."
    write(*,*) "# With velocity: ", mcs_c

    do T = T_max, T_min, - dT
        call measure(                       &
            T, norm_B, dir_B, mcs_c,        &
            s, smp,                         &
            ek, norm_k,                     &
            n, nnb, nbh,                    &
            energy, magnetization           &
        )
        
        write (*,*) T, sum(energy) / mcs_c,  & 
            sum(sqrt(magnetization(:,1) ** 2 + magnetization(:,2) ** 2 + magnetization(:,3) ** 2)) / mcs_c
    end do

else

    write(*,*) "# Sorry unknown work:", work

end if

end program core
