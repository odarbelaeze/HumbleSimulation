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

module core_utils

use rdn_utils

implicit none

contains


	subroutine estabilize(              &
		T, norm_B, dir_B, mcs,          &
		s, smp,                         &
	    ek, norm_k,                     &
		n, nnb, nbh,                    &
		energy                          &
	)
	
	real, intent(in) :: T, norm_B, norm_k
	real, intent(in), dimension(3) :: dir_B, ek
	real, intent(in), dimension(:), allocatable :: smp
	real, intent(inout), dimension(:,:), allocatable :: s
	integer, intent(in) :: mcs, n
	integer, intent(in), dimension(:,:), allocatable :: nbh
	integer, intent(in), dimension(:), allocatable :: nnb
	
	real, intent (out), dimension(:), allocatable :: energy
	
	real :: dE, E1, E2, Et, p
	real, dimension(3) :: magt, mag, st
	integer :: i, j, k
	
	allocate (energy(mcs))
	
	do k = 1, mcs
        Et = 0.0
        magt = 0.0
        
        do i = 1, n
            call cheap_rand_changed_vec (st, s(i,:))
            E1 = - dot_product (norm_B * dir_B, s(i,:))
            E2 = - dot_product (norm_B * dir_B, st)
            
            E1 = - norm_k * dot_product (ek, s(i,:)) ** 2
            E2 = - norm_k * dot_product (ek, st) ** 2
            
            do j = 1, nnb(i)
                E1 = E1 - (smp(i) + smp(j)) * dot_product (s(i,:), s(nbh(i,j),:))
                E2 = E2 - (smp(i) + smp(j)) * dot_product (st,     s(nbh(i,j),:))
            end do
            
            dE = E2 - E1
            if (dE .le. 0) then
                s(i,:) = st
                Et = Et + E2
            else
                call random_number (p)
                if (p .le. exp(- dE / T)) then
                    s(i,:) = st
                    Et = Et + E2
                else
                    Et = Et + E1
                end if
            end if
            magt = magt + s(i,:)
        end do

        
        energy(k) = Et
        mag = mag + (magt / n)
        
    end do
	
	end subroutine estabilize
	
	subroutine measure(                 &
	    T, norm_B, dir_B, mcs,          &
	    s, smp,                         &
	    ek, norm_k,                     &
	    n, nnb, nbh,                    &
	    energy, mag                     &
	)
	
	real, intent(in) :: T, norm_B, norm_k
	real, intent(in), dimension(3) :: dir_B, ek
	real, intent(in), dimension(:), allocatable :: smp
	real, intent(inout), dimension(:,:), allocatable :: s
	integer, intent(in) :: mcs, n
	integer, intent(in), dimension(:,:), allocatable :: nbh
	integer, intent(in), dimension(:), allocatable :: nnb
	
	real, intent (out), dimension(:), allocatable :: energy
	real, intent (out), dimension(:,:), allocatable :: mag
	
	real :: dE, E1, E2, Et, p
	real, dimension(3) :: magt, st
	integer :: i, j, k
	
	allocate (energy(mcs), mag(mcs, 3))
	
	do k = 1, mcs
        Et = 0.0
        magt = 0.0
        
        do i = 1, n
            call cheap_rand_changed_vec (st, s(i,:))
            E1 = - dot_product (norm_B * dir_B, s(i,:))
            E2 = - dot_product (norm_B * dir_B, st)
            
            E1 = E1 - norm_k * dot_product (ek, s(i,:)) ** 2
            E2 = E2 - norm_k * dot_product (ek, st) ** 2
            
            do j = 1, nnb(i)
                E1 = E1 - (smp(i) + smp(j)) * dot_product (s(i,:), s(nbh(i,j),:))
                E2 = E2 - (smp(i) + smp(j)) * dot_product (st,     s(nbh(i,j),:))
            end do
            
            dE = E2 - E1
            if (dE .le. 0) then
                s(i,:) = st
                Et = Et + E2
            else
                call random_number (p)
                if (p .le. exp(- dE / T)) then
                    s(i,:) = st
                    Et = Et + E2
                else
                    Et = Et + E1
                end if
            end if
            magt = magt + s(i,:)
        end do

        
        energy(k) = Et
        mag(k,:) = magt
        
    end do
	
	end subroutine measure
	
	

    subroutine mag_curves (     &
        T_max, T_min, dT,       &
        n, nnb, nbh,            &
        J_ex, norm_s,           &
        norm_B, dir_B,          &
        mcs_max, mcs_c, k_B,    &
        fname                   &
    )
    
    real, intent(in) :: T_max, T_min, dT
    
    integer, intent(in) :: n
    integer, dimension(:,:), allocatable, intent(in) :: nbh
    integer, dimension(:), allocatable, intent(in) :: nnb
    
    real, intent(in) :: J_ex, norm_s
    
    real, intent(in) :: norm_B
    real, dimension(3), intent(in) :: dir_B
    
    integer, intent(in) :: mcs_max, mcs_c
    real, intent(in) :: k_B
    
    character (len = 50) :: fname
    
    integer :: i, j, k
    real :: E, Et, E1, E2, dE
    real, dimension(3) :: mag, magt, st
    real, dimension(:,:), allocatable :: s
    real :: p, T
    
    call init_random_seed()
    
    allocate (s(n,3))
    do i = 1, n
        call rdn_vec (s(i,:), norm_s)
    end do

    open(1, file = trim(fname), status = 'new')
    
    do T = T_max, T_min, - dT   
        E = 0.0
        mag = 0.0
        
        do k = 1, mcs_max
            Et = 0.0
            magt = 0.0
            
            do i = 1, n
                call cheap_rand_changed_vec (st, s(i,:))
                E1 = - dot_product (norm_B * dir_B, s(i,:))
                E2 = - dot_product (norm_B * dir_B, st)
                
                do j = 1, nnb(i)
                    E1 = E1 - J_ex * dot_product (s(i,:), s(nbh(i,j),:))
                    E2 = E2 - J_ex * dot_product (st,     s(nbh(i,j),:))
                end do
                
                dE = E2 - E1
                if (dE .le. 0) then
                    s(i,:) = st
                    Et = Et + E2
                else
                    call random_number (p)
                    if (p .le. exp(- dE / (k_B * T))) then
                        s(i,:) = st
                        Et = Et + E2
                    else
                        Et = Et + E1
                    end if
                end if
                magt = magt + s(i,:)
            end do

            if (k > mcs_c) then
                E = E + Et / n
                mag = mag + (magt / n)
            end if
        end do

        write (1,*) T, mag/(mcs_max - mcs_c), E/(mcs_max - mcs_c)
        write (*,*) T, mag/(mcs_max - mcs_c), E/(mcs_max - mcs_c)
    end do
    
    close (1)
    end subroutine mag_curves

    subroutine hyst_loop (                          &
        T,                                          &
        n, nnb, nbh,                                &
        J_ex, norm_s,                               &
        norm_B_max, dnorm_B, dir_B,                 &
        mcs_max, mcs_c, k_B,                        &
        fname                                       &
    )
    
    real, intent(in) :: T
    
    integer, intent(in) :: n
    integer, dimension(:,:), allocatable, intent(in) :: nbh
    integer, dimension(:), allocatable, intent(in) :: nnb
    
    real, intent(in) :: J_ex, norm_s
    
    real, intent(in) :: norm_B_max
    real :: dnorm_B
    real, dimension(3), intent(in) :: dir_B
    
    integer, intent(in) :: mcs_max, mcs_c
    real, intent(in) :: k_B
    
    character (len = 50) :: fname
    
    integer :: i, j, k, id = 0
    real :: E, Et, E1, E2, dE
    real, dimension(3) :: mag, magt, st
    real, dimension(:,:), allocatable :: s
    real :: p, norm_B
    
    call init_random_seed()
    
    allocate (s(n,3))
    do i = 1, n
        call rdn_vec (s(i,:), norm_s)
    end do

    open(1, file = trim(fname), status = 'new')
 
    do while (id .le. 2)   
        E = 0.0
        mag = 0.0
        
        do k = 1, mcs_max
            Et = 0.0
            magt = 0.0
            
            do i = 1, n
                call cheap_rand_changed_vec (st, s(i,:))
                E1 = - dot_product (norm_B * dir_B, s(i,:))
                E2 = - dot_product (norm_B * dir_B, st)
                
                do j = 1, nnb(i)
                    E1 = E1 - J_ex * dot_product (s(i,:), s(nbh(i,j),:))
                    E2 = E2 - J_ex * dot_product (st,     s(nbh(i,j),:))
                end do
                
                dE = E2 - E1
                if (dE .le. 0) then
                    s(i,:) = st
                    Et = Et + E2
                else
                    call random_number (p)
                    if (p .le. exp(- dE / (k_B * T))) then
                        s(i,:) = st
                        Et = Et + E2
                    else
                        Et = Et + E1
                    end if
                end if
                magt = magt + s(i,:)
            end do

            if (k > mcs_c) then
                E = E + Et / n
                mag = mag + (magt / n)
            end if
        end do

        write (1,*) T, mag/(mcs_max - mcs_c), E/(mcs_max - mcs_c)
        write (*,*) T, mag/(mcs_max - mcs_c), E/(mcs_max - mcs_c)
        
        if (abs(norm_B) .gt. norm_B_max) then
            dnorm_B = - dnorm_B
            id = id + 1
        end if
        
        norm_B = norm_B + dnorm_B
        
    end do
    
    close (1)
    end subroutine hyst_loop

end module core_utils
