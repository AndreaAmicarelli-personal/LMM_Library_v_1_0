!-----------------------------------------------------------------------
! "LMM Library v.1.0" is a Fortran library, which may be useful for a  
! generic stochastic Lagrangian Micro-Mixing (LMM) numerical model for 
! pollutant/scalar dispersion. 
! LMM Library v.1.0 Copyright 2008-2015 Andrea Amicarelli
! email contact: Andrea.Amicarelli@gmail.com
!-----------------------------------------------------------------------
! This file is part of LMM Library v.1.0. 
! LMM Library v.1.0 is free software: you can redistribute it and/or 
! modify it under the terms of the GNU Lesser General Public License as 
! published by the Free Software Foundation, either version 3 of the 
! License, or (at your option) any later version.
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU Lesser General Public License for more details.
! You should have received a copy of the GNU Lesser General Public 
! License along with this library. If not, see 
! <http://www.gnu.org/licenses/>.
!-----------------------------------------------------------------------
! Subroutine name: rect_source_100
! Subroutine description: to compute the matrix "quad_source_100_mat",  
! which indicates the presence/absence of a rectangular source 
! (whose normal is aligned with x-axis) in a domain cell, provided the
! source minima and maxima values of y- and z-coordinate. To  
! compute the number of particles per source cell ("np1_per_cell"). To 
! update the total number of particles in the conservative phase of a 
! LMM stationary code run ("np" or "np1") to have a uniform initial 
! distribution of particles within the source section. 
!-----------------------------------------------------------------------
      subroutine rect_source_100(dy,dz,np,rect_source_100_mat,
     &y_qs_100_min,y_qs_100_max,z_qs_100_min,z_qs_100_max,n_cel_qs_100,
     &ndy,ndz,np1_per_cell)
!-----------------------------------------------------------------------
! Declarations
!-----------------------------------------------------------------------
      integer :: np,ndy,ndz
      integer(4) :: np1_per_cell,n_cel_qs_100
      real :: dy
      double precision :: dz,y_qs_100_min,y_qs_100_max,z_qs_100_min
      double precision :: z_qs_100_max
      integer(4) :: rect_source_100_mat(ndy,ndz)
! Local variables           
      double precision :: y_cel,z_cel
!-----------------------------------------------------------------------
! Initializations
!-----------------------------------------------------------------------
      n_cel_qs_100 = 0
!-----------------------------------------------------------------------
! Statements
!-----------------------------------------------------------------------
      do j=1,ndz 
         z_cel = (j - 0.5) * dz 
         if ((z_cel>=z_qs_100_min).and.(z_cel<=z_qs_100_max)) then
            do i=1,ndy  
               y_cel = (i - 0.5) * dy
               if ((y_cel>=y_qs_100_min).and.(y_cel<=y_qs_100_max)) then
                  rect_source_100_mat(i,j) = 1
                  n_cel_qs_100 =  n_cel_qs_100 + 1
               endif
            enddo
         endif
      enddo
      np = np - mod(np,n_cel_qs_100)
      np1_per_cell = np / n_cel_qs_100
      return
      end 

