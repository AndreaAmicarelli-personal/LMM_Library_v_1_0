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
! Subroutine name: rect_source_100_xy_part
! Subroutine description: provided a particle ID, it initializes y- and 
! z-coordinate of the particle position (ypr,zpr) for a rectangular 
! source, whose normal is aligned with x-axis (conservative phase of a 
! LMM stationary code run).
!-----------------------------------------------------------------------
      subroutine rect_source_100_xy_part(i_part,dy,dz,ypr,zpr,
     &rect_source_100_mat,np1_per_cell,ndy,ndz)
!-----------------------------------------------------------------------
! Declarations
!-----------------------------------------------------------------------
      integer :: i_part,ndy,ndz
      integer(4) :: np1_per_cell
      real :: dy,ypr,zpr
      double precision :: dz
      integer(4) :: rect_source_100_mat(ndy,ndz)
! Local variables
      integer :: i_cel_aux,i_cel_part
!-----------------------------------------------------------------------
! Initializations
!-----------------------------------------------------------------------
      i_cel_aux = 0
!-----------------------------------------------------------------------
! Statements
!-----------------------------------------------------------------------
      i_cel_part = int((i_part - 1) / np1_per_cell) + 1
      do j=1,ndz
         do i=1,ndy
            if (rect_source_100_mat(i,j)==1) then
               i_cel_aux = i_cel_aux + 1
               if (i_cel_part==i_cel_aux) then
                  ypr = (i - 0.5) * dy
                  zpr = (j - 0.5) * dz
! In Fortran 77 one cannot name a cycle
                  exit 
               endif
            endif
         enddo
! In Fortran 77 one cannot name a cycle
         if (i_cel_part==i_cel_aux) then
            exit
         endif
      enddo
      return
      end

