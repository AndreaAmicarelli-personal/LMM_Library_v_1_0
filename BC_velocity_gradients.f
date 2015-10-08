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
! Subroutine name: BC_velocity_gradients
! Subroutine description: to zero velocit gradients at wall frontiers.
!-----------------------------------------------------------------------
      subroutine BC_velocity_gradients(BC_lat_x_min,BC_lat_x_max,
     &BC_lat_y_min,BC_lat_y_max,ndx,ndy,ndz,derxvel,derxvar)
!----------------------------------------------------------------------
! Declarations
!----------------------------------------------------------------------
      integer(4) :: BC_lat_x_min,BC_lat_x_max,BC_lat_y_min,BC_lat_y_max
      integer :: ndx,ndy,ndz
      real derxvel(ndx,ndy,ndz,3)
      real derxvar(ndx,ndy,ndz,3,3)
      integer :: i,j,k
!-----------------------------------------------------------------------
! Initializations
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Statements
!-----------------------------------------------------------------------
      if (BC_lat_x_min==1) then
         do k=1,ndz
            do j=1,ndy
               derxvel(1,j,k,1) = 0.
               derxvar(1,j,k,1,1) = 0.
               derxvar(1,j,k,2,1) = 0.
               derxvar(1,j,k,3,1) = 0.
            enddo
         enddo   
      endif
      if (BC_lat_x_max==1) then
         do k=1,ndz
            do j=1,ndy
               derxvel(ndx,j,k,1) = 0.
               derxvar(ndx,j,k,1,1) = 0.
               derxvar(ndx,j,k,2,1) = 0.
               derxvar(ndx,j,k,3,1) = 0.
            enddo
         enddo
      endif
      if (BC_lat_y_min==1) then
         do k=1,ndz
            do i=1,ndx
               derxvel(i,1,k,2) = 0.
               derxvar(i,1,k,1,2) = 0.
               derxvar(i,1,k,2,2) = 0.
               derxvar(i,1,k,3,2) = 0.
            enddo
         enddo  
      endif
      if (BC_lat_y_max==1) then
         do k=1,ndz
            do i=1,ndx
               derxvel(i,ndy,k,2) = 0.
               derxvar(i,ndy,k,1,2) = 0.
               derxvar(i,ndy,k,2,2) = 0.
               derxvar(i,ndy,k,3,2) = 0.
            enddo
         enddo
      endif
      do j=1,ndy
         do i=1,ndx
            derxvel(i,j,1,3) = 0.
            derxvar(i,j,1,1,3) = 0.
            derxvar(i,j,1,2,3) = 0.
            derxvar(i,j,1,3,3) = 0.
            derxvel(i,j,ndz,3) = 0.
            derxvar(i,j,ndz,1,3) = 0.
            derxvar(i,j,ndz,2,3) = 0.
            derxvar(i,j,ndz,3,3) = 0.
         enddo
      enddo
      return
      end

