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
! Subroutine name: BC_lateral_boundaries
! Subroutine description: to impose wall Boundary Conditions at lateral 
! boundaries, if required in input.
!-----------------------------------------------------------------------
      subroutine BC_lateral_boundaries(x_min,dx,x_max,y_min,dy,y_max,
     &BC_lat_x_min,BC_lat_x_max,BC_lat_y_min,BC_lat_y_max,X,U_prime,Y,
     &V_prime) 
!----------------------------------------------------------------------
! Declarations
!----------------------------------------------------------------------
      integer(4) :: BC_lat_x_min,BC_lat_x_max,BC_lat_y_min,BC_lat_y_max
      real x_min,dx,x_max,y_min,dy,y_max
      real X,U_prime,Y,V_prime
!----------------------------------------------------------------------
! Initializations
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Statements
!-----------------------------------------------------------------------
      if ((BC_lat_x_min==1).and.(X.le.x_min)) then
         X = -X
         if (X==0.) X = 0.001 * dx
         U_prime = -U_prime
      endif
      if ((BC_lat_x_max==1).and.(X.ge.x_max)) then
         X = 2. * x_max - X
         if (X==x_max) X = x_max - 0.001 * dx
         U_prime = -U_prime
      endif
      if ((BC_lat_y_min==1).and.(Y.le.y_min)) then
         Y = -Y
         if (Y==0.) Y = 0.001 * dy
         V_prime = -V_prime
      endif
      if ((BC_lat_y_max==1).and.(Y.ge.y_max)) then
         Y = 2. * y_max - Y
         if (Y==y_max) Y = y_max - 0.001 * dy
         V_prime = -V_prime
      endif
      return
      end

