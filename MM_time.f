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
! Subroutine name: MM_time
! Subroutine description: to compute the Micro-Mixing time, depending on 
! the chosen formulation:
! MM_time_flag=1: Equation (42) of Amicarelli et al. (2012,AmeJouEnvSci) 
!-----------------------------------------------------------------------
      subroutine MM_time(mu_MM,sigma_z0,eps,ratio_CR_C0,C0,t_fly,
     &MM_time_flag,tm) 
!-----------------------------------------------------------------------
! Declarations
!-----------------------------------------------------------------------
      integer :: MM_time_flag
      real :: mu_MM,eps,C0
      double precision :: tm,sigma_z0,t_fly,ratio_CR_C0
      double precision :: CR
!-----------------------------------------------------------------------
! Initializations
!-----------------------------------------------------------------------
      CR = C0 * ratio_CR_C0
!-----------------------------------------------------------------------
! Statements
!-----------------------------------------------------------------------
      if (MM_time_flag.eq.1) then
         tm = mu_MM * ((3.d0/2.d0) ** (0.5d0) * sigma_z0 ** (2.d0 / 
     &        3.d0) / eps ** (1.d0/3.d0) + dsqrt(3.d0/2.d0) * CR **  
     &        (1.d0/3.d0) * t_fly)
      endif
      return
      end

