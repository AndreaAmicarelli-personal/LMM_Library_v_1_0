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
! Subroutine name: BC_C_Dirichlet
! Subroutine description: to impose Dirichlet Boundary Conditions at the
! selected domain frontiers. If (phase==1), then the concentration means 
! and conditonal means are modified. If (phase==2), then the particle 
! instantaneous concentration is modified.
!-----------------------------------------------------------------------
      subroutine BC_C_Dirichlet(phase,Cmc,icel,jcel,kcel,
     &ndx,ndy,ndz,nbu,nbv,nbw,Cm,BC_C_Dirichlet_flag_x_min,
     &BC_C_Dirichlet_flag_x_max,BC_C_Dirichlet_flag_y_min,
     &BC_C_Dirichlet_flag_y_max,BC_C_Dirichlet_flag_z_min,
     &BC_C_Dirichlet_flag_z_max,BC_C_Dirichlet_value_x_min,
     &BC_C_Dirichlet_value_x_max,BC_C_Dirichlet_value_y_min,
     &BC_C_Dirichlet_value_y_max,BC_C_Dirichlet_value_z_min,
     &BC_C_Dirichlet_value_z_max,Concentration)
!----------------------------------------------------------------------
! Declarations
!----------------------------------------------------------------------
      logical :: BC_C_Dirichlet_flag_x_min,BC_C_Dirichlet_flag_x_max
      logical :: BC_C_Dirichlet_flag_y_min,BC_C_Dirichlet_flag_y_max
      logical :: BC_C_Dirichlet_flag_z_min,BC_C_Dirichlet_flag_z_max
      integer :: phase,icel,jcel,kcel,ndx,ndy,ndz,nbu,nbv,nbw
      double precision :: BC_C_Dirichlet_value_x_min
      double precision :: BC_C_Dirichlet_value_x_max
      double precision :: BC_C_Dirichlet_value_y_min
      double precision :: BC_C_Dirichlet_value_y_max
      double precision :: BC_C_Dirichlet_value_z_min
      double precision :: BC_C_Dirichlet_value_z_max
      double precision :: Cm(ndx,ndy,ndz)
      double precision :: Cmc(ndx,ndy,ndz,nbu,nbv,nbw)
      double precision :: Concentration
      integer :: i,j,k,i2,j2,k2
!----------------------------------------------------------------------
! Initializations
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Statements
!-----------------------------------------------------------------------
      if (phase.eq.1) then
         do k=1,ndz
            do j=1,ndy
               if (BC_C_Dirichlet_flag_x_min.eqv.(.true.)) then
                  Cm(1,j,k) = BC_C_Dirichlet_value_x_min
               endif
               if (BC_C_Dirichlet_flag_x_max.eqv.(.true.)) then
                  Cm(ndx,j,k) = BC_C_Dirichlet_value_x_max
               endif
               do k2=1,nbw
                  do j2=1,nbv
                     do i2=1,nbu
                        if (BC_C_Dirichlet_flag_x_min.eqv.(.true.)) then
                           Cmc(1,j,k,i2,j2,k2) = 
     &                        BC_C_Dirichlet_value_x_min
                        endif
                        if (BC_C_Dirichlet_flag_x_max.eqv.(.true.)) then
                           Cmc(ndx,j,k,i2,j2,k2) = 
     &                        BC_C_Dirichlet_value_x_max
                        endif
                     enddo
                  enddo
               enddo      
            enddo
         enddo
         do k=1,ndz
            do i=1,ndx
               if (BC_C_Dirichlet_flag_y_min.eqv.(.true.)) then
                  Cm(i,1,k) = BC_C_Dirichlet_value_y_min
               endif
               if (BC_C_Dirichlet_flag_y_max.eqv.(.true.)) then
                  Cm(i,ndy,k) = BC_C_Dirichlet_value_y_max
               endif
               do k2=1,nbw
                  do j2=1,nbv
                     do i2=1,nbu
                        if (BC_C_Dirichlet_flag_y_min.eqv.(.true.)) then
                           Cmc(i,1,k,i2,j2,k2) = 
     &                        BC_C_Dirichlet_value_y_min
                        endif
                        if (BC_C_Dirichlet_flag_y_max.eqv.(.true.)) then
                           Cmc(i,ndy,k,i2,j2,k2) = 
     &                        BC_C_Dirichlet_value_y_max
                        endif
                     enddo
                  enddo
               enddo      
            enddo
         enddo
         do j=1,ndy
            do i=1,ndx
               if (BC_C_Dirichlet_flag_z_min.eqv.(.true.)) then
                  Cm(i,j,1) = BC_C_Dirichlet_value_z_min
               endif
               if (BC_C_Dirichlet_flag_z_max.eqv.(.true.)) then
                  Cm(i,j,ndz) = BC_C_Dirichlet_value_z_max
               endif
               do k2=1,nbw
                  do j2=1,nbv
                     do i2=1,nbu
                        if (BC_C_Dirichlet_flag_z_min.eqv.(.true.)) then
                           Cmc(i,j,1,i2,j2,k2) = 
     &                        BC_C_Dirichlet_value_z_min
                        endif
                        if (BC_C_Dirichlet_flag_z_max.eqv.(.true.)) then
                           Cmc(i,j,ndz,i2,j2,k2) = 
     &                        BC_C_Dirichlet_value_z_max
                        endif
                     enddo
                  enddo
               enddo      
            enddo
         enddo
         elseif (phase.eq.2) then
            if ((icel.eq.1).and.(BC_C_Dirichlet_flag_x_min.eqv.(.true.))
     &         ) then
               Concentration = BC_C_Dirichlet_value_x_min
            endif
            if ((icel.eq.ndx).and.(BC_C_Dirichlet_flag_x_max.eqv.(.true.
     &         ))) then
               Concentration = BC_C_Dirichlet_value_x_max
            endif
            if ((jcel.eq.1).and.(BC_C_Dirichlet_flag_y_min.eqv.(.true.))
     &         ) then
               Concentration = BC_C_Dirichlet_value_y_min
            endif
            if ((jcel.eq.ndy).and.(BC_C_Dirichlet_flag_y_max.eqv.(.true.
     &         ))) then
               Concentration = BC_C_Dirichlet_value_y_max
            endif
            if ((kcel.eq.1).and.(BC_C_Dirichlet_flag_z_min.eqv.(.true.))
     &         ) then
               Concentration = BC_C_Dirichlet_value_z_min
            endif
            if ((kcel.eq.ndz).and.(BC_C_Dirichlet_flag_z_max.eqv.(.true.
     &         ))) then
               Concentration = BC_C_Dirichlet_value_z_max
            endif
      endif
      return
      end

