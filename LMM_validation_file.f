!-----------------------------------------------------------------------
! "LMM Library" is a Fortran library, which may be useful for a generic 
! stochastic Lagrangian Micro-Mixing (LMM) numerical model for 
! pollutant/scalar dispersion. 
! LMM Library Copyright 2008-2015 Andrea Amicarelli
! email contact: Andrea.Amicarelli@gmail.com
!-----------------------------------------------------------------------
! This file is part of LMM Library. 
! LMM_library is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as 
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
! Subroutine name: LMM_validation_file 
! Subroutine description: to compute the validation statistics on scalar
! Renolds' statistics modelling, provided measures (input file) vs. 
! numerical estimations of the associated LMM model. 
!-----------------------------------------------------------------------
      subroutine LMM_validation_file(conc,sigc,c_skew,c_kurt,
     &conc_react1,conc_react2,dx,dy,dz,reactive,ndx,ndy,ndz,nmeas)
!-----------------------------------------------------------------------
! Declarations
!-----------------------------------------------------------------------
      integer :: reactive,ndx,ndy,ndz,nmeas
      real :: dx,dy
      double precision :: dz
      double precision :: conc(ndx,ndy,ndz),sigc(ndx,ndy,ndz)
      double precision :: c_skew(ndx,ndy,ndz),c_kurt(ndx,ndy,ndz)
      double precision :: conc_react1(ndx,ndy,ndz)
      double precision :: conc_react2(ndx,ndy,ndz)
! Local variables
      integer :: ix,iy,iz
      real :: test(nmeas,11)
!-----------------------------------------------------------------------
! Initializations
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Statements
!-----------------------------------------------------------------------
      open(21,file='Cmeas.txt')
1114  format(/)
1112  format(7(F12.6))
      read(21,1114)
! Measured values are read from input file (just for validation 
! purposes) and written in the array "test".
      do i=1,nmeas
         read(21,1112,IOSTAT=IOstatus) (test(i,j),j=1,7)
      enddo
      close (21)
      do i=1,nmeas
         ix = test(i,1) / dx + 1
         iy = test(i,2) / dy + 1
         iz = test(i,3) / dz + 1
         if (ix.gt.ndx) ix = ndx
! Simulated values in the array "test" 
         if (reactive==0) then
            test(i,8) = conc(ix,iy,iz)
            test(i,9) = sigc(ix,iy,iz)
            test(i,10) = c_skew(ix,iy,iz)
            test(i,11) = c_kurt(ix,iy,iz)
            else
               test(i,8) = conc(ix,iy,iz)
               test(i,9) = conc_react1(ix,iy,iz)
               test(i,10) = conc_react2(ix,iy,iz)
               test(i,11) = conc_react2(ix,iy,iz)
         endif
      enddo
      open(22,file="test")
1113  format(11(F12.6))
1115  format(11a12)
      if (reactive==0) then
         write(22,1115) 'X(m)','Y(m)','Z(m)','Cmmeas','sigCmeas'
     &      ,'c_skew_meas','c_kurt_meas','Cmsim','sigCsim','c_skew_sim',
     &      'c_kurt_sim'
         else
            write(22,1115) 'X(m)','Y(m)','Z(m)','Cm_Fm_meas',
     &         'Cm_r1_meas','Cm_r2_meas','Cm_r2_meas','Cm_Fm_sim',
     &         'Cm_r1_sim','Cm_r2__sim','Cm_r2__sim'
      endif
      do i=1,nmeas
         write(22,1113) (test(i,j),j=1,11)
      enddo
      close(22)
      return
      end

