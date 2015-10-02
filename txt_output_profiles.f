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
! Subroutine name: txt_output_profiles
! Subroutine description: to write the output profiles in ".txt" file 
! format. 
!-----------------------------------------------------------------------
      subroutine txt_output_profiles(conc,sigc,c_skew,c_kurt,
     &conc_react1,conc_react2,sigc_react1,c_skew_react1,c_kurt_react1,
     &sigc_react2,c_skew_react2,c_kurt_react2,conc_product,sigc_product,
     &c_skew_product,c_kurt_product,Iseg,sigz,x_profiles,y_profiles,
     &z_profiles,dx,dy,dz,C0_react1,C0_react2,reactive,ndx,ndy,ndz,
     &n_profiles)
!-----------------------------------------------------------------------
! Declarations
!-----------------------------------------------------------------------
      integer :: reactive,ndx,ndy,ndz,n_profiles
      real :: dx,dy
      double precision :: dz,C0_react1,C0_react2
      double precision :: sigz(ndx),x_profiles(n_profiles)
      double precision :: y_profiles(n_profiles),z_profiles(n_profiles)
      double precision :: conc(ndx,ndy,ndz),sigc(ndx,ndy,ndz)
      double precision :: c_skew(ndx,ndy,ndz),c_kurt(ndx,ndy,ndz)
      double precision :: conc_react1(ndx,ndy,ndz)
      double precision :: conc_react2(ndx,ndy,ndz)
      double precision :: sigc_react1(ndx,ndy,ndz)
      double precision :: c_skew_react1(ndx,ndy,ndz)
      double precision :: c_kurt_react1(ndx,ndy,ndz)
      double precision :: sigc_react2(ndx,ndy,ndz)
      double precision :: c_skew_react2(ndx,ndy,ndz)
      double precision :: c_kurt_react2(ndx,ndy,ndz)
      double precision :: conc_product(ndx,ndy,ndz)
      double precision :: sigc_product(ndx,ndy,ndz)
      double precision :: c_skew_product(ndx,ndy,ndz)
      double precision :: c_kurt_product(ndx,ndy,ndz),Iseg(ndx,ndy,ndz)
! Local variables      
      integer :: i_prof,ix,iy,iz,i,j,k
      double precision :: x,y,z,sigmaz,Cm_r1_FL,Cm_r2_FL,sigmaC_r1_FL
      double precision :: sigmaC_r2_FL,SkC_r2_FL,Cm_r1,Cm_r2,Cm_pr
      double precision :: sigmaC_r1,sigmaC_r2,sigmaC_pr,SkC_r1,SkC_r2 
      double precision :: SkC_pr,KuC_r1,KuC_r2,KuC_pr,Is
!-----------------------------------------------------------------------
! Initializations
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Statements
!-----------------------------------------------------------------------
! Writing output profiles for concentration (CA,CB,CC) statistics, Iseg 
! and sigz    
      open(11,file='profiles.txt')
      do i_prof=1,n_profiles
         write(11,'(/,a,i6)') "profile: ",i_prof
         write(11,'(3(14x,a),(9x,a),(13x,a),2(7x,a),3(10x,a),
     &      (9x,a),2(3x,a),3(6x,a),(12x,a),(6x,a),3(9x,a),(12x,a),
     &      3(9x,a),(11x,a))') "x","y","z","sigmaz","Cm","Cm_r1_FL",
     &      "Cm_r2_FL","Cm_r1","Cm_r2","Cm_pr","sigmaC",
     &      "sigmaC_r1_FL","sigmaC_r2_FL","sigmaC_r1","sigmaC_r2",
     &      "sigmaC_pr","SkC","SkC_r2_FL","SkC_r1","SkC_r2","SkC_pr",
     &      "KuC","KuC_r1","KuC_r2","KuC_pr","Iseg"
         do i=1,ndx
            if (x_profiles(i_prof).ne.-999.d0) then
               ix = int(x_profiles(i_prof) / dx) + 1  
               x = x_profiles(i_prof)
               else
                  ix = i
                  x = dx * (ix - 1 + 0.5)
            endif
            do j=1,ndy
               if (y_profiles(i_prof).ne.-999.d0) then
                  iy = int(y_profiles(i_prof) / dy) + 1  
                  y = y_profiles(i_prof)
                  else
                     iy = j
                     y = dy * (iy - 1 + 0.5)
               endif
               do k=1,ndz
                  if (z_profiles(i_prof).ne.-999.d0) then
                     iz = int(z_profiles(i_prof) / dz) + 1 
                     z = z_profiles(i_prof)
                     else
                        iz = k
                        z = dz * (iz - 1 + 0.5)
                  endif
                  if (x_profiles(i_prof).eq.-999.d0) then
                     sigmaz = sigz(ix)
                     else
                        sigmaz = -999.d0
                  endif
                  if (reactive==1) then
                     Cm_r1_FL = conc(ix,iy,iz) * C0_react1
                     Cm_r2_FL = (1.d0 - conc(ix,iy,iz)) * C0_react2
                     sigmaC_r1_FL = sigc(ix,iy,iz) * C0_react1
                     sigmaC_r2_FL = sigc(ix,iy,iz) * C0_react2
                     SkC_r2_FL = - c_skew(ix,iy,iz)
                     Cm_r1 = conc_react1(ix,iy,iz)
                     Cm_r2 = conc_react2(ix,iy,iz)
                     Cm_pr = conc_product(ix,iy,iz)
                     sigmaC_r1 = sigc_react1(ix,iy,iz)
                     sigmaC_r2 = sigc_react2(ix,iy,iz)
                     sigmaC_pr = sigc_product(ix,iy,iz)
                     SkC_r1 = c_skew_react1(ix,iy,iz)
                     SkC_r2 = c_skew_react2(ix,iy,iz)
                     SkC_pr = c_skew_product(ix,iy,iz)
                     KuC_r1 = c_kurt_react1(ix,iy,iz)
                     KuC_r2 = c_kurt_react2(ix,iy,iz)
                     KuC_pr = c_kurt_product(ix,iy,iz)
                     Is = Iseg(ix,iy,iz)
                     else
                        conc_react1_FL = -999999.d0
                        conc_react2_FL = -999999.d0
                        sigc_react1_FL = -999999.d0
                        sigc_react2_FL = -999999.d0
                        c_skew_react2_FL = -999999.d0
                        Cm_r1 = -999999.d0
                        Cm_r2 = -999999.d0
                        Cm_pr = -999999.d0
                        sigmaC_r1 = -999999.d0
                        sigmaC_r2 = -999999.d0
                        sigmaC_pr = -999999.d0
                        SkC_r1 = -999999.d0
                        SkC_r2 = -999999.d0
                        SkC_pr = -999999.d0
                        KuC_r1 = -999999.d0
                        KuC_r2 = -999999.d0
                        KuC_pr = -999999.d0
                        Is = -999999.d0
                  endif
                  write(11,'(27(f15.5))') x,y,z,sigmaz,
     &               conc(ix,iy,iz),Cm_r1_FL,Cm_r2_FL,Cm_r1,Cm_r2,Cm_pr,
     &               sigc(ix,iy,iz),sigmaC_r1_FL,sigmaC_r2_FL,sigmaC_r1,
     &               sigmaC_r2,sigmaC_pr,c_skew(ix,iy,iz),SkC_r2_FL,
     &               SkC_r1,SkC_r2,SkC_pr,c_kurt(ix,iy,iz),KuC_r1,
     &               KuC_r2,KuC_pr,Is
                  if (z_profiles(i_prof).ne.-999.d0) exit
               enddo
               if (y_profiles(i_prof).ne.-999.d0) exit
            enddo
            if (x_profiles(i_prof).ne.-999.d0) exit
         enddo
      enddo       
      close(11)
      return
      end 

