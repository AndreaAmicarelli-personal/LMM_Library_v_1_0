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
! Subroutine name: txt_output_fields
! Subroutine description: to write the output fields in ".txt" file 
! format. 
!-----------------------------------------------------------------------
      subroutine txt_output_fields(conc,sigc,c_skew,c_kurt,conc_react1,
     &conc_react2,sigc_react1,c_skew_react1,c_kurt_react1,sigc_react2,
     &c_skew_react2,c_kurt_react2,conc_product,sigc_product,
     &c_skew_product,c_kurt_product,Iseg,sigz,u,v,w,varu,varv,varw,eps,
     &dx,dy,dz,reactive,ndx,ndy,ndz)
!-----------------------------------------------------------------------
! Declarations
!-----------------------------------------------------------------------
      integer :: reactive,ndx,ndy,ndz
      real :: dx,dy
      double precision :: dz
      double precision :: sigz(ndx)
      real :: u(ndx,ndy,ndz),v(ndx,ndy,ndz),w(ndx,ndy,ndz)
      real :: varu(ndx,ndy,ndz),varv(ndx,ndy,ndz),varw(ndx,ndy,ndz)
      real :: eps(ndx,ndy,ndz)
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
      double precision :: x,y,z
!-----------------------------------------------------------------------
! Initializations
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Statements
!-----------------------------------------------------------------------
! Writing output for concentration statistics (for Fm,CA,CB,CC), Iseg 
! and sigz    
      open(11,file='Cm.txt')
      open(14,file='sigz.txt')
      open(22,file='sigmaC.txt')
      open(23,file='Sk_C.txt')
      open(24,file='Ku_C.txt')
      if (reactive==1) then
         open(12,file='Cm_react1.txt')
         open(13,file='Cm_react2.txt')
         open(25,file='Cm_product.txt')
         open(15,file='sigma_C_react1.txt')
         open(16,file='sigma_C_react2.txt')
         open(26,file='sigma_C_product.txt')
         open(17,file='Sk_C_react1.txt')
         open(18,file='Sk_C_react2.txt')
         open(27,file='Sk_C_product.txt')
         open(19,file='Ku_C_react1.txt')
         open(20,file='Ku_C_react2.txt')
         open(28,file='Ku_C_product.txt')
         open(21,file='Iseg.txt')
! AA!!! test: start
         open(31,file='um.txt')
         open(32,file='vm.txt')
         open(33,file='wm.txt')
         open(34,file='varu.txt')
         open(35,file='varv.txt')
         open(36,file='varw.txt')
         open(37,file='eps.txt')
! AA!!! test: end
      endif
      do k=1,ndz
         z = dz * (k - 1 + 0.5d0)
         do j=1,ndy
            y = dy * (j - 1 + 0.5d0)
            do i=1,ndx
               x = dx * (i - 1 + 0.5d0)
               write(11,'(4(f15.5))') x,y,z,conc(i,j,k)
               write(22,'(4(f15.5))') x,y,z,sigc(i,j,k)
               write(23,'(4(f15.5))') x,y,z,c_skew(i,j,k)
               write(24,'(4(f15.5))') x,y,z,c_kurt(i,j,k)
               if (reactive==1) then
                  write(12,'(4(f15.5))') x,y,z,conc_react1(i,j,k)
                  write(13,'(4(f15.5))') x,y,z,conc_react2(i,j,k)
                  write(25,'(4(f15.5))') x,y,z,conc_product(i,j,k)
                  write(15,'(4(f15.5))') x,y,z,sigc_react1(i,j,k)
                  write(16,'(4(f15.5))') x,y,z,sigc_react2(i,j,k)
                  write(26,'(4(f15.5))') x,y,z,sigc_product(i,j,k)
                  write(17,'(4(f15.5))') x,y,z,c_skew_react1(i,j,k)
                  write(18,'(4(f15.5))') x,y,z,c_skew_react2(i,j,k)
                  write(27,'(4(f15.5))') x,y,z,c_skew_product(i,j,k)
! Post-processing limiters for c_kurt_react1 and c_kurt_react2 
                  if ((c_kurt_react1(i,j,k).gt.-99999999.d0).and.
     &               (c_kurt_react1(i,j,k).lt.(99999999.d0))) then 
                     write(19,'(4(f15.5))') x,y,z,c_kurt_react1(i,j,k)
                     else
                        write(19,'(3(f15.5),a)') x,y,z,"     -999.00000"
                  endif
                  if ((c_kurt_react2(i,j,k).gt.-99999999.d0).and.
     &               (c_kurt_react2(i,j,k).lt.(99999999.d0))) then 
                     write(20,'(4(f15.5))') x,y,z,c_kurt_react2(i,j,k)
                     else
                        write(20,'(3(f15.5),a)') x,y,z,"     -999.00000"
                  endif
                  write(28,'(4(f15.5))') x,y,z,c_kurt_product(i,j,k)
                  write(21,'(4(f15.5))') x,y,z,Iseg(i,j,k)
! AA!!! test: start
                  write(31,'(4(f15.5))') x,y,z,u(i,j,k)
                  write(32,'(4(f15.5))') x,y,z,v(i,j,k)
                  write(33,'(4(f15.5))') x,y,z,w(i,j,k)
                  write(34,'(4(f15.5))') x,y,z,varu(i,j,k)
                  write(35,'(4(f15.5))') x,y,z,varv(i,j,k)
                  write(36,'(4(f15.5))') x,y,z,varw(i,j,k)
                  write(37,'(4(f15.5))') x,y,z,eps(i,j,k)
! AA!!! test: end
               endif
            enddo
         enddo
      enddo
      do i=1,ndx
         x = dx * (i - 1 + 0.5d0)
         write(14,'(2(f15.5))') x,sigz(i)
      enddo
      close(11)
      close(14)
      close(22)
      close(23)
      close(24)
      if (reactive==1) then
         close(12)
         close(13)
         close(15)
         close(16)
         close(17)
         close(18)
         close(19)
         close(20)
         close(21)
         close(25)
         close(26)
         close(27)
         close(28)
      endif
! AA!!! test: start
      close(31)
      close(32)
      close(33)
      close(34)
      close(35)
      close(36)
      close(37)
! AA!!! test: end      
      return
      end 

