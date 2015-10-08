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
! Subroutine name: LMM_sections
! Subroutine description: post-processing the fields of Reynolds' scalar
! stastistics to extract 2D fields (sections). 
!-----------------------------------------------------------------------
      subroutine LMM_sections(conc,sigc,c_skew,c_kurt,conc_react1,
     &conc_react2,dz,dx,dy,reactive,tm,ndx,ndy,ndz)
!-----------------------------------------------------------------------
! Declarations
!-----------------------------------------------------------------------
      integer :: reactive,ndx,ndy,ndz
      real :: dx,dy
      double precision :: dz
      double precision :: conc(ndx,ndy,ndz),sigc(ndx,ndy,ndz)
      double precision :: c_skew(ndx,ndy,ndz),c_kurt(ndx,ndy,ndz)
      double precision :: conc_react1(ndx,ndy,ndz)
      double precision :: conc_react2(ndx,ndy,ndz),tm(ndx,ndy,ndz)
! Local variables      
      integer :: iy1,iy2,iy3,ix1,ix2,ix3 
      real :: x,y,z
!-----------------------------------------------------------------------
! Initializations
!-----------------------------------------------------------------------
! Indices for plotting vertical sections      
      iy1 = 1
      iy2 = 1
      iy3 = 1
      ix1 = 1
      ix2 = 1
      ix3 = 1
!-----------------------------------------------------------------------
! Statements
!-----------------------------------------------------------------------
      open(25,file='Cm_iy1')
      open(26,file='Cm_iy2')
      open(27,file='Cm_iy3')
      open(28,file='sigmaC_iy1')
      open(29,file='sigmaC_iy2')
      open(30,file='sigmaC_iy3')
      open(38,file='tm_iy3')
      open(44,file='c_skew_iy1')
      open(45,file='c_skew_iy2')
      open(46,file='c_skew_iy3')
      open(47,file='c_kurt_iy1')
      open(48,file='c_kurt_iy2')
      open(49,file='c_kurt_iy3')
      do j=1,ndz
         do i=1,ndx
            x = (i - 1) * dx + dx / 2.
            z = (j - 1) * dz + dz / 2.
300   format(1x,3(f12.4,1x))
            write(25,fmt=300) x,z,conc(i,iy1,j)
            write(26,fmt=300) x,z,conc(i,iy2,j)
            write(27,fmt=300) x,z,conc(i,iy3,j)
            write(28,fmt=300) x,z,sigc(i,iy1,j)
            write(29,fmt=300) x,z,sigc(i,iy2,j)
            write(30,fmt=300) x,z,sigc(i,iy3,j)
            write(38,fmt=300) x,z,tm(i,iy3,j)
            write(44,fmt=300) x,z,c_skew(i,iy1,j)
            write(45,fmt=300) x,z,c_skew(i,iy2,j)
            write(46,fmt=300) x,z,c_skew(i,iy3,j)
            write(47,fmt=300) x,z,c_kurt(i,iy1,j)
            write(48,fmt=300) x,z,c_kurt(i,iy2,j)
            write(49,fmt=300) x,z,c_kurt(i,iy3,j)
	 enddo
      enddo
      close(25)
      close(26)
      close(27)
      close(28)
      close(29)
      close(30)
      close(38)
      close(44)
      close(45)
      close(46)
      close(47)
      close(48)
      close(49)
      if (reactive==1) then
         open(25,file='Cm_r1_iy1')
         open(26,file='Cm_r1_iy2')
         open(27,file='Cm_r1_iy3')
         do j=1,ndz
            do i=1,ndx
               x = (i - 1) * dx + dx / 2.
               z = (j - 1) * dz + dz / 2.
               write(25,fmt=300) x,z,conc_react1(i,iy1,j)
               write(26,fmt=300) x,z,conc_react1(i,iy2,j)
               write(27,fmt=300) x,z,conc_react1(i,iy3,j)
            enddo
         enddo
         close(25)
         close(26)
         close(27)
         open(25,file='Cm_r2_iy1')
         open(26,file='Cm_r2_iy2')
         open(27,file='Cm_r2_iy3')
         do j=1,ndz
            do i=1,ndx
               x = (i - 1) * dx + dx / 2.
               z = (j - 1) * dz + dz / 2.
               write(25,fmt=300) x,z,conc_react2(i,iy1,j)
               write(26,fmt=300) x,z,conc_react2(i,iy2,j)
               write(27,fmt=300) x,z,conc_react2(i,iy3,j)
            enddo
         enddo
         close(25)
         close(26)
         close(27)
         endif
         open(31,file='Cm_ix1')
         open(32,file='Cm_ix2')
         open(33,file='Cm_ix3')
         open(34,file='sigmaC_ix1')
         open(35,file='sigmaC_ix2')
         open(36,file='sigmaC_ix3')
         open(39,file='tm_ix1')
         open(50,file='c_skew_ix1')
         open(51,file='c_skew_ix2')
         open(52,file='c_skew_ix3')
         open(53,file='c_kurt_ix1')
         open(54,file='c_kurt_ix2')
         open(55,file='c_kurt_ix3')
         do i=1,ndy
            do j=1,ndz
               y = (i - 1) * dy + dy / 2.
               z = (j - 1) * dz + dz / 2.
               write(31,fmt=300) y,z,conc(ix1,i,j)
               write(32,fmt=300) y,z,conc(ix2,i,j)
               write(33,fmt=300) y,z,conc(ix3,i,j)
               write(34,fmt=300) y,z,sigc(ix1,i,j)
               write(35,fmt=300) y,z,sigc(ix2,i,j)
               write(36,fmt=300) y,z,sigc(ix3,i,j)
               write(39,fmt=300) y,z,tm(ix1,i,j)
               write(50,fmt=300) y,z,c_skew(ix1,i,j)
               write(51,fmt=300) y,z,c_skew(ix2,i,j)
               write(52,fmt=300) y,z,c_skew(ix3,i,j)
               write(53,fmt=300) y,z,c_kurt(ix1,i,j)
               write(54,fmt=300) y,z,c_kurt(ix2,i,j)
               write(55,fmt=300) y,z,c_kurt(ix3,i,j)
            enddo
         enddo
         close(31)
         close(32)
         close(33)
         close(34)
         close(35)
         close(36)
         close(39)
         close(50)
         close(51)
         close(52)
         close(53)
         close(54)
         close(55)
         if (reactive==1) then
            open(31,file='Cm_r1_ix1')
            open(32,file='Cm_r1_ix2')
            open(33,file='Cm_r1_ix3')
            do j=1,ndz
               do i=1,ndy
                  y = (i - 1) * dy + dy / 2.
                  z = (j - 1) * dz + dz / 2.
                  write (31,fmt=300) y,z,conc_react1(ix1,i,j)
                  write (32,fmt=300) y,z,conc_react1(ix2,i,j)
                  write (33,fmt=300) y,z,conc_react1(ix3,i,j)
            enddo
         enddo
         close(31)
         close(32)
         close(33)
         open(31,file='Cm_r2_ix1')
         open(32,file='Cm_r2_ix2')
         open(33,file='Cm_r2_ix3')
         do j=1,ndz
            do i=1,ndy
               y = (i - 1) * dy + dy / 2.
               z = (j - 1) * dz + dz / 2.
               write(31,fmt=300) y,z,conc_react2(ix1,i,j)
               write(32,fmt=300) y,z,conc_react2(ix2,i,j)
               write(33,fmt=300) y,z,conc_react2(ix3,i,j)
            enddo
         enddo
         close(31)
         close(32)
         close(33)
      endif
      return
      end

