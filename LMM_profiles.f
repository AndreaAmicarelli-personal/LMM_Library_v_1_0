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
! Subroutine name: LMM_profiles
! Subroutine description. Post-processing the fields of Reynolds' scalar
! stastistics to extract horizontal and vertical profiles, which are 
! parallel to the reference axis. 
!-----------------------------------------------------------------------
      subroutine LMM_profiles(conc,sigc,c_skew,c_kurt,conc_react1,
     &conc_react2,dz,dx,dy,case2D,reactive,ndx,ndy,ndz,n_prof,
     &n_prof_vert)
!-----------------------------------------------------------------------
! Declarations
!-----------------------------------------------------------------------
      integer :: case2D,reactive,ndx,ndy,ndz,n_prof,n_prof_vert
      real :: dx,dy
      double precision :: dz
      double precision :: conc(ndx,ndy,ndz),sigc(ndx,ndy,ndz)
      double precision :: c_skew(ndx,ndy,ndz),c_kurt(ndx,ndy,ndz)
      double precision :: conc_react1(ndx,ndy,ndz)
      double precision :: conc_react2(ndx,ndy,ndz)
! Local variables
      integer :: IOstatus
      integer :: x_prof_vert_int(n_prof_vert),y_prof_int(n_prof)
      real :: y_prof(n_prof),x_prof_vert(n_prof_vert)
      real :: Cm_prof(ndx,n_prof+1),sigC_prof(ndx,n_prof+1)
      real :: iC_prof(ndx,n_prof+1),SkC_prof(ndx,n_prof+1)
      real :: KuC_prof(ndx,n_prof+1),Cm_react1_prof(ndx,n_prof+1)
      real :: Cm_react2_prof(ndx,n_prof+1)
      real :: Cm_prof_vert(ndz,n_prof_vert+1)
      real :: sigC_prof_vert(ndz,n_prof_vert+1)
      real :: iC_prof_vert(ndz,n_prof_vert+1)
      real :: SkC_prof_vert(ndz,n_prof_vert+1)
      real :: KuC_prof_vert(ndz,n_prof_vert+1)
      real :: Cm_react1_prof_vert(ndz,n_prof_vert+1)
      real :: Cm_react2_prof_vert(ndz,n_prof_vert+1)
!-----------------------------------------------------------------------
! Initializations
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Statements
!-----------------------------------------------------------------------
! Horizontal profiles: start
      if (case2D.ne.1) then 
         open(56,file='y_profiles.txt')
1120     format(/)
         read(56,1120)
         read(56,*) (y_prof(i),i=1,n_prof)
         do i=1,n_prof
            y_prof_int(i) = int(y_prof(i) / dy) + 1
         enddo
         do i=1,ndx
            Cm_prof(i,1) = i * dx - dx / 2.
            sigC_prof(i,1) = i * dx - dx / 2.
            iC_prof(i,1) = i * dx - dx / 2.
            SkC_prof(i,1) = i * dx - dx / 2.
            KuC_prof(i,1) = i * dx - dx / 2.
            if (reactive==1) then
               Cm_react1_prof(i,1) = i * dx - dx / 2.
               Cm_react2_prof(i,1) = i * dx - dx / 2.
            endif
         enddo
         do j=2,(n_prof+1)
            do i=1,ndx
               Cm_prof(i,j) = conc(i,y_prof_int(j-1),3)
               sigC_prof(i,j) = sigc(i,y_prof_int(j-1),3)
               if ((conc(i,y_prof_int(j-1),3).eq.-999.d0).or.
     &            (sigc(i,y_prof_int(j-1),3).eq.-999.d0).or.
     &            (conc(i,y_prof_int(j-1),3).lt.(0.000001d0)).or.
     &            ((conc(i,y_prof_int(j-1),3).eq.(0.d0)).and.
     &            (sigc(i,y_prof_int(j-1),3).eq.(0.d0)))) then
                  iC_prof(i,j) = -999.d0
                  else 
                     iC_prof(i,j) = sigc(i,y_prof_int(j-1),3) / 
     &                              conc(i,y_prof_int(j-1),3)
               endif
               if ((c_skew(i,y_prof_int(j-1),3).eq.-999999.d0).or.
     &            (c_skew(i,y_prof_int(j-1),3).eq.(100000.d0)).or.
     &            (c_skew(i,y_prof_int(j-1),3).eq.-100000d0)) then
                  SkC_prof(i,j) = -999.d0
                  else 
                     SkC_prof(i,j) = c_skew(i,y_prof_int(j-1),3)
               endif
               if ((c_kurt(i,y_prof_int(j-1),3).eq.-999.d0).or.
     &            (c_kurt(i,y_prof_int(j-1),3).eq.(100000.d0)).or.
     &            (c_kurt(i,y_prof_int(j-1),3).eq.(-1.d0))) then
                  KuC_prof(i,j) = -999.d0
                  else 
                     KuC_prof(i,j) = c_kurt(i,y_prof_int(j-1),3)
               endif
               if (reactive==1) then
                  Cm_react1_prof(i,j) = conc_react1(i,y_prof_int(j-1),3)
                  Cm_react2_prof(i,j) = conc_react2(i,y_prof_int(j-1),3)
               endif
            enddo
         enddo
         close(56)
         open(57,file='Cm_profile.txt')
1122     format(12(F12.6))
1123     format(11(a12))
         write(57,1123) 'x(m)','Cm_prof1','Cm_prof2','Cm_prof3',
     &      'Cm_prof4','Cm_prof5','Cm_prof6','Cm_prof7','Cm_prof8',
     &      'Cm_prof9','Cm_prof10','Cm_prof11'
         do i=1,ndx
            write (57,*) (Cm_prof(i,j),j=1,(n_prof+1))
         enddo
         close (57)
   ! sigma_C profiles
         open(58,file='sigmaC_profile.txt')
         write (58,1123) 'x(m)','sigmaC_prof1','sigmaC_prof2',
     &      'sigmaC_prof3','sigmaC_prof4', 'sigmaC_prof5',
     &      'sigmaC_prof6','sigmaC_prof7','sigmaC_prof8',
     &      'sigmaC_prof9','sigmaC_prof10','sigmaC_prof11'
         do i=1,ndx
            write (58,*) (sigC_prof(i,j),j=1,(n_prof+1))
         enddo
         close (58)
   ! i_C profiles
         open(59,file='iC_profile.txt')
         write (59,1123) 'x(m)','i_prof1','i_prof2','i_prof3','i_prof4',
     &      'i_prof5','i_prof6','i_prof7','i_prof8','i_prof9','i_prof10'
     &      ,'i_prof11'
         do i=1,ndx
            write (59,*) (iC_prof(i,j),j=1,(n_prof+1))
         enddo
         close (59)
   ! Sk_C profiles
         open(60,file='SkC_profile.txt')
         write (60,1123) 'x(m)','SkC_prof1','SkC_prof2','SkC_prof3',
     &      'SkC_prof4', 'SkC_prof5','SkC_prof6','SkC_prof7','SkC_prof8'
     &      ,'SkC_prof9','SkC_prof10','SkC_prof11'
         do i=1,ndx
            write (60,*) (SkC_prof(i,j),j=1,(n_prof+1))
         enddo
         close (60)
   ! Ku_C profiles
         open(61,file='KuC_profile.txt')
         write (61,1123) 'x(m)','KuC_prof1','KuC_prof2','KuC_prof3',
     &      'KuC_prof4', 'KuC_prof5','KuC_prof6','KuC_prof7','KuC_prof8'
     &      ,'KuC_prof9','KuC_prof10','KuC_prof11'
         do i=1,ndx
            write (61,*) (KuC_prof(i,j),j=1,(n_prof+1))
         enddo
         close (61)
         if (reactive==1) then
   ! Profiles of reactant 1
            open(57,file='Cm_r1_profile.txt')
            write (57,1123) 'x(m)','Cm_r1_prof1','Cm_r1_prof2',
     &         'Cm_r1_prof3','Cm_r1_prof4','Cm_r1_prof5','Cm_r1_prof6',
     &         'Cm_r1_prof7','Cm_r1_prof8','Cm_r1_prof9','Cm_r1_prof10',
     &         'Cm_r1_prof11'
            do i=1,ndx
               write (57,*) (Cm_react1_prof(i,j),j=1,(n_prof+1))
            enddo
            close (57)
   ! Profiles of reactant 2
            open(57,file='Cm_r2_profile.txt')
            write (57,1123) 'x(m)','Cm_r2_prof1','Cm_r2_prof2',
     &         'Cm_r2_prof3','Cm_r2_prof4','Cm_r2_prof5','Cm_r2_prof6',
     &         'Cm_r2_prof7','Cm_r2_prof8','Cm_r2_prof9','Cm_r2_prof10',
     &         'Cm_r2_prof11'
            do i=1,ndx
               write (57,*) (Cm_react2_prof(i,j),j=1,(n_prof+1))
            enddo
            close (57)
         endif
      endif
! Horizontal profiles: end 
! Vertical profiles: start 
      open(62,file='x_vertical_profiles.txt')
1124  format (11(F12.6))
      read (62,1120)
      read (62,*,IOSTAT=IOstatus) (x_prof_vert(i),i=1,n_prof_vert)
      do i=1,n_prof_vert
         x_prof_vert_int(i) = int(x_prof_vert(i) / dx) + 1
      enddo
      do i=1,ndz
         Cm_prof_vert(i,1) = i * dz - dz / 2.
         sigC_prof_vert(i,1) = i * dz - dz / 2.
         iC_prof_vert(i,1) = i * dz - dz / 2.
         SkC_prof_vert(i,1) = i * dz - dz / 2.
         KuC_prof_vert(i,1) = i * dz - dz / 2.
         if (reactive==1) then
            Cm_react1_prof_vert(i,1) = i * dz - dz / 2.
            Cm_react2_prof_vert(i,1) = i * dz - dz / 2.
         endif
      enddo
      do j=2,(n_prof_vert+1)
         do i=1,ndz
            Cm_prof_vert(i,j) = conc(x_prof_vert_int(j-1),1,i)
            sigC_prof_vert(i,j) = sigc(x_prof_vert_int(j-1),1,i)
            if ((conc(x_prof_vert_int(j-1),1,i).eq.-999.d0).or.
     &         (sigc(x_prof_vert_int(j-1),1,i).eq.-999.d0).or.
     &         (conc(x_prof_vert_int(j-1),1,i).lt.(0.000001d0)).or.
     &         ((conc(x_prof_vert_int(j-1),1,i).eq.(0.d0)).and.
     &         (sigc(x_prof_vert_int(j-1),1,i).eq.(0.d0)))) then
               iC_prof(i,j)=-999.d0
               else 
                  iC_prof(i,j) = sigc(x_prof_vert_int(j-1),1,i) /
     &                           conc(x_prof_vert_int(j-1),1,i)
            endif
            if ((c_skew(x_prof_vert_int(j-1),1,i).eq.-999999.d0).or.
     &         (c_skew(x_prof_vert_int(j-1),1,i).eq.(100000.d0)).or.
     &         (c_skew(x_prof_vert_int(j-1),1,i).eq.(-100000.d0))) then
               SkC_prof(i,j)=-999.d0
               else 
                  SkC_prof(i,j) = c_skew(x_prof_vert_int(j-1),1,i)
            endif
            if ((c_kurt(x_prof_vert_int(j-1),1,i).eq.-999.d0).or.
     &         (c_kurt(x_prof_vert_int(j-1),1,i).eq.(100000.d0)).or.
     &         (c_kurt(x_prof_vert_int(j-1),1,i).eq.(-1.d0))) then
               KuC_prof(i,j)=-999.d0
               else 
                  KuC_prof(i,j)=c_kurt(x_prof_vert_int(j-1),1,i)
            endif
            if (reactive==1) then
               Cm_react1_prof_vert(i,j) =
     &conc_react1(x_prof_vert_int(j-1),1,i)
               Cm_react2_prof_vert(i,j)=
     &conc_react2(x_prof_vert_int(j-1),1,i)
            endif
         enddo
      enddo
      close(62)
1125  format(12(F12.6))
1126  format(11(a12))
! C_m profiles
      open(63,file='Cm_vertical_profile.txt')
      write(63,1126) 'z(m)','Cm_prov1','Cm_prov2','Cm_prov3',
     &   'Cm_prov4','Cm_prov5','Cm_prov6','Cm_prov7','Cm_prov8',
     &   'Cm_prov9','Cm_prov10','Cm_prov11'
      do i=1,ndz
         write(63,*) (Cm_prof_vert(i,j),j=1,(n_prof_vert+1))
      enddo
      close (63)
! sigma_C profiles
      open(64,file='sigmaC_vertical_profile.txt')
      write (64,1126) 'z(m)','sigC_prov1','sigC_prov2','sigC_prov3',
     &   'sigC_prov4', 'sigC_prov5','sigC_prov6','sigC_prov7',
     &   'sigC_prov8','sigC_prov9','sigC_prov10','sigC_prov11'
      do i=1,ndz
         write(64,*) (sigC_prof_vert(i,j),j=1,(n_prof_vert+1))
      enddo
      close(64)
! i_C profiles
      open(65,file='iC_vertical_profile.txt')
      write (65,1126) 'z(m)','iC_prov1','iC_prov2','iC_prov3',
     &'iC_prov4', 'iC_prov5','iC_prov6','iC_prov7','iC_prov8',
     &'iC_prov9','iC_prov10','iC_prov11'
      do i=1,ndz
         write(65,*) (iC_prof_vert(i,j),j=1,(n_prof_vert+1))
      enddo
      close(65)
! Sk_C profiles
      open(66,file='SkC_vertical_profile.txt')
      write(66,1126) 'z(m)','SkC_prov1','SkC_prov2','SkC_prov3',
     &   'SkC_prov4', 'SkC_prov5','SkC_prov6','SkC_prov7','SkC_prov8',
     &   'SkC_prov9','SkC_prov10','SkC_prov11'
      do i=1,ndz
         write(66,*) (SkC_prof_vert(i,j),j=1,(n_prof_vert+1))
      enddo
      close(66)
! Ku_C profiles
      open(67,file='KuC_vertical_profile.txt')
      write (67,1126) 'z(m)','KuC_prov1','KuC_prov2','KuC_prov3',
     &   'KuC_prov4', 'KuC_prov5','KuC_prov6','KuC_prov7','KuC_prov8',
     &   'KuC_prov9','KuC_prov10','KuC_prov11'
      do i=1,ndz
         write(67,*) (KuC_prof_vert(i,j),j=1,(n_prof_vert+1))
      enddo
      close(67)
      if (reactive==1) then
! C_m_react1 vertical profiles
         open(63,file='Cm_r1_vertical_profile.txt')
         write(63,1126) 'z(m)','Cm_r1_prov1','Cm_r1_prov2',
     &      'Cm_r1_prov3','Cm_r1_prov4','Cm_r1_prov5','Cm_r1_prov6',
     &      'Cm_r1_prov7','Cm_r1_prov8','Cm_r1_prov9','Cm_r1_prov10',
     &      'Cm_r1_prov11'
         do i=1,ndz
            write(63,*) (Cm_react1_prof_vert(i,j),j=1,(n_prof_vert+1))
         enddo
         close(63)
! C_m_react2 vertical profiles
         open(63,file='Cm_r2_vertical_profile.txt')
         write(63,1126) 'z(m)','Cm_r2_prov1','Cm_r2_prov2',
     &      'Cm_r2_prov3','Cm_r2_prov4','Cm_r2_prov5','Cm_r2_prov6',
     &      'Cm_r2_prov7','Cm_r2_prov8','Cm_r2_prov9','Cm_r2_prov10',
     &      'Cm_r2_prov11'
         do i=1,ndz
            write(63,*) (Cm_react2_prof_vert(i,j),j=1,(n_prof_vert+1))
         enddo
         close(63)
      endif
! Vertical profiles: end
      return
      end
      
