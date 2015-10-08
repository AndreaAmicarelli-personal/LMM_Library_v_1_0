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
! Subroutine name: concentration_statistics
! Subroutine description: to compute the fields of Reynolds' mean,   
! standard deviation, skewness and kurtosis of reactive scalars (2 
! reactants, 1 product species), provided their first four moments (LMM 
! non-conservative phase). This library computes the fields of Reynolds'
! skewness and kurtosis of a passive scalar, provided the sum of the  
! differences (with respect to the mean) to the 3rd/4th power and the  
! number of samples (LMM non-conservative phase).
!----------------------------------------------------------------------
      subroutine concentration_statistics(ndx,ndy,ndz,isig,reactive,
     &conc_react1,conc_react2,C0_react2,conc_product,sigc_react1,
     &sigc_react2,sigc_product,c_skew_react1,c_kurt_react1,
     &c_skew_react2,c_kurt_react2,c_skew_product,c_kurt_product,Iseg,
     &c_skew,c_kurt,sigc)
!----------------------------------------------------------------------
! Declarations
!----------------------------------------------------------------------
      integer :: ndx,ndy,ndz,reactive
      double precision :: C0_react2
      integer :: isig(ndx,ndy,ndz)
      double precision :: conc_react1(ndx,ndy,ndz)
      double precision :: conc_react2(ndx,ndy,ndz)
      double precision :: conc_product(ndx,ndy,ndz)
      double precision :: sigc_react1(ndx,ndy,ndz)
      double precision :: sigc_react2(ndx,ndy,ndz)
      double precision :: sigc_product(ndx,ndy,ndz)
      double precision :: c_skew_react1(ndx,ndy,ndz)
      double precision :: c_kurt_react1(ndx,ndy,ndz)
      double precision :: c_skew_react2(ndx,ndy,ndz)
      double precision :: c_kurt_react2(ndx,ndy,ndz)
      double precision :: c_skew_product(ndx,ndy,ndz)
      double precision :: c_kurt_product(ndx,ndy,ndz)
      double precision :: Iseg(ndx,ndy,ndz),c_skew(ndx,ndy,ndz)
      double precision :: c_kurt(ndx,ndy,ndz),sigc(ndx,ndy,ndz)
      integer :: i,j,k
!-----------------------------------------------------------------------
! Initializations
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Statements
!-----------------------------------------------------------------------
      do k=1,ndz
         do j=1,ndy
            do i=1,ndx
               if (isig(i,j,k).eq.0)then
                  if (reactive==1) then
                     conc_react1(i,j,k) = 0.0d0 
                     conc_react2(i,j,k) = C0_react2
                     conc_product(i,j,k) = 0.0d0
                     sigc_react1(i,j,k) = 0.0d0
                     sigc_react2(i,j,k) = 0.0d0
                     sigc_product(i,j,k) = 0.0d0
                  endif
                  else
                     if (reactive==1) then
                        conc_react1(i,j,k) = conc_react1(i,j,k) /
     &                                       isig(i,j,k)
                        conc_react2(i,j,k) = conc_react2(i,j,k) / 
     &                                       isig(i,j,k)
                        conc_product(i,j,k) = conc_product(i,j,k) / 
     &                                        isig(i,j,k)
                        sigc_react1(i,j,k) = (sigc_react1(i,j,k) / 
     &                                       isig(i,j,k) - 
     &                                       conc_react1(i,j,k) ** 2)
                        sigc_react2(i,j,k) = (sigc_react2(i,j,k) / 
     &                                       isig(i,j,k) - 
     &                                       conc_react2(i,j,k) ** 2)
                        sigc_product(i,j,k) = (sigc_product(i,j,k) / 
     &                                        isig(i,j,k) - 
     &                                        conc_product(i,j,k) ** 2)
                        if(sigc_react1(i,j,k).ne.0.) sigc_react1(i,j,k) 
     &                     = sqrt(sigc_react1(i,j,k))
                        if(sigc_react2(i,j,k).ne.0.) sigc_react2(i,j,k) 
     &                     = sqrt(sigc_react2(i,j,k))
                        if(sigc_product(i,j,k).ne.0.) 
     &                     sigc_product(i,j,k) = 
     &                     sqrt(sigc_product(i,j,k))
                        if (sigc_react1(i,j,k).ne.0.d0) then
                           c_skew_react1(i,j,k) = (c_skew_react1(i,j,k)
     &                                            / isig(i,j,k) - 
     &                                            conc_react1(i,j,k) **
     &                                            3) / 
     &                                            sigc_react1(i,j,k) **
     &                                            3 - 3.d0 * 
     &                                            conc_react1(i,j,k) / 
     &                                            sigc_react1(i,j,k)
                           c_kurt_react1(i,j,k) = (c_kurt_react1(i,j,k)
     &                                            / isig(i,j,k) - 
     &                                            conc_react1(i,j,k) **
     &                                            4) / 
     &                                            sigc_react1(i,j,k) **
     &                                            4 - 6.d0 * 
     &                                            (conc_react1(i,j,k) / 
     &                                            sigc_react1(i,j,k)) **
     &                                            2 - 4.d0 * 
     &                                            c_skew_react1(i,j,k) *
     &                                            conc_react1(i,j,k) / 
     &                                            sigc_react1(i,j,k)
                        endif
                        if (sigc_react2(i,j,k).ne.0.d0) then
                           c_skew_react2(i,j,k) = (c_skew_react2(i,j,k)
     &                                            / isig(i,j,k) - 
     &                                            conc_react2(i,j,k) **
     &                                            3) / 
     &                                            sigc_react2(i,j,k) **
     &                                            3 - 3.d0 * 
     &                                            conc_react2(i,j,k) / 
     &                                            sigc_react2(i,j,k)
                           c_kurt_react2(i,j,k) = (c_kurt_react2(i,j,k)
     &                                            / isig(i,j,k) -
     &                                            conc_react2(i,j,k) **
     &                                            4) / 
     &                                            sigc_react2(i,j,k) **
     &                                            4 - 6.d0 * 
     &                                            (conc_react2(i,j,k) /
     &                                            sigc_react2(i,j,k)) **
     &                                            2 - 4.d0 * 
     &                                            c_skew_react2(i,j,k) *
     &                                            conc_react2(i,j,k) /
     &                                            sigc_react2(i,j,k)
                        endif
                        if (sigc_product(i,j,k).ne.0.d0) then
                           c_skew_product(i,j,k) = 
     &(c_skew_product(i,j,k) / isig(i,j,k) - conc_product(i,j,k) ** 3) /
     &                                             sigc_product(i,j,k)
     &                                             ** 3 - 3.d0 * 
     &                                             conc_product(i,j,k) /
     &                                             sigc_product(i,j,k)
                           c_kurt_product(i,j,k) = 
     &(c_kurt_product(i,j,k) / isig(i,j,k) - conc_product(i,j,k) ** 4) /
     &                                             sigc_product(i,j,k)
     &                                             ** 4 - 6.d0 * 
     &                                             (conc_product(i,j,k)
     &                                             / sigc_product(i,j,k)
     &                                             )**2 - 4.d0 * 
     &                                             c_skew_product(i,j,k)
     &                                             * conc_product(i,j,k)
     &                                             / sigc_product(i,j,k)
                        endif
                        if ((conc_react1(i,j,k).ne.0.d0).and.
     &                     (conc_react2(i,j,k).ne.0.d0)) then
                           Iseg(i,j,k) = Iseg(i,j,k) / (isig(i,j,k) * 
     &                                   conc_react1(i,j,k) * 
     &                                   conc_react2(i,j,k)) - 1.d0
                        endif
                     endif
! Mean cubic difference (for skewness) and mean quartic difference (for
! kurtosis) where non-conservative particles had been detected
                     c_skew(i,j,k) = c_skew(i,j,k) / isig(i,j,k) 
                     c_kurt(i,j,k) = c_kurt(i,j,k) / isig(i,j,k)
               endif
! Computation of skewness and kurtosis (passive scalar) where 
! non-conservative particles had been detected: normalization of the 3rd
! and 4th moment of the scalar pdf. To invalidate the field values where  
! non-conservative particles are absent. 
               if (isig(i,j,k).ne.0.) then
                  if (sigc(i,j,k).ne.0.) then
                     c_skew(i,j,k) = c_skew(i,j,k) / (sigc(i,j,k) ** 3)
                     c_kurt(i,j,k) = c_kurt(i,j,k) / (sigc(i,j,k) ** 4)
	             else
                        c_skew(i,j,k) = 100000.d0
                        c_kurt(i,j,k) = 100000.d0
	          endif
	          else
	             c_skew(i,j,k) = -999999.d0
	             c_kurt(i,j,k) = -1.d0
	       endif
            enddo
         enddo
      enddo
      return
      end
  

