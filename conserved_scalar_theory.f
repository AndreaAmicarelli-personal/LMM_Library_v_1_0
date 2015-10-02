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
! Subroutine name: conserved_scalar_theory
! Subroutine description. To compute the limits of the conserved scalar 
! theory, provided the particle instantaneous value of the mixture 
! fraction for a 2nd-order kinetics reaction. To update the 
! contributions to the fields of the associated concentration statistics
! (mean, standard deviation, skewness, kurtosis, segregation 
! coefficient).
!-----------------------------------------------------------------------
      subroutine conserved_scalar_theory(time,ipcont,cpr,ipu,jpv,kpw,
     &concc,cnew,conc,sigc,isig,c_skew,c_kurt,conc_react1,
     &conc_react2,sigc_react1,c_skew_react1,c_kurt_react1,
     &sigc_react2,c_skew_react2,c_kurt_react2,conc_product,
     &sigc_product,c_skew_product,c_kurt_product,Iseg,sigz,isigz,cmax,
     &C_A1,C_B2,reaction_rate,U_scale,dt,dx,t_re_flag,react_lim,
     &sum_C_A0_C_B0,ndx,ndy,ndz,nbu,nbv,nbw,inew,jnew,knew,xpr,
     &n_warn_Fm_max,n_warn_Fm_min,n_warn_C_A_max,n_warn_C_A_min,
     &n_warn_C_C_min,n_warn_C_B_max,n_warn_C_B_min)
!-----------------------------------------------------------------------
! Declarations
!-----------------------------------------------------------------------
      integer :: t_re_flag,react_lim,ndx,ndy,ndz,nbu,nbv,nbw,inew,jnew
      integer :: n_warn_Fm_max,n_warn_Fm_min,n_warn_C_A_max
      integer :: n_warn_C_B_max,n_warn_C_B_min,n_warn_C_C_min,ipcont,ipu
      integer :: jpv,kpw,knew,n_warn_C_A_min 
      real :: cmax,dt,dx,xpr,time,cpr
      double precision :: reaction_rate,U_scale,sum_C_A0_C_B0,C_A1,C_B2
      double precision :: cnew
      integer :: isigz(ndx)
      double precision :: sigz(ndx)
      integer :: isig(ndx,ndy,ndz)
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
      double precision :: concc(ndx,ndy,ndz,nbu,nbv,nbw)
! Local variables
      double precision :: F_m,F_ms,C_A,C_B,C_C,C_A0,C_B0
      double precision :: time_react_eff,x_cell
!-----------------------------------------------------------------------
! Initializations
!-----------------------------------------------------------------------
      F_m = cnew
!-----------------------------------------------------------------------
! Statements
!-----------------------------------------------------------------------
      if (F_m.gt.1.d0) then
         if ((F_m.gt.1.05d0).and.(n_warn_Fm_max.le.1000)) then 
            write(*,*) 
     &         "Warning:F_m>1.05; F_m, time, ix_cel, iy_cel, iz_cel,",
     &         " ID_part, Fm_m_c, Fm_old: "
            write(*,*) F_m,time,inew,jnew,knew,ipcont,
     &         concc(inew,jnew,knew,ipu,jpv,kpw),cpr
            n_warn_Fm_max = n_warn_Fm_max + 1
            if (n_warn_Fm_max.eq.1000) then
               write(*,*) 
     &            "Reached the maximum number (1000) of warnings ",
     &            "on Fm_max (non-conservative phase)"
            endif
         endif
         F_m = 1.d0
      endif
      if (F_m.lt.0.d0) then
         if ((F_m.lt.(-(0.05d0))).and.(n_warn_Fm_min.le.1000)) then 
            write(*,*) 
     &         "Warning:F_m<0.05; F_m, time, ix_cel, iy_cel, iz_cel,",
     &         " ID_part, Fm_m_c, Fm_old: "
            write(*,*) F_m,time,inew,jnew,knew,ipcont,
     &         concc(inew,jnew,knew,ipu,jpv,kpw),cpr
            n_warn_Fm_min = n_warn_Fm_min + 1
            if (n_warn_Fm_min.eq.1000) then
               write(*,*) 
     &            "Reached the maximum number (1000) of warnings ",
     &            "on Fm_min (non-conservative phase)"
            endif
         endif
         F_m = 0.d0
      endif
      F_ms = C_B2 / (C_A1 + C_B2)
! Choosing the reaction limit:
! 0: Main Reaction-Dominated Limits (RDL, NHRDL, RDL_plume)
! 1: Frozen Limit (FL)
! 2: Equilibrium Limit (EL)
! 3: Reaction Dominated Limit based on the reactant means 
! (NFRDL, NI-NHRDL) 
      if (react_lim.eq.3) F_m = conc(inew,jnew,knew)
      if (react_lim.eq.2) then
         if (F_m.gt.F_ms) then          
            C_A = C_A1 * (F_m - F_ms) / (1.d0 - F_ms)
            C_B = 0.d0
            else
               C_A = 0.d0
               C_B = C_B2 * (F_ms - F_m) / F_ms
         endif
         elseif (react_lim.eq.1) then 
            C_A = F_m * C_A1
            C_B = (1.d0 - F_m) * C_B2
            elseif ((react_lim.eq.0).or.(react_lim.eq.3)) then
! Modified formulation for the reaction/contact time (EL and FL are not 
! affected by this part):
! t_re_flag=0: contact time=fly time (since the numerical release -from 
!              the plume edge-) (RDL_plume)
! t_re_flag=1: new formulation for contact time weighted on FL 
!              concentrations along the trajectory (NHRDL, NI-NHRDL)
! t_re_flag=2: contact time=fly time (since a virtual release from the 
!              inlet section) with inlet section at x=0 (RDL) 
! t_re_flag=3: contact time=fly time (since a virtual release from the 
!              inlet section) with inlet section at x=0, computed at 
!              the cell centre (no fluctuation) (NFRDL)
               if (t_re_flag.eq.1) then
! Frozen limit formulas
                  C_A0 = F_m * C_A1
                  C_B0 = (1.d0 - F_m) * C_B2
                  sum_C_A0_C_B0 = sum_C_A0_C_B0 + C_A0 * C_B0
                  if ((C_A0.ne.0.d0).and.(C_B0.ne.0.d0)) then
                     time_react_eff = sum_C_A0_C_B0 * dt / (C_A0 * C_B0)
                     else
                        time_react_eff = time
                  endif
                  elseif (t_re_flag.eq.0) then
                     time_react_eff = time
                     elseif (t_re_flag.eq.2) then
                        time_react_eff = xpr / U_scale
                        elseif (t_re_flag.eq.3) then
                           x_cell = dx * (inew - 1 + 0.5d0)
                           time_react_eff = x_cell / U_scale
               endif
! Other formulas of the conserved scalar theory
               if(F_m.eq.0.d0)then
                  C_A = 0.d0
                  C_B = C_B2
                  elseif(F_m.eq.1.d0) then
                     C_A = C_A1
                     C_B = 0.d0
                     else
                        if (F_m.ne.F_ms) then
                           C_A = C_A1 * F_m * (F_ms - F_m) / (F_ms * 
     &                           (1.d0 - F_m) * dexp((F_ms - F_m) * 
     &                           reaction_rate * (C_A1 + C_B2) * 
     &                           time_react_eff) - F_m * (1.d0 - F_ms))
                           C_B = C_B2 * ((F_ms - F_m) / F_ms + (1.d0 - 
     &                           F_ms) * C_A / (F_ms * C_A1))
                           else
                              C_A = (C_A1 * F_ms) / (1.d0 + F_ms * (1.d0
     &                              - F_ms) * reaction_rate * (C_A1 + 
     &                              C_B2) * time_react_eff)
                              C_B = C_A
                        endif
               endif
      endif
! Check concentration values
      if (C_A.gt.C_A1) then
         if ((C_A.gt.(1.05d0*C_A1)).and.(n_warn_C_A_max.le.1000)) then
            write(*,*) 
     &         "Warning: C_A>1.05*C_A1; F_m, time, ix_cel, iy_cel,",
     &         " iz_cel, ID_part, Fm_m_c, Fm_old, C_A, C_B:"
            write(*,*) F_m,time,inew,jnew,knew,ipcont,
     &         concc(inew,jnew,knew,ipu,jpv,kpw),cpr,C_A,C_B
            n_warn_C_A_max = n_warn_C_A_max + 1
            if (n_warn_C_A_max.eq.1000) then
               write(*,*) 
     &            "Reached the maximum number (1000) of warnings ",
     &            "on C_A_max"
            endif
         endif
         C_A = C_A1
      endif
      if (C_A.lt.0.0d0) then
         if ((C_A.lt.(-(0.05d0)*C_A1)).and.(n_warn_C_A_min.le.1000)) 
     &      then
            write(*,*) 
     &         "Warning: C_A<-0.05C_A1; F_m, time, ix_cel, iy_cel,",
     &         " iz_cel, ID_part, Fm_m_c, Fm_old, C_A, C_B:"
            write(*,*) F_m,time,inew,jnew,knew,ipcont,
     &         concc(inew,jnew,knew,ipu,jpv,kpw),cpr,C_A,C_B
            n_warn_C_A_min = n_warn_C_A_min + 1
            if (n_warn_C_A_min.eq.1000) then
               write(*,*) 
     &            "Reached the maximum number (1000) of warnings ",
     &            "on C_A_min"
            endif
         endif
         C_A = 0.d0
      endif
      if (C_B.gt.C_B2) then
         if ((C_B.gt.(1.05d0*C_B2)).and.(n_warn_C_B_max.le.1000)) then
            write(*,*) 
     &         "Warning: C_B>1.05*C_B2; F_m, time, ix_cel, iy_cel,",
     &         " iz_cel, ID_part, Fm_m_c, Fm_old, C_A, C_B:"
            write(*,*) F_m,time,inew,jnew,knew,ipcont,
     &         concc(inew,jnew,knew,ipu,jpv,kpw),cpr,C_A,C_B
            n_warn_C_B_max = n_warn_C_B_max + 1
            if (n_warn_C_B_max.eq.1000) then
               write(*,*) 
     &            "Reached the maximum number (1000) of warnings ",
     &            "on C_B_max"
            endif
         endif
         C_B = C_B2
      endif
      if (C_B.lt.0.0) then
         if ((C_B.lt.(-(0.05d0)*C_B2)).and.(n_warn_C_B_min.le.1000)) 
     &      then 
            write(*,*) 
     &         "Warning: C_B<-0.05*C_B2; F_m, time, ix_cel, iy_cel,",
     &         " iz_cel, ID_part, Fm_m_c, Fm_old, C_A, C_B:"
            write(*,*) F_m,time,inew,jnew,knew,ipcont,
     &         concc(inew,jnew,knew,ipu,jpv,kpw),cpr,C_A,C_B
            n_warn_C_B_min = n_warn_C_B_min + 1
            if (n_warn_C_B_min.eq.1000) then
               write(*,*) 
     &            "Reached the maximum number (1000) of warnings ",
     &            "on C_B_min"
            endif
         endif
         C_B = 0.d0
      endif
! Estimation of the instantaneous concentration of the product 
! species
      C_C = F_m * C_A1 - C_A
      if (C_C.lt.0.0d0) then
         if ((C_C.lt.(-(0.05d0)*C_B2)).and.(n_warn_C_C_min.le.1000)) 
     &      then
            write(*,*) 
     &         "Warning: C_C<-0.05C_B2; F_m, time, ix_cel, iy_cel,",
     &         " iz_cel, ID_part, Fm_m_c, Fm_old, C_A, C_B, C_C: "
            write(*,*) F_m,time,inew,jnew,knew,ipcont,
     &         concc(inew,jnew,knew,ipu,jpv,kpw),cpr,C_A,C_B,C_C
            n_warn_C_C_min = n_warn_C_C_min + 1
            if (n_warn_C_C_min.eq.1000) then
               write(*,*) 
     &            "Reached the maximum number (1000) of warnings ",
     &            "on C_C_min"
            endif
         endif
         C_C = 0.d0
      endif
! Updating summations for the reactive concentration statistics
      conc_react1(inew,jnew,knew) = conc_react1(inew,jnew,knew) + C_A
      conc_react2(inew,jnew,knew) = conc_react2(inew,jnew,knew) + C_B
      conc_product(inew,jnew,knew) = conc_product(inew,jnew,knew) + C_C
      sigc_react1(inew,jnew,knew) = sigc_react1(inew,jnew,knew) + C_A **
     &                              2
      sigc_react2(inew,jnew,knew) = sigc_react2(inew,jnew,knew) + C_B **
     &                              2
      sigc_product(inew,jnew,knew) = sigc_product(inew,jnew,knew) + C_C
     &                               ** 2
      c_skew_react1(inew,jnew,knew) = c_skew_react1(inew,jnew,knew) + 
     &                                C_A ** 3
      c_skew_react2(inew,jnew,knew) = c_skew_react2(inew,jnew,knew) + 
     &                                C_B ** 3 
      c_skew_product(inew,jnew,knew) = c_skew_product(inew,jnew,knew) + 
     &                                 C_C ** 3
      c_kurt_react1(inew,jnew,knew) = c_kurt_react1(inew,jnew,knew) +
     &                                C_A ** 4
      c_kurt_react2(inew,jnew,knew) = c_kurt_react2(inew,jnew,knew) +
     &                                C_B ** 4
      c_kurt_product(inew,jnew,knew) = c_kurt_product(inew,jnew,knew) +
     &                                 C_C ** 4
      Iseg(inew,jnew,knew) = Iseg(inew,jnew,knew) + C_A * C_B
      return
      end 
      
