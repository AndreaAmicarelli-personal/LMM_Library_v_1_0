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
! Subroutine name: meteo_input
! Subroutine description: Pre-processing of grid turbulence main flow 
! data. Year: 2014 
!-----------------------------------------------------------------------
PROGRAM meteo_input
!Declarations
implicit none
integer(4) :: ndx,ndy,ndz,i,j,k,i_thresh
!Don't use double precision as the meteo arrays in LAGFLUM are real
real :: U,c,x,y,z,dx,M,A,x_g,c1,c2,A1,A2
real, allocatable, dimension(:,:,:) :: mu,mv,mw
real, allocatable, dimension(:,:,:) :: varu,varv,varw,eps
!Initializations (hard coding)
ndx=220
ndy=5
ndz=68
U = 0.5
dx = 0.032
M = 10.* dx
x_g = -0.41*M
i_thresh = 126
 c1 = 1.7
 c2 = 1.27
A1 = 3.014
A2 = 1.
! Allocations
allocate(mu(ndx,ndy,ndz))
allocate(mv(ndx,ndy,ndz))
allocate(mw(ndx,ndy,ndz))
allocate(varu(ndx,ndy,ndz))
allocate(varv(ndx,ndy,ndz))
allocate(varw(ndx,ndy,ndz))
allocate(eps(ndx,ndy,ndz))
mu = U
mv = 0.
mw = 0.
! Formatted files
open(1,file='u.txt')
open(2,file='v.txt')
open(3,file='w.txt')
open(4,file='sigmau.txt')
open(5,file='sigmav.txt')
open(6,file='sigmaw.txt')
open(7,file='eps.txt')
do i=1,ndx
   if (i<=i_thresh) then
      c=c1
      A=A1
      else
         c=c2
         A=A2
   endif
   x = dx * (i-1+0.5)
   do j=1,ndy
      y = dx * (j-1+0.5)
      do k=1,ndz
         z = dx * (k-1+0.5)
         varu(i,j,k) = (0.05*U)**2
         varv(i,j,k) = A * 0.041 * U**2 * ((x-x_g)/M)**(-c)
         varw(i,j,k) = varv(i,j,k)
         eps(i,j,k) = varw(i,j,k) * 1.5 * c * U / (x-x_g)
         write(1,'(4(e12.5))') x,y,z,mu(i,j,k)
         write(2,'(4(e12.5))') x,y,z,mv(i,j,k)
         write(3,'(4(e12.5))') x,y,z,mw(i,j,k)
         write(4,'(4(e12.5))') x,y,z,sqrt(varu(i,j,k))
         write(5,'(4(e12.5))') x,y,z,sqrt(varv(i,j,k))
         write(6,'(4(e12.5))') x,y,z,sqrt(varw(i,j,k))
         write(7,'(4(e12.5))') x,y,z,eps(i,j,k)
      enddo
   enddo
enddo
close(1)
close(2)
close(3)
close(4)
close(5)
close(6)
close(7)
! Unformatted files
open(1,file='u.dat',form='unformatted')
open(2,file='v.dat',form='unformatted')
open(3,file='w.dat',form='unformatted')
open(4,file='varu.dat',form='unformatted')
open(5,file='varv.dat',form='unformatted')
open(6,file='varw.dat',form='unformatted')
open(7,file='eps.dat',form='unformatted')
write(1) mu
write(2) mv
write(3) mw
write(4) varu
write(5) varv
write(6) varw
write(7) eps
close(1)
close(2)
close(3)
close(4)
close(5)
close(6)
close(7)
!Deallocations
deallocate(mu)
deallocate(mv)
deallocate(mw)
deallocate(varu)
deallocate(varv)
deallocate(varw)
deallocate(eps)
end program meteo_input
