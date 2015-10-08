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
! Program unit name: WieneromogeneoCBL
! Program unit description: processo di Wiener. Come Wiener.f, ma con 
! in CBL, dt 10 volte + piccolo e steptmax 3 volte + grande. Year: 2006.
!-----------------------------------------------------------------------
      PROGRAM WieneromogeneoCBL

	  REAL U, D, n
      REAL zs, ztop, xend, deltaz, stepmax, dt, alfa
      REAL i,j,k
	  REAL C (40,3000)
!     stepzmax=40, steptmax=1000
	  REAL Yu, dmu, dZ, X, out, Z,stepzmax
	  REAL Zmat(10,3000)
!     vedo 10 traiettorie,steptmax=1000
!	  WRITE (*,*) 'Inserisci il numero di particelle (n): '
!	  WRITE (*,*) 'Inserisci la velocità media orizzontale del fluido (U): '
!	  WRITE (*,*) 'Inserisci la dispersività verticale (D), unica presente: '
!	  WRITE (*,*) 'Inserisci la quota di rilascio delle particelle (zs)'
!	  WRITE (*,*) 'Inserisci l''altezza dello strato mescolato scalare (ztop): '
!	  WRITE (*,*) 'Inserisci la lunghezza dello strato mescolato scalare (xend): '
!	  WRITE (*,*) 'Inserisci il numero di step temporali (steptmax) da simulare: '
!	  WRITE (*,*) 'Inserisci lo step temporale (dt): '
!	  WRITE (*,*) 'Inserisci il numero di step spaziali verticali(stepzmax): '
!	  READ (*,*) n, U, D, zs, ztop, xend, steptmax, dt,stepzmax
      n=100000
	  U=1
	  D=1110
	  zs=1000
	  ztop=2000
	  xend=1000001
	  steptmax=3000
	  dt=10
	  stepzmax=40
	       	 	  
	  alfa=SQRT(2*D)
	  CALL SEED (1)
      deltaz=ztop/stepzmax
	  DO 500, i=1,stepzmax
	         C(i,1)=0
500   CONTINUE
      i=stepzmax/2  
	  C(i,1)=n
	  out=0
	  DO 100, j=1,n
	     Z=zs
		 IF (j.LE.10) Zmat(j,1)=Z
	     DO 110, i=2,steptmax
		    IF (out.EQ.0) THEN
			   CALL RANDOM (ranval)
			   Yu=ranval
			   dmu=(Yu-0.5)*SQRT(12*dt)
			   dZ=alfa*dmu
			   Z=Z+dZ
			   X=U*dt*(i-1)
			END IF
600     	IF (Z.gt.ztop) then
			   Z=ztop-(Z-ztop)
			   go to 600
			end if
			IF (Z.lt.0) then 
			   Z=-Z
			   go to 600
			end if
			IF (X>xend) out=1
			IF (j.LE.10) Zmat(j,i)=Z
			DO 130, k=1,40
		   	   IF ((Z.GT.deltaz*(k-1)) .AND. (Z.LE.deltaz*k)) C(k,i)=C(k,i)+1
130   CONTINUE
110   CONTINUE
100   CONTINUE
      
      OPEN (1,file='generale')
	  OPEN (2,FILE='C')
	  OPEN (3,FILE='Zmat')
	  WRITE (1,FMT=400) U,D,n,zs,ztop,xend,dt,alfa,deltaz,out,steptmax,stepzmax
	  WRITE (2,FMT=410) C
	  WRITE (3,FMT=420) Zmat
400   FORMAT(1X,12(F11.3,1X))
410   FORMAT(1X,40(F10.3,1X))
420   FORMAT(1X,10(F10.3,1X))
      CLOSE (1,STATUS='KEEP')
	  CLOSE (2,STATUS='KEEP')
	  CLOSE (3,STATUS='KEEP')
	  	  	  
	  END

