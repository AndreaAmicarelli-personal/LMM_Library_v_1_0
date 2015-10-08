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
! Program unit name: WienersenzCBL
! Program unit description: processo di Wiener. Deviazione standard 
! della velocità verticale è sinusoidale (CBL), dt=5s, rilascio puntuale
! e istantaneo. Year: 2006.
!-----------------------------------------------------------------------
      PROGRAM WienersenzCBL
 
	  REAL U, D, n, sigma2wmax,TL,sigma2w
      REAL zs, ztop, xend, deltaz, stepmax, dt, alfa
      REAL i,j,k
	  REAL C (40,3000)
!     stepzmax=40, steptmax=3000
	  REAL Yu, dmu, dZ, X, out, Z,stepzmax
	  REAL Zmat(10,3000)
!     vedo 10 traiettorie,steptmax=3000
!	  WRITE (*,*) 'Inserisci il numero di particelle (n): '
!	  WRITE (*,*) 'Inserisci la velocità media orizzontale del fluido (U): '
!	  WRITE (*,*) 'Inserisci la varianza massima delle velocità verticali (sigma2wmax): '
!	  WRITE (*,*) 'Inserisci la varianza minima delle velocità verticali (sigma2wmin): '
!	  WRITE (*,*) 'Inserisci la quota di rilascio delle particelle (zs)'
!	  WRITE (*,*) 'Inserisci l''altezza dello strato mescolato scalare (ztop): '
!	  WRITE (*,*) 'Inserisci la lunghezza dello strato mescolato scalare (xend): '
!	  WRITE (*,*) 'Inserisci il numero di step temporali (steptmax) da simulare: '
!	  WRITE (*,*) 'Inserisci lo step temporale (dt): '
!	  WRITE (*,*) 'Inserisci il numero di step spaziali verticali(stepzmax): '
!	  WRITE (*,*) 'Inserisci la scala di decorrelazione delle velocità verticali (TL): '
!	  READ (*,*) n, U, sigma2wmax, zs, ztop, xend, steptmax, dt,stepzmax,TL
      n=100000
	  U=1
	  sigma2wmax=1.60
	  sigma2wmin=0.24
	  zs=1000
	  ztop=2000
	  xend=1000001
	  steptmax=3000
	  dt=5
	  stepzmax=40
	  TL=1000
	       	 	  
	  deltaz=ztop/stepzmax
	  CALL Ciniziali(C,n)
	  out=0
	  DO 100, j=1,n
	     Z=zs
		 if (j.LE.10) Zmat(j,1)=Z
		 out=0
	     DO 110, i=2,steptmax
		    IF (out.EQ.0) THEN
               call calcolodZ(Yu,dmu,dt,alfa,D,dZ,Z,X,TL,U,i,ztop,sigma2wmax,sigma2wmin)
			   CALL calcoloC (ztop,xend,deltaz,C,Zmat,stepzmax,Z,X,out,j,i)
		    END IF
110   CONTINUE
100   CONTINUE
      call risultati (U,n,zs,ztop,xend,dt,deltaz,out,steptmax,stepzmax,sigma2wmax,sigma2wmin,TL,C,Zmat)
	  END

		    

      SUBROUTINE Ciniziali(C,n)
      REAL C(40,3000)
	  real n
	  real stepzmax
	  stepzmax=40
	  DO 200, i=1,stepzmax
	         C(i,1)=0
200   CONTINUE
      i=stepzmax/2  
	  C(i,1)=n
	  return
	  END
	    


      subroutine calcolodZ(Yu,dmu,dt,alfa,D,dZ,Z,X,TL,U,i,ztop,sigma2wmax,sigma2wmin)
	  real Yu,dmu,dt,alfa,D,dZ,Z,X,i,ztop,sigma2wmax,sigma2wmin
	  CALL RANDOM (ranval)
			   Yu=ranval
			   dmu=(Yu-0.5)*SQRT(12*dt)
			   sigma2w=sigma2wfun(Z,sigma2wmax,ztop,sigma2wmin)
			   D=sigma2w*TL
			   alfa=SQRT(abs(2*D))
			   dZ=alfa*dmu
			   Z=Z+dZ
			   X=U*dt*(i-1)
	  return
	  end


      real function sigma2wfun(Z,sigma2wmax,ztop,sigma2wmin)
	  real Z,sigma2wmax,ztop
	  PARAMETER (pgreco=3.14159265)
	  b=(pgreco/ztop)*Z
	  a=sin(b)
	  sigma2wfun=(sigma2wmax-sigma2wmin)*a+sigma2wmin
	  return
	  end



      subroutine calcoloC (ztop,xend,deltaz,C,Zmat,stepzmax,Z,X,out,j,i)
	  real j,i
	  REAL C(40,3000)
	  REAL Zmat(10,3000)
500   IF (Z>ztop) then
	     Z=ztop-(Z-ztop)
	     go to 500
	  end if
	  IF (Z<0) then
	     Z=-Z
	     go to 500
	  end if
	  IF (X>xend) out=1
	  IF (j.LE.10) Zmat(j,i)=Z
	  DO 300, k=1,stepzmax
	  IF ((Z.GT.deltaz*(k-1)) .AND. (Z.LE.deltaz*k)) C(k,i)=C(k,i)+1
300   CONTINUE
      return
      end


	  
      subroutine risultati (U,n,zs,ztop,xend,dt,deltaz,out,steptmax,stepzmax,sigma2wmax,sigma2wmin,TL,C,Zmat)
	  REAL U,n,zs,ztop,xend,dt,deltaz,out,steptmax,stepzmax,sigma2wmax,TL
	  REAL C(40,3000)
	  REAL Zmat(10,3000)
	  open (1,file='generale')
	  OPEN (2,FILE='C')
	  OPEN (3,FILE='Zmat')
	  WRITE (1,FMT=400) U,n,zs,ztop,xend,dt,deltaz,out,steptmax,sigma2wmax,sigma2wmin,TL,stepzmax
	  WRITE (2,FMT=410) C
	  WRITE (3,FMT=420) Zmat
400   FORMAT(1X,13(F11.3,1X))
410   FORMAT(1X,40(F10.3,1X))
420   FORMAT(1X,10(F10.3,1X))
      CLOSE (1,STATUS='KEEP')
	  CLOSE (2,STATUS='KEEP')
	  CLOSE (3,STATUS='KEEP')
      return
      end
      
