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
! Program unit name: Thomsoncos
! Program unit description: 1D Lagrangian Macro-Mixing. Applicazione  
! dello schema di Thomson (1984), con deviazione standard della  
! componente verticale (lungo asse z) della velocitÓ dipende da cos(z). 
! Year: 2006.
!-----------------------------------------------------------------------
      PROGRAM Thomsoncos
!     applicazione di Thomson del 1984, turbolenza verticale dipende da cos(z)
      
	  REAL U, n, zs, ztop, xend, deltaz, stepmax, dt, TL, pgreco,stepzmax, sigma2wmax,gamma
      REAL C (40,3000)
!     stepzmax=40, steptmax=3000
	  REAL dZ, X, out, Z, W, dW, mdeps,mW,alea1W,pdfW,alea2maxW,alea2W,alea1deps,pdfdeps,alea2maxdeps,alea2deps,deps
	  REAL Zmat(10,3000), depsmat(10,3000)
!     osservo 10 particelle, stetmax=3000
	  REAL Wvet(1000)
!     vedo la velocitÓ iniziale di 1000 particelle 
!	  WRITE (*,*) 'Inserisci il numero di particelle (n): '
!	  WRITE (*,*) 'Inserisci la velocitÓ media orizzontale del fluido (U): '
!	  WRITE (*,*) 'Inserisci la quota di rilascio delle particelle (zs)'
!	  WRITE (*,*) 'Inserisci l''altezza dello strato mescolato scalare (ztop): '
!	  WRITE (*,*) 'Inserisci la lunghezza dello strato mescolato scalare (xend): '
!	  WRITE (*,*) 'Inserisci il numero di step temporali (steptmax) da simulare: '
!	  WRITE (*,*) 'Inserisci lo step temporale (dt): '
!	  WRITE (*,*) 'Inserisci la varianza massima della velocitÓ verticale (sigma2wmax): '
!	  WRITE (*,*) 'Inserisci la varianza minima della velocitÓ verticale (sigma2wmin): '
!	  WRITE (*,*) 'Inserisci la scala di decorrelazione delle velocitÓ verticali (TL, valore di riferimento: 1000s): '
!	  WRITE (*,*) 'Inserisci il numero di step spaziali verticali (stepzmax): '
!	  WRITE (*,*) 'Inserisci il fattore di disomogeneitÓ gamma: '
!	  READ (*,*) n, U, zs, ztop, xend, stepmax, dt, sigma2w, TL, stepzmax,gamma
	  n=100000
	  U=1
	  zs=1000
	  ztop=2000
      xend=1000001
	  steptmax=3000
	  dt=5
	  TL=1000
	  stepzmax=40
	  sigma2wmax=1.96
	  gamma=1
	  
	  mW=0
	  deltaz=ztop/stepzmax
	  mdeps=0
	  CALL Ciniziali(C,n)
!	  open (6,file='verificaj')
!	  open (7,file='verificai')
!     verifico a che punto si ferma la simulazione
700         format(f6.0)
	  DO 100, j=1,n
	     Z=zs
		 if (j.LE.10) Zmat(j,1)=Z
		 CALL estragauss(alea1W,mW,sigma2wmax,pdfW,alea2maxW,alea2W,W)
		 if (j.LE.1000) Wvet(j)=W
		 out=0
!		 if (j.GT.1000) write (6,fmt=700) j
		 DO 110, i=2,steptmax
			IF (out.EQ.0) THEN
!			   if (j.EQ.9105) write (7,fmt=700) i
!			   if ((Z.gt.2000).or.(z.lt.0)) then
!			   pause
!			   end if
!			   if ((j.EQ.9105).and.(i.eq.1837)) then
!			   pause
!			   endif
               CALL estragauss (alea1deps,mdeps,dt,pdfdeps,alea2maxdeps,alea2deps,deps)
    		   CALL calcoloC (U,dW,W,dt,TL,deps,dZ,Z,X,i,j,ztop,xend,deltaz,C,depsmat,Zmat,stepzmax,gamma)
		    END IF
110   CONTINUE
100   CONTINUE
      call risultati (U,n,zs,ztop,xend,dt,deltaz,out,steptmax,sigma2wmax,TL,alea2maxdeps,C,Zmat,depsmat,Wvet,stepzmax,gamma)
!	  close (6,status='keep')
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



      real function gauss(a,sigma2)
      PARAMETER (pgreco=3.14159265)
	  gauss=(1/(SQRT(2*sigma2*pgreco)))*(EXP(-(0.5)*(((ABS(a))**2)/sigma2)))
	  return
	  end



      subroutine estragauss (alea1,m,sigma2,pdf,alea2max,alea2,ris)
	  REAL alea1,m,sigma2,pdf,alea2max,alea2,ris
	  alea2max=gauss(m,sigma2)
300   CALL RANDOM (ranval)
		   alea1=(ranval-0.5)*2*(3*SQRT(sigma2))
           pdf=gauss(alea1,sigma2)
		   CALL RANDOM (ranval)
		   alea2=ranval*alea2max
		   IF (alea2.LE.pdf) THEN 
		   ris=alea1
		   ELSE
		   GO TO 300
           END IF
	  return
	  END			    
		    


      subroutine calcoloC (U,dW,W,dt,TL,deps,dZ,Z,X,i,j,ztop,xend,deltaz,C,depsmat,Zmat,stepzmax,gamma)
	  REAL C(40,3000)
	  REAL Zmat(10,3000)
	  REAL depsmat(10,3000)
      real sigma2w,Z
	  PARAMETER (pgreco=3.14159265)
	  sigma2w=sigma2wThomsoncos(Z,ztop,TL,gamma)
	  dW=-W*dt/TL+(SQRT(abs(2*sigma2w/TL)))*deps-(2*((1.2)+cos(z*pgreco/ztop)*gamma*ztop/(10*TL))*(gamma*pgreco/(10*TL))*sin(z*pgreco/ztop))*(sigma2w+(W**2))*(dt)/(2*sigma2w)
	  W=W+dW
	  dZ=W*dt
	  Z=Z+dZ
	  X=U*dt*(i-1)
600	  IF (Z.gt.ztop) then
	     Z=ztop-(Z-ztop)
		 W=-W
 		 go to 600
	  end if
	  IF (Z.lt.0) then
	     Z=-Z
		 W=-W
		 go to 600
	  end if
	  IF (X.gt.xend) out=1
	  IF (j.LE.10) then
	     Zmat(j,i)=Z
		 depsmat(j,i)=deps
	  end if
	  DO 400, k=1,stepzmax
	  IF ((Z.GT.deltaz*(k-1)) .AND. (Z.LE.deltaz*k)) C(k,i)=C(k,i)+1
400   CONTINUE
      return
      end



      real function sigma2wThomsoncos(Z,ztop,TL,gamma)
	  real Z,ztop,gamma
	  PARAMETER (pgreco=3.14159265)
	  b=(pgreco/ztop)*Z
	  a=cos(b)
	  sigma2wThomsoncos=(1.2+a*gamma*ztop/(10*TL))**2
	  return
	  end


      subroutine risultati (U,n,zs,ztop,xend,dt,deltaz,out,stepmax,sigma2wmax,TL,alea2maxdeps,C,Zmat,depsmat,Wvet,stepzmax,gamma)
	  REAL n
	  REAL C(40,3000)
	  REAL Zmat(10,3000)
	  REAL depsmat(10,3000)
	  REAL Wvet(1000)
	  OPEN (1,file='generale')
	  OPEN (2,FILE='C')
	  OPEN (3,FILE='Zmat')
	  OPEN (4,FILE='depsmat')
	  OPEN (5,FILE='Wvet')
	  WRITE (1,FMT=500) U,n,zs,ztop,xend,dt,deltaz,out,stepmax,sigma2wmax,TL,alea2maxdeps,stepzmax,gamma
	  WRITE (2,FMT=510) C
	  WRITE (3,FMT=520) Zmat
	  WRITE (4,FMT=520) depsmat
	  WRITE (5,FMT=530) Wvet
500   FORMAT(1X,14(F11.3,1X))
510   FORMAT(1X,40(F10.3,1X))
520   FORMAT(1X,10(F10.3,1X))
530   FORMAT (1X,F10.3,1X)
      CLOSE (1,STATUS='KEEP')
	  CLOSE (2,STATUS='KEEP')
	  CLOSE (3,STATUS='KEEP')
	  CLOSE (4,STATUS='KEEP')
	  CLOSE (5,STATUS='KEEP')
      return
      end






