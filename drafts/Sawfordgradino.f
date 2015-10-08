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
! Subroutine name: Sawfordgradino
! Subroutine description: 1D Lagrangian Macro-Mixing. Modello di Thomson
! 1987 applicato al caso Sawford 2006, in decaying grid turbulence. 
! Year: 2006.
!-----------------------------------------------------------------------
      PROGRAM Sawfordgradino
!     dati
      REAL U, n, zs, ztop, zbot, xend, deltaz, steptmax, dt, C0, pgreco,stepzmax, sigma2wsor,i,j,mSaw,mp,Volcella
      REAL np (40,3000)
!     np è il numero di particelle inquinate in una cella, stepzmax=40, steptmax=3000
      REAL C (40,3000)
!     C è la media d'insieme della concentrazione di inquinante in una cella, stepzmax=40, steptmax=3000
	  REAL dZ, X, out, Z, W, dW, mdeps,mW,alea1W,pdfW,alea2maxW,alea2W,alea1deps,pdfdeps,alea2maxdeps,alea2deps,deps
	  REAL Zmat(10,3000), depsmat(10,3000)
!     osservo 10 particelle, stetmax=3000
	  REAL Wvet(1000)
!     vedo la velocità iniziale di 1000 particelle
!	  WRITE (*,*) 'Inserisci il numero di particelle (n): '
!	  WRITE (*,*) 'Inserisci la velocità media orizzontale del fluido (U): '
!	  WRITE (*,*) 'Inserisci la quota di rilascio delle particelle (zs)'
!	  WRITE (*,*) 'Inserisci la quota massima (ztop): '
!	  WRITE (*,*) 'Inserisci la quota minima (zbot): '
!	  WRITE (*,*) 'Inserisci la lunghezza dello strato mescolato scalare (xend): '
!	  WRITE (*,*) 'Inserisci il numero di step temporali (steptmax) da simulare: '
!	  WRITE (*,*) 'Inserisci lo step temporale (dt): '
!	  WRITE (*,*) 'Inserisci la varianza della velocità verticale alla sorgente (sigma2wsor): '
!	  WRITE (*,*) 'Inserisci il numero di step spaziali verticali (stepzmax): '
!	  WRITE (*,*) 'Inserisci la costante C0: '
!	  WRITE (*,*) 'Inserisci il tempo t0: '
!	  WRITE (*,*) 'Inserisci il parametro di decadimento mSaw: '
!	  WRITE (*,*) 'Inserisci la massa di inquinante in microgrammi al metro cubo in una particella inquinata (mp): '
!	  READ (*,*) n, U, zs, ztop, xend, steptmax, dt, sigma2w,stepzmax,C0,t0,mSaw,mp
	  n=500000
	  U=6.2
	  zs=0
	  ztop=0.10
	  zbot=-0.10
      xend=1000001
	  steptmax=3000
	  dt=0.00403
	  sigma2wsor=0.049
	  stepzmax=40
	  C0=3
	  t0=0.168
	  mSaw=1.26
	  mp=1
	  
!     simulazione	  
	  mW=0
	  deltaz=(ztop-zbot)/stepzmax
	  mdeps=0
	  Volcella=U*dt*deltaZ
	  CALL npinizialigradino(np,n)
	  DO 100, j=1,n
	     Z=-(j-1)*((zs-zbot)/n)
		 if (j.LE.10) Zmat(j,1)=Z
		 CALL estragauss(alea1W,mW,sigma2wsor,pdfW,alea2maxW,alea2W,W)
		 if (j.LE.1000) Wvet(j)=W
		 out=0
	     DO 110, i=2,steptmax
		    IF (out.EQ.0) THEN
               CALL estragauss (alea1deps,mdeps,dt,pdfdeps,alea2maxdeps,alea2deps,deps)
    		   CALL calcoloW (U,dW,W,dt,sigma2wsor,deps,dZ,Z,X,i,j,ztop,zbot,xend,deltaz,np,depsmat,Zmat,stepzmax,C0,mSaw,t0,mp,Volcella)
		    END IF
110   CONTINUE
100   CONTINUE
      call calcoloC(C,mp,np,Volcella,n,stepzmax,steptmax)
      call risultati (U,n,zs,ztop,zbot,xend,dt,deltaz,out,steptmax,sigma2wsor,C,Zmat,depsmat,Wvet,stepzmax,mSaw,C0,t0,mp,Volcella)
	  END


      SUBROUTINE npinizialigradino(np,n)
!     definisce il numero di particelle per cella iniziali
      REAL np(40,3000)
	  real n
	  real stepzmax
	  stepzmax=40
	  DO 200, i=1,stepzmax
             if (i.le.20) then
			    np(i,1)=n/20
				else
				np(i,1)=0
		     end if
200   CONTINUE
      return
	  END



      real function gauss(a,sigma2)
!     restituisce la pdf gaussiana di media "a" e varianza "sigma2"
      PARAMETER (pgreco=3.14159265)
	  gauss=(1/(SQRT(2*sigma2*pgreco)))*(EXP(-(0.5)*(((ABS(a))**2)/sigma2)))
	  return
	  end



      subroutine estragauss (alea1,m,sigma2,pdf,alea2max,alea2,ris)
!     estrae una variabile aleatoria (ris) con pdf gaussiana di media "m" e varianza "sigma2"
	  REAL m
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
		    


      subroutine calcoloW (U,dW,W,dt,sigma2wsor,deps,dZ,Z,X,i,j,ztop,zbot,xend,deltaz,np,depsmat,Zmat,stepzmax,C0,mSaw,t0,mp,Volcella)
!     calcola le traiettorie, con riflessioni geometriche al top e bottom, e il numero di particelle per cella 
	  REAL np(40,3000)
	  REAL Zmat(10,3000)
	  REAL depsmat(10,3000)
	  real U,dW,W,dt,sigma2wsor,deps,dZ,Z,X,i,j,ztop,zbot,xend,deltaz,stepzmax,C0,mSaw,t0,mp,Volcella
	  call calcolodW (dW,W,U,dt,C0,deps,mSaw,t0,i,sigma2wsor)
	  W=W+dW
	  dZ=W*dt
	  Z=Z+dZ
	  X=U*dt*(i-1)
600	  IF (Z.gt.ztop) then
	     Z=2*ztop-Z
		 W=-W
		 go to 600
	  end if
	  IF (Z.lt.zbot) then
	     Z=2*zbot-Z
		 W=-W
		 go to 600
	  end if
	  IF (X.gt.xend) out=1
	  IF (j.LE.10) then
	     Zmat(j,i)=Z
		 depsmat(j,i)=deps
	  end if
	  DO 400, k=1,stepzmax
	  IF ((Z.GT.(zbot+deltaz*(k-1))) .AND. (Z.LE.(zbot+deltaz*k))) np(k,i)=np(k,i)+1
400   CONTINUE
      return
      end



      subroutine calcolodW (dW,W,U,dt,C0,deps,mSaw,t0,i,sigma2wsor)
!     calcola le velocità verticali secondo Thomson 1987 1D con disomogeneità orizzontale
	  real dW,W,U,dt,C0,deps,mSaw,t0,dsigma2wdt,i,sigma2wsor
	  sigma2w=sigma2wSaw(sigma2wsor,t0,mSaw,i,dt)
	  dsigma2wdt=dsigma2wdtSaw(mSaw,t0,sigma2wsor,i,dt)
	  dW=-C0*(-3/2)*dsigma2wdt/(2*sigma2w)*W*dt+W*dsigma2wdt*dt/(2*sigma2w)+sqrt(C0*(-3/2)*dsigma2wdt)*deps
      return
	  end


      real function sigma2wSaw(sigma2wsor,t0,mSaw,i,dt)
!     dà la turbolenza secondo la legge di potenza di Sawford 2006
	  real sigma2wsor,t0,mSaw,i,dt
	  sigma2wSaw=sigma2wsor*((1+dt*i/t0)**(-mSaw))
	  return
	  end



      real function dsigma2wdtSaw(mSaw,t0,sigma2wsor,i,dt)
!     calcola la derivata sostanziale della turbolenza nel tempo
	  real sigma2wsor,t0,mSaw,i,dt
	  dsigma2wdtSaw=(-mSaw/t0)*sigma2wsor*((1+dt*i/t0)**(-mSaw-1))
	  return
	  end


      subroutine calcoloC(C,mp,np,Volcella,n,stepzmax,steptmax)
!     calcola le concentrazioni
	  real mp,Volcella,n
      real np(40,3000)
	  real C(40,3000)
	  real i,k
      do 700, i=1,steptmax
	     do 800 k=1,stepzmax
		    C(k,i)=(mp*np(k,i))/(Volcella*n)
800      continue
700   continue 
      return
      end


      subroutine risultati (U,n,zs,ztop,zbot,xend,dt,deltaz,out,steptmax,sigma2wsor,C,Zmat,depsmat,Wvet,stepzmax,mSaw,C0,t0,mp,Volcella)
!     output
	  REAL U,n,zs,ztop,zbot,xend,dt,deltaz,out,steptmax,sigma2wsor,stepzmax,mSaw,C0,t0,mp,Volcella
	  REAL C(40,3000)
	  REAL Zmat(10,3000)
	  REAL depsmat(10,3000)
	  REAL Wvet(1000)
	  OPEN (1,file='generale')
	  OPEN (2,FILE='C')
	  OPEN (3,FILE='Zmat')
	  OPEN (4,FILE='depsmat')
	  OPEN (5,FILE='Wvet')
	  WRITE (1,FMT=500) U,n,zs,ztop,zbot,xend,dt,deltaz,out,steptmax,sigma2wsor,stepzmax,mSaw,C0,t0,mp,Volcella
	  WRITE (2,FMT=510) C
	  WRITE (3,FMT=520) Zmat
	  WRITE (4,FMT=520) depsmat
	  WRITE (5,FMT=530) Wvet
500   FORMAT(1X,17(F11.3,1X))
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






