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
! Subroutine name: Sawfordfull
! Subroutine description: 1D Lagrangian Macro-Mixing. E' equivalente a 
! "Sawford2006", con tutte le particelle inquinate, per verificare la 
! condizione di completo mescolamento. Year: 2006.
!-----------------------------------------------------------------------
      PROGRAM Sawfordfull
	  REAL U, n, zs, ztop, zbot, xend, deltaz, steptmax, dt, C0, pgreco,stepzmax, sigma2wsor,i,j,mSaw,deltaW,classeWmax,Wclassemin,sigma2w22,sigma2wt0,r,s,TL0,mp,Volcella
      REAL C (100,22)
!     stepzmax=100, steptmax=22
	  REAL dZ, X, out, Z, W, dW, mdeps,mW,alea1W,pdfW,alea2maxW,alea2W,alea1deps,pdfdeps,alea2maxdeps,alea2deps,deps
	  REAL Zmat(10,22), depsmat(10,22)
!     osservo 10 particelle, stetmax=22
	  REAL Wvet(1000)
!     vedo la velocità iniziale di 1000 particelle
	  real sigmazvet(22)
!     vettore delle sigmaz numeriche, stept=22
	  real mZvet(22)
!     vettore della media delle quote, stept=22
	  real Zquadrovet(22)
!     vettore della somma degli scarti quadratici medi delle quote, stept=22
      real nW(100,40)
!     matrice della media delle C in funzione di W e Z, (stepz=100,stepwmax=40) 
      real pdfWmat(100,40)
!     matrice della pdf numerica di W in al variare di Z, (stepz=100,stepwmax=40)  
      real CdataW(100,40)
!     matrice della media delle C data una classe di W, (stepz=100,stepwmax=40) 
      real argerf(100,40)
!     parametro di verifica analitica per CdataW (stepz=100,stepwmax=40) 
      real  fN(40)
!     fN è il vettore della frazione di realizzazioni ad una certa W, dato il tempo, classeWmax=40
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
!	  WRITE (*,*) 'Inserisci il numero di classi di velocità (classeWmax): '
!	  WRITE (*,*) 'Inserisci l'intervallo di una classe di velocità (deltaW): '
!	  WRITE (*,*) 'Inserisci la velocità minima classificata (Wclassemin): '
!	  WRITE (*,*) 'Inserisci la massa di inquinante in microgrammi al metro cubo in una particella inquinata (mp): '
!	  READ (*,*) n, U, zs, ztop, xend, steptmax, dt, sigma2w,stepzmax,C0,t0,mSaw,classeWmax,deltaW,Wclassemin,mp
	  n=1000000
	  U=6.2
	  zs=0.10
	  ztop=0.10
	  zbot=-0.10
      xend=1000001
	  steptmax=22
	  dt=0.00403
	  sigma2wsor=0.049
	  stepzmax=100
	  C0=3
	  t0=0.168
	  mSaw=1.26
	  classeWmax=40
	  deltaW=0.034
	  Wclassemin=-0.68
	  mp=1
	  
	  
	  mW=0
	  deltaz=(ztop-zbot)/stepzmax
	  mdeps=0
	  Volcella=U*dt*deltaZ
	  sigma2w22=sigma2wSaw(sigma2wsor,t0,mSaw,steptmax,dt)
	  open (12,file='i')
	  open (13,file='j')
	  CALL Cinizialifull(C,n,stepzmax,mp,Volcella)
	  DO 100, j=1,n
	     WRITE (13,*) j
         i=1
	     Z=zbot+(j-1)*((zs-zbot)/n)
		 if (j.LE.10) Zmat(j,1)=Z
		 CALL estragauss(alea1W,mW,sigma2wsor,pdfW,alea2maxW,alea2W,W)
		 if (j.LE.1000) Wvet(j)=W
		 out=0
	     DO 110, i=2,steptmax
		    WRITE (12,*) i 
		    IF (out.EQ.0) THEN
               CALL estragauss (alea1deps,mdeps,dt,pdfdeps,alea2maxdeps,alea2deps,deps)
			   CALL calcoloC(U,dW,W,dt,sigma2wsor,deps,dZ,Z,X,i,j,ztop,zbot,xend,deltaz,C,depsmat,Zmat,stepzmax,C0,mSaw,t0,mp,Volcella,n)
               if (i.eq.22) then
			   call calcolonW(stepzmax,classeWmax,Z,zbot,deltaz,W,Wclassemin,deltaW,nW)
			   end if
			END IF
110   CONTINUE
100   CONTINUE
	  call calcolopdfWmat(stepzmax,classeWmax,pdfWmat,nW,C,Volcella,n,mp)
	  call calcoloCdataW(stepzmax,classeWmax,CdataW,mp,nW,Volcella,n,fN,zbot,Wclassemin,deltaW,sigma2w22,sigma2wsor,t0,mSaw,dt)
	  call risultati(U,n,zs,ztop,zbot,xend,dt,deltaz,out,steptmax,sigma2wsor,C,Zmat,depsmat,Wvet,stepzmax,mSaw,C0,t0,deltaW,classeWmax,Wclassemin,nW,sigma2w22,TL0,CdataW,fN,mp,Volcella,pdfWmat)
	  CLOSE (12,STATUS='KEEP')
	  CLOSE (13,STATUS='KEEP')
	  END



      SUBROUTINE Cinizialifull(C,n,stepzmax,mp,Volcella)
      REAL C(100,22)
	  real n,stepzmax,mp,Volcella
	  DO 200, i=1,stepzmax
              C(i,1)=((n/stepzmax)*mp)/(n*Volcella)
200   CONTINUE
      return
	  END




      subroutine estragauss (alea1,m,sigma2,pdf,alea2max,alea2,ris)
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
		    


      subroutine calcoloC(U,dW,W,dt,sigma2wsor,deps,dZ,Z,X,i,j,ztop,zbot,xend,deltaz,C,depsmat,Zmat,stepzmax,C0,mSaw,t0,mp,Volcella,n)
	  REAL C(100,22)
	  REAL Zmat(10,22)
	  REAL depsmat(10,22)
	  real U,dW,W,dt,sigma2wsor,deps,dZ,Z,X,i,j,ztop,zbot,xend,deltaz,stepzmax,C0,mSaw,t0,n,mp,Volcella
	  call calcolodW(dW,W,U,dt,C0,deps,mSaw,t0,i,sigma2wsor)
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
	     IF ((Z.GT.(zbot+deltaz*(k-1))) .AND. (Z.LE.(zbot+deltaz*k))) C(k,i)=C(k,i)+mp/(Volcella*n)
400      CONTINUE
	  return
      end



      subroutine calcolodW (dW,W,U,dt,C0,deps,mSaw,t0,i,sigma2wsor)
	  real dW,W,U,dt,C0,deps,mSaw,t0,dsigma2wdt,i,sigma2wsor
	  sigma2w=sigma2wSaw(sigma2wsor,t0,mSaw,i,dt)
	  dsigma2wdt=dsigma2wdtSaw(mSaw,t0,sigma2wsor,i,dt)
	  dW=-C0*(-3/2)*dsigma2wdt/(2*sigma2w)*W*dt+W*dsigma2wdt*dt/(2*sigma2w)+sqrt(C0*(-3/2)*dsigma2wdt)*deps
      return
	  end



      subroutine calcolonW(stepzmax,classeWmax,Z,zbot,deltaz,W,Wclassemin,deltaW,nW)
	  real stepzmax,classeWmax,Z,zbot,deltaz,W,Wclassemin,deltaW
	  real nW(100,40)
	  real a,b
      do 800, a=1,stepzmax
	     do 900, b=1,classeWmax
		    if (((Z.GT.(zbot+deltaz*(a-1))) .AND. (Z.LE.(zbot+deltaz*a))).and.((W.GT.(Wclassemin+deltaW*(b-1))) .AND. (W.LE.(Wclassemin+deltaW*b)))) then
			nW(a,b)=nW(a,b)+1
		    end if
900   continue
800   continue
      return
	  end


      
	  subroutine calcolopdfWmat(stepzmax,classeWmax,pdfWmat,nW,C,Volcella,n,mp)
      real stepzmax,classeWmax,Volcella,n,mp
	  real pdfWmat(100,40)
	  real nW(100,40)
	  real C(100,22)
	  real a,b
	  do 1100, a=1,stepzmax
	     do 1110, b=1,classeWmax
		    pdfWmat(a,b)=nW(a,b)/(C(a,22)*Volcella*n/mp)
1110     continue
1100  continue
      return
      end   
	  
	  
	  
	  subroutine calcoloCdataW(stepzmax,classeWmax,CdataW,mp,nW,Volcella,n,fN,zbot,Wclassemin,deltaW,sigma2w22,sigma2wsor,t0,mSaw,dt)
      real stepzmax,classeWmax,mp,Volcella,n,zbot,Wclassemin,sigma2w22,sigma2wsor,t0,mSaw,dt
	  real CdataW(100,40)
	  real nW(100,40)
	  real fN(40)
	  real a,b
	  do 1000, b=1,classeWmax
	     WfN=Wclassemin+deltaW*(b-1)+deltaW/2
		 fN(b)=gauss(WfN,sigma2w22)*deltaW
1000  continue
	  do 1010, a=1,stepzmax
	     do 1020, b=1,classeWmax
            CdataW(a,b)=(mp*nW(a,b))/(Volcella*n*fN(b))
1020     continue
1010  continue
      return
	  end



      subroutine risultati (U,n,zs,ztop,zbot,xend,dt,deltaz,out,steptmax,sigma2wsor,C,Zmat,depsmat,Wvet,stepzmax,mSaw,C0,t0,deltaW,classeWmax,Wclassemin,nW,sigma2w22,TL0,CdataW,fN,mp,Volcella,pdfWmat)
	  REAL U,n,zs,ztop,zbot,xend,dt,deltaz,out,steptmax,sigma2wsor,stepzmax,mSaw,C0,t0,deltaW,classeWmax,Wclassemin,sigma2w22,TL0,mp,Volcella
	  REAL C(100,22)
	  REAL Zmat(10,22)
	  REAL depsmat(10,22)
	  REAL Wvet(1000)
	  real nW(100,40)
	  real CdataW(100,40)
	  real pdfWmat(100,40)
	  real fN(40)
	  OPEN (1,file='generale')
	  OPEN (2,FILE='C')
	  OPEN (3,FILE='Zmat')
	  OPEN (4,FILE='depsmat')
	  OPEN (5,FILE='Wvet')
	  OPEN (6,FILE='nW')
	  open (7,file='CdataW')
	  open (8,file='fN')
	  open (9,file='pdfWmat')
	  WRITE (1,FMT=500) U,n,zs,ztop,zbot,xend,dt,deltaz,out,steptmax,sigma2wsor,stepzmax,mSaw,C0,t0,deltaW,classeWmax,Wclassemin,sigma2w22,TL0,mp,Volcella
	  WRITE (2,FMT=510) C
	  WRITE (3,FMT=520) Zmat
	  WRITE (4,FMT=520) depsmat
	  WRITE (5,FMT=530) Wvet
	  WRITE (6,FMT=510) nW
	  WRITE (7,FMT=510) CdataW
	  WRITE (8,FMT=530) fN
	  WRITE (9,FMT=510) pdfWmat
500   FORMAT(1X,24(F11.6,1X))
510   FORMAT(1X,100(F10.3,1X))
520   FORMAT(1X,10(F10.3,1X))
530   FORMAT (1X,F10.6,1X)
      CLOSE (1,STATUS='KEEP')
	  CLOSE (2,STATUS='KEEP')
	  CLOSE (3,STATUS='KEEP')
	  CLOSE (4,STATUS='KEEP')
	  CLOSE (5,STATUS='KEEP')
	  CLOSE (6,STATUS='KEEP')
	  CLOSE (7,STATUS='KEEP')
	  CLOSE (8,STATUS='KEEP')
	  CLOSE (9,STATUS='KEEP')
      return
      end




      real function sigma2wSaw(sigma2wsor,t0,mSaw,i,dt)
	  real sigma2wsor,t0,mSaw,i,dt
	  sigma2wSaw=sigma2wsor*((1+dt*i/t0)**(-mSaw))
	  return
	  end



      real function gauss(a,sigma2)
      PARAMETER (pgreco=3.14159265)
	  gauss=(1/(SQRT(2*sigma2*pgreco)))*(EXP(-(0.5)*(((ABS(a))**2)/sigma2)))
	  return
	  end



      real function dsigma2wdtSaw(mSaw,t0,sigma2wsor,i,dt)
	  real sigma2wsor,t0,mSaw,i,dt
	  dsigma2wdtSaw=(-mSaw/t0)*sigma2wsor*((1+dt*i/t0)**(-mSaw-1))
	  return
	  end
