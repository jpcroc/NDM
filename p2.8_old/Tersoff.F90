!     nvar=3*nb d'atome
!     ndim = dimension des tableaux d'atomes
!     xvar(3*nbd'at)=positions
!     epot=energie potentielle 
!     fvar(ndim) =forces
!     ivoisi =nb des voisins  de i
!     listvoi(i,1->IVOIS(i)) = indices des voisins de i 
!     
!     *******************************************************************
!     ********************  FORCES ET POTENTIEL DE TERSOFF **************
!     **** (La structuration de cette subroutine provient d'une *********
!     **** version 2D-periodique de Pietro Ballone 10 octobre 91) *******
!     ****  corrigee par Olivier HARDOUIN DUPARC  05 janvier 1993  ****** 
!     ****  simplifiee pour cluster Francois Willaime 21 Juin 1996 ******
!     *******************************************************************
!     ******************************************************************** 
      SUBROUTINE Tersoff(nvar,ndim,xvar,epot,fvar,ivoisi, &
          listvoi,zl,zls2)
!     
!     ***** NDIM=Nombre max d'atomes *****
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
!     implicit none
      dimension DZEKX(NDIM),DZEKY(NDIM),DZEKZ(NDIM),EAT(NDIM), &
          ivoisi(ndim),listvoi(ndim,15)
      dimension xvar(ndim*3),fvar(ndim*3)
      real*8 zl(3),zls2(3)
!     
!     ***** Constantes du potentiel ***** 
!     parametres pour le carbone:
!     J.Tersoff PRL vol.25 (1988) 2879
      DATA AA/1.3936D3/,BB/3.4674D2/ !
      DATA RL1/3.4879D0/,RL2/2.2119D0/ ! A-1
      DATA RR/1.95D0/,DD/0.15D0/ ! A
      DATA BETA/1.5724D-7/,ENNE/0.72751D0/
      DATA SC/3.8049D4/,SD/4.3484/,SH/-5.7058D-1/

!     ***** Initialisations *****
!     
      PI=4.d0*datan(1.d0)
      nat=nvar/3
!     print *,'natom=',nat
      UN=1.d0
      ZERO=0.d0
      HALF=0.5d0
      DEUX=2.d0
      RMD    = RR - DD
      RPD    = RR + DD
      RPD2   = RPD**2
      SC2    = SC**2
      SD2    = SD**2
      AUS1   = UN + (SC2/SD2)
      TWOSC2 = DEUX*SC2
      EXPON  = - UN / (DEUX*ENNE)
      PIMDD  = HALF*PI/DD
      do i=1,nvar
         fvar(i)=0.d0
      enddo
!     
!     ***** Debut boucle sur tous les atomes I *****
!     

!     boucle sur les atomes
      DO 1 I = 1, NAT  
!     
         EAT(I) = ZERO
         xi = xvar(3*i-2)
         yi = xvar(3*i-1)
         zi = xvar(3*i)
!     
!     ***** Debut boucle sur les voisins de I *****
!     boucle sur les voisins j<>i

         DO 2 JJ = 1, IVOISI(I)
!     do 2 j=1,nat
!     if(j.eq.i) goto 2
            J = LISTVOI(I,JJ)
!     
            DZIJ = xvar(3*j) - ZI
!     rc	  if(dabs(dzij).gt.rpd) goto 2
            if((dzij-ZLS2(3)).gt.0.0) dzij=dzij-ZL(3)
            if((dzij+ZLS2(3)).lt.0.0) dzij=dzij+ZL(3)

            DXIJ = xvar(3*j-2) - XI
            if((dxij-ZLS2(1)).gt.0.0) dxij=dxij-ZL(1)
            if((dxij+ZLS2(1)).lt.0.0) dxij=dxij+ZL(1)

!     rc	  if(dabs(dxij).gt.rpd) goto 2
            DYIJ = xvar(3*j-1) - YI
            if((dyij-ZLS2(2)).gt.0.0) dyij=dyij-ZL(2)
            if((dyij+ZLS2(2)).lt.0.0) dyij=dyij+ZL(2)

!     rc	  if(dabs(dyij).gt.rpd) goto 2
            DR2IJ = DXIJ**2+DYIJ**2+DZIJ**2
!     rpd2 rayon de coupure sur rij

                                !          IF (DR2IJ.GT.RPD2) then
                                !             write(6,*)'ceci ne doit pas arriver !',dr2ij,rpd2
                                !             stop
                                !          endif
!     

            DRIJ  = DSQRT(DR2IJ)
!     print *,i,j,drij
            IF (DRIJ.LT.RMD) THEN
               FC    =  UN
               DFCDR =  ZERO
            ELSE
               ARG   =  PIMDD*(DRIJ-RR)
               FC    =  HALF * ( UN-DSIN(ARG) )
               DFCDR = -HALF * PIMDD * DCOS(ARG)
            END IF
            FR = +AA*DEXP(-RL1*DRIJ)
            FA = -BB*DEXP(-RL2*DRIJ)
!     
            ZETAIJ = ZERO
            DZEDXI = ZERO
            DZEDYI = ZERO
            DZEDZI = ZERO
            DZEDXJ = ZERO
            DZEDYJ = ZERO
            DZEDZJ = ZERO
            DO 29 KK = 1, IVOISI(I)
               DZEKX(KK) = ZERO
               DZEKY(KK) = ZERO
               DZEKZ(KK) = ZERO
 29         CONTINUE 
!     *************************** loop over k **********************
            DO 3 KK = 1, IVOISI(I)
               K = LISTVOI(I,KK)
               IF (K.EQ.J) GO TO 3
!     
               DZIK = xvar(3*k) - ZI
               if((dzik-ZLS2(3)).gt.0.0) dzik=dzik-ZL(3)
               if((dzik+ZLS2(3)).lt.0.0) dzik=dzik+ZL(3)

               DXIK = xvar(3*k-2) - XI
               if((dxik-ZLS2(1)).gt.0.0) dxik=dxik-ZL(1)
               if((dxik+ZLS2(1)).lt.0.0) dxik=dxik+ZL(1)

               DYIK = xvar(3*k-1) - YI
               if((dyik-ZLS2(2)).gt.0.0) dyik=dyik-ZL(2)
               if((dyik+ZLS2(2)).lt.0.0) dyik=dyik+ZL(2)

               DR2IK = DXIK**2+DYIK**2+DZIK**2
                                !            IF (DR2IK.GT.RPD2) then
                                !               write(6,*)'ceci ne doit pas arriver K!',dr2ik,rpd2
                                !               stop
                                !            endif

!     
               DRIK = DSQRT(DR2IK)
!     
               DRIJIK = UN/(DRIJ*DRIK)
               COST   = (DXIJ*DXIK+DYIJ*DYIK+DZIJ*DZIK)*(DRIJIK)
               AUSIL  = SD2+(SH-COST)**2
               GTETA  = AUS1-SC2/AUSIL
               IF (DRIK.LT.RMD) THEN
                  FCP  =  UN
                  DFCP =  ZERO
               ELSE
                  ARG  =  PIMDD*(DRIK-RR)
                  FCP  =  HALF * ( UN-DSIN(ARG) )
                  DFCP = -HALF * PIMDD * DCOS(ARG)
               END IF
               DGDCOS = -TWOSC2*(SH-COST)/AUSIL**2
               CRT    = DFCP*GTETA/DRIK
!     Contribution en vue de F(I) :
               DCDRIJ=(COST/DRIJ - UN/DRIK)/DRIJ
               DCDRIK=(COST/DRIK - UN/DRIJ)/DRIK
               DGTDXI=DGDCOS*(DCDRIJ*DXIJ+DCDRIK*DXIK)
               DGTDYI=DGDCOS*(DCDRIJ*DYIJ+DCDRIK*DYIK)
               DGTDZI=DGDCOS*(DCDRIJ*DZIJ+DCDRIK*DZIK)
               DZEDXI=DZEDXI+FCP*DGTDXI-CRT*DXIK
               DZEDYI=DZEDYI+FCP*DGTDYI-CRT*DYIK
               DZEDZI=DZEDZI+FCP*DGTDZI-CRT*DZIK
!     Contribution en vue de F(J) :
               COR2IJ=COST/DR2IJ
               DGTDXJ=DGDCOS*(COR2IJ*(-DXIJ)+DRIJIK*DXIK)
               DGTDYJ=DGDCOS*(COR2IJ*(-DYIJ)+DRIJIK*DYIK)
               DGTDZJ=DGDCOS*(COR2IJ*(-DZIJ)+DRIJIK*DZIK)
               DZEDXJ=DZEDXJ+FCP*DGTDXJ
               DZEDYJ=DZEDYJ+FCP*DGTDYJ
               DZEDZJ=DZEDZJ+FCP*DGTDZJ
!     Contribution en vue de F(K) :
               COR2IK=COST/DR2IK
               DGTDXK=DGDCOS*(COR2IK*(-DXIK)+DRIJIK*DXIJ)
               DGTDYK=DGDCOS*(COR2IK*(-DYIK)+DRIJIK*DYIJ)
               DGTDZK=DGDCOS*(COR2IK*(-DZIK)+DRIJIK*DZIJ)
               DZEKX(KK)=FCP*DGTDXK-CRT*(-DXIK)
               DZEKY(KK)=FCP*DGTDYK-CRT*(-DYIK)
               DZEKZ(KK)=FCP*DGTDZK-CRT*(-DZIK)
!     
               ZETAIJ=ZETAIJ+FCP*GTETA
!     
!     ************** end of loop over k ******************************
 3          CONTINUE
            IF(ZETAIJ.GT.1.D-10) THEN
               BZENNE = (BETA*ZETAIJ)**ENNE
               BODY   = UN + BZENNE
               BIJ    = (BODY)**(EXPON)
               DBIJ   = -HALF*BZENNE*BIJ/(ZETAIJ*BODY)
            ELSE
               BIJ    = 1.D0
               DBIJ   = 0.D0
            ENDIF
            FRBFA = FR + BIJ*FA
            EAT(I)= EAT(I) + FC*FRBFA
            RLRBA = RL1*FR + RL2*BIJ*FA 
            DVDR  = (DFCDR*FRBFA - FC*RLRBA)/DRIJ
            FCADB = FC*FA*DBIJ
            DVDXI = - DVDR*DXIJ + FCADB*DZEDXI
            DVDYI = - DVDR*DYIJ + FCADB*DZEDYI
            DVDZI = - DVDR*DZIJ + FCADB*DZEDZI
            fvar(3*i-2) = fvar(3*i-2) - DVDXI
            fvar(3*i-1) = fvar(3*i-1) - DVDYI
            fvar(3*i)   = fvar(3*i)   - DVDZI
            DVDXJ = - DVDR*(-DXIJ) + FCADB*DZEDXJ
            DVDYJ = - DVDR*(-DYIJ) + FCADB*DZEDYJ
            DVDZJ = - DVDR*(-DZIJ) + FCADB*DZEDZJ
            fvar(3*j-2) = fvar(3*j-2) - DVDXJ
            fvar(3*j-1) = fvar(3*j-1) - DVDYJ
            fvar(3*j)   = fvar(3*j)   - DVDZJ
!     do 33 k=1,nat
!     if (k.eq.i) goto 33
            DO 33 KK = 1, IVOISI(I)
               K = LISTVOI(I,KK)
               IF (K.EQ.J) GO TO 33
               fvar(3*k-2) = fvar(3*k-2) - FCADB*DZEKX(KK)
               fvar(3*k-1) = fvar(3*k-1) - FCADB*DZEKY(KK)
               fvar(3*k)   = fvar(3*k)   - FCADB*DZEKZ(KK)
 33         CONTINUE
 2       CONTINUE
 1    CONTINUE
!     
!     
!     ***** Valeur finale de l'energie potentielle (en eV) *****
!     ***** et valeurs finales des forces en eV/A *****
      EPOT = 0.D0
      DO  I = 1, NAT
         EPOT  = EPOT+EAT(I)
      ENDDO
      do i=1,nvar
         fvar(i)=0.5d0*fvar(i)
      enddo
      EPOT = EPOT/2.D0
      end
!
!                                   ***** Fin routine TSPOTENTIEL *****
