module calfo3c_mod
        use cryst_to_cart_mod
        use notperiod_mod
        implicit none 
        contains
! *****************************************************************
subroutine calfo3c(xp,  vp,  fp, ielat, iwmax, ityp)
  !version du 20.11.2001
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use var_pot
  implicit none
  !-----------------------------------------------
  !   G l o b a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  integer , intent(in) :: ielat(imm)
  integer  :: iwmax(imm)
  integer , intent(in) :: ityp(imm)
  real(double) , intent(inout) :: xp(3,imm)
  real(double)  :: vp(3,imm)
  real(double)  :: ax(3,imm)
  real(double) , intent(inout) :: fp(3,imm)
  !-----------------------------------------------
  !   L o c a l   P a r a m e t e r s
  !-----------------------------------------------





  integer :: i,koo,i1,ko1,ko3,j,i2,k,i3, &
       iti,itj,itk,itrip,lj,lk,ic

  real(double) :: &
       x1,x2,x3,c1,c2,c3,r2,        &
       pscal,pscal2,pror2,pror,tetjik,cosi,cosi2,c2osi, &
       tcp31,tcp32, &
       tcp1,tcp2,tcp3,tcp4,tcp5,tcpx1,tcpx2,tcpy1, &
       tcpy2,tcpz1,tcpz2,cv(1,3), &
       Rayij,Rayik,inv_Rij,inv_Rik, &
       INTER1j,INTER1k,INTER2exp,INTER2,INTER3j,INTER3k

  real (double), dimension(26*natperc) :: &
       rtc,xtc, ytc, ztc, indic



  real(double) &
       w1,w1c,w2,w2c,w3,w3c, &
       uxij,uyij,uzij,uxik,uyik,uzik, &
       rij,riji,rij2,rik,riki,rik2, &
       rrijk(6,3)
  integer :: Deb
  integer :: Fin,ncelvois


  real(double) :: xpnp(3,imm)
  REAL(double), dimension(1:3) :: dxp


  Deb=1
  Fin=im
  potcp=0




  !C-----------------------------------------------------
  !C                 Triplets j-i-k
  !C-----------------------------------------------------
  !C --- ouverture de la boucle sur i

  if (noxyz==1) then
     call cryst_to_cart (imm, xp, bg, -1)    !cart vers cryst

  else
     if (lperiod) then
        xpnp(:,:)=xp(:,:)
     else 
        call notperiod(xp,xpnp)
     end if

  end if

      if (any(free).NEQV..true.)then
         write(6,*)'free +SW =pas code'
         stop
      end if



  DO  I=Deb,Fin
!write(6,*)'debut',i
     ITI=ityp(i)
     I3=0

     if (noxyz==1)then
        do j=1,im
           if (j==i) cycle

           dxp(1:3) = xp(1:3,i) - xp(1:3,j)

           WHERE ( (dxp.GT.0.5d0).OR.(dxp.LT.-0.5d0) )
              dxp(1:3) = dxp(1:3) - Dble(Nint(dxp(1:3)))
           END WHERE

           dxp = MatMul(at,dxp)

           r2 = Sum( dxp(1:3)**2 )


           if(r2.le.r3cm2) then
              I3=I3+1
              RTC(I3)=R2            !distance rij**2
              XTC(I3)=dxp(1)          !Xi-Xj
              YTC(I3)=dxp(2)            !Yi-Yj
              ZTC(I3)=dxp(3)            !Zi-Zj
              INDIC(I3)=J           !# de l'atome j
!           write(6,*)i,j,r2
           endif

        end do
     else
        X1=XPnp(1,I)
        X2=XPnp(2,I)
        X3=XPnp(3,I)
        !C --- calcul du # de cellule KOO de l'atome i
        KOO=ielat(i)

        !C --- calcul du tableau MERKEN des voisins du i considere
        ncelvois = min(noxyz,27)-1
        ! pour chaque cel. voisine
        do i1 = 0, ncelvois

           !     DO  I1=0,26
           KO1=NCEL(KOO,I1)
           DO  I2=1,NATO(KO1)
              j=last(i2,ko1)

              if(i.eq.j) cycle
              !C --- calcul de la distance (=>rij)
              C1=X1-XPnp(1,J)
              C2=X2-XPnp(2,J)
              C3=X3-XPnp(3,J)
              !           !c--- conditions aux limites

              do ic=1,3
                 c1=c1+at(1,ic)*deltadist(ic,i1,koo)
                 c2=c2+at(2,ic)*deltadist(ic,i1,koo)
                 c3=c3+at(3,ic)*deltadist(ic,i1,koo)
              end do

              if (noxyz.ne.1) then
                 if (abs(c1)>r3cm2) cycle
                 if (abs(c2)>r3cm2) cycle
                 if (abs(c3)>r3cm2) cycle 
              else
                 cv(1,1) = c1
                 cv(1,2) = c2
                 cv(1,3) = c3
                 call cryst_to_cart (1, cv, bg, -1) !cryst vers cart sur cv
                 WHERE ( (cv.GT.0.5d0).OR.(cv.LT.-0.5d0) )
                    cv(:,1:3) = cv(:,1:3) - Dble(Nint(cv(:,1:3)))
                 END WHERE
                 call cryst_to_cart (1, cv, at, 1) !cryst vers cart sur cv
                 c1=cv(1,1)
                 c2=cv(1,2)
                 c3=cv(1,3)
              end if

!              if (abs(C1).gt.r3cm) cycle
!              if (abs(C2).gt.r3cm) cycle
!              if (abs(C3).gt.r3cm) cycle




              R2=C1*C1+C2*C2+C3*C3
              !C --- stockage des valeurs de distance pour l'atome # j
              if(r2.le.r3cm2) then
                 I3=I3+1
                 RTC(I3)=R2            !distance rij**2
                 XTC(I3)=C1            !Xi-Xj
                 YTC(I3)=C2            !Yi-Yj
                 ZTC(I3)=C3            !Zi-Zj
                 INDIC(I3)=J           !# de l'atome j
              endif
           end DO
        end DO

     end if

!write(6,*)'milieu',i

     !C --- fin de la boucle sur les paires

     !C --- Traitement du triplet JIK
     !C --- ouverture de la boucle sur les j

     DO  I1=1,I3-1
        !C --- ouverture de la boucle sur les k
        DO  I2=I1+1,I3

           J=INDIC(I1)
           K=INDIC(I2)
           !write(6,*)i1,i2,i3,j,k

           ITJ=ITYP(J)
           ITK=ITYP(K)
           ITRIP=IPO3C(ITI,ITJ,ITK)
           LJ=IPO(ITI,ITJ)
           LK=IPO(ITI,ITK)
           LJ=IPO(ITI,ITJ)
           if(RTC(I1).GT.COUP3C2(ITRIP,LJ)) cycle
           LK=IPO(ITI,ITK)
           if(RTC(I2).GT.COUP3C2(ITRIP,LK)) cycle

           !C --- calcul du produit scalaire rij.rik
           PSCAL=XTC(I1)*XTC(I2)+YTC(I1)*YTC(I2)+ZTC(I1)*ZTC(I2)
           PSCAL2=PSCAL*PSCAL
           PROR2=RTC(I1)*RTC(I2)
           !write(6,*)'tota'
           PROR=SQRT(PROR2)
           !write(6,*)'tota2'
           !write(6,*)rtc(i1),rtc(i2),lj,lk
           !C --- calcul de l'angle jik
           TETJIK=PSCAL/PROR     ! tetijk=cos_angle ijk
           if (TETJIK.lt.-1) TETJIK=-1.0
           !C --- calcul de termes intermediaires
           Rayij=SQRT(RTC(i1))
           !write(6,*)'tota3'
           Rayik=SQRT(RTC(i2))
           !write(6,*)'tota4'
           inv_Rij=1/(Rayij-coup3c(itrip,lj))
           !write(6,*)'tota5'
           inv_Rik=1/(Rayik-coup3c(itrip,lk))
           !write(6,*)'tota6'
           INTER1j=gam(itrip,lj)*inv_Rij**2/Rayij
           !write(6,*)'tota7'
           INTER1k=gam(itrip,lk)*inv_Rik**2/Rayik
           !write(6,*)'tota8'
           INTER2exp=gam(itrip,lj)*inv_Rij+gam(itrip,lk)*inv_Rik
           INTER2=0.
           !write(6,*)'toto'
           if(INTER2exp.gt.(-1./precexp))&
                INTER2=lamb(itrip)* &
                EXP(INTER2exp)
           INTER3j=TETJIK*Rayij**2
           INTER3k=TETJIK*Rayik**2

           !C --- calcul des termes cosiaires pour la force
           COSI=TETJIK-cangle(itrip) !difference d'angle
           !crc           COSI2=COSI*COSI
           !crc          C2OSI=2.*COSI

           COSI2=COSI*COSI/(1+c3c(itrip)*COSI*COSI)
           C2OSI=2.*(COSI/(1+c3c(itrip)*COSI*COSI))*(1.-c3c(itrip)*cosi*cosi/(1+c3c(itrip)*cosi*cosi))


           !C --- calcul de l'energie potentiel du terme a 3 corps
           POTCP=POTCP+INTER2*COSI2
!           if ((i==1).or.(j==1).or.(k==1))then
              !write(6,'(3I5,2F15.5,D15.5)')i,j,k,rayij*1d8,rayik*1d8,inter2*cosi2*erg2ev
!           end if

           !C --- calcul des differents terme intervenant dans le calcul des forces
           TCP31=INTER3j/PROR2
           TCP32=INTER3k/PROR2
           PROR=PROR/PROR2
           TCP1=INTER2*C2OSI
           TCP2=INTER2*INTER1j
           TCP3=INTER2*INTER1k
           TCP4=(-TCP31+PROR)
           TCP5=(-TCP32+PROR)
           TCPX1=XTC(I1)*TCP2*COSI2
           TCPX2=XTC(I2)*TCP3*COSI2
           TCPY1=YTC(I1)*TCP2*COSI2
           TCPY2=YTC(I2)*TCP3*COSI2
           TCPZ1=ZTC(I1)*TCP2*COSI2
           TCPZ2=ZTC(I2)*TCP3*COSI2

           !C             Force sur i

           FP(1,I)=FP(1,I)+TCPX1+TCPX2-TCP1* &
                (XTC(I1)*TCP5+XTC(I2)*TCP4)
           FP(2,I)=FP(2,I)+TCPY1+TCPY2-TCP1* &
                (YTC(I1)*TCP5+YTC(I2)*TCP4)
           FP(3,I)=FP(3,I)+TCPZ1+TCPZ2-TCP1* &
                (ZTC(I1)*TCP5+ZTC(I2)*TCP4)

           !C             Force sur j et k

           FP(1,J)=FP(1,J)-TCPX1+TCP1* &
                (-XTC(I1)*TCP32+XTC(I2)*PROR)
           FP(2,J)=FP(2,J)-TCPY1+TCP1* &
                (-YTC(I1)*TCP32+YTC(I2)*PROR)
           FP(3,J)=FP(3,J)-TCPZ1+TCP1* &
                (-ZTC(I1)*TCP32+ZTC(I2)*PROR)
           FP(1,K)=FP(1,K)-TCPX2+TCP1* &
                (XTC(I1)*PROR-XTC(I2)*TCP31)
           FP(2,K)=FP(2,K)-TCPY2+TCP1* &
                (YTC(I1)*PROR-YTC(I2)*TCP31)
           FP(3,K)=FP(3,K)-TCPZ2+TCP1* &
                (ZTC(I1)*PROR-ZTC(I2)*TCP31)

           !c-------------Calculate 3body contribution to stress tensor, sig
           if (itesigma.gt.0) then
              if(mod(it,itesigma).eq.0) then
                 uxij=XTC(I1)
                 uyij=YTC(I1)   ! avant: uyij=XTC(I1)
                 uzij=ZTC(I1)   ! avant: uzij=XTC(I1)
                 rij2=RTC(I1)
                 rij=SQRT(rij2)
                 riji=1.0/rij
                 uxik=XTC(I2)
                 uyik=YTC(I2)   ! uyik=XTC(I2)
                 uzik=ZTC(I2)   ! uzik=XTC(I2)
                 rik2=RTC(I2)
                 rik=SQRT(rik2)
                 riki=1.0/rik
                 w1=(-TCP2*rij2*COSI2-TCP1*tetjik)/volu
                 w2=TCP1/volu
                 w3=(-TCP3*rik2*COSI2-TCP1*tetjik)/volu
                 w1c=(-TCP2*rij2*COSI2-TCP1*tetjik)*noxyz/volu
                 w2c=TCP1*noxyz/volu
                 w3c=(-TCP3*rik2*COSI2-TCP1*tetjik)*noxyz/volu

                 rrijk(1,1)=uxij*riji*uxij*riji
                 rrijk(1,2)=2.e0*uxij*riji*uxik*riki
                 rrijk(1,3)=uxik*riki*uxik*riki
                 rrijk(2,1)=uyij*riji*uyij*riji
                 rrijk(2,2)=2.e0*uyij*riji*uyik*riki
                 rrijk(2,3)=uyik*riki*uyik*riki
                 rrijk(3,1)=uzij*riji*uzij*riji
                 rrijk(3,2)=2.e0*uzij*riji*uzik*riki
                 rrijk(3,3)=uzik*riki*uzik*riki
                 rrijk(4,1)=uyij*riji*uzij*riji
                 rrijk(4,2)=uyij*riji*uzik*riki+uyik*riki*uzij*riji
                 rrijk(4,3)=uyik*riki*uzik*riki
                 rrijk(5,1)=uzij*riji*uxij*riji
                 rrijk(5,2)=uzij*riji*uxik*riki+uzik*riki*uxij*riji
                 rrijk(5,3)=uzik*riki*uxik*riki
                 rrijk(6,1)=uxij*riji*uyij*riji
                 rrijk(6,2)=uxij*riji*uyik*riki+uxik*riki*uyij*riji
                 rrijk(6,3)=uxik*riki*uyik*riki

                 sig(1,1)=sig(1,1)-rrijk(1,1)*w1-rrijk(1,2)*w2 &
                      &                         -rrijk(1,3)*w3
                 sig(1,2)=sig(1,2)-rrijk(6,1)*w1-rrijk(6,2)*w2 &
                      &                         -rrijk(6,3)*w3
                 sig(1,3)=sig(1,3)-rrijk(5,1)*w1-rrijk(5,2)*w2 &
                      &                         -rrijk(5,3)*w3
                 sig(2,1)=sig(2,1)-rrijk(6,1)*w1-rrijk(6,2)*w2 &
                      &                         -rrijk(6,3)*w3
                 sig(2,2)=sig(2,2)-rrijk(2,1)*w1-rrijk(2,2)*w2 &
                      &                         -rrijk(2,3)*w3
                 sig(2,3)=sig(2,3)-rrijk(4,1)*w1-rrijk(4,2)*w2 &
                      &                         -rrijk(4,3)*w3
                 sig(3,1)=sig(3,1)-rrijk(5,1)*w1-rrijk(5,2)*w2 &
                      &                         -rrijk(5,3)*w3
                 sig(3,2)=sig(3,2)-rrijk(4,1)*w1-rrijk(4,2)*w2 &
                      &                         -rrijk(4,3)*w3
                 sig(3,3)=sig(3,3)-rrijk(3,1)*w1-rrijk(3,2)*w2 &
                      &                         -rrijk(3,3)*w3

                 ko3=ielat(i)
                 if (lTPcel.EQV..true.) then
                    sigc(1,1,ko3)=sigc(1,1,ko3)-rrijk(1,1)*w1c-rrijk(1,2)*w2c &
                         &                           -rrijk(1,3)*w3c
                    sigc(1,2,ko3)=sigc(1,2,ko3)-rrijk(6,1)*w1c-rrijk(6,2)*w2c &
                         &                           -rrijk(6,3)*w3c
                    sigc(1,3,ko3)=sigc(1,3,ko3)-rrijk(5,1)*w1c-rrijk(5,2)*w2c &
                         &                           -rrijk(5,3)*w3c
                    sigc(2,1,ko3)=sigc(2,1,ko3)-rrijk(6,1)*w1c-rrijk(6,2)*w2c &
                         &                           -rrijk(6,3)*w3c
                    sigc(2,2,ko3)=sigc(2,2,ko3)-rrijk(2,1)*w1c-rrijk(2,2)*w2c &
                         &                           -rrijk(2,3)*w3c
                    sigc(2,3,ko3)=sigc(2,3,ko3)-rrijk(4,1)*w1c-rrijk(4,2)*w2c &
                         &                           -rrijk(4,3)*w3c
                    sigc(3,1,ko3)=sigc(3,1,ko3)-rrijk(5,1)*w1c-rrijk(5,2)*w2c &
                         &                           -rrijk(5,3)*w3c
                    sigc(3,2,ko3)=sigc(3,2,ko3)-rrijk(4,1)*w1c-rrijk(4,2)*w2c &
                         &                           -rrijk(4,3)*w3c
                    sigc(3,3,ko3)=sigc(3,3,ko3)-rrijk(3,1)*w1c-rrijk(3,2)*w2c &
                         &                           -rrijk(3,3)*w3c
                 endif
              endif
           endif
           !C --- fin de la boucle sur les k
        end DO
        !C --- fin de la boucle sur les j
     end DO
     !C --- fin de la boucle sur les i

  end DO

if(noxyz==1)  call cryst_to_cart (imm, xp, at, 1)     !cryst vers cart



  potist=potist+potcp
  !  write(6,*)'sortie 3c'
  return
end subroutine calfo3c

end module
