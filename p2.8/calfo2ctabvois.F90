module calfo2ctabvois_mod
        USE cryst_to_cart_mod,only: cryst_to_cart
        USE calerf_mod,only: calerf
        USE potrep_mod,only: potrep
        USE calfocommon
        implicit none
        contains
! **********************************************************
subroutine calfo2ctabvois(im,imm,xp,  vp, fp,  iwmax, ityp,indi )
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  USE gen_com_m, ONLY:at,bg,&
       &lcalcjq,pi,potis1,potis2,volu

  USE var_pot, ONLY:alpha,csive,ipo,zz,rue_pair,ipo,pot
  USE jqmod
  ! **********************************************************

  implicit none
  !-----------------------------------------------
  !   G l o b a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  integer,intent(in)::im,imm
  integer , intent(in),allocatable :: iwmax(:),ityp(:),indi(:)
  real(double),intent(in),allocatable  :: vp(:,:)
  real(double),intent(inout),allocatable  :: xp(:,:)
  real(double) , intent(inout),allocatable :: fp(:,:)
  !-----------------------------------------------
  !   L o c a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: iw2, iti, l, iw1, i, j, itj, k, &
       iw,  ic

  real(double) :: aux, alp, a1, a2, a3, f1, f2, f3,  &
       c1, c2, c3, ddsq,  sk, phu, ra(3),eatcomp
  real(double) ::  dr, r,r2,partsig,deltaepot
  real(double), dimension(1,3) :: cv
  REAL(double), dimension(1:3) :: dxp, aCell, gradij

  aux = 23.06134575D-20
  alp = alpha/sqrt(pi)*aux
  iw2 = 0
  do i = 1, im

     iti = ityp(i)
     l = ipo(iti,iti)
     ! --- Calcul du second potentiel de la somme d'Ewald ---
     potis2 = potis2-zz(l)*alp
  end do

  ! --------------------------
  !   OUVERTURE BOUCLE SUR I
  ! --------------------------
  call cryst_to_cart (imm, xp, bg, -1)    !cart vers cryst
  do i = 1, im-1

     iti = ityp(i)
     iw1 = iw2+1
     iw2 = iwmax(i)
     do iw = iw1, iw2
        j = indi(iw)
        itj=ityp(j);l=ipo(iti,itj)
        dxp(1:3) = xp(1:3,i) - xp(1:3,j)

        WHERE ( (dxp.GT.0.5d0).OR.(dxp.LT.-0.5d0) )
           dxp(1:3) = dxp(1:3) - Dble(Nint(dxp(1:3)))
        END WHERE

        ! Application des conditions aux limites pÃ©riodiques
        !           IF (dxp(1)>0.5d0) THEN
        !              dxp(1) = dxp(1)-1.d0
        !           ELSE IF (dxp(1)<-0.5d0) THEN
        !              dxp(1) = dxp(1)+1.d0
        !           END IF
        !           IF (dxp(2)>0.5d0) THEN
        !              dxp(2) = dxp(2)-1.d0
        !           ELSE IF (dxp(2)<-0.5d0) THEN
        !              dxp(2) = dxp(2)+1.d0
        !           END IF
        !           IF (dxp(3)>0.5d0) THEN
        !              dxp(3) = dxp(3)-1.d0
        !           ELSE IF (dxp(3)<-0.5d0) THEN
        !              dxp(3) = dxp(3)+1.d0
        !           END IF

        dxp = MatMul(at,dxp)



        r2 = Sum( dxp(1:3)**2 )

        if (r2>rue_pair(l)**2) cycle
        r=sqrt(r2)

        itj = ityp(j)
        l = ipo(iti,itj)
        sk = r/csive
        k = sk
        ! spline
        dr = r-float(k)*csive
        
        deltaepot=0.5*(pot(1,l,k)+r*dr*(pot(2,l,k)+dr*(pot(3,l,k)+dr*pot(4,l,k))))
        phu = -1.0*(pot(2,l,k)+dr*(2.0*pot(3,l,k)+dr*(3.0*pot(4,l,k))))

!        write(6,*)'i,j,r,sk,k'
!        write(6,*)i,j,r,sk,k
!        write(6,*)pot(1,l,k),pot(2,l,k),pot(3,l,k),pot(4,l,k)

!        if (allocated(free)) then
!           if (free(i).EQV..true.)potis1 = potis1+deltaepot
!           if (free(j).EQV..true.)potis1 = potis1+deltaepot
!        else
        potis1 = potis1+2*deltaepot
        
!        end if

        ! --- Fin du calcul ---
        fp(:,i)=fp(:,i)+phu*dxp(:)
        fp(:,j)=fp(:,j)-phu*dxp(:)

        !        ra(1) = phu*cv(1,1)
        !        f1 = f1+ra(1)
        !        fp(1,j) = fp(1,j)-ra(1)
        !        ra(2) = phu*cv(1,2)
        !        f2 = f2+ra(2)
        !        fp(2,j) = fp(2,j)-ra(2)
        !        ra(3) = phu*cv(1,3)
        !        f3 = f3+ra(3)
        !        fp(3,j) = fp(3,j)-ra(3)
        !
        ! A commenter qd lcalcjq=false pour ne pas perdre de temps dans le test
        !ra(3)=Force de j sur i
           if (lprteat) then
!              if (allocated(free)) then
!                 if (free(i).EQV..true.)eat(i) = eat(i)+deltaepot
!                 if (free(j).EQV..true.)eat(j) = eat(j)+deltaepot
!              else
                 eat(i) = eat(i)+deltaepot
                 eat(j) = eat(j)+deltaepot
!              end if
           end if
        if (lcalcjq) then
           jqf=0.0
           eat(i) = eat(i)+deltaepot
           eat(j) = eat(j)+deltaepot
           
           do ic=1,3
              jqf=jqf-0.5*(ra(ic)*(vp(ic,i)+vp(ic,j)))
           end do
           do ic=1,3
              jq(ic)=jq(ic)-jqf*dxp(ic)
           end do
        end if

        if (test_sigma) then
           partsig=phu/volu
           ! calcul de sigma contrainte
           sig(1,1) = sig(1,1)+partsig*dxp(1)*dxp(1)
           sig(1,2) = sig(1,2)+partsig*dxp(1)*dxp(2)
           sig(1,3) = sig(1,3)+partsig*dxp(1)*dxp(3)
           sig(2,1) = sig(2,1)+partsig*dxp(2)*dxp(1)
           sig(2,2) = sig(2,2)+partsig*dxp(2)*dxp(2)
           sig(2,3) = sig(2,3)+partsig*dxp(2)*dxp(3)
           sig(3,1) = sig(3,1)+partsig*dxp(3)*dxp(1)
           sig(3,2) = sig(3,2)+partsig*dxp(3)*dxp(2)
           sig(3,3) = sig(3,3)+partsig*dxp(3)*dxp(3)
           if(lsigat)then
              sigat(1,1,i) = sigat(1,1,i)+0.5*partsig*dxp(1)*dxp(1)
              sigat(1,2,i) = sigat(1,2,i)+0.5*partsig*dxp(1)*dxp(2)
              sigat(1,3,i) = sigat(1,3,i)+0.5*partsig*dxp(1)*dxp(3)
              sigat(2,1,i) = sigat(2,1,i)+0.5*partsig*dxp(2)*dxp(1)
              sigat(2,2,i) = sigat(2,2,i)+0.5*partsig*dxp(2)*dxp(2)
              sigat(2,3,i) = sigat(2,3,i)+0.5*partsig*dxp(2)*dxp(3)
              sigat(3,1,i) = sigat(3,1,i)+0.5*partsig*dxp(3)*dxp(1)
              sigat(3,2,i) = sigat(3,2,i)+0.5*partsig*dxp(3)*dxp(2)
              sigat(3,3,i) = sigat(3,3,i)+0.5*partsig*dxp(3)*dxp(3)
              sigat(1,1,j) = sigat(1,1,j)+0.5*partsig*dxp(1)*dxp(1)
              sigat(1,2,j) = sigat(1,2,j)+0.5*partsig*dxp(1)*dxp(2)
              sigat(1,3,j) = sigat(1,3,j)+0.5*partsig*dxp(1)*dxp(3)
              sigat(2,1,j) = sigat(2,1,j)+0.5*partsig*dxp(2)*dxp(1)
              sigat(2,2,j) = sigat(2,2,j)+0.5*partsig*dxp(2)*dxp(2)
              sigat(2,3,j) = sigat(2,3,j)+0.5*partsig*dxp(2)*dxp(3)
              sigat(3,1,j) = sigat(3,1,j)+0.5*partsig*dxp(3)*dxp(1)
              sigat(3,2,j) = sigat(3,2,j)+0.5*partsig*dxp(3)*dxp(2)
              sigat(3,3,j) = sigat(3,3,j)+0.5*partsig*dxp(3)*dxp(3)
           end if
        end if
     end do

  end do

!        do i=1,im,100
!           write(6,*)i,fp(1,i),fp(2,i),fp(3,i)
!        end do


  call cryst_to_cart (imm, xp, at, 1)     !cryst vers cart

  return
end subroutine calfo2ctabvois
end module
