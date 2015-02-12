!----------------------------------------------------------------------
SUBROUTINE calfoeamtabvois(xp, vp,  fp,  ielat, iwmax, ityp)
  !tentative de calfoeam avec une seule grande boucle sur i
  USE T_kind_param_m
  use gen_com_m
  use var_pot
  use jqmod
  implicit none

  !           version du 4 juin 2010, 15h20 - last chaged by MCM (xpnp sa mere)
  ! *****************************************************************
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  ! eam variables


  integer  :: ielat(imm)
  integer  :: iwmax(imm)
  integer  :: ityp(imm)
  real(double)  :: xp(3,imm)
  real(double)  :: vp(3,imm)
  real(double)  :: fp(3,imm)

  !local variables
  integer :: i,j !atomes
  integer ::iti,itj !types
  integer :: l !paires
  integer::  iw1,iw2, iw,k,izero,icccc,ic!vois




  REAL(double), dimension(1:3) :: dxp, aCell, gradij
  real(double) :: r,r2 !distance i-j
  real(double) :: Erep,dErep ! potentiel et gradient de la repulsion de paire ij
  real(double) :: dEembi, Eembi ! potentiel et gradient de l'immersion
  real(double) :: rhoi,rhoj, drhoi, drhoj ! densite de i sur j et j sur i et leurs derivees radiales
  REAL(double) :: Femb, dFemb
  real(double):: fpnemd(3,imm),fpnemdmoy(3), XijdotF
  real(double) :: drk, ktor, inv_ktor, ktorho, inv_ktorho
  real(double), dimension(3) :: fij
  real(double) :: inv_volu, inv_atomic_volu

  real(double) :: densityi,tabdensity(imm)
  LOGICAL :: test_sigma

  real(double)::rue,rue2
  
  real(double), dimension(:,:), allocatable :: xpnp
 
!  write(6,*)'eamtabvois'
  rue=rue_pot(ipotentiel)
!  if (lprteat.EQV..true.) then
!     eatom(:)=0.
!  end if
  ktor=rue/ngrid
  inv_ktor=1.d0/ktor
  ktorho=(rhomax-rhomin)/ngrid
  inv_ktorho = 1.d0/ktorho
  tabdensity(:)=0.
!  fp(:,:) = 0.0
!  sig(:,:)=0.
!  potist = zero
  potisrep=0.
  potisglue=0.
  rue2=rue**2
  test_sigma=(mod(it,itesigma)==0)
  jq=0.
  fpnemd(:,:)=0.
  inv_volu = 1.d0/volu
  inv_atomic_volu = dble(im)/volu

  iw2=0
  ALLOCATE(xpnp(3,imm))
  if (lperiod) then
   xpnp(:,:)=xp(:,:)
  else
   call notperiod(xp,xpnp)
  end if
   
  call cryst_to_cart (imm, xpnp, bg, -1)    !cart vers cryst

  ! ===============================
  loop1at1: do i=1,im
     !       densityi=tabdensity(i)
     !       nvi=0
     ! --- Calcul de la densite sur i ---    
     Eembi=0.0; dEembi=0.0


     iti = ityp(i)
     iw1 = iw2+1
     iw2 = iwmax(i)
     loopvois :do iw = iw1, iw2
        j = indi(iw)

        dxp(1:3) = xpnp(1:3,i) - xpnp(1:3,j)
        WHERE ( (dxp.GT.0.5d0).OR.(dxp.LT.-0.5d0) )
           dxp(1:3) = dxp(1:3) - Dble(Nint(dxp(1:3)))
        END WHERE
        ! Transformation des coordonnees reduites en cartesiennes
        dxp = MatMul(at,dxp)

        ! Calcul du carre de la distance
        do izero=1,3
           if (dabs(dxp(izero)).lt.low_limit) then
              dxp(izero) = zero
           end if
        end do
        r2 = Sum( dxp(1:3)**2 )

        if (r2>rue2) cycle

        itj=ityp(j)
        r=sqrt(r2)             
        if (r.eq.zero) & 
             write(*,*) '1. WARNING IN calfoeamtabvois TWO ATOMS VERY CLOSE i ,j , dist(angst)', i ,j , r*angst
        k=Int(r*inv_ktor)
        drk=r-k*ktor

        !densites sur i et j
        rhoj = eamrho(1,itj,k) + drk*( eamrho(2,itj,k) + drk*( eamrho(3,itj,k) + drk*eamrho(4,itj,k) ) )  !rho de j sur i
        tabdensity(i) = tabdensity(i) + rhoj
        IF (ldemitab) THEN
                rhoi = eamrho(1,iti,k) + drk*( eamrho(2,iti,k) + drk*( eamrho(3,iti,k) + drk*eamrho(4,iti,k) ) )  !rho de i sur j
                tabdensity(j)=tabdensity(j)+rhoi
        END IF

     end do loopvois
     !       tabdensity(i)=densityi
  end do loop1at1

  ! ===============================
  iw2=0
  ! calcul et stockage de Eembi et dEembi
  loop2at1: do i=1,im
     iti=ityp(i)
     k=Int((tabdensity(i)-rhomin)*inv_ktorho)
     if(k.gt.ngrid) then
        write(6,*)k, ngrid, 'k> ngrid ; augmenter le facteur multiplicatif de rhomax dans calpo'
        write(6,*)'densityi',k,ngrid,tabdensity(i), rhomin, inv_ktorho, densityi
        stop
     end if
     drk=tabdensity(i)-(rhomin+k*ktorho)
     Eembi = eamglue(1,iti,k) + drk*( eamglue(2,iti,k) + drk*( eamglue(3,iti,k) + drk*eamglue(4,iti,k) ) )
     if( associated (free)) then
        if( ( (lprteat.EQV..true.).or.(lcalcjq==.true.) ).and.( free(i).EQV..true.)) eatom(i)=eatom(i)+Eembi
        if( free(i).EQV..true.)potisglue = potisglue+Eembi
     else
        if((lprteat.EQV..true.).or.(lcalcjq==.true.)) eatom(i)=eatom(i)+Eembi
        potisglue = potisglue+Eembi
     end if
    if (lforcetabulate) then
     tabdensity(i)= eamglue_d(1,iti,k) + drk*( eamglue_d(2,iti,k) + drk*( eamglue_d(3,iti,k) + drk*eamglue_d(4,iti,k) ) )
    else 
     tabdensity(i) = eamglue(2,iti,k) + drk*( 2.0*eamglue(3,iti,k) + 3.0*drk*eamglue(4,iti,k) )
    end if 
 end do loop2at1


  ! ===============================
  !boucle des forces
  loop3at1: do i=1,im
     iti = ityp(i)
     iw1 = iw2+1
     iw2 = iwmax(i)
     loopvois2 :do iw = iw1, iw2
        j = indi(iw)

        dxp(1:3) = xpnp(1:3,i) - xpnp(1:3,j)
        WHERE ( (dxp(:).GT.0.5d0).OR.(dxp(:).LT.-0.5d0) )
           dxp(:) = dxp(:) - Dble(Nint(dxp(:)))
        END WHERE

        ! Transformation des coordonnees reduites en cartesiennes
        dxp = MatMul(at,dxp)

        do izero=1,3
          if (dxp(izero).eq.zero) cycle
           if (dabs(dxp(izero)).lt.low_limit) then
               write(*,*) 'WARNING low_limit'
               dxp(izero) = zero
           end if
        end do
        r2 = Sum( dxp(1:3)**2 )

        if (r2>rue2) cycle

        itj=ityp(j)
        r=sqrt(r2)      

        if (r.eq.zero) & 
             write(*,*) '2. WARNING IN calfoeamtabvois TWO ATOMS VERY CLOSE i ,j , dist(angst)', i ,j , r*angst

        gradij(1:3) = dxp(1:3)/r
        k=Int(r*inv_ktor)
        drk=r-k*ktor

        !terme de repulsion 
        l = ipo(iti,itj)
        Erep = eamrep(1,l,k) + drk*( eamrep(2,l,k) + drk*( eamrep(3,l,k) + drk*eamrep(4,l,k) ) )
        if (lforcetabulate) then
          dErep = eamrep_d(1,l,k) + drk*( eamrep_d(2,l,k) + drk*( eamrep_d(3,l,k) + drk*eamrep_d(4,l,k) ) )
          drhoj = eamrho_d(1,itj,k) + drk*( eamrho_d(2,itj,k) + drk*( eamrho_d(3,itj,k) + drk*eamrho_d(4,itj,k) ) )  !rho_d de j sur i
          drhoi = eamrho_d(1,iti,k) + drk*( eamrho_d(2,iti,k) + drk*( eamrho_d(3,iti,k) + drk*eamrho_d(4,itj,k) ) )  !rho_d de j sur i
        else
          dErep = eamrep(2,l,k) + drk*( 2.0*eamrep(3,l,k) + 3.0*drk*eamrep(4,l,k) )
          drhoj = eamrho(2,itj,k) + drk*( 2.0*eamrho(3,itj,k) + 3.0*drk*eamrho(4,itj,k) )
          drhoi = eamrho(2,iti,k) + drk*( 2.0*eamrho(3,iti,k) + 3.0*drk*eamrho(4,iti,k) )
        end if
        !dFemb = ( tabdensity(i) + tabdensity(j))*drhoj
        dFemb = tabdensity(i)*drhoj + tabdensity(j)*drhoi  ! THIS is WRONG in my SENSE

        if((lprteat.EQV..true.).or.(lcalcjq==.true.))then
           if (associated (free)) then
              if( free(i).EQV..true.)                eatom(i)=eatom(i) + 0.5d0*Erep
              if ((free(j).EQV..true.).and.ldemitab) eatom(j)=eatom(j) + 0.5d0*Erep
           else
              eatom(i)=eatom(i) + 0.5d0*Erep
              if (ldemitab)  eatom(j)=eatom(j) + 0.5d0*Erep
           end if
        end if
        if (associated (free)) then
           if( free(i).EQV..true.)                potisrep = potisrep + 0.5*Erep
           if(( free(j).EQV..true.).and.ldemitab) potisrep = potisrep + 0.5*Erep
        else
                IF (ldemitab) THEN
                        potisrep = potisrep + Erep
                ELSE
                        potisrep = potisrep + 0.5d0*Erep
                END IF
        end if
        
        fij(:) = - (dFemb+dErep)*gradij(1:3)
        fp(1:3,i) = fp(1:3,i) + fij(:)
        if (ldemitab) fp(1:3,j) = fp(1:3,j) - fij(:)

        if (lnemd) then
           XijdotF=dxp(1)*Fnemd
           do ic=1,3
              fpnemd(ic,i)=fpnemd(ic,i) +0.5*(dFemb+dErep)*gradij(ic)*XijdotF
              fpnemd(ic,j)=fpnemd(ic,j) +0.5*(dFemb+dErep)*gradij(ic)*XijdotF
           end do
        end if
        if(lcalcjq) then
           jqf=0.0
           do ic=1,3
              jqf = jqf + ( dFemb*vp(ic,j) + 0.5*dErep*(vp(ic,i)+vp(ic,j)) )*gradij(ic)
           end do
           do ic=1,3
              jq(ic) = jq(ic) - jqf*gradij(ic)*r
           end do
        end if

        if (test_sigma) then                   
                IF (ldemitab) THEN
                        sig(1:3,1) = sig(1:3,1) + inv_volu*fij(1:3)*dxp(1)
                        sig(1:3,2) = sig(1:3,2) + inv_volu*fij(1:3)*dxp(2)
                        sig(1:3,3) = sig(1:3,3) + inv_volu*fij(1:3)*dxp(3)
                        IF (lPrtSigat) THEN
                                sigat(1:3,1,i) = sigat(1:3,1,i) + 0.5d0*inv_atomic_volu*fij(1:3)*dxp(1)
                                sigat(1:3,2,i) = sigat(1:3,2,i) + 0.5d0*inv_atomic_volu*fij(1:3)*dxp(2)
                                sigat(1:3,3,i) = sigat(1:3,3,i) + 0.5d0*inv_atomic_volu*fij(1:3)*dxp(3)
                                sigat(1:3,1,j) = sigat(1:3,1,j) + 0.5d0*inv_atomic_volu*fij(1:3)*dxp(1)
                                sigat(1:3,2,j) = sigat(1:3,2,j) + 0.5d0*inv_atomic_volu*fij(1:3)*dxp(2)
                                sigat(1:3,3,j) = sigat(1:3,3,j) + 0.5d0*inv_atomic_volu*fij(1:3)*dxp(3)
                        END IF
                ELSE
                        sig(1:3,1) = sig(1:3,1) - 0.5d0*inv_volu*fij(1:3)*dxp(1)
                        sig(1:3,2) = sig(1:3,2) - 0.5d0*inv_volu*fij(1:3)*dxp(2)
                        sig(1:3,3) = sig(1:3,3) - 0.5d0*inv_volu*fij(1:3)*dxp(3)
                        IF (lPrtSigat) THEN
                                sigat(1:3,1,i) = sigat(1:3,1,i) + 0.5d0*inv_atomic_volu*fij(1:3)*dxp(1)
                                sigat(1:3,2,i) = sigat(1:3,2,i) + 0.5d0*inv_atomic_volu*fij(1:3)*dxp(2)
                                sigat(1:3,3,i) = sigat(1:3,3,i) + 0.5d0*inv_atomic_volu*fij(1:3)*dxp(3)
                        END IF
                END IF
        end if
     end do loopvois2


  end do loop3at1

  ! Energie potentielle totale
  potist = potisglue + potisrep

  if (lnemd) then
     fpnemdmoy=0
     do i=1,imd   
        do l=1,3
           fpnemdmoy(l)=fpnemdmoy(l)+fpnemd(l,i)/float(imd)
        enddo
     end do

     do i=1,imd
!        write(6,*)'A',i,fp(:,i)
        do l=1,3
           fp(l,i)=fp(l,i)-fpnemdmoy(l)
           fp(l,i)=fp(l,i)+fpnemd(l,i)
        enddo
!        write(6,*)'B',i,fp(:,i)
     end do
  end if


  !call cryst_to_cart (imm, xp, at, 1)     !cryst vers cart
   DEALLOCATE (xpnp)
!  write(6,*)'eamtabvois'
  return
end SUBROUTINE calfoeamtabvois
