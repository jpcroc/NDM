module calfoeamtabvois_mod
  USE arret_ndm_mod,only:arret_ndm
  USE gen_com_m, only:uwrt,lwrt,angst,fnemd,lcalcjq,ldemitab,&
       &lnemd,low_limit,lperiod,zero,pi
  USE calfocommon
  USE atomconfig,only : atom_config,atom_config_d,atom_config_e
  USE cellconfig, only : cell_config
  use boxconfig,only: box_config
  use vect_dist_mod,only:vect_dist

  implicit none
contains
  !----------------------------------------------------------------------
  SUBROUTINE calfoeamtabvois(atcf,celcf,boxcf)! 
    !tentative de calfoeam avec une seule grande boucle sur i
    USE T_kind_param_m
    USE var_pot, ONLY:ipotentiel,lforcetabulate,ngrid,potisglue,potisrep,rhomax,rhomin,eamrho,eamglue_d,&
         &eamglue,ipo,eamrep,eamrep_d,eamrho_d,rue_pot,alpha,ntyp

    USE jqmod
    implicit none
  class(atom_config),intent(inout)::atcf
  type(cell_config),intent(in)::celcf
  type(box_config),intent(in)::boxcf



  !local variables
    integer :: i,j !atomes
    integer ::iti,itj !types
    integer :: l !paires
    integer::  iw1,iw2, iw,k,ic!vois
  REAL(double), dimension(1:3) :: dxp, gradij
    real(double) :: r !distance i-j
    real(double) :: Erep,dErep ! potentiel et gradient de la repulsion de paire ij
    real(double) :: dEembi, Eembi ! potentiel et gradient de l'immersion
    real(double) :: rhoi,rhoj, drhoi, drhoj ! densite de i sur j et j sur i et leurs derivees radiales
    REAL(double) ::  dFemb
    real(double):: fpnemd(3,atcf%im),fpnemdmoy(3), XijdotF
    real(double) :: drk, ktor, inv_ktor
    real(double),dimension(:),allocatable::ktorho(:), inv_ktorho(:)
    real(double), dimension(3) :: fij
    real(double) :: inv_volu, inv_atomic_volu
    real(double) :: densityi,tabdensity(atcf%imm)
    real(double)::rue

    real(double)::aux,alp
    logical::linter
    aux = 23.06134575D-20
    alp = alpha/sqrt(pi)*aux
    allocate(ktorho(ntyp))
    allocate(inv_ktorho(ntyp))

    !  write(6,*)'eamtabvois'
    rue=rue_pot(ipotentiel)
    !  if (lprteat.EQV..true.) then
    !     eat(:)=0.
    !  end if
    ktor=rue/ngrid
    inv_ktor=1.d0/ktor
    ktorho(:)=(rhomax(:)-rhomin(:))/ngrid
    inv_ktorho(:) = 1.d0/ktorho(:)
    tabdensity(:)=0.
    !  fp(:,:) = 0.0
    !  sig(:,:)=0.
    !  potist = zero
    potisrep=0.
    potisglue=0.
!    rue2=rue**2
    jq=0.
    fpnemd(:,:)=0.
    inv_volu = 1.d0/boxcf%volu
    inv_atomic_volu = dble(atcf%im)/boxcf%volu

    iw2=0



    ! ===============================
    loop1at1: do i=1,atcf%im
       !       densityi=tabdensity(i)
       !       nvi=0
       ! --- Calcul de la densite sur i ---    
       Eembi=0.0; dEembi=0.0


       iti = atcf%ityp(i)
             
       iw1 = iw2+1
       iw2 = atcf%iwmax(i)
       loopvois :do iw = iw1, iw2
          j = atcf%indi(iw)
          itj=atcf%ityp(j)
          

          call vect_dist(atcf,celcf,boxcf,i,j,VJI=dxp, lperiod=boxcf%lperiod,rum=rue&
               &,linter=linter,dist=r)
          if (.not.linter) cycle

          if (r.eq.zero) & 
               write(*,*) '1. WARNING IN calfoeamtabvois TWO ATOMS VERY CLOSE i ,j , dist(angst)', i ,j , r*angst
!          k=Int(r*inv_ktor)
          k=min(ngrid,Int(r*inv_ktor))
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
    loop2at1: do i=1,atcf%im
       iti=atcf%ityp(i)
!       k=Int((tabdensity(i)-rhomin(iti))*inv_ktorho(iti))
       k=min(ngrid,Int((tabdensity(i)-rhomin(iti))*inv_ktorho(iti)))
       if(k.gt.ngrid) then
          write(uwrt,*)k, ngrid, 'k> ngrid ; augmenter le facteur multiplicatif de rhomax dans calpo'
          write(uwrt,*)'densityi',k,ngrid,tabdensity(i), rhomin(iti), inv_ktorho(iti), densityi
          call arret_ndm(.true.)
       end if
       drk=tabdensity(i)-(rhomin(iti)+k*ktorho(iti))
       Eembi = eamglue(1,iti,k) + drk*( eamglue(2,iti,k) + drk*( eamglue(3,iti,k) + drk*eamglue(4,iti,k) ) )
       if(lprteat.EQV..true.) then
          select type (atcf)
          class is (atom_config_e)
             atcf%eat(i)=atcf%eat(i)+Eembi
          end select
       end if
          potisglue = potisglue+Eembi
       if (lforcetabulate) then
          tabdensity(i)= eamglue_d(1,iti,k) + drk*( eamglue_d(2,iti,k) + drk*( eamglue_d(3,iti,k) + drk*eamglue_d(4,iti,k) ) )
       else 
          tabdensity(i) = eamglue(2,iti,k) + drk*( 2.0*eamglue(3,iti,k) + 3.0*drk*eamglue(4,iti,k) )
       end if
    end do loop2at1


    ! ===============================
    !boucle des forces
    loop3at1: do i=1,atcf%im
       iti = atcf%ityp(i)
       iw1 = iw2+1
       iw2 = atcf%iwmax(i)
       loopvois2 :do iw = iw1, iw2
          j = atcf%indi(iw)
          itj=atcf%ityp(j)
          call vect_dist(atcf,celcf,boxcf,i,j,VJI=dxp, lperiod=boxcf%lperiod,rum=rue&
               &,linter=linter,dist=r)
          if (.not.linter) cycle
!          r=sqrt(r2)      
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

          if((lprteat.EQV..true.).or.(lcalcjq.EQV..true.))then
             select type (atcf)
             class is (atom_config_e)
                atcf%eat(i)=atcf%eat(i) + 0.5d0*Erep
                if (ldemitab)  then
                   atcf%eat(j)=atcf%eat(j) + 0.5d0*Erep
!             end if
                end if
             end select
          end if

             IF (ldemitab) THEN
                potisrep = potisrep + Erep
             ELSE
                potisrep = potisrep + 0.5d0*Erep
             END IF

          fij(:) = - (dFemb+dErep)*gradij(1:3)
          atcf%fp(1:3,i) = atcf%fp(1:3,i) + fij(:)
          if (ldemitab) atcf%fp(1:3,j) = atcf%fp(1:3,j) - fij(:)

          if (lnemd) then
             XijdotF=dxp(1)*Fnemd
             do ic=1,3
                fpnemd(ic,i)=fpnemd(ic,i) +0.5*(dFemb+dErep)*gradij(ic)*XijdotF
                fpnemd(ic,j)=fpnemd(ic,j) +0.5*(dFemb+dErep)*gradij(ic)*XijdotF
             end do
          end if
!!$          if(lcalcjq) then
!!$             jqf=0.0
!!$             do ic=1,3
!!$                jqf = jqf + ( dFemb*vp(ic,j) + 0.5*dErep*(vp(ic,i)+vp(ic,j)) )*gradij(ic)
!!$             end do
!!$             do ic=1,3
!!$                jq(ic) = jq(ic) - jqf*gradij(ic)*r
!!$             end do
!!$          end if
!!$
          if (test_sigma) then                   
             IF (ldemitab) THEN
                sigcalfo(1:3,1) = sigcalfo(1:3,1) + inv_volu*fij(1:3)*dxp(1)
                sigcalfo(1:3,2) = sigcalfo(1:3,2) + inv_volu*fij(1:3)*dxp(2)
                sigcalfo(1:3,3) = sigcalfo(1:3,3) + inv_volu*fij(1:3)*dxp(3)
                IF (lSigat) THEN
                   select type (atcf)
                   class is (atom_config_e)
                   atcf%sigat(1:3,1,i) = atcf%sigat(1:3,1,i) + 0.5d0*inv_atomic_volu*fij(1:3)*dxp(1)
                   atcf%sigat(1:3,2,i) = atcf%sigat(1:3,2,i) + 0.5d0*inv_atomic_volu*fij(1:3)*dxp(2)
                   atcf%sigat(1:3,3,i) = atcf%sigat(1:3,3,i) + 0.5d0*inv_atomic_volu*fij(1:3)*dxp(3)
                   atcf%sigat(1:3,1,j) = atcf%sigat(1:3,1,j) + 0.5d0*inv_atomic_volu*fij(1:3)*dxp(1)
                   atcf%sigat(1:3,2,j) = atcf%sigat(1:3,2,j) + 0.5d0*inv_atomic_volu*fij(1:3)*dxp(2)
                   atcf%sigat(1:3,3,j) = atcf%sigat(1:3,3,j) + 0.5d0*inv_atomic_volu*fij(1:3)*dxp(3)
                end select
                END IF
             ELSE
                sigcalfo(1:3,1) = sigcalfo(1:3,1) + 0.5d0*inv_volu*fij(1:3)*dxp(1)
                sigcalfo(1:3,2) = sigcalfo(1:3,2) + 0.5d0*inv_volu*fij(1:3)*dxp(2)
                sigcalfo(1:3,3) = sigcalfo(1:3,3) + 0.5d0*inv_volu*fij(1:3)*dxp(3)
                IF (lSigat) THEN
                   select type (atcf)
                   class is (atom_config_e)
                   atcf%sigat(1:3,1,i) = atcf%sigat(1:3,1,i) + 0.5d0*inv_atomic_volu*fij(1:3)*dxp(1)
                   atcf%sigat(1:3,2,i) = atcf%sigat(1:3,2,i) + 0.5d0*inv_atomic_volu*fij(1:3)*dxp(2)
                   atcf%sigat(1:3,3,i) = atcf%sigat(1:3,3,i) + 0.5d0*inv_atomic_volu*fij(1:3)*dxp(3)
                end select
                END IF
             END IF
          end if
       end do loopvois2


    end do loop3at1

    ! Energie potentielle totale
    potistcalfo = potisglue + potisrep

    if (lnemd) then
       fpnemdmoy=0
       do i=1,atcf%im   
          do l=1,3
             fpnemdmoy(l)=fpnemdmoy(l)+fpnemd(l,i)/float(atcf%im)
          enddo
       end do

       do i=1,atcf%im
          !        write(6,*)'A',i,fp(:,i)
          do l=1,3
             atcf%fp(l,i)=atcf%fp(l,i)-fpnemdmoy(l)
             atcf%fp(l,i)=atcf%fp(l,i)+fpnemd(l,i)
          enddo
          !        write(6,*)'B',i,fp(:,i)
       end do
    end if


    !  write(6,*)'eamtabvois'
    return
  end SUBROUTINE calfoeamtabvois
end module calfoeamtabvois_mod
