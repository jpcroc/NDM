module calfoeamcel_mod
  USE arret_ndm_mod,only:arret_ndm
  USE gen_com_m, ONLY:angst,nvat,low_limit,lperiod,zero,pi
  USE calfocommon
  use vect_dist_mod,only:vect_dist
  USE atomconfig,only : atom_config,atom_config_d,atom_config_e
  USE cellconfig, only : cell_config
  use boxconfig,only: box_config

  implicit none
contains
  !----------------------------------------------------------------------
  SUBROUTINE calfoeamcel(atcf,celcf,boxcf,psc)
    USE T_kind_param_m

    USE var_pot, ONLY:ipotentiel,ngrid,potiseam,potisglue,potisrep,rue_pot,&
         &typ_and_pot,typ_pot_pair,ipotentiel,ngrid,potiseam,potisglue,potisrep,rhomax,rhomin,eamrho,ipo,eamrep,&
         &eamglue,alpha,zz,ntyp

#ifdef PARA
    use Tpara,only:nprocspace,para_space_config,ierr,comm_space
    USE mod_para,only:maj_tabdensity_ftm


#else
    USE Tpara,only:nprocspace,para_space_config
#endif

    class(atom_config),intent(inout)::atcf
    type(cell_config),intent(in)::celcf
    type(box_config),intent(in)::boxcf

    type(para_space_config)::psc
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    ! eam variables
    !local variables
    integer :: i,j !atomes
    integer ::iti,itj !types
    integer :: l !paires
    integer:: koo,ko1,i2,i1 !cel.
    integer :: k ! aux pour splines

!    real(double) :: rue2 !coupure**2
    REAL(double), dimension(1:3) :: cp, dxp, gradij
    real(double) :: r!distance i-j
    real(double) :: Fij,cv(1,3)
    real(double) :: Eembi,dEembi ! potentiel et gradient de l'immersion
    real(double) :: rhoi,rhoj ! densite de i sur j et j sur i
    real(double) :: Erep,dErep ! potentiel et gradient de la repulsion de paire ij
    REAL(double) :: Femb
    real(double) :: rk, drk,ktor, inv_ktor
    real(double),dimension(:),allocatable::ktorho(:), inv_ktorho(:)
    real(double) :: densityi !densite totale sur i
    integer :: izero
    real(double) :: tabdensity(atcf%imm)
    real(double)::rue,alp,aux
    logical ::linter
    
    rue=rue_pot(ipotentiel)
    aux = 23.06134575D-20
    alp = alpha/sqrt(pi)*aux
    allocate(ktorho(ntyp))
    allocate(inv_ktorho(ntyp))
    ktor=rue/ngrid
    inv_ktor=1.d0/ktor
    ktorho(:)=(rhomax(:)-rhomin(:))/ngrid
    inv_ktorho(:) = 1.d0/ktorho(:)
    tabdensity(:)=0.
    potisrep=0.
    potisglue=0.
!    rue2=rue**2



    loop1at1: do i=1,atcf%im
       if (typ_and_pot(atcf%ityp(i),ipotentiel).eqv..false.)cycle

       ! --- Calcul de la densite sur i ---    

       !     nvi=0
       densityi=0.0 ; dEembi=0.0
       koo = atcf%ielat(i)                          ! Numero de la cellule
       iti = atcf%ityp(i)
!!$       if (ipotentiel==16) then
!!$          l=ipo(iti,iti)
!!$          potis2=potis2-zz(l)*alp
!!$       end if

       ! pour chaque cel. voisine
       loop1cel:   do i1 = 0, celcf%ncelvois(koo)
          ko1 = celcf%ncel(koo,i1)
          ! pour chaque atome ds la cel. voisine
          loop1at2: do i2 = 1, celcf%nato(ko1)
             j = celcf%atincel(i2,ko1)
!             write(6,*)i,koo,i1, celcf%ncelvois(koo),i2,j
             if (typ_pot_pair(ipo(atcf%ityp(i),atcf%ityp(j))).ne.ipotentiel) cycle

             itj=atcf%ityp(j)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             !CRC             if(i.eq.j) cycle
#ifdef PARA

             ! Methode pour ne prendre qu'une seule fois en compte
             ! le couple i,j en paralle :
             ! - i est necesairement local (boucle i<=im)
             ! - si j est local on ne retient que le couple i<j
             ! - si j n'est pas local, le couple n'est par definition
             !   pris qu'une fois puisque i est local
             if (nprocspace.gt.1) then

                if (j.le.atcf%im) then

                   ! les deux atomes sont locaux
                   if (atcf%num_at_glob(i).ge.atcf%num_at_glob(j)) cycle !terme deja calcule
                else
                   ! j n'est pas local, on ne fait le calcul normal           
                endif
             else
                if (atcf%num_at_glob(i).ge.atcf%num_at_glob(j)) cycle !terme dÃ£Â©ja calculÃ£Â©
             end if
                
#else
             if (atcf%num_at_glob(i).ge.atcf%num_at_glob(j)) cycle !terme dÃ£Â©ja calculÃ£Â©
#endif

           call vect_dist(atcf,celcf,boxcf,i,j,VJI=dxp,indcv=i1, lperiod=boxcf%lperiod,rum=rue,linter=linter,dist=r)
           if(.not.linter) cycle
             k=Int(r*inv_ktor)
             gradij(1:3) = dxp(1:3)/r
             drk=r-k*ktor
             !             rhoj = eamrho(1,itj,k) + drk*( eamrho(2,itj,k) + drk*( eamrho(3,itj,k) + drk*eamrho(4,itj,k) ) )  !rho de j sur i
             !             densityi=densityi+rhoj

             rhoj = eamrho(1,itj,k) + drk*( eamrho(2,itj,k) + drk*( eamrho(3,itj,k) + drk*eamrho(4,itj,k) ) )  !rho de j sur i
             tabdensity(i)=tabdensity(i)+rhoj
             rhoi = eamrho(1,iti,k) + drk*( eamrho(2,iti,k) + drk*( eamrho(3,iti,k) + drk*eamrho(4,iti,k) ) )  !rho de i sur j
             tabdensity(j)=tabdensity(j)+rhoi
!!$             write(100,'(2I3,3G17.8)')i,j,r,tabdensity(i),tabdensity(j)
!!$             write(100,'(2I3,4G17.8)')i,j,eamrho(1,iti,k) , eamrho(2,iti,k),eamrho(3,iti,k),eamrho(4,iti,k)
!!$             write(100,'(2I3,4G17.8)')i,j,eamrho(1,itj,k) , eamrho(2,itj,k),eamrho(3,itj,k),eamrho(4,itj,k)
             !           nvi=nvi+1
             !           dxpij(1:3,nvi)=dxp(1:3)
             !           jvi(nvi)=j
             !           rij(nvi)=r            

             !terme de repulsion 
             l = ipo(iti,itj)
             Erep = eamrep(1,l,k) + drk*( eamrep(2,l,k) + drk*( eamrep(3,l,k) + drk*eamrep(4,l,k) ) )
             if(lprteat.EQV..true.)then
                select type (atcf)
                class is (atom_config_e)
                   atcf%eat(i)=atcf%eat(i)+Erep/2.d0
                   if (j.le.atcf%im) atcf%eat(j)=atcf%eat(j)+Erep/2.d0
                end select
             end if
             dErep = eamrep(2,l,k) + drk*( 2.0*eamrep(3,l,k) + 3.0*drk*eamrep(4,l,k) )

             if (atcf%num_at_glob(i).lt.atcf%num_at_glob(j)) then
                potisrep = potisrep+Erep
             endif

             atcf%fp(1:3,i)=atcf%fp(1:3,i)-dErep*gradij(1:3)
             atcf%fp(1:3,j)=atcf%fp(1:3,j)+dErep*gradij(1:3)

             if (test_sigma) then        
                if (atcf%num_at_glob(i).lt.atcf%num_at_glob(j)) then          
                   sig(1:3,1) = sig(1:3,1)-dErep*gradij(1:3)*dxp(1)/boxcf%volu
                   sig(1:3,2) = sig(1:3,2)-dErep*gradij(1:3)*dxp(2)/boxcf%volu
                   sig(1:3,3) = sig(1:3,3)-dErep*gradij(1:3)*dxp(3)/boxcf%volu
                   if (lTPcel.EQV..true.) then
                      sigc(1:3,1,koo) =sigc(1:3,1,koo) -0.5*dErep*gradij(1:3)*dxp(1)*celcf%noxyz/boxcf%volu
                      sigc(1:3,2,koo) =sigc(1:3,2,koo) -0.5*dErep*gradij(1:3)*dxp(2)*celcf%noxyz/boxcf%volu
                      sigc(1:3,3,koo) =sigc(1:3,3,koo) -0.5*dErep*gradij(1:3)*dxp(3)*celcf%noxyz/boxcf%volu
                      sigc(1:3,1,ko1) =sigc(1:3,1,ko1) -0.5*dErep*gradij(1:3)*dxp(1)*celcf%noxyz/boxcf%volu
                      sigc(1:3,2,ko1) =sigc(1:3,2,ko1) -0.5*dErep*gradij(1:3)*dxp(2)*celcf%noxyz/boxcf%volu
                      sigc(1:3,3,ko1) =sigc(1:3,3,ko1) -0.5*dErep*gradij(1:3)*dxp(3)*celcf%noxyz/boxcf%volu
                   end if

                endif
             end if



          end do loop1at2
       end do loop1cel
    end do loop1at1


    ! calcul et stockage de Eembi et dEembi
    loop2at1: do i=1,atcf%im
       if (typ_and_pot(atcf%ityp(i),ipotentiel).eqv..false.)cycle
       iti=atcf%ityp(i)
       k=Int((tabdensity(i)-rhomin(iti))*inv_ktorho(iti))
!       write(6,*)i,iti,k,tabdensity(i),rhomin(iti),inv_ktorho(iti)
!       write(110,'(2I8,3G17.8)')i,k,tabdensity(i),rhomin(iti),inv_ktorho(iti)
       if(k.gt.ngrid) then
          write(6,*)k, ngrid, 'k> ngrid ; augmenter le facteur multiplicatif de rhomax dans calpo'
          write(6,*)'densityi',k,ngrid,densityi
          call arret_ndm
       end if
       drk=tabdensity(i)-(rhomin(iti)+k*ktorho(iti))
       Eembi = eamglue(1,iti,k) + drk*( eamglue(2,iti,k) + drk*( eamglue(3,iti,k) + drk*eamglue(4,iti,k) ) )
!       write(120,'(2I8,4G17.8)')i,k, eamglue(1,iti,k) , eamglue(2,iti,k),eamglue(3,iti,k),eamglue(4,iti,k)

       if(lprteat.EQV..true.)then
          select type (atcf)
          class is (atom_config_e)
             atcf%eat(i)=atcf%eat(i)+Eembi
          end select
       end if
       potisglue = potisglue+Eembi

       tabdensity(i) = eamglue(2,iti,k) + drk*( 2.0*eamglue(3,iti,k) + 3.0*drk*eamglue(4,iti,k) )
    end do loop2at1
    


#ifdef PARA

    if (nprocspace.gt.1) then
       call maj_tabdensity_ftm(tabdensity,atcf%imm,celcf%nato,atcf%num_at_glob,psc,atcf%im)
    end if

    !    write(3000+i,*)it
!    do i=1,im
!       write(6,*)i,num_at_glob(i),tabdensity(i)
!    end do
#endif

!    tabdensity=0
    !boucle des forces

    loop3at1: do i=1,atcf%im
       if (typ_and_pot(atcf%ityp(i),ipotentiel).eqv..false.)cycle

       ! --- Calcul de la densite sur i ---    

       koo = atcf%ielat(i)                          ! Numero de la cellule
       iti = atcf%ityp(i)
       ! pour chaque cel. voisine
       loop2cel:   do i1 = 0, celcf%ncelvois(koo)
          ko1 = celcf%ncel(koo,i1)
          ! pour chaque atome ds la cel. voisine
          loop2at2: do i2 = 1, celcf%nato(ko1)
             j = celcf%atincel(i2,ko1)
             if (typ_pot_pair(ipo(atcf%ityp(i),atcf%ityp(j))).ne.ipotentiel) cycle
             itj=atcf%ityp(j)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             !CRC             if(i.eq.j) cycle
#ifdef PARA

             ! Methode pour ne prendre qu'une seule fois en compte
             ! le couple i,j en paralle :
             ! - i est necesairement local (boucle i<=im)
             ! - si j est local on ne retient que le couple i<j
             ! - si j n'est pas local, le couple n'est par definition
             !   pris qu'une fois puisque i est local
             if (nprocspace.gt.1) then

                if (j.le.atcf%im) then

                   ! les deux atomes sont locaux
                   if (atcf%num_at_glob(i).ge.atcf%num_at_glob(j)) cycle !terme deja calcule
                else
                   ! j n'est pas local, on ne fait le calcul normal           
                endif
             else
                if (atcf%num_at_glob(i).ge.atcf%num_at_glob(j)) cycle !terme dÃ£Â©ja calculÃ£Â©
             end if
                
#else
             if (atcf%num_at_glob(i).ge.atcf%num_at_glob(j)) cycle !terme dÃ£Â©ja calculÃ£Â©
#endif

           call vect_dist(atcf,celcf,boxcf,i,j,VJI=dxp,indcv=i1, lperiod=boxcf%lperiod,rum=rue,linter=linter,dist=r)
           if(.not.linter) cycle

             if (r.eq.zero) & 
                  write(*,*) '2. WARNING IN calfoeamcell TWO ATOMS VERY CLOSE i ,j , dist(angst)', i ,j , r*angst

             k=Int(r*inv_ktor)
             drk=r-k*ktor
             gradij(1:3) = dxp(1:3)/r
             Femb = ( eamrho(2,itj,k) + drk*( 2.0*eamrho(3,itj,k) + 3.0*drk*eamrho(4,itj,k) ) )*tabdensity(i) &
                  + ( eamrho(2,iti,k) + drk*( 2.0*eamrho(3,iti,k) + 3.0*drk*eamrho(4,iti,k) ) )*tabdensity(j)
             atcf%fp(1:3,i) = atcf%fp(1:3,i) - Femb*gradij(1:3)
             atcf%fp(1:3,j) = atcf%fp(1:3,j) + Femb*gradij(1:3)

             if (test_sigma) then                   
                if (atcf%num_at_glob(i).lt.atcf%num_at_glob(j)) then
                   sig(1:3,1) = sig(1:3,1) - Femb*gradij(1:3)*dxp(1)/boxcf%volu
                   sig(1:3,2) = sig(1:3,2) - Femb*gradij(1:3)*dxp(2)/boxcf%volu
                   sig(1:3,3) = sig(1:3,3) - Femb*gradij(1:3)*dxp(3)/boxcf%volu
                   if (lTPcel.EQV..true.) then
                      sigc(1:3,1,koo) =sigc(1:3,1,koo) - 0.5*Femb*gradij(1:3)*dxp(1)*celcf%noxyz/boxcf%volu
                      sigc(1:3,2,koo) =sigc(1:3,2,koo) - 0.5*Femb*gradij(1:3)*dxp(2)*celcf%noxyz/boxcf%volu
                      sigc(1:3,3,koo) =sigc(1:3,3,koo) - 0.5*Femb*gradij(1:3)*dxp(3)*celcf%noxyz/boxcf%volu
                      sigc(1:3,1,ko1) =sigc(1:3,1,ko1) - 0.5*Femb*gradij(1:3)*dxp(1)*celcf%noxyz/boxcf%volu
                      sigc(1:3,2,ko1) =sigc(1:3,2,ko1) - 0.5*Femb*gradij(1:3)*dxp(2)*celcf%noxyz/boxcf%volu
                      sigc(1:3,3,ko1) =sigc(1:3,3,ko1) - 0.5*Femb*gradij(1:3)*dxp(3)*celcf%noxyz/boxcf%volu
                   end if

                endif
             end if

          end do loop2at2
       end do loop2cel
    end do loop3at1



#ifdef PARA

    if (nprocspace.gt.1) then
       call comm_space%sum(potisrep)
       call comm_space%sum(potisglue)
       if (test_sigma) then 
       call comm_space%sum(sig)
          if (associated(sigc)) then
             call comm_space%sum(sigc)
          endif
       endif
    end if
#endif

    potiseam=potisglue+potisrep

    return
  end SUBROUTINE calfoeamcel
end module calfoeamcel_mod
