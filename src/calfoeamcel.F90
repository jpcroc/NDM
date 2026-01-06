module calfoeamcel_mod
  USE arret_ndm_mod,only:arret_ndm
  USE gen_com_m, ONLY:angst,nvat,low_limit,lperiod,zero,pi,rang
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

    USE var_pot, ONLY:ipotentiel,ngrid,rue_pot,&
         &typ_and_pot,typ_pot_pair,ipotentiel,ngrid,potiseam,potisglue,potisrep,rhomax,rhomin,eamrho,ipo,eamrep,&
         &eamglue,alpha,ntyp

#ifdef PARA
    use Tpara,only:nprocspace,para_space_config,comm_space
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
    integer:: koo,ko1,i2,i1,nvi !cel.
    integer :: k ! aux pour splines

    !    real(double) :: rue2 !coupure**2
    REAL(double), dimension(1:3) :: dxp, gradij
    real(double) :: r!distance i-j
    real(double) :: Eembi,dEembi ! potentiel et gradient de l'immersion
    real(double) :: rhoi,rhoj ! densite de i sur j et j sur i
    real(double) :: Erep,dErep ! potentiel et gradient de la repulsion de paire ij
    REAL(double) :: Femb
    real(double) ::  drk,ktor, inv_ktor
    real(double),dimension(:),allocatable::ktorho(:), inv_ktorho(:)
    real(double) :: densityi !densite totale sur i
    real(double) :: tabdensity(atcf%imm)
    real(double)::rue,alp,aux
    logical ::linter
    real(double)::sig2p(3,3),sigem(3,3)!,dxptab(1000,3),xptab(1000,3),deltadist(1000,3)
!    integer::cellv(1000),indv(1000),m,n
    !    write(6,*)'INNNI CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
    sig2p=0;sigem=0
    rue=rue_pot(ipotentiel)
    aux = 23.06134575D-20
    alp = alpha/sqrt(pi)*aux
    allocate(ktorho(ntyp))
    allocate(inv_ktorho(ntyp))
    ktor=rue/ngrid
    inv_ktor=1.d0/ktor
    do iti=1,ntyp
       if (typ_and_pot(iti,ipotentiel).eqv..false.)then
          ktorho(iti)=1
       else
          ktorho(iti)=(rhomax(iti)-rhomin(iti))/ngrid
       end if
    end do
    inv_ktorho(:) = 1.d0/ktorho(:)
    tabdensity(:)=0.
    potisrep=0.
    potisglue=0.
    !    rue2=rue**2



    loop1at1: do i=1,atcf%im

       if (typ_and_pot(atcf%ityp(i),ipotentiel).eqv..false.)cycle

       ! --- Calcul de la densite sur i ---    

       nvi=0
       densityi=0.0 ; dEembi=0.0
       koo = atcf%ielat(i)                          ! Numero de la cellule
       iti = atcf%ityp(i)
       !       write(6,*)'IXP, koo', atcf%xp(:,1), koo
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
             if (celcf%isghost(ko1)) then
                !CRC les interactions des ghost doivent toujouts être calculées

             else
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
             end if
             call vect_dist(atcf,celcf,boxcf,i,j,VJI=dxp,indcv=i1, lperiod=boxcf%lperiod,rum=rue,linter=linter,dist=r)
             if(.not.linter) cycle
!!$             nvi=nvi+1
!!$             indv(nvi)=j
!!$             cellv(nvi)=ko1
!!$             dxptab(nvi,:)=dxp
!!$             xptab(nvi,:)=atcf%xp(:,j)
!!$             deltadist(nvi,:)=celcf%deltadist(:,i1,koo)
!!$           do m=1,nvi-1
!!$              if (all(dxptab(m,:)-dxptab(nvi,:)==0)) then
!!$                 write(6,*) 'm==nvi',m,nvi
!!$                 write(6,*)'INDV',indv(m),cellv(m)
!!$                 write(6,'(A,2I4,3E15.7)')'DXP',indv(m),cellv(m),dxptab(m,1:3)
!!$                 write(6,*)'Xptab',xptab(m,1:3)
!!$                 write(6,*)'celtadist',deltadist(m,:)
!!$                 
!!$                 write(6,*)'INDV',indv(nvi),cellv(nvi)
!!$                 write(6,'(A,2I4,3E15.7)')'DXP',indv(nvi),cellv(nvi),dxptab(nvi,1:3)
!!$                 write(6,*)'Xptab',xptab(nvi,1:3)
!!$                 write(6,*)'celtadist',deltadist(nvi,:)
!!$                 stop
!!$              end if
!!$           end do
             k=min(ngrid,Int(r*inv_ktor))
             gradij(1:3) = dxp(1:3)/r
             drk=r-k*ktor
             !             rhoj = eamrho(1,itj,k) + drk*( eamrho(2,itj,k) + drk*( eamrho(3,itj,k) + drk*eamrho(4,itj,k) ) )  !rho de j sur i
             !             densityi=densityi+rhoj

             rhoj = eamrho(1,itj,k) + drk*( eamrho(2,itj,k) + drk*( eamrho(3,itj,k) + drk*eamrho(4,itj,k) ) )  !rho de j sur i
             tabdensity(i)=tabdensity(i)+rhoj
             rhoi = eamrho(1,iti,k) + drk*( eamrho(2,iti,k) + drk*( eamrho(3,iti,k) + drk*eamrho(4,iti,k) ) )  !rho de i sur j
             if (.not.(celcf%isghost(ko1)))  tabdensity(j)=tabdensity(j)+rhoi
!!$             write(100,'(2I3,3G17.8)')i,j,r,tabdensity(i),tabdensity(j)
!!$             write(100,'(2I3,4G17.8)')i,j,eamrho(1,iti,k) , eamrho(2,iti,k),eamrho(3,iti,k),eamrho(4,iti,k)
!!$             write(100,'(2I3,4G17.8)')i,j,eamrho(1,itj,k) , eamrho(2,itj,k),eamrho(3,itj,k),eamrho(4,itj,k)
             !           nvi=nvi+1
             !                        dxpij(1:3,nvi)=dxp(1:3)
             !                        jvi(nvi)=j
             !                        rij(nvi)=r            

             !terme de repulsion 
             l = ipo(iti,itj)
             Erep = eamrep(1,l,k) + drk*( eamrep(2,l,k) + drk*( eamrep(3,l,k) + drk*eamrep(4,l,k) ) )
             if(lprteat.EQV..true.)then
                select type (atcf)
                class is (atom_config_e)
                   atcf%eat(i)=atcf%eat(i)+Erep/2.d0
                   if (.not.(celcf%isghost(ko1)))then
                      if (j.le.atcf%im) atcf%eat(j)=atcf%eat(j)+Erep/2.d0
                   end if
                end select
             end if
             dErep = eamrep(2,l,k) + drk*( 2.0*eamrep(3,l,k) + 3.0*drk*eamrep(4,l,k) )

             if (celcf%isghost(ko1)) then
                potisrep = potisrep+0.5*Erep                
             else
                if (atcf%num_at_glob(i).lt.atcf%num_at_glob(j)) then
                   potisrep = potisrep+Erep
                endif
             end if
             atcf%fp(1:3,i)=atcf%fp(1:3,i)-dErep*gradij(1:3)
             if (.not.(celcf%isghost(ko1)))             atcf%fp(1:3,j)=atcf%fp(1:3,j)+dErep*gradij(1:3)

             if (test_sigma) then        
                if (celcf%isghost(ko1)) then
                   sig2p(1:3,1) = sig2p(1:3,1)-0.5*dErep*gradij(1:3)*dxp(1)/boxcf%volu
                   sig2p(1:3,2) = sig2p(1:3,2)-0.5*dErep*gradij(1:3)*dxp(2)/boxcf%volu
                   sig2p(1:3,3) = sig2p(1:3,3)-0.5*dErep*gradij(1:3)*dxp(3)/boxcf%volu
                   if (lcalcsigc.EQV..true.) then
                      sigc(1:3,1,koo) =sigc(1:3,1,koo) -0.5*dErep*gradij(1:3)*dxp(1)*celcf%noxyzact/boxcf%volu
                      sigc(1:3,2,koo) =sigc(1:3,2,koo) -0.5*dErep*gradij(1:3)*dxp(2)*celcf%noxyzact/boxcf%volu
                      sigc(1:3,3,koo) =sigc(1:3,3,koo) -0.5*dErep*gradij(1:3)*dxp(3)*celcf%noxyzact/boxcf%volu
                   end if
                   
                else
                   if (atcf%num_at_glob(i).lt.atcf%num_at_glob(j)) then
                      sig2p(1:3,1) = sig2p(1:3,1)-dErep*gradij(1:3)*dxp(1)/boxcf%volu
                      sig2p(1:3,2) = sig2p(1:3,2)-dErep*gradij(1:3)*dxp(2)/boxcf%volu
                      sig2p(1:3,3) = sig2p(1:3,3)-dErep*gradij(1:3)*dxp(3)/boxcf%volu
                      if (lcalcsigc.EQV..true.) then
                         sigc(1:3,1,koo) =sigc(1:3,1,koo) -0.5*dErep*gradij(1:3)*dxp(1)*celcf%noxyzact/boxcf%volu
                         sigc(1:3,2,koo) =sigc(1:3,2,koo) -0.5*dErep*gradij(1:3)*dxp(2)*celcf%noxyzact/boxcf%volu
                         sigc(1:3,3,koo) =sigc(1:3,3,koo) -0.5*dErep*gradij(1:3)*dxp(3)*celcf%noxyzact/boxcf%volu
                         sigc(1:3,1,ko1) =sigc(1:3,1,ko1) -0.5*dErep*gradij(1:3)*dxp(1)*celcf%noxyzact/boxcf%volu
                         sigc(1:3,2,ko1) =sigc(1:3,2,ko1) -0.5*dErep*gradij(1:3)*dxp(2)*celcf%noxyzact/boxcf%volu
                         sigc(1:3,3,ko1) =sigc(1:3,3,ko1) -0.5*dErep*gradij(1:3)*dxp(3)*celcf%noxyzact/boxcf%volu
                      end if
                   endif
                end if
             end if

          end do loop1at2
       end do loop1cel
!!$       write(16,*)'NVI',i
!!$       write(16,*)i,nvi
!!$       do j=1,nvi
!!$          write(16,*)
!!$          write(16,*)'INDV',indv(j),cellv(j)
!!$          write(16,'(A,2I4,3E15.7)')'DXP',indv(j),cellv(j),dxptab(j,1:3)
!!$          write(16,*)'Xptab',xptab(j,1:3)
!!$
!!$       end do
!!$       stop
    end do loop1at1


    ! calcul et stockage de Eembi et dEembi
    loop2at1: do i=1,atcf%im
       if (typ_and_pot(atcf%ityp(i),ipotentiel).eqv..false.)cycle
       iti=atcf%ityp(i)

       !       k=Int((tabdensity(i)-rhomin(iti))*inv_ktorho(iti))
       k=min(ngrid,Int((tabdensity(i)-rhomin(iti))*inv_ktorho(iti)))
       !       if (tabdensity(i).gt.rhmax) then
       !          rhmax=tabdensity(i)
       !          imax=i; kmax=k ; itimax=iti
       !       end if

       !       write(6,*)i,iti,k,tabdensity(i),rhomin(iti),inv_ktorho(iti)
       !       write(110,'(2I8,3G17.8)')i,k,tabdensity(i),rhomin(iti),inv_ktorho(iti)
       !       if(k.gt.ngrid) then
       !          write(6,*)k, ngrid, 'k> ngrid ; augmenter le facteur multiplicatif de rhomax dans calpo OU MAX DENS POUR CRG '
       !          write(6,*)'densityi',k,ngrid,tabdensity(i)
       !          call arret_ndm
       !       end if
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
    !    write(6,'(A,2I2,I4,G17.5,I7)')'rhm', rang, itimax,imax,rhmax,kmax
    !  end block


#ifdef PARA

    if (nprocspace.gt.1) then
       call maj_tabdensity_ftm(tabdensity,atcf%imm,atcf%num_at_glob,psc,atcf%im)
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
             if (celcf%isghost(ko1)) then
                !CRC les interactions des ghost doivent toujouts être calculées

             else

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
             end if
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
             if (.not.(celcf%isghost(ko1)))             atcf%fp(1:3,j) = atcf%fp(1:3,j) + Femb*gradij(1:3)

             if (test_sigma) then        
                if (celcf%isghost(ko1)) then
                   sigem(1:3,1) = sigem(1:3,1) - 0.5*Femb*gradij(1:3)*dxp(1)/boxcf%volu
                   sigem(1:3,2) = sigem(1:3,2) - 0.5*Femb*gradij(1:3)*dxp(2)/boxcf%volu
                   sigem(1:3,3) = sigem(1:3,3) - 0.5*Femb*gradij(1:3)*dxp(3)/boxcf%volu
                   if (lcalcsigc.EQV..true.) then
                      sigc(1:3,1,koo) =sigc(1:3,1,koo) - 0.5*Femb*gradij(1:3)*dxp(1)*celcf%noxyzact/boxcf%volu
                      sigc(1:3,2,koo) =sigc(1:3,2,koo) - 0.5*Femb*gradij(1:3)*dxp(2)*celcf%noxyzact/boxcf%volu
                      sigc(1:3,3,koo) =sigc(1:3,3,koo) - 0.5*Femb*gradij(1:3)*dxp(3)*celcf%noxyzact/boxcf%volu
                   end if
                   
                else
                   if (atcf%num_at_glob(i).lt.atcf%num_at_glob(j)) then
                      sigem(1:3,1) = sigem(1:3,1) - Femb*gradij(1:3)*dxp(1)/boxcf%volu
                      sigem(1:3,2) = sigem(1:3,2) - Femb*gradij(1:3)*dxp(2)/boxcf%volu
                      sigem(1:3,3) = sigem(1:3,3) - Femb*gradij(1:3)*dxp(3)/boxcf%volu
                      if (lcalcsigc.EQV..true.) then
                         sigc(1:3,1,koo) =sigc(1:3,1,koo) - 0.5*Femb*gradij(1:3)*dxp(1)*celcf%noxyzact/boxcf%volu
                         sigc(1:3,2,koo) =sigc(1:3,2,koo) - 0.5*Femb*gradij(1:3)*dxp(2)*celcf%noxyzact/boxcf%volu
                         sigc(1:3,3,koo) =sigc(1:3,3,koo) - 0.5*Femb*gradij(1:3)*dxp(3)*celcf%noxyzact/boxcf%volu
                         sigc(1:3,1,ko1) =sigc(1:3,1,ko1) - 0.5*Femb*gradij(1:3)*dxp(1)*celcf%noxyzact/boxcf%volu
                         sigc(1:3,2,ko1) =sigc(1:3,2,ko1) - 0.5*Femb*gradij(1:3)*dxp(2)*celcf%noxyzact/boxcf%volu
                         sigc(1:3,3,ko1) =sigc(1:3,3,ko1) - 0.5*Femb*gradij(1:3)*dxp(3)*celcf%noxyzact/boxcf%volu
                      end if
                   endif
                end if
             end if
             
          end do loop2at2
       end do loop2cel
    end do loop3at1



#ifdef PARA

    if (nprocspace.gt.1) then
       call comm_space%sum(potisrep)
       call comm_space%sum(potisglue)
       if (test_sigma) then 
          !          call comm_space%sum(sig)
          call comm_space%sum(sig2p)
          call comm_space%sum(sigem)

          if (associated(sigc)) then
             call comm_space%sum(sigc)
          endif
       endif
    end if
#endif
    if (test_sigma)sigcalfo=sigcalfo+sig2p+sigem
!!$    if (rang==0)    write(6,*)
!!$    if (rang==0)    write(6,*)'sig2p',sig2p
!!$    if (rang==0)    write(6,*)
!!$    if (rang==0)    write(6,*)'sigem',sigem
    potiseam=potisglue+potisrep
!    call atcf%print
    return
  end SUBROUTINE calfoeamcel
end module calfoeamcel_mod
