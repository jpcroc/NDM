module calfo2ccel_mod
  USE var_pot, ONLY:alpha,csive,ipotentiel,ipo,zz,ipo,rue_pair,pot,typ_and_pot,typ_pot_pair,potis1
  USE calfocommon
  use vect_dist_mod,only:vect_dist
  USE atomconfig,only : atom_config,atom_config_d,atom_config_e,atom_config_arps
  USE cellconfig, only : cell_config
  use boxconfig,only: box_config

  implicit none
contains
  ! ***************************************************************
  subroutine calfo2ccel(atcf,celcf,boxcf)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double
    USE gen_com_m , ONLY:pi,rang,dmtype
    USE jqmod
#ifdef PARA
    USE Tpara,only:COMM_space,nprocspace
#else
    USE Tpara,only:nprocspace
#endif
    implicit none
    class(atom_config),intent(inout)::atcf
    type(cell_config),intent(in)::celcf
    type(box_config),intent(in)::boxcf

    integer :: iti, l, i, koo, i1, ko1, j, i2, itj, k
    real(double) :: aux, alp, f1, f2, f3, sk, r, phu
    real(double) :: dr,deltaepot,fcontr
    logical ::linter
    real(double)::dxp(3),gradij(3)

    sig2p=0
    ! Declarations de constantes
    aux = 23.06134575D-20
    alp = alpha/sqrt(pi)*aux

    ! Initialisation des contraintes pour le systeme global et les process
    ! terme de paires
    !pour chaque atome
    do i = 1, atcf%im
       if(typ_and_pot(atcf%ityp(i),ipotentiel).eqv..false.) cycle

       koo = atcf%ielat(i)                          ! Numero de la cellule
       iti = atcf%ityp(i)
       ! --- Calcul du second potentiel de la somme d'Ewald ---

       ! pour chaque cel. voisine
       do i1 = 0, celcf%ncelvois(koo)
          ko1 = celcf%ncel(koo,i1)
          ! pour chaque atome ds la cel. voisine
          do i2 = 1, celcf%nato(ko1)
             j = celcf%atincel(i2,ko1)
             itj = atcf%ityp(j)
             l = ipo(iti,itj)
             if (typ_pot_pair(l).ne.ipotentiel) cycle


!#ifdef PARA
             ! Methode pour ne prendre qu'une seule fois en compte
             ! le couple i,j en paralle :
             ! - i est necesairement local (boucle i<=im)
             ! - si j est local on ne retient que le couple i<j
             ! - si j n'est pas local, le couple n'est par definition
             !   pris qu'une fois puisque i est local
             if (nprocspace.gt.1) then
                if (j.le.atcf%im) then
                   ! les deux atomes sont locaux
                   if (atcf%num_at_glob(i).ge.atcf%num_at_glob(j)) cycle !terme déja calculé
                else
                   ! j n'est pas local, on fait le calcul normal           
                endif
             else
                if (atcf%num_at_glob(i).ge.atcf%num_at_glob(j)) cycle !terme déja calculé
             end if
             
             select type(atcf) ! ARPS dynamics see arps.F90
             type is (atom_config_arps)
                if (dmtype==41) then
                   if((.not.atcf%lgul(i)).and.(.not.atcf%lgul(j)))cycle
                end if
             end select
             
!#else
!             if (num_at_glob(i).ge.num_at_glob(j)) cycle !terme déja calculé
!#endif
           call vect_dist(atcf,celcf,boxcf,i,j,VJI=dxp,indcv=i1, lperiod=boxcf%lperiod,rum=rue_pair(l),linter=linter,dist=r)
           if(.not.linter) cycle
           
           sk = r/csive
           k = int(sk)
           gradij(1:3) = dxp(1:3)/r

             ! spline
             dr = r-float(k)*csive
             phu = -1.0*(pot(2,l,k)+(2.0*pot(3,l,k)+3.0*pot(4,l,k)*dr)*dr)
             deltaepot=0.5*(pot(1,l,k)+pot(2,l,k)*dr+pot(3,l,k)*dr**2+pot(4,l,k)*dr**3)
             f1 = phu*gradij(1)
             f2 = phu*gradij(2)
             f3 = phu*gradij(3)
             fcontr=sqrt(f1**2+f2**2+f3**2)
             atcf%fp(1,i) = atcf%fp(1,i)+f1
             atcf%fp(2,i) = atcf%fp(2,i)+f2
             atcf%fp(3,i) = atcf%fp(3,i)+f3
             potis1 = potis1+deltaepot
#ifdef PARA
             if (j.le.atcf%im) then
#endif
                potis1 = potis1+deltaepot
#ifdef PARA
             endif
#endif
             atcf%fp(1,j) = atcf%fp(1,j)-f1
             atcf%fp(2,j) = atcf%fp(2,j)-f2
             atcf%fp(3,j) = atcf%fp(3,j)-f3
!!
             if (lprteat) then
                select type (atcf)
                class is (atom_config_e)
                   atcf%eat(i) = atcf%eat(i)+deltaepot
                   atcf%eat(j) = atcf%eat(j)+deltaepot
                end select
                !              end if
             end if

             ! calcul des contraintes

             if (test_sigma) then
                if (atcf%num_at_glob(i).lt.atcf%num_at_glob(j)) then
                   sig2p(1,1) = sig2p(1,1)+phu*gradij(1)*dxp(1)/boxcf%volu
                   sig2p(1,2) = sig2p(1,2)+phu*gradij(1)*dxp(2)/boxcf%volu
                   sig2p(1,3) = sig2p(1,3)+phu*gradij(1)*dxp(3)/boxcf%volu
                   sig2p(2,1) = sig2p(2,1)+phu*gradij(2)*dxp(1)/boxcf%volu
                   sig2p(2,2) = sig2p(2,2)+phu*gradij(2)*dxp(2)/boxcf%volu
                   sig2p(2,3) = sig2p(2,3)+phu*gradij(2)*dxp(3)/boxcf%volu
                   sig2p(3,1) = sig2p(3,1)+phu*gradij(3)*dxp(1)/boxcf%volu
                   sig2p(3,2) = sig2p(3,2)+phu*gradij(3)*dxp(2)/boxcf%volu
                   sig2p(3,3) = sig2p(3,3)+phu*gradij(3)*dxp(3)/boxcf%volu
                   if (lcalcsigc.EQV..true.) then
                      sigc(1,1,koo) = sigc(1,1,koo)+0.5*phu*gradij(1)*dxp(1)*celcf%noxyz/boxcf%volu
                      sigc(1,2,koo) = sigc(1,2,koo)+0.5*phu*gradij(1)*dxp(2)*celcf%noxyz/boxcf%volu
                      sigc(1,3,koo) = sigc(1,3,koo)+0.5*phu*gradij(1)*dxp(3)*celcf%noxyz/boxcf%volu
                      sigc(2,1,koo) = sigc(2,1,koo)+0.5*phu*gradij(2)*dxp(1)*celcf%noxyz/boxcf%volu
                      sigc(2,2,koo) = sigc(2,2,koo)+0.5*phu*gradij(2)*dxp(2)*celcf%noxyz/boxcf%volu
                      sigc(2,3,koo) = sigc(2,3,koo)+0.5*phu*gradij(2)*dxp(3)*celcf%noxyz/boxcf%volu
                      sigc(3,1,koo) = sigc(3,1,koo)+0.5*phu*gradij(3)*dxp(1)*celcf%noxyz/boxcf%volu
                      sigc(3,2,koo) = sigc(3,2,koo)+0.5*phu*gradij(3)*dxp(2)*celcf%noxyz/boxcf%volu
                      sigc(3,3,koo) = sigc(3,3,koo)+0.5*phu*gradij(3)*dxp(3)*celcf%noxyz/boxcf%volu
                      !                      if (j.le.atcf%im) then
                      sigc(1,1,ko1) = sigc(1,1,ko1)+0.5*phu*gradij(1)*dxp(1)*celcf%noxyz/boxcf%volu
                      sigc(1,2,ko1) = sigc(1,2,ko1)+0.5*phu*gradij(1)*dxp(2)*celcf%noxyz/boxcf%volu
                      sigc(1,3,ko1) = sigc(1,3,ko1)+0.5*phu*gradij(1)*dxp(3)*celcf%noxyz/boxcf%volu
                      sigc(2,1,ko1) = sigc(2,1,ko1)+0.5*phu*gradij(2)*dxp(1)*celcf%noxyz/boxcf%volu
                      sigc(2,2,ko1) = sigc(2,2,ko1)+0.5*phu*gradij(2)*dxp(2)*celcf%noxyz/boxcf%volu
                      sigc(2,3,ko1) = sigc(2,3,ko1)+0.5*phu*gradij(2)*dxp(3)*celcf%noxyz/boxcf%volu
                      sigc(3,1,ko1) = sigc(3,1,ko1)+0.5*phu*gradij(3)*dxp(1)*celcf%noxyz/boxcf%volu
                      sigc(3,2,ko1) = sigc(3,2,ko1)+0.5*phu*gradij(3)*dxp(2)*celcf%noxyz/boxcf%volu
                      sigc(3,3,ko1) = sigc(3,3,ko1)+0.5*phu*gradij(3)*dxp(3)*celcf%noxyz/boxcf%volu
                   end if
                endif
             end if
          end do  ! fin i2=j
    end do ! fin i1=koo


    end do ! fin i

#ifdef PARA

    if (nprocspace.gt.1) then
       call comm_space%sum(potis1)
       if (test_sigma)then
       call comm_space%sum(sig2p)
       if (associated(sigc)) then
          call comm_space%sum(sigc)
          
       endif
    end if
 end if
#endif


    ! fin du calcul du terme de paire dans l'espace direct
    return
  end subroutine calfo2ccel

end module calfo2ccel_mod
