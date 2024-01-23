module force_tersoff_cel_mod
  USE arret_ndm_mod,only:arret_ndm
  USE cryst_to_cart_mod,only: cryst_to_cart
  USE gen_com_m, ONLY:lcalcjq,potistersoff,potiszbl
  USE calfocommon
  use vect_dist_mod,only:vect_dist
  USE atomconfig,only : atom_config,atom_config_d,atom_config_e
  USE cellconfig, only : cell_config
  use boxconfig,only: box_config

  implicit none
contains
! ***************************************************************
  subroutine force_tersoff_cel(atcf,celcf,boxcf,psc)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double
    USE var_pot, ONLY:npair,csive,typ_and_pot,ipo,typ_pot_pair,roff2,pot
    USE jqmod
    USE force_tersoff_facteurs
#ifdef PARA
    USE mod_para,only:maj_fp_frt
    use Tpara,only:COMM_space,nprocspace,para_space_config
#else
    USE Tpara,only:nprocspace,para_space_config
#endif


    implicit none
    class(atom_config),intent(inout)::atcf
    type(cell_config),intent(in)::celcf
    type(box_config),intent(in)::boxcf
    type(para_space_config)::psc
    !-----------------------------------------------
    integer :: i,j,k,nk,n_voisin,l,ij,ik, idv, m, moi
    integer , dimension(32) :: indice   ! Recense le nombre de voisins
    integer :: icelnumber,jcelvois,jcelnumber,jnumber, kcelvois,kcelnumber,knumber
    real(double) ::  rij,  rik, sui_ij, bij, n, v_ij
    real(double) :: fc_rij, dfc_rij, fr_rij, fa_rij, fc_rik, dfc_rik, paire_ij, triplet_ij, triplet_ik 
    real(double) :: exponentiel, cos_theta, g_cos, dg_cos, ER1, ER2, ER3
    real(double) :: Scal_FiVi, Scal_FjVj, Scal_FijVi, Scal_FijVj, Scal_FikVi, Scal_FikVk

    real(double) , dimension(15,6) :: tmp
    real(double) , dimension(15,3) :: tmp1
    real(double) , dimension(1,3) :: cvij, cvik
    real(double):: coupR(npair),sigT(3,3)

    real(double) :: phu,sk,dr
    integer::kk
    logical::linter,linterik

    select case (ipotentiel) 
    case(13)
       coupR(:)=ster(:)
    case(14,15)
       coupR(:)=rter(:)+ster(:)
    end select

    !-----------------------------------------------
    !INIALISATION
    moi =0
    !  do i=1,im
    !     fp(1,i)=0; fp(2,i)=0; fp(3,i)=0
    !  end do
    !  jq(:)=0; sigcalfo(:,:)=0

    ER1=0. ;  ER2=0. ;  ER3=0.
    !         write(6,*)'boite quelc'
!    call cryst_to_cart(imm,xp,bg,-1)

    !  write(6,*)sigcalfo
    !  write(6,*)
    do i=1,atcf%im
       if(typ_and_pot(atcf%ityp(i),ipotentiel).eqv..false.) cycle
       v_ij = 0

       icelnumber = atcf%ielat(i)

       do jcelvois=0, celcf%ncelvois(icelnumber)
          jcelnumber = celcf%ncel(icelnumber,jcelvois)
          do jnumber =1, celcf%nato(jcelnumber)
             j = celcf%atincel(jnumber,jcelnumber)
             if (j==i) then  !Cette condition n'est pas necessaire si JP. fais correctement sa table
                cycle
             else
                ij=ipo(atcf%ityp(i),atcf%ityp(j))
                !              write(6,*)ij,typ_pot_pair(ij),ipotentiel
!!$              if (typ_pot_pair(ij).ne.ipotentiel) cycle
!!$              do l=1,3
!!$                 Xij(l)= xp(l,i)-xp(l,j)
!!$              end do
!!$
!!$              do l=1,3
!!$                 !                        if (abs(Xij(l))>0.5) Xij(l)= Xij(l)-sign(1.0d0,xp(l,i))
!!$                 if (Xij(l)>0.5) Xij(l)=Xij(l)-1.
!!$                 if (Xij(l)<(-0.5)) Xij(l)=Xij(l)+1.
!!$                 cvij(1,l) = Xij(l)
!!$              end do
!!$              call cryst_to_cart(1,cvij,at,1)
!!$              rij2=cvij(1,1)**2+cvij(1,2)**2+cvij(1,3)**2
!!$              if (rij2>Coupr(ij)**2) then         !borne sup de Lisa Porter 89
!!$                 !if (rij2>(Rter(ij)+Coupr(ij))**2) then   !borne sup de Tersoff 88
!!$                 cycle
!!$              else
!!$
                call vect_dist(atcf,celcf,boxcf,i,j,VJI=cvij(1,:), lperiod=boxcf%lperiod,rum=CoupR(ij)&
                     &,linter=linter,dist=rij)
                if (.not.linter) cycle

                n=nter(atcf%ityp(i))
                !                 rij=sqrt(rij2)
                sui_ij = 0
                n_voisin = 0

                do kcelvois=0, celcf%ncelvois(icelnumber)
                   kcelnumber = celcf%ncel(icelnumber,kcelvois)
                   do knumber =1, celcf%nato(kcelnumber)
                      k = celcf%atincel(knumber,kcelnumber)

                      if (k==j .or. k==i) then
                         cycle
                      else
                         ik=ipo(atcf%ityp(i),atcf%ityp(k))
                         if (typ_pot_pair(ik).ne.ipotentiel) cycle
                         call vect_dist(atcf,celcf,boxcf,i,k,VJI=cvik(1,:), lperiod=boxcf%lperiod,rum=CoupR(ik)&
                              &,linter=linterik,dist=rik)
                         if (.not.linterik) cycle

                         ! CONDITIONS PERIODIQUES
                         !write(6,*)'XP i j k',i,j,k
                         !write(6,*)xp(:,i)
                         !write(6,*)xp(:,j)
                         !write(6,*)xp(:,k)
!!$                          do l=1,3
!!$                             Xik(l)= xp(l,i)-xp(l,k)
!!$                          end do
!!$                          do l=1,3
!!$                             if (XiK(l)>0.5) Xik(l)=Xik(l)-1.
!!$                             if (XiK(l)<(-0.5)) Xik(l)=Xik(l)+1.
!!$                             cvik(1,l) = Xik(l)
!!$                          end do
!!$                          !                       write(6,*)'cvik'
!!$                          !                       write(6,*)cvik(1,:)
!!$                          call cryst_to_cart(1,cvik,at,1)
!!$                          rik2=cvik(1,1)**2+cvik(1,2)**2+cvik(1,3)**2
!!$                          if (rik2>Coupr(ik)**2) then               !borne sup de Lisa Porter 89
!!$                             !if (rik2>(Rter(ik)+Coupr(ik))**2) then     !borne sup de Tersoff 88
!!$                             cycle
!!$                          else
                         ! ne prend en compte que les k voisin de la paire ij
                         n_voisin = n_voisin + 1
                         indice(n_voisin) = k
                         !                             rik = sqrt(rik2)
                         tmp(n_voisin,1) = rik                                                   ! rik
                         call facteur_amortissement (rik, ik, fc= tmp(n_voisin,2) )              ! fc(rik)
                         tmp(n_voisin,3) = exp(lambda3(ik)**3*(rij-rik)**3)                      ! exponentiel
                         tmp(n_voisin,4) = (cvij(1,1)*cvik(1,1)+cvij(1,2)*cvik(1,2)+cvij(1,3)*cvik(1,3))/(rij*rik)
                         ! cos theta ijk
                         call facteur_angulaire(tmp(n_voisin,4), atcf%ityp(i), g = tmp(n_voisin,5))   ! g(cos_theta)
                         tmp1(n_voisin,1)=cvik(1,1) ; tmp1(n_voisin,2)=cvik(1,2) ; tmp1(n_voisin,3)=cvik(1,3)

                         sui_ij = sui_ij + tmp(n_voisin,2)*tmp(n_voisin,5)*tmp(n_voisin,3)       ! sui ij

                         !                          end if
                      end if
                   end do
                end do

                ! Tersoff standard                        bij = psi(ij)*(1.+beta(ityp(i))**n*sui_ij**n)**(-1./(2.*n)) 
                !Brenner
                bij = psi(ij)*(1.+beta(atcf%ityp(i))**n*sui_ij**n)**(-deltater(atcf%ityp(i)))

                call facteur_amortissement (rij, ij, fc= fc_rij )
                call facteur_amortissement (rij, ij, dfc= dfc_rij)
                fr_rij = fr(rij,Ater(ij),lambda1(ij))
                fa_rij = fa(rij,Bter(ij),lambda2(ij))

                v_ij = v_ij + fc_rij*(fr_rij+bij*fa_rij)
                ! DEBUT CALCUL DE  FORCE $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  
                if(n_voisin==0) then
                   Scal_FiVi=0. ; Scal_FjVj=0.
                   do l=1,3
                      paire_ij = -0.5*dfc_rij*(fr_rij+bij*fa_rij)/rij*cvij(1,l) & 
                           +0.5*fc_rij*(lambda1(ij)*fr_rij+lambda2(ij)*bij*fa_rij)/rij*cvij(1,l)

                      atcf%fp(l,i) = atcf%fp(l,i) + paire_ij
                      atcf%fp(l,j) = atcf%fp(l,j) - paire_ij

!!$                       if (lcalcjq) then
!!$                          Scal_FjVj=Scal_FjVj - paire_ij*vp(l,j)
!!$                       end if
                      !Contrainte
                      if (test_sigma) then 
                         do m=1,3
                            sigT(l,m)=sigT(l,m) + paire_ij*cvij(1,m)/boxcf%volu
                            if (lTPcel.EQV..true.) then
                               sigc(l,m,icelnumber) = sigc(l,m,icelnumber) + 0.5*paire_ij*cvij(1,m)*celcf%noxyz/boxcf%volu
                               sigc(l,m,jcelnumber) = sigc(l,m,jcelnumber) + 0.5*paire_ij*cvij(1,m)*celcf%noxyz/boxcf%volu
                            end if
                         end do



                      end if
                   end do
                   if (lcalcjq) then
                      do l=1,3
                         jq(l) = jq(l) - Scal_FjVj*(cvij(1,l))
                      end do
                   end if

                else
                   Scal_FiVi=0. ; Scal_FjVj=0.
                   do l=1,3
                      paire_ij = -0.5*dfc_rij*(fr_rij+bij*fa_rij)/rij*cvij(1,l) & 
                           +0.5*fc_rij*(lambda1(ij)*fr_rij+lambda2(ij)*bij*fa_rij)/rij*cvij(1,l)

                      atcf%fp(l,i) = atcf%fp(l,i) + paire_ij
                      atcf%fp(l,j) = atcf%fp(l,j) - paire_ij

!!$                       if (lcalcjq) then
!!$                          Scal_FjVj=Scal_FjVj - paire_ij*vp(l,j)
!!$                       end if
                      !Contrainte
                      if (test_sigma) then 
                         do m=1,3
                            sigT(l,m)=sigT(l,m) + paire_ij*cvij(1,m)/boxcf%volu
                            if (lTPcel.EQV..true.) then
                               sigc(l,m,icelnumber) = sigc(l,m,icelnumber) + paire_ij*cvij(1,m)*celcf%noxyz/boxcf%volu
                               sigc(l,m,jcelnumber) = sigc(l,m,jcelnumber) + paire_ij*cvij(1,m)*celcf%noxyz/boxcf%volu
                            end if

                         end do
                      end if
                   end do
                   !write(6,*)i,j,'P'
                   !write(6,*)sig
                   !write(6,*)

                   if (lcalcjq) then
                      do l=1,3
                         jq(l) = jq(l) - Scal_FjVj*(cvij(1,l))
                      end do
                   end if
                   !*********************************************************************************************
                   do nk=1,n_voisin
                      !write(6,*)'NV',n_voisin
                      k=indice(nk) 
                      !write(6,*)k,nk
                      rik = tmp(nk,1)
                      fc_rik = tmp(nk,2)
                      ik=ipo(atcf%ityp(i),atcf%ityp(k))
                      call facteur_amortissement (rik, ik, dfc= dfc_rik)
                      exponentiel = tmp(nk,3)
                      cos_theta = tmp(nk,4)
                      g_cos = tmp(nk,5)
                      call facteur_angulaire (cos_theta, atcf%ityp(i), dg= dg_cos)
                      Scal_FijVi=0. ; Scal_FijVj=0. ; Scal_FikVi=0. ; Scal_FikVk=0.
                      cvik(1,:)=tmp1(nk,:)
                      do l=1,3
                         !                          cvik(1,l) = tmp1(nk,l)
                         !                          write(6,*)'cvik l',cvik(1,l),l

                         !                                 triplet_ij = 0.5*fc_rij*fa_rij*0.5*psi(ij)*beta(ityp(i))**n*sui_ij**(n-1)*&
                         !                                      (1.+beta(ityp(i))**n*sui_ij**n)**(-(1.+2.*n)/(2.*n)) * &   
                         !                                      (3.*lambda3(ik)**3*fc_rik*(rij-rik)**2*g_cos*exponentiel/rij*cvij(1,l) &
                         !                                      +fc_rik*dg_cos*exponentiel* &
                         !                                      (cvik(1,l)/(rij*rik)-cos_theta*cvij(1,l)/rij**2))


                         triplet_ij = 0.5*fc_rij*fa_rij*deltater(atcf%ityp(i))*n*psi(ij)*beta(atcf%ityp(i))**n*sui_ij**(n-1)*&
                              (1.+beta(atcf%ityp(i))**n*sui_ij**n)**(-deltater(atcf%ityp(i))-1.) * &   
                              (3.*lambda3(ik)**3*fc_rik*(rij-rik)**2*g_cos*exponentiel/rij*cvij(1,l) &
                              +fc_rik*dg_cos*exponentiel* &
                              (cvik(1,l)/(rij*rik)-cos_theta*cvij(1,l)/rij**2))


                         !                                 triplet_ik = 0.5*fc_rij*fa_rij*0.5*psi(ij)*beta(atcf%ityp(i))**n*sui_ij**(n-1)*&
                         !                                      (1.+beta(atcf%ityp(i))**n*sui_ij**n)**(-(1.+2.*n)/(2.*n)) * &
                         !                                      (-3.*lambda3(ik)**3*fc_rik*(rij-rik)**2*g_cos*exponentiel/rik*cvik(1,l) &
                         !                                      +dfc_rik*g_cos*exponentiel/rik*cvik(1,l) &
                         !                                      +fc_rik*dg_cos*exponentiel* &
                         !                                      (cvij(1,l)/(rij*rik)-cos_theta*cvik(1,l)/rik**2))

                         triplet_ik = 0.5*fc_rij*fa_rij*deltater(atcf%ityp(i))*n*psi(ij)*beta(atcf%ityp(i))**n*sui_ij**(n-1)*&
                              (1.+beta(atcf%ityp(i))**n*sui_ij**n)**(-deltater(atcf%ityp(i))-1.)* &
                              (-3.*lambda3(ik)**3*fc_rik*(rij-rik)**2*g_cos*exponentiel/rik*cvik(1,l) &
                              +dfc_rik*g_cos*exponentiel/rik*cvik(1,l) &
                              +fc_rik*dg_cos*exponentiel* &
                              (cvij(1,l)/(rij*rik)-cos_theta*cvik(1,l)/rik**2))

                         atcf%fp(l,i) = atcf%fp(l,i) + triplet_ij
                         atcf%fp(l,j) = atcf%fp(l,j) - triplet_ij
                         atcf%fp(l,i) = atcf%fp(l,i) + triplet_ik
                         atcf%fp(l,k) = atcf%fp(l,k) - triplet_ik

                         !Flux
!!$                          if (lcalcjq) then
!!$                             Scal_FijVj=Scal_FijVj - triplet_ij*vp(l,j)
!!$                             Scal_FikVk=Scal_FikVk - triplet_ik*vp(l,k)
!!$                          end if
                         !Contrainte
                         !write(6,*)'sig AV l',l
                         !write(6,*)sig
                         if (test_sigma) then 
                            do m=1,3
                               sigT(l,m)=sigT(l,m) + triplet_ij*cvij(1,m)/boxcf%volu
                               sigT(l,m)=sigT(l,m) + triplet_ik*cvik(1,m)/boxcf%volu
                               if (lTPcel.EQV..true.) then
                                  sigc(l,m,icelnumber) = sigc(l,m,icelnumber) + celcf%noxyz*0.5*(triplet_ij*cvij(1,m)/boxcf%volu  &
                                       + triplet_ik*cvik(1,m)/boxcf%volu)
                                  sigc(l,m,jcelnumber) = sigc(l,m,jcelnumber) + celcf%noxyz*0.5*triplet_ij*cvij(1,m)/boxcf%volu 
                                  sigc(l,m,kcelnumber) = sigc(l,m,kcelnumber) + celcf%noxyz*0.5*triplet_ik*cvik(1,m)/boxcf%volu
                               end if

                            end do
                         end if
                         !write(6,*)i,j,k,'T'
                         !write(6,*)sig
                         !write(6,*)'tik tik'
                         !write(6,*)triplet_ij,triplet_ik
                         !write(6,*)'cvij cvik'
                         !                       do m=1,3
                         !write(6,*)cvij(1,m),cvik(1,m)
                         !                       end do

                      end do
                      if (lcalcjq) then
                         do l=1,3
                            jq(l) = jq(l) - Scal_FijVj*(cvij(1,l)) - Scal_FikVk*(cvik(1,l))
                         end do
                      end if
                   end do
                end if
                if (rij.le.roff2(ij)) then
                   !                    write(6,*)'BINGO'
                   sk = rij/csive
                   kk = sk
                   !                    ! spline
                   dr = rij-float(kk)*csive
                   !                    write(6,*)'TZBL',rij,dr,kk,csive
                   potiszbl = potiszbl+0.5*(pot(1,ij,kk)+ rij*(dr*(pot(2,ij,kk)+dr*(pot(3,ij,kk) +dr*(pot(4,ij,kk))))))
                   phu = -1.0*(pot(2,ij,kk)+dr*(2.0*pot(3,ij,kk)+dr*(3.0*pot(4,ij,kk))))
                   atcf%fp(:,i)=atcf%fp(:,i)+0.5*phu*cvij(1,:)/rij
                   atcf%fp(:,j)=atcf%fp(:,j)-0.5*phu*cvij(1,:)/rij
                end if

                !              end if

             end if
          end do
       end do
       potisTersoff = potisTersoff + 0.5*v_ij
       if(lprteat.EQV..true.)then
          select type (atcf)
          class is (atom_config_e)
             atcf%eat(i) = atcf%eat(i)+atcf%eat(i)+0.5*v_ij
          end select
       end if


    end do

#ifdef PARA
    if (nprocspace.gt.1) then
       call comm_space%sum(potisTersoff)
       call comm_space%sum(sigT)
       if (associated(sigc)) then
          call comm_space%sum(sigc)
       endif
    end if
#endif
    sigcalfo=sigcalfo+sigT

#ifdef PARA
    if (nprocspace.gt.1) then
       write(6,*)'tersoff para ne fonctionne pas (envoi de "fp" non définis)'
!       call arret_ndm
       call maj_fp_frt(psc,atcf,celcf)
    end if
#endif
    !  write(6,*)sig
    !  write(6,*)

    return

  end subroutine force_tersoff_cel

end module force_tersoff_cel_mod
