module force_tersoff_cel_mod
        use cryst_to_cart_mod
        implicit none
        contains
! ***************************************************************
subroutine force_tersoff_cel
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use var_pot
  use tab_imm_m
  use jqmod
  use force_tersoff_facteurs
#ifdef PARA
  use mod_para
#endif
  ! **************************************************************
  ! Programme par NGUYEN Quoc Hoang
  ! CEA-Saclay/DEN/DMN/SRMP
  ! Printemps 2004 version du 27 septembre 2004

  ! Calcul des forces en utilisant la methode des cellules
  ! Test des contraintes -> OK (warning: convention de signe speciale)
  ! Probleme de flux
  ! ***************************************************************

  implicit none
  !-----------------------------------------------
  !   G l o b a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: i,j,k,nk,n_voisin,l,ij,ik, ipv,idv,ivj,ivk, m, moi
  integer , dimension(32) :: indice   ! Recense le nombre de voisins
  integer :: icelnumber,ncelvois,jcelvois,jcelnumber,jnumber, kcelvois,kcelnumber,knumber
  real(double) :: rij2, rij, rik2, rik, sui_ij, bij, n, v_ij, energie_i
  real(double) :: fc_rij, dfc_rij, fr_rij, fa_rij, fc_rik, dfc_rik, paire_ij, triplet_ij, triplet_ik 
  real(double) :: exponentiel, cos_theta, g_cos, dg_cos, flux, pression,ER1, ER2, ER3
  real(double) :: Scal_FiVi, Scal_FjVj, Scal_FijVi, Scal_FijVj, Scal_FikVi, Scal_FikVk
  real(double) , dimension(3) :: Xij, Xik, ai, Eflow, Flux1
  real(double) , dimension(15,6) :: tmp
  real(double) , dimension(15,3) :: tmp1
  real(double) , dimension(1,3) :: cvij, cvik
#ifdef PARA
  real(double) :: potist_tot, ER1_tot, ER2_tot, ER3_tot 
  real(double), dimension(3)   :: jq_tot
  real(double), dimension(3,3) :: sig_tot
  real(double) :: potisTersoff_tot
  real(double), dimension(3,3,noxyz) :: sigc_tot
#endif
  real(double):: coupR(npair)

  real(double) :: phu,sk,dr
  integer::kk


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
  !  jq(:)=0; sig(:,:)=0

  ER1=0. ;  ER2=0. ;  ER3=0.
  !         write(6,*)'boite quelc'
  call cryst_to_cart(imm,xp,bg,-1)

  !  write(6,*)sig
  !  write(6,*)
  do i=1,im
     if(typ_and_pot(ityp(i),ipotentiel).eqv..false.) cycle
     v_ij = 0

     icelnumber = ielat(i)
     ncelvois = min(noxyz,27)-1

     do jcelvois=0, ncelvois
        jcelnumber = ncel(icelnumber,jcelvois)
        do jnumber =1, nato(jcelnumber)
           j = last(jnumber,jcelnumber)
           if (j==i) then  !Cette condition n'est pas necessaire si JP. fais correctement sa table
              cycle
           else
              ij=ipo(ityp(i),ityp(j))
              !              write(6,*)ij,typ_pot_pair(ij),ipotentiel
              if (typ_pot_pair(ij).ne.ipotentiel) cycle
              do l=1,3
                 Xij(l)= xp(l,i)-xp(l,j)
              end do

              do l=1,3
                 !                        if (abs(Xij(l))>0.5) Xij(l)= Xij(l)-sign(1.0d0,xp(l,i))
                 if (Xij(l)>0.5) Xij(l)=Xij(l)-1.
                 if (Xij(l)<(-0.5)) Xij(l)=Xij(l)+1.
                 cvij(1,l) = Xij(l)
              end do
              call cryst_to_cart(1,cvij,at,1)
              rij2=cvij(1,1)**2+cvij(1,2)**2+cvij(1,3)**2
              if (rij2>Coupr(ij)**2) then         !borne sup de Lisa Porter 89
                 !if (rij2>(Rter(ij)+Coupr(ij))**2) then   !borne sup de Tersoff 88
                 cycle
              else
                 n=nter(ityp(i))
                 rij=sqrt(rij2)
                 sui_ij = 0
                 n_voisin = 0

                 do kcelvois=0, ncelvois
                    kcelnumber = ncel(icelnumber,kcelvois)
                    do knumber =1, nato(kcelnumber)
                       k = last(knumber,kcelnumber)

                       if (k==j .or. k==i) then
                          cycle
                       else
                          ik=ipo(ityp(i),ityp(k))
                          if (typ_pot_pair(ik).ne.ipotentiel) cycle
                          ! CONDITIONS PERIODIQUES
                          !write(6,*)'XP i j k',i,j,k
                          !write(6,*)xp(:,i)
                          !write(6,*)xp(:,j)
                          !write(6,*)xp(:,k)
                          do l=1,3
                             Xik(l)= xp(l,i)-xp(l,k)
                          end do
                          do l=1,3
                             if (XiK(l)>0.5) Xik(l)=Xik(l)-1.
                             if (XiK(l)<(-0.5)) Xik(l)=Xik(l)+1.
                             cvik(1,l) = Xik(l)
                          end do
                          !                       write(6,*)'cvik'
                          !                       write(6,*)cvik(1,:)
                          call cryst_to_cart(1,cvik,at,1)
                          rik2=cvik(1,1)**2+cvik(1,2)**2+cvik(1,3)**2
                          if (rik2>Coupr(ik)**2) then               !borne sup de Lisa Porter 89
                             !if (rik2>(Rter(ik)+Coupr(ik))**2) then     !borne sup de Tersoff 88
                             cycle
                          else
                             ! ne prend en compte que les k voisin de la paire ij
                             n_voisin = n_voisin + 1
                             indice(n_voisin) = k
                             rik = sqrt(rik2)
                             tmp(n_voisin,1) = rik                                                   ! rik
                             call facteur_amortissement (rik, ik, fc= tmp(n_voisin,2) )              ! fc(rik)
                             tmp(n_voisin,3) = exp(lambda3(ik)**3*(rij-rik)**3)                      ! exponentiel
                             tmp(n_voisin,4) = (cvij(1,1)*cvik(1,1)+cvij(1,2)*cvik(1,2)+cvij(1,3)*cvik(1,3))/(rij*rik)
                             ! cos theta ijk
                             call facteur_angulaire(tmp(n_voisin,4), ityp(i), g = tmp(n_voisin,5))   ! g(cos_theta)
                             tmp1(n_voisin,1)=cvik(1,1) ; tmp1(n_voisin,2)=cvik(1,2) ; tmp1(n_voisin,3)=cvik(1,3)

                             sui_ij = sui_ij + tmp(n_voisin,2)*tmp(n_voisin,5)*tmp(n_voisin,3)       ! sui ij

                          end if
                       end if
                    end do
                 end do

                 ! Tersoff standard                        bij = psi(ij)*(1.+beta(ityp(i))**n*sui_ij**n)**(-1./(2.*n)) 
                 !Brenner
                 bij = psi(ij)*(1.+beta(ityp(i))**n*sui_ij**n)**(-deltater(ityp(i)))

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

                       fp(l,i) = fp(l,i) + paire_ij
                       fp(l,j) = fp(l,j) - paire_ij

                       if (lcalcjq) then
                          Scal_FjVj=Scal_FjVj - paire_ij*vp(l,j)
                       end if
                       !Contrainte
                       if (mod(it,itesigma)==0) then 
                          do m=1,3
                             sig(l,m)=sig(l,m) + paire_ij*cvij(1,m)/volu
                             if (lTPcel.EQV..true.) then
                                sigc(l,m,icelnumber) = sigc(l,m,icelnumber) + 0.5*paire_ij*cvij(1,m)*noxyz/volu
                                sigc(l,m,jcelnumber) = sigc(l,m,jcelnumber) + 0.5*paire_ij*cvij(1,m)*noxyz/volu
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

                       fp(l,i) = fp(l,i) + paire_ij
                       fp(l,j) = fp(l,j) - paire_ij

                       if (lcalcjq) then
                          Scal_FjVj=Scal_FjVj - paire_ij*vp(l,j)
                       end if
                       !Contrainte
                       if (mod(it,itesigma)==0) then 
                          do m=1,3
                             sig(l,m)=sig(l,m) + paire_ij*cvij(1,m)/volu
                             if (lTPcel.EQV..true.) then
                                sigc(l,m,icelnumber) = sigc(l,m,icelnumber) + paire_ij*cvij(1,m)*noxyz/volu
                                sigc(l,m,jcelnumber) = sigc(l,m,jcelnumber) + paire_ij*cvij(1,m)*noxyz/volu
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
                       ik=ipo(ityp(i),ityp(k))
                       call facteur_amortissement (rik, ik, dfc= dfc_rik)
                       exponentiel = tmp(nk,3)
                       cos_theta = tmp(nk,4)
                       g_cos = tmp(nk,5)
                       call facteur_angulaire (cos_theta, ityp(i), dg= dg_cos)
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


                          triplet_ij = 0.5*fc_rij*fa_rij*deltater(ityp(i))*n*psi(ij)*beta(ityp(i))**n*sui_ij**(n-1)*&
                               (1.+beta(ityp(i))**n*sui_ij**n)**(-deltater(ityp(i))-1.) * &   
                               (3.*lambda3(ik)**3*fc_rik*(rij-rik)**2*g_cos*exponentiel/rij*cvij(1,l) &
                               +fc_rik*dg_cos*exponentiel* &
                               (cvik(1,l)/(rij*rik)-cos_theta*cvij(1,l)/rij**2))


                          !                                 triplet_ik = 0.5*fc_rij*fa_rij*0.5*psi(ij)*beta(ityp(i))**n*sui_ij**(n-1)*&
                          !                                      (1.+beta(ityp(i))**n*sui_ij**n)**(-(1.+2.*n)/(2.*n)) * &
                          !                                      (-3.*lambda3(ik)**3*fc_rik*(rij-rik)**2*g_cos*exponentiel/rik*cvik(1,l) &
                          !                                      +dfc_rik*g_cos*exponentiel/rik*cvik(1,l) &
                          !                                      +fc_rik*dg_cos*exponentiel* &
                          !                                      (cvij(1,l)/(rij*rik)-cos_theta*cvik(1,l)/rik**2))

                          triplet_ik = 0.5*fc_rij*fa_rij*deltater(ityp(i))*n*psi(ij)*beta(ityp(i))**n*sui_ij**(n-1)*&
                               (1.+beta(ityp(i))**n*sui_ij**n)**(-deltater(ityp(i))-1.)* &
                               (-3.*lambda3(ik)**3*fc_rik*(rij-rik)**2*g_cos*exponentiel/rik*cvik(1,l) &
                               +dfc_rik*g_cos*exponentiel/rik*cvik(1,l) &
                               +fc_rik*dg_cos*exponentiel* &
                               (cvij(1,l)/(rij*rik)-cos_theta*cvik(1,l)/rik**2))

                          fp(l,i) = fp(l,i) + triplet_ij
                          fp(l,j) = fp(l,j) - triplet_ij
                          fp(l,i) = fp(l,i) + triplet_ik
                          fp(l,k) = fp(l,k) - triplet_ik

                          !Flux
                          if (lcalcjq) then
                             Scal_FijVj=Scal_FijVj - triplet_ij*vp(l,j)
                             Scal_FikVk=Scal_FikVk - triplet_ik*vp(l,k)
                          end if
                          !Contrainte
                          !write(6,*)'sig AV l',l
                          !write(6,*)sig
                          if (mod(it,itesigma)==0) then 
                             do m=1,3
                                sig(l,m)=sig(l,m) + triplet_ij*cvij(1,m)/volu
                                sig(l,m)=sig(l,m) + triplet_ik*cvik(1,m)/volu
                             if (lTPcel.EQV..true.) then
                                sigc(l,m,icelnumber) = sigc(l,m,icelnumber) + noxyz*0.5*(triplet_ij*cvij(1,m)/volu  &
                                                       + triplet_ik*cvik(1,m)/volu)
                                sigc(l,m,jcelnumber) = sigc(l,m,jcelnumber) + noxyz*0.5*triplet_ij*cvij(1,m)/volu 
                                sigc(l,m,kcelnumber) = sigc(l,m,kcelnumber) + noxyz*0.5*triplet_ik*cvik(1,m)/volu
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
                    fp(:,i)=fp(:,i)+0.5*phu*cvij(1,:)/rij
                    fp(:,j)=fp(:,j)-0.5*phu*cvij(1,:)/rij
                 end if

              end if

           end if
        end do
     end do
     if (associated (free)) then
        if (free(i).EQV..true.)potisTersoff = potisTersoff + 0.5*v_ij
     else
        potisTersoff = potisTersoff + 0.5*v_ij
     end if
     if (associated (free)) then
        if ((associated(eatom)).and.(free(i).EQV..true.)) eatom(i) = eatom(i)+eatom(i)+0.5*v_ij
     else
        if (associated(eatom)) eatom(i) = eatom(i)+eatom(i)+0.5*v_ij
     end if


  end do
  call cryst_to_cart(imm,xp,at,1)

#ifdef PARA
  call MPI_ALLREDUCE(potisTersoff,potisTersoff_tot,1,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
  potisTersoff=potisTersoff_tot
  !     call MPI_ALLREDUCE(jq,    jq_tot,    3,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
  !     jq=jq_tot
  call MPI_ALLREDUCE(sig,   sig_tot,   9,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
  sig=sig_tot  
  if (associated(sigc)) then
     call MPI_ALLREDUCE(sigc,      sigc_tot,      9*noxyz,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
     sigc=sigc_tot
  endif

#endif


#ifdef PARA
  call maj_fp_frt
#endif
  !  write(6,*)sig
  !  write(6,*)

  return

end subroutine force_tersoff_cel

end module
