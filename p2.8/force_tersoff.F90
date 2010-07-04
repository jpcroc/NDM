! ***************************************************************
subroutine force_tersoff (xp,  vp,  fp,  iwmax, ityp)
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use jqmod
  use force_tersoff_facteurs
  ! **************************************************************
  ! Programme par NGUYEN Quoc Hoang, J-P Crocombette
  ! CEA-Saclay/DEN/DMN/SRMP
  ! Printemps 2004 version du 27 septembre 2004

  ! Calcul des forces en utilisant la table des voisins
  ! ***************************************************************

  implicit none
  !-----------------------------------------------
  !   G l o b a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  integer :: iwmax(imm)
  integer , intent(in) :: ityp(imm)
  real(double) , intent(inout) :: xp(3,imm)
  real(double) :: vp(3,imm)
  real(double) , intent(inout) :: fp(3,imm)
  !-----------------------------------------------
  !   L o c a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: i,j,k,nk,n_voisin,l,ij,ik, ipv,idv,ivj,ivk, m, moi,nvij
  integer , dimension(32) :: indice   ! Recense le nombre de voisins
  real(double) :: rij2, rij, rik2, rik, sui_ij, bij, n, v_ij, energie_i
  real(double) :: fc_rij, dfc_rij, fr_rij, fa_rij, fc_rik, dfc_rik, paire_ij, triplet_ij, triplet_ik 
  real(double) :: exponentiel, cos_theta, g_cos, dg_cos, flux, pression,ER1, ER2, ER3
  real(double) :: Scal_FiVi, Scal_FjVj, Scal_FijVi, Scal_FijVj, Scal_FikVi, Scal_FikVk
  real(double) :: XijdotF, XikdotF,fpnemd(3),fpnemdmoy(3)
  real(double) , dimension(3) :: Xij, Xik, ai
  real(double) , dimension(15,6) :: tmp
  real(double) , dimension(15,3) :: tmp1
  real(double) , dimension(1,3) :: cvij, cvik
  logical::test_sigma
  real(double):: coupR(npair)
  real(double) :: phu,sk,dr
  integer::kk
  !-----------------------------------------------
  !INIALISATION

#ifdef paraTersoff



  ! declarations supplementaires pour MPI
  real(double), dimension (3,imm) :: fpTemp
  real(double) ::  potistTemp
  real(double) :: sigTemp (3,3)
  !  real(double) :: eatomtemp(imm)
  integer :: Fin,Deb
  real(double) :: div
#endif

  select case (ipotentiel) 
  case(13)
     coupR(:)=ster(:)
  case(14,15)
     coupR(:)=rter(:)+ster(:)
  end select
  !write(6,*)'coupR',coupR

  if (mod(it,itesigma)==0) then
     test_sigma=.true.
  else
     test_sigma=.false.
  end if

!  potist = 0
  moi =0
!  do i=1,im
!     fp(1,i)=0; fp(2,i)=0; fp(3,i)=0
!  end do
!  jq(:)=0; sig(:,:)=0

  fpnemdmoy(:)=0
  ER1=0. ;  ER2=0. ;  ER3=0.
  call cryst_to_cart(imm,xp,bg,-1)

  idv = 0

#ifdef paraTersoff
  ! MPI

  deb=1+Int(im*sqrt(float(rang)/nb_procs))
  fin=Int(im*sqrt(float(rang+1)/nb_procs))
  if (it==1) write(6,*)'rang deb Fin nb d_at  ',rang,deb,fin,fin-deb+1
  Tloop1at1:  do i=Deb,Fin


#else
     ! Sequentiel
     Tloop1at1:  do i=1,im
#endif
! Sequentiel

        v_ij = 0

        ipv = idv+1
        idv = iwmax(i)
        Tloop1at2 :do ivj=ipv, idv
           j = indi (ivj)
           !               write(6,*)i,j
           if (j==i) then  !Cette condition n'est pas necessaire si JP. fais correctement sa table
              cycle
           else
              ij=ipo(ityp(i),ityp(j))

              do l=1,3
                 Xij(l)= xp(l,i)-xp(l,j)
              end do

              do l=1,3

                 if (Xij(l)>0.5) Xij(l)=Xij(l)-1.
                 if (Xij(l)<(-0.5)) Xij(l)=Xij(l)+1.
                 cvij(1,l) = Xij(l)
              end do
              call cryst_to_cart(1,cvij,at,1)
              rij2=cvij(1,1)**2+cvij(1,2)**2+cvij(1,3)**2
              !                  write(6,*)i,j,rij2
              if (rij2>CoupR(ij)**2) then         !borne sup de Lisa Porter 89
                 !if (rij2>(Rter(ij)+CoupR(ij))**2) then   !borne sup de Tersoff 88
                 cycle
              else
                 n=nter(ityp(i))
                 rij=sqrt(rij2)
                 sui_ij = 0
                 n_voisin = 0
                 Tloop1at3:    do ivk=ipv, idv
                    k = indi (ivk)
                    if (k==j .or. k==i) then
                       cycle
                    else
                       ik=ipo(ityp(i),ityp(k))

                       ! CONDITIONS PERIODIQUES
                       do l=1,3
                          Xik(l)= xp(l,i)-xp(l,k)
                       end do
                          !write(6,*)'XP i j k',i,j,k
                          !write(6,*)xp(:,i)
                          !write(6,*)xp(:,j)
                          !write(6,*)xp(:,k)
 
                       do l=1,3
                          if (XiK(l)>0.5) Xik(l)=Xik(l)-1.
                          if (XiK(l)<(-0.5)) Xik(l)=Xik(l)+1.
                          cvik(1,l) = Xik(l)
                       end do
!                       write(6,*)'cvik'
!                       write(6,*)cvik(1,:)
                       call cryst_to_cart(1,cvik,at,1)
                       rik2=cvik(1,1)**2+cvik(1,2)**2+cvik(1,3)**2
                       if (rik2>CoupR(ik)**2) then               !borne sup de Lisa Porter 89
                          !if (rik2>(Rter(ik)+CoupR(ik))**2) then     !borne sup de Tersoff 88
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
                 end do Tloop1at3
                 !
                 !Tersoff OR                     bij = psi(ij)*(1.+beta(ityp(i))**n*sui_ij**n)**(-1./(2.*n)) 
                 ! Brenner
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
                       if (test_sigma.EQV..true.) then 
                          do m=1,3
                             sig(l,m)=sig(l,m) + paire_ij*cvij(1,m)/volu
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
                       if (test_sigma.EQV..true.) then  
                          do m=1,3
                             sig(l,m)=sig(l,m) + paire_ij*cvij(1,m)/volu
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
                    Tloop2at3:        do nk=1,n_voisin
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

                          !                              triplet_ij = 0.5*fc_rij*fa_rij*0.5*psi(ij)*beta(ityp(i))**n*sui_ij**(n-1)*&
                          !                                   (1.+beta(ityp(i))**n*sui_ij**n)**(-(1.+2.*n)/(2.*n)) * &   
                          !                                   (3.*lambda3(ik)**3*fc_rik*(rij-rik)**2*g_cos*exponentiel/rij*cvij(1,l) &
                          !                                   +fc_rik*dg_cos*exponentiel* &
                          !                                   (cvik(1,l)/(rij*rik)-cos_theta*cvij(1,l)/rij**2))
                          !                              
                          !                              
                          !                              triplet_ik = 0.5*fc_rij*fa_rij*0.5*psi(ij)*beta(ityp(i))**n*sui_ij**(n-1)*&
                          !                                   (1.+beta(ityp(i))**n*sui_ij**n)**(-(1.+2.*n)/(2.*n)) * &
                          !                                   (-3.*lambda3(ik)**3*fc_rik*(rij-rik)**2*g_cos*exponentiel/rik*cvik(1,l) &
                          !                                   +dfc_rik*g_cos*exponentiel/rik*cvik(1,l) &
                          !                                   +fc_rik*dg_cos*exponentiel* &
                          !                                   (cvij(1,l)/(rij*rik)-cos_theta*cvik(1,l)/rik**2))

                          triplet_ij = 0.5*fc_rij*fa_rij*deltater(ityp(i))*n*psi(ij)*beta(ityp(i))**n*sui_ij**(n-1)*&
                               (1.+beta(ityp(i))**n*sui_ij**n)**(-deltater(ityp(i))-1.) * &   
                               (3.*lambda3(ik)**3*fc_rik*(rij-rik)**2*g_cos*exponentiel/rij*cvij(1,l) &
                               +fc_rik*dg_cos*exponentiel* &
                               (cvik(1,l)/(rij*rik)-cos_theta*cvij(1,l)/rij**2))


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
                          if (test_sigma.EQV..true.) then  
                             do m=1,3
                                sig(l,m)=sig(l,m) + triplet_ij*cvij(1,m)/volu
                                sig(l,m)=sig(l,m) + triplet_ik*cvik(1,m)/volu
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
                    end do Tloop2at3
                 end if
!                 write(6,*)'rij roff',rij,roff2(ij)
!                 if (rij.le.roff2(ij)) then
!                       sk = rij/csive
!                       kk = sk
!                       ! spline
!                       dr = rij-float(kk)*csive
!!                       write(6,*)'rij roff',rij,roff2(ij),pot(1,l,kk)
!                       if (free(i)==.true)potist = potist+pot(1,l,kk)+ rij*(dr*(pot(2,l,kk)+dr*(pot(3,l,kk) +dr*(pot(4,l,kk)))))                       
!                       phu = -1.0*(pot(2,l,kk)+dr*(2.0*pot(3,l,kk)+dr*(3.0*pot(4,l,kk))))
!                       fp(:,i)=fp(:,i)+phu*cvij(1,:)
!                       fp(:,j)=fp(:,j)-phu*cvij(1,:)
!
!                 end if

              end if
           end if

        end do Tloop1at2
        if (associated (free)) then
           if (free(i).EQV..true.)potist = potist + 0.5*v_ij
        else
           potist = potist + 0.5*v_ij
        end if


        !energie_i =  0.5*v_i
        if (associated (free)) then
           if ((lprteat.or.lcalcjq.or.lnemd).and.(free(i).EQV..true.))eatom(i) = eatom(i)+0.5*v_ij
        else
           if (lprteat.or.lcalcjq.or.lnemd)eatom(i) = eatom(i)+0.5*v_ij
        end if

     end do Tloop1at1


#ifdef paraTersoff
     CALL MPI_BARRIER(MPI_COMM_WORLD,code)

     CALL MPI_ALLREDUCE(fp,fpTemp,3*im,MPI_DOUBLE_PRECISION,&
          MPI_SUM,MPI_COMM_WORLD,code)
     !      fp(:,1:im)=fpTemp(:,1:im)
     fp=fpTemp

     CALL MPI_ALLREDUCE(potist,potistTemp,1,MPI_DOUBLE_PRECISION,&
          MPI_SUM,MPI_COMM_WORLD,code)

     potist=potistTemp

     ! MPI : collecte generale et somme des contraintes calcules par les process
     if(test_sigma)then
        CALL MPI_ALLREDUCE(sig,sigTemp,9,MPI_DOUBLE_PRECISION,&
             MPI_SUM,MPI_COMM_WORLD,code)
        sig=sigTemp
     end if

     !  if (lcalcjq) then
     !     CALL MPI_ALLREDUCE(jq,jqTemp,3,MPI_DOUBLE_PRECISION,&
     !          MPI_SUM,MPI_COMM_WORLD,code)
     !     jq=jqTemp
     !     CALL MPI_ALLREDUCE(eatom,eatomTemp,imm,MPI_DOUBLE_PRECISION,&
     !          MPI_SUM,MPI_COMM_WORLD,code)
     !     eatom=eatomTemp
     !  end if

#endif



     call cryst_to_cart(imm,xp,at,1)



     return

   end subroutine force_tersoff


