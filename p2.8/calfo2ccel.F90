module calfo2ccel_mod
  USE notperiod_mod,only: notperiod
  USE cryst_to_cart_mod,only: cryst_to_cart
  USE var_pot, ONLY:alpha,csive,ipotentiel,ipo,zz,ipo,rue_pair,pot,typ_and_pot,typ_pot_pair
  USE calfocommon
  implicit none
contains
  ! ***************************************************************
  subroutine calfo2ccel(im,imm,xp,   fp,  ityp,ielat,num_at_glob,noxyz,natperc,atincel,nato,ncel,deltadist,at,bg,volu)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double
    USE gen_com_m , ONLY:lcalcjq,lperiod,pi,potis1,potis2
    USE jqmod
#ifdef PARA
    use mpi
    USE Tpara,only:MPI_COMM_space,NDM_MPI_REAL_DOUBLE,nprocspace
#else
    USE Tpara,only:nprocspace
#endif
    implicit none
    !-----------------------------------------------
    !   G l o b a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer,intent(in)::im,imm
    integer , intent(in),allocatable :: ielat(:),ityp(:),num_at_glob(:)
    real(double),intent(in),allocatable  :: xp(:,:)
    real(double) , intent(inout),allocatable :: fp(:,:)

    integer,intent(in)::noxyz,natperc
    integer, intent(in), allocatable::nato(:),ncel(:,:),atincel(:,:),deltadist(:,:,:)
    real(double),intent(in),dimension(3,3)::at,bg
    real(double),intent(in)::volu
    !-----------------------------------------------
    !   L o c a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: iti, l, i, koo, i1, ko1, j, i2, itj, k, &
         ic, ncelvois,itimin,itimax
    real(double) :: aux, alp, f1, f2, f3,  c1, c2&
         , c3, c1p,c2p,c3p, sk, r, phu, c1abs,c2abs,c3abs, ra(3),cv(1,3)
    real(double) :: dr,deltaepot,fcontr
#ifdef PARA
    real(double) :: potis1_tot, potis2_tot
    real(double) :: deltaF_tot,deltaEpot_tot,deltaEspr_tot,Espr_tot,deltafcomp
    real(double), dimension(3,3) :: sig_tot
    real(double), dimension(3,3,noxyz) :: sigc_tot
#endif

    !-----------------------------------------------
    !
    !
    ! *** Initialisations ***
    real(double),allocatable :: xpnp(:,:)

    allocate(xpnp(3,imm))
    ! Initialisation des termes du potentiel
    !  potis2 = zero
    !  potis0 = zero
    !  potis3 = zero
    !  potis1 = zero
    !  potist = zero

    if (lperiod) then
       xpnp(:,:)=xp(:,:)
    else 
       call notperiod(imm,xp,xpnp,at,bg)
    end if


    ! Declarations de constantes
    aux = 23.06134575D-20
    alp = alpha/sqrt(pi)*aux

    ! Initialisation des contraintes pour le systeme global et les process

    ! terme de paires

    !pour chaque atome





    do i = 1, im
       if(typ_and_pot(ityp(i),ipotentiel).eqv..false.) cycle
       koo = ielat(i)                          ! Numero de la cellule
       iti = ityp(i)
       ! --- Calcul du second potentiel de la somme d'Ewald ---
       l = ipo(iti,iti)

       potis2 = potis2-zz(l)*alp

       ncelvois = min(noxyz,27)-1

       ! pour chaque cel. voisine
       do i1 = 0, ncelvois

          ko1 = ncel(koo,i1)
          c1p = xpnp(1,i)+sum(at(1,:)*deltadist(:,i1,koo))
          c2p = xpnp(2,i)+sum(at(2,:)*deltadist(:,i1,koo))
          c3p = xpnp(3,i)+sum(at(3,:)*deltadist(:,i1,koo))

          ! pour chaque atome ds la cel. voisine
          do i2 = 1, nato(ko1)
             j = atincel(i2,ko1)
             itj = ityp(j)
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
                if (j.le.im) then
                   ! les deux atomes sont locaux
                   if (num_at_glob(i).ge.num_at_glob(j)) cycle !terme déja calculé
                else
                   ! j n'est pas local, on fait le calcul normal           
                endif
             else
                if (num_at_glob(i).ge.num_at_glob(j)) cycle !terme déja calculé
             end if
             
!#else
!             if (num_at_glob(i).ge.num_at_glob(j)) cycle !terme déja calculé
!#endif

             if (noxyz.ne.1) then
                c1 = c1p-xpnp(1,j)
                c1abs=abs(c1)
                if(c1abs>rue_pair(l)) cycle

                c2 = c2p-xpnp(2,j)
                c2abs=abs(c2)
                if(c2abs>rue_pair(l)) cycle
                c3 = c3p-xpnp(3,j)
                c3abs=abs(c3)
                if(c3abs>rue_pair(l)) cycle
                cv(1,1) = c1
                cv(1,2) = c2
                cv(1,3) = c3

             else
                c1 = c1p-xpnp(1,j)
                c2 = c2p-xpnp(2,j)
                c3 = c3p-xpnp(3,j)
                cv(1,1) = c1
                cv(1,2) = c2
                cv(1,3) = c3
                call cryst_to_cart (1, cv, bg, -1) !cart vers cryst sur cv
                WHERE ( (cv.GT.0.5d0).OR.(cv.LT.-0.5d0) )
                   cv(:,1:3) = cv(:,1:3) - Dble(Nint(cv(:,1:3)))
                END WHERE
                call cryst_to_cart (1, cv, at, 1) !cryst vers cart sur cv
                c1=cv(1,1)
                c2=cv(1,2)
                c3=cv(1,3)

             end if


             r = c1*c1+c2*c2+c3*c3
             if (r>rue_pair(l)**2) cycle
             r=sqrt(r)
             sk = r/csive
             k = sk
             ! spline
             dr = r-float(k)*csive
             phu = -1.0*(pot(2,l,k)+(2.0*pot(3,l,k)+3.0*pot(4,l,k)*dr)*dr)
             deltaepot=0.5*(pot(1,l,k)+(pot(2,l,k)*dr+pot(3,l,k)*dr**2+pot(4,l,k)*dr**3)*r)




             f1 = phu*c1
             f2 = phu*c2
             f3 = phu*c3
             ra(1)=f1 ; ra(2)=f2 ; ra(3)=f3
             fcontr=sqrt(f1**2+f2**2+f3**2)
             fp(1,i) = fp(1,i)+f1
             fp(2,i) = fp(2,i)+f2
             fp(3,i) = fp(3,i)+f3
             potis1 = potis1+deltaepot
#ifdef PARA
             if (j.le.im) then
#endif
                potis1 = potis1+deltaepot
#ifdef PARA
             endif
#endif


             !           endif


             fp(1,j) = fp(1,j)-f1
             fp(2,j) = fp(2,j)-f2
             fp(3,j) = fp(3,j)-f3
!!$
!!$             if (lcalcjq) then
!!$                jqf=0.0
!!$                eat(i) = eat(i)+deltaepot
!!$#ifdef PARA
!!$                if (nprocspace.gt.1) then
!!$                   if (j.le.im) then
!!$                      eat(j) = eat(j)+deltaepot
!!$                      do ic=1,3
!!$                         jqf=jqf-0.5*(ra(ic)*(vp(ic,i)+vp(ic,j)))
!!$                      end do
!!$                   else
!!$                      do ic=1,3
!!$                         jqf=jqf-0.5*(ra(ic)*(vp(ic,i)))
!!$                      end do
!!$                   end if
!!$                end if
!!$#else
!!$                eat(j) = eat(j)+deltaepot
!!$                do ic=1,3
!!$                   jqf=jqf-0.5*(ra(ic)*(vp(ic,i)+vp(ic,j)))
!!$                end do
!!$                do ic=1,3
!!$                   jq(ic)=jq(ic)-jqf*cv(1,ic)
!!$                end do
!!$#endif
!!$             end if

             if (lprteat) then
                !              if (allocated (free)) then
                !                 if (free(i).EQV..true.)eat(i) = eat(i)+deltaepot
                !                 if (free(j).EQV..true.)eat(j) = eat(j)+deltaepot
                !              else
                eat(i) = eat(i)+deltaepot
                eat(j) = eat(j)+deltaepot
                !              end if
             end if

             ! calcul des contraintes

             if (test_sigma) then
                if (num_at_glob(i).lt.num_at_glob(j)) then
                   sig(1,1) = sig(1,1)+phu*c1*c1/volu
                   sig(1,2) = sig(1,2)+phu*c1*c2/volu
                   sig(1,3) = sig(1,3)+phu*c1*c3/volu
                   sig(2,1) = sig(2,1)+phu*c2*c1/volu
                   sig(2,2) = sig(2,2)+phu*c2*c2/volu
                   sig(2,3) = sig(2,3)+phu*c2*c3/volu
                   sig(3,1) = sig(3,1)+phu*c3*c1/volu
                   sig(3,2) = sig(3,2)+phu*c3*c2/volu
                   sig(3,3) = sig(3,3)+phu*c3*c3/volu
                endif
                if (lTPcel.EQV..true.) then
                   sigc(1,1,koo) = sigc(1,1,koo)+0.5*phu*c1*c1*noxyz/volu
                   sigc(1,2,koo) = sigc(1,2,koo)+0.5*phu*c1*c2*noxyz/volu
                   sigc(1,3,koo) = sigc(1,3,koo)+0.5*phu*c1*c3*noxyz/volu
                   sigc(2,1,koo) = sigc(2,1,koo)+0.5*phu*c2*c1*noxyz/volu
                   sigc(2,2,koo) = sigc(2,2,koo)+0.5*phu*c2*c2*noxyz/volu
                   sigc(2,3,koo) = sigc(2,3,koo)+0.5*phu*c2*c3*noxyz/volu
                   sigc(3,1,koo) = sigc(3,1,koo)+0.5*phu*c3*c1*noxyz/volu
                   sigc(3,2,koo) = sigc(3,2,koo)+0.5*phu*c3*c2*noxyz/volu
                   sigc(3,3,koo) = sigc(3,3,koo)+0.5*phu*c3*c3*noxyz/volu
                   sigc(1,1,ko1) = sigc(1,1,ko1)+0.5*phu*c1*c1*noxyz/volu
                   sigc(1,2,ko1) = sigc(1,2,ko1)+0.5*phu*c1*c2*noxyz/volu
                   sigc(1,3,ko1) = sigc(1,3,ko1)+0.5*phu*c1*c3*noxyz/volu
                   sigc(2,1,ko1) = sigc(2,1,ko1)+0.5*phu*c2*c1*noxyz/volu
                   sigc(2,2,ko1) = sigc(2,2,ko1)+0.5*phu*c2*c2*noxyz/volu
                   sigc(2,3,ko1) = sigc(2,3,ko1)+0.5*phu*c2*c3*noxyz/volu
                   sigc(3,1,ko1) = sigc(3,1,ko1)+0.5*phu*c3*c1*noxyz/volu
                   sigc(3,2,ko1) = sigc(3,2,ko1)+0.5*phu*c3*c2*noxyz/volu
                   sigc(3,3,ko1) = sigc(3,3,ko1)+0.5*phu*c3*c3*noxyz/volu
                end if
             endif

          end do  ! fin i2=j
       end do ! fin i1=koo


    end do ! fin i

#ifdef PARA

    if (nprocspace.gt.1) then

       call MPI_ALLREDUCE(potis1,potis1_tot,1,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_space,ierr)
       potis1=potis1_tot
       call MPI_ALLREDUCE(potis2,potis2_tot,1,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_space,ierr)
       potis2=potis2_tot
       call MPI_ALLREDUCE(sig,sig_tot,9,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_space,ierr)
       sig=sig_tot
       if (associated(sigc)) then
          call MPI_ALLREDUCE(sigc,      sigc_tot,      9*noxyz,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_space,ierr)
          sigc=sigc_tot
       endif
    end if

#endif


    ! fin du calcul du terme de paire dans l'espace direct
    deallocate ( xpnp)
    return
  end subroutine calfo2ccel

end module calfo2ccel_mod
