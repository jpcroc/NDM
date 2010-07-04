! ***************************************************************
subroutine calfo2ccel
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use jqmod
  use tab_imm_m
#if(PARA)
  use mod_mpi
#endif
  !...Translated by PSUITE Trans90                  4.3ZH 16:03:53   7/03/ 1
  !...Switches: -nqp -rl -xf -xhm -x
  !      version para-seq du 22 fevrier 2001
  ! **************************************************************

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
  integer :: iti, l, i, koo, i1, ko1, j, i2, itj, k, &
       ic, ncelvois,itimin,itimax
  real(double) :: aux, alp, f1, f2, f3,  c1, c2&
       , c3, c1p,c2p,c3p, sk, r, phu, c1abs,c2abs,c3abs, ra(3),cv(1,3)
  real(double) :: dr,deltaepot,fcontr
#if(PARA)
  real(double) :: potis1_tot, potis2_tot
  real(double) :: deltaF_tot,deltaEpot_tot,deltaEspr_tot,Espr_tot,deltafcomp
  real(double), dimension(3,3) :: sig_tot
#endif

  !-----------------------------------------------
  !
  !
  ! *** Initialisations ***

  real(double),pointer :: xpnp(:,:)
    allocate(xpnp(3,imm))
  ! Initialisation des termes du potentiel
  !  potis2 = zero
  !  potis0 = zero
  !  potis3 = zero
  !  potis1 = zero
  !  potist = zero

  ! write(6,*)'PARA-T R I I ',rang,im,imm
  if (lperiod) then
     xpnp(:,:)=xp(:,:)
  else 
     call notperiod(xp,xpnp)
  end if


  ! Declarations de constantes
  aux = 23.06134575D-20
  alp = alpha/sqrt(pi)*aux

  ! Initialisation des contraintes pour le systeme global et les process

  ! terme de paires

  !pour chaque atome

           if (ldesinteg) then
#if(PARA)
	deltafcomp=0.
#endif
              if (pm1des==-1) then
                 lambdades=1.-float(itdes)/float(nstepdes)           
                else
                 lambdades=float(itdes)/float(nstepdes)           
              end if
           endif





  do i = 1, im
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
           j = last(i2,ko1)
           itj = ityp(j)
           l = ipo(iti,itj)
           if (typ_pot_pair(l).ne.ipotentiel) cycle

#if(PARA)
	   ! Methode pour ne prendre qu'une seule fois en compte
           ! le couple i,j en paralle :
	   ! - i est necesairement local (boucle i<=im)
           ! - si j est local on ne retient que le couple i<j
           ! - si j n'est pas local, le couple n'est par definition
	   !   pris qu'une fois puisque i est local
	   if (j.le.im) then
       ! les deux atomes sont locaux
              if (num_at_glob(i).ge.num_at_glob(j)) cycle !terme déja calculé
	   else
       ! j n'est pas local, on fait le calcul normal           
	   endif
#else
           if (num_at_glob(i).ge.num_at_glob(j)) cycle !terme déja calculé
#endif

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



           if ((ldesinteg).and.(num_at_glob(i)==1)) then

#if(PARA)
	deltafcomp=deltafcomp-deltaepot*1./float(nstepdes)
#else
                 deltaF=deltaF-deltaepot*1./float(nstepdes)
#endif
                 phu=phu*lambdades
                 deltaepot=deltaepot*lambdades
!#if(PARA)
!	write(6,*)'des',rang,i,deltafcomp
!#endif
           endif




           f1 = phu*c1
           f2 = phu*c2
           f3 = phu*c3
           ra(1)=f1 ; ra(2)=f2 ; ra(3)=f3
           fcontr=sqrt(f1**2+f2**2+f3**2)
           fp(1,i) = fp(1,i)+f1
           fp(2,i) = fp(2,i)+f2
           fp(3,i) = fp(3,i)+f3




              if (associated (free)) then
                 if (free(i).EQV..true.)potis1 = potis1+deltaepot
#if(PARA)
                 if (j.le.im) then
#endif
                    if (free(j).EQV..true.)potis1 = potis1+deltaepot
#if(PARA)
                 endif
#endif

              else
                 potis1 = potis1+deltaepot
#if(PARA)
                 if (j.le.im) then
#endif
                    potis1 = potis1+deltaepot
#if(PARA)
                 endif
#endif

                 
              endif


           fp(1,j) = fp(1,j)-f1
           fp(2,j) = fp(2,j)-f2
           fp(3,j) = fp(3,j)-f3

           if (lcalcjq) then
              jqf=0.0
              eatom(i) = eatom(i)+deltaepot
#if(PARA)
              if (j.le.im) then
                 eatom(j) = eatom(j)+deltaepot
                 do ic=1,3
                    jqf=jqf-0.5*(ra(ic)*(vp(ic,i)+vp(ic,j)))
                 end do
              else
                 do ic=1,3
                    jqf=jqf-0.5*(ra(ic)*(vp(ic,i)))
                 end do
              end if
#else
              eatom(j) = eatom(j)+deltaepot
              do ic=1,3
                 jqf=jqf-0.5*(ra(ic)*(vp(ic,i)+vp(ic,j)))
              end do
              do ic=1,3
                 jq(ic)=jq(ic)-jqf*cv(1,ic)
              end do
#endif
           end if

           if (lprteat) then
              if (associated (free)) then
                 if (free(i).EQV..true.)eatom(i) = eatom(i)+deltaepot
                 if (free(j).EQV..true.)eatom(j) = eatom(j)+deltaepot
              else
                 eatom(i) = eatom(i)+deltaepot
                 eatom(j) = eatom(j)+deltaepot
              end if
           end if

           ! calcul des contraintes

           if (itesigma>0) then
              if (mod(it,itesigma)==0) then

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
                    sigc(1,1,koo) = sigc(1,1,koo)+phu*c1*c1*noxyz/volu
                    sigc(1,2,koo) = sigc(1,2,koo)+phu*c1*c2*noxyz/volu
                    sigc(1,3,koo) = sigc(1,3,koo)+phu*c1*c3*noxyz/volu
                    sigc(2,1,koo) = sigc(2,1,koo)+phu*c2*c1*noxyz/volu
                    sigc(2,2,koo) = sigc(2,2,koo)+phu*c2*c2*noxyz/volu
                    sigc(2,3,koo) = sigc(2,3,koo)+phu*c2*c3*noxyz/volu
                    sigc(3,1,koo) = sigc(3,1,koo)+phu*c3*c1*noxyz/volu
                    sigc(3,2,koo) = sigc(3,2,koo)+phu*c3*c2*noxyz/volu
                    sigc(3,3,koo) = sigc(3,3,koo)+phu*c3*c3*noxyz/volu
                 end if
                 if (lsigtyp.EQV..true.) then
                    sigtyp(1,1,ityp(i)) = sigtyp(1,1,ityp(i))+phu*c1*c1*0.5
                    sigtyp(1,2,ityp(i)) = sigtyp(1,2,ityp(i))+phu*c1*c2*0.5
                    sigtyp(1,3,ityp(i)) = sigtyp(1,3,ityp(i))+phu*c1*c3*0.5
                    sigtyp(2,1,ityp(i)) = sigtyp(2,1,ityp(i))+phu*c2*c1*0.5
                    sigtyp(2,2,ityp(i)) = sigtyp(2,2,ityp(i))+phu*c2*c2*0.5
                    sigtyp(2,3,ityp(i)) = sigtyp(2,3,ityp(i))+phu*c2*c3*0.5
                    sigtyp(3,1,ityp(i)) = sigtyp(3,1,ityp(i))+phu*c3*c1*0.5
                    sigtyp(3,2,ityp(i)) = sigtyp(3,2,ityp(i))+phu*c3*c2*0.5
                    sigtyp(3,3,ityp(i)) = sigtyp(3,3,ityp(i))+phu*c3*c3*0.5
                    sigtyp(1,1,ityp(j)) = sigtyp(1,1,ityp(j))+phu*c1*c1*0.5
                    sigtyp(1,2,ityp(j)) = sigtyp(1,2,ityp(j))+phu*c1*c2*0.5
                    sigtyp(1,3,ityp(j)) = sigtyp(1,3,ityp(j))+phu*c1*c3*0.5
                    sigtyp(2,1,ityp(j)) = sigtyp(2,1,ityp(j))+phu*c2*c1*0.5
                    sigtyp(2,2,ityp(j)) = sigtyp(2,2,ityp(j))+phu*c2*c2*0.5
                    sigtyp(2,3,ityp(j)) = sigtyp(2,3,ityp(j))+phu*c2*c3*0.5
                    sigtyp(3,1,ityp(j)) = sigtyp(3,1,ityp(j))+phu*c3*c1*0.5
                    sigtyp(3,2,ityp(j)) = sigtyp(3,2,ityp(j))+phu*c3*c2*0.5
                    sigtyp(3,3,ityp(j)) = sigtyp(3,3,ityp(j))+phu*c3*c3*0.5
                    itimin=min(ityp(i),ityp(j))
                    itimax=max(ityp(i),ityp(j))
                    sigtyptyp(1,1,itimin,itimax) = sigtyptyp(1,1,itimin,itimax)+phu*c1*c1
                    sigtyptyp(1,2,itimin,itimax) = sigtyptyp(1,2,itimin,itimax)+phu*c1*c2
                    sigtyptyp(1,3,itimin,itimax) = sigtyptyp(1,3,itimin,itimax)+phu*c1*c3
                    sigtyptyp(2,1,itimin,itimax) = sigtyptyp(2,1,itimin,itimax)+phu*c2*c1
                    sigtyptyp(2,2,itimin,itimax) = sigtyptyp(2,2,itimin,itimax)+phu*c2*c2
                    sigtyptyp(2,3,itimin,itimax) = sigtyptyp(2,3,itimin,itimax)+phu*c2*c3
                    sigtyptyp(3,1,itimin,itimax) = sigtyptyp(3,1,itimin,itimax)+phu*c3*c1
                    sigtyptyp(3,2,itimin,itimax) = sigtyptyp(3,2,itimin,itimax)+phu*c3*c2
                    sigtyptyp(3,3,itimin,itimax) = sigtyptyp(3,3,itimin,itimax)+phu*c3*c3
                 end if
              endif
           endif

        end do  ! fin i2=j
     end do ! fin i1=koo
     if ((ldesinteg).and.(num_at_glob(i)==1)) then
        cv(1,1) = xpnp(1,i)-xpspr(1)
        cv(1,2) = xpnp(2,i)-xpspr(2)
        cv(1,3) = xpnp(3,i)-xpspr(3)
        call cryst_to_cart (1, cv, bg, -1) !cart vers cryst sur cv
        WHERE ( (cv.GT.0.5d0).OR.(cv.LT.-0.5d0) )
           cv(:,1:3) = cv(:,1:3) - Dble(Nint(cv(:,1:3)))
        END WHERE
        call cryst_to_cart (1, cv, at, 1) !cryst vers cart sur cv
        r = cv(1,1)*cv(1,1)+cv(1,2)*cv(1,2)+cv(1,3)*cv(1,3)
!        deltaEspr=deltaEspr+kspr*r/float(nstepdes)
        deltaEspr=0
        Espr=(1-lambdades)*kspr*r
        Espr=0
        fp(:,i)=fp(:,i)-2*kspr*cv(1,:)
     end if


  end do ! fin i

#if(PARA)
  call MPI_ALLREDUCE(potis1,potis1_tot,1,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
  potis1=potis1_tot
  call MPI_ALLREDUCE(potis2,potis2_tot,1,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
  potis2=potis2_tot
  call MPI_ALLREDUCE(sig,sig_tot,9,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
  sig=sig_tot
if(ldesinteg) then
  call MPI_ALLREDUCE(Espr,Espr_tot,1,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
  Espr=Espr_tot
  call MPI_ALLREDUCE(deltafcomp,deltaF_tot,1,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
  deltaF=deltaF+deltaF_tot
!	write(6,*)'deltaf',deltaf
!  call MPI_ALLREDUCE(deltaEspr,deltaEspr_tot,1,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
!  deltaEspr=deltaEspr_tot
endif	

#endif


  ! fin du calcul du terme de paire dans l'espace direct

  deallocate ( xpnp)
  return
end subroutine calfo2ccel


