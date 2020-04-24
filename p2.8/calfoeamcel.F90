module calfoeamcel_mod
        USE notperiod_mod,only: notperiod
        USE cryst_to_cart_mod,only: cryst_to_cart
        USE gen_com_m, ONLY:angst,at,bg,nvat,it,itesigma,low_limit,lperiod,lprteat,&
             &ltpcel,noxyz,zero,ncel,last,nato,ncel,last,nato,sigc,eatom,&
             &deltadist,nox,noy,noz,volu
        USE calfocommon

        implicit none
        contains
!----------------------------------------------------------------------
SUBROUTINE calfoeamcel(im,xp,  vp,  fp, ielat, ityp,num_at_glob)
  USE T_kind_param_m

  USE var_pot, ONLY:ipotentiel,ngrid,potiseam,potisglue,potisrep,rhomax,rhomin,eamrho,ipo,eamrep,eamglue,eamrho,rue_pot,&
       &typ_and_pot,typ_pot_pair,ipotentiel,ngrid,potiseam,potisglue,potisrep,rhomax,rhomin,eamrho,eamrho,ipo,eamrep,eamrep,&
  &eamglue,eamglue,eamrho
!  USE tab_imm_m
#ifdef PARA
  USE mod_para
#endif
  implicit none

  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  ! eam variables
  integer,intent(in)::im
  real(double),intent(inout),allocatable,dimension(:,:)::xp,vp,fp
  integer,intent(in),allocatable,dimension(:)::ityp,ielat,num_at_glob
  !local variables
  integer :: i,j !atomes
  integer ::iti,itj !types
  integer :: l !paires
  integer:: koo,ko1,ncelvois,i2,i1 !cel.
  integer :: k ! aux pour splines

  real(double) :: rue2 !coupure**2
  REAL(double), dimension(1:3) :: cp, dxp, gradij
  real(double) :: r,r2 !distance i-j
  real(double) :: Fij,cv(1,3)
  real(double) :: Eembi,dEembi ! potentiel et gradient de l'immersion
  real(double) :: rhoi,rhoj ! densite de i sur j et j sur i
  real(double) :: Erep,dErep ! potentiel et gradient de la repulsion de paire ij
  REAL(double) :: Femb
  real(double) :: rk, drk,ktor, inv_ktor, ktorho, inv_ktorho
  real(double) :: densityi !densite totale sur i

  integer ::nvi,iw
  integer, dimension(nvat) :: jvi
  real(double), dimension (nvat) ::rij
  real(double), dimension (1:3,nvat) ::dxpij
  LOGICAL :: test_sigma
  integer :: izero

  real(double) :: tabdensity(im)

#ifdef PARA
  ! declarations supplementaires pour MPI
  real(double) ::  potisglue_tot
  real(double) ::  potisrep_tot
  real(double), dimension(3,3) :: sig_tot
  real(double), dimension(3,3,noxyz) :: sigc_tot

#endif

  real(double) :: xpnp(3,im)
  real(double)::rue
  rue=rue_pot(ipotentiel)
  test_sigma=(mod(it,itesigma)==0)
  !  if (it.le.1)   write(6,*)'ROUTINE CALFOEAMCEL :boucle (i,j) complete et non j>i car parallelise *'




  ktor=rue/ngrid
  inv_ktor=1.d0/ktor
  ktorho=(rhomax-rhomin)/ngrid
  inv_ktorho = 1.d0/ktorho

  tabdensity(:)=0.
  potisrep=0.
  potisglue=0.
  rue2=rue**2


  if (lperiod) then
     xpnp(:,:)=xp(:,:)
  else 
     call notperiod(im,xp,xpnp)
  end if

 !       open(unit=806, file='CALFOGMT.csv', form='formatted', &
 !            status='unknown')
 !   do i=1,im_glob
!       if (it.ge.4)write(806,'(A,2I6,3G22.13)')'CF1 ', it,i,xp(:,i)
 !   end do
 !   do i=1,im_glob
! if (it.ge.4)       write(806,'(A,2I6,3G22.13)')'CF2 ',it,i,xpnp(:,i)
 !   end do


  loop1at1: do i=1,im
     if (typ_and_pot(ityp(i),ipotentiel).eqv..false.)cycle

     ! --- Calcul de la densite sur i ---    

     nvi=0
     densityi=0.0 ; dEembi=0.0
     koo = ielat(i)                          ! Numero de la cellule
 !if (it.ge.4)       write(806,'(A,2I6,G22.13)')'CF2B ' ,it,i,ielat(i)


     iti = ityp(i)
     ncelvois = min(noxyz,27)-1
     ! pour chaque cel. voisine
     loop1cel:   do i1 = 0, ncelvois
        ko1 = ncel(koo,i1)
 !if (it.ge.4)               write(806,'(A,4I6)')'CF3 ', it,i,ko1,nato(ko1)
        cp(1:3) = xpnp(1:3,i) + MatMul(at(1:3,:),deltadist(:,i1,koo))
        ! pour chaque atome ds la cel. voisine

        loop1at2: do i2 = 1, nato(ko1)
           j = last(i2,ko1)
 !if (it.ge.4)           write(806,'(A,4I6)')'CF4 ' ,it,ko1,i2,j
 !if (it.ge.4)         write(806,'(A,2I6,3G22.13)')'CF5 ', it,j,xp(:,j)
 !if (it.ge.4)         write(806,'(A,2I6,3G22.13)')'CF6 ', it,j,xpnp(:,j)
       
           if (typ_pot_pair(ipo(ityp(i),ityp(j))).ne.ipotentiel) cycle


           

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !CRC             if(i.eq.j) cycle
#ifdef PARA
	   ! Methode pour ne prendre qu'une seule fois en compte
           ! le couple i,j en paralle :
	   ! - i est necesairement local (boucle i<=im)
           ! - si j est local on ne retient que le couple i<j
           ! - si j n'est pas local, le couple n'est par definition
	   !   pris qu'une fois puisque i est local

	   if (j.le.im) then

       ! les deux atomes sont locaux
              if (num_at_glob(i).ge.num_at_glob(j)) cycle !terme deja calcule
	   else
       ! j n'est pas local, on ne fait le calcul normal           
	   endif
#else
           if (num_at_glob(i).ge.num_at_glob(j)) cycle !terme dÃ£Â©ja calculÃ£Â©
#endif


           if (noxyz.ne.1) then
              dxp(1) = cp(1)-xpnp(1,j)
              IF ( (dxp(1)>rue).OR.(dxp(1)<-rue) ) Cycle
              dxp(2) = cp(2)-xpnp(2,j)
              IF ( (dxp(2)>rue).OR.(dxp(2)<-rue) ) Cycle
              dxp(3) = cp(3)-xpnp(3,j)
              IF ( (dxp(3)>rue).OR.(dxp(3)<-rue) ) Cycle
           else
              dxp(1) = cp(1)-xpnp(1,j)
              dxp(2) = cp(2)-xpnp(2,j)
              dxp(3) = cp(3)-xpnp(3,j)
              cv(1,1) = dxp(1)
              cv(1,2) = dxp(2)
              cv(1,3) = dxp(3)
              call cryst_to_cart (1, cv, bg, -1) !cryst vers cart sur cv
              WHERE ( (cv.GT.0.5d0).OR.(cv.LT.-0.5d0) )
                 cv(:,1:3) = cv(:,1:3) - Dble(Nint(cv(:,1:3)))
              END WHERE
              call cryst_to_cart (1, cv, at, 1) !cryst vers cart sur cv
              dxp(1)=cv(1,1)
              dxp(2)=cv(1,2)
              dxp(3)=cv(1,3)
           end if

           do izero=1,3
	    if (dabs(dxp(izero)).lt.low_limit) then
	     dxp(izero) = zero
	    end if
           end do

           r2 = Sum(dxp(1:3)**2)
             if (r2.eq.zero*zero) & 
	      write(*,*) '1. WARNING IN calfoeamcell TWO ATOMS VERY CLOSE i ,j , dist(angst)', i ,j , sqrt(r2)*angst

           if (r2>rue2) cycle

           itj=ityp(j)
           r=sqrt(r2)
           k=Int(r*inv_ktor)
           gradij(1:3) = dxp(1:3)/r
           drk=r-k*ktor
           !             rhoj = eamrho(1,itj,k) + drk*( eamrho(2,itj,k) + drk*( eamrho(3,itj,k) + drk*eamrho(4,itj,k) ) )  !rho de j sur i
           !             densityi=densityi+rhoj

           rhoj = eamrho(1,itj,k) + drk*( eamrho(2,itj,k) + drk*( eamrho(3,itj,k) + drk*eamrho(4,itj,k) ) )  !rho de j sur i
           tabdensity(i)=tabdensity(i)+rhoj
           rhoi = eamrho(1,iti,k) + drk*( eamrho(2,iti,k) + drk*( eamrho(3,iti,k) + drk*eamrho(4,iti,k) ) )  !rho de i sur j
           tabdensity(j)=tabdensity(j)+rhoi

           !           nvi=nvi+1
           !           dxpij(1:3,nvi)=dxp(1:3)
           !           jvi(nvi)=j
           !           rij(nvi)=r            

           !terme de repulsion 
           l = ipo(iti,itj)
           Erep = eamrep(1,l,k) + drk*( eamrep(2,l,k) + drk*( eamrep(3,l,k) + drk*eamrep(4,l,k) ) )
           if(lprteat.EQV..true.)then
!              if (allocated (free)) then              
!                 if( free(i).EQV..true.) eatom(i)=eatom(i)+Erep/2.d0
!                 if( free(j).EQV..true.) eatom(j)=eatom(j)+Erep/2.d0
!              else
                  eatom(i)=eatom(i)+Erep/2.d0
                  eatom(j)=eatom(j)+Erep/2.d0
!               end if
           end if
           dErep = eamrep(2,l,k) + drk*( 2.0*eamrep(3,l,k) + 3.0*drk*eamrep(4,l,k) )

         if (num_at_glob(i).lt.num_at_glob(j)) then
!              if (allocated (free)) then              
!                 if( free(i).EQV..true.) potisrep = potisrep+Erep
!              else
                 potisrep = potisrep+Erep
!              end if
           endif

           fp(1:3,i)=fp(1:3,i)-dErep*gradij(1:3)
           fp(1:3,j)=fp(1:3,j)+dErep*gradij(1:3)

           if (test_sigma) then        
	      if (num_at_glob(i).lt.num_at_glob(j)) then          
                 sig(1:3,1) = sig(1:3,1)-dErep*gradij(1:3)*dxp(1)/volu
                 sig(1:3,2) = sig(1:3,2)-dErep*gradij(1:3)*dxp(2)/volu
                 sig(1:3,3) = sig(1:3,3)-dErep*gradij(1:3)*dxp(3)/volu
                 if (lTPcel.EQV..true.) then
                    sigc(1:3,1,koo) =sigc(1:3,1,koo) -0.5*dErep*gradij(1:3)*dxp(1)*nox*noy*noz/volu
                    sigc(1:3,2,koo) =sigc(1:3,2,koo) -0.5*dErep*gradij(1:3)*dxp(2)*nox*noy*noz/volu
                    sigc(1:3,3,koo) =sigc(1:3,3,koo) -0.5*dErep*gradij(1:3)*dxp(3)*nox*noy*noz/volu
                    sigc(1:3,1,ko1) =sigc(1:3,1,ko1) -0.5*dErep*gradij(1:3)*dxp(1)*nox*noy*noz/volu
                    sigc(1:3,2,ko1) =sigc(1:3,2,ko1) -0.5*dErep*gradij(1:3)*dxp(2)*nox*noy*noz/volu
                    sigc(1:3,3,ko1) =sigc(1:3,3,ko1) -0.5*dErep*gradij(1:3)*dxp(3)*nox*noy*noz/volu
                 end if

	      endif
           end if



        end do loop1at2
     end do loop1cel
  end do loop1at1

  
  ! calcul et stockage de Eembi et dEembi
  loop2at1: do i=1,im
     if (typ_and_pot(ityp(i),ipotentiel).eqv..false.)cycle
     iti=ityp(i)
     k=Int((tabdensity(i)-rhomin)*inv_ktorho)
     if(k.gt.ngrid) then
        write(6,*)k, ngrid, 'k> ngrid ; augmenter le facteur multiplicatif de rhomax dans calpo'
        write(6,*)'densityi',k,ngrid,densityi
        stop
     end if
     drk=tabdensity(i)-(rhomin+k*ktorho)
     Eembi = eamglue(1,iti,k) + drk*( eamglue(2,iti,k) + drk*( eamglue(3,iti,k) + drk*eamglue(4,iti,k) ) )

!     if (allocated (free)) then
!        if((lprteat.EQV..true.).and.( free(i).EQV..true.)) eatom(i)=eatom(i)+Eembi
!        if( free(i).EQV..true.)   potisglue = potisglue+Eembi
!     else
        if(lprteat.EQV..true.) eatom(i)=eatom(i)+Eembi
        potisglue = potisglue+Eembi
!     end if

     tabdensity(i) = eamglue(2,iti,k) + drk*( 2.0*eamglue(3,iti,k) + 3.0*drk*eamglue(4,iti,k) )
  end do loop2at1



#ifdef PARA
  call maj_tabdensity_ftm(tabdensity)
#endif


  !boucle des forces

  loop3at1: do i=1,im
     if (typ_and_pot(ityp(i),ipotentiel).eqv..false.)cycle

     ! --- Calcul de la densite sur i ---    

     koo = ielat(i)                          ! Numero de la cellule
     iti = ityp(i)
     ncelvois = min(noxyz,27)-1
     ! pour chaque cel. voisine
     loop2cel:   do i1 = 0, ncelvois
        ko1 = ncel(koo,i1)
        cp(1:3) = xpnp(1:3,i) + MatMul(at(1:3,:),deltadist(:,i1,koo))
        ! pour chaque atome ds la cel. voisine
        loop2at2: do i2 = 1, nato(ko1)
           j = last(i2,ko1)
           if (typ_pot_pair(ipo(ityp(i),ityp(j))).ne.ipotentiel) cycle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !CRC             if(i.eq.j) cycle
#ifdef PARA
	   ! Methode pour ne prendre qu'une seule fois en compte
           ! le couple i,j en paralle :
	   ! - i est necesairement local (boucle i<=im)
           ! - si j est local on ne retient que le couple i<j
           ! - si j n'est pas local, le couple n'est par definition
	   !   pris qu'une fois puisque i est local
	   if (j.le.im) then
       ! les deux atomes sont locaux
              if (num_at_glob(i).ge.num_at_glob(j)) cycle !terme dÃ£Â©ja calculÃ£Â©
!	   else
       ! j n'est pas local, on ne fait le calcul normal           
	   endif
#else
           if (num_at_glob(i).ge.num_at_glob(j)) cycle !terme dÃ£Â©ja calculÃ£Â©
#endif

           if (noxyz.ne.1) then
              dxp(1) = cp(1)-xpnp(1,j)
              IF ( (dxp(1)>rue).OR.(dxp(1)<-rue) ) Cycle
              dxp(2) = cp(2)-xpnp(2,j)
              IF ( (dxp(2)>rue).OR.(dxp(2)<-rue) ) Cycle
              dxp(3) = cp(3)-xpnp(3,j)
              IF ( (dxp(3)>rue).OR.(dxp(3)<-rue) ) Cycle
           else
              dxp(1) = cp(1)-xpnp(1,j)
              dxp(2) = cp(2)-xpnp(2,j)
              dxp(3) = cp(3)-xpnp(3,j)
              cv(1,1) = dxp(1)
              cv(1,2) = dxp(2)
              cv(1,3) = dxp(3)
              call cryst_to_cart (1, cv, bg, -1) !cryst vers cart sur cv
              WHERE ( (cv.GT.0.5d0).OR.(cv.LT.-0.5d0) )
                 cv(:,1:3) = cv(:,1:3) - Dble(Nint(cv(:,1:3)))
              END WHERE
              call cryst_to_cart (1, cv, at, 1) !cryst vers cart sur cv
              dxp(1)=cv(1,1)
              dxp(2)=cv(1,2)
              dxp(3)=cv(1,3)
           end if

           do izero=1,3
	    if (dabs(dxp(izero)).lt.low_limit) then
	     dxp(izero) = zero
	    end if
           end do


           r2 = Sum(dxp(1:3)**2)
           if (r2>rue2) cycle

           itj=ityp(j)
           r=sqrt(r2)
           if (r.eq.zero) & 
	     write(*,*) '2. WARNING IN calfoeamcell TWO ATOMS VERY CLOSE i ,j , dist(angst)', i ,j , r*angst

           k=Int(r*inv_ktor)
           drk=r-k*ktor
           gradij(1:3) = dxp(1:3)/r
           Femb = ( eamrho(2,itj,k) + drk*( 2.0*eamrho(3,itj,k) + 3.0*drk*eamrho(4,itj,k) ) )*tabdensity(i) &
                + ( eamrho(2,iti,k) + drk*( 2.0*eamrho(3,iti,k) + 3.0*drk*eamrho(4,iti,k) ) )*tabdensity(j)
           fp(1:3,i) = fp(1:3,i) - Femb*gradij(1:3)
           fp(1:3,j) = fp(1:3,j) + Femb*gradij(1:3)

           if (test_sigma) then                   
              if (num_at_glob(i).lt.num_at_glob(j)) then
                 sig(1:3,1) = sig(1:3,1) - Femb*gradij(1:3)*dxp(1)/volu
                 sig(1:3,2) = sig(1:3,2) - Femb*gradij(1:3)*dxp(2)/volu
                 sig(1:3,3) = sig(1:3,3) - Femb*gradij(1:3)*dxp(3)/volu
                 if (lTPcel.EQV..true.) then
                    sigc(1:3,1,koo) =sigc(1:3,1,koo) - 0.5*Femb*gradij(1:3)*dxp(1)*nox*noy*noz/volu
                    sigc(1:3,2,koo) =sigc(1:3,2,koo) - 0.5*Femb*gradij(1:3)*dxp(2)*nox*noy*noz/volu
                    sigc(1:3,3,koo) =sigc(1:3,3,koo) - 0.5*Femb*gradij(1:3)*dxp(3)*nox*noy*noz/volu
                    sigc(1:3,1,ko1) =sigc(1:3,1,ko1) - 0.5*Femb*gradij(1:3)*dxp(1)*nox*noy*noz/volu
                    sigc(1:3,2,ko1) =sigc(1:3,2,ko1) - 0.5*Femb*gradij(1:3)*dxp(2)*nox*noy*noz/volu
                    sigc(1:3,3,ko1) =sigc(1:3,3,ko1) - 0.5*Femb*gradij(1:3)*dxp(3)*nox*noy*noz/volu
                 end if

	      endif
           end if

        end do loop2at2
     end do loop2cel
  end do loop3at1


!  if (test_sigma) sig(1:3,1:3) = sig(1:3,1:3)/volu

#ifdef PARA
  CALL MPI_ALLREDUCE(potisrep, potisrep_tot, 1,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
  potisrep=potisrep_tot
  CALL MPI_ALLREDUCE(potisglue,potisglue_tot,1,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
  potisglue=potisglue_tot
!  write(6,'(I3,9G15.7)')rang,sig
 if (test_sigma) then 
     call MPI_ALLREDUCE(sig,      sig_tot,      9,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
     sig=sig_tot
       if (allocated(sigc)) then
     call MPI_ALLREDUCE(sigc,      sigc_tot,      9*noxyz,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
     sigc=sigc_tot
  endif
endif
#endif

  potiseam=potisglue+potisrep

  return
end SUBROUTINE calfoeamcel
end module
