module calpo_mod
  USE spline_mod,only: cspline
  USE zieg2_mod,only: zieg2

  USE dervbeest_mod,only: deriVBEEST,maxVBEEST,potvbeest

  USE arret_ndm_mod,only: arret_ndm
  USE potrep_mod,only: potrep
  USE calerf_mod,only: calerf
  USE gen_com_m, ONLY:ecgs,half,one,precexp,rang,pi
  USE var_pot, ONLY:bspg,ngrid,catom,csive,csive_g,cspg,dip,dspg,gd,gm2,gm3,gm4,gm5,gr,gz,ipotentiel,ipotrep,&
       &npair,lprtpot,lu_roff_pair,ngr,ntyp,pot_d,roff1,roff2,rue_pair,sigmawat,ro,lue_paire,ipo,pot,rawat,qwat,&
       &ipo_2_pair_tab,zz,ipo,capdij,caphij,capwij,gm1,ietaij,lambda,rbp5,rp3c,rp5p3,xsi,poly5,poly3,r8p,pwat,&
       &typ_pot_pair,pot_pair_tab,ray,a_factor,fcr,potw,bspw,cspw,dspw,bspf,dspf,cspf,shel,pm,bwat,bm,awat,&
       &alpha,auxe,iewald,q,Afd,Bfd,r0fd, Aig,big,r0ig,dmorse,remorse,amorse

  
  implicit none
contains

  subroutine calpo
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m


    implicit none
    !-----------------------------------------------
    !   G l o b a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   L o c a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------

    integer :: i, l, k, j
    real(double), dimension(npair) :: sigp
    real(double) :: bmh, r, r2, r3, r4, r5, r6, r8
    real(double), dimension(npair) :: pau
    real(double) :: factor, ar, ar2, damp, ddam

    real(double), dimension(ngrid) ::  kxsp, &
         potpart ,fcpart
    real(double), dimension (ngrid) ::potpartw
    real(double), dimension(ngrid) :: bsppart, csppart, dsppart


    !real(double),external  :: derfc


    real(double) :: &
         f2,f2exp, & ! intermediaires de calcul
         gexp1,gexp2,fcr1,&  ! intermediaires de calcul
         kgz(0:ngrid+1),z  ! valeur de z, variable z


    real(double), parameter:: maxSiO=8.0       ! 8 c'est deja beaucoup
    ! repulsion polynomiale

    integer :: i1,i2,lw

    integer, dimension(npair) :: irrep
    integer :: convrep
    real(double), dimension(npair) :: r0rep, V0rep
    real(double) :: rrep, potV0


    real(double),allocatable::potrc(:),dpotrc(:)  ! redressement en rc


    real(double):: drkp,skp
    integer::kp,lpt

    factor = (2.0D0*alpha)/sqrt(pi)
555 format(1x,'Q =',f5.1,3x,'RAY =',f6.2,3x,'BM =',f7.4,3x,'N=',f4.1)
557 format(1x,4f10.3)                          ! utile ???

    ! ****************************************************************
    !    CALCUL DU POTENTIEL D'INTERACTION ENTRE 2 TYPES DE PAIRE
    !    CHOIX ENTRE : 0. Born-Mayer-Huggins et 1. Buckingham
    ! ****************************************************************
    write(6,*)'ipotentiel ',ipotentiel
    select case (ipotentiel)
    case(0,1,3,4,5,8)  ! FORMULES ANALYTIQUES
       sigp=0.0
       ! +++++++++ 0. POTENTIEL DE BORN-MAYER-HUGGINS +++++++++++++
       !         ip2: select case (ipotentiel)
       !         case(0)
       select case (ipotentiel)
       case(0)
          l = 0
          do i = 1, ntyp
             if (ntyp-i+1>0) then
                sigp(l+1:ntyp-i+1+l) = ray(i)+ray(i:ntyp)
                pau(l+1:ntyp-i+1+l) = (one+q(i)/shel(i)+q(i:ntyp)/shel(i:ntyp))*&
                     (bm(i)+bm(i:ntyp))*half
                l = ntyp-i+1+l
             endif
          end do


          !                                                ! erg
          pau(:npair) = pau(:npair)*ecgs*pm(:npair)
          sigp(:npair) = sigp(:npair)*1.D-08      ! cm
          !       dip(l)=dip(l)*65.4D-60*1.6021892    ! conversion deja faite dans input.f

          ! exponential term
          do k = 1, ngrid
             r = float(k)*csive
             r2 = r*r                             ! utile ????
             kxsp(k) = r
             r3 = r2*r                            ! utile ????
             r6 = r3*r3
             r8 = r6*r2                           ! utile ????
             do l = 1, npair
                !                                                ! erg
                bmh = pau(l)*exp((sigp(l)-r)/ro(l))
                pot(1,l,k) = bmh-dip(l)/r6          ! erg
             end do
          end do
          ! ++++++++++++ Fin de Born-Mayer-Huggins +++++++++++++

          ! +++++++++++ 1. POTENTIEL UO2 +++++++++++++
          !case(4)
       case(4)
          pau(:npair) = a_factor(:npair)          ! erg
          ! exponential term
          if (rang==0) then
             write (6,*)
             write (6,*)
             write (6,*)
             write (6,*)
             write(6,*)'ATTENTION!!! UO2 === Oxygène = TYPE 2!!!!'
             write (6,*)
             write (6,*)
             write (6,*)
             write (6,*)
          endif

          do k = 1, ngrid
             r = float(k)*csive
             r2 = r*r
             kxsp(k) = r
             r3 = r2*r
             r4 = r2*r2
             r5 = r4*r
             r6 = r3*r3
             do l=1,npair
                if((typ_pot_pair(l)==4).and.(lue_paire(l).eqv..true.)) then
                   pot(1,l,k) = pau(l)*exp((-r)/ro(l))-dip(l)/r6 
                end if
             end do

             !           pot(1,:npair,k) = pau(:npair)*exp((-r)/ro(:npair))-dip(:npair)/r6
             ! cas particulier de l'interaction O-O
             l=ipo(2,2)

             if (r < rbp5) pot(1,l,k) = pau(l)*exp((-r)/ro(l))
             if (r < rp5p3.and.r >= rbp5) pot(1,l,k) = poly5(1)+poly5(2)*r+ &
                  poly5(3)*r2+poly5(4)*r3+poly5(5)*r4+poly5(6)*r5
             if (r < rp3c.and.r >= rp5p3) pot(1,l,k) = poly3(1)+poly3(2)*r+ &
                  poly3(3)*r2+poly3(4)*r3
             if (r > rp3c) pot(1,l,k) = -dip(l)/r6

             !               write (6,*) k,r,pot(1,l,k) 
          end do
          ! ++++++++++++ Fin du potentiel UO2 +++++++++++++
       case default
          !         case(1,3,5,8)

          ! +++++++++++ 1. POTENTIEL DE BUCKINGHAM +++++++++++++

          pau(:npair) = a_factor(:npair)          
          !       write(6,*)'L, Aij, ROij, Cij',l,pau(l),
          !    &            ro(l),dip(l)

          ! exponential term
          do k = 1, ngrid
             r = float(k)*csive
             r2 = r*r                             ! utile ????
             kxsp(k) = r
             r3 = r2*r                            ! utile ????
             r6 = r3*r3
             r8 = r6*r2                           ! utile ????
             !                                                ! erg
             if (ipotentiel==3) then
                do l=1,npair
                   if((typ_pot_pair(l)==3).and.(lue_paire(l).eqv..true.)) then
                      pot(1,l,k) = pau(l)*exp((-r)/ro(l))-dip(l)/r6 +r8p(l)/r8
                   end if
                end do
             else
                do l=1,npair
                   if((typ_pot_pair(l)==5).and.(lue_paire(l).eqv..true.)) then
                      pot(1,l,k) = pau(l)*exp((-r)/ro(l))-dip(l)/r6 
                   end if
                end do
                do l=1,npair
                   if((typ_pot_pair(l)==1).and.(lue_paire(l).eqv..true.)) then
                      pot(1,l,k) = pau(l)*exp((-r)/ro(l))-dip(l)/r6

                   end if
                end do

                !              where ((typ_pot_pair==5).or.(typ_pot_pair==1)) &
                !                   &          pot(1,:,k) = pau(:)*exp((-r)/ro(:))-dip(:)/r6
             endif
             if (ipotentiel==5) then
                do l=1,npair
                   if((typ_pot_pair(l)==5).and.(lue_paire(l).eqv..true.)) then
                      pot(1,l,k) = pau(l)*exp((-r)/ro(l))-dip(l)/r6 
                   end if
                end do

                              where (typ_pot_pair==5)& 
                                   &              pot(1,:,k)=pot(1,:,k)+dmorse(:)*((1.-exp(-1.*amorse(:)*(r-remorse(:))))**2 -1.)
             end if
!8888888888888             
             if (ipotentiel==8) then 
                do l=1,npair
                   if((typ_pot_pair(l)==8).and.(lue_paire(l).eqv..true.)) then
                      pot(1,l,k) = pau(l)*exp((-r)/ro(l))-dip(l)/r6
                      pot(1,l,k) =  pot(1,l,k)+ dmorse(l)*((1.-exp(-1.*amorse(l)*(r-remorse(l))))**2 -1.)
                      pot(1,l,k) =  pot(1,l,k) -afd(l)/(1+exp(bfd(l)*(r-r0fd(l))))
                      pot(1,l,k) =  pot(1,l,k) -aig(l)*exp(-big(l)*((r-r0ig(l))**2))
                   end if
                end do
             end if
                
          end do


       end select
       !        end select ip2


       ! ++++++++++++++++++++ Fin de Buckingham ++++++++++++++++++++

       ! ********* FIN POTENTIEL INTERACTION DE TYPES DE PAIRES **********


       ! ************************************************************
       !               CALCUL DU POTENTIEL COULOMBIEN
       ! ************************************************************
       ! Produits des charges entre 2 types (pour terme coulombien)
       l = 0
       do i = 1, ntyp
          if (ntyp-i+1>0) then
             zz(l+1:ntyp-i+1+l) = q(i)*q(i:ntyp)
             l = ntyp-i+1+l
          endif
       end do

       ! - Tableau des potentiels et forces correspondant aux interactions coulombiennes
       ! --- Le premier terme du potentiel et de la force est calcule ---
       do k = 1, ngrid
          r = k*csive
          r2 = r*r
          r3 = r2*r
          ar = alpha*r
          ar2 = ar*ar
          damp = 0.0
          ! calcul de derfc par sous routine exterieure
          damp = derfc(ar)
          ddam = factor*r*exp((-ar2))
          !      write(6,*)'r k dampam ',r, k ,damp+ddam

          !    interaction de paire + interaction couenne
          pot(1,:npair,k) = pot(1,:npair,k)+auxe*zz(:npair)*damp/r

       end do
       do l=1,npair
          potpart(:ngrid) = pot(1,l,1:ngrid)
          call cspline (ngrid, kxsp, potpart, bsppart, csppart, dsppart)
          pot(2,l,1:ngrid) = bsppart(:ngrid)
          pot(3,l,1:ngrid) = csppart(:ngrid)
          pot(4,l,1:ngrid) = dsppart(:ngrid)
          !write(6,*)pot(1,l,k),pot(2,l,k),pot(3,l,k),pot(4,l,k)
       end do


       !*********************************cas : rpulsion par Pot polynomial******************************
       select case (ipotrep)
       case(1)
          ! Calcul du premier maximum local

          if (rang==0) write(6,*)'calcul du max loc du pot VBEEST'
          irrep(:npair)=0
          l=0
          do i1=1,ntyp
             do i2=1,ntyp
                l=ipo(i1,i2)
                if (lue_paire(l)) then
                   if(irrep(l)==20) cycle
                   !           write(6,*)'entre maxVBEEST'
                   call maxVBEEST(rrep,csive,l,auxe,alpha,ngrid,ntyp, &
                        npair,pau,dip,ro,zz,convrep)
                   !           write(6,*)'sortie max VBEEST'
                   irrep(l)=convrep
                   r0rep(l)=rrep
                   call potVBEEST(potV0,rrep,l,auxe,alpha,ngrid, &
                        ntyp,npair,pau,dip,ro,zz)
                   V0rep(l)=potV0
                endif
             enddo
          enddo

          ! Terme repulsif a courte distance
          call potrep(csive,r0rep,V0rep,ngrid,ntyp,npair)



          !***********************************************************************************


       case(2)  !(lpotrep)

          ! ********************************Cas : Pot  de Ziegler ****************************lprtpo

          ! **** terme de Ziegler *******

          call zieg2 (pot,pot_d, csive,ngrid, ntyp,npair,catom,roff1,roff2,lu_roff_pair,ipotentiel,typ_pot_pair,ipo)

       case default    !(lpotrep)
       end select


       ! affichage potentiel de chaque paire
       !          do l=1,npair
       !             if (typ_pot_pair(l)==ipotentiel)then
       !                write(6,*)'l,k,r,pot(1,l,k)'
       !                do k=10,ngrid,10
       !                   r=float(k)*csive*1.0D8
       !                   write(6,'(2I6,3D17.6)')l,k,r,pot(1,l,k),pot(2,l,k)
       !                enddo
       !             end if
       !          endd

    case(2)
       if (rang==0) write (6, *) '----------- POTENTIEL WATANABE --------------'

       ! -----Terme a 2 corps de base
       do l=1,npair
          do k=1,ngrid
             kxsp(k)=float(k)*csive
             !unites reduites
             r=float(k)*csive
             if ((r-rawat(l)).gt.-precexp*sigmawat) then
                f2exp = 0.0
             else
                f2exp = EXP(sigmawat/(r-rawat(l)))
             endif
             f2=Awat(l)*(Bwat(l)*(r/sigmawat)**(-pwat(l))-(r/sigmawat)**(-qwat(l)))
             potw(l,k)=f2*f2exp
          enddo
       enddo

       !-------Fonction g(z)
       csive_g= maxSiO/float(ngrid)
       do k=0,ngrid+1     ! +1 pas vraiment necessaire...
          z=float(k)*csive_g
          kgz(k)= z
          gexp1= EXP((gm2-z)/gm3)
          gexp2= EXP(gm4*(z-gm5)**2)
          gz(k)= gm1*gexp2/(gexp1+1.0)
       enddo

       !------Fonction fc(r) pour determiner z
       do k=1,ngrid
          !on repasse en unites reduites
          r= float(k)*csive

          if(r.lt.(gR-gD)) then
             fcr(k)=1.0
          elseif(r.lt.(gR+gD))then
             fcr1= (r-gR+gD)/gD
             fcr(k)= 1.0-fcr1/2.0+ SIN(pi*fcr1)/(2*pi)
          else
             fcr(k)=0.0
          endif
       enddo


    case(7)
       ! calculer pot par le spline de  pot_pair_tab
       ! puis resplinner
 !      write(6,*)'csive',csive
       loopk:     do k=1,ngrid
          r= float(k)*csive
          kxsp(k) = r
          do l=1,npair

             if (typ_pot_pair(l)==ipotentiel)then
                lpt=ipo_2_pair_tab(l)
!                write(6,*)'l',l,k, r,pot_pair_tab(ngr,0,lpt)
                if (r.gt.pot_pair_tab(ngr,0,lpt)) then
                   if (rang==0) write(6,*)'pot tab pair trop court',r,k,pot_pair_tab(ngr,0,lpt),l,lpt
                   call arret_ndm
                end if
                !              if (r.lt.pot_pair_tab(1,0,lpt)) cycle loopk



                skp=r*ngr/pot_pair_tab(ngr,0,lpt)
                kp=skp
                drkp = (skp-kp)/(ngr/pot_pair_tab(ngr,0,lpt))
                !              write(6,*)k,r,kp,l              
                !              pot(1,l,k) = pot_pair_tab(kp,1,lpt)+ r*(drkp*(pot_pair_tab(kp,2,lpt)+drkp*(pot_pair_tab(kp,3,lpt) +drkp*(pot_pair_tab(kp,4,lpt)))))
                pot(1,l,k) = pot_pair_tab(kp,1,lpt)+ drkp*(pot_pair_tab(kp,2,lpt)+drkp*(pot_pair_tab(kp,3,lpt) &
                     & +drkp*(pot_pair_tab(kp,4,lpt))))
                !                           write(6,*)r,pot(1,l,k)              

                if (iewald.ne.0) then
                   pot(1,l,k) = pot_pair_tab(kp,1,lpt)+ drkp*(pot_pair_tab(kp,2,lpt)+drkp*(pot_pair_tab(kp,3,lpt) &
                        & +drkp*(pot_pair_tab(kp,4,lpt))))

                end if
             end if
          end do
       end do loopk


       ! Si il existe des espèces chargées
       ! - Tableau des potentiels et forces correspondant aux interactions coulombiennes
       if (iewald.ne.0) then
          do i=1,ntyp
             do j=i,ntyp
                l=ipo(i,j)
                if (typ_pot_pair(l)==ipotentiel)then
                   zz(l)=q(i)*q(j)
                end if
             end do
          end do

          do k = 1, ngrid
             r = k*csive
             r2 = r*r
             r3 = r2*r
             ar = alpha*r
             ar2 = ar*ar
             damp = 0.0
             ! calcul de derfc par sous routine exterieure
             damp = derfc(ar)
             ddam = factor*r*exp((-ar2))
             !      write(6,*)'r k dampam ',r, k ,damp+ddam

             !    interaction de paire + interaction couenne
             do l=1,npair
                if (typ_pot_pair(l)==ipotentiel)then
                   if (zz(l).ne.0) then
                      pot(1,l,k) = pot(1,l,k)+auxe*zz(l)*damp/r
                   end if
                end if
             end do
          end do
       end if

       if (ipotrep==2) then
          call zieg2 (pot,pot_d, csive,ngrid, ntyp,npair,catom,roff1,roff2,lu_roff_pair,ipotentiel,typ_pot_pair,ipo)
       end if


    case(6)
       allocate(potrc(npair));allocate(dpotrc(npair))
       do i=1,ntyp
          do j=i,ntyp
             l=ipo(i,j)
             zz(l)=q(i)*q(j)
          end do
       end do
       r=maxval(rue_pair)
       potrc(:)=capHij(:)/(r**ietaij(:)) +auxe*zz(:npair)*exp(-r/lambda)/r &
            & -capDij(:)*exp(-r/xsi)/(2*r**4) -capWij(:)/r**6
       !ecriture vectorielle implicite
       dpotrc=-ietaij*capHij/r**ietaij -zz*auxe*exp(-r/lambda)*(1./r**2+1./(lambda*r))&
            &+capDij*exp(-r/xsi)*(2.0/r**5+1./(2*xsi*r**4)) +6*capWij/r**7




       do k = 1, ngrid
          r = float(k)*csive
          kxsp(k) = r

          !                                                ! erg
          pot(1,:npair,k)=capHij(:)/(r**ietaij(:)) +auxe*zz(:npair)*exp(-r/lambda)/r &
               & -capDij(:)*exp(-r/xsi)/(2*r**4) -capWij(:)/r**6 &
               &    -potrc(:)-(r-rue_pair(:))*dpotrc(:)
       end do

       !     do l=1,npair
       !        do k = 1, ngrid
       !           r = float(k)*csive
       !           write (l,*) k,r,pot(1,l,k) 
       !        end do
       !     end do


    case default
       write (6, *) rang,'Bienvenue dans le cote obscur de la force : pas de potentiel ?'
       call arret_ndm
    end select



    ! spline
    if(ipotentiel.eq.2)then
       do l = 1, npair
          potpartw(:ngrid) = potw(l,1:ngrid)

          call cspline (ngrid, kxsp, potpartw, bsppart, csppart, dsppart)
          bspw(l,1:ngrid)=bsppart(1:ngrid)
          cspw(l,1:ngrid)=csppart(1:ngrid)
          dspw(l,1:ngrid)=dsppart(1:ngrid)

          bspw(l,0) = 0.0
          cspw(l,0) = 0.0
          dspw(l,0) = 0.0
          bspw(l,ngrid+1) = 0.0
          cspw(l,ngrid+1) = 0.0
          dspw(l,ngrid+1) = 0.0

       enddo

    else

       do l = 1, npair
          if (typ_pot_pair(l)==ipotentiel)then
             potpart(1:ngrid) = pot(1,l,1:ngrid)
             !           write(6,*)'uuuuuuu'
             !           write(6,*)potpart
             call cspline (ngrid, kxsp, potpart, bsppart, csppart, dsppart)
             pot(2,l,1:ngrid) = bsppart(1:ngrid)/kxsp(1:ngrid)
             pot(3,l,1:ngrid) = csppart(1:ngrid)/kxsp(1:ngrid)
             pot(4,l,1:ngrid) = dsppart(1:ngrid)/kxsp(1:ngrid)

             pot(2,l,0) = 0.0
             pot(3,l,0) = 0.0
             pot(4,l,0) = 0.0
             pot(2,l,ngrid+1) = 0.0
             pot(3,l,ngrid+1) = 0.0
             pot(4,l,ngrid+1) = 0.0
          end if

       enddo
    endif

    if (lprtpot.EQV..true.) then
       do l=1,npair
          if (typ_pot_pair(l)==ipotentiel)then
             write(6,*)'l,k,r,pot(1,l,k)'
             do k=1,ngrid
                r=float(k)*csive*1.0D8
                lw=320+l
                write(lw,'(2I6,5D15.6)')l,k,r,pot(1,l,k),pot(2,l,k),pot(3,l,k),pot(4,l,k)
             enddo
          end if
       enddo
    end if



    if(ipotentiel.eq.2) then ! ben non, pas fini
       call cspline(ngrid+2,kgz,gz, bspg,cspg,dspg) ! pour g


       fcpart(1:ngrid)=fcr(1:ngrid)

       call cspline(ngrid,kxsp,fcpart,bsppart,csppart,dsppart) !pour f

       bspf(1:ngrid)=bsppart(1:ngrid)
       cspf(1:ngrid)=csppart(1:ngrid)
       dspf(1:ngrid)=dsppart(1:ngrid)
       bspf(0)=0.0
       cspf(0)=0.0
       dspf(0)=0.0
       bspf(ngrid+1)=0.0
       cspf(ngrid+1)=0.0
       dspf(ngrid+1)=0.0
    endif




    return
  end subroutine calpo


  !******************************************************************
  real(kind(0.0d0)) function fac (ll)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double
    !******************************************************************
    implicit none
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer , intent(in) :: ll
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: i
    !-----------------------------------------------
    fac = 1
    do i = 1, ll
       fac = fac*i
    end do
    return
  end function fac
end module calpo_mod
