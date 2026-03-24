module Parrinello_Rahman

  ! Algorithme de Parrinello-Rahman [1] combiné ou non à un thermostat de Nosé-Hoover [4]
  ! basé sur les article de Ray et Rahman [2,3]
  ! Si le thermostat de Nosé-Hoover n'est pas branché (lTHoover=.false.), la
  ! viscosité zHoover restera nulle.
  !
  ! Quelques remarques
  ! * La vitesse des particules ne tient pas compte de la dérivée du tenseur h:
  ! * On impose une tension thermodynamique plutôt qu'une contrainte constante:
  !   eq. (2.22) de [2]
  !   Du coup, la matrice h0R définissant l'état de référence n'apparaît plus
  !   nulle part dans l'algorithme. Cette matrice est nécessaire seulement si à
  !   un instant donné on souhaite calculer la déformation epsilon de la boîte
  !   et également son énergie potentielle Ucell, mais ces 2 quantités ne sont
  !   pas nécessaires à l'algorithme.
  ! * Le thermostat est imposé sur l'énergie cinétique des atomes et de la
  !   cellule.
  ! * L'algorithme utilisé pour intégrer les équations de mouvement est
  !   "velocity Verlet". Les positions et les vitesses sont donc correctes
  !   jusqu'à respectivement l'ordre 3 et 2 inclus.
  ! * Pour que l'algorithme soit correcte, les forces ne doivent pas faire
  !   intervenir la vitesse des particules. Des thermostats du type lTcon=.true.
  !   ou lTBerendsen=.true. ne sont donc pas autorisés.
  !
  !
  ! [1] Parrinello, M. & Rahman,
  !     A. Polymorphic Transitions in Single rcystals: A New Molecular Dynamics Method
  !     J. Appl. Phys., 1981, 52, 7182-7190
  ! [2] Ray, J.R. & Rahman, A.controleT
  !     Statistical Ensembles and Molecular Dynamics Studies of Anisotropic Solids
  !     J. Chem. Phys., 198bulk.crcin4, 80, 4423-4428
  ! [3] Ray, J.R. & Rahman, A.
  !     Statistical Ensembles and Molecular Dynamics Studies of Anisotropic Solids. II
  !     J. Chem. Phys., 1985, 82, 4243-4247
  ! [4] Nosé, S.
  !     A Molecular Dynamics Method for Simulations in the Canonical Ensemble
  !     Mol. Phys., 1984, 52, 255-268tabv
  USE T_kind_param_m
  USE gen_com_m, only:uwrt,lwrt,kine,lpcon2,lthoover,sigext,erg2ev,&
       &leev,lucell,nhoover,timel,wboxf, ihbox0 ,tbox, bk,&
       &potist,sigtot,text,tstep,iteration,potist,rang,sig,text,sigkine,lpcube,&
       &pi,l2t,ltberendsen,lperiod,lspaceNDM,h0R,dmtype,usdh,llangevin,gamlg,gamprfact,unitP,&
       & lmaxvp,vplim,astarsig,sig0dir,thsig,lpconxyz

  use FireModule,only:alph_start,f_alph,fdec,finc,nstepmin,tstep_mm,tstep0,init_trempe_fire

  USE var_pot, ONLY:cm,ntyp,gamlt
  USE recips_mod,only: recips,calcvol
#ifdef PARA
  use Tpara, only:nprocspace,ierr,comm_space,myidsp
#else
  use Tpara, only:nprocspace,ierr,comm_space,myidsp
#endif
  USE calfo_mod,only: calfo
  USE scalebox_mod,only: scalebox
  USE atomconfig,only : atom_config_d,atom_config_e
  USE cellconfig, only:cell_config,caltabtc
  USE boxconfig,only:box_config_lpr,box_config,periodbox,updatebox
  USE eloss, ONLY : calceloss,ibrake !, tcalfcelec,ecelec,ibrake,elstopforce,elosselectot,elosselectot1,elosselec1,ngrdel,elosselec
  USE elec_cell, ONLY :i2t
  USE calfoberend_mod,only:calfoberend
  use Tpara,only:para_space_config
  USE calpo_ew_mod,only: calpo_ew
  USE calctemp_mod,only: calctemp
  USE calpo_ew_mod,only: calpo_ew
  use sigkinetot_mod,only:sigkinetot
  USE tempinstT_mod,only: tempinstT
  use mat_utils_mod

  implicit none
!!$  ! Vecteurs de la boîte et leurs dérivées
!!$  real(double), dimension(3,3), save , private :: h, hDot
!!$  real(double), dimension(3,3), save , private :: trh, invh, invtrh, Gmat, invGmat, Gdot
!!$  real(double), save, private :: invVolu

  ! Coordonnées réduites des atomes et leurs dérivées
  real(double), allocatable :: sp(:,:), sdot(:,:), sdot_new(:,:),sfp(:,:),spp(:,:)



  

  ! Nombre de degrés de liberté

  !  real(double)::wbox
  ! Variables uniquement nécessaires au calcul de l'énergie potentielle de la
  ! boîte
  REAL(double) ::  fire_alph
  INTEGER :: fire_nstep,ic

  real(double)::TInitBox

contains

  subroutine initlpr (atpr,celndm,boxndm,psc,lwrtprR,uwrtprR)

    implicit none
    type(para_space_config)::psc
    class(box_config_lpr)::boxndm
    class(atom_config_d)::atpr
    type(cell_config):: celndm
    logical, optional :: lwrtprR
    integer,optional ::uwrtprR
    
    logical :: lwrtpr
    integer ::uwrtpr


    INTEGER :: ia, i, j
    !real(double), external :: calcvol
    real(double):: unitE,v1,z1,z2
    character*5 :: cunitE
    real(double)::pint
    if (present(lwrtprR))then
       lwrtpr=lwrtprR
    else
       lwrtpr=lwrt
    end if
    if (present(uwrtprR)) then
       uwrtpr=uwrtprR
    else
       uwrtpr=uwrt
    end if

    if (dmtype==24) then
       CALL init_trempe_fire(tstep, fire_nstep, fire_alph)
    END IF

    IF(lwrtpr) WRITE(uwrtpr,*)
    if (dmtype.ne.15) then
       if (llangevin) then 
          IF(lwrtpr) WRITE(uwrtpr,*) 'Algorithme deP cst Langevin Parrinello-Rahman '
       else
          IF(lwrtpr) WRITE(uwrtpr,*) 'Algorithme de Parrinello-Rahman (V2)'
       end if
       IF(lwrtpr) WRITE(uwrtpr,'(a)') '  -> la vitesse de la boîte ne prend pas en compte la dérivée du tenseur h à t=0'
       IF(lwrtpr) WRITE(uwrtpr,*)
    end if

    IF (lUcell) THEN
       IF(lwrtpr) WRITE(uwrtpr,'(a)') "Repère de référence pour Parrinello-Rahman  (A):"
       IF(lwrtpr) WRITE(uwrtpr,'(a,3(f0.5,1x))') ' h0(1:3,1) = ', 1e8*h0R(1:3,1)
       IF(lwrtpr) WRITE(uwrtpr,'(a,3(f0.5,1x))') ' h0(1:3,2) = ', 1e8*h0R(1:3,2)
       IF(lwrtpr) WRITE(uwrtpr,'(a,3(f0.5,1x))') ' h0(1:3,3) = ', 1e8*h0R(1:3,3)
       IF(lwrtpr) WRITE(uwrtpr,*)
       boxndm%h0 = h0R
    ELSE
       boxndm%h0 = boxndm%at
    END IF
    IF(lwrtpr) WRITE(uwrtpr,'(a)') "Repère actuel  (A):"
    IF(lwrtpr) WRITE(uwrtpr,'(a,3(f0.5,1x))') ' h (1:3,1) = ', 1e8*boxndm%at(1:3,1)
    IF(lwrtpr) WRITE(uwrtpr,'(a,3(f0.5,1x))') ' h (1:3,2) = ', 1e8*boxndm%at(1:3,2)
    IF(lwrtpr) WRITE(uwrtpr,'(a,3(f0.5,1x))') ' h (1:3,3) = ', 1e8*boxndm%at(1:3,3)
    boxndm%wbox =wboxf*sum(0.5*cm(atpr%ityp(:atpr%im)))       ! La moitié de la masse totale des atomes
#ifdef PARA
    if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
       call comm_space%sum(boxndm%wbox)
    end if
#endif


    IF(lwrtpr) WRITE(uwrtpr,'(a,g20.12)')'Masse de la boîte pour Parrinello-Rahman: wbox=',boxndm%wbox

    ! État de référence défini par la matrice h0
    !   Cet état de référence doit correspondre à un tenseur de contrainte nul.
    !   Il n'est utile que pour calculer la déformation et l'énergie potentielle
    !   de la boîte.
    boxndm%volu0 = calcvol(boxndm%h0(1:3,1),boxndm%h0(1:3,2),boxndm%h0(1:3,3))
    boxndm%invVolu0 = 1.d0/boxndm%volu0
    boxndm%trh0=Transpose(boxndm%h0)
    CALL MatInv(boxndm%h0,boxndm%invh0)
    boxndm%invtrh0=Transpose(boxndm%invh0)

    ! Vecteurs de la boîte et grandeurs associées à l'instant initial
!!$    h(:,:)=boxndm%at(:,:)
!!$    trh=Transpose(h)
!!$    Gmat = MatMul(trh,h)
!!$    CALL MatInv(Gmat,invGmat)
!!$    call MatInv(h,invh)
!!$    invtrh = Transpose(invh)
!!$    boxndm%volu = calcvol(h(1:3,1),h(1:3,2),h(1:3,3))
!!$    invVolu = 1.d0/boxndm%volu

    ! Initialisation de la vitesse de la boîte

    IF(lwrtpr) WRITE(uwrtpr,'(a,f0.3,a)') 'Initialisation de la vitesse de la boîte pour la température ', TinitBox, ' K'
    boxndm%hdot(:,:) = 0.d0
    !#ifdef PARA
    if (all(ihbox0==1))then 
       if (myidsp==0) then
          !#endif
          DO i=1, 3
             DO j=1, 3
                call random_number(z1)
                call random_number(z2)
                if(z1.eq.0.d0) z1=0.000000001d0
                if(z2.eq.0.d0) z2=0.000000001d0

                v1 =  sqrt(-2*log(z1))*cos(2*pi*z2)

                boxndm%hdot(1:3,1:3)=sqrt(2*bk*Tinitbox/boxndm%Wbox)*v1
             end DO
          end DO
          !#ifdef PARA
       endif
    end if
    !   boxndm%hdot=0
    call comm_space%bcast(0,boxndm%hdot)

    !#endif

    boxndm%Kcell = 0.5d0*boxndm%wbox*Sum( boxndm%hDot(1:3,1:3)**2 )
    boxndm%Tempcell=boxndm%Kcell*2./(9.*bk)
    DO i=1, 3
       DO j=1, 3

          ! Time derivative of the metric tensor Gmat
          boxndm%hdot(1:3,1:3)=boxndm%hdot(1:3,1:3)*ihbox0(1:3,1:3)
          boxndm%Gdot(i,j) = Sum( boxndm%h(1:3,i)*boxndm%hdot(1:3,j) + boxndm%hdot(1:3,i)*boxndm%h(1:3,j) )
       END DO
    END DO

    ! Kinetic energy of the cell (Eq. 2.14 of Ref. [2])

    Boxndm%Kcell = 0.5d0*boxndm%wbox*Sum( boxndm%hDot(1:3,1:3)**2 )


    if(lEev) then
       unitE=erg2eV
       cunitE='  eV'
    else
       unitE=1.0
       cunitE=' erg'
    end if
    if(lwrtpr) write(uwrtpr,'(I7,D10.3,A,D21.12,A,a,f0.3,a)') 0,0.d0,'*Kcell = ',Boxndm%Kcell*unitE,cunitE, &
         '  (', 2.d0*Boxndm%Kcell/(9.d0*bk), ' K)'

    ! Nombre de thermostats de Hoover
    IF (nHoover.LT.0) nHoover=0
    IF (.NOT.lTHoover) nHoover=0
    if (.not.(allocated(boxndm%zhoover)))ALLOCATE(boxndm%zHoover(1:nHoover+1))


    IF (lTHoover) THEN

       ! Allocation des tableaux nécessaires au thermostat de Hoover
       IF (nHoover.LE.0) THEN
          WRITE(0,'(a)') 'Si le thermostat de Nosé-Hoover est utilisé (lTHoover=.true.), '
          WRITE(0,'(a)') 'le nombre de thermostats doit être supérieur à 0'
          WRITE(0,'(a,i0)') '  nHoover = ', nHoover
          STOP '< init_lpr >'
       END IF
       ALLOCATE(boxndm%zOld(1:nHoover), boxndm%zNew(1:nHoover), &
            &boxndm%zDot(1:nHoover), boxndm%wHoover(1:nHoover), &
            Boxndm%Khoover(1:nHoover))
       !ALLOCATE(UHoover(1:nHoover), UHoover_new(1:nHoover), UHoover_old(1:nHoover))

       ! Nombre de degrés de liberté pour le thermostat de Nosé-Hoover
       boxndm%gnose=dble(3*atpr%im_glob)

       ! Masse de chaque thermostat
       IF (boxndm%wNose.EQ.0) THEN
          ! On veut qu'une variation de la température de 10K corresponde à
          ! une variation de fNose de 1%
          boxndm%wNose = boxndm%gnose*bk*10.d0*tstep**2/1.d-2**2
       END IF
       boxndm%whoover(1)=boxndm%wNose
       DO i=2, nHoover
          boxndm%whoover(i)=boxndm%wNose/boxndm%gnose
       END DO

       DO i=1, nHoover
          Boxndm%Khoover(i) = 0.5d0*bk*Text
          boxndm%zHoover(i) = Sqrt( 2.d0*Boxndm%Khoover(i)/boxndm%whoover(i) )
          boxndm%zDot(i) = 0.d0
          boxndm%zOld(i) = boxndm%zHoover(i) - tstep*boxndm%zdot(i)
          !UHoover(i) = 0.d0
          !UHoover_old(i) = UHoover(i) - tstep*bk*Text*boxndm%zhoover(i)
       END DO
       !UHoover(1) = gNose*UHoover(1)
       boxndm%KNose = Sum( Boxndm%Khoover(1:nHoover) )
       !UNose = Sum( UHoover(1:nHoover) )
       boxndm%zHoover(nHoover+1)=0.d0

       IF(lwrtpr) WRITE(uwrtpr,*)
       IF(lwrtpr) WRITE(uwrtpr,'(a)') 'Thermostat de Nosé-Hoover (V2)'
       IF(lwrtpr) WRITE(uwrtpr,'(a)') "  -> l'énergie cinétique des atomes et de la boîte est thermalisée"
       IF(lwrtpr) WRITE(uwrtpr,'(a,g20.12)')'Masse de la boîte pour thermostat de Nosé-Hoover: wHoover=',boxndm%wNose
       IF(lwrtpr) WRITE(uwrtpr,'(a,g20.12)')'Nombre de degrés de liberté: gNose=',  boxndm%gNose
       IF(lwrtpr) WRITE(uwrtpr,'(a,i0)')    'Nombre de thermostats: nHoover=', nHoover
       IF(lwrtpr) WRITE(uwrtpr,'(a,f0.3,a)') 'Initialisation du thermostat de Nosé-Hoover pour la température ', &
            2.d0*boxndm%KNose/(bk*dble(nHoover)), ' K'
       IF(lwrtpr) WRITE(uwrtpr,'(I7,D10.3,A,D21.12,A,a,f0.3,a)') iteration,timel,'*KNose = ',boxndm%KNose*unitE,cunitE, &
            '  (', 2.d0*boxndm%KNose/(bk*nHoover), ' K)'
       IF(lwrtpr) WRITE(uwrtpr,*)
    ELSE
       boxndm%gNose=0.d0; boxndm%zHoover(1)=0.d0
    END IF
    if (.not.(allocated(sp)))then
       ALLOCATE(sp(1:3,1:atpr%imm), sdot(1:3,1:atpr%imm), sdot_new(1:3,1:atpr%imm))
       if (dmtype==24)     ALLOCATE(sfp(1:3,1:atpr%imm),spp(1:3,1:atpr%imm))
    end if

    ! Coordonnées réduites des atomes et leurs dérivées à l'instant initial
    if (dmtype.ne.15) then
       sp(:,1:atpr%im) = MatMul(boxndm%invh(:,:), atpr%xp(:,1:atpr%im) )
       sdot(:,1:atpr%im) = MatMul(boxndm%invh(:,:), atpr%vp(:,1:atpr%im) )

       ! Forces à l'instant initial
       CALL CalFo(sig,potist,atpr,celndm,boxndm%box_config,t_sigma=.true.,psc=psc)
       if (l2t)then
          if (i2t==1)  call calceloss (celndm,atpr)
       else
          if(ibrake.gt.0) call calceloss(celndm,atpr)
       end if


       !  Contrainte thermique à l'instant initial
       sigkine(:,:)=0.d0
       do ia = 1, atpr%im
          do j = 1,3
             sigkine(1:3,j) = sigkine(1:3,j) + cm(atpr%ityp(ia))*atpr%vp(1:3,ia)*atpr%vp(j,ia)
          enddo
       enddo
       sigkine(1:3,1:3) = boxndm%invVolu*sigkine(1:3,1:3)

#ifdef PARA
       if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
          call comm_space%sum(sigkine)
       end if

#endif
       if (lpcube) then
          pint=0.33333333333*(sigkine(1,1)+sigkine(2,2)+sigkine(3,3))
          sigkine=0
          do ic=1,3
             sigkine(ic,ic)=pint
          end do
          boxndm%ppot=(sig(1,1)+sig(2,2)+sig(3,3))/3.0
       end if

       ! Énergie cinétique des atomes à l'instant initial
       !       kine = 0.5d0*boxndm%volu*( sigKine(1,1) + sigKine(2,2) + sigKine(3,3) )

       ! Contrainte totale à l'instant initial
       sigtot = 0.5d0*(sigkine + Transpose(sigkine) + sig + Transpose(sig) )
    end if

    RETURN

  end subroutine initlpr

  !-----------------------------------------------
  !-----------------------------------------------
  !-----------------------------------------------

  subroutine pr1 (atpr,celndm,boxndm,psc)

    implicit none
    class(atom_config_d),target::atpr
    type(cell_config),target:: celndm
    type(box_config_lpr)::boxndm
    type(para_space_config)::psc

    real(double),dimension(3,3)::mf,mfi, grsig, hdot_new, hdot_last,forcebox,glanh
    REAL(double) :: diff, tdiff
    integer:: i,j,ia, iter,ic,im,ic2
    !real(double) , external ::  calcvol
    ! Parameter for Parrinello-Rahman self consistency loop
    REAL(double), parameter :: tol=1.0d-12        ! Tolerance for h convergency
    INTEGER, parameter :: max_Iter=1000            ! Maximal number of iterations in self-consistency loop
    real(double)::T1,kin1,tstepN,u1,u2,rga,rgah
    real(double), dimension(1:3) :: xprov
    real(double):: norme_de_fp, norme_de_vp, pscal,pint,gs3
    real(double), dimension(ntyp) :: aux
    integer,save::nstep=0
    real(double)::maxvploc,maxvp,nrmvp,sigrel(3,3)
    integer::ib
    real(double)::fbox

!!$    if (lpconxyz) then
!!$       sigtot(2,1)=0
!!$       sigtot(1,2)=0
!!$       sigtot(3,1)=0
!!$       sigtot(1,3)=0
!!$       sigtot(2,3)=0
!!$       sigtot(3,2)=0
!!$    end if
!!$
!!$    sigrel=sigtot
!!$    if (any(astarsig.eqv..true.)) then
!!$       call set_MP(boxndm,astarsig,sigtot,sigrel)
!!$    end if
!!$    write(uwrt,*) 'sigrel' ,sigrel(:,1)
!!$    write(uwrt,*) 'sigrel' ,sigrel(:,2)p
!!$    write(uwrt,*) 'sigrel' ,sigrel(:,3)
    select case(dmtype)
    case(24)
       if (lpcube) then
          fbox=(1/3.0)*boxndm%volu*boxndm%ppot*(1/(boxndm%h(1,1)**2+boxndm%h(2,2)**2+boxndm%h(3,3)**2))
       end if

       if (lpconxyz) then
          sig(2,1)=0
          sig(1,2)=0
          sig(3,1)=0
          sig(1,3)=0
          sig(2,3)=0
          sig(3,2)=0
       end if

       sigrel=sig
       if (any(astarsig.eqv..true.)) then
          call set_MP(boxndm,astarsig,sig,sigrel)
       end if

       !       write(uwrt,*)'IN',atpr%xp(1,1)
       ! Coordonnées réduites des atomes (au cas où elles ont été modifiées à l'extérieur)
       sp(:,1:atpr%im) = MatMul(boxndm%invh(:,:), atpr%xp(:,1:atpr%im) )
       select type (atpr)
       class is (atom_config_e)
          if (atpr%lxpp)       spp(:,1:atpr%im) = MatMul(boxndm%invh(:,:), atpr%xpp(:,1:atpr%im) )
       end select
       ! De même pour les vitesses au cas où, par exemple, on utilise le thermostat
       sdot(:,1:atpr%im) = MatMul(boxndm%invh(:,:), atpr%vp(:,1:atpr%im) )
       !       write(uwrt,*)'VEL',sdot(1,1),atpr%vp(1,1)
       DO ia=1, atpr%im
          sfp(:,ia)= MatMul( boxndm%invh(:,:), atpr%fp(:,ia) )
       END DO
       im=atpr%im
       ! 1/ Intégration de l'équation de mouvement
       aux(:ntyp) = tstep**2/(2.d0*cm(:ntyp))
       usdh = 1.d0/(2.d0*tstep)
       DO i=1, im
          xprov(:) = sp(:,i) + sdot(:,i)*tstep + sfp(:,i)*aux(atpr%iTyp(i))
          !            if (i==1) write(uwrt,*)'VPROV',(xprov(1) - spp(1,i))*usdh
          sdot(:,i) = (xprov(:) - spp(:,i))*usdh
          spp(:,i) = sp(:,i)
          sp(:,i) = xprov(:)
       END DO

       ! 2/ Renormalisation des vitesses par l'algorithme fire
       ! Puissance dissipée
       pScal = Sum( sdot(:,1:im)*sfp(:,1:im) )
#ifdef PARA
       if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
          call comm_space%sum(pscal)
       end if
#endif
       !         write(uwrt,*)'PSCAL',pscal,hdot(1,1)
       ! Modification du vecteur vitesse
       !         write(uwrt,*)'PSCAL',iteration,nstep,pscal,tstep
       if (pScal.gt.0) then
          ! Norme du vecteur force
          norme_de_fp = Sqrt( Sum( sfp(:,1:im)**2 ) )
          ! Norme du vecteur vitesse
          norme_de_vp = Sqrt( Sum( sdot(:,1:im)**2 ) )
          ! Nouveau vecteur vitesse
          sdot(:,1:im) = (1.d0-fire_alph)*sdot(:,1:im) + fire_alph*norme_de_vp/norme_de_fp*sfp(:,1:im)
          nStep = nStep + 1
          if (nStep.gt.nStepMin) then
             tstepN=min(tstep*finc,tstep_MM*tstep0)
             fire_alph=fire_alph*f_alph
          else
             tstepN=tstep
          end if
       else
          sdot(:,:)=0.
          tstepN=tstep*fdec
          fire_alph=alph_start
          nstep=0
       end if
       if (lpcube) then
          boxndm%hdot(2,1)=0
          boxndm%hdot(1,2)=0
          boxndm%hdot(3,1)=0
          boxndm%hdot(1,3)=0
          boxndm%hdot(2,3)=0
          boxndm%hdot(3,2)=0
          if (boxndm%hdot(1,1)*fbox<0) boxndm%hdot(:,:)=0
       else
	  forcebox(:,:)=MatMul( sigrel(:,:) - sigext(:,:), boxndm%invtrh(:,:) )*boxndm%volu
          do i = 1, 3
             do ic = 1, 3
                if (boxndm%hdot(ic,i)*forcebox(ic,i)<0) then
                   boxndm%hdot(ic,i)=0.
                end if
             end do
          end do
       end if

!!$       do i = 1, 3
!!$          do ic = 1, 3
!!$             if (boxndm%hdot(ic,i)*forcebox(ic,i)<0) then
!!$                boxndm%hdot(ic,i)=0.
!!$             end if
!!$          end do
!!$       end do

       !          write(uwrt,*)'hdot',hdot(1,1),tstep/(1*wBox)*forcebox(1,1)*ihbox0(1,1),invtrh(1,1),sigtot(1,1),tstep,tstepN
       if (lpcube) then
          do ic=1,3
             boxndm%hdot(ic,ic) = boxndm%h(ic,ic)*(( 1.d0  )* boxndm%hdot(ic,ic)  + tstep/(1.d0*boxndm%wBox)*fbox)
          end do
       else
          boxndm%hdot(:,:) = boxndm%hdot(:,:)*ihbox0(:,:)+ tstep/(1*boxndm%wBox)*forcebox(:,:)*ihbox0(:,:)          
       end if


!!$          write(uwrt,*)'hdot', boxndm%hdot(:,1)
!!$          write(uwrt,*)'hdot', boxndm%hdot(:,2)
!!$          write(uwrt,*)'hdot', boxndm%hdot(:,3)
       ! Tenseur h à l'instant t+dt
       boxndm%h(:,:) = boxndm%h(:,:) + boxndm%hdot(:,:)*tstep*ihbox0(:,:)
       !          write(uwrt,*)'FHdot',hdot(1,1),tstep,(1*wBox),forcebox(1,1)*ihbox0(1,1),h(1,1)

       !       write(uwrt,*)'MED',atpr%xp(1,1),sp(1,1),sfp(1,1)
       ! Coordonnées réelles à l'instant t+dt
       select type (atpr)
       class is (atom_config_e)
          if (atpr%lxpp)       atpr%xpp(:,1:atpr%im) = MatMul( boxndm%h, spp(:,1:atpr%im) )
       end select
       atpr%xp(:,1:atpr%im) = MatMul( boxndm%h, sp(:,1:atpr%im) )
       call updatebox(boxndm,boxndm%h)

       !#ifdef PARA
       atpr%vp(:,1:atpr%im) = MatMul( boxndm%h(:,:), sdot(:,1:atpr%im) )
       !#endif

       CALL ScaleBox(atpr,celndm,boxndm,psc)

       !#ifdef PARA
       !    atpr%vp(:,1:atpr%im) = MatMul( h(:,:), sdot(:,1:atpr%im) )
       sdot(:,1:atpr%im) = MatMul(boxndm%invh(:,:), atpr%vp(:,1:atpr%im) )
       !#endif
       !    write(uwrt,*)'MED2',atpr%xp(1,1),sp(1,1),sfp(1,1)
       ! Calcul des forces et des contraintes à l'instant t+dt
       CALL CalFo(sig,potist,atpr,celndm,boxndm%box_config,t_sigma=.true.,psc=psc)
       if (dmtype==24) tstep=tstepN
       if (lpcube) then
          pint=0.33333333333*(sigkine(1,1)+sigkine(2,2)+sigkine(3,3))
          sigkine=0
          do ic=1,3
             sigkine(ic,ic)=pint
          end do
          boxndm%ppot=(sig(1,1)+sig(2,2)+sig(3,3))/3.0
          fbox=(1/3.0)*boxndm%volu*boxndm%ppot*(1/(boxndm%h(1,1)**2+boxndm%h(2,2)**2+boxndm%h(3,3)**2))

       end if


       sigkine(:,:)=0.d0
       do ia = 1, atpr%im
          do j = 1,3
             sigkine(1:3,j) = sigkine(1:3,j) + cm(atpr%ityp(ia))*atpr%vp(1:3,ia)*atpr%vp(j,ia)
          enddo
       enddo
       sigkine(1:3,1:3) = boxndm%invVolu*sigkine(1:3,1:3)
#ifdef PARA
       if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
          call comm_space%sum(sigkine)
       end if
#endif
!!$          if (lpcube) then
!!$             pint=0.33333333333*(sigkine(1,1)+sigkine(2,2)+sigkine(3,3))
!!$             sigkine=0
!!$             do ic=1,3
!!$                sigkine(ic,ic)=pint
!!$             end do
!!$          end if

       sigtot = 0.5d0*(sigkine + Transpose(sigkine) + sig + Transpose(sig) )
!!$    write(uwrt,*) 'sigtot' ,sigtot(:,1)
!!$    write(uwrt,*) 'sigtot' ,sigtot(:,2)
!!$    write(uwrt,*) 'sigtot' ,sigtot(:,3)

       if ((any(astarsig.eqv..true.)).or.(any(sig0dir.ne.0))) then
       else
          thsig=maxval(abs(sig-sigext))   
       end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! LPR +LANGEVIN : evolution de VP en corrdonnées réelles pour éviter de se tromper dans les dimensions       
    case(88,16)
       if (lpcube) then
          fbox=(1/3.0)*boxndm%volu*boxndm%ppot*(1/(boxndm%h(1,1)**2+boxndm%h(2,2)**2+boxndm%h(3,3)**2))
       end if

       if (lpconxyz) then
          sigtot(2,1)=0
          sigtot(1,2)=0
          sigtot(3,1)=0
          sigtot(1,3)=0
          sigtot(2,3)=0
          sigtot(3,2)=0
       end if

       sigrel=sigtot
       if (any(astarsig.eqv..true.)) then
          call set_MP(boxndm,astarsig,sigtot,sigrel)
       end if
       call comm_space%barrier
       select type(atpr)
       class is (atom_config_e)
          boxndm%Kcell = 0.5d0*boxndm%wbox*Sum( boxndm%hDot(1:3,1:3)**2 )
          boxndm%Tempcell=Boxndm%Kcell*2./(sum(ihbox0)*bk)
          DO i=1, atpr%im

             rga=exp(-gamlt(atpr%ityp(i))*tstep/2)
             do ic=1,3
                call random_number(u1)
                call random_number(u2)
                atpr%Glangv(ic,i)=sqrt(-2.*log(u1))*cos(2.*pi*u2)   
                atpr%vp(ic,i) = atpr%vp(ic,i)*rga+ atpr%fp(ic,i)*tstep/(cm(atpr%ityp(i))*2)&
                     &+atpr%Glangv(ic,i)*sqrt(cm(atpr%ityp(i))*bk*text*(1-rga))/cm(atpr%ityp(i))
             end do


             ! Coordonnées réduites des atomes (au cas où elles ont été modifiées à l'extérieur)
             sp(:,1:atpr%im) = MatMul(boxndm%invh(:,:), atpr%xp(:,1:atpr%im) )
             ! De même pour les vitesses au cas où, par exemple, on utilise le thermostat
             sdot(:,1:atpr%im) = MatMul(boxndm%invh(:,:), atpr%vp(:,1:atpr%im) )

          END DO
          call sigkinetot(atpr,boxndm,sig,sigkine,sigtot)
          rgah=exp(-gamlg*gamprfact*tstep/2)
          if (myidsp==0) then
             do ic=1,3
                do ic2=1,3
                   call random_number(u1)
                   call random_number(u2)
                   glanh(ic,ic2)=sqrt(-2.*log(u1))*cos(2.*pi*u2)   
                end do
             end do
             if (lpcube) then
!!$                gs3=(glanh(1,1)+glanh(2,2)+glanh(3,3))/3.
!!$                glanh=0.
!!$                do ic=1,3
!!$                   glanh(ic,ic)=gs3
!!$                end do
!!$                boxndm%hdot(2,1)=0
!!$                boxndm%hdot(1,2)=0
!!$                boxndm%hdot(3,1)=0
!!$                boxndm%hdot(1,3)=0
!!$                boxndm%hdot(2,3)=0
!!$                boxndm%hdot(3,2)=0
!!$                do ic=1,3
!!$                   boxndm%hdot(ic,ic) =( boxndm%h(ic,ic)* rgah&
!!$                        &+ tstep/(2.d0*boxndm%wBox)*fbox &
!!$                        + (glanh(ic,ic)/boxndm%wbox)*sqrt(boxndm%wbox*bk*text*(1-rgah))  )*boxndm%h(ic,ic)*ihbox0(ic,ic)
!!$                end do
                write(uwrt,*)'langevin +lpcube =STOP'
                stop

             else


                boxndm%hdot(:,:) = (  boxndm%hdot(:,:)*rgah  &
                     + tstep/(2.d0*boxndm%wBox)*boxndm%volu*MatMul( sigtot(:,:) - sigext(:,:),boxndm%invtrh(:,:) ) &
                     + (glanh(:,:)/boxndm%wbox)*sqrt(boxndm%wbox*bk*text*(1-rgah))  )*ihbox0(:,:)
             end if
          end if

          call comm_space%bcast(0,boxndm%hdot)
          call comm_space%bcast(0,glanh)
          boxndm%Kcell = 0.5d0*boxndm%wbox*Sum( boxndm%hDot(1:3,1:3)**2 )
          boxndm%Tempcell=Boxndm%Kcell*2./(sum(ihbox0)*bk)

          sp(:,1:atpr%im) = sp(:,1:atpr%im) + sdot(:,1:atpr%im)*tstep

          ! Tenseur h à l'instant t+dt
          boxndm%h(:,:) = boxndm%h(:,:) + boxndm%hdot(:,:)*tstep*ihbox0(:,:)

          ! Coordonnées réelles à l'instant t+dt
          select type (atpr)
          class is (atom_config_e)
             if (atpr%lxpp)          atpr%xpp(:,1:atpr%im) = atpr%xp(:,1:atpr%im)
          end select
          call updatebox(boxndm,boxndm%h)
          atpr%xp(:,1:atpr%im) = MatMul( boxndm%h, sp(:,1:atpr%im) )

          atpr%vp(:,1:atpr%im) = MatMul( boxndm%h(:,:), sdot(:,1:atpr%im) ) ! retour à vp car transfert d'atomes  dans scalebox en PARA
          Boxndm%Kcell = 0.5d0*boxndm%wbox*Sum( boxndm%hDot(1:3,1:3)**2 )
          boxndm%Tempcell=Boxndm%Kcell*2./(sum(ihbox0)*bk)


          CALL ScaleBox(atpr,celndm,boxndm,psc)
          call caltabtC(celndm,atpr,lperiod,boxndm,lchktrav=.true.)

          ! Calcul des forces et des contraintes à l'instant t+dt
          CALL CalFo(sig,potist,atpr,celndm,boxndm%box_config,t_sigma=.true.,psc=psc)
          call sigkinetot(atpr,boxndm,sig,sigkine,sigtot)

          do i=1,atpr%im
             rga=exp(-gamlt(atpr%ityp(i))*tstep/2)
             do ic=1,3
                atpr%vp(ic,i) = atpr%vp(ic,i)*rga+ atpr%fp(ic,i)*tstep/(cm(atpr%ityp(i))*2)&
                     &+atpr%Glangv(ic,i)*sqrt(cm(atpr%ityp(i))*bk*text*(1-rga))/cm(atpr%ityp(i))
             end do
          end do
          !    sdot(:,1:atpr%im) = MatMul(boxndm%invh(:,:), atpr%vp(:,1:atpr%im) )

          call sigkinetot(atpr,boxndm,sig,sigkine,sigtot)
          if (lpcube) then
             !                gs3=(glanh(1,1)+glanh(2,2)+glanh(3,3))/3.
             !                glanh=0.
             !                do ic=1,3
             !                   glanh(ic,ic)=gs3
             !                end do
             boxndm%hdot(2,1)=0
             boxndm%hdot(1,2)=0
             boxndm%hdot(3,1)=0
             boxndm%hdot(1,3)=0
             boxndm%hdot(2,3)=0
             boxndm%hdot(3,2)=0
             do ic=1,3
                boxndm%hdot(ic,ic) =( boxndm%h(ic,ic)* rgah&
                     &+ tstep/(2.d0*boxndm%wBox)*fbox &
                     + (glanh(ic,ic)/boxndm%wbox)*sqrt(boxndm%wbox*bk*text*(1-rgah))  )*boxndm%h(ic,ic)*ihbox0(ic,ic)
             end do

          else


             boxndm%hdot(:,:) = (  boxndm%hdot(:,:)*rgah  &
                  + tstep/(2.d0*boxndm%wBox)*boxndm%volu*MatMul( sigtot(:,:) - sigext(:,:),boxndm%invtrh(:,:) ) &
                  + (glanh(:,:)/boxndm%wbox)*sqrt(boxndm%wbox*bk*text*(1-rgah))  )*ihbox0(:,:)
          end if

          !!          boxndm%hdot(:,:) = (boxndm%hdot(:,:)*rgah  &
          !!               + tstep/(2.d0*boxndm%wBox)*boxndm%volu*MatMul( sigtot(:,:) - sigext(:,:), boxndm%invtrh(:,:) ) &
          !!              + (glanh(:,:)/boxndm%wbox)*sqrt(boxndm%wbox*bk*text*(1-rgah))  )*ihbox0(:,:)
          ! Estimation de la dérivée du tenseur Gmat à l'instant t+dt
          DO i=1, 3
             DO j=1, 3
                boxndm%Gdot(i,j) = Sum( boxndm%h(1:3,i)*boxndm%hdot(1:3,j) + boxndm%hdot(1:3,i)*boxndm%h(1:3,j) )
             END DO
          END DO

          grsig = boxndm%volu * MatMul(boxndm%invh, MatMul( sigext, boxndm%invtrh) )
          boxndm%tension = boxndm%invVolu0*MatMul( MatMul( boxndm%h0, grsig), boxndm%trh0 )

          ! Énergie potentielle de la cellule (Eq. 2.25, Ref.2)
          boxndm%Ucell = boxndm%volu0*Sum( boxndm%tension(1:3,1:3) * boxndm%epsi(1:3,1:3) )

          ! Énergie cinétique de la cellule (Eq. 2.14, Ref.2)
          boxndm%Kcell = 0.5d0*boxndm%wbox*Sum( boxndm%hDot(1:3,1:3)**2 )
          boxndm%Tempcell=Boxndm%Kcell*2./(sum(ihbox0)*bk)
          boxndm%EcellPR = boxndm%Kcell + boxndm%Ucell

       end select
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    case(22,8)
       if (lpcube) then
          fbox=(1/3.0)*boxndm%volu*boxndm%ppot*(1/(boxndm%h(1,1)**2+boxndm%h(2,2)**2+boxndm%h(3,3)**2))
       end if


       if (dmtype==22) then
          if (lpconxyz) then
             sig(2,1)=0
             sig(1,2)=0
             sig(3,1)=0
             sig(1,3)=0
             sig(2,3)=0
             sig(3,2)=0
          end if

          sigrel=sig
          if (any(astarsig.eqv..true.)) then
             call set_MP(boxndm,astarsig,sig,sigrel)
          end if

          do i = 1, atpr%im
             do ic = 1, 3
                if (atpr%vp(ic,i)*atpr%fp(ic,i)<0) then
                   atpr%vp(ic,i)=0.
                end if
             end do
          end do
          forcebox(:,:)=MatMul( sigrel(:,:) - sigext(:,:), boxndm%invtrh(:,:) )
          if (lpcube) then
             boxndm%hdot(2,1)=0
             boxndm%hdot(1,2)=0
             boxndm%hdot(3,1)=0
             boxndm%hdot(1,3)=0
             boxndm%hdot(2,3)=0
             boxndm%hdot(3,2)=0
             if (boxndm%hdot(1,1)*fbox<0) boxndm%hdot(:,:)=0
          else

             do i = 1, 3
                do ic = 1, 3
                   if (boxndm%hdot(ic,i)*forcebox(ic,i)<0) then
                      boxndm%hdot(ic,i)=0.
                   end if
                end do
             end do
          end if
       else
          if (lpconxyz) then
             sigtot(2,1)=0
             sigtot(1,2)=0
             sigtot(3,1)=0
             sigtot(1,3)=0
             sigtot(2,3)=0
             sigtot(3,2)=0
          end if
          if (lpcube) then
             boxndm%hdot(2,1)=0
             boxndm%hdot(1,2)=0
             boxndm%hdot(3,1)=0
             boxndm%hdot(1,3)=0
             boxndm%hdot(2,3)=0
             boxndm%hdot(3,2)=0
          end if


          sigrel=sigtot
          if (any(astarsig.eqv..true.)) then
             call set_MP(boxndm,astarsig,sigtot,sigrel)
          end if

       end if

       !    case(8)
       ! Coordonnées réduites des atomes (au cas où elles ont été modifiées à l'extérieur)
       sp(:,1:atpr%im) = MatMul(boxndm%invh(:,:), atpr%xp(:,1:atpr%im) )
       ! De même pour les vitesses au cas où, par exemple, on utilise le thermostat
       sdot(:,1:atpr%im) = MatMul(boxndm%invh(:,:), atpr%vp(:,1:atpr%im) )
       ! Dérivée des coordonnées réduites des atomes à l'instant t+dt/2
       mf(:,:) = -0.5d0*tstep*MatMul(boxndm%invGmat,boxndm%Gdot)
       DO i=1, 3
          mf(i,i) = 1.d0 - 0.5d0*tstep*boxndm%zHoover(1) + mf(i,i)
       END DO
       DO ia=1, atpr%im
          sdot(:,ia) = MatMul( mf(:,:), sdot(:,ia) ) &
               + tstep/(2.d0*cm(atpr%ityp(ia))) * MatMul( boxndm%invh(:,:), atpr%fp(:,ia) )
       END DO
       ! Dérivée du tenseur h à l'instant t+dt/2
       if (lpcube) then
          IF (lpcon2.EQV..true.) THEN    ! On ajoute une force de friction
             do ic=1,3
                boxndm%hdot(ic,ic) = boxndm%h(ic,ic)*(( 1.d0 - 0.5d0*tstep*boxndm%zHoover(1)- 0.5d0/tbox )* boxndm%hdot(ic,ic) &
                     + tstep/(2.d0*boxndm%wBox)*fbox)
             end do
          ELSE        ! Équation sans force de friction supplémentaire
             do ic=1,3
                boxndm%hdot(ic,ic) = boxndm%h(ic,ic)*(( 1.d0 - 0.5d0*tstep*boxndm%zHoover(1) )* boxndm%hdot(ic,ic) &
                     + tstep/(2.d0*boxndm%wBox)*fbox)
             end do
          END IF

       else

          IF (lpcon2.EQV..true.) THEN    ! On ajoute une force de friction
             boxndm%hdot(:,:) = ( 1.d0 - 0.5d0*tstep*boxndm%zHoover(1) - 0.5d0/tbox )* boxndm%hdot(:,:)*ihbox0(:,:) &
                  + tstep/(2.d0*boxndm%wBox)*boxndm%volu*MatMul( &
                  &sigrel(:,:) - sigext(:,:), boxndm%invtrh(:,:) )*ihbox0(:,:)
          ELSE        ! Équation sans force de friction supplémentaire
             boxndm%hdot(:,:) = ( 1.d0 - 0.5d0*tstep*boxndm%zHoover(1) )* boxndm%hdot(:,:)*ihbox0(:,:) &
                  + tstep/(2.d0*boxndm%wBox)*boxndm%volu*MatMul(&
                  &sigrel(:,:) - sigext(:,:), boxndm%invtrh(:,:) )*ihbox0(:,:)
          END IF
       end if
       !          write(uwrt,*)'hdot',hdot
       ! Coordonnées réduites des atomes à l'instant t+dt
       sp(:,1:atpr%im) = sp(:,1:atpr%im) + sdot(:,1:atpr%im)*tstep

       ! Tenseur h à l'instant t+dt
       boxndm%h(:,:) = boxndm%h(:,:) + boxndm%hdot(:,:)*tstep*ihbox0(:,:)

       ! Coordonnées réelles à l'instant t+dt
       select type (atpr)
       class is (atom_config_e)
          if (atpr%lxpp)       atpr%xpp(:,1:atpr%im) = atpr%xp(:,1:atpr%im)
       end select
       atpr%xp(:,1:atpr%im) = MatMul( boxndm%h, sp(:,1:atpr%im) )
       call updatebox(boxndm,boxndm%h)
!!$    invVolu = 1.d0/boxndm%volu
!!$    trh=Transpose(h)                            ! Matrices associées à h
!!$    Gmat = MatMul(trh,h)
!!$    CALL MatInv(Gmat,invGmat)
!!$    call MatInv(h,invh)
!!$    invtrh = Transpose(invh)

#ifdef PARA
       atpr%vp(:,1:atpr%im) = MatMul( boxndm%h(:,:), sdot(:,1:atpr%im) )
       !    sdot(:,1:atpr%im) = MatMul(invh(:,:), atpr%vp(:,1:atpr%im) )
#endif

       CALL ScaleBox(atpr,celndm,boxndm,psc)

#ifdef PARA
       !    atpr%vp(:,1:atpr%im) = MatMul( h(:,:), sdot(:,1:atpr%im) )
       sdot(:,1:atpr%im) = MatMul(boxndm%invh(:,:), atpr%vp(:,1:atpr%im) )
#endif
       !    write(uwrt,*)'MED2',atpr%xp(1,1),sp(1,1),sfp(1,1)
       ! Calcul des forces et des contraintes à l'instant t+dt
       !           block
       !             real(double),allocatable::normvp (:)
       !             allocate(normvp(atpr%im))
       if (lmaxvp) then
          do ib=1,atpr%im
             nrmvp=norm2(atpr%vp(:,ib))
             if(nrmvp.ge.vplim) then
                atpr%vp(:,ib)=atpr%vp(:,ib)*5d6/nrmvp
             end if
          end do
       end if
!!$             maxvp=maxval(normvp)
!!$             call comm_space%max(maxvp)
!!$             if (rang==0) write(uwrt,*)'MAXVP', maxvp
       !           end block

       CALL CalFo(sig,potist,atpr,celndm,boxndm%box_config,t_sigma=.true.,psc=psc)
       if (lpcube) then
          pint=0.33333333333*(sigkine(1,1)+sigkine(2,2)+sigkine(3,3))
          sigkine=0
          do ic=1,3
             sigkine(ic,ic)=pint
          end do
          boxndm%ppot=(sig(1,1)+sig(2,2)+sig(3,3))/3.0
          fbox=(1/3.0)*boxndm%volu*boxndm%ppot*(1/(boxndm%h(1,1)**2+boxndm%h(2,2)**2+boxndm%h(3,3)**2))

       end if

       if (l2t)then
          if (i2t==1)  call calceloss(celndm,atpr)
       else
          if(ibrake.gt.0) call calceloss(celndm,atpr)
       end if
       if (lTberendsen) call calfoberend(atpr)

       ! Calcul de la viscosité à l'instant ...
       DO i=1, nHoover
          boxndm%zNew(i) = boxndm%zOld(i) + 2.d0*boxndm%zDot(i)*tstep    ! ... t+dt
          boxndm%zOld(i) = boxndm%zHoover(i)                      ! ... t
          boxndm%zHoover(i) = boxndm%zNew(i)                      ! ... t+dt
       END DO

       ! Estimation de la contrainte totale à l'instant t+dt
       ! (la contrainte cinétique est calculée à l'instant t
       ! et la contrainte potentielle à l'instant t+dt)
       sigtot = 0.5d0*(sigkine + Transpose(sigkine) + sig + Transpose(sig) )
       if (lpconxyz) then
          sigtot(2,1)=0
          sigtot(1,2)=0
          sigtot(3,1)=0
          sigtot(1,3)=0
          sigtot(2,3)=0
          sigtot(3,2)=0
       end if

       sigrel=sigtot
       if (any(astarsig.eqv..true.)) then
          call set_MP(boxndm,astarsig,sigrel,sigrel)
       end if

       ! Estimation de la dérivée du tenseur h à l'instant t+dt
       if (lpcube) then
          IF (lpcon2.EQV..true.) THEN    ! On ajoute une force de friction
             hdot_new=0.
             do ic=1,3
                hdot_new(ic,ic) = boxndm%h(ic,ic)*(( 1.d0 - 0.5d0*tstep*boxndm%zHoover(1)- 0.5d0/tbox )* boxndm%hdot(ic,ic) &
                     + tstep/(2.d0*boxndm%wBox)*fbox)
             end do
          ELSE        ! Équation sans force de friction supplémentaire
             hdot_new=0.
             do ic=1,3
                hdot_new(ic,ic) = boxndm%h(ic,ic)*(( 1.d0 - 0.5d0*tstep*boxndm%zHoover(1) )* boxndm%hdot(ic,ic) &
                     + tstep/(2.d0*boxndm%wBox)*fbox)
             end do
          END IF

       else
          IF (lpcon2.EQV..true.) THEN    ! On ajoute une force de friction
             hdot_new(:,:) = 1.d0/(1.d0+0.5d0*tstep*boxndm%zHoover(1)+ 0.5d0/tbox )*( boxndm%hdot(:,:)*ihbox0(:,:) &
                  + tstep/(2.d0*boxndm%wBox)*boxndm%volu*MatMul(sigrel(:,:)-sigext(:,:),boxndm%invtrh(:,:)))*ihbox0(:,:)
          ELSE        ! Équation sans force de friction supplémentaire
             hdot_new(:,:) = 1.d0/( 1.d0 + 0.5d0*tstep*boxndm%zHoover(1) )*( boxndm%hdot(:,:)*ihbox0(:,:) &
                  + tstep/(2.d0*boxndm%wBox)*boxndm%volu*MatMul(sigrel(:,:)-sigext(:,:),boxndm%invtrh(:,:)))*ihbox0(:,:)
          END IF
       end if

!!$       
!!$       ! Estimation de la dérivée du tenseur Gmat à l'instant t+dt
       DO i=1, 3
          DO j=1, 3
             boxndm%Gdot(i,j) = Sum(boxndm%h(1:3,i)*hdot_new(1:3,j) + hdot_new(1:3,i)*boxndm%h(1:3,j) )
          END DO
       END DO

       ! Cycle autocohérent (hdot -> sdot -> sigkine -> hdot)
       DO iter=1, Max_Iter
          ! Valeurs de la dernière itération du cycle d'autocohérence
          hdot_last(:,:) = hdot_new(:,:)*ihbox0(:,:)

          ! Dérivée des coordonnées réduites des atomes à l'instant t+dt
          mf(:,:) = 0.5d0*tstep*MatMul(boxndm%invGmat,boxndm%Gdot)
          DO i=1, 3
             mf(i,i) = 1.d0 + 0.5d0*tstep*boxndm%zHoover(1) + mf(i,i)
          END DO
          CALL MatInv(mf, mfi)
          DO ia=1, atpr%im
             sdot_new(:,ia) = MatMul( mfi(:,:), sdot(:,ia)  &
                  + tstep/(2.d0*cm(atpr%ityp(ia))) * MatMul(boxndm%invh(:,:), atpr%fp(:,ia) ) )
          END DO

          ! Vitesse des atomes à l'instant t+dt
          atpr%vp(:,1:atpr%im) = MatMul(boxndm%h(:,:), sdot_new(:,1:atpr%im) )
          if (lmaxvp) then
             do ib=1,atpr%im
                nrmvp=norm2(atpr%vp(:,ib))
                if(nrmvp.ge.5d6) then
                   atpr%vp(:,ib)=atpr%vp(:,ib)*5d6/nrmvp
                end if
             end do
          end if

          !  Contrainte thermique à l'instant t+dt
          sigkine(:,:)=0.d0
          do ia = 1, atpr%im
             do j = 1,3
                sigkine(1:3,j) = sigkine(1:3,j) + cm(atpr%ityp(ia))*atpr%vp(1:3,ia)*atpr%vp(j,ia)
             enddo
          enddo
          sigkine(1:3,1:3) = boxndm%invVolu*sigkine(1:3,1:3)
#ifdef PARA

          if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
             call comm_space%sum(sigkine)
          end if
#endif
          if (lpcube) then
             pint=0.33333333333*(sigkine(1,1)+sigkine(2,2)+sigkine(3,3))
             sigkine=0
             do ic=1,3
                sigkine(ic,ic)=pint
             end do
          end if

          ! Contrainte totale à l'instant t+dt
          sigtot = 0.5d0*(sigkine + Transpose(sigkine) + sig + Transpose(sig) )
          if (lpcube) then
             IF (lpcon2.EQV..true.) THEN    ! On ajoute une force de friction
                hdot_new=0
                do ic=1,3
                   hdot_new(ic,ic) = boxndm%h(ic,ic)*(( 1.d0 - 0.5d0*tstep*boxndm%zHoover(1)- 0.5d0/tbox )* boxndm%hdot(ic,ic) &
                        + tstep/(2.d0*boxndm%wBox)*fbox)
                end do

             ELSE        ! Équation sans force de friction supplémentaire
                hdot_new=0
                do ic=1,3
                   hdot_new(ic,ic) = boxndm%h(ic,ic)*(( 1.d0 - 0.5d0*tstep*boxndm%zHoover(1) )* boxndm%hdot(ic,ic) &
                        + tstep/(2.d0*boxndm%wBox)*fbox)
                end do
             END IF

          else
             ! Dérivée du tenseur h à l'instant t+dt
             IF (lpcon2.EQV..true.) THEN    ! On ajoute une force de friction
                hdot_new(:,:) = 1.d0/( 1.d0 + 0.5d0*tstep*boxndm%zHoover(1) + 0.5d0/tbox )*(boxndm%hdot(:,:)*ihbox0(:,:) &
                     + tstep/(2.d0*boxndm%wBox)*boxndm%volu*MatMul( sigtot(:,:)-sigext(:,:),boxndm%invtrh(:,:)))*ihbox0(:,:)
             ELSE        ! Équation sans force de friction supplémentaire
                hdot_new(:,:) = 1.d0/( 1.d0 + 0.5d0*tstep*boxndm%zHoover(1) )*(boxndm%hdot(:,:) &
                     + tstep/(2.d0*boxndm%wBox)*boxndm%volu*MatMul(sigtot(:,:)-sigext(:,:),boxndm%invtrh(:,:)))*ihbox0(:,:)
             END IF
          end if


          ! Dérivée du tenseur Gmat à l'instant t+dt
          DO i=1, 3
             DO j=1, 3
                boxndm%Gdot(i,j) = Sum( boxndm%h(1:3,i)*hdot_new(1:3,j) + hdot_new(1:3,i)*boxndm%h(1:3,j) )
             END DO
          END DO

          ! Vérifie l'autocohérence de h
          diff = Sum( abs( hdot_new(1:3,1:3) - hdot_last(1:3,1:3) ) )
          tdiff = Sum( abs( hdot_last(1:3,1:3) ) )
          !          write(uwrt,*)'hdotnew',hdot_new
          !          write(uwrt,*)'hdotlast',hdot_last
          !          write(uwrt,*)'conv',diff/tdiff,hdot_new(1,1),hdot_last(1,1)
          !          if (tdiff .eq. 0.0d0) then
          !             if (diff .LE. tol) exit
          !          else
          if (abs(diff/tdiff) .LE. tol) exit
          !          endif

       END DO

       ! Valeurs convergées de hdot et sdot à l'instant t+dt
       boxndm%hdot(:,:) = hdot_new(:,:)*ihbox0(:,:)
       sdot(:,:) = sdot_new(:,:)
       IF (iter.GE.Max_Iter) THEN
          WRITE(0,'(a,i0,a)') 'Maximal number of iterations (', Max_Iter, &
               ') in Parrinello-Rahman self consistency loop has been reached'
          WRITE(0,*)
          WRITE(0,'(a)') 'Last hdot proposed:'
          WRITE(0,'(a,3(f0.5,1x))') ' hdot(1,1:3) = ', 1e8*hdot_last(1,1:3)
          WRITE(0,'(a,3(f0.5,1x))') ' hdot(2,1:3) = ', 1e8*hdot_last(2,1:3)
          WRITE(0,'(a,3(f0.5,1x))') ' hdot(3,1:3) = ', 1e8*hdot_last(3,1:3)
          WRITE(0,*)
          WRITE(0,'(a)') 'New hdot proposed:'
          WRITE(0,'(a,3(f0.5,1x))') ' hdot(1,1:3) = ', 1e8*hdot_new(1,1:3)
          WRITE(0,'(a,3(f0.5,1x))') ' hdot(2,1:3) = ', 1e8*hdot_new(2,1:3)
          WRITE(0,'(a,3(f0.5,1x))') ' hdot(3,1:3) = ', 1e8*hdot_new(3,1:3)
          STOP '< PRNose >'
       END IF

       ! Énergie cinétique des atomes à l'instant t+dt
       kine = 0.5d0*boxndm%volu*( sigKine(1,1) + sigKine(2,2) + sigKine(3,3) )

       ! Déformation (Eq. 2.16, Ref.2)
       boxndm%epsi=0.5d0*MatMul( MatMul( boxndm%invtrh0, boxndm%Gmat ), boxndm%invh0 )
       DO i=1, 3
          boxndm%epsi(i,i) = boxndm%epsi(i,i) - 1.d0
       END DO
       ! Tension thermodynamique (Eq. 2.22 et 2.26, Ref.2)
       grsig = boxndm%volu * MatMul(boxndm%invh, MatMul( sigext, boxndm%invtrh) )
       boxndm%tension = boxndm%invVolu0*MatMul( MatMul( boxndm%h0, grsig), boxndm%trh0 )

       ! Énergie potentielle de la cellule (Eq. 2.25, Ref.2)
       boxndm%Ucell = boxndm%volu0*Sum( boxndm%tension(1:3,1:3) * boxndm%epsi(1:3,1:3) )

       ! Énergie cinétique de la cellule (Eq. 2.14, Ref.2)
       Boxndm%Kcell = 0.5d0*boxndm%wbox*Sum(boxndm%hDot(1:3,1:3)**2 )
       boxndm%EcellPR = boxndm%Kcell + boxndm%Ucell
       IF (lTHoover) THEN
          ! Dérivée de la viscosité et énergie cinétique du thermostat
          boxndm%zDot(1) = (2.d0*(kine + Boxndm%Kcell) - boxndm%gNose*bk*Text)/boxndm%wHoover(1) &
               - boxndm%zHoover(2)*boxndm%zHoover(1)
          Boxndm%Khoover(1) = 0.5d0*boxndm%whoover(1)*boxndm%zHoover(1)**2         ! t+dt
          DO i=2, nHoover
             boxndm%zDot(i) = ( 2.d0*Boxndm%Khoover(i-1) - bk*Text )/boxndm%whoover(i) &
                  - boxndm%zhoover(i+1)*boxndm%zhoover(i)
             Boxndm%Khoover(i) = 0.5d0*boxndm%whoover(i)*boxndm%zhoover(i)**2
          END DO
          boxndm%KNose = Sum( Boxndm%Khoover(1:nHoover) )
          ! Énergie potentielle du thermostat
          !DO i=1, nHoover
          !   UHoover_new(i) = UHoover_old(i) + 2.d0*tstep*bk*Text*boxndm%zhoover(i)   ! t+dt
          !   UHoover_old(i) = UHoover(i)              ! t
          !   UHoover(i) = UHoover_new(i)              ! t+dt
          !END DO
          !UHoover(1) = gNose*UHoover(1)               ! t+dt
          !UNose = Sum( UHoover(1:nHoover) )           ! t+dt
          !ENose = KNose + UNose                       ! t+dt
       END IF
    end select
    !    write(uwrt,*)'OUT',atpr%xp(1,1),tstep
    call calctemp(T1,kin1,atpr,celndm)
  end subroutine pr1



  subroutine set_MP(box,astarsig,sig,sigr)
    class(box_config)::box
    logical::astarsig(3)
    real(double)::sig(3,3),sigr(3,3)
    real(double)::sigt(3,3),sigt0(3,3)
    real(double)::MP(3,3),tMP(3,3),atp(3,3),vp(3)
    integer::ichk,ic,i
    ichk=0
    !    do i=1,3
    !       write(uwrt,*)'SIG',sig(1:3,i)
    !    end do

    if (astarsig(1)) then
       ichk=ichk+1
       atp(:,1)=box%as(:,1)/norm2(box%as(:,1))
       atp(:,2)=box%at(:,2)/norm2(box%at(:,2))
    end if
    if (astarsig(2)) then
       ichk=ichk+1
       atp(:,1)=box%as(:,2)/norm2(box%as(:,2))
       atp(:,2)=box%at(:,3)/norm2(box%at(:,3))
    end if
    if (astarsig(3)) then
       ichk=ichk+1
       atp(:,1)=box%as(:,3)/norm2(box%as(:,3))
       atp(:,2)=box%at(:,1)/norm2(box%at(:,1))
    end if
    if (ichk.ne.1) then
       write(uwrt,*)'ichk',ichk
       call arret_ndm
    end if
    call vectprod(atp(:,1),atp(:,2),atp(:,3))
    atp(:,3)=atp(:,3)/norm2(atp(:,3))
    !   write(uwrt,*)'NORMatp ',norm2(atp(:,1)),norm2(atp(:,2)),norm2(atp(:,3))

    call mattrp(atp,tMP)
    sigt=matmul(tMP,sig)
    sigt=matmul(sigt,atp)
    sigt0(:,:)=0
    sigt0(1,1)=sigt(1,1)
    thsig=abs(sigt(1,1))
    !   write(uwrt,*)'thsig',thsig*unitP
    sigr=matmul(atp,sigt0)
    sigr=matmul(sigr,tMP)
    !   do i=1,3
    !      write(uwrt,*)'SIGR',sigr(1:3,i)
    !   end do
    return
  end subroutine set_MP
!!$  subroutine set_MP2(box,sig0dir,sigt,sigr)
!!$    class(box_config_lpr)::box
!!$    real(double)::sig(3,3),sigr(3,3),sig0dir(3)
!!$    real(double)::sigt(3,3),sigt0(3,3)
!!$    real(double)::MP(3,3),tMP(3,3),atp(3,3),vp(3)
!!$    integer::ichk=0,ic,i
!!$    do i=1,3
!!$       write(uwrt,*)'SIG',sig(1:3,i)
!!$    end do
!!$    
!!$    atp(:,1)=sig0dir(:)/norm2(sig0dir)
!!$    call vectprod(atp(:,1),box%at(:,3),atp(:,2))
!!$    atp(:,2)=atp(:,2)/norm2(atp(:,2))
!!$    call vectprod(atp(:,1),atp(:,2),atp(:,3))
!!$    atp(:,3)=atp(:,3)/norm2(atp(:,3))
!!$    write(uwrt,*)'NORMatp ',norm2(atp(:,1)),norm2(atp(:,2)),norm2(atp(:,3))
!!$
!!$    call mattrp(atp,tMP)
!!$    sigt=matmul(tMP,sig)
!!$    sigt=matmul(sigt,atp)
!!$    sigt0(:,:)=0
!!$    sigt0(1,1)=sigt(1,1)
!!$    write(uwrt,*)'thsig',thsig*unitP
!!$    sigr=matmul(atp,sigt0)
!!$    sigr=matmul(sigr,tMP)
!!$    do i=1,3
!!$       write(uwrt,*)'SIGR',sigr(1:3,i)
!!$    end do
!!$    return
!!$    
!!$  end subroutine set_MP2
!!$    
!!$       

end module Parrinello_Rahman !Parrinello_Rahman





!******************
