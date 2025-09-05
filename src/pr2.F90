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
  !   Du coup, la matrice h0 définissant l'état de référence n'apparaît plus
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
  USE gen_com_m, ONLY:ecellpr,kcell,kine,knose,lpcon2,lthoover,nhoover,sigext,ucell,erg2ev,&
       &kcell,kine,knose,leev,lthoover,lucell,nhoover,timel,wboxf,wnose,zhoover, ihbox0,tbox, bk,&
       &potist,sig,sigtot,text,tstep,iteration,potist,rang,sig,text,sigkine,lpcube,&
       &pi,l2t,ltberendsen,lperiod,lspaceNDM,h0,dmtype,usdh,llangevin,gamlg,gamprfact,unitP,&
       & lmaxvp,vplim
  
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
  USE Mat_utils_mod,only:  MatInv
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

  implicit none
!!$  ! Vecteurs de la boîte et leurs dérivées
!!$  real(double), dimension(3,3), save , private :: h, hDot
!!$  real(double), dimension(3,3), save , private :: trh, invh, invtrh, Gmat, invGmat, Gdot
!!$  real(double), save, private :: invVolu

  ! Coordonnées réduites des atomes et leurs dérivées
  real(double), allocatable :: sp(:,:), sdot(:,:), sdot_new(:,:),sfp(:,:),spp(:,:)

  real(double)::kinx,tempx,pre,pint

  ! Variable associée au thermostat de Nosé-Hoover
  !  (zHoover est défini dans gen_com_m.F90)
  REAL(double), dimension(:), allocatable, save, private :: wHoover   ! Poids associé au thermostat de Hoover
  REAL(double), dimension(:), allocatable, save, private :: zOld, zNew, zDot !  Viscosité et dérivée
  REAL(double), dimension(:), allocatable, save, private :: KHoover
  !REAL(double), dimension(:), allocatable, save, private :: UHoover, UHoover_new, UHoover_old Énergies

  ! Nombre de degrés de liberté
  REAL(double), save, private :: gNose
  !  real(double)::wbox
  ! Variables uniquement nécessaires au calcul de l'énergie potentielle de la
  ! boîte
  real(double), dimension(3,3) ::trh0,invh0,invtrh0,epsi, tension
  real(double) ::volu0, invVolu0
  REAL(double) ::  fire_alph
  INTEGER :: fire_nstep,ic

  real(double)::TInitBox

contains

  subroutine initlpr (atpr,celndm,boxndm,psc)

    implicit none
    type(para_space_config)::psc
    type(box_config_lpr)::boxndm
    class(atom_config_d)::atpr
    type(cell_config):: celndm


    INTEGER :: ia, i, j
    !real(double), external :: calcvol
    real(double):: unitE,v1,z1,z2,tempcell
    character*5 :: cunitE

    if (dmtype==24) then
       CALL init_trempe_fire(tstep, fire_nstep, fire_alph)
    END IF

    IF(RANG==0) WRITE(6,*)
    if (dmtype.ne.15) then
       if (llangevin) then 
          IF(RANG==0) WRITE(6,*) 'Algorithme deP cst Langevin Parrinello-Rahman '
       else
          IF(RANG==0) WRITE(6,*) 'Algorithme de Parrinello-Rahman (V2)'
       end if
       IF(RANG==0) WRITE(6,'(a)') '  -> la vitesse de la boîte ne prend pas en compte la dérivée du tenseur h à t=0'
       IF(RANG==0) WRITE(6,*)
    end if

    IF (lUcell) THEN
       IF(RANG==0) WRITE(6,'(a)') "Repère de référence pour Parrinello-Rahman  (A):"
       IF(RANG==0) WRITE(6,'(a,3(f0.5,1x))') ' h0(1:3,1) = ', 1e8*h0(1:3,1)
       IF(RANG==0) WRITE(6,'(a,3(f0.5,1x))') ' h0(1:3,2) = ', 1e8*h0(1:3,2)
       IF(RANG==0) WRITE(6,'(a,3(f0.5,1x))') ' h0(1:3,3) = ', 1e8*h0(1:3,3)
       IF(RANG==0) WRITE(6,*)
    ELSE
       h0 = boxndm%at
    END IF
    IF(RANG==0) WRITE(6,'(a)') "Repère actuel  (A):"
    IF(RANG==0) WRITE(6,'(a,3(f0.5,1x))') ' h (1:3,1) = ', 1e8*boxndm%at(1:3,1)
    IF(RANG==0) WRITE(6,'(a,3(f0.5,1x))') ' h (1:3,2) = ', 1e8*boxndm%at(1:3,2)
    IF(RANG==0) WRITE(6,'(a,3(f0.5,1x))') ' h (1:3,3) = ', 1e8*boxndm%at(1:3,3)
    boxndm%wbox =wboxf*sum(0.5*cm(atpr%ityp(:atpr%im)))       ! La moitié de la masse totale des atomes
#ifdef PARA
    if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
       call comm_space%sum(boxndm%wbox)
    end if
#endif


    IF(RANG==0) WRITE(6,'(a,g20.12)')'Masse de la boîte pour Parrinello-Rahman: wbox=',boxndm%wbox

    ! État de référence défini par la matrice h0
    !   Cet état de référence doit correspondre à un tenseur de contrainte nul.
    !   Il n'est utile que pour calculer la déformation et l'énergie potentielle
    !   de la boîte.
    volu0 = calcvol(h0(1:3,1),h0(1:3,2),h0(1:3,3))
    invVolu0 = 1.d0/volu0
    trh0=Transpose(h0)
    CALL MatInv(h0,invh0)
    invtrh0=Transpose(invh0)

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
    IF(RANG==0) WRITE(6,'(a,f0.3,a)') 'Initialisation de la vitesse de la boîte pour la température ', TinitBox, ' K'
    boxndm%hdot(:,:) = 0.d0
    !#ifdef PARA
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
 !   boxndm%hdot=0
    call comm_space%bcast(0,boxndm%hdot)

    !#endif

    Kcell = 0.5d0*boxndm%wbox*Sum( boxndm%hDot(1:3,1:3)**2 )
    Tempcell=Kcell*2./(9.*bk)
!    if (Tinitbox.gt.0)     boxndm%hdot(1:3,1:3)=sqrt(tinitbox/Tempcell)*boxndm%hdot(1:3,1:3)
    DO i=1, 3
       DO j=1, 3

          ! Time derivative of the metric tensor Gmat
          boxndm%hdot(1:3,1:3)=boxndm%hdot(1:3,1:3)*ihbox0(1:3,1:3)
          boxndm%Gdot(i,j) = Sum( boxndm%h(1:3,i)*boxndm%hdot(1:3,j) + boxndm%hdot(1:3,i)*boxndm%h(1:3,j) )
       END DO
    END DO

    ! Kinetic energy of the cell (Eq. 2.14 of Ref. [2])

    Kcell = 0.5d0*boxndm%wbox*Sum( boxndm%hDot(1:3,1:3)**2 )


    if(lEev) then
       unitE=erg2eV
       cunitE='  eV'
    else
       unitE=1.0
       cunitE=' erg'
    end if
    if(rang==0) write(6,'(I7,D10.3,A,D21.12,A,a,f0.3,a)') 0,0.d0,'*Kcell = ',Kcell*unitE,cunitE, &
         '  (', 2.d0*Kcell/(9.d0*bk), ' K)'

    ! Nombre de thermostats de Hoover
    IF (nHoover.LT.0) nHoover=0
    IF (.NOT.lTHoover) nHoover=0
    if (.not.(allocated(zhoover)))ALLOCATE(zHoover(1:nHoover+1))


    IF (lTHoover) THEN

       ! Allocation des tableaux nécessaires au thermostat de Hoover
       IF (nHoover.LE.0) THEN
          WRITE(0,'(a)') 'Si le thermostat de Nosé-Hoover est utilisé (lTHoover=.true.), '
          WRITE(0,'(a)') 'le nombre de thermostats doit être supérieur à 0'
          WRITE(0,'(a,i0)') '  nHoover = ', nHoover
          STOP '< init_lpr >'
       END IF
       ALLOCATE(zOld(1:nHoover), zNew(1:nHoover), zDot(1:nHoover), wHoover(1:nHoover), &
            KHoover(1:nHoover))
       !ALLOCATE(UHoover(1:nHoover), UHoover_new(1:nHoover), UHoover_old(1:nHoover))

       ! Nombre de degrés de liberté pour le thermostat de Nosé-Hoover
       gNose=dble(3*atpr%im_glob)

       ! Masse de chaque thermostat
       IF (wNose.EQ.0) THEN
          ! On veut qu'une variation de la température de 10K corresponde à
          ! une variation de fNose de 1%
          wNose = gNose*bk*10.d0*tstep**2/1.d-2**2
       END IF
       wHoover(1)=wNose
       DO i=2, nHoover
          wHoover(i)=wNose/gNose
       END DO

       DO i=1, nHoover
          KHoover(i) = 0.5d0*bk*Text
          zHoover(i) = Sqrt( 2.d0*KHoover(i)/wHoover(i) )
          zDot(i) = 0.d0
          zOld(i) = zHoover(i) - tstep*zdot(i)
          !UHoover(i) = 0.d0
          !UHoover_old(i) = UHoover(i) - tstep*bk*Text*zHoover(i)
       END DO
       !UHoover(1) = gNose*UHoover(1)
       KNose = Sum( KHoover(1:nHoover) )
       !UNose = Sum( UHoover(1:nHoover) )
       zHoover(nHoover+1)=0.d0

       IF(RANG==0) WRITE(6,*)
       IF(RANG==0) WRITE(6,'(a)') 'Thermostat de Nosé-Hoover (V2)'
       IF(RANG==0) WRITE(6,'(a)') "  -> l'énergie cinétique des atomes et de la boîte est thermalisée"
       IF(RANG==0) WRITE(6,'(a,g20.12)')'Masse de la boîte pour thermostat de Nosé-Hoover: wHoover=',wNose
       IF(RANG==0) WRITE(6,'(a,g20.12)')'Nombre de degrés de liberté: gNose=',gNose
       IF(RANG==0) WRITE(6,'(a,i0)')    'Nombre de thermostats: nHoover=', nHoover
       IF(RANG==0) WRITE(6,'(a,f0.3,a)') 'Initialisation du thermostat de Nosé-Hoover pour la température ', &
            2.d0*KNose/(bk*dble(nHoover)), ' K'
       IF(RANG==0) WRITE(6,'(I7,D10.3,A,D21.12,A,a,f0.3,a)') iteration,timel,'*KNose = ',KNose*unitE,cunitE, &
            '  (', 2.d0*KNose/(bk*nHoover), ' K)'
       IF(RANG==0) WRITE(6,*)
    ELSE
       gNose=0.d0; zHoover(1)=0.d0
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
    real(double):: norme_de_fp, norme_de_vp, pscal,tempcell,pint,gs3
    real(double), dimension(ntyp) :: aux
    integer,save::nstep=0
    real(double)::maxvploc,maxvp,nrmvp
    integer::ib

    select case(dmtype)
    case(24)
       !       write(6,*)'IN',atpr%xp(1,1)
       ! Coordonnées réduites des atomes (au cas où elles ont été modifiées à l'extérieur)
       sp(:,1:atpr%im) = MatMul(boxndm%invh(:,:), atpr%xp(:,1:atpr%im) )
       select type (atpr)
       class is (atom_config_e)
          if (atpr%lxpp)       spp(:,1:atpr%im) = MatMul(boxndm%invh(:,:), atpr%xpp(:,1:atpr%im) )
       end select
       ! De même pour les vitesses au cas où, par exemple, on utilise le thermostat
       sdot(:,1:atpr%im) = MatMul(boxndm%invh(:,:), atpr%vp(:,1:atpr%im) )
       !       write(6,*)'VEL',sdot(1,1),atpr%vp(1,1)
       DO ia=1, atpr%im
          sfp(:,ia)= MatMul( boxndm%invh(:,:), atpr%fp(:,ia) )
       END DO
       im=atpr%im
       ! 1/ Intégration de l'équation de mouvement
       aux(:ntyp) = tstep**2/(2.d0*cm(:ntyp))
       usdh = 1.d0/(2.d0*tstep)
       DO i=1, im
          xprov(:) = sp(:,i) + sdot(:,i)*tstep + sfp(:,i)*aux(atpr%iTyp(i))
          !            if (i==1) write(6,*)'VPROV',(xprov(1) - spp(1,i))*usdh
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
       !         write(6,*)'PSCAL',pscal,hdot(1,1)
       ! Modification du vecteur vitesse
       !         write(6,*)'PSCAL',iteration,nstep,pscal,tstep
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
       forcebox(:,:)=MatMul( sigtot(:,:) - sigext(:,:), boxndm%invtrh(:,:) )*boxndm%volu
       do i = 1, 3
          do ic = 1, 3
             if (boxndm%hdot(ic,i)*forcebox(ic,i)<0) then
                boxndm%hdot(ic,i)=0.
             end if
          end do
       end do

       !          write(6,*)'hdot',hdot(1,1),tstep/(1*wBox)*forcebox(1,1)*ihbox0(1,1),invtrh(1,1),sigtot(1,1),tstep,tstepN
       boxndm%hdot(:,:) = boxndm%hdot(:,:)*ihbox0(:,:)+ tstep/(1*boxndm%wBox)*forcebox(:,:)*ihbox0(:,:)

       ! Tenseur h à l'instant t+dt
       boxndm%h(:,:) = boxndm%h(:,:) + boxndm%hdot(:,:)*tstep*ihbox0(:,:)
       !          write(6,*)'FHdot',hdot(1,1),tstep,(1*wBox),forcebox(1,1)*ihbox0(1,1),h(1,1)

       !       write(6,*)'MED',atpr%xp(1,1),sp(1,1),sfp(1,1)
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
       !    write(6,*)'MED2',atpr%xp(1,1),sp(1,1),sfp(1,1)
       ! Calcul des forces et des contraintes à l'instant t+dt
       CALL CalFo(sig,potist,atpr,celndm,boxndm%box_config,t_sigma=.true.,psc=psc)
       if (dmtype==24) tstep=tstepN


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

       sigtot = 0.5d0*(sigkine + Transpose(sigkine) + sig + Transpose(sig) )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! LPR +LANGEVIN : evolution de VP en corrdonnées réelles pour éviter de se tromper dans les dimensions       
    case(88)
       call comm_space%barrier
       select type(atpr)
       class is (atom_config_e)
          Kcell = 0.5d0*boxndm%wbox*Sum( boxndm%hDot(1:3,1:3)**2 )
          Tempcell=Kcell*2./(sum(ihbox0)*bk)
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
                gs3=(glanh(1,1)+glanh(2,2)+glanh(3,3))/3.
                glanh=0.
                do ic=1,3
                   glanh(ic,ic)=gs3
                end do
             end if
                
             boxndm%hdot(:,:) = (  boxndm%hdot(:,:)*rgah  &
                  + tstep/(2.d0*boxndm%wBox)*boxndm%volu*MatMul( sigtot(:,:) - sigext(:,:),boxndm%invtrh(:,:) ) &
                  + (glanh(:,:)/boxndm%wbox)*sqrt(boxndm%wbox*bk*text*(1-rgah))  )*ihbox0(:,:)
          end if
          call comm_space%bcast(0,boxndm%hdot)
          call comm_space%bcast(0,glanh)
         tempx= tempinstT(atpr)
          Kcell = 0.5d0*boxndm%wbox*Sum( boxndm%hDot(1:3,1:3)**2 )
          Tempcell=Kcell*2./(sum(ihbox0)*bk)

          sp(:,1:atpr%im) = sp(:,1:atpr%im) + sdot(:,1:atpr%im)*tstep

          ! Tenseur h à l'instant t+dt
          boxndm%h(:,:) = boxndm%h(:,:) + boxndm%hdot(:,:)*tstep*ihbox0(:,:)

          ! Coordonnées réelles à l'instant t+dt
          select type (atpr)
          class is (atom_config_e)
             if (atpr%lxpp)          atpr%xpp(:,1:atpr%im) = atpr%xp(:,1:atpr%im)
          end select
          atpr%xp(:,1:atpr%im) = MatMul( boxndm%h, sp(:,1:atpr%im) )

          call updatebox(boxndm,boxndm%h)
          atpr%vp(:,1:atpr%im) = MatMul( boxndm%h(:,:), sdot(:,1:atpr%im) ) ! retour à vp car transfert d'atomes  dans scalebox en PARA
          Kcell = 0.5d0*boxndm%wbox*Sum( boxndm%hDot(1:3,1:3)**2 )
          Tempcell=Kcell*2./(sum(ihbox0)*bk)

          tempx= tempinstT(atpr)

          CALL ScaleBox(atpr,celndm,boxndm,psc)
          call caltabtC(celndm,atpr,lperiod,boxndm)

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

          boxndm%hdot(:,:) = (boxndm%hdot(:,:)*rgah  &
               + tstep/(2.d0*boxndm%wBox)*boxndm%volu*MatMul( sigtot(:,:) - sigext(:,:), boxndm%invtrh(:,:) ) &
               + (glanh(:,:)/boxndm%wbox)*sqrt(boxndm%wbox*bk*text*(1-rgah))  )*ihbox0(:,:)
          ! Estimation de la dérivée du tenseur Gmat à l'instant t+dt
          DO i=1, 3
             DO j=1, 3
                boxndm%Gdot(i,j) = Sum( boxndm%h(1:3,i)*boxndm%hdot(1:3,j) + boxndm%hdot(1:3,i)*boxndm%h(1:3,j) )
             END DO
          END DO

          grsig = boxndm%volu * MatMul(boxndm%invh, MatMul( sigext, boxndm%invtrh) )
          tension = invVolu0*MatMul( MatMul( h0, grsig), trh0 )

          ! Énergie potentielle de la cellule (Eq. 2.25, Ref.2)
          Ucell = volu0*Sum( tension(1:3,1:3) * epsi(1:3,1:3) )

          ! Énergie cinétique de la cellule (Eq. 2.14, Ref.2)
          Kcell = 0.5d0*boxndm%wbox*Sum( boxndm%hDot(1:3,1:3)**2 )
          Tempcell=Kcell*2./(sum(ihbox0)*bk)
          EcellPR = Kcell + Ucell

       end select
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    case(22,8)
       if (dmtype==22) then

          do i = 1, atpr%im
             do ic = 1, 3
                if (atpr%vp(ic,i)*atpr%fp(ic,i)<0) then
                   atpr%vp(ic,i)=0.
                end if
             end do
          end do
          forcebox(:,:)=MatMul( sigtot(:,:) - sigext(:,:), boxndm%invtrh(:,:) )
          do i = 1, 3
             do ic = 1, 3
                if (boxndm%hdot(ic,i)*forcebox(ic,i)<0) then
                   boxndm%hdot(ic,i)=0.
                end if
             end do
          end do
       end if

       !    case(8)
       ! Coordonnées réduites des atomes (au cas où elles ont été modifiées à l'extérieur)
       sp(:,1:atpr%im) = MatMul(boxndm%invh(:,:), atpr%xp(:,1:atpr%im) )
       ! De même pour les vitesses au cas où, par exemple, on utilise le thermostat
       sdot(:,1:atpr%im) = MatMul(boxndm%invh(:,:), atpr%vp(:,1:atpr%im) )
       ! Dérivée des coordonnées réduites des atomes à l'instant t+dt/2
       mf(:,:) = -0.5d0*tstep*MatMul(boxndm%invGmat,boxndm%Gdot)
       DO i=1, 3
          mf(i,i) = 1.d0 - 0.5d0*tstep*zHoover(1) + mf(i,i)
       END DO
       DO ia=1, atpr%im
          sdot(:,ia) = MatMul( mf(:,:), sdot(:,ia) ) &
               + tstep/(2.d0*cm(atpr%ityp(ia))) * MatMul( boxndm%invh(:,:), atpr%fp(:,ia) )
       END DO
       ! Dérivée du tenseur h à l'instant t+dt/2
       IF (lpcon2.EQV..true.) THEN    ! On ajoute une force de friction
          boxndm%hdot(:,:) = ( 1.d0 - 0.5d0*tstep*zHoover(1) - 0.5d0/tbox )* boxndm%hdot(:,:)*ihbox0(:,:) &
               + tstep/(2.d0*boxndm%wBox)*boxndm%volu*MatMul( &
               &sigtot(:,:) - sigext(:,:), boxndm%invtrh(:,:) )*ihbox0(:,:)
       ELSE        ! Équation sans force de friction supplémentaire
          boxndm%hdot(:,:) = ( 1.d0 - 0.5d0*tstep*zHoover(1) )* boxndm%hdot(:,:)*ihbox0(:,:) &
               + tstep/(2.d0*boxndm%wBox)*boxndm%volu*MatMul(&
               &sigtot(:,:) - sigext(:,:), boxndm%invtrh(:,:) )*ihbox0(:,:)
       END IF
       !          write(6,*)'hdot',hdot
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
       !    write(6,*)'MED2',atpr%xp(1,1),sp(1,1),sfp(1,1)
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
!!$             if (rang==0) write(6,*)'MAXVP', maxvp
!           end block
    
       CALL CalFo(sig,potist,atpr,celndm,boxndm%box_config,t_sigma=.true.,psc=psc)
       if (l2t)then
          if (i2t==1)  call calceloss(celndm,atpr)
       else
          if(ibrake.gt.0) call calceloss(celndm,atpr)
       end if
       if (lTberendsen) call calfoberend(atpr)

       ! Calcul de la viscosité à l'instant ...
       DO i=1, nHoover
          zNew(i) = zOld(i) + 2.d0*zDot(i)*tstep    ! ... t+dt
          zOld(i) = zHoover(i)                      ! ... t
          zHoover(i) = zNew(i)                      ! ... t+dt
       END DO

       ! Estimation de la contrainte totale à l'instant t+dt
       ! (la contrainte cinétique est calculée à l'instant t
       ! et la contrainte potentielle à l'instant t+dt)
       sigtot = 0.5d0*(sigkine + Transpose(sigkine) + sig + Transpose(sig) )

       ! Estimation de la dérivée du tenseur h à l'instant t+dt
       IF (lpcon2.EQV..true.) THEN    ! On ajoute une force de friction
          hdot_new(:,:) = 1.d0/(1.d0+0.5d0*tstep*zHoover(1)+ 0.5d0/tbox )*( boxndm%hdot(:,:)*ihbox0(:,:) &
               + tstep/(2.d0*boxndm%wBox)*boxndm%volu*MatMul(sigtot(:,:)-sigext(:,:),boxndm%invtrh(:,:)))*ihbox0(:,:)
       ELSE        ! Équation sans force de friction supplémentaire
          hdot_new(:,:) = 1.d0/( 1.d0 + 0.5d0*tstep*zHoover(1) )*( boxndm%hdot(:,:)*ihbox0(:,:) &
               + tstep/(2.d0*boxndm%wBox)*boxndm%volu*MatMul(sigtot(:,:)-sigext(:,:),boxndm%invtrh(:,:)))*ihbox0(:,:)
       END IF
       ! Estimation de la dérivée du tenseur Gmat à l'instant t+dt
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
             mf(i,i) = 1.d0 + 0.5d0*tstep*zHoover(1) + mf(i,i)
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
          ! Dérivée du tenseur h à l'instant t+dt
          IF (lpcon2.EQV..true.) THEN    ! On ajoute une force de friction
             hdot_new(:,:) = 1.d0/( 1.d0 + 0.5d0*tstep*zHoover(1) + 0.5d0/tbox )*(boxndm%hdot(:,:)*ihbox0(:,:) &
                  + tstep/(2.d0*boxndm%wBox)*boxndm%volu*MatMul( sigtot(:,:)-sigext(:,:),boxndm%invtrh(:,:)))*ihbox0(:,:)
          ELSE        ! Équation sans force de friction supplémentaire
             hdot_new(:,:) = 1.d0/( 1.d0 + 0.5d0*tstep*zHoover(1) )*(boxndm%hdot(:,:) &
                  + tstep/(2.d0*boxndm%wBox)*boxndm%volu*MatMul(sigtot(:,:)-sigext(:,:),boxndm%invtrh(:,:)))*ihbox0(:,:)
          END IF

          ! Dérivée du tenseur Gmat à l'instant t+dt
          DO i=1, 3
             DO j=1, 3
                boxndm%Gdot(i,j) = Sum( boxndm%h(1:3,i)*hdot_new(1:3,j) + hdot_new(1:3,i)*boxndm%h(1:3,j) )
             END DO
          END DO

          ! Vérifie l'autocohérence de h
          diff = Sum( abs( hdot_new(1:3,1:3) - hdot_last(1:3,1:3) ) )
          tdiff = Sum( abs( hdot_last(1:3,1:3) ) )
          if (tdiff .eq. 0.0d0) then
             if (diff .LE. tol) exit
          else
             if (diff/tdiff .LE. tol) exit
          endif

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
       epsi=0.5d0*MatMul( MatMul( invtrh0, boxndm%Gmat ), invh0 )
       DO i=1, 3
          epsi(i,i) = epsi(i,i) - 1.d0
       END DO
       ! Tension thermodynamique (Eq. 2.22 et 2.26, Ref.2)
       grsig = boxndm%volu * MatMul(boxndm%invh, MatMul( sigext, boxndm%invtrh) )
       tension = invVolu0*MatMul( MatMul( h0, grsig), trh0 )

       ! Énergie potentielle de la cellule (Eq. 2.25, Ref.2)
       Ucell = volu0*Sum( tension(1:3,1:3) * epsi(1:3,1:3) )

       ! Énergie cinétique de la cellule (Eq. 2.14, Ref.2)
       Kcell = 0.5d0*boxndm%wbox*Sum(boxndm%hDot(1:3,1:3)**2 )
       EcellPR = Kcell + Ucell
       IF (lTHoover) THEN
          ! Dérivée de la viscosité et énergie cinétique du thermostat
          zDot(1) = (2.d0*(kine + Kcell) - gNose*bk*Text)/wHoover(1) &
               - zHoover(2)*zHoover(1)
          KHoover(1) = 0.5d0*wHoover(1)*zHoover(1)**2         ! t+dt
          DO i=2, nHoover
             zDot(i) = ( 2.d0*KHoover(i-1) - bk*Text )/wHoover(i) &
                  - zHoover(i+1)*zHoover(i)
             KHoover(i) = 0.5d0*wHoover(i)*zHoover(i)**2
          END DO
          KNose = Sum( KHoover(1:nHoover) )
          ! Énergie potentielle du thermostat
          !DO i=1, nHoover
          !   UHoover_new(i) = UHoover_old(i) + 2.d0*tstep*bk*Text*zHoover(i)   ! t+dt
          !   UHoover_old(i) = UHoover(i)              ! t
          !   UHoover(i) = UHoover_new(i)              ! t+dt
          !END DO
          !UHoover(1) = gNose*UHoover(1)               ! t+dt
          !UNose = Sum( UHoover(1:nHoover) )           ! t+dt
          !ENose = KNose + UNose                       ! t+dt
       END IF
    end select
    !    write(6,*)'OUT',atpr%xp(1,1),tstep

    call calctemp(T1,kin1,atpr,celndm)

  end subroutine pr1

end module Parrinello_Rahman !Parrinello_Rahman





!******************
