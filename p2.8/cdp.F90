module cdp_mod
  USE T_kind_param_m, ONLY:  double
  USE gen_com_m, ONLY: iseed_glob=>iseed,rang,dmtype,itmax,lspacendm,it,lperiod,itloopmax,ivisu
  !  use temp_com,only:im
  USE arret_ndm_mod,only: arret_ndm
  USE var_pot,only:ntyp
  USE atomconfig,only : atom_config,atom_config_d
  USE boxconfig,only:box_config
  USE cellconfig, only:cell_config,caltabtC
  use rasmolT_mod,only:rasmolT
  use vect_dist_mod,only:closest_at
  USE cryst_to_cart_mod,only: cryst_to_cart

  USE Tpara,only:COMM_space,myidsp,para_space_config
#ifdef PARA  
  USE mod_para,only:maj_atomes_frt_ftm
#endif
  use constrconf_mod,only:coord_to_cell
!#else
!  USE Tpara,only:myidsp,para_space_config
!#endif
  !  use endrunT_mod,only:endrunT
  use NGC_mod,only:ngc
  use dmloop_pilot_mod,only:dmloop_pilot

  implicit none

  integer :: &
       itecdp, &     ! introduction de DP tout les itecdp pas
       nposI,&        ! nombre de positions interstitielles
       iseed, &      ! racine des nombres aléatoires
       ideftyp, &      ! racine des nombres aléatoires
       typint ! type d'introduction des Intestitiels : 0 dans les sites prédéfinis, 1 aléatoirement
  integer,allocatable::nvac(:),nbint(:)
  real(double), dimension(:,:), allocatable :: xposint ! positions des interstitiels POSSIBLES
  real(double), dimension(:,:), allocatable :: xposI ! positions des interstitiels réalisés
  real(double)::maxposint(3),minposint(3)
  real (double) :: dminins,rsphdef,centresphdef(3)
  integer:: ioxdef,itprep


contains
  ! **************************************************************
  subroutine initcdp
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    implicit none

    !-----------------------------------------------
    !   L o c a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: i,itapp,nfp,iti
    !-----------------------------------------------
    namelist /inputcdp/itecdp,nfp,nposI,iseed,dminins,itecdp,itprep,maxposint,minposint,nvac,nbint,typint

    allocate(nvac(ntyp));allocate(nbint(ntyp))

    itprep=1
    itecdp=-1      ! introduction de DP tout les itecdp pas
    nvac(:)=0 ! number of vacancies 
    nbint(:)=0 ! number of interstitials 
    nposI=-1      ! nombre de positions interstitielles
    iseed=-1      ! graine pour la generation aleatoire si <0 tirage avec SECNDS
    dminins=1.0   ! distance minimum entre nouvel interstitiel et atomes deja present
    typint=1     ! type d'introduction des Intestitiels : 0 dans les sites prédéfinis, 1 aléatoirement
    minposint(1:3) =0;maxposint(1:3)=1
    nfp=-1

    open(unit=73, file='creaDPin', status='unknown')
    read (73, nml=inputcdp)
    if (nfp.gt.0) then
       if (rang==0) write(6,*)'NFP VAC INT for all types'
       nvac=nfp
       nbint=nvac
    end if
    dminins=dminins*1d-8
    if(all(nvac==-1).and.all(nbint==-1)) then
       if (rang==0)write(6,*)'what defects ?'
       call arret_ndm
    end if
    if(itecdp==-1) then
       if (rang==0)write(6,*)'when defects ?'
       call arret_ndm
    end if
    if((typint.lt.0).or.(typint.GT.1)) then
       write(6,*)'mauvaise introduction des interstitiels stop'
       call arret_ndm
    end if


    if (typint==0) then
       allocate(xposint(3,nposI))
       do i=1,nposI
          read(73,*)xposint(1,i),xposint(2,i),xposint(3,i)
       end do
       where (xposint(:,:)<0.0)
          xposint(:,:)=xposint(:,:)+1.
       end where
       where (xposint(:,:)>1.0)
          xposint(:,:)=xposint(:,:)-1.
       end where
    end if
    ! call    random_seed (put=iseedt)

    return
  end subroutine initcdp
  !**********************************************************
  subroutine creadp(atdml,celndm,boxndm,psc)
    USE var_pot, ONLY:ntyp,ty
    implicit none
    type atomvac_typ
       integer,allocatable,dimension(:)::iproc,natg,iloc,ityp
       real(double),allocatable::pos(:,:)
    end type atomvac_typ
    type atomint_typ
       integer,allocatable,dimension(:)::natg,iatpos,ityp
       real(double),allocatable::pos(:,:)
    end type atomint_typ

    integer, dimension(2) :: iseedt
    type(para_space_config)::psc
    type(box_config)::boxndm
    class(atom_config)::atdml
    type(cell_config):: celndm

    type(atomvac_typ)::atomvac
    type(atomint_typ)::atomint
    integer :: ic,j,ntry,iti,i,natyp,nvactot,iat,ivac,ivacloc,ivactot,jvac,ninttot,iproc
    integer :: itapp,npp
    integer :: idep,itinser
    integer :: iposI,iint,natgm
    real(double) :: a1,a2,a3,c1,c2,c3,z1,r2,rd,z3,z2, edt
    real(double),dimension(3):: xdec, xavant,xapres,xpositest
    real(double), dimension(1,3) :: cv
    real(double),dimension(:),allocatable:: edrat
    integer,allocatable::nb_at_typ(:),last_at_typ(:),iatvac(:)
    logical::l2close
    integer::numproc,iatint
    integer::jint,iinttot,numcell,imt,ntry2
    if (myidsp==0) then
       if (iseed.le.0) then
          call system_clock (iseed) 
       end if
       write(6,*)'rang iseed CREADP',rang,iseed
       iseedt(1)=iseed
       call random_seed (iseedt(1))
    end if


    !    call atdml%print
    if (itprep.gT.0) then 
       itloopmax=itprep
       itinser=0
       select type(atdml)
       type is (atom_config)
          select case (dmtype)
          case(32,33,34)
             call NGC (atdml,celndm,boxndm,psc)
          case default
             stop
          end select
       class is (atom_config_d)

          select case (dmtype)
          case(32,33,34)
             call NGC (atdml,celndm,boxndm,psc)
          case(4,10,8,1,21,22)
             call dmloop_pilot(atdml,celndm,boxndm,psc)
          case default
             stop
          end select
       end select
    end if

    nvactot=sum(nvac(1:ntyp)); ninttot=sum(nbint(1:ntyp))
    allocate(atomvac%iproc(nvactot))
    allocate(atomvac%ityp(nvactot))
    allocate(atomvac%natg(nvactot))
    allocate(atomvac%iloc(nvactot))
    allocate(atomvac%pos(3,nvactot))
    atomvac%natg=0
    atomvac%ityp=0
    atomvac%iproc=0
    atomvac%iloc=0
    allocate(atomint%iatpos(ninttot))
    allocate(atomint%ityp(ninttot))
    allocate(atomint%pos(3,ninttot))

!#ifdef PARA
    npp=comm_space%nproc
!#else
!    npp=1
!#endif

    allocate(nb_at_typ(0:npp-1))
    allocate(last_at_typ(-1:npp-1))
    open (unit=121,file='vac_int')
    do while (it.le.itmax)
       last_at_typ=0
       natyp=0
       nb_at_typ=0
       itinser=itinser+1
       call rasmolT(atdml,boxndm,itinser,'PRE_INSER',latcomp=.false.,ivisumol=ivisu)

       itloopmax=min(it+itecdp,itmax)
       if (myidsp==0) then
          write(6,*)'****************************************'
          write(6,*)'POINT DEFECT CREATION ',nvactot, ninttot
          write(6,*)'****************************************'
       end if
       natgm=maxval(atdml%num_at_glob(1:atdml%im))
!#ifdef PARA
       call comm_space%max(natgm)
!#endif


       ivactot=0
       atomvac%natg=0 ! on remet à 0 les indices
       atomvac%iproc=0
       atomvac%iloc=0
       atomvac%ityp=0
       atomint%iatpos=0
       atomint%ityp=0
       atomint%pos=0
       if (nvactot.ne.0) then
          ! insérer les lacunes
          do iti=1,ntyp
             allocate(iatvac(nvac(iti)))
             iatvac=0
             last_at_typ(:)=0
             natyp=count(atdml%ityp(1:atdml%im)==iti)
!#ifdef PARA
             nb_at_typ(comm_space%rank)=natyp
             if (lspacendm) then
!                call comm_space%build(natyp,nb_at_typ,torank=0)
                call comm_space%sum(nb_at_typ)
                call comm_space%sum(natyp)
             end if
             do iproc=0,comm_space%nproc-1
                last_at_typ(iproc)=last_at_typ(iproc-1)+nb_at_typ(iproc)
             end do
!#else
!             last_at_typ(0)=natyp
!#endif
             if (natyp.lt.nvac(iti)) then
                write(6,*)'impossible to delete that many atoms of type ',iti,nvac(iti),natyp
                call arret_ndm
             end if
             do ivac=1,nvac(iti)
                ivactot=ivactot+1 ! indice l'ensemble des lacunes (inter-types)

                if (myidsp==0) then
                   ntry=0
1                  continue
                   ntry=ntry+1
                   if (ntry==100) then
                      write(6,*)'VAC NTRY exceeded'
                      call arret_ndm
                   end if
                   call random_number(z1)
                   iatvac(ivac)=1+int(z1*natyp)
                   do jvac=1,ivac-1
                      if (iatvac(jvac)==iatvac(ivac)) goto 1
                   end do
!#ifdef PARA
                   if (lspacendm) then
                      looppr:do iproc=0,comm_space%nproc-1
                         if (last_at_typ(iproc).ge.iatvac(ivac)) then
!                         write(6,*)'CHOIX',iproc, last_at_typ(iproc),iatvac(ivac)
                            ivacloc=iatvac(ivac)-last_at_typ(iproc-1)
                            exit looppr !iproc est l'indice du proc qui contient iatvac et ivacloc est le rang  de l'atome de cette lacune dans les atomes de ce type
                         end if
                      end do looppr
                   else
                      iproc=0
                      ivacloc=iatvac(ivac)
                   endif
!#else
!                   iproc=0
!                   ivacloc=iatvac(ivac)
!#endif
                end if

!#ifdef PARA
                call comm_space%bcast(0,iproc)
                call comm_space%bcast(0,ivacloc)
                if (myidsp==iproc) then
!#endif
                   iat=0
                   loopi: do i=1,atdml%im
                      if (atdml%ityp(i)==iti) then
                         iat=iat+1
                         if (iat==ivacloc) then
                            atomvac%iproc(ivactot)=myidsp
                            atomvac%pos(:,ivactot)=atdml%xp(:,i)
                            atomvac%ityp(ivactot)= atdml%ityp(i)
                            atomvac%natg(ivactot)=atdml%num_at_glob(i)
                            atomvac%iloc(ivactot)=i  ! %iloc est le numéro de l'atome de la alcune (tous types confondus) <> ivacloc
                            exit loopi
                         end if
                      end if
                   end do loopi
!#ifdef PARA
                end if
!#endif

             end do
!#ifdef PARA
             call comm_space%sum(atomvac%iproc) ! avant ça seul le proc iproc connaissait ces chiffres
             call comm_space%sum(atomvac%natg)
             call comm_space%sum(atomvac%iloc)
             call comm_space%sum(atomvac%pos)
             call comm_space%sum(atomvac%ityp)
!#endif
             deallocate (iatvac)
          end do

          do ivactot=1,nvactot
!#ifdef PARA          
             if (comm_space%rank==atomvac%iproc(ivactot)) then
!#endif
                call atdml%switch_atom(atomvac%iloc(ivactot),atdml%im)
                atdml%im=atdml%im-1
!#ifdef PARA          
             end if
!#endif
          end do
          if (myidsp==0) then
             write(121,*)itinser, 'VAC'
             do i=1,nvactot
                write(121,'(3G15.6,2I6)')atomvac%pos(:,i),atomvac%ityp(i),atomvac%natg(i)
             end do
             flush(121)
          end if
          
       end if
       
       ! insérer les interstitiels       
       if (ninttot.ne.0) then
          iinttot=0 
          do iti=1,ntyp
             do iint=1,nbint(iti)
                l2close=.true. ! le do while doit être fait au mins une fois
                ntry=0; ntry2=0
                iinttot=iinttot+1 ! indice l'ensemble des interstitiels (inter-types)
                do while (l2close)
                   if (myidsp==0) then
                      select case(typint)
                      case(0)
2                        continue
                         ntry=ntry+1
                         if (ntry==100) then
                            write(6,*)' INT typ0 NTRY exceeded'
                            call arret_ndm
                         end if
                         call random_number(z1)
                         iatint=1+int(z1*nposI)
                         do jint=1,iinttot-1
                            if(atomint%iatpos(jint)==iatint) goto 2
                         end do
                         atomint%iatpos(iinttot)=iatint
                         xpositest(:)=xposint(:,iatint)
                      case(1)
                         z1=-1.0
                         do while ((z1.lt.minposint(1)).or.z1.gt.maxposint(1))
                            call random_number(z1)
                         end do
                         xposItest(1)=z1

                         z1=-1.0
                         do while ((z1.lt.minposint(2)).or.z1.gt.maxposint(2))
                            call random_number(z1)
                         end do
                         xpositest(2)=z1

                         z1=-1.0
                         do while ((z1.lt.minposint(3)).or.z1.gt.maxposint(3))
                            call random_number(z1)
                         end do
                         xposItest(3)=z1
                      end select
                      call cryst_to_cart (1, xpositest, boxndm%at, 1) !cryst vers cart
                   end if

#ifdef PARA
                   call comm_space%bcast(0,xpositest)
                   if (lspacendm) then
                      call coord_to_cell(xposItest,numcell,boxndm%bg,celndm%nox,celndm%noy,celndm%noz)
                      numproc=celndm%proc_cell(numcell)
                   else
                      numproc=0
                   end if
!                   if (myidsp==0)write(6,*)'PROCint',numcell,numproc,xpositest
#endif
                   if (numproc == myidsp) then
                      ntry2=ntry2+1
                      if (ntry2==100) then
                         call arret_ndm
                      end if

                      l2close=.false.             
                      if (dminins.gT.0) then
                         call closest_at(xpositest,atdml,celndm,boxndm,lperiod,rumin=dminins,lclose=l2close)
                      end if
!#ifdef PARA
                   end if
                   call comm_space%bcast(numproc,l2close)
!#endif

                   if (.not.l2close) then ! not too close ==> intertsitiel+1
                      atomint%ityp(iinttot)=iti
                      atomint%pos(:,iinttot)=xpositest(:)
                      natgM=maxval(atdml%num_at_glob(1:atdml%im))
!#ifdef PARA
                      if (lspacendm) call comm_space%max(natgM)
                      if (myidsp==numproc) then
!#endif                   
                         atdml%im=atdml%im+1
                         atdml%xp(:,atdml%im)=xpositest(:)
                         atdml%ityp(atdml%im)=iti
                         atdml%fp(:,atdml%im)=0
                         atdml%num_at_glob(atdml%im)=natgM+1
                      

                         select type(atdml)
                         class is (atom_config_d)
                            atdml%vp(:,atdml%im)=0
                            atdml%xpp(:,atdml%im)=atdml%xp(:,atdml%im)
                         end select
!#ifdef PARA
                      end if
!#endif
                   end if
                end do
             end do
          end do

          if (iinttot.ne.ninttot) then
             write(6,*)'pb nombre de int',iinttot,ninttot
             call arret_ndm
          end if
          if (myidsp==0) then
             write(121,*)itinser, 'INT'
             do i=1,ninttot
                write(121,'(3G15.6,I6)')atomint%pos(:,i),atomint%ityp(i)
             end do
          end if
          
       end if
       atdml%im_glob=atdml%im_glob-nvactot+ninttot
!#ifdef PARA

       imt=atdml%im
       if (lspacendm)call comm_space%sum(imt)
       if(imt.ne.atdml%im_glob) then
          write(6,*)'imt <> %im_glob'
          call arret_ndm
       end if
!#endif
       call caltabtC(celndm,atdml,lperiod,boxndm)
#ifdef PARA
       if (lspacendm) call maj_atomes_frt_ftm(atdml,celndm,boxndm,psc)
#endif

       if (lspacendm) then
          call rasmolT(atdml,boxndm,itinser,'POST_INSER',latcomp=.false.,ivisumol=ivisu)
       else
          call rasmolT(atdml,boxndm,itinser,'POST_INSER',latcomp=.true.,ivisumol=ivisu)
       end if
       select type(atdml)
       type is (atom_config)
          select case (dmtype)
          case(32,33,34)
             call NGC (atdml,celndm,boxndm,psc)
          case default
             write(6,*)'WTFCDP1'
             call arret_ndm
          end select
       class is (atom_config_d)

          select case (dmtype)
          case(32,33,34)
             call NGC (atdml,celndm,boxndm,psc)
          case(4,10,8,1,21,22)
             call dmloop_pilot(atdml,celndm,boxndm,psc)
          case default
             write(6,*)'WTFCDP2'
             call arret_ndm
          end select
       end select

    end do

  end subroutine creadp

end module cdp_mod


