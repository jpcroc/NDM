module cdp_mod
  USE T_kind_param_m, ONLY:  double
  USE gen_com_m, ONLY: iseed_glob=>iseed,rang,dmtype,itmax,lspacendm,it,lperiod
  !  use temp_com,only:im
  USE arret_ndm_mod,only: arret_ndm
  USE var_pot,only:ntyp
  USE atomconfig,only : atom_config,atom_config_d
  USE boxconfig,only:box_config
  USE cellconfig, only:cell_config,caltabtC

#ifdef PARA  
  USE Tpara,only:COMM_space,myidsp,para_space_config
  USE mod_para,only:maj_atomes_frt_ftm

#else
  USE Tpara,only:myidsp,para_space_config
#endif
!  use endrunT_mod,only:endrunT
  use NGC_mod,only:ngc
  use dmloop_pilot_mod,only:dmloop_pilot

  implicit none

  integer :: &
       itecdp, &     ! introduction de DP tout les itecdp pas
       nposI,&        ! nombre de positions interstitielles
       iseed, &      ! racine des nombres aléatoires
       ideftyp, &      ! racine des nombres aléatoires
       typint,&! type d'introduction des Intestitiels : 0 dans les sites prédéfinis, 1 aléatoirement
       ittot ! nombre total d'itérations (va remplacer itmax)
  integer,allocatable::nvac(:),nbint(:)

  real(double), dimension(:,:), allocatable :: xposint ! positions des interstitiels
  real (double) :: dminins,rsphdef,centresphdef(3)
  integer:: ioxdef


contains
  ! **************************************************************
  subroutine initcdp
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    implicit none
    integer, dimension(2) :: iseedt
    !-----------------------------------------------
    !   L o c a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: i,itapp,nfp
    !-----------------------------------------------
    namelist /inputcdp/itecdp,nfp,nposI,iseed,dminins,itecdp

    allocate(nvac(ntyp));allocate(nbint(ntyp))

    itecdp=-1      ! introduction de DP tout les itecdp pas
    nvac(:)=-1 ! number of vacancies 
    nbint(:)=-1 ! number of interstitials 
    nposI=-1      ! nombre de positions interstitielles
    iseed=-1      ! graine pour la generation aleatoire si <0 tirage avec SECNDS
    dminins=1.3   ! distance minimum entre nouvel interstitiel et atomes deja present
    typint=0     ! type d'introduction des Intestitiels : 0 dans les sites prédéfinis, 1 aléatoirement

    open(unit=73, file='creaDPin', status='unknown')
    read (73, nml=inputcdp)

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

    ittot=itmax
    
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

       !         do i=1,nposI
       !            xposint(:,i)=(xposint(:,i)-0.5)*zl(:)
       !         end do
    end if
    if (iseed.le.0) then
       call system_clock (iseed) 
       write(6,*)'rang iseed ',rang,iseed
    end if
    iseedt(1)=iseed
    call random_seed (iseedt(1))
    ! call    random_seed (put=iseedt)

    return
  end subroutine initcdp

  subroutine creadp(atdml,celndm,boxndm,psc)
    USE var_pot, ONLY:ntyp,ty
    implicit none
    type atomvac_typ
       integer,allocatable,dimension(:)::iproc,natg,iloc
    end type atomvac_typ

    type(para_space_config)::psc
    type(box_config)::boxndm
    class(atom_config)::atdml
    type(cell_config):: celndm

    type(atomvac_typ),allocatable::atomvac
    integer :: ic,j,ntry,iti,i,natyp,nvactot,iat,iproc,ivac,ivacloc,ivactot,jvac
    integer :: itapp,npp
    integer :: idep
    integer :: iposI,iint,natgm
    real(double) :: a1,a2,a3,c1,c2,c3,z1,r2,rd,z3,z2, edt
    real(double),dimension(3):: xdec, xavant,xapres,xposinttest
    real(double), dimension(1,3) :: cv
    real(double),dimension(:),allocatable:: edrat
    integer,allocatable::nb_at_typ(:),last_at_typ(:),iatvac(:)
    
    nvactot=sum(nvac(1:ntyp))
    allocate(atomvac%iproc(nvactot))
    allocate(atomvac%natg(nvactot))
    allocate(atomvac%iloc(nvactot))
    atomvac%natg=0
    atomvac%iproc=0
    atomvac%iloc=0
#ifdef PARA
    npp=comm_space%nproc
#else
    npp=1
#endif
    
    allocate(nb_at_typ(0:npp-1))
    allocate(last_at_typ(-1:npp-1))

    do while (it.le.ittot)
    
       
       itmax=it+itecdp

       select type(atdml)
       type is (atom_config)
          select case (dmtype)
          case(32,33,34)
             call NGC (atdml,celndm,boxndm,psc)
          case default
             write(6,*)'WTFCDP'
             stop
          end select
       class is (atom_config_d)
          
          select case (dmtype)
          case(32,33,34)
             call NGC (atdml,celndm,boxndm,psc)
          case(4,10,8,1,21,22)
             call dmloop_pilot(atdml,celndm,boxndm,psc)
          case default
             write(6,*)'WTFCDP'
             stop
          end select
       end select
       
       natgm=maxval(atdml%num_at_glob(1:atdml%im))
#ifdef PARA
       call comm_space%max(natgm)
#endif
    
       
       ivactot=0
       ! insérer les défauts
       do iti=1,ntyp
          allocate(iatvac(nvac(iti)))
          last_at_typ(:)=0
          natyp=count(atdml%ityp(1:atdml%im)==iti)
#ifdef PARA
          if (lspacendm) then
             call comm_space%build(natyp,nb_at_typ,torank=0)
             call comm_space%sum(natyp,torank=0)
          end if

          do iproc=0,comm_space%nproc-1
             last_at_typ(iproc)=last_at_typ(iproc-1)+nb_at_typ(iproc)
          end do
#endif
          if (natyp.lt.nvac(iti)) then
             write(6,*)'impossible to delete that many atoms of type ',iti
             stop
          end if
          do ivac=1,nvac(iti)
             ivactot=ivactot+1
             if (myidsp==0) then
                1 continue
                call random_number(z1)
                iatvac(ivac)=1+int(z1*natyp)
                do jvac=1,ivac-1
                   if (iatvac(jvac)==iatvac(ivac)) goto 1
                end do
#ifdef PARA
                if (lspacendm) then
                   looppr:do iproc=0,comm_space%nproc-1
                      if (last_at_typ(iproc).ge.iatvac(ivac)) then
                         ivacloc=iatvac(ivac)-last_at_typ(iproc-1)
                         exit looppr !iproc est l'indice du proc qui contient iatvac et ivacloc est le numéro de l'atome de cette lacune
                      end if
                   end do looppr
                end if
#else
                iproc=0
#endif
             end if
#ifdef PARA
             call comm_space%bcast(0,iproc)
             call comm_space%bcast(0,ivacloc)
             if (myidsp==iproc) then
#endif
                iat=0
                do i=1,atdml%im
                   if (atdml%ityp(i)==iti) then
                      iat=iat+1
                      if (iat==ivacloc) then
                         atomvac%iproc(ivactot)=myidsp 
                         atomvac%natg(ivactot)=atdml%num_at_glob(i)
                         atomvac%iloc(ivactot)=i
                         exit
                      end if
                   end if
                end do
#ifdef PARA
             end if
#endif
             
          end do
#ifdef PARA
          call comm_space%sum(atomvac%iproc) ! avant ça seul le proc iproc connaissait ces chiffres
          call comm_space%sum(atomvac%natg)
          call comm_space%sum(atomvac%iloc)
#endif
          deallocate (iatvac)
       end do
                
       do ivactot=1,nvactot
#ifdef PARA          
          if (comm_space%rank==atomvac%iproc(ivactot)) then
#endif
             call atdml%switch_atom(atomvac%iloc(ivactot),atdml%im)
             atdml%im=atdml%im-1
#ifdef PARA          
          end if
#endif
       end do
       atdml%im_glob=atdml%im_glob-nvactot
       call caltabtC(celndm,atdml,lperiod,boxndm)
#ifdef PARA
       call maj_atomes_frt_ftm(atdml,celndm,boxndm,psc)
#endif

    end do

  end subroutine creadp

!!$    ntry=0
!!$
!!$1      continue
!!$       ntry=ntry+1
!!$       ! tirer une position d'insertion
!!$
!!$
!!$
!!$       select case (typint)
!!$       case(0)
!!$
!!$          call random_number(z1)
!!$          iposI=1+Int(z1*nposI)
!!$          !        write(6,*)iposI
!!$          if (iposI.gt.nposI) idep=nposI
!!$          !iposI est l'indice de la position int.
!!$          !xposinttest est la position effective de l'int.
!!$
!!$          xposinttest(:)=xposint(:,iposI)
!!$          call cryst_to_cart (1, xposinttest(:), bg, -1) !cryst vers cart sur cv
!!$
!!$       case(1)
!!$          call random_number(z1)
!!$          call random_number(z2)
!!$          call random_number(z3)
!!$          xposinttest(1)=z1
!!$          xposinttest(2)=z2
!!$          xposinttest(3)=z3
!!$       end select
!!$
!!$       ! verifier qu'elle est loin de tout atome
!!$
!!$       call cryst_to_cart (imm, xp, bg, -1)    !cart vers cryst
!!$
!!$
!!$
!!$       do j=1,im
!!$          c1 = xposinttest(1)-xp(1,j)
!!$          c2 = xposinttest(2)-xp(2,j)
!!$          c3 = xposinttest(3)-xp(3,j)
!!$          if (c1>0.5) c1 = c1-1.
!!$          if (c1<(-0.5)) c1 = c1+1.
!!$          if (c2>0.5) c2 = c2-1.
!!$          if (c2<(-0.5)) c2 = c2+1.
!!$          if (c3>0.5) c3 = c3-1.
!!$          if (c3<(-0.5)) c3 = c3+1.
!!$
!!$          cv(1,1) = c1
!!$          cv(1,2) = c2
!!$          cv(1,3) = c3
!!$          call cryst_to_cart (1, cv, at, 1) !cryst vers cart sur cv
!!$          r2 = cv(1,1)*cv(1,1)+cv(1,2)*cv(1,2)+cv(1,3)*cv(1,3)
!!$          !           write(6,*)'toto',sqrt(r2),j
!!$
!!$
!!$          rd=sqrt(r2)
!!$          if (rd<dminins) then 
!!$             !           write(6,*)'rate', rd, dminins
!!$             call cryst_to_cart (imm, xp, at, 1)     !cryst vers cart
!!$             goto 1                ! position proche d'un atome
!!$          end if
!!$
!!$       end do
!!$       call cryst_to_cart (imm, xp, at, 1)     !cryst vers cart
!!$       call cryst_to_cart (1, xposinttest(:), at, 1) !cryst vers cart sur cv
!!$
!!$
!!$       ! deplacement acceptable
!!$       ! tirer un atome
!!$3      continue
!!$       call random_number(z1)   
!!$       idep=imin+Int(z1*(imax-imin-1))
!!$       if (ioxdef.ne.0) then
!!$          if (ioxdef.gt.0) then
!!$             if (ityp(idep).ne.ioxdef) goto 3
!!$          else
!!$             if (ityp(idep).eq.ioxdef) goto 3
!!$          end if
!!$       end if
!!$       !     select case (ioxdef)
!!$       !     case(0)
!!$       !     case(-2)
!!$       !        if (ityp(idep)==2) goto 3
!!$       !     case(2)
!!$       !        if (ityp(idep).ne.2) goto 3
!!$       !     case(-4)
!!$       !        if (ityp(idep)==4) goto 3
!!$       !     case(4)
!!$       !        if (ityp(idep).ne.4) goto 3
!!$       !     case(-5)
!!$       !        if (ityp(idep)==5) goto 3
!!$       !    case(5)
!!$       !       if (ityp(idep).ne.5) goto 3
!!$
!!$
!!$
!!$
!!$       xavant(:)=xp(:,idep)
!!$       xdec(:)=xp(:,idep)-xpp(:,idep)
!!$       xp(:,idep)=xposinttest(:)
!!$       xapres(:)=xposinttest(:)
!!$       xpp(:,idep)=xposinttest(:)-xdec(:)
!!$
!!$
!!$       write(6,*)'INTRDUCTION PF de type ', ty(ityp(idep))
!!$       write(6,*)'ntry',ntry
!!$       !     write(6,*)
!!$       write(6,*)'indice lac int', idep, iposI
!!$       !     write(6,*)
!!$       write(6,'(A,3F12.5)')'pos. lac.', xavant(1)*1.d8,xavant(2)*1.d8,xavant(3)*1.d8
!!$       write(6,'(A,3F12.5)')'pos. int.', xapres(1)*1.d8,xapres(2)*1.d8,xapres(3)*1.d8
!!$       write(6,*)
!!$
!!$
!!$
!!$
!!$    itapp=it
!!$    write(6,*)'outcdp'
!!$    im_glob=im
!!$
!!$    return
!!$  end subroutine creadp


end module cdp_mod
