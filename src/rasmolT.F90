! ****************************************************************
module rasmolT_mod
  USE arret_ndm_mod,only:arret_ndm
  USE cryst_to_cart_mod,only: cryst_to_cart
  USE gen_com_m, ONLY:rang,ivisu,lpkbar,lspaceNDM,&
       &cunitP,iteration,lcasca,timel,unitP,fnam,erg2ev,lenfnam,umass,rang
  USE var_pot, ONLY:ntyp,ntyp_buffer,ty,ty_buffer,cm_buffer,cm

  use atomconfig,only: atom_config,atom_config_d,atom_config_e,atom_config_arps
  use paraconfig,only:para_config
  use boxconfig,only:box_config
  USE T_kind_param_m, ONLY:  double
  implicit none
  integer, dimension(:), allocatable       :: ityp_buffer   ! temp/iorary store the types buffer when

contains

  subroutine rasmolT(atmol,boxmol,itapp,namefr,rty,latcomp,ivisumol,naux,charaux,vaux,lappend)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    !atmol et boxmol sont les configurations atomiques et de boite
    !itapp est l'itération (apprente) en cours
    !namefr est la racine nom du fichier (par défaut celui de name.in
    !rty est un tableau     character*3,intent(in), dimension(1:atmol%im),optional  :: rty qui donne les symboles des atomes. utile pour utiliser d'autres symboles que les symboles chimiques associés aux types des atomes. En l'absence de rty, on utilise les symboles des types des atomes.
    !latcomp= en PARA latcomp=.true.=> atmol est une cofiguration complète/latcomp=false=>atmol est distributé sur comm_space
    !ivisu dans gen_com_m : 1 :.mol, 4=.cfg ; 2=.xred ; 5 =.gin
    !    naux=nb de carac auxiliaires,characaux string de description des carac ,vaux valeurs des auxiliaires

#ifdef PARA
    USE Tpara,only:COMM_space,nprocspace,myidsp,nprocs
#else
    USE Tpara,only:myidsp,nprocs
#endif

    implicit none
    integer,intent(in),optional  :: itapp
    class(atom_config)::atmol
    class(box_config),intent(in)::boxmol
    character*3,intent(in), dimension(1:atmol%im),optional,target  :: rty
    character(len=*), optional ::namefr
    logical::latcomp ! true= pas besoinde rapatrier atdml, false= il faut rapatrier atdml sur les masters
    integer, optional:: ivisumol
    logical,optional::lappend
    real(double)::at(3,3),bg(3,3)
    logical::lappendF

    integer,optional,intent(in)::naux
    character(len=*),optional::charaux(:)
    real(double),optional,intent(in)::vaux(:,:) 
    logical::laux


    integer::nauxw,nauxv,im,imm,im_glob,im_proc
    real(double),allocatable::vauxw(:,:) 
    character(len=10),allocatable::charauxw(:)

    integer::ivisum
    character*80::nameo,end_name
    integer :: rgloc,j,ic

    character*3, dimension(:), pointer  :: tyw

#ifdef PARA

    real(double),allocatable::xp_proc(:,:),aux_proc(:,:)
    integer,allocatable::ityp_proc(:),natg_proc(:)
    character*3, dimension(:), allocatable  :: ty_proc


#endif

#ifdef DKIO
    logical :: lvelocities=.false.
    character(len=6) :: format
    real(double),allocatable::xposg(:,:),vpg(:,:),masses_g(:)
    character*3, dimension(:),allocatable :: tywg
#endif

    !    class(atom_config),allocatable::atcomp
    integer :: i, luvisu, luvisu2,lenfn2,iaux,iaux2,i_proc,proc_source
    real(double) :: xp1, xp2, xp3,pat
    character :: extension*9
    integer::iax
    logical::latc

    latc=latcomp
    if (nprocs==1) latc=.true.
    if (present(lappend)) then
       lappendF=lappend
    else
       lappendF=.false.
    end if

    if((.not.latc).and.(atmol%im_glob==0)) then
       write(6,*)'rasmolT im_glob stop',latc,atmol%im_glob
       call arret_ndm
    end if

    if (present(ivisumol)) then
       ivisum=ivisumol
    else
       ivisum=ivisu
    end if
    laux=.false.
    nauxw=0 ; nauxv=0
#ifdef DKIO
    if ((ivisum==41).or.(ivisum==61).or.(ivisum==21).or.(ivisum==22).or.(ivisum==23).or.(ivisum==25))then
#else
    if ((ivisum==41).or.(ivisum==61))then
#endif
       laux=.true.
       nauxw=nauxw+3
    end if
    if (present(naux)) then
       if (naux.gt.0) then
          laux=.true.
          nauxw=nauxw+naux
          nauxv=nauxv+naux
       end if
    end if
    select type (atmol)
    class is (atom_config_e)
       if (atmol%lsigat) then
          laux=.true.
          nauxw=nauxw+1
          nauxv=nauxv+1
       end if
       if (atmol%lprteat) then
          laux=.true.
          nauxw=nauxw+1
          nauxv=nauxv+1
       end if
    end select
    block
      integer::jaux,jaux2
      if (laux) then 
         allocate (vauxw(nauxw,atmol%im))
         allocate(charauxw(nauxw))
#ifdef DKIO
         if ((ivisum==41).or.(ivisum==61).or.(ivisum==21).or.(ivisum==22).or.(ivisum==23).or.(ivisum==25))then
#else
         if ((ivisum==41).or.(ivisum==61))then
#endif
            charauxw(1)='VPX';          charauxw(2)='VPX';          charauxw(3)='VPZ'
            jaux=3
         else
            jaux=0
         end if
         if (present(naux)) then
            if (naux.gt.0) then
               do jaux2=1,naux
                  jaux=jaux+1
                  charauxw(jaux)=charaux(jaux2)
               end do
            end if
         end if
         select type (atmol)
         class is (atom_config_e)
            if (atmol%lsigat) then
               jaux=jaux+1
               charauxw(jaux)='sigat'
            end if
            if (atmol%lprteat) then
               jaux=jaux+1
               charauxw(jaux)='Eat'
            end if
         end select
      end if

    end block

    iaux=0
#ifdef DKIO
    if ((ivisum==41).or.(ivisum==61).or.(ivisum==21).or.(ivisum==22).or.(ivisum==23).or.(ivisum==25))then
#else
    if ((ivisum==41).or.(ivisum==61))then
#endif
       select type (atmol)
       type is(atom_config)
          write(6,*)' no velocity in atom-config and cfg with velocities stop'
          call arret_ndm
       class is (atom_config_d)
          laux = .true.
          !          nauxw = nauxw + 3
          do ic = 1, 3
             iaux = iaux + 1
             do i = 1, atmol%im
                vauxw(iaux, i) = atmol%vp(ic, i)*1d8*1d-12
             end do
          end do

       end select
    end if
    if (present(naux)) then
       if (naux.gt.0) then
          do iaux2=1,naux
             iaux=iaux+1
             do i=1,atmol%im
                vauxw(iaux,i)=vaux(iaux2,i)
             end do
          end do
       end if
    end if
    select type (atmol)
    class is (atom_config_e)
       if (atmol%lsigat) then
          iaux=iaux+1
          do i=1,atmol%im
             pat=unitP*(atmol%sigat(1,1,i)+atmol%sigat(2,2,i)+atmol%sigat(3,3,i))/3.
             vauxw(iaux,i)=pat
          end do
       end if
       if (atmol%lprteat) then
          iaux=iaux+1
          do i=1,atmol%im
             vauxw(iaux,i)=atmol%eat(i)*erg2ev
          end do
       end if
    end select



#ifdef PARA
    !    latcin=latc
    !    if (latcin) then 
    !       if (rang==0) then
    !          latcin=.true.
    !       else
    !          latcin=.false. !latcin intègre lw0 et rang=0
    !       end if
    !    end if
    rgloc=myidsp
    !    call atmol%deftype(atcomp)
    if (latc.eqv..false.) then
       if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
          im_glob =atmol%im_glob
          rgloc=myidsp
       else
          write(6,*)'latc=false et (nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) ??? stop'
          write(6,*)latc,nprocspace,lspaceNDM
          call arret_ndm

       end if
    else
       !          call atcomp%init(atmol%im)
       !          call atmol%copy_config(atcomp,lrescl=.true.)
       im_glob =atmol%im
       imm=atmol%imm
       rgloc=myidsp
       !im=atmol%im
       !imm=atmol%imm
    end if
#else
    !   call atmol%deftype(atcomp)
    rgloc=0
    !   call atcomp%init(atmol%im,ltabvois=.false.,nvois=0)
    !   call atmol%copy_config(atcomp,lrescl=.true.)
    im_glob=atmol%im
    imm=atmol%imm
#endif



    if(rgloc==0) then
       !*****************PPPPPPPPPAAAAAAAAASSSSSSSSAAAAAAAAAAAGGGGGGGGEEEEEEEEEE en AngSTROMS!!!!!!!!!!!!!!!!
       at =boxmol%at*1d8 ; bg=boxmol%bg*1d-8
       !       atcomp%xp(1:3,1:atcomp%im)=atcomp%xp(1:3,1:atcomp%im)*1d8


       if(lPkbar) then
          unitP=1.0d-9
          cunitP='kbar'
       else
          unitP=1.0
          cunitP='d/cm2'
       endif


       luvisu = 86
       luvisu2 = 87
#ifdef PARA   
       luvisu = 86+rang
       luvisu2 = 87+rang
#endif    

       ! TJ: change the formatting so that files are well listed.

       lenfn2 = 9
       if (present(itapp))then
          select case(itapp)
          case(-1)
             extension='iiiiiiiii' 
             !          if (itapp < 0  )
          case(999999999)
             extension='fffffffff'
          case default
             write(extension,'(i9.9)') itapp
          end select
       else
          extension='ooooooooo'
       end if

       if (present(namefr)) then
          nameo=namefr(1:len(namefr))
       else
          nameo=fnam(1:lenfnam)
       end if
       select case (ivisum)
       case(5)
          end_name='.newgin'
       case (1)
          end_name='.mol'
       case(3)
          end_name='.xred'
       case(40,41)
          end_name='.cfg'
       case(60,61)
          end_name='.xfg'
       case(7)
          end_name='.xyz'
#ifdef DKIO
!       case(10) ----------------------------------- Soon availible
!          ! Dk_io Abinit format
!          end_name='.'
       case(11,21)
          ! Dk_io Atomeye's CFG format
          end_name='.cfg'
          format='xfg'
          if (ivisum==21) lvelocities=.true.
       case(12,22)
          ! Dk_io Atomeye's extended CFG format
          end_name='.xfg'
          format='xfg'
          if (ivisum==22) lvelocities=.true.
       case(13,23)
          ! Dk_io CASTEP format
          end_name='.cell'
          format='castep'
          if (ivisum==23) lvelocities=.true.
       case(14)
          ! Dk_io CIF format
          end_name='.cif'
          format='cif'
       case(15,25)
          ! Dk_io DL_POLY format
          end_name='.CONFIG'
          format='dlpoly'
          if (ivisum==25) lvelocities=.true.
       case(16)
          ! Dk_io GULP format
          end_name='.gulp'
          format='gulp'
       case(17)
          ! Dk_io LAMMPS format
          end_name='.lmp'
          format='lammps'
       case(18)
          ! Dk_io VASP format
          end_name='.POSCAR'
          format='vasp'
       case(19)
          ! Dk_io XYZ format
          end_name='.xyz'
          format='xyz'
#endif
       case default
          write(6,*)'wrong ivisu',ivisum,ivisu
          call arret_ndm
       end select
       
       if (ivisum < 10 .or. ivisum > 30)then
          if (present(itapp))then
             call openfilemol( luvisu,nameo,end_name,extension)
          else
             call openfilemol( luvisu,nameo,end_name)
          end if

         !       if (ivisum.ne.5) then
         !          tyw='000'
         !    do i=1,im
         !       write(6,*)i,atmol%ityp(i),ty(atmol%ityp(i))
         !    end do
         !       if (naux.gt.0)then
         call write_header(ivisum,at,im_glob,luvisu, itapp,atmol,nauxv,nauxw,laux,charauxw)
         !       else
         !          call write_header(ivisum,at,im_glob,luvisu, itapp,atmol,nauxv,nauxw,laux)
         !       end if
       end if

    end if
    if (present (rty))then
       tyw=>rty
    else
       allocate(tyw(atmol%imm))
       tyw(1:atmol%im)=ty(atmol%ityp(1:atmol%im))
    end if
    !      end if
    select type (atmol)
    class is (atom_config_arps)
       do i=1,atmol%im
          select case (atmol%mov(i))
          case(0)
             tyw(i)=' Re'
          case(1)
             tyw(i)=' In'             
          case(2)
             tyw(i)=' Mo'
          end select
       end do
    end select
    
    if (latc) then
       if (ivisum < 10 .or. ivisum > 30) then
         call writepos(ivisum, im_glob,atmol%xp,tyw,atmol%ityp,atmol%num_at_glob,luvisu,boxmol%at,boxmol%bg,laux,nauxw,vauxw)
#ifdef DKIO
       else
         call dk_io_write(nameo,end_name,boxmol%at,boxmol%bg,trim(format),lvelocities,extension,atmol,tyw)
#endif
       end if
    else

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!§§LATC FALSE!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef PARA
       call comm_space%maxval(atmol%im,imm)
       allocate (xp_proc(3,imm))
       if (laux) then
          allocate (aux_proc(nauxw,imm))
       end if
       allocate (ityp_proc(imm))
       allocate (natg_proc(imm))
       allocate (ty_proc(imm))



       if (rgloc==0) then
          im_proc=atmol%im
          xp_proc(:,1:im_proc)=atmol%xp(:,1:im_proc)
          ityp_proc(1:im_proc)=atmol%ityp(1:im_proc)
          ty_proc(1:im_proc)=tyw(1:im_proc)
          if (laux) then
             aux_proc(1:nauxw,1:im_proc)=vauxw(1:nauxw,1:im_proc)
          end if
          natg_proc(1:im_proc)=atmol%num_at_glob(1:im_proc)
          
          if (ivisum < 10 .or. ivisum > 30) then
             call writepos(ivisum,im_proc,xp_proc,ty_proc,ityp_proc,natg_proc, luvisu,boxmol%at,boxmol%bg,laux,nauxw,aux_proc)
#ifdef DKIO
          else
             call reconstruction(atmol%im_glob,im_proc,natg_proc,xp_proc,ty_proc,ityp_proc,xposg,tywg,masses_g,vpg,lvelocities,aux_proc)
#endif
          end if

          do i_proc=1,nprocspace-1
             !             call comm_space%probe(11001,sourceout=proc_source)
             proc_source=i_proc
             call comm_space%recv(im_proc,proc_source,11001)
             call comm_space%recv(xp_proc(1:3,1:im_proc),proc_source,11002)
             call comm_space%recv(ityp_proc(1:im_proc),proc_source,11003)
             call comm_space%recv(natg_proc(1:im_proc),proc_source,11004)
             call comm_space%recv(ty_proc(1:im_proc),proc_source,11005)
             if (laux) then
                call comm_space%recv(aux_proc(1:nauxw,1:im_proc),proc_source,11006)
             end if
             if (ivisum < 10 .or. ivisum > 30) then
                call writepos(ivisum,im_proc,xp_proc,ty_proc,ityp_proc,natg_proc, luvisu,boxmol%at,boxmol%bg,laux,nauxw,aux_proc)
#ifdef DKIO
             else
                call reconstruction(atmol%im_glob,im_proc,natg_proc,xp_proc,ty_proc,ityp_proc,xposg,tywg,masses_g,vpg,lvelocities,aux_proc)
#endif
             end if
          enddo
#ifdef DKIO
          if (ivisum > 9 .and. ivisum < 31) then
             call dk_io_write(nameo,end_name,boxmol%at,boxmol%bg,trim(format),lvelocities,extension,xposg=xposg,masses_g=masses_g,vpg=vpg,tywg=tywg,img=atmol%im_glob)
             deallocate(xposg)
             deallocate(tywg)
             deallocate(masses_g)
             if (allocated(vpg)) deallocate(vpg)
          end if
#endif

       else ! myidsp different de 0 :
          im_proc=atmol%im
          call comm_space%send(im_proc,0,11001)
          call comm_space%send(atmol%xp(1:3,1:im_proc),0,11002)
          call comm_space%send(atmol%ityp(1:im_proc),0,11003)
          call comm_space%send(atmol%num_at_glob(1:im_proc),0,11004)
          call comm_space%send(tyw(1:im_proc),0,11005)
          if (laux) then
             call comm_space%send(vauxw(1:nauxw,1:im_proc),0,11006)
          end if

       endif

#endif
    end if
    !    write(6,*)'OUT Rasmol',rang

    !if(itapp==0) open (file='filmtot.mol',unit=47)

    select case (ivisum)
    case(7) !xyz
       write(luvisu,*)
       write (luvisu,'(3F15.9)')at(1,1),at(2,1),at(3,1)
       write (luvisu,'(3F15.9)')at(1,2),at(2,2),at(3,2)
       write (luvisu,'(3F15.9)')at(1,3),at(2,3),at(3,3)
       write (luvisu,*) 
    end select
    if ((.not.lappendF).and.(rgloc==0)) close(luvisu)
    return

    ! -------------------------------------------------------------
    !                          formats
    ! -------------------------------------------------------------
!!$134 format(i6)
!!$135 format(A,3f10.4,I9)
!!$136 format(A,3f10.4,D14.5,I9)
!!$
!!$200 format(i2)
!!$300 format(i3)
!!$400 format(i4)
!!$500 format(i5)
!!$600 format(i6)
!!$700 format(i7)
!!$800 format(i8)
!!$900 format(i9)
!!$101 format(a1)
!!$201 format(a2)
!!$301 format(a3)
!!$401 format(a4)
!!$501 format(a5)
!!$601 format(a6)
!!$701 format(a7)
!!$801 format(a8)
!!$901 format(a9)



  end subroutine rasmolT

  subroutine openfilemol(luvisu,nameo,end_name,ext)
    integer,intent(in)::luvisu
    character(len=*),intent(in)::nameo,end_name
    character(len=9),optional::ext
    logical::lopen
    character*80::namef
    !    write(6,*)'name o end_name ',nameo, ' ; ',end_name
    if (present(ext)) then
       namef=trim(nameo)//'.'//trim(ext)//trim(end_name)
       !              namef=nameo(1:len(nameo))//'.'//ext//'.'//end_name
    else
       namef=trim(nameo)//trim(end_name)
    end if
    !    write(6,*)'atomic output file name ',namef
    inquire(FILE=namef,opened=lopen)
    if (.not.lopen)open(luvisu, file=namef, form='formatted', &
         &         status='unknown')
  end subroutine openfilemol

#ifdef DKIO
  subroutine dk_io_write(nameo,end_name,box,invbox,format,lvelocities,ext,atcomp,tyw,xposg,masses_g,vpg,tywg,img)
    !-----------------------------------------------------
    !  Subroutine for interfacing with the dk_io library
    !-----------------------------------------------------
    use dk_structure_io, only: write_structure, TAG_LENGTH, Element
    use T_kind_param_m, only:  double
    
    character(len=*), intent(in) :: nameo,end_name,format
    real(double), intent(in) :: box(3,3), invbox(3,3)
    class(atom_config), optional :: atcomp
    character(len=3), dimension(:), intent(in), optional :: tyw
    real(double),allocatable, optional ::xposg(:,:),vpg(:,:), masses_g(:)
    character*3, dimension(:),allocatable, optional :: tywg
    character(len=9), intent(in), optional :: ext
    logical, intent(in) :: lvelocities
    character(len=80) :: namef
    real(double), allocatable, dimension(:) :: masses
    character(TAG_LENGTH), dimension(:), allocatable :: tags
    integer :: i
    integer, intent(in), optional :: img
    type(Element), dimension(:), allocatable :: elements

    if (present(ext)) then
       namef=trim(nameo)//'.'//trim(ext)//trim(end_name)
    else
       namef=trim(nameo)//trim(end_name)
    end if

    allocate(elements(size(ty)))
    do i=1, size(ty)
       elements(i)%tag=ty(i)
    end do

    if (present(atcomp)) then ! atcomp et tyw => en séquentiel
       allocate(tags(atcomp%im))
       allocate(masses(atcomp%im))

       ! Convert to fractional coordinates, and get atoms tag
       call cryst_to_cart(atcomp%im, atcomp%xp,  invbox,  -1) !cart vers cryst
       do i=1, atcomp%im
          tags(i) = tyw(i)
          masses(i) = cm(atcomp%ityp(i))/umass        ! Mass (g/mol)
       end do

       if(lvelocities) then
          select type (atcomp)
          type is(atom_config)
             write(6,*)'no velocity in atom-config and export asked with velocities stop'
             call arret_ndm
          class is (atom_config_d)
             call write_structure(trim(namef), box*1d8, atcomp%xp(:,1:atcomp%im), tags, format=format, velocities=atcomp%vp(:,1:atcomp%im)*1d8*1d-12, masses=masses, elements=elements)
          end select
       else
          call write_structure(trim(namef), box*1d8, atcomp%xp(:,1:atcomp%im), tags, format=format, masses=masses, elements=elements)
       end if

       call cryst_to_cart (atcomp%im, atcomp%xp,  box,  1) !cryst vers cart
       deallocate(masses)

    else ! xposg, masses_g, tywg, img et vpg => tableaux reconstruits sur proc 0
       allocate(tags(img))
       ! Convert to fractional coordinates, and get atoms tag
       call cryst_to_cart(img, xposg,  invbox,  -1) !cart vers cryst
       do i=1, img
          tags(i) = tywg(i)
       end do
       
       if(lvelocities) then
          call write_structure(trim(namef), box*1d8, xposg, tags, format=format, velocities=vpg, masses=masses_g, elements=elements)
       else
          call write_structure(trim(namef), box*1d8, xposg, tags, format=format, masses=masses_g, elements=elements)
       end if
       
    end if
    deallocate(tags)
    deallocate(elements)
  end subroutine dk_io_write

  subroutine reconstruction(img,im_proc,natg_proc,xpos_proc,tyw_proc,ityp_proc,xposg,tywg,masses_g,vpg,lvelocities,aux_proc)
    integer,intent(in)::img, im_proc, natg_proc(:), ityp_proc(:)
    real(double),allocatable::xposg(:,:),vpg(:,:), masses_g(:)
    character*3, dimension(:),allocatable :: tywg
    real(double), intent(in) :: xpos_proc(:,:)
    real(double), intent(in), optional :: aux_proc(:,:)
    character*3, dimension(:),intent(in) :: tyw_proc
    logical, intent(in) :: lvelocities
    integer :: i

    if (.not. allocated(xposg)) then
       write(6,*) "Attention : utilisation de dk-io en parallèle => reconstruction du tableau de position complet sur le proc 0"
       allocate(xposg(3,img))
       allocate(tywg(img))
       allocate(masses_g(img))
       if (lvelocities) then
         write(6,*) "reconstruction du tableau des vitesses complet sur le proc 0"
         allocate(vpg(3,img))
       end if
    end if

    do i=1, im_proc
       xposg(:,natg_proc(i))=xpos_proc(:,i)
       tywg(natg_proc(i))=tyw_proc(i)
       masses_g(natg_proc(i))=cm(ityp_proc(i))/umass        ! Mass (g/mol)
       if(lvelocities) then
          vpg(1,natg_proc(i))= aux_proc(1,i)
          vpg(2,natg_proc(i))= aux_proc(2,i)
          vpg(3,natg_proc(i))= aux_proc(3,i)
       end if
    end do
  end subroutine reconstruction
#endif

  subroutine writexred(im,xpos,tyw,luvisu,atw,bgw)
    integer,intent(in)::im
    real(double)::xpos(:,:)
    integer::nvar,luvisu
    character*3, dimension(:)  :: tyw(:)
    real(double)::atw(3,3),bgw(3,3)
    real(double) :: xp1, xp2, xp3

    integer::i
    call cryst_to_cart (im, xpos,  bgw,  -1) !cart vers cryst
    do i = 1, im
       xp1 = xpos(1,i)
       xp2 = xpos(2,i)
       xp3 = xpos(3,i)
       write (luvisu,'(3es15.6,2x,2a)') xp1, xp2, xp3, ' ! ', tyw(i)
    end do
    call cryst_to_cart (im, xpos,  atw,  1) !cart vers cryst
  end subroutine writexred

  subroutine writegin(im,xpos,itypw,natg,luvisu,atw,bgw,laux,nvar,auxvar)
    integer,intent(in)::im
    real(double)::xpos(:,:)
    integer,intent(in),optional::nvar
    integer,intent(in)::luvisu
    real(double),intent(in),optional::auxvar(:,:)
    logical,optional::laux
    integer,intent(in)::itypw(:)
    integer,intent(in)::natg(:)
    real(double),intent(in)::atw(3,3),bgw(3,3)

    logical::lvar=.false.
    integer::i,iax
    real(double) :: xp1, xp2, xp3
    if (present(laux))lvar=laux
    call cryst_to_cart (im, xpos,  bgw,  -1) !cart vers cryst
    do i = 1, im
       xp1 = xpos(1,i)
       xp2 = xpos(2,i)
       xp3 = xpos(3,i)
       write (luvisu,'(3es15.6,I3, I9)',advance='no') xp1, xp2, xp3, itypw(i),natg(i)
       if (lvar) then
          if (nvar.gt.0) then
             do iax=1,nvar
                write(luvisu,'(G20.12)',advance='no')auxvar(iax,i)
             end do
          end if
       end if
       write(luvisu,*)
    end do
    call cryst_to_cart (im, xpos,  atw,  1) !cart vers cryst
  end subroutine writegin


  subroutine writemol(im,xpos,tyw,natg,luvisu,atw,bgw,laux,nvar,auxvar)
    integer,intent(in)::im
    real(double)::xpos(:,:)
    integer,intent(in),optional::nvar
    integer,intent(in)::luvisu
    real(double),intent(in),optional::auxvar(:,:)
    logical,optional::laux
    character*3, dimension(:),intent(in)  :: tyw(:)
    integer,intent(in)::natg(:)
    real(double),intent(in)::atw(3,3),bgw(3,3)

    logical::lvar=.false.
    integer::i,iax
    real(double) :: xp1, xp2, xp3
    if (present(laux))lvar=laux
    do i = 1, im
       xp1 = xpos(1,i)*1d8
       xp2 = xpos(2,i)*1d8
       xp3 = xpos(3,i)*1d8

       write (luvisu, '(A,3G18.4)',advance='no') tyw(i),xp1, xp2, xp3
       write (luvisu, '(I9)',advance='no')  natg(i)
       if (lvar) then
          if (nvar.gt.0) then
             do iax=1,nvar
                write(luvisu,'(G20.12)',advance='no')auxvar(iax,i)
             end do
          end if
       end if
       write(luvisu,*)
    end do
  end subroutine writemol

  subroutine writexyz(im,xpos,tyw,luvisu,atw,bgw,laux,nvar,auxvar)

    integer,intent(in)::im
    real(double)::xpos(:,:)
    integer,intent(in),optional::nvar
    integer,intent(in)::luvisu
    real(double),intent(in),optional::auxvar(:,:)
    logical,optional::laux
    character*3, dimension(:),intent(in)  :: tyw(:)
    real(double),intent(in)::atw(3,3),bgw(3,3)

    logical::lvar=.false.
    integer::i,iax
    real(double) :: xp1, xp2, xp3
    if (present(laux))lvar=laux
    do i = 1, im
       xp1 = xpos(1,i)*1d8
       xp2 = xpos(2,i)*1d8
       xp3 = xpos(3,i)*1d8

       write (luvisu, '(A,3G18.4)',advance='no') tyw(i),xp1, xp2, xp3
       if (lvar) then
          if (nvar.gt.0) then
             do iax=1,nvar
                write(luvisu,'(G20.12)',advance='no')auxvar(iax,i)
             end do
          end if
       end if
       write(luvisu,*)
    end do
  end subroutine writexyz

  subroutine write40(im,xpos,tyw,itypw,natg,luvisu,atw,bgw,laux,nvar,auxvar)
    integer,intent(in)::im
    real(double)::xpos(:,:)
    integer,intent(in),optional::nvar
    integer,intent(in)::luvisu
    real(double),intent(in),optional::auxvar(:,:)
    logical,optional::laux
    character*3, dimension(:),intent(in)  :: tyw(:)
    integer,intent(in)::itypw(:)
    integer,intent(in)::natg(:)
    real(double),intent(in)::atw(3,3),bgw(3,3)

    logical::lvar=.false.
    integer::i,iax
    if (present(laux))lvar=laux

    call cryst_to_cart (im, xpos,  bgw,  -1) !cart vers cryst
    do i = 1, im
       WRITE(luvisu,'(f0.3)') cm(itypw(i))/umass        ! Mass (g/mol)
       WRITE(luvisu,'(a,A,I9)') tyw(i),' ! ',natg(i)                 ! Atom type
       write(luvisu, '(3(g20.12,1x))',advance='no') xpos(:,i)
       if (lvar) then
          if (nvar.gt.0) then
             do iax=1,nvar
                write(luvisu,'(G20.12)',advance='no')auxvar(iax,i)
             end do
          end if
       end if
       write(luvisu,*)


    end do

    call cryst_to_cart (im, xpos,  atw,  1) !cart vers cryst
  end subroutine write40

  subroutine write60(im,xpos,tyw,itypw,luvisu,atw,bgw,laux,nvar,auxvar)    
    integer,intent(in)::im
    real(double)::xpos(:,:)
    integer,intent(in),optional::nvar
    integer,intent(in)::luvisu
    real(double),intent(in),optional::auxvar(:,:)
    logical,optional::laux
    character*3, dimension(:),intent(in)  :: tyw(:)
    integer,intent(in)::itypw(:)
    real(double),intent(in)::atw(3,3),bgw(3,3)

    logical::lvar=.false.
    integer::i,iax
    if (present(laux))lvar=laux

    call cryst_to_cart (im, xpos,  bgw,  -1) !cart vers cryst
    do i = 1, im

       if (i==1) then
          WRITE(luvisu,'(f0.3)')cm(itypw(i))/umass        ! Mass (g/mol)
          WRITE(luvisu,'(a)') tyw(i)
       else
          if (itypw(i).ne.itypw(i-1)) then
             WRITE(luvisu,'(f0.3)') cm(itypw(i))/umass        ! Mass (g/mol)      ! Mass (g/mol)
             WRITE(luvisu,'(a)') tyw(i)
          end if
       end if
       write(luvisu, '(3(g18.10,1x))',advance='no') xpos(:,i)
       if (lvar) then
          if (nvar.gt.0) then
             do iax=1,nvar
                write(luvisu,'(E15.7)',advance='no')auxvar(iax,i)
             end do
          end if
       end if

       write(luvisu,*)


    end do
    call cryst_to_cart (im, xpos,  atw,  1) !cart vers cryst

  end subroutine write60

  subroutine write_header(ivisum,at,im_g,luvisu,itapp,atmol,nauxv,nauxw,laux,charaux)
    class (atom_config),intent(in)::atmol
    integer,intent(in)::ivisum,luvisu,im_g,nauxv,nauxw
    integer,optional::itapp
    real(double),intent(in)::at(3,3)
    character(len=*),optional::charaux(:)
    logical,intent(in)::laux

    integer::iax,j,ic,iax2

    select case (ivisum)

    case(5) !newgin
       write (luvisu,'(A)',advance='no')' 1 1 1 !'
       if (laux) then
          if (present(charaux)) then
             do iax=1,nauxw
                write(luvisu,'(A)',advance='no')trim(charaux(iax))
             end do
          end if
       end if
       write(luvisu,*)
       write (luvisu,'(3F15.9)')at(1,1),at(2,1),at(3,1)
       write (luvisu,'(3F15.9)')at(1,2),at(2,2),at(3,2)
       write (luvisu,'(3F15.9)')at(1,3),at(2,3),at(3,3)
       write (luvisu,*) im_g
    case (1) !mol

       if (present(itapp))then
          write (luvisu, '(I9,A,I12,A,F12.6)',advance='no') im_g, ' IT =', itapp, ' Time = ', timel
       else
          write (luvisu, '(I9,A,I7,A,F12.6)',advance='no') im_g
       end if
       if (laux) then
          if (present(charaux)) then
             do iax=1,nauxw
                write(luvisu,'(A)',advance='no')trim(charaux(iax))
             end do
          end if
       end if
       write(luvisu,*)

       write (luvisu,'(9F12.6)')at(1,1),at(2,1),at(3,1),at(1,2),at(2,2),at(3,2),at(1,3),at(2,3),at(3,3)
    case(3) 
       write (luvisu,'(9F12.6)')at(1,1),at(2,1),at(3,1),at(1,2),at(2,2),at(3,2),at(1,3),at(2,3),at(3,3)

    case(40,41,60,61) 

       write(luvisu,'(a,i0)')'Number of particles = ',  im_g
       write(luvisu,'(a)')'A = 1.000 Angstrom (basic length-scale)'
       do j=1,3
          write(luvisu,'(a,i0)') '# Unit cell vector #', j    
          do ic=1,3
             write(luvisu,'(A,I1,A,I1,A,g16.8,A)')'H0(',j,',',ic,') = ',at(ic,j),' A'
          end do
       end do
       if ((ivisum==40).or.(ivisum==60)) then
          write(luvisu,'(A)')'.NO_VELOCITY.'
          write(luvisu,'(A,I0)')'entry_count = ', 3+nauxw
          if ((laux).and.(present(charaux))) then
             do iax=0,nauxw-1
                iax2=iax
                write(luvisu,'(A,I0,A,A)')'auxiliary[',iax2,'] = ',trim(charaux(iax+1))
             end do
          end if
       else
          write(luvisu,'(A,I0)')'entry_count = ', 3+nauxw
          if ((laux).and.(present(charaux))) then
             do iax=0,nauxw-1
                iax2=iax
                write(luvisu,'(A,I0,A,A)')'auxiliary[',iax2,'] = ',trim(charaux(iax+1))
             end do
          end if
       end if
       iax=nauxV
    end select
  end subroutine write_header

  subroutine writepos(ivisum,im,xpos,tyw,itypw,natg,luvisu,atw,bgw,laux,nvar,auxvar)
    integer,intent(in)::im,ivisum
    real(double)::xpos(:,:)
    integer,intent(in),optional::nvar
    integer,intent(in)::luvisu
    real(double),intent(in),optional::auxvar(:,:)
    logical,optional::laux
    character*3, dimension(:),intent(in)  :: tyw(:)
    integer,intent(in)::itypw(:)
    integer,intent(in)::natg(:)
    real(double),intent(in)::atw(3,3),bgw(3,3)
    select case (ivisum)
    case(3)
       call writexred(im,xpos,tyw,luvisu,atw,bgw)
    case(5)
       call writegin(im,xpos,itypw,natg,luvisu,atw,bgw,laux,nvar,auxvar)
    case(1)
       call writemol(im,xpos,tyw,natg,luvisu,atw,bgw,laux,nvar,auxvar)
    case(7)
       call writexyz(im,xpos,tyw,luvisu,atw,bgw,laux,nvar,auxvar)
    case(40,41)
       call write40(im,xpos,tyw,itypw,natg,luvisu,atw,bgw,laux,nvar,auxvar)
    case(60,61)
       call write60(im,xpos,tyw,itypw,luvisu,atw,bgw,laux,nvar,auxvar)
    end select
  end subroutine writepos

end module rasmolT_mod
