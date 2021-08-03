! ****************************************************************
module rasmolT_mod
  USE cryst_to_cart_mod,only: cryst_to_cart
  USE gen_com_m, ONLY:rang,ivisu,ldesinteg,lpkbar,im_glob,lspaceNDM,&
       &cunitP,it,lcasca,timel,unitP,fnam,erg2ev,lenfnam,dmtype,umass,rang
  USE var_pot, ONLY:ntyp,ntyp_buffer,ty,ty_buffer,cm_buffer,cm

  use atomconfig,only: atom_config,atom_config_d,atom_config_e
  use paraconfig,only:para_config
  use boxconfig,only:box_config
  implicit none
  integer, dimension(:), allocatable       :: ityp_buffer   ! temp/iorary store the types buffer when we
contains

  subroutine rasmolT(atmol,boxmol,itapp,namefr,rty,latcomp)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    !atmol et boxmol sont les configurations atomiques et de boite
    !itapp est l'itération (apprente) en cours
    !namefr est la racine nom du fichier (par défaut celui de name.in
    !rty est un tableau     character*3,intent(in), dimension(1:atmol%im),optional  :: rty qui donne les symboles des atomes. utile pour utiliser d'autres symboles que les symboles chimiques associés aux types des atomes. En l'absence de rty, on utilise les symboles des types des atomes.
    !latcomp= en PARA latcomp=.true.=> atmol est une cofiguration complète/latcomp=false=>atmol est distributé sur comm_space
    !lw0= .true. supprimé seul le proc 0 écrit la configuration
    !ivisu dans gen_com_m : 1 :.mol, 4=.cfg ; 2=.xred ; 5 =.gin



    USE T_kind_param_m, ONLY:  double
#ifdef PARA
    USE Tpara,only:COMM_space,nprocspace,myidsp
#else
    USE Tpara,only:myidsp
#endif
    ! ****************************************************************

    implicit none



    integer,intent(in),optional  :: itapp
    class(atom_config),intent(in)::atmol
    type(box_config),intent(in)::boxmol
    character*3,intent(in), dimension(1:atmol%im),optional  :: rty
    character(len=*), optional ::namefr
    logical,intent(in)::latcomp ! true= pas besoinde rapatrier atdml, false= il faut rapatrier atdml sur les masters

    character*80::namef,nameo,end_name
    integer :: rgloc,im,imm,j,ic,e_c,e_c0

    character*3, dimension(:), allocatable  :: tyw
    real(double),allocatable::xp(:,:)
    real(double),allocatable::sigat(:,:,:),eat(:)
    integer,allocatable:: num_at_glob(:)
    integer,allocatable::ityp(:)

#ifdef PARA
    integer :: iproc
    real(double), allocatable :: xp_loc(:,:)
    integer, allocatable      :: ityp_loc(:)
    integer, allocatable      :: num_at_glob_loc(:)
    character*3, allocatable      :: tyw_loc(:)
    real(double),allocatable::sigat_loc(:,:,:),eat_loc(:)
    integer :: im_loc
    integer :: proc_source
    type(para_config)::div

#endif
    type(atom_config_e)::atcomp
    integer :: i, luvisu, luvisu2, iti,lenfn2
    real(double) :: xp1, xp2, xp3,at(3,3),bg(3,3),pat
    character :: extension*9

#ifdef PARA
    !    latcompin=latcomp
    !    if (latcompin) then 
    !       if (rang==0) then
    !          latcompin=.true.
    !       else
    !          latcompin=.false. !latcompin intègre lw0 et rang=0
    !       end if
    !    end if
    rgloc=myidsp

    if (latcomp.eqv..false.) then

       if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
          call atcomp%init(im_glob)
          div%mpi_image%rank=myidsp
          div%mpi_image%nproc=nprocspace
          div%mpi_image%comm=COMM_space%comm
          call atmol%vers_master(atcomp,div,'ixnlusv')
          im =atcomp%im
          imm=atcomp%im
          rgloc=myidsp
       else
          write(6,*)'latcomp=false et (nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) ??? stop'
          write(6,*)latcomp,nprocspace,lspaceNDM
          stop
       end if
    else
       if (rgloc==0) then
          call atcomp%init(atmol%im)
          call atmol%copy_config(atcomp,lrescl=.true.)
          im=atmol%im
          imm=atmol%imm
       end if
    end if
#else
    rgloc=0
    call atcomp%init(atmol%im,ltabvois=.false.,nvois=0)
    call atmol%copy_config(atcomp,lrescl=.true.)
    im=atmol%im
    imm=atmol%imm
#endif  

    if(rgloc==0) then

       at =boxmol%at*1d8 ; bg=boxmol%bg*1d-8
       allocate(xp(3,im));allocate(ityp(im));allocate(num_at_glob(imm))
       xp(1:3,1:atmol%im)=atcomp%xp(1:3,1:atmol%im)*1d8
       num_at_glob(1:atmol%im)=atcomp%num_at_glob(1:atmol%im)
       ityp(1:im)=atcomp%ityp(1:im)
       !       call atcomp%print
       allocate(tyw(imm))
       tyw='000'
       !    do i=1,im
       !       write(6,*)i,atmol%ityp(i),ty(atmol%ityp(i))
       !    end do
       if (present (rty))then
          tyw(1:im)=rty(1:im)
       else
          !       do i=1,im
          !          write(6,*)i, atmol%ityp(i),ty(atmol%ityp(i))
          !       end do
          tyw(1:im)=ty(ityp(1:im))
       end if


       select type (atmol)
       type is (atom_config_e)
          if(atmol%lsigat) then
             allocate (sigat(3,3,imm))
             sigat=0
             sigat(:,:,1:im)=atmol%sigat(:,:,:1:im)
          end if

          if(atmol%lprteat) then

             allocate (eat(im)) ; eat=0; eat(1:im)=atmol%eat(1:im)
          end if
       end select



       ! Notes about V_sim:
       ! * works if at(:,:) "encompasses" all the system (no duplication of lattice cells)
       ! * at(:,1) must be along x and at(:,2) must have no component along z.
       !   Otherwise a rotation matrix should be coded.
       !-----------------------------------------------
       !
       !
       !    if (dmtype==17) then 
       !       call redefine_ty() 
       !    end if
       !
       if(lPkbar) then
          unitP=1.0d-9
          cunitP='kbar'
       else
          unitP=1.0
          cunitP='d/cm2'
       endif

       !    iksp=-1
       !    if (lcasca) iksp=iko
       !    if (ldesinteg)iksp=1

       !      write (*, *) 'entree dans rasmol.f',im

       ! ****ouverture fichier sortie pour traitement images rasmol****
       !if(rgloc==0)       OPEN(luvisu,file='donnrasmol',form='formatted', &
       !       status='unknown')
       !   la fonction char ne marche que pour de faibles valeurs de it!!!!
       !     open(luvisu,file='donn_rasmolit.'//char(48+it),form='formatted',status='unknown')

       ! conversion entier-->alphanumerique par transfert du nombre
       ! de l'iteration vers fichier tampon relu sous format caractere.

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
       select case (ivisu)
       case(5)
          end_name='.newgin'
       case (1)
          end_name='.mol'
       case(3)
          end_name='.xred'
       case(4)
          end_name='.cfg'
       case(40,41,42)
          end_name='.xfg'
       end select
       if (present(itapp))then
          call openfilemol( luvisu,nameo,end_name,extension)
       else
          call openfilemol( luvisu,nameo,end_name)
       end if

       select case (ivisu)
       case(5)
          write (luvisu,*)' 1 1 1 '
          write (luvisu,'(3F12.6)')at(1,1),at(2,1),at(3,1)
          write (luvisu,'(3F12.6)')at(1,2),at(2,2),at(3,2)
          write (luvisu,'(3F12.6)')at(1,3),at(2,3),at(3,3)
          write (luvisu,*) atcomp%im
          call cryst_to_cart (im, xp,  bg,  -1) !cart vers cryst
          do i = 1, im
             xp1 = xp(1,i)
             xp2 = xp(2,i)
             xp3 = xp(3,i)
             write (luvisu,'(3es15.6,I3)') xp1, xp2, xp3, ityp(i)
          end do
       case (1)

          if (present(itapp))then
             write (luvisu, '(I9,A,I7,A,F12.6)') im, ' IT =', itapp, ' Time = ', timel
          else
             write (luvisu, '(I9,A,I7,A,F12.6)') im
          end if

          write (luvisu,'(9F12.6)')at(1,1),at(2,1),at(3,1),at(1,2),at(2,2),at(3,2),at(1,3),at(2,3),at(3,3)
          !at=at/1.d8
          do i = 1, im
             xp1 = xp(1,i)
             xp2 = xp(2,i)
             xp3 = xp(3,i)
             !                write (6,*) 't',tyw(i)
             !                write(6,*)'x', xp1,xp2, xp3
             write (luvisu, '(A,3f10.4)',advance='no') tyw(i),xp1, xp2, xp3
             select type (atmol)
             class is (atom_config_e)
                if (atmol%lsigat) then
                   if(it.eq.0)then
                      pat=0.0
                   else
                      pat=unitP*(sigat(1,1,i)+sigat(2,2,i)+sigat(3,3,i))/3.
                   end if
                   write (luvisu, '(D14.5)',advance='no') pat
                end if
                if (atmol%lprteat) write (luvisu, '(D14.5)',advance='no') eat(i)*erg2ev
             end select
#ifdef PARA
             write (luvisu, '(I9)')  num_at_glob(i)
#else
             write (luvisu, '(I9)')  i
#endif
          end do
!       case (2)
       case(3) 
          write (luvisu,'(9F12.6)')at(1,1),at(2,1),at(3,1),at(1,2),at(2,2),at(3,2),at(1,3),at(2,3),at(3,3)
          call cryst_to_cart (imm, xp,  bg,  -1) !cart vers cryst
          do i = 1, im
             xp1 = xp(1,i)
             xp2 = xp(2,i)
             xp3 = xp(3,i)
             write (luvisu,'(3es15.6,2x,2a)') xp1, xp2, xp3, ' ! ', tyw(i)
          end do
          
       case(4,40,41,42) 

          write(luvisu,'(a,i0)')'Number of particles = ', im
          write(luvisu,'(a)')'A = 1.000 Angstrom (basic length-scale)'
          do j=1,3
             write(luvisu,'(a,i0)') '# Unit cell vector #', j    
             do ic=1,3
                write(luvisu,'(A,I1,A,I1,A,g16.8,A)')'H0(',j,',',ic,') = ',at(ic,j),' A'
             end do
          end do
          if (ivisu.le.40) then
             write(luvisu,'(A)')'.NO_VELOCITY.'
             write(luvisu,'(A,I0)')'entry_count = ', 3
          else
             e_c=6
             select case(ivisu)
             case(41) !RAS
             case(42)
                select type (atmol)
                class is (atom_config_e)
                   if (atmol%lprteat) e_c=e_c+1
                   if (atmol%lsigat) e_c=e_c+9
                end select
             end select
             write(luvisu,'(A,I0)')'entry_count = ', e_c
             select case(ivisu)
             case(41) !RAS
             case(42)
                e_c0=-1
                select type (atmol)
                class is (atom_config_e)
                   if (atmol%lprteat) then
                      e_c0=e_c0+1
                      write(luvisu,'(A,I0,A)')'auxiliary[',e_c0,'] = eat'
                   end if
                   if (atmol%lsigat)then
                      e_c0=e_c0+1
                      write(luvisu,'(A,I0,A)')'auxiliary[',e_c0,'] = sigxx'
                      e_c0=e_c0+1
                      write(luvisu,'(A,I0,A)')'auxiliary[',e_c0,'] = sigxy'
                      e_c0=e_c0+1
                      write(luvisu,'(A,I0,A)')'auxiliary[',e_c0,'] = sigxz'
                      e_c0=e_c0+1
                      write(luvisu,'(A,I0,A)')'auxiliary[',e_c0,'] = sigyx'
                      e_c0=e_c0+1
                      write(luvisu,'(A,I0,A)')'auxiliary[',e_c0,'] = sigyy'
                      e_c0=e_c0+1
                      write(luvisu,'(A,I0,A)')'auxiliary[',e_c0,'] = sigyz'
                      e_c0=e_c0+1
                      write(luvisu,'(A,I0,A)')'auxiliary[',e_c0,'] = sigzx'
                      e_c0=e_c0+1
                      write(luvisu,'(A,I0,A)')'auxiliary[',e_c0,'] = sigzy'
                      e_c0=e_c0+1
                      write(luvisu,'(A,I0,A)')'auxiliary[',e_c0,'] = sigzz'
                   end if
                end select
             end select
          end if
          call cryst_to_cart (im, xp,  bg,  -1) !cart vers cryst
          do i=1,im
             select case(ivisu)
             case(4)
                WRITE(luvisu,'(f0.3)') cm(ityp(i))/umass        ! Mass (g/mol)
                WRITE(luvisu,'(a)') tyw(i)                 ! Atom type
                write(luvisu, '(3(g20.12,1x))') xp(:,i)
             case(40,41,42)
                if (i==1) then
                   WRITE(luvisu,'(f0.3)') cm(ityp(i))/umass        ! Mass (g/mol)
                   WRITE(luvisu,'(a)') tyw(i)                 ! Atom type
                else
                   if (ityp(i).ne.ityp(i-1)) then
                      WRITE(luvisu,'(f0.3)') cm(ityp(i))/umass        ! Mass (g/mol)
                      WRITE(luvisu,'(a)') tyw(i)                 ! Atom type
                   end if
                end if
                write (luvisu, '(3(g20.12,1x))',advance='no') xp(:,i)
                select case(ivisu)
                case(41)
                   write (luvisu, '(3(g20.12,1x))',advance='no') atcomp%vp(:,i)*1d8*1d-12
                case(42)
                   write (luvisu, '(3(g20.12,1x))',advance='no') atcomp%vp(:,i)*1d8*1d-12
                   select type (atmol)
                   class is (atom_config_e)
                      if (atmol%lprteat) then
                         write(luvisu,'(g20.12)',advance='no')atcomp%eat*erg2ev
                      end if
                      if (atmol%lsigat)then
                         write(luvisu,'(9g20.12)',advance='no')atcomp%sigat(1,1,i),atcomp%sigat(1,2,i),atcomp%sigat(1,3,i),&
                              &atcomp%sigat(2,1,i),atcomp%sigat(2,2,i),atcomp%sigat(2,3,i),&
                              &atcomp%sigat(3,1,i),atcomp%sigat(3,2,i),atcomp%sigat(3,3,i)
                      end if
                   end select
                end select
                write(luvisu,'(A)')' '
             end select
          end do
       end select
       
       close(luvisu)
    end if


    !if(itapp==0) open (file='filmtot.mol',unit=47)


    ! -------------------------------------------------------------
    !                          formats
    ! -------------------------------------------------------------
134 format(i6)
135 format(A,3f10.4,I9)
136 format(A,3f10.4,D14.5,I9)

200 format(i2)
300 format(i3)
400 format(i4)
500 format(i5)
600 format(i6)
700 format(i7)
800 format(i8)
900 format(i9)
101 format(a1)
201 format(a2)
301 format(a3)
401 format(a4)
501 format(a5)
601 format(a6)
701 format(a7)
801 format(a8)
901 format(a9)



    return
  end subroutine rasmolT

  subroutine openfilemol(luvisu,nameo,end_name,ext)
    integer,intent(in)::luvisu
    integer::lenfn2
    character(len=*),intent(in)::nameo,end_name
    character(len=9),optional::ext

    character*80::namef
    if (present(ext)) then
       namef=trim(nameo)//'.'//trim(ext)//trim(end_name)
!              namef=nameo(1:len(nameo))//'.'//ext//'.'//end_name
    else
       namef=trim(nameo)//'.'//trim(end_name)
    end if
    write(6,*)'atom config file name ',namef
    open(luvisu, file=namef, form='formatted', &
         status='unknown')
  end subroutine openfilemol
end module rasmolT_mod
