! ****************************************************************
module rasmol_mod
  USE cryst_to_cart_mod,only: cryst_to_cart
  USE gen_com_m, ONLY:rang,ivisu,sigat,ldesinteg,lpkbar,im_glob,&
       &cunitP,it,lcasca,timel,unitP,at,fnam,bg,erg2ev,lenfnam,eatom,dmtype,umass
  USE var_pot, ONLY:ntyp,ntyp_buffer,ty,ty_buffer,cm_buffer,cm

!    USE paraneb_mod
    !USE tab_imm_m,only:num_at_glob,ityp,xp
    use atomconfig,only: atom_config,atom_config_d,atom_config_e
    use boxconfig,only:box_config
    implicit none
    integer, dimension(:), allocatable       :: ityp_buffer   ! temp/iorary store the types buffer when we
contains

  subroutine rasmol(atmol,boxmol,itapp,namefr,rty)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    !atmol et boxmol sont les configurations atomiques et de boite
    !itapp est l'itération (apprente) en cours
    !namefr est la racine nom du fichier (par défaut celui de name.in
    !rty est un tableau     character*3,intent(in), dimension(1:atmol%im),optional  :: rty qui donne les symboles des atomes. utile pour utiliser d'autres symboles que les symboles chimiques associés aux types des atomes. En l'absence de rty, on utilise les symboles des types des atomes.
    !ivisu dans gen_com_m : 1 :.mol, 4=.cfg ; 2=.xred
   
 

    USE T_kind_param_m, ONLY:  double
#ifdef PARA
    USE mpi
    USE mod_para,only:MPI_COMM_space,status,ierr,nprocs,myid,NDM_MPI_REAl_DOUBLE
#else
    USE mod_para,only:myid
#endif
    ! ****************************************************************

    implicit none



    integer,intent(in),optional  :: itapp
    class(atom_config),intent(in)::atmol
    type(box_config),intent(in)::boxmol
    character*3,intent(in), dimension(1:atmol%im),optional  :: rty
    character(len=*), optional ::namefr

    character*80::namef
    integer :: rgloc,im,imm,j,ic

    character*3, dimension(:), allocatable  :: tyw
    real(double)::xp(3,atmol%im)
    real(double),allocatable::sigat(:,:,:),eat(:)
    integer:: num_at_glob(atmol%im)
    integer::ityp(atmol%im)

#ifdef PARA
    integer :: iproc
    real(double), allocatable :: xp_loc(:,:)
    integer, allocatable      :: ityp_loc(:)
    integer, allocatable      :: num_at_glob_loc(:)
    character*3, allocatable      :: tyw_loc(:)
    real(double),allocatable::sigat_loc(:,:,:),eat_loc(:)
    integer :: im_loc
    integer :: proc_source
#endif

    integer :: i, luvisu, luvisu2, iti,lenfn2
    real(double) :: xp1, xp2, xp3,at(3,3),bg(3,3),pat
    character :: extension*9


    
    im=atmol%im; imm=atmol%imm
    at =boxmol%at*1d8 ; bg=boxmol%bg*1d-8
    xp(1:3,1:atmol%im)=atmol%xp(1:3,1:atmol%im)*1d8
    num_at_glob(1:atmol%im)=atmol%num_at_glob(1:atmol%im)
    ityp(1:im)=atmol%ityp(1:im)
    allocate(tyw(1:atmol%imm))
    tyw='000'
    !    do i=1,im
    !       write(6,*)i,atmol%ityp(i),ty(atmol%ityp(i))
    !    end do
    if (present (rty))then
       tyw(1:atmol%im)=rty(1:atmol%im)
    else
!       do i=1,im
!          write(6,*)i, atmol%ityp(i),ty(atmol%ityp(i))
!       end do
       tyw(1:atmol%im)=ty(atmol%ityp(1:atmol%im))
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


#ifdef PARA
    rgloc=myid
#else
    rgloc=0
#endif    
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
    luvisu = 86
    luvisu2 = 87
#ifdef PARA   
    luvisu = 86+rang
    luvisu2 = 87+rang
#endif    

    !      write (*, *) 'entree dans rasmol.f',im

    ! ****ouverture fichier sortie pour traitement images rasmol****
    !if(rgloc==0)       OPEN(luvisu,file='donnrasmol',form='formatted', &
    !       status='unknown')
    !   la fonction char ne marche que pour de faibles valeurs de it!!!!
    !     open(luvisu,file='donn_rasmolit.'//char(48+it),form='formatted',status='unknown')

    ! conversion entier-->alphanumerique par transfert du nombre
    ! de l'iteration vers fichier tampon relu sous format caractere.

    if(rgloc==0) then

       ! TJ: change the formatting so that files are well listed.

       lenfn2 = 9
       if (present(itapp))then
          if (itapp < 0  )   extension='iiiiiiiii' 
          if (itapp >= 0 )   write(extension,'(i9.9)') itapp
       end if



       select case (ivisu)
       case (1)
          if (present(namefr))then
             if (present(itapp)) then
                namef=namefr(1:len(namefr))//'.'//extension(1:lenfn2)//'.mol'
             else
                namef=namefr(1:len(namefr))//'.mol'
             end if
          else
             if (present(itapp)) then
                namef=fnam(1:lenfnam)//'.'//extension(1:lenfn2)//'.mol'
             else
                namef=fnam(1:lenfnam)//'.mol'
             end if
             open(luvisu, file=namef, form='formatted', &
                  status='unknown')

             write (luvisu, '(I9,A,I7,A,F12.6)') im_glob, ' IT =', itapp, ' Time = ', timel

             write (luvisu,'(9F12.6)')at(1,1),at(2,1),at(3,1),at(1,2),at(2,2),at(3,2),at(1,3),at(2,3),at(3,3)
             !at=at/1.d8
          end if
       case (2)
       case(3) 
          if (present(namefr))then
             if (present(itapp)) then
                namef=namefr(1:len(namefr))//'.'//extension(1:lenfn2)//'.xred'
             else
                namef=namefr(1:len(namefr))//'.xred'
             end if
          else
             if (present(itapp)) then
                namef=fnam(1:lenfnam)//'.'//extension(1:lenfn2)//'.xred'
             else
                namef=fnam(1:lenfnam)//'.xred'
             end if
          end if
          open(luvisu, file=namef, form='formatted', status='unknown')


          write (luvisu,'(9F12.6)')at(1,1),at(2,1),at(3,1),at(1,2),at(2,2),at(3,2),at(1,3),at(2,3),at(3,3)
          call cryst_to_cart (imm, xp,  bg,  -1) !cart vers cryst
       case(4) 
          if (present(namefr))then
             if (present(itapp)) then
                namef=namefr(1:len(namefr))//'.'//extension(1:lenfn2)//'.cfg'
             else
                namef=namefr(1:len(namefr))//'.cfg'
             end if
          else
             if (present(itapp)) then
                namef=fnam(1:lenfnam)//'.'//extension(1:lenfn2)//'.cfg'
             else
                namef=fnam(1:lenfnam)//'.cfg'
             end if
          end if
          open(luvisu, file=namef, form='formatted', status='unknown')
#ifdef PARA
          write(luvisu,'(a,i0)')'Number of particles = ', im_glob
#else
          write(luvisu,'(a,i0)')'Number of particles = ', im

#endif



          write(luvisu,'(a)')'A = 1.000 Angstrom (basic length-scale)'
          do j=1,3
             write(luvisu,'(a,i0)') '# Unit cell vector #', j    
             do ic=1,3
                write(luvisu,'(A,I1,A,I1,A,g16.8,A)')'H0(',j,',',ic,') = ',at(ic,j),' A'
             end do
          end do
          write(luvisu,'(A)')'.NO_VELOCITY.'              
          write(luvisu,'(A,I0)')'entry_count = ', 3
          call cryst_to_cart (im, xp,  bg,  -1) !cart vers cryst

       end select

    end if
    !if(itapp==0) open (file='filmtot.mol',unit=47)
    !if(itapp==0) open (file='filmtot.mol',unit=47)
    !      write (47, 134) im
    !      write (47, *) 'IT =', itapp, '    Time = ', timel



#ifdef PARA
    ! Le processeur maitre recoit les information des autres processeurs pour les ecrire sur fichier



    if (myid==0) then
       ! Copie des tableaux xp,num_at_glob et ityp locaux 
       allocate(xp_loc(3,imm))
       allocate(ityp_loc(imm))
       allocate(num_at_glob_loc(imm))
       allocate(tyw_loc(imm))
       xp_loc = xp
       ityp_loc = ityp
       num_at_glob_loc = num_at_glob
       im_loc = im
       tyw_loc=tyw
       select type (atmol)
       type is (atom_config_e)
          if (atmol%lprteat) then
             allocate (eat_loc(imm))
             eat_loc=eat
          end if
          if (atmol%lsigat)then
             allocate (sigat_loc(3,3,imm))
             sigat_loc=sigat
          end if
       end select
       ! Boucle sur les processeurs
       do iproc=0,nprocs-1
          ! Pour le processeur maitre il n'y a rien a faire
          ! reception des donnees des autres processeurs
          if (iproc.ne.0) then
             call MPI_RECV(im,               1,    MPI_INTEGER,      MPI_ANY_SOURCE, 10001, MPI_COMM_space, status, ierr)
             proc_source = status(MPI_SOURCE)
             call MPI_RECV(xp(1:3,1:im),     3*im, NDM_MPI_REAL_DOUBLE, proc_source, 10002, MPI_COMM_space, status, ierr)
             call MPI_RECV(ityp(1:im),       im,   MPI_INTEGER,         proc_source, 10003, MPI_COMM_space, status, ierr)
             call MPI_RECV(num_at_glob(1:im),im,   MPI_INTEGER,         proc_source, 10004, MPI_COMM_space, status, ierr)
             call MPI_RECV(tyw(1:im),3*im,   MPI_CHARACTER,         proc_source, 10004, MPI_COMM_space, status, ierr)
             select type (atmol)
             type is (atom_config_e)
                if (atmol%lprteat) call MPI_RECV(eat(1:im),im,  NDM_MPI_REAL_DOUBLE  ,  proc_source, 10005, MPI_COMM_space, status, ierr)
                if (atmol%lsigat) call MPI_RECV(sigat(:,:,1:im),9*im,  NDM_MPI_REAL_DOUBLE,proc_source, 10006, MPI_COMM_space, status, ierr)
             end select
!       select type (atmol)
!       type is (atom_config_e)
!          if (atmol%lprteat) call MPI_send(eat(1:im),im,  NDM_MPI_REAL_DOUBLE  ,         proc_source, 10005, MPI_COMM_space, status, ierr)
!          if (atmol%lsigat) call MPI_send(sigat(:,:,1:im),9*im,  NDM_MPI_REAL_DOUBLE  ,         proc_source, 10006, MPI_COMM_space, status, ierr)
!       end select

          endif
#endif
          do i = 1, im
             xp1 = xp(1,i)
             xp2 = xp(2,i)
             xp3 = xp(3,i)
             select case (ivisu)
             case(1)
!                write (6,*) 't',tyw(i)
!                write(6,*)'x', xp1,xp2, xp3
                write (luvisu, '(A,3f10.4)',advance='no') tyw(i),xp1, xp2, xp3
                select type (atmol)
                type is (atom_config_e)
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
             case(3)                    
                write (luvisu,'(3es15.6,2x,2a)') xp1, xp2, xp3, ' ! ', tyw(i)
             case(4)
                WRITE(luvisu,'(f0.3)') cm(ityp(i))/umass        ! Mass (g/mol)
                WRITE(luvisu,'(a)') tyw(i)                 ! Atom type
                write(luvisu, '(3(g24.16,1x))') xp(:,i)

             end select
             !               write (47, 135) tyw(i),xp1, xp2, xp3
             !end if


          end do

#ifdef PARA
       enddo  ! fin de boucle sur les processeurs
       ! Le processeur maitre recupere ses donnees locales
       xp = xp_loc
       ityp = ityp_loc
       num_at_glob = num_at_glob_loc
       im = im_loc
       deallocate(xp_loc)
       deallocate(ityp_loc)
       deallocate(num_at_glob_loc)
       select type (atmol)
       type is (atom_config_e)
          
          if (atmol%lprteat) then
             eat=eat_loc
             deallocate (eat_loc)
          end if
          if (atmol%lsigat)then
             sigat=sigat_loc
             deallocate (sigat_loc)
          end if
       end select
    else ! Les autres processeurs envoient leurs donnees locales
       call MPI_SEND(im,               1,   MPI_INTEGER,        0,10001,MPI_COMM_space,ierr)
       call MPI_SEND(xp(1:3,1:im),     3*im,NDM_MPI_REAL_DOUBLE,0,10002,MPI_COMM_space,ierr)
       call MPI_SEND(ityp(1:im),       im,  MPI_INTEGER,        0,10003,MPI_COMM_space,ierr)
       call MPI_SEND(num_at_glob(1:im),im,  MPI_INTEGER,        0,10004,MPI_COMM_space,ierr)
       select type (atmol)
       type is (atom_config_e)
          if (atmol%lprteat) call MPI_SEND(eat(1:im),im,  NDM_MPI_REAL_DOUBLE  ,    proc_source, 10005, MPI_COMM_space,  ierr)
          if (atmol%lsigat) call MPI_SEND(sigat(:,:,1:im),9*im,  NDM_MPI_REAL_DOUBLE  ,         proc_source, 10006, MPI_COMM_space,  ierr)
       end select
    endif
#endif


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

    if(rgloc==0)      close(luvisu)

    return
  end subroutine rasmol




end module rasmol_mod
