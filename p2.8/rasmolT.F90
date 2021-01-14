! ****************************************************************
module rasmolT_mod
  USE cryst_to_cart_mod,only: cryst_to_cart
  USE gen_com_m, ONLY:rang,ivisu,ldesinteg,lpkbar,im_glob,lspaceNDM,&
       &cunitP,it,lcasca,timel,unitP,fnam,erg2ev,lenfnam,dmtype,umass
  USE var_pot, ONLY:ntyp,ntyp_buffer,ty,ty_buffer,cm_buffer,cm

  use atomconfig,only: atom_config,atom_config_d,atom_config_e
  use paraconfig,only:para_config
  use boxconfig,only:box_config
  implicit none
  integer, dimension(:), allocatable       :: ityp_buffer   ! temp/iorary store the types buffer when we
contains

  subroutine rasmolT(atmol,boxmol,itapp,namefr,rty,latcomp,lw0)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    !atmol et boxmol sont les configurations atomiques et de boite
    !itapp est l'itération (apprente) en cours
    !namefr est la racine nom du fichier (par défaut celui de name.in
    !rty est un tableau     character*3,intent(in), dimension(1:atmol%im),optional  :: rty qui donne les symboles des atomes. utile pour utiliser d'autres symboles que les symboles chimiques associés aux types des atomes. En l'absence de rty, on utilise les symboles des types des atomes.
    !latcomp= en PARA latcomp=.true.=> atmol est une cofiguration complète/latcomp=false=>atmol est distributé sur comm_space
    !lw0= .true. seul le proc 0 écrit la configuration
    !ivisu dans gen_com_m : 1 :.mol, 4=.cfg ; 2=.xred ; 5 =.gin



    USE T_kind_param_m, ONLY:  double
#ifdef PARA
    USE mpi
    USE mod_para,only:MPI_COMM_space,status,ierr,nprocspace,myidsp,NDM_MPI_REAl_DOUBLE,rang
#else
    USE mod_para,only:myidsp
#endif
    ! ****************************************************************

    implicit none



    integer,intent(in),optional  :: itapp
    class(atom_config),intent(in)::atmol
    type(box_config),intent(in)::boxmol
    character*3,intent(in), dimension(1:atmol%im),optional  :: rty
    character(len=*), optional ::namefr
    logical,intent(in)::latcomp ! true= pas besoinde rapatrier atdml, false= il faut rapatrier atdml sur les masters
    logical, optional,intent(in):: lw0 ! seul le rang=0 écrit (implique latcomp=.true.)
    
    character*80::namef
    integer :: rgloc,im,imm,j,ic

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
    logical::latcompin=.false.
    type(para_config)::div
    logical :: lw0in=.false.

#endif
    type(atom_config)::atcomp
    integer :: i, luvisu, luvisu2, iti,lenfn2
    real(double) :: xp1, xp2, xp3,at(3,3),bg(3,3),pat
    character :: extension*9

#ifdef PARA
    write(6,*)'IN rasmol'
    latcompin=latcomp

    if (present (lw0))lw0in=lw0
    if (lw0in) then
       if(latcompin.eqv..false.) then
          write(6,*)'comment sauvegarder seulement rang0 si latcomp=.false. ?'
          stop
       end if
       if (rang==0) then
          latcompin=.true.
       else
          latcompin=.false. !latcompin intègre lw0 et rang=0
       end if
    end if

    
    if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.).and.(latcompin.eqv..false.)) then
       call atcomp%init(im_glob)
       div%rgim=myidsp
       div%npim=nprocspace
       div%comm_image=MPI_COMM_space
       call atmol%vers_master(atcomp,div)
       im =atcomp%im
       imm=atcomp%im
       rgloc=myidsp
    else
       rgloc=0
       call atcomp%init(atmol%im)
       call atmol%copy_config(atcomp,lrescl=.true.)
       im=atmol%im
       imm=atmol%imm
    end if
#else
    rgloc=0
    call atcomp%init(atmol%im,ltabvois=.false.,nvois=0)
    call atmol%copy_config(atcomp,lrescl=.true.)
    im=atmol%im
    imm=atmol%imm
#endif  


    at =boxmol%at*1d8 ; bg=boxmol%bg*1d-8
    allocate(xp(3,im));allocate(ityp(im));allocate(num_at_glob(imm))
    xp(1:3,1:atmol%im)=atcomp%xp(1:3,1:atmol%im)*1d8
    num_at_glob(1:atmol%im)=atcomp%num_at_glob(1:atmol%im)
    ityp(1:im)=atcomp%ityp(1:im)
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
    if(rgloc==0) then
       luvisu = 86
       luvisu2 = 87
#ifdef PARA   
       luvisu = 86+rang
       luvisu2 = 87+rang
#endif    

       ! TJ: change the formatting so that files are well listed.

       lenfn2 = 9
       if (present(itapp))then
          if (itapp < 0  )   extension='iiiiiiiii' 
          if (itapp >= 0 )   write(extension,'(i9.9)') itapp
       end if
       



       select case (ivisu)
       case(5)
          if (present(namefr))then
             if (present(itapp)) then
                namef=namefr(1:len(namefr))//'.'//extension(1:lenfn2)//'.newgin'
             else
                !                write(6,*)'BINGO'
                namef=namefr(1:len(namefr))//'.newgin'
             end if
          else
             if (present(itapp)) then
                namef=fnam(1:lenfnam)//'.'//extension(1:lenfn2)//'.newgin'
             else
                namef=fnam(1:lenfnam)//'.newgin'
             end if
          end if
          open(luvisu, file=namef, form='formatted', &
               status='unknown')
          write (luvisu,*)' 1 1 1 '
          write (luvisu,'(3F12.6)')at(1,1),at(2,1),at(3,1)
          write (luvisu,'(3F12.6)')at(1,2),at(2,2),at(3,2)
          write (luvisu,'(3F12.6)')at(1,3),at(2,3),at(3,3)
          write (luvisu,*) atcomp%im
          call cryst_to_cart (im, xp,  bg,  -1) !cart vers cryst          
       case (1)
          if (present(namefr))then
             if (present(itapp)) then
                namef=namefr(1:len(namefr))//'.'//extension(1:lenfn2)//'.mol'
             else
!                write(6,*)'BINGO'
                namef=namefr(1:len(namefr))//'.mol'
             end if
          else
             if (present(itapp)) then
                namef=fnam(1:lenfnam)//'.'//extension(1:lenfn2)//'.mol'
             else
                namef=fnam(1:lenfnam)//'.mol'
             end if
          end if

          open(luvisu, file=namef, form='formatted', &
               status='unknown')
          
          if (present(itapp))then
             write (luvisu, '(I9,A,I7,A,F12.6)') im, ' IT =', itapp, ' Time = ', timel
          else
             write (luvisu, '(I9,A,I7,A,F12.6)') im
          end if
          
          write (luvisu,'(9F12.6)')at(1,1),at(2,1),at(3,1),at(1,2),at(2,2),at(3,2),at(1,3),at(2,3),at(3,3)
          !at=at/1.d8

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
          write(luvisu,'(a,i0)')'Number of particles = ', im




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
          case(5)                    
             write (luvisu,'(3es15.6,I3)') xp1, xp2, xp3, ityp(i)
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




end module rasmolT_mod
