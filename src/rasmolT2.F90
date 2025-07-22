! ****************************************************************
module rasmolT2_mod
  USE arret_ndm_mod,only:arret_ndm
  USE cryst_to_cart_mod,only: cryst_to_cart
  USE gen_com_m, ONLY:rang,ivisu,lpkbar,lspaceNDM,&
       &cunitP,iteration,lcasca,timel,unitP,fnam,erg2ev,lenfnam,umass,rang
  USE var_pot, ONLY:ntyp,ntyp_buffer,ty,ty_buffer,cm_buffer,cm

  use atomconfig,only: atom_config,atom_config_d,atom_config_e,atom_config_arps
  use paraconfig,only:para_config
  use boxconfig,only:box_config
  implicit none
  integer, dimension(:), allocatable       :: ityp_buffer   ! temp/iorary store the types buffer when 
contains

  subroutine rasmolT2(atmol,boxmol,itapp,namefr,rty,latcomp,ivisumol,naux,charaux,vaux,lappend)
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
    USE T_kind_param_m, ONLY:  double
#ifdef PARA
    USE Tpara,only:COMM_space,nprocspace,myidsp,nprocs
#else
    USE Tpara,only:myidsp,nprocs
#endif

    implicit none
    integer,intent(in),optional  :: itapp
    class(atom_config),intent(in)::atmol
    class(box_config),intent(in)::boxmol
    character*3,intent(in), dimension(1:atmol%im),optional,target  :: rty
    character(len=*), optional ::namefr
    logical::latcomp ! true= pas besoinde rapatrier atdml, false= il faut rapatrier atdml sur les masters
    integer, optional:: ivisumol
    logical,optional::lappend

    logical::lappendF

    integer,optional::naux
    character(len=*),optional::charaux(:)
    real(double),optional::vaux(:,:) 
    logical::laux
    integer::nauxV,nauxtot
    integer::ivisum
    character*80::nameo,end_name
    integer :: rgloc,j,ic

    character*3, dimension(:), pointer  :: tyw

#ifdef PARA
    type(para_config)::div

#endif
    class(atom_config),allocatable::atcomp
    integer :: i, luvisu, luvisu2,lenfn2
    real(double) :: xp1, xp2, xp3,at(3,3),bg(3,3),pat
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
    nauxv=0
    if (present(naux)) then
       if (naux.gt.0) then
          laux=.true.
          nauxV=naux
       end if
    end if
    nauxtot=nauxv
    select type (atmol)
    class is (atom_config_e)
       if (atmol%lsigat) then
          nauxtot=nauxtot+1
       end if
       if (atmol%lprteat) nauxtot=nauxtot+1
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
    call atmol%deftype(atcomp)
    if (latc.eqv..false.) then
       if (laux) then
          write(6,*)'laux TRUE et latc FAUX  stop (FLEMME)'
          call arret_ndm
       end if
       if (present(rty)) then
          write(6,*)'RTY TRUE and latc FAUX  stop (FLEMME) use a type extension'
          call arret_ndm
       end if
       if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
          call atcomp%init(atmol%im_glob,im_glob=atmol%im_glob,imm_glob=atmol%imm_glob)
!          call atcomp%print
          div%mpi_image%rank=myidsp
          div%mpi_image%nproc=nprocspace
          div%mpi_image%comm=COMM_space%comm
         ! write(6,*)
          !          call atmol%vers_master(atcomp,div,'ixnlusv')
          call atmol%vers_master(atcomp,div)
          !im =atcomp%im
          !imm=atcomp%im
          rgloc=myidsp
       else
          write(6,*)'latc=false et (nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) ??? stop'
          write(6,*)latc,nprocspace,lspaceNDM
          call arret_ndm
       end if
    else
       if (rgloc==0) then
          call atcomp%init(atmol%im)
          call atmol%copy_config(atcomp,lrescl=.true.)
          !im=atmol%im
          !imm=atmol%imm
       end if
    end if
#else
    call atmol%deftype(atcomp)
    rgloc=0
    call atcomp%init(atmol%im,ltabvois=.false.,nvois=0)
    call atmol%copy_config(atcomp,lrescl=.true.)
!    im=atmol%im
    !imm=atmol%imm
#endif  
    if(rgloc==0) then
!*****************PPPPPPPPPAAAAAAAAASSSSSSSSAAAAAAAAAAAGGGGGGGGEEEEEEEEEE en AngSTROMS!!!!!!!!!!!!!!!!
       at =boxmol%at*1d8 ; bg=boxmol%bg*1d-8
       atcomp%xp(1:3,1:atcomp%im)=atcomp%xp(1:3,1:atcomp%im)*1d8
       if (ivisum.ne.5) then
!          tyw='000'
          !    do i=1,im
          !       write(6,*)i,atmol%ityp(i),ty(atmol%ityp(i))
          !    end do
          if (present (rty))then
             tyw=>rty
          else
             allocate(tyw(atcomp%im))
             tyw(1:atcomp%im)=ty(atcomp%ityp(1:atcomp%im))
          end if
       end if
       select type (atcomp)
       class is (atom_config_arps)
          do i=1,atcomp%im
             select case (atcomp%mov(i))
             case(0)
                tyw(i)=' Re'
             case(1)
                tyw(i)=' In'             
             case(2)
                tyw(i)=' Mo'
             end select
          end do
       end select


!!$       select type (atmol)
!!$       type is (atom_config_e)
!!$          if(atmol%lsigat) then
!!$             allocate (sigat(3,3,imm))
!!$             sigat=0
!!$             sigat(:,:,1:im)=atmol%sigat(:,:,:1:im)
!!$          end if
!!$
!!$          if(atmol%lprteat) then
!!$
!!$             allocate (eat(im)) ; eat=0; eat(1:im)=atmol%eat(1:im)
!!$          end if
!!$       end select

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
       case default
          write(6,*)'wrong ivisu',ivisum,ivisu
          call arret_ndm
       end select
       
       if (present(itapp))then
          call openfilemol( luvisu,nameo,end_name,extension)
       else
          call openfilemol( luvisu,nameo,end_name)
       end if

       select case (ivisum)
       case(7) !xyz
          write (luvisu,*)atcomp%im 
          write (luvisu,'(A)',advance='no')nameo
          if (laux) then
             if (present(charaux)) then
                do iax=1,naux
                   write(luvisu,'(A)',advance='no')trim(charaux(iax))
                end do
             end if
          end if
          write(luvisu,*)' '
         do i = 1, atcomp%im
             xp1 = atcomp%xp(1,i)
             xp2 = atcomp%xp(2,i)
             xp3 = atcomp%xp(3,i)
             if (laux) then
                write (luvisu,'(A,3es15.6)',advance='no')ty(atcomp%ityp(i)), xp1, xp2, xp3
                do iax=1,naux
                   write(luvisu,'(G20.12)',advance='no')vaux(iax,i)
                end do
                write(luvisu,*)' '
             else
                write (luvisu,'(A,3es15.6,I9)')ty(atcomp%ityp(i)), xp1, xp2, xp3
             end if
          end do
          write(luvisu,*)
          write (luvisu,'(3F15.9)')at(1,1),at(2,1),at(3,1)
          write (luvisu,'(3F15.9)')at(1,2),at(2,2),at(3,2)
          write (luvisu,'(3F15.9)')at(1,3),at(2,3),at(3,3)
          write (luvisu,*) 
       case(5) !newgin
          write (luvisu,'(A)',advance='no')' 1 1 1 !'
          if (laux) then
             if (present(charaux)) then
                do iax=1,naux
                   write(luvisu,'(A)',advance='no')trim(charaux(iax))
                end do
             end if
          end if
          write(luvisu,*)' '
          write (luvisu,'(3F15.9)')at(1,1),at(2,1),at(3,1)
          write (luvisu,'(3F15.9)')at(1,2),at(2,2),at(3,2)
          write (luvisu,'(3F15.9)')at(1,3),at(2,3),at(3,3)
          write (luvisu,*) atcomp%im
          call cryst_to_cart (atcomp%im, atcomp%xp,  bg,  -1) !cart vers cryst
          do i = 1, atcomp%im
             xp1 = atcomp%xp(1,i)
             xp2 = atcomp%xp(2,i)
             xp3 = atcomp%xp(3,i)
             write (luvisu,'(3es15.6,I3, I9)',advance='no') xp1, xp2, xp3, atcomp%ityp(i),atcomp%num_at_glob(i)
!!$             if (laux) then
!!$                do iax=1,naux
!!$                   write(luvisu,'(G20.12)',advance='no')vaux(iax,i)
!!$                end do
!!$             end if
             select type (atcomp)
             class is (atom_config_e)
                if (atcomp%lsigat) then
                   if(iteration.eq.0)then
                      pat=0.0
                   else
                      pat=unitP*(atcomp%sigat(1,1,i)+atcomp%sigat(2,2,i)+atcomp%sigat(3,3,i))/3.
                   end if
                   write (luvisu, '(G20.12)',advance='no') pat*1d-9
                end if
                if (atcomp%lprteat) write (luvisu, '(G20.12)',advance='no') atcomp%eat(i)*erg2ev
             end select
             if (laux) then
                do iax=1,naux
                   write(luvisu,'(G20.12)',advance='no')vaux(iax,i)
                end do
             end if

             write(luvisu,*)' '
!             else
!                write (luvisu,'(3es15.6,I3)') xp1, xp2, xp3, atcomp%ityp(i)
!             end if
          end do
       case (1) !mol
          
          if (present(itapp))then
             write (luvisu, '(I9,A,I12,A,F12.6)',advance='no') atcomp%im, ' IT =', itapp, ' Time = ', timel
          else
             write (luvisu, '(I9,A,I7,A,F12.6)',advance='no') atcomp%im
          end if
          if (laux) then
             if (present(charaux)) then
                do iax=1,naux
                   write(luvisu,'(A)',advance='no')trim(charaux(iax))
                end do
             end if
          end if
          write(luvisu,*)' '

          write (luvisu,'(9F12.6)')at(1,1),at(2,1),at(3,1),at(1,2),at(2,2),at(3,2),at(1,3),at(2,3),at(3,3)
          !at=at/1.d8
          do i = 1, atcomp%im
             xp1 = atcomp%xp(1,i)
             xp2 = atcomp%xp(2,i)
             xp3 = atcomp%xp(3,i)
             !                write (6,*) 't',tyw(i)
             !                write(6,*)'x', xp1,xp2, xp3
             write (luvisu, '(A,3G18.4)',advance='no') tyw(i),xp1, xp2, xp3
#ifdef PARA
             write (luvisu, '(I9)',advance='no')  atcomp%num_at_glob(i)
#else
             write (luvisu, '(I9)',advance='no')  i
#endif
             if (laux) then
                do iax=1,naux
                   write(luvisu,'(G20.12)',advance='no')vaux(iax,i)
                end do
             end if
             select type (atcomp)
             class is (atom_config_e)
                if (atcomp%lsigat) then
                   if(iteration.eq.0)then
                      pat=0.0
                   else
                      pat=unitP*(atcomp%sigat(1,1,i)+atcomp%sigat(2,2,i)+atcomp%sigat(3,3,i))/3.
                   end if
                   write (luvisu, '(D14.5)',advance='no') pat*1d-9
                end if
                if (atcomp%lprteat) write (luvisu, '(D14.5)',advance='no') atcomp%eat(i)*erg2ev
             end select
             write(luvisu,*)' '


          end do
          !       case (2)
       case(3) 
          write (luvisu,'(9F12.6)')at(1,1),at(2,1),at(3,1),at(1,2),at(2,2),at(3,2),at(1,3),at(2,3),at(3,3)
          call cryst_to_cart (atcomp%im, atcomp%xp,  bg,  -1) !cart vers cryst
          do i = 1, atcomp%im
             xp1 = atcomp%xp(1,i)
             xp2 = atcomp%xp(2,i)
             xp3 = atcomp%xp(3,i)
             write (luvisu,'(3es15.6,2x,2a)') xp1, xp2, xp3, ' ! ', tyw(i)
          end do

       case(40,41,60,61) 

          write(luvisu,'(a,i0)')'Number of particles = ',  atcomp%im
          write(luvisu,'(a)')'A = 1.000 Angstrom (basic length-scale)'
          do j=1,3
             write(luvisu,'(a,i0)') '# Unit cell vector #', j    
             do ic=1,3
                write(luvisu,'(A,I1,A,I1,A,g16.8,A)')'H0(',j,',',ic,') = ',at(ic,j),' A'
             end do
          end do
          if ((ivisum==40).or.(ivisum==60)) then
             write(luvisu,'(A)')'.NO_VELOCITY.'
             write(luvisu,'(A,I0)')'entry_count = ', 3+nauxV

          else
             write(luvisu,'(A,I0)')'entry_count = ', 6+nauxV
             if ((laux).and.(present(charaux))) then
                do iax=1,naux
                   write(luvisu,'(A,I0,A,A)')'auxiliary[',iax,'] = ',trim(charaux(iax))
                end do
             end if
          end if
!          write(6,*)'nuaxv nauxtot',nauxv,nauxtot
          iax=nauxV
          select type (atcomp)
          class is (atom_config_e)
             if (atcomp%lsigat) then
                iax=iax+1
                write(luvisu,'(A,I0,A)')'auxiliary[',iax,'] = at_pressure'
             end if
             if (atcomp%lprteat) then
                iax=iax+1
                write(luvisu,'(A,I0,A)')'auxiliary[',iax,'] = at_energy'
             end if
          end select
          call cryst_to_cart (atcomp%im, atcomp%xp,  bg,  -1) !cart vers cryst
          do i=1, atcomp%im
             select case(ivisum)
             case(40,41)
                WRITE(luvisu,'(f0.3)') cm(atcomp%ityp(i))/umass        ! Mass (g/mol)
                WRITE(luvisu,'(a,A,I9)') tyw(i),' ! ',atcomp%num_at_glob(i)                 ! Atom type
                write(luvisu, '(3(g20.12,1x))',advance='no') atcomp%xp(:,i)
                if (ivisum==41) then
                   select type (atcomp)
                   class is (atom_config_d)
                      write (luvisu, '(3(g20.12,1x))',advance='no') atcomp%vp(:,i)*1d8*1d-12
                   end select
                end if
                if (laux) then
                   do iax=1,naux
                      write(luvisu,'(G20.12)',advance='no')vaux(iax,i)
                   end do
                end if
                select type (atcomp)
                class is (atom_config_e)
                   if (atcomp%lsigat) then
                      if(iteration.eq.0)then
                         pat=0.0
                      else
                         pat=unitP*(atcomp%sigat(1,1,i)+atcomp%sigat(2,2,i)+atcomp%sigat(3,3,i))/3.
                      end if
                      write (luvisu, '(D14.5)',advance='no') pat*1d-9
                   end if
                   if (atcomp%lprteat) write (luvisu, '(D14.5)',advance='no') atcomp%eat(i)*erg2ev
                end select

                write(luvisu,*)' '

             case(60,61)
                if (i==1) then
                   WRITE(luvisu,'(f0.3)') cm(atcomp%ityp(i))/umass        ! Mass (g/mol)
                WRITE(luvisu,'(a,A,I9)') tyw(i),' ! ',atcomp%num_at_glob(i)
                else
                   if (atcomp%ityp(i).ne.atcomp%ityp(i-1)) then
                      WRITE(luvisu,'(f0.3)') cm(atcomp%ityp(i))/umass        ! Mass (g/mol)
                      WRITE(luvisu,'(a,A,I9)') tyw(i),' ! ',atcomp%num_at_glob(i)
                   end if
                end if

                write(luvisu, '(3(g20.12,1x))',advance='no') atcomp%xp(:,i)
                if (ivisum==61) then
                   select type (atcomp)
                   class is (atom_config_d)
                      write (luvisu, '(3(g20.12,1x))',advance='no') atcomp%vp(:,i)*1d8*1d-12
                   end select
                end if
                select type (atcomp)
                class is (atom_config_e)
                   if (atcomp%lsigat) then
                      if(iteration.eq.0)then
                         pat=0.0
                      else
                         pat=unitP*(atcomp%sigat(1,1,i)+atcomp%sigat(2,2,i)+atcomp%sigat(3,3,i))/3.
                      end if
                      write (luvisu, '(D14.5)',advance='no') pat*1d-9
                   end if
                   if (atcomp%lprteat) write (luvisu, '(D14.5)',advance='no') atcomp%eat(i)*erg2ev
                end select

                if (laux) then
                   do iax=1,naux
                      write(luvisu,'(G20.12)',advance='no')vaux(iax,i)
                   end do
                end if
                write(luvisu,*)' '
             end select
          end do
       end select
       if (.not.lappendF) close(luvisu)

    end if
!    write(6,*)'OUT Rasmol',rang

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
  end subroutine rasmolT2

  subroutine openfilemol2(luvisu,nameo,end_name,ext)
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
  end subroutine openfilemol2
end module  rasmolT2_mod
