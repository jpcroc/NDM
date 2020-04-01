! ****************************************************************
subroutine rasmol(itapp)
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use var_pot
  use tab_imm_m
#if(PARA)
  use mod_mpi
#endif
  ! ****************************************************************

  implicit none
  !-----------------------------------------------
  !   G l o b a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  integer  :: itapp,iksp
#if(PARA)
  integer :: iproc
  real(double), allocatable :: xp_loc(:,:)
  integer, allocatable      :: ityp_loc(:)
  integer, allocatable      :: num_at_glob_loc(:)
  integer :: im_loc
  integer :: proc_source
#endif
  !-----------------------------------------------
  !   L o c a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: i, luvisu, luvisu2, iti,lenfn2
  real(double) :: xp1, xp2, xp3,pat
  character :: extension*9

  ! Notes about V_sim:
  ! * works if at(:,:) "encompasses" all the system (no duplication of lattice cells)
  ! * at(:,1) must be along x and at(:,2) must have no component along z.
  !   Otherwise a rotation matrix should be coded.
  !-----------------------------------------------
  !
  !
  !
  iksp=-1
  if (lcasca) iksp=iko
  if (ldesinteg)iksp=1
  luvisu = 86
  luvisu2 = 87
  !      write (*, *) 'entree dans rasmol.f',im

  ! ****ouverture fichier sortie pour traitement images rasmol****
  !if(rang==0)       OPEN(luvisu,file='donnrasmol',form='formatted', &
  !       status='unknown')
  !   la fonction char ne marche que pour de faibles valeurs de it!!!!
  !     open(luvisu,file='donn_rasmolit.'//char(48+it),form='formatted',status='unknown')

  ! conversion entier-->alphanumerique par transfert du nombre
  ! de l'iteration vers fichier tampon relu sous format caractere.

  if(rang==0) then

     ! TJ: change the formatting so that files are well listed.
     lenfn2 = 9
     if (itapp < 0  )   extension='iiiiiiiii' 
     if (itapp >= 0 )   write(extension,'(i9.9)') itapp
     ! -------------------------------------------------------------
     !     creation du  fichier positions pour le logiciel de visulation
     ! -------------------------------------------------------------
     select case (ivisu)
     case (1)
       open(luvisu, file=fnam(1:lenfnam)//'.'//extension(1:lenfn2)//'.mol', form='formatted', &
            status='unknown')
       write (luvisu, '(I9,A,I7,A,F12.6)') im_glob, ' IT =', itapp, ' Time = ', timel
       at=at*1.d8
       write (luvisu,'(9F12.6)')at(1,1),at(2,1),at(3,1),at(1,2),at(2,2),at(3,2),at(1,3),at(2,3),at(3,3)
       at=at/1.d8
     case (2)
       open(luvisu, file=fnam(1:lenfnam)//'.'//extension(1:lenfn2)//'.ascii', form='formatted', &
            status='unknown')
       write (luvisu, '(A,I9,A,I7,A,F12.6)') '#',im_glob, ' IT =', itapp, ' Time = ', timel
       if (abs(at(2,1)) > 1e-13 .or.abs(at(3,1)) >1e-13) then
         write(*,*) "Error, at(:,1) should only have a component along x with V_Sim"
         write(*,*) at(:,1)
       end if
       if (abs(at(3,2)) > 1e-13) then
         write(*,*) "Error, at(:,2) should have no component along z with V_Sim"
         write(*,*) at(:,2)
       end if
       at=at*1.d8
       write (luvisu,'(3F12.6)') at(1,1), at(1,2), at(2,2)
       write (luvisu,'(3F12.6)') at(1,3), at(2,3), at(3,3)
       at=at/1.d8
       open(luvisu2,file=fnam(1:lenfnam)//'.'//extension(1:lenfn2)//'.dat', form='formatted', &
            status='unknown')
       write(luvisu2,'(a)') '# Data file associated with an ascii file (V_Sim)'
    case(3) 
       open(luvisu, file=fnam(1:lenfnam)//'.'//extension(1:lenfn2)//'.xred', form='formatted', &
            status='unknown')
       at=at*1.d8
       write (luvisu,'(9F12.6)')at(1,1),at(2,1),at(3,1),at(1,2),at(2,2),at(3,2),at(1,3),at(2,3),at(3,3)
       at=at/1.d8
     end select
     
     
     !if(itapp==0) open (file='filmtot.mol',unit=47)
     !if(itapp==0) open (file='filmtot.mol',unit=47)
     !      write (47, 134) im
     !      write (47, *) 'IT =', itapp, '    Time = ', timel
  end if
   if (ivisu==3)  call cryst_to_cart (imm, xp,  bg,  -1) !cart vers cryst
#if(PARA)
  ! Le processeur maitre recoit les information des autres processeurs pour les ecrire sur fichier

  

  if (myid==0) then
     ! Copie des tableaux xp,num_at_glob et ityp locaux 
     allocate(xp_loc(3,imm))
     allocate(ityp_loc(imm))
     allocate(num_at_glob_loc(imm))
     xp_loc = xp
     ityp_loc = ityp
     num_at_glob_loc = num_at_glob
     im_loc = im

     ! Boucle sur les processeurs
     do iproc=0,nprocs-1
        ! Pour le processeur maitre il n'y a rien a faire
        ! reception des donnees des autres processeurs
        if (iproc.ne.0) then
           call MPI_RECV(im,               1,    MPI_INTEGER,      MPI_ANY_SOURCE, 10001, MPI_COMM_WORLD, status, ierr)
           proc_source = status(MPI_SOURCE)
           call MPI_RECV(xp(1:3,1:im),     3*im, NDM_MPI_REAL_DOUBLE, proc_source, 10002, MPI_COMM_WORLD, status, ierr)
           call MPI_RECV(ityp(1:im),       im,   MPI_INTEGER,         proc_source, 10003, MPI_COMM_WORLD, status, ierr)
           call MPI_RECV(num_at_glob(1:im),im,   MPI_INTEGER,         proc_source, 10004, MPI_COMM_WORLD, status, ierr)
        endif
#endif
        do i = 1, im
           xp1 = xp(1,i)*1D+08
           xp2 = xp(2,i)*1D+08
           xp3 = xp(3,i)*1D+08
           if (lPrtSigat) then
              if(it.eq.0)then
                 pat=0.0
              else
                 pat=unitP*(sigat(1,1,i)+sigat(2,2,i)+sigat(3,3,i))/3.
              end if
              select case (ivisu)
              case (1)
                write (luvisu, 136) ty(ityp(i)),xp1, xp2, xp3,pat,num_at_glob(i)
              case (2)                
                write (luvisu,'(3es15.6,2x,a)') xp1, xp2, xp3, ty(ityp(i))
                write (luvisu2,'(es15.6)') pat
              end select
           else if ((lprteat).and.(it.ne.0)) then
              select case (ivisu)
              case (1)
                write (luvisu, 136) ty(ityp(i)),xp1, xp2, xp3,eatom(i)*erg2eV ,num_at_glob(i)
              case (2)
                ! cubic case
!!$                xp1 = xp1 + 0.5
!!$                xp2 = xp2 + 0.5
!!$                xp3 = xp3 + 0.5
!!$                if (xp1 > at(1,1)*1e8) xp1 = xp1 - at(1,1)*1e8
!!$                if (xp2 > at(2,2)*1e8) xp2 = xp2 - at(2,2)*1e8
!!$                if (xp3 > at(3,3)*1e8) xp3 = xp3 - at(3,3)*1e8
                write (luvisu,'(3es15.6,2x,a)') xp1, xp2, xp3, ty(ityp(i))
                write (luvisu2,'(es15.6)') eatom(i)*erg2eV
              end select
           else  
              if (num_at_glob(i)==iksp) then
#if(PARA)
                 write (luvisu, 138) xp1, xp2, xp3,num_at_glob(i)
#else
                 write (luvisu, 138) xp1, xp2, xp3,i
#endif
               else
                 !write (luvisu, 135) ty(ityp(i)),xp1, xp2, xp3
                 select case (ivisu)
                 case (1)
#if(PARA)
                 write (luvisu, 135)  ty(ityp(i)),xp1, xp2, xp3,num_at_glob(i)
#else
                 write (luvisu, 135)  ty(ityp(i)),xp1, xp2, xp3,i
#endif

!                   write (luvisu, 135) ty(ityp(i)),xp1, xp2, xp3
                 case (2)
                   ! cubic case
!!$                   xp1 = xp1 + 0.5
!!$                   xp2 = xp2 + 0.5
!!$                   xp3 = xp3 + 0.5
!!$                   if (xp1 > at(1,1)*1e8) xp1 = xp1 - at(1,1)*1e8
!!$                   if (xp2 > at(2,2)*1e8) xp2 = xp2 - at(2,2)*1e8
!!$                   if (xp3 > at(3,3)*1e8) xp3 = xp3 - at(3,3)*1e8
                   write (luvisu,'(3es15.6,2x,a)') xp1, xp2, xp3, ty(ityp(i))
                case(3)                    
                   write (luvisu,'(3es15.6,2x,2a)') xp1/1d8, xp2/1d8, xp3/1d8, ' ! ', ty(ityp(i))
                end select
                 !               write (47, 135) ty(ityp(i)),xp1, xp2, xp3
              end if
           end if
        end do
#if(PARA)
     enddo  ! fin de boucle sur les processeurs
     ! Le processeur maitre recupere ses donnees locales
     xp = xp_loc
     ityp = ityp_loc
     num_at_glob = num_at_glob_loc
     im = im_loc
     deallocate(xp_loc)
     deallocate(ityp_loc)
     deallocate(num_at_glob_loc)

  else ! Les autres processeurs envoient leurs donnees locales
     call MPI_SEND(im,               1,   MPI_INTEGER,        0,10001,MPI_COMM_WORLD,ierr)
     call MPI_SEND(xp(1:3,1:im),     3*im,NDM_MPI_REAL_DOUBLE,0,10002,MPI_COMM_WORLD,ierr)
     call MPI_SEND(ityp(1:im),       im,  MPI_INTEGER,        0,10003,MPI_COMM_WORLD,ierr)
     call MPI_SEND(num_at_glob(1:im),im,  MPI_INTEGER,        0,10004,MPI_COMM_WORLD,ierr)
  endif
#endif
if (ivisu==3)    call cryst_to_cart (imm, xp , at,  1)  !cryst vers cart

  ! -------------------------------------------------------------
  !                          formats
  ! -------------------------------------------------------------
134 format(i6)
135 format(A,3f10.4,I9)
136 format(A,3f10.4,D14.5,I9)
138 format('Pb ',3f10.4,I9)

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

  if(rang==0)      close(luvisu)
  if ((rang == 0).and.(ivisu == 2)) then
    close(luvisu2)
  end if

  return
end subroutine rasmol


 subroutine redefine_ty
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use var_pot
  use tab_imm_m
#if(PARA)
  use mod_mpi
#endif
implicit none
ntyp_buffer=ntyp
allocate(ityp_buffer(imm),ty_buffer(ntyp),cm_buffer(ntyp))
 ityp_buffer(:)=ityp(:)
 ty_buffer(:)=ty(:)
 cm_buffer(:)=cm(:)
 
 ntyp=3
 ityp(1:6)=2
 ityp(8:15)=2
 ityp(7)=3
 
 
 deallocate(ty,cm)
 allocate(ty(ntyp),cm(ntyp))

 ty(1)='Fe' 
 ty(2)='Cu'
 ty(3)='O ' 

 cm(1:ntyp)=cm(1)

 return
 end subroutine redefine_ty


 subroutine refix_ty
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use var_pot
  use tab_imm_m
#if(PARA)
  use mod_mpi
#endif
implicit none 

 ntyp=ntyp_buffer
 deallocate(ty,ityp)
 allocate(ty(ntyp),cm(ntyp),ityp(imm))
 ityp(:)=ityp_buffer(:)
 ty(:)=ty_buffer(:)
 cm(:)=cm_buffer(:)

 return
 end subroutine refix_ty







