module gcII_mod
  USE endrun_mod,only: endrun
  USE analyse_mod,only: analyse
  USE initspeed_mod,only: initspeed,bruit_xp
  USE gcmodII_mod,only: ZXCGRII
  USE tab_imm_m,only:xp,ityp,bruitmd,num_at_glob
  implicit none

contains
  ! *************************************************************
  subroutine  gcII   ! (xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double
    USE gen_com_m, ONLY:bg,itetemp2,imm_glob,dmtype,rang,im,it,itmax,mdcg_noise,&
         &angst,erg2ev,imm,potist
    USE var_pot, ONLY:nad,na,ntyp
    USE work_cgII,only: funct
#ifdef PARA

use mod_para
#endif
    ! *************************************************************
    ! xp positions des atomes
    ! xpp previous positions
    ! vp  velocities
    ! ax starting positions
    ! fp forces
    !iwmax =  indice du dernier voisin de chaque atome
    !ityp  tableau des types
    implicit none
    !-----------------------------------------------
    !   G l o b a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    !integer  :: ielat(imm)
    !integer  :: iwmax(imm)
    !integer  :: ityp(imm)
    !real(double)  :: xp(3,imm)
    !real(double)  :: xpp(3,imm)
    !real(double)  :: vp(3,imm)
    !real(double)  :: ax(3,imm)
    !real(double)  :: fp(3,imm)
    !-----------------------------------------------
    !   L o c a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: n,  i, igc,ims
    real(double) :: efinal
    !GC settings
    integer  :: NGC,criterion,NCALLS,IER
    double precision, dimension(:), allocatable:: X,G,W
    double precision     :: ACC,F
    character :: extension*2
    integer::lenfn2
    !-----------------------------------------------
    !
    !
    real(double), allocatable :: xp_all(:,:),fp_all(:,:)
    integer, allocatable      :: ityp_all(:)

#ifdef PARA
    integer :: iproc
    integer, allocatable      :: num_at_glob_all(:)
    !  integer :: im_loc
    integer :: proc_source
    integer::nag(imm)
#endif


!    stop

    nad(:ntyp) = na(:ntyp)

    if (mdcg_noise /= 0 ) then
       call bruit_xp
    end if



    
#ifdef PARA
!    write(6,*)'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT,rang'
!    write(6,*)rang,im,imm,im_glob,imm_glob


    NGC=3*imm_glob
    allocate (X(NGC),G(NGC),W(6*NGC))
    allocate(xp_all(3,imm_glob))
    allocate(fp_all(3,imm_glob))
    allocate(ityp_all(imm_glob))
    xp_all=0
    if (rang==0) then
       do iproc=0,nprocs-1
          ! Pour le processeur maitre il n'y a rien a faire
          ! reception des donnees des autres processeurs
          if (iproc.ne.0) then
             call MPI_RECV(im,               1,    MPI_INTEGER,      MPI_ANY_SOURCE, 10001, MPI_COMM_WORLD, status, ierr)
             proc_source = status(MPI_SOURCE)
             call MPI_RECV(xp(1:3,1:im),     3*im, NDM_MPI_REAL_DOUBLE, proc_source, 10002, MPI_COMM_WORLD, status, ierr)
             call MPI_RECV(ityp(1:im),       im,   MPI_INTEGER,         proc_source, 10003, MPI_COMM_WORLD, status, ierr)
             call MPI_RECV(nag(1:im),im,   MPI_INTEGER,         proc_source, 10004, MPI_COMM_WORLD, status, ierr)
          else
             nag(1:im)=num_at_glob(1:im)
          end if
!          write(6,*)'proc im',iproc,im,nag(1:10)
          do i=1,im
             IF (dmtype.EQ.30) THEN
                ! Variables = reduced coordinates
                if (mdcg_noise==0) then

                   X(3*nag(i)-2:3*nag(i)) = MatMul( xp(1:3,i), bg)
                else
                   xp(1:3,i)= xp(1:3,i)+bruitmd(1:3,i)
                   X(3*nag(i)-2:3*nag(i)) = MatMul( xp(1:3,i), bg)
                end if
             ELSE
                ! Variables = cartesian coordinates (in A)
                if (mdcg_noise==0) then
                   X(3*nag(i)-2:3*nag(i)) = xp(1:3,i)*angst
                else
                   xp(1:3,i)= xp(1:3,i)+bruitmd(1:3,i)
                  X(3*nag(i)-2:3*nag(i)) = (xp(1:3,i)+bruitmd(1:3,i))*angst
                end if
             END IF
             if (any(xp_all(:,nag(i)).ne.0)) then
                write(6,*)'iproc,i,nag'
                write(6,*)iproc,i,nag(i)
             end if
             xp_all(:,nag(i))=xp(:,i)
             ityp_all(nag(i))=ityp(i)
          end do
       end do

    else
       call MPI_SEND(im,               1,   MPI_INTEGER,        0,10001,MPI_COMM_WORLD,ierr)
       call MPI_SEND(xp(1:3,1:im),     3*im,NDM_MPI_REAL_DOUBLE,0,10002,MPI_COMM_WORLD,ierr)
       call MPI_SEND(ityp(1:im),       im,  MPI_INTEGER,        0,10003,MPI_COMM_WORLD,ierr)
       call MPI_SEND(num_at_glob(1:im),im,  MPI_INTEGER,        0,10004,MPI_COMM_WORLD,ierr)
       xp_all=0; fp_all=0;ityp_all=0
       X=0
    end if

    X(3*im_glob+1:3*imm_glob)=0.d0
    ims=imm_glob
!    write(6,*)rang, 'POST REPART',im_glob
!    write(extension,'(i2.2)') rang
!    lenfn2=2
!    open(unit=616, file='posall.'//extension(1:lenfn2)//'.csv', form='formatted', &
!             status='unknown')
!        do i=1,im_glob
!       write(6,*)rang,i,xp_all(:,i)
!       write(616,'(I2,I6,3G18.9)') rang,i,xp_all(:,i)
!    end do
!    call mpi_barrier(MPI_COMM_WORLD,ierr)
!    stop
#else
    allocate(xp_all(3,imm))
    allocate(fp_all(3,imm))
    allocate(ityp_all(imm))

!    write(6,*)'preGC0A',xp(:,1),xp_all(:,1),it,'dmtype ',dmtype
    !New GC settings ....:
    NGC=3*imm

    allocate (X(NGC),G(NGC),W(6*NGC))
    do i=1,im
       IF (dmtype.EQ.30) THEN
          ! Variables = reduced coordinates
          if (mdcg_noise==0) then

             X(3*i-2:3*i) = MatMul( xp(1:3,i), bg)
          else
             xp(1:3,i)= xp(1:3,i)+bruitmd(1:3,i)
             X(3*i-2:3*i) = MatMul( xp(1:3,i), bg)
          end if
       ELSE
          ! Variables = cartesian coordinates (in A)
          if (mdcg_noise==0) then
             X(3*i-2:3*i) = xp(1:3,i)*angst
          else
             xp(1:3,i)= xp(1:3,i)+bruitmd(1:3,i)
             X(3*i-2:3*i) = (xp(1:3,i)+bruitmd(1:3,i))*angst
          end if
       END IF
          xp_all(:,i)=xp(:,i)
          ityp_all(i)=ityp(i)

    end do
!    write(6,*)'preGC0B',xp(:,1),xp_all(:,1),it,'dmtype ',dmtype
    X(3*im+1:3*imm)=0.d0
    ims=imm
#endif

    !
    ! internal units
    !      epsilon_force=epsilon_force/6.251d3
    criterion=0
!!$ACC=1.d-8
    ACC=0.d0


    it = 0
    efinal = 0.0
!    write(6,*)'preGC',xp_all(:,1),it
    CALL ZXCGRII(FUNCT,NGC,ACC,itmax,X,G,F,W,IER,criterion,NCALLS, &
         xp_all,fp_all, ityp_all,ims)       

    deallocate (X,G,W)
    deallocate(xp_all,ityp_all)

    if (rang==0)write(*,*) 'DEBUG after GCII ... just SAY HALLO'


#ifdef ART
    write (*, *) 'energie ', potist,'erg',potist*erg2eV,'eV'
    return

#else
    it=0
    itetemp2=0
    if (rang==0)     write (*, *) 'energie ',potist*erg2eV,'   eV', potist,'erg'



    call analyse
    call endrun 
    return
#endif


  endsubroutine gcII

end module gcII_mod
