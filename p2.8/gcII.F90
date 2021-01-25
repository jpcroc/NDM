module gcII_mod
  USE endrun_mod,only: endrun
  USE analyse_mod,only: analyse
  USE initspeed_mod,only: initspeed,bruit_xp
  USE gcmodII_mod,only: ZXCGRII
  USE work_cgII,only: atcgcomp,cellcgcomp,atcgloc,cellcgloc,cellcible,atcible,boxcg,gcpara
  USE atomconfig,only:atom_config
  USE cellconfig,only:cell_config
  USE boxconfig,only:box_config
  USE endrunT_mod,only:endrunT
  implicit none

  
contains
  ! *************************************************************
  subroutine  gcII  (atcgin,celcgin,boxndm)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double
    USE gen_com_m, ONLY:itetemp2,imm_glob,dmtype,rang,it,itmax,mdcg_noise,&
         &angst,erg2ev,potist,im_glob,lperiod,lspacendm,latcomp
    USE var_pot, ONLY:ntyp
    USE work_cgII,only: funct
#ifdef PARA
    use paraconfig,only:para_config,initparapuresp
    USE parautils,only:initcomp
    use mod_para,only:nprocspace,MPI_COMM_space ,ierr,status,NDM_MPI_REAL_DOUBLE
#else
    use mod_para,only:nprocspace
#endif
#ifdef PARA


    include "mpif.h" 
#endif

    type(atom_config),target::atcgin
    type(cell_config),target::celcgin
    type(box_config)::boxndm

    integer :: n,  i, igc
    real(double) :: efinal
    !GC settings
    integer  :: NGC,criterion,NCALLS,IER
    double precision, dimension(:), allocatable:: X,G,W
    double precision     :: ACC,F
    character :: extension*2
    integer::lenfn2,i1
    !-----------------------------------------------
    !
    !
    real(double), allocatable :: bruitmd(:,:)
    integer, allocatable      :: ityp_all(:)

#ifdef PARA
    integer :: iproc
    integer, allocatable      :: num_at_glob_all(:)
    !  integer :: im_loc
    integer :: proc_source


#endif
!    PROBLEME AVEC NUM_AT_GLOB DANS INITCOMP (VERS_MASTER ?)
    cellcgloc=>celcgin
    atcgloc=>atcgin
    boxcg=boxndm

    !    stop



#ifdef PARA
if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
    call atcgcomp%init(im_glob,imm_glob)
    call initparapuresp(gcpara,rang,mpi_comm_space,nprocspace)
    call initcomp(atcgcomp,cellcgcomp,atcgin,celcgin,boxcg,gcpara,lperiod)
 else
    call initparapuresp(gcpara,rang,mpi_comm_space,nprocspace)
    atcgcomp=atcgin
    cellcgcomp=celcgin
 end if

#else
    atcgcomp=atcgin
    cellcgcomp=celcgin
#endif
    if (mdcg_noise /= 0 ) then
       call bruit_xp (atcgcomp%xp,bruitmd,atcgcomp%im)
    end if

       NGC=3*atcgcomp%im
       allocate (X(NGC),G(NGC),W(6*NGC))
       X=0;G=0;W=0
       if (gcpara%rgim==0)  then   
          do i=1,atcgcomp%im
             i1=atcgcomp%num_at_glob(i)
             
             
             IF (dmtype.EQ.30) THEN
                ! Variables = reduced coordinates
                if (mdcg_noise==0) then
                   
                   X(3*i1-2:3*i1) = MatMul( atcgcomp%xp(1:3,i), boxcg%bg)
                else
                   atcgcomp%xp(1:3,i)= atcgcomp%xp(1:3,i)+bruitmd(1:3,i)
                   X(3*i1-2:3*i1) = MatMul( atcgcomp%xp(1:3,i), boxcg%bg)
                end if
             ELSE
                ! Variables = cartesian coordinates (in A)
                if (mdcg_noise==0) then
                   X(3*i1-2:3*i1) = atcgcomp%xp(1:3,i)*angst
                else
                   atcgcomp%xp(1:3,i)= atcgcomp%xp(1:3,i)+bruitmd(1:3,i)
                X(3*i1-2:3*i1) = (atcgcomp%xp(1:3,i)+bruitmd(1:3,i))*angst
             end if
          END IF
       end do
    end if
       criterion=0
!!$ACC=1.d-8
       ACC=0.d0
       
       it = 0
       efinal = 0.0
       CALL ZXCGRII(FUNCT,NGC,ACC,itmax,X,G,F,W,IER,criterion,NCALLS)       
       
       deallocate (X,G,W)

    if (rang==0)write(*,*) 'DEBUG after GCII ... just SAY HALLO'

    latcomp=.true.
    call endrunT(atcgcomp,cellcgcomp,boxcg,latcomp) 
    return


  endsubroutine gcII

end module gcII_mod
