! This is the main program for PHONDY version june 2017
!
!
!/ PHONDY method, Copyright Mihai-Cosmin Marinica, July 2010



subroutine phondy
!-----------------------------------------------

      USE T_kind_param_m, ONLY:  double
      use gen_com_m
      use var_pot
      use tab_imm_m

!#if (ML)
!   use ml_in_ndm_module, only: md_iconf
!   use derived_types, only: config_real
!#endif


#if(PARAPH)
! use mkl_service
 use mpi
 use mod_mpi_phondy
#if(LAMMPS_VERSION)
 use LAMMPS
 use vars_lammps
#endif
#endif
      use phondy_in_ndm_module
      use derived_types_ph, only: config_real
      implicit none
#if(PARAPH)
    !  integer(4) :: status
    !  integer :: numthreads
      integer :: Nin
#endif

  call init_phondy()
  call allocate_phondy()
  call read_phondy_file()


  if (rangph==0) write(6,*)
  if (rangph==0) write(6,*)
  if (rangph==0) write(6,*)'************ DEBUT DE PHONDY ****************'
  if (rangph==0) write(6,*)
#if(LAMMPS_VERSION)
  if (rangph==0) write(6,*)'---------------LAMMPS FORCES-----------------'
#endif

#if(ML || PARAML)
  if (rangph==0) write(6,*)'-----------------ML FORCES-------------------'
#endif

  if (rangph==0) write(6,*)


  call print_phondy(rangph)



!#if (ML)
  if (lq_points) then
   call md_init_config_ph
   ! put the initial gin config into config_real(1) object
   call put_ndm_into_ph_config(iconf_ini)
   ! Test the size of the  config_real(1) object
   call test_if_config_is_small(iconf_ini)
   ! Replicate the config_real(1) object into config_real(2) object
   call ph_build_box
   ! Compute the index for phase
   call allocate_object_ph
   !>begin ... neighbours subroutines for the big box ....
   !prepare  the neighbours for iconf_big
   call neighbours_ndm_layer(iconf_big)
   ! put evrything in ndm
   call put_ph_config_into_ndm(iconf_big)
   ! compute the neighbours ...
   call calc_neighbours_ph(iconf_big)
   !>end  ... neighbours
   call allocate_phondy
   nmat = 3*config_real(iconf_ini)%nat
#ifdef LAMMPS_VERSION
   call dump_lammps(config_real(iconf_big)%imm, config_real(iconf_big)%im, config_real(iconf_big)%itype, config_real(iconf_big)%cell, config_real(iconf_big)%pos_cart)
   call MPI_BARRIER(MPI_COMM_WORLD,codeph)

   if (rangph==0) write(6,*) 'init communicators for lammps...'
   call define_communicators()
!v1   call MPI_BARRIER(orig_group,codeph)
   call init_potential_lammps()
!v1   call MPI_BARRIER(orig_group,codeph)
   if (rangph==0) write(6,*) '...init LAMMPS forces'
#endif

  end if

#ifdef LAMMPS_VERSION
if (.not.lq_points) then
   if (rangph==0) write(6,*) 'init communicators for lammps...'
   call define_communicators()
!v1   call MPI_BARRIER(orig_group,codeph)
   call init_potential_lammps()
!v1   call MPI_BARRIER(orig_group,codeph)
   call dump_lammps(imm, im, ityp(:), at(:,:)/A2cm, xp(:,:)/A2cm)
   if (rangph==0) write(6,*) '...init LAMMPS forces'
end if
#endif

!#endif
  call test_minimum
  call force_constant
!end if

  ! if this is uncomment everything blows up.
  !call lammps_close(lmp)


#if PARAPH && PHONDY


#if(LAMMPS_VERSION)
!v1 call MPI_BARRIER(orig_group,codeph)
!v1 all MPI_COMM_FREE (lammps_comm(w_rang), codeph)
call MPI_BARRIER(MPI_COMM_WORLD,codeph)
call MPI_COMM_SIZE(MPI_COMM_WORLD,w_size,codeph)
call MPI_COMM_RANK(MPI_COMM_WORLD,rangph,codeph)
#endif
 if (isave==0) then
  Nin=nmat
  call diago_scalapack (Nin,W)
 end if
 !
 if (isave==2) then
  if (rangph==0 ) then

    if (lq_points) then
      call diago_q ()
    else
      call diago_threading_real ()
    end if

  end if
!debug write(6,*) 'rangph 1', rangph
end if

!PARAPH && PHONDY
#else
 if (isave==0) then
  call diago_threading_real ()
 end if
#endif


#if(PARAPH && PHONDY)
! if I uncomment that evertything blows up ... all the rangph becomes 0 !!!!
!  call MPI_BARRIER(MPI_COMM_WORLD, rangph)
#endif

!write(6,*) 'rangph 2', rangph

 if ((isave==0).or.(isave==2)) then
  call write_dos_ldos ()
 end if

  if (isave==1) then

   if (rangph==0) write(6,*) 'PHONDY: the DM matrix have been saved for the subsequent diagonalization '

  end if

!write(6,*) 'rangph 2', rangph

  if (rangph==0) write(6,*)
  if (rangph==0) write(6,*)
  if (rangph==0) write(6,*)'************  FIN  DE PHONDY ****************'
  if (rangph==0) write(6,*)
  if (rangph==0) write(6,*)

#if(PARAPH && PHONDY)
call MPI_FINALIZE(codeph)
#endif

  return


  end subroutine phondy
