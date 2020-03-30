MODULE vars_lammps
  use LAMMPS
  type (C_ptr) :: lmp
  !old not workingversion
  integer :: comm_lammps, orig_group
  integer, dimension(:), allocatable :: lammps_comm, lammps_group
  !Characteristic of each  group
  integer, dimension(:),allocatable :: lammps_size, lammps_rank, all_rang
  integer :: no_of_lammps_group
  integer :: no_procs_of_lammps_group = 1
  integer :: ranks(1),new_group, new_comm
  integer :: w_rang, w_size
end MODULE vars_lammps


subroutine define_communicators ()
  use LAMMPS
  use vars_lammps
  use mpi
  use mod_mpi_phondy
  use gen_com_m , ONLY: rangph
  implicit none
  integer :: i,  istat, beggin,  endding
     if (mod(nb_procsph,  no_procs_of_lammps_group).ne.0) then
         if (rangph == 0) then
           write(6,*) 'The no of procs in each lammps group is not  divisor of the total number of procs'
           write(6,*) 'total number of procs, nb_procsph  :',  nb_procsph
           write(6,*) 'the  number of procs in each lammps group:',  no_procs_of_lammps_group
           stop
         end  if
     endif

     no_of_lammps_group= nb_procsph / no_procs_of_lammps_group

     if (allocated(lammps_group)) deallocate(lammps_group)
     if (allocated(lammps_comm))  deallocate(lammps_comm)
     if (allocated(lammps_size))  deallocate(lammps_size)
     if (allocated(lammps_rank))  deallocate(lammps_rank)
     allocate(lammps_group(0:no_of_lammps_group-1)      ,STAT=istat)
     allocate(lammps_comm (0:no_of_lammps_group-1)      ,STAT=istat)
     allocate(lammps_size (0:no_of_lammps_group-1)      ,STAT=istat)
     allocate(lammps_rank (0:no_of_lammps_group-1)      ,STAT=istat)

   call MPI_COMM_DUP(MPI_COMM_WORLD,orig_group,codeph)
   call MPI_COMM_GROUP(MPI_COMM_WORLD,orig_group,codeph)
!v1   orig_group=MPI_COMM_WORLD
!v1   call MPI_COMM_SIZE(orig_group,w_size,codeph)
!v1  call MPI_COMM_RANK(orig_group,w_rang,codeph)

!v1   call MPI_COMM_SPLIT(orig_group,w_rang,0,lammps_comm(w_rang),codeph)
!v1        call MPI_COMM_SIZE(lammps_comm(w_rang),lammps_size(w_rang),codeph)
!v1        call MPI_COMM_RANK(lammps_comm(w_rang),lammps_rank(w_rang),codeph)

     do i = 0 , no_of_lammps_group-1
        lammps_comm(i)= 0
     enddo

    if (allocated(all_rang)) deallocate(all_rang)
    allocate (all_rang(0:nb_procsph-1), STAT=istat)

    do i=0,nb_procsph-1
       all_rang(i)=i
    end do

   call MPI_COMM_GROUP(MPI_COMM_WORLD,orig_group,codeph)
   do i=0, no_of_lammps_group-1
     beggin  = no_procs_of_lammps_group*i
     endding = no_procs_of_lammps_group*(i+1)-1
     call MPI_GROUP_INCL(orig_group,no_procs_of_lammps_group, all_rang(beggin:endding),lammps_group(i),codeph)
     call MPI_COMM_CREATE(MPI_COMM_WORLD,lammps_group(i),lammps_comm(i),codeph)
     if (any(all_rang(beggin:endding) == rangph)) then
        call MPI_GROUP_SIZE(lammps_group(i),lammps_size(i),codeph)
        call MPI_GROUP_RANK(lammps_group(i),lammps_rank(i),codeph)
     endif
   enddo


  end  subroutine  define_communicators





subroutine init_potential_lammps()
  use LAMMPS
  use vars_lammps
  use gen_com_m, ONLY: rangph
  use mpi
  use mod_mpi_phondy
  implicit none
  integer :: i,error, beggin,  endding
  logical :: firsttime_lammps
  character*128 :: INPUT_LAMMPS_FILE

   INPUT_LAMMPS_FILE='in.lammps'
   if (rangph==0) write(*,*) "init potential lammps"


    do i=0, no_of_lammps_group-1
     beggin  = no_procs_of_lammps_group*i
     endding = no_procs_of_lammps_group*(i+1)-1
     if (any(all_rang(beggin:endding) == rangph)) then
        call lammps_open('lmp -log none -screen none', lammps_comm(i), lmp)
         if (rangph==0) write(*,*) "before init potential lammps"
         call lammps_file (lmp, INPUT_LAMMPS_FILE)
         if (rangph==0) write(*,*) "after init potential lammps"
        !debug the communicators:  write(*,*) "this is lammps_comm(",i,") for rangh :", rangph
     endif
   enddo

!v2     ranks(1)=rangph
!v2     call mpi_group_incl(orig_group,1,ranks,new_group,codeph)
!v2     call mpi_comm_create(MPI_COMM_WORLD,new_group,new_comm,codeph)
!v2     call lammps_open('lmp -log none -screen none', new_comm, lmp)


!v1   call lammps_open('lmp -log none -screen none', lammps_comm(w_rang), lmp)
!v1   write(*,*) rangph, lmp
!v1   call lammps_file (lmp, INPUT_LAMMPS_FILE)


     !call lammps_open('lmp -log none -screen none',  lmp)
     !call lammps_open_no_mpi('lmp -log none -screen none', lmp)
     !call lammps_file (lmp, INPUT_LAMMPS_FILE)
     !call lammps_open('lmp',MPI_COMM_WORLD,lmp)
   if (rangph==0) write(*,'("PHONDY: reading INPUT_LAMMPS_FILE file  :", (a))') INPUT_LAMMPS_FILE
   firsttime_lammps= .TRUE.

   if (rangph==0) write(*,'("PHONDY: LAMMPS force field init done")')

end subroutine init_potential_lammps



subroutine calcforce_lammps(NATOMS,posART,boxl,tmp_force,tmp_pos,pot_energy)
  use LAMMPS
  use vars_lammps
  use tab_imm_m, ONLY: ityp
  use gen_com_m, ONLY : rangph,im,umass
  use mod_mpi_phondy
  use mpi
  use phondy_in_ndm_module, ONLY: VECSIZE, firsttime_lammps, position_conversion_lammps,energy_conversion_lammps, cm_phondy_at, passage

  implicit none
  integer i,k,num,ierr
!  type (C_ptr) :: lmp
  integer, save :: cnt
  integer, intent(in)                               :: NATOMS
  real(kind=8), intent(in), dimension(3*NATOMS) :: posART
  real(8),dimension(3), intent(inout)  :: boxl
  real(8), intent(out), dimension(3*NATOMS) :: tmp_force, tmp_pos
  real (C_double), pointer :: energy => NULL()
  real(8), intent(out) :: pot_energy
  !real(8), intent(out), dimensiMPI_COM_WORLD can be changedon(VECSIZE), target:: tmp_force
  !real(8), intent(out), dimension(3*NATOMS), target:: tmp_pos
  real(8),  dimension(:), allocatable:: g_pos
  real(C_double),  dimension(:,:), pointer:: for_tmp => NULL()
  real(C_double),  dimension(:), pointer:: lammps_mass_per_type => NULL()
  !real(C_double),  dimension(:), allocatable:: lammps_mass
  integer, dimension(:), allocatable :: lammps_types
  real(8) :: xdummy

  double precision, dimension(:), allocatable :: pos_lammps
  double precision, dimension(:), allocatable :: force_lammps
  real(kind=8),dimension(3) :: box
  integer :: itemp

  box(:) = boxl(:)

  if (allocated(pos_lammps)) deallocate (pos_lammps) ; allocate(pos_lammps(VECSIZE), stat=ierr)

  if (firsttime_lammps) then
     if (rangph==0) write(*,'("PHONDY: Number of atoms   :",i7)') NATOMS
     num=lammps_get_natoms(lmp)
     if ((NATOMS /= im).or.(num /= im)) then
       write(*,*) 'Big problem: gin and lammps files contain different number of atoms'
     end if
     !if(allocated(lammps_mass)) deallocate(lammps_mass)
     !lammps_mass=>NULL()
     call lammps_extract_atom (lammps_mass_per_type, lmp, 'mass')
     if ( .not. associated(lammps_mass_per_type)) write(*,*) 'WARNING mass not associated'
     cnt=0
     !call lammps_gather_atoms (lmp, 'mass', 1, lammps_mass)

     call lammps_gather_atoms(lmp, 'type', 1, lammps_types)

     if (num /= size(lammps_types)) then
        write(*,*) 'WARNING:  the atoms type is not correctly read in the LAMMPS wrapper ndm_lammps'
     end if

     if (allocated(cm_phondy_at)) deallocate(cm_phondy_at)
     allocate (cm_phondy_at(num))

     do i=1,num
       cm_phondy_at (i) = lammps_mass_per_type(lammps_types(i))*umass
     end do

     do i=1,im
      if (lammps_types(i) - ityp(i) > 0) then
        write(*,*) 'Inconsistencies between *.lmp and *.gin file. Fix it! No any other choice!'
        stop
      end if
     end do

     firsttime_lammps = .FALSE.
  endif
  cnt=cnt+1

  do i=1, NATOMS
	pos_lammps(3*i-2) = posART(i)
	pos_lammps(3*i-1) = posART(i+NATOMS)
	pos_lammps(3*i  ) = posART(i+NATOMS*2)
  enddo


  ! Put the coordinates to LAMMPS
  call lammps_scatter_atoms (lmp, 'x',  pos_lammps)

  ! Call LAMMPS to compute energy and forces
  call lammps_command (lmp, 'run 0')

  ! Extract energy from LAMMPS
  call lammps_extract_compute (energy, lmp, 'thermo_pe',0,0)

  ! The energy  - we transform in eV

  pot_energy = energy*energy_conversion_lammps

  ! Extract forces from LAMMPS
  !v call lammps_gather_atoms (lmp, 'f', 3, force_lammps)
  call lammps_extract_atom (for_tmp, lmp, 'f')

  if (allocated(force_lammps)) deallocate(force_lammps)
  allocate(force_lammps(VECSIZE))

  do i=1,NATOMS
     force_lammps(3*i-2)=for_tmp(1,i)
     force_lammps(3*i-1)=for_tmp(2,i)
     force_lammps(3*i  )=for_tmp(3,i)
  end do


  do i=1, NATOMS
     tmp_force(i)          = force_lammps(3*i-2)*energy_conversion_lammps*position_conversion_lammps
     tmp_force(i+NATOMS)   = force_lammps(3*i-1)*energy_conversion_lammps*position_conversion_lammps
     tmp_force(i+NATOMS*2) = force_lammps(3*i  )*energy_conversion_lammps*position_conversion_lammps
  enddo
  ! call mpi_barrier(MPI_COMM_WORLD,codeph)
  return
end subroutine calcforce_lammps





