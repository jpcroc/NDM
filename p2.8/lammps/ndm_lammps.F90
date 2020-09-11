#ifdef LAMMPS_VERSION
MODULE vars_lammps
  USE T_kind_param_m, ONLY:  double
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
end MODULE !vars_lammps


module lammps_util_mod
        use gen_com_m
        use LAMMPS
        use vars_lammps
        implicit none
        contains
<<<<<<< HEAD

=======
!
>>>>>>> f48eb679ca3c0a28dbdd70a92b7134f26bc0b4ca
subroutine read_lammps
  use gen_com_m, ONLY: rang,firsttime_lammps
  use LAMMPS
  use vars_lammps
  character*128 :: INPUT_LAMMPS_FILE
!  type (C_ptr) :: lmp 


   INPUT_LAMMPS_FILE='in.lammps' 
   call lammps_open_no_mpi('lmp -log none -screen none', lmp)
   if (rang==0) write(*,*) "before init potential lammps"
   call lammps_file (lmp, INPUT_LAMMPS_FILE)
   if (rang==0) write(*,*) "after init potential lammps"

   if (rang==0) write(*,*) "init potential lammps"



   if (rang==0)   write(*,'("NDM: reading INPUT_LAMMPS_FILE file  :", (a))') INPUT_LAMMPS_FILE
   firsttime_lammps= .TRUE.

   write(*,'("NDM: LAMMPS force field init done")')

 end subroutine read_lammps



subroutine calcforce_lammps2 (im,xp,ityp,fp,potislammps)
  use vars_lammps
  use LAMMPS
!  use tab_imm_m, ONLY: ityp,xp,fp
  use var_pot,only:ntyp,cm
  use gen_com_m, ONLY : umass,firsttime_lammps,energy_conversion_lammps,position_conversion_lammps,sig,it,itesigma,rskin
!  use mod_para_phondy
!  use mpi


  implicit none

  integer,intent(in)::im
  integer, intent (in),allocatable:: ityp(:)
  real(double),intent(in),allocatable::xp(:,:)
  real(double), intent(inout),allocatable::fp(:,:)
  real(double),intent(out)::potislammps
  integer i,k,num,ierr
  real (C_double), pointer :: energy => NULL()
  real(C_double), dimension(:), pointer :: p_tensor=>NULL()
  real(C_double),  dimension(:,:), pointer:: for_tmp => NULL()
  real(C_double),  dimension(:), pointer:: lammps_mass_per_type => NULL()
  integer, dimension(:), allocatable :: lammps_types

  double precision, dimension(:), allocatable :: pos_lammps
  double precision, dimension(:), allocatable :: force_lammps
  double precision, dimension(:,:), allocatable,save :: axlmp
  real(kind=8),dimension(3) :: box
  integer :: itemp,iti,ic
  logical::lrun0
  real*8 :: rdiff
  
!  box(:) = boxl(:)

  if (allocated(pos_lammps)) deallocate (pos_lammps)
  allocate(pos_lammps(3*im), stat=ierr)
  if (firsttime_lammps) then
!     write(6,*)'calfolammps0'
     num=lammps_get_natoms(lmp)
     if (num /= im) then
        write(*,*) 'Big problem: gin and lammps files contain different number of atoms'
        stop
     end if 
     call lammps_gather_atoms(lmp, 'type', 1, lammps_types)
!     write(6,*)'calfolammps1'
     if (num /= size(lammps_types)) then
        write(*,*) 'WARNING:  the atoms type is not correctly read in the LAMMPS wrapper ndm_lammps'
     end if 
!     write(6,*)'calfolammps1.1'
     do i=1,im
        if (ityp(i).ne.lammps_types(i)) then
           write(6,*)'erreur de transmission de type atome ',i,' type ndm lammps ',ityp(i),lammps_types(i)
           stop
        end if
     end do
!     write(6,*)'calfolammps1.2',im
     if(.not.allocated(axlmp))allocate(axlmp(3,im))
!     write(6,*)'calfolammps1.2.1'
     axlmp(:,1:im)=xp(:,1:im)
  endif
!     write(6,*)'calfolammps1.3'
  do i=1, im
	pos_lammps(3*i-2) = xp(1,i)/position_conversion_lammps
	pos_lammps(3*i-1) = xp(2,i)/position_conversion_lammps
	pos_lammps(3*i  ) = xp(3,i)/position_conversion_lammps
  enddo
!     write(6,*)'calfolammps2'

  ! Put the coordinates to LAMMPS

  call lammps_scatter_atoms (lmp, 'x',  pos_lammps)
!     write(6,*)'calfolammps3'

  ! Call LAMMPS to compute energy and forces
  if (firsttime_lammps) then
     lrun0=.true.
  else
     rdiff=0
     do i=1,im
        do ic=1,3
           rdiff=max(abs(axlmp(ic,i)-xp(ic,i))*1d8,rdiff)
        end do
     end do
     if (rdiff.ge.rskin) then
        axlmp(:,:)=xp(:,:)
        write(6,*)'LRUN0'
        lrun0=.true.
     end if
  end if
  
 if (lrun0.eqv..true.)then
     call lammps_command (lmp, 'run 0')
     lrun0=.false.
 else
    call lammps_command (lmp, 'run 1 pre no post yes')
 end if
!  call lammps_command (lmp, 'run 0')
!     write(6,*)'calfolammps4'


  ! Extract energy from LAMMPS
  call lammps_extract_compute (energy, lmp, 'thermo_pe',0,0)
!     write(6,*)'calfolammps5'
  potislammps=energy*energy_conversion_lammps 
  if (mod(it,itesigma)==0) then
     call lammps_extract_compute (p_tensor, lmp, 'thermo_press',0,1)
!     write(6,*)'calfolammps6'
!  pot_energy = energy*energy_conversion_lammps
!     write (6,*)p_tensor
     sig(1,1)=p_tensor(1)*energy_conversion_lammps/(position_conversion_lammps**3)
     sig(2,2)=p_tensor(2)*energy_conversion_lammps/(position_conversion_lammps**3)
     sig(3,3)=p_tensor(3)*energy_conversion_lammps/(position_conversion_lammps**3)
     sig(1,2)=p_tensor(4)*energy_conversion_lammps/(position_conversion_lammps**3)
     sig(1,3)=p_tensor(5)*energy_conversion_lammps/(position_conversion_lammps**3)
     sig(2,3)=p_tensor(6)*energy_conversion_lammps/(position_conversion_lammps**3)
     sig(2,1)=sig(1,2)
     sig(3,1)=sig(1,3)
     sig(3,2)=sig(2,3)
  end if
  
  ! Extract forces from LAMMPS  
  !v call lammps_gather_atoms (lmp, 'f', 3, force_lammps)
  call lammps_extract_atom (for_tmp, lmp, 'f')

!       write(6,*)'calfolammps7'
!    call lammps_extract_atom (vel_tmp, lmp, 'v')
!  write(6,*) vel_tmp

 ! if (allocated(force_lammps)) deallocate(force_lammps)
 ! allocate(force_lammps(3*im))

  do i=1,im
     fp(1,i)=for_tmp(1,i)*energy_conversion_lammps/position_conversion_lammps ! / (A2cm*erg2ev)
     fp(2,i)=for_tmp(2,i)*energy_conversion_lammps/position_conversion_lammps ! / (A2cm*erg2ev)
     fp(3,i)=for_tmp(3,i)*energy_conversion_lammps/position_conversion_lammps ! / (A2cm*erg2ev)
!     force_lammps(3*i-2)=for_tmp(1,i)
!     force_lammps(3*i-1)=for_tmp(2,i)
!     force_lammps(3*i  )=for_tmp(3,i)
  end do
!  write(6,*)
!  write(6,*)'fp1',fp(:,1), 'Z'
!  write(6,*)'fp2',fp(:,2)
!  write(6,*)'fp3',fp(:,3)
!  do i=1, NATOMS    
!     tmp_force(i)          = force_lammps(3*i-2)*energy_conversion_lammps*position_conversion_lammps
!     tmp_force(i+NATOMS)   = force_lammps(3*i-1)*energy_conversion_lammps*position_conversion_lammps
!     tmp_force(i+NATOMS*2) = force_lammps(3*i  )*energy_conversion_lammps*position_conversion_lammps
!  enddo
  ! call mpi_barrier(MPI_COMM_WORLD,codeph) 
  return
end subroutine calcforce_lammps2


end module


#endif
