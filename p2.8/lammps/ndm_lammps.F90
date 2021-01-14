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
  integer :: w_rang, w_size
  integer::MPI_COMM_lammps
end MODULE !vars_lammps


module lammps_util_mod
  use gen_com_m, ONLY: rang,firsttime_lammps
  !use gen_com_m,only:
        use LAMMPS
        use vars_lammps
#ifdef PARA
  use mod_para, only: nprocspace,mpi_comm_space,ierr,mpi_comm_world
#else
  use mod_para, only: nprocspace
#endif

        implicit none

      contains

subroutine read_lammps(inplammps)

  use LAMMPS
  use vars_lammps
  character(*),optional :: inplammps
    character*128 :: INPUT_LAMMPS_FILE
!  type (C_ptr) :: lmp

    INPUT_LAMMPS_FILE='in.lammps'
    if (present(inplammps))  INPUT_LAMMPS_FILE=inplammps
#ifdef PARA
   if (nprocspace==1) then
    call lammps_open_no_mpi('lmp -log none -screen none', lmp)
    write(*,*) "LAMMPS OPEN_NO_MPI_",rang
    if (rang==0) write(*,*) "before init potential lammps"
    call lammps_file (lmp, INPUT_LAMMPS_FILE)
    if (rang==0) write(*,*) "after init potential lammps"
    if (rang==0) write(*,*) "init potential lammps"
    if (rang==0)   write(*,'("NDM: reading INPUT_LAMMPS_FILE file  :", (a))') INPUT_LAMMPS_FILE
  else
     !call define_communicators_lammps
     call MPI_COMM_DUP(MPI_COMM_SPACE,MPI_COMM_lammps,ierr)     
     call lammps_open('lmp -log none -screen none', MPI_COMM_lammps, lmp)
     write(*,*) "LAMMPS OPEN_MPI_",rang
    call lammps_file (lmp, INPUT_LAMMPS_FILE)
  end if
#else

    call lammps_open_no_mpi('lmp -log none -screen none', lmp)
    write(*,*) "LAMMPS OPEN_NO_MPI-SEQ"
    if (rang==0) write(*,*) "before init potential lammps SEQ"
    call lammps_file (lmp, INPUT_LAMMPS_FILE)
    if (rang==0) write(*,*) "after init potential lammps"
    if (rang==0) write(*,*) "init potential lammps"
    if (rang==0)   write(*,'("NDM: reading INPUT_LAMMPS_FILE file  :", (a))') INPUT_LAMMPS_FILE
#endif

  firsttime_lammps= .TRUE.
  if (rang==0)write(*,'("NDM: LAMMPS force field init done")')
!  stop

 end subroutine read_lammps


subroutine calcforce_lammps2 (im,imm,xp,ityp,fp,potislammps)
  use vars_lammps
  use LAMMPS
!  use tab_imm_m, ONLY: ityp,xp,fp
  use var_pot,only:ntyp,cm
  use gen_com_m, ONLY : umass,firsttime_lammps,energy_conversion_lammps,position_conversion_lammps,sig,it,itesigma,rskin&
       &,pressure_conversion_lammps
!  use mod_para_phondy
!  use mpi


  implicit none

  integer,intent(in)::im,imm
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
!     write(6,*)'calfolammps0',rang
     num=lammps_get_natoms(lmp)
     if (num /= im) then
        write(*,*) 'Big problem: gin and lammps files contain different number of atoms'
        stop
     end if
     call lammps_gather_atoms(lmp, 'type', 1, lammps_types)
!     write(6,*)'calfolammps1',rang
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
     firsttime_lammps=.false.
  endif
!     write(6,*)'calfolammps1.3'
  do i=1, im
	pos_lammps(3*i-2) = xp(1,i)/position_conversion_lammps
	pos_lammps(3*i-1) = xp(2,i)/position_conversion_lammps
	pos_lammps(3*i  ) = xp(3,i)/position_conversion_lammps
!        write(500+rang,*) rang,i,pos_lammps(3*i-2),pos_lammps(3*i-1),pos_lammps(3*i  )
enddo
!     write(6,*)'calfolammps2',rang

  ! Put the coordinates to LAMMPS
!     write(100+rang,*)'X',pos_lammps
  call lammps_scatter_atoms (lmp, 'x',  pos_lammps)
 !    write(6,*)'calfolammps3',rang

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
!        write(6,*)'LRUN0'
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
!     write(6,*)'calfolammps4',rang


  ! Extract energy from LAMMPS
  call lammps_extract_compute (energy, lmp, 'thermo_pe',0,0)
!     write(6,*)'calfolammps5'
  potislammps=energy*energy_conversion_lammps
  if (mod(it,itesigma)==0) then
     call lammps_extract_compute (p_tensor, lmp, 'thermo_press',0,1)
!     write(6,*)'calfolammps6'
!  pot_energy = energy*energy_conversion_lammps
!     write (6,*)p_tensor
      sig(1,1)=p_tensor(1)*pressure_conversion_lammps
     sig(2,2)=p_tensor(2)*pressure_conversion_lammps
     sig(3,3)=p_tensor(3)*pressure_conversion_lammps
     sig(1,2)=p_tensor(4)*pressure_conversion_lammps
     sig(1,3)=p_tensor(5)*pressure_conversion_lammps
     sig(2,3)=p_tensor(6)*pressure_conversion_lammps
     sig(2,1)=sig(1,2)
     sig(3,1)=sig(1,3)
     sig(3,2)=sig(2,3)
  end if

  ! Extract forces from LAMMPS
  call lammps_gather_atoms (lmp, 'f', 3, force_lammps)
 ! call lammps_extract_atom (for_tmp, lmp, 'f')
!!$En séquentiel, gather_atoms et extract_atom donnent la même chose.
!!$En parallèle :
!!$1/extract_atom donne des choses différentes sur chaque proc
!!$2/gather_atoms donne des choses égales sur tous les procs
!!$3/gather_atoms donne des choses égales au gather ou extract du séquentiel
!!$4/Il semble que ce qui change dans les différents extract_atoms soit l'ordre des atomes (on dirait, il y a des nombres qui se ressemblent). Il faut peut-être els réarranger selon un indice interproc inconnu.

do i=1,im
   
     fp(1,i)=force_lammps(3*i-2)*energy_conversion_lammps/position_conversion_lammps ! / (A2cm*erg2ev)
     fp(2,i)=force_lammps(3*i-1)*energy_conversion_lammps/position_conversion_lammps ! / (A2cm*erg2ev)
     fp(3,i)=force_lammps(3*i)*energy_conversion_lammps/position_conversion_lammps ! / (A2cm*erg2ev)
  end do
  return
end subroutine calcforce_lammps2


end module


#endif
