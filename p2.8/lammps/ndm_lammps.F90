#ifdef LAMMPS_VERSION
MODULE vars_lammps
  use LAMMPS
  USE T_kind_param_m, ONLY:  double
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




subroutine read_lammps
  use gen_com_m, ONLY: rang,firsttime_lammps,iverbose
  use LAMMPS
  use vars_lammps
  character*128 :: INPUT_LAMMPS_FILE
!  type (C_ptr) :: lmp 


   INPUT_LAMMPS_FILE='in.lammps' 
   if (iverbose ==0) then
      call lammps_open_no_mpi('lmp -log none -screen none', lmp)
   else
      call lammps_open_no_mpi('lmp  -screen none', lmp)
   endif
   if (rang==0) write(*,*) "before init potential lammps"
   call lammps_file (lmp, INPUT_LAMMPS_FILE)
   if (rang==0) write(*,*) "after init potential lammps"

   if (rang==0) write(*,*) "init potential lammps"



   if (rang==0)   write(*,'("NDM: reading INPUT_LAMMPS_FILE file  :", (a))') INPUT_LAMMPS_FILE
   firsttime_lammps= .TRUE.

   write(*,'("NDM: LAMMPS force field init done")')

 end subroutine read_lammps



 subroutine calcforce_lammps2
   use vars_lammps
   use LAMMPS
   use tab_imm_m, ONLY: ityp,xp,fp
   use var_pot,only:ntyp,cm
   use gen_com_m, ONLY : im,umass,firsttime_lammps,potist,energy_conversion_lammps,position_conversion_lammps,&
        &sig,it,itesigma,rskin,lpr,at,unitP,pressure_conversion_lammps
   use mat_util
   !  use mod_mpi_phondy
   !  use mpi


   implicit none
   integer i,k,num,ierr
   real (C_double), pointer :: energy => NULL()
   real(C_double), dimension(:), pointer :: p_tensor=>NULL()
   real(C_double),  dimension(:,:), pointer:: for_tmp => NULL()
   real(C_double),  dimension(:), pointer:: lammps_mass_per_type => NULL()
   integer, dimension(:), allocatable :: lammps_types

   double precision, dimension(:), allocatable :: pos_lammps
   !  double precision, dimension(:), allocatable :: force_lammps
   double precision, dimension(:,:), allocatable,save :: axlmp
   real(kind=8),dimension(3) :: box
   integer :: itemp,iti,ic
   logical::lrun0
   real*8 :: rdiff

   logical lrotated,upper
   real(double)::xhi,yhi,zhi,xy,xz,yz,xlo,ylo,zlo
   real(double), dimension(3,3)::at_lammps,passage,passage_inv
   character :: lmpcom*300
   real(double), dimension(3) :: tmp_coord_i,new_tmp_coord_i
   real(double), dimension(3,3) ::   sigboxlmp,siginter



   !  box(:) = boxl(:)

   if (firsttime_lammps) then
      if (allocated(pos_lammps)) deallocate (pos_lammps)
      allocate(pos_lammps(3*im), stat=ierr)

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

   !IF variable box, box should be changed into triclin lammps format, positions must be transformed 
   if (lpr.eqv..true.) then
      !     call transformat2lammps

      call is_upper_triangular(at,upper)
      if (.not.upper) then
         call convert_cell (at,at_lammps,passage)
         call matinv(passage, passage_inv)
      else
         at_lammps=at
         passage(:,:)=0
         do ic=1,3
            passage(ic,ic)=1
         end do
         passage_inv(:,:)=passage(:,:)
      end if

      !building  lammps header
      xlo = 0.d0
      ylo = 0.d0
      zlo = 0.d0
      xhi = at_lammps(1,1)/position_conversion_lammps
      yhi = at_lammps(2,2)/position_conversion_lammps
      zhi = at_lammps(3,3)/position_conversion_lammps
      xy = at_lammps(1,2)/position_conversion_lammps
      xz = at_lammps(1,3)/position_conversion_lammps
      yz = at_lammps(2,3)/position_conversion_lammps

      write(lmpcom,'(A,2F17.12,A,2F17.12,A,2F17.12,A,F17.12,A,F17.12,A,F17.12,A)')'change_box all x final',&
           &xlo,xhi,' y final', ylo,yhi,' z final',zlo,zhi, ' xy final',xy,' xz final',xz,' yz final',yz, ' remap '
!      write(6,*)lmpcom
      !  stop
!      write(6,*)'OHOH'

      call lammps_command (lmp, lmpcom)
      !    call lammps_command (lmp,'change_box all x final 0.0 10.0 y final 0. 10.0 z &
      !&final 0. 10.0 xy final 2.0 xz final 2.0 yz final 2.0 remap ' )
!      write(6,*)passage
      !  stop
      do i=1, im
         tmp_coord_i = xp(:,i)
         new_tmp_coord_i = matmul(passage,tmp_coord_i)/position_conversion_lammps
         pos_lammps(3*i-2) = new_tmp_coord_i(1)
         pos_lammps(3*i-1) = new_tmp_coord_i(2)
         pos_lammps(3*i  ) = new_tmp_coord_i(3)
      enddo


   else  
      do i=1, im
         pos_lammps(3*i-2) = xp(1,i)/position_conversion_lammps
         pos_lammps(3*i-1) = xp(2,i)/position_conversion_lammps
         pos_lammps(3*i  ) = xp(3,i)/position_conversion_lammps
      enddo
   endif
   !     write(6,*)'calfolammps2'

   ! Put the coordinates to LAMMPS
!   write(6,*)pos_lammps


   call lammps_scatter_atoms (lmp, 'x',  pos_lammps)
   !     write(6,*)'calfolammps3'

   ! Call LAMMPS to compute energy and forces
   if ((firsttime_lammps).or.(lpr)) then
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
   potist=energy*energy_conversion_lammps 

   if (lpr) then
      call lammps_extract_compute (p_tensor, lmp, 'thermo_press',0,1)
      sigboxlmp(1,1)=p_tensor(1)*pressure_conversion_lammps
      sigboxlmp(2,2)=p_tensor(2)*pressure_conversion_lammps
      sigboxlmp(3,3)=p_tensor(3)*pressure_conversion_lammps
      sigboxlmp(1,2)=p_tensor(4)*pressure_conversion_lammps
      sigboxlmp(1,3)=p_tensor(5)*pressure_conversion_lammps
      sigboxlmp(2,3)=p_tensor(6)*pressure_conversion_lammps
      sigboxlmp(2,1)=sigboxlmp(1,2)
      sigboxlmp(3,1)=sigboxlmp(1,3)
      sigboxlmp(3,2)=sigboxlmp(2,3)
      siginter=matmul(sigboxlmp,passage)
      sig=matmul(passage_inv,siginter)
!         do ic = 1, 3
!            write (6, '(I1,3(A,I1),A,3G18.10)') ic,' sigma POST LAMMPS potentiel (1,', ic, ') (2,', ic, &
!                 ') (3,', ic, ') =',sig(1:3,ic)*unitP
!         end do
     
      call lammps_extract_atom (for_tmp, lmp, 'f')
      do i=1,im
         tmp_coord_i = for_tmp(:,i)*energy_conversion_lammps/position_conversion_lammps 
         fp(:,i)= matmul(passage_inv,tmp_coord_i)
         
!         fp(1,i)=for_tmp(1,i)*energy_conversion_lammps/position_conversion_lammps ! / (A2cm*erg2ev)
!         fp(2,i)=for_tmp(2,i)*energy_conversion_lammps/position_conversion_lammps ! / (A2cm*erg2ev)
!         fp(3,i)=for_tmp(3,i)*energy_conversion_lammps/position_conversion_lammps ! / (A2cm*erg2ev)
      end do

      
   else
      
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
      call lammps_extract_atom (for_tmp, lmp, 'f')
      
      
      do i=1,im
         fp(1,i)=for_tmp(1,i)*energy_conversion_lammps/position_conversion_lammps ! / (A2cm*erg2ev)
         fp(2,i)=for_tmp(2,i)*energy_conversion_lammps/position_conversion_lammps ! / (A2cm*erg2ev)
         fp(3,i)=for_tmp(3,i)*energy_conversion_lammps/position_conversion_lammps ! / (A2cm*erg2ev)
      end do
   end if

   ! Extract forces from LAMMPS  
   !v call lammps_gather_atoms (lmp, 'f', 3, force_lammps)
   return
 end subroutine calcforce_lammps2





#endif
