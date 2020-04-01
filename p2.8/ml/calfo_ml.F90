module calfo_ml_mod
        use alloc_typ_ml_mod 
        implicit none
        contains
!!this code is copyrighted @mihai-cosmin.marinica@cea.fr
!$-------------------------------------------------------------
subroutine md_init_potential_ml
!$-------------------------------------------------------------
! Here the ml potential is initialized. This subroutine is used
! by MD program and is called in NDM's init.F90 sub
use gen_com_m, only: A2cm, rvois, umass, dmtype
use ml_in_ndm_module, only: prepare_factorial,r_cut, rangml, periodic_table_element, fix_no_of_elements, fix_type_to_periodic
use var_pot
use time_measure, only: temps_energy, temps_force, temps_descripteurs, temps_neigh, temps_stress
use read_ml_file_mod
use compute_descriptors_mod
use snap
implicit none
integer :: i
!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ntyp and cm(:ntyp) not read from input ................!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!HERE ntyp should be changed ....
call read_ml_file()
ntyp=fix_no_of_elements
call alloc_typ_ml

!HERE ntyp and cm set by hand ...
!MAYBE cm should be initialised from ml file ...
!W mass
do i=1,ntyp
  cm(:ntyp)=periodic_table_element(fix_type_to_periodic(i))%mass
end do
!Fe mass
!cm(:ntyp)=periodic_table_element(26)%mass
if (rangml==0) then
   write(6,'("ML: the initial type and masses in uam: ", i7)') ntyp
   !write(6,*) "ML: mass setting", cm(:ntyp)
end if
cm(:ntyp)=cm(:ntyp)*umass
!here im not sure. Maybe rvois should be read in * .din and rue_pot fixed to t_cut of descipors
!like that rue_pot(:)=r_cut*A2cm
!
rvois=r_cut*A2cm
rue_pot(:)=rvois
do i=1,npair
    if (typ_pot_pair(i)==ipotentiel) rue_pair(i)=rue_pot(ipotentiel)
end do

iewald=0; l3c=.false. ; r3cm=0.d0
allocate(typ_and_pot(ntyp,npotmax))
typ_and_pot(1:ntyp,ipotentiel)=.true.
typ_pot_pair(:)=ipotentiel

if (ipotentiel /= 20) then
      if (rangml==0) write(6,'("ML: Fatal Error ipotentiel should be 20 in din file")')
      stop "ipotentiel error in md_init_potential_ml"
end if
lpotentiel(ipotentiel)=.true.
rue_pot(ipotentiel)=rvois
rue_pair(:) = rue_pot(ipotentiel)
rumax=rvois
!initialisation ... very general
call prepare_factorial()
call init_descriptors

!initialisation ... only for snap
call md_allocate_snap_params()
if (dmtype /= 18) then
   call read_parameters_snap
end if

!initialization of time counter ...
temps_force=0.d0
temps_neigh=0.d0
temps_energy=0.d0
temps_descripteurs=0.d0
temps_stress=0.d0

return
end subroutine md_init_potential_ml
!<-------------------------------------------------------------



!$-------------------------------------------------------------
subroutine md_init_config_ml
!$-------------------------------------------------------------
! Here the ml potential is initialized. This subroutine is used
! by MD program and is called in NDM's init.F90 sub after that configuration was read
use ml_in_ndm_module, only: iconf_data, iconf_data_train, iconf_data_test, md_iconf, write_desc
use gen_com_m, only: im, imm, volu, A2cm, dmtype
use derived_types, only:config_real, config_desc
use var_pot, only: ntyp
use snap
implicit none


!C_DEBUG if (ntyp > 1) then
!C_DEBUG  write(6,*) 'The driver MD - Milady not (yet) implemented for ntyp >1', ntyp
!C_DEBUG  stop 'fatal ntyp in md_init_config_ml'
!C_DEBUG end if

!initialisation ... very general
iconf_data=1
iconf_data_train=1
iconf_data_test=0
! config_real and condif_desc objects
if (allocated(config_real)) deallocate(config_real) ; allocate(config_real(iconf_data))
if (allocated(config_desc)) deallocate(config_desc) ; allocate(config_desc(iconf_data))

md_iconf=iconf_data
if (imm .gt. im) then
  imm=im
  write(6,'("ML: WARNING imm should be resized to im. MiLaDy cannot work otherwise")')
  stop 'md_init_potential_ml imm and im sizes are diffrent'
end if

call md_allocate_snap_desc()

config_real(md_iconf)%volume = volu/A2cm**3
config_real(md_iconf)%has_energy = .true.
config_real(md_iconf)%has_force = .true.
config_real(md_iconf)%has_stress = .true.

if (write_desc.and.(dmtype /= 18)) then
   write(6,'("ML: writting descriptors using MD is not possible ... yet!")')
   stop "md_init_config_ml writting descriptor during MD"
end if


if (allocated(config_real(md_iconf)%itype))    deallocate(config_real(md_iconf)%itype)    ;  allocate(config_real(md_iconf)%itype(im))
if (allocated(config_real(md_iconf)%pos_cart)) deallocate(config_real(md_iconf)%pos_cart) ;  allocate(config_real(md_iconf)%pos_cart(3,im))
if (allocated(config_real(md_iconf)%pos_crst)) deallocate(config_real(md_iconf)%pos_crst) ;  allocate(config_real(md_iconf)%pos_crst(3,im))
if (allocated(config_real(md_iconf)%force))    deallocate(config_real(md_iconf)%force)    ;  allocate(config_real(md_iconf)%force(3,im))
if (allocated(config_real(md_iconf)%mass_per_type))    &
                                               deallocate(config_real(md_iconf)%mass_per_type)  &
                                                                                          ;  allocate(config_real(md_iconf)%mass_per_type(ntyp))



return
end subroutine md_init_config_ml
!<-------------------------------------------------------------



!$-------------------------------------------------------------
subroutine put_ndm_into_ml_config(iconf)
!$-------------------------------------------------------------
use gen_com_m, only: im, at, bg, imm, A2cm, erg2ev, umass, volu
use tab_imm_m, only : xp, fp, ityp
use var_pot, only: ntyp, cm
use derived_types, only:config_real
implicit none
integer, intent(in) :: iconf


config_real(iconf)%volume = volu/A2cm**3
config_real(iconf)%ntypes=ntyp
config_real(iconf)%nat=imm
im=imm
config_real(iconf)%itype(1:imm) = ityp(1:imm)
config_real(iconf)%pos_cart(1:3,1:imm) = xp(1:3,1:imm)/A2cm
!config_real(iconf)%pos_crst(3,1:imm) = xc(3,1:imm)/A2cm
config_real(iconf)%force(1:3,1:imm) = fp(1:3,1:imm)*(A2cm*erg2ev)
config_real(iconf)%cell = at(:,:)/A2cm
config_real(iconf)%bg_cell = bg(:,:)*A2cm
config_real(iconf)%mass_per_type(:) = cm(:)/umass
!future config_real(iconf)%prev_pos_cart(1:3,1:imm)=xpp(1:3,1:imm)/A2cm

return
end subroutine put_ndm_into_ml_config
!<-------------------------------------------------------------


!$-------------------------------------------------------------
subroutine put_ml_config_into_ndm (iconf)
!$-------------------------------------------------------------
use gen_com_m, only: at, bg , imm, A2cm, erg2ev, umass, volu
use tab_imm_m, only : xp,  fp,  ityp
use derived_types, only:config_real
use var_pot, only: ntyp, cm
implicit none
integer, intent(in) :: iconf

imm = config_real(iconf)%nat
ityp(1:imm)=config_real(iconf)%itype(1:imm)
xp(1:3,1:imm)=A2cm*config_real(iconf)%pos_cart(1:3,1:imm)
!xc(3,1:imm)=A2cm*config_real(iconf)%pos_crst(3,1:imm)
fp(1:3,1:imm)= config_real(iconf)%force(1:3,1:imm) / (A2cm*erg2ev)
at(:,:)=config_real(iconf)%cell(:,:) * A2cm
bg(:,:)=config_real(iconf)%bg_cell(:,:) / A2cm
cm(:)=config_real(iconf)%mass_per_type(:)*umass
volu=config_real(iconf)%volume*A2cm**3
ntyp=config_real(iconf)%ntypes
!future xpp(1:3,1:imm)=config_real(iconf)%prev_pos_cart(1:3,1:imm)*A2cm
return
end subroutine put_ml_config_into_ndm
!<-------------------------------------------------------------


!$-------------------------------------------------------------
subroutine md_calfo_ml
!$-------------------------------------------------------------
!This subroutine provides the energy and forces for NDM MD main
!program.
!This subroutine is called in the main NDM's calfo.F90
!Output: (to be completed)
!         potist and fp ? through NDM module
#ifdef PARAML
   use mpi
   use mod_mpi_ml
#endif

USE T_kind_param_m, ONLY:  double
use gen_com_m, only: potist, sig, ev2erg, erg2ev, A2cm, evA2dyn, rangml, dmtype
use tab_imm_m, only: fp
use ml_in_ndm_module, only: ml_type, ml_type_basis, &
                            md_iconf, &
                            prepare_factorial, allocate_ml, deallocate_ml
use compute_descriptors_mod
use snap
use derived_types, only: config_real
!use var_pot, only: rumax
use time_measure, only: temps_energy, temps_force, temps_descripteurs, temps_neigh, temps_stress

implicit none
double precision,dimension(:,:),allocatable :: xdesc_i
real(kind(0.d0)) :: temps1, temps2, temps3, tmp_val
integer :: ia

    !if (debug) write(6,*) 'init md_calfo_ml'
if (ml_type==ml_type_basis) then
    !ph if (rangml==0) write(*,*) 'here0', md_iconf, nox, noxyz
    call put_ndm_into_ml_config(md_iconf)
    !if (debug) write(6,*) 'was put_ml_config_into_ndm'

    !ph if (rangml==0) write(*,*) 'here1', nox, noxyz
    call test_if_config_is_small(md_iconf)
    !ph if (rangml==0) write(*,*) 'here2', nox,noxyz
    !if (debug) write(6,*) 'was test_if_config_is_small'
#ifdef PARAML
    temps1=MPI_Wtime()
#endif
    call calc_neighbours(md_iconf)
#ifdef PARAML
    temps2=MPI_Wtime()
#endif
    temps_neigh = temps_neigh + (temps2 - temps1)
    !ph if (rangml==0) write(*,*) 'here3', nox, noxyz
    !if (debug) write(6,*) 'was calc_neighbours'
    call allocate_ml
    !if (debug)
    !ph if (rangml==0) write(6,*) 'here4 was allocate_ml'

#ifdef PARAML
    temps1=MPI_Wtime()
#endif
    call compute_descriptors(xdesc_i, md_iconf)
#ifdef PARAML
    temps2=MPI_Wtime()
#endif
    temps_descripteurs = temps_descripteurs + (temps2 - temps1)
    !ph if (rangml==0) write(*,*) 'here5',  nox, noxyz
    call deallocate_ml
    !ph if (rangml==0) write(*,*) 'here6 ....'
    !call MPI_BARRIER

#ifdef PARAML
    temps1=MPI_Wtime()
#endif
    !debug write(6,*) 'ALL THAT', config_real(md_iconf)%has_energy, config_real(md_iconf)%has_force, config_real(md_iconf)%has_stress
    call md_snap_compute_energy(md_iconf)
#ifdef PARAML
    temps2=MPI_Wtime()
#endif
    call md_snap_compute_force(md_iconf)
#ifdef PARAML
    temps3=MPI_Wtime()
#endif
    temps_energy=temps_energy + (temps2 - temps1)
    temps_force=temps_force+ (temps3 - temps2)


#ifdef PARAML
    temps1=MPI_Wtime()
#endif
    call md_snap_compute_stress(md_iconf)
#ifdef PARAML
    temps2=MPI_Wtime()
#endif

    temps_stress=temps_stress+ (temps2 - temps1)
    !exwrite(*,*) temps_energy, temps_force
    !call put_ml_config_into_ndm(md_iconf)
    !put energy in NDM units
    potist = ene_snap*ev2erg
    !put energy in NDM units
    fp(:,:)=fp_snap(:,:) /(A2cm*erg2ev)
    sig(1,1) = stress_snap(1)
    sig(2,2) = stress_snap(2)
    sig(3,3) = stress_snap(3)
    sig(2,3) = stress_snap(4)
    sig(1,3) = stress_snap(5)
    sig(1,2) = stress_snap(6)
    sig(2,1) =sig(1,2)
    sig(3,1) =sig(1,3)
    sig(3,2) =sig(2,3)

    !put stres in NDM units
    sig(:,:) = sig(:,:)*1.d+09/evA2dyn
    !write(*,*) 'ENE_SNAP', ene_snap, maxval(fp_snap)
    !write(*,*) 'here5',  nox, noxyz, rumax
    !stop "testing clafo_ml"

end if
if (dmtype==5)  then
    if (rangml==0) then

     open(unit=531, file='COORD', status='unknown')
     open(unit=532, file='FORCE', status='unknown')

     ! Writting COORD ....
     write(531,'(a)') 'ITEM: TIMESTEP'
     write(531,'(a)') '0'
     write(531,'(a)') 'ITEM: NUMBER OF ATOMS'
     write(531,'(i6)') config_real(md_iconf)%nat
     write(531,'(a)') 'ITEM: BOX BOUNDS xy xz yz pp pp pp'
     write(531,'(3e23.15)') 0.d0, config_real(md_iconf)%cell(1,1), 0.d0
     write(531,'(3e23.15)') 0.d0, config_real(md_iconf)%cell(2,2), 0.d0
     write(531,'(3e23.15)') 0.d0, config_real(md_iconf)%cell(3,3), 0.d0
     tmp_val= config_real(md_iconf)%cell(1,2)**2 + config_real(md_iconf)%cell(1,3)**2 + &
     config_real(md_iconf)%cell(2,3)**2 + config_real(md_iconf)%cell(2,1)**2 + &
     config_real(md_iconf)%cell(3,1)**2 + config_real(md_iconf)%cell(3,2)**2
     if (tmp_val >=1.d-15) then
       if (rangml==0) write(6,*) 'This save type is not implemented for triclinic box. Stop in calfo_ml'
       stop 'calfo_ml writting triclinic box error'
     end if
     write(531,'(a)') 'ITEM: ATOMS id xu yu zu'
    do ia=1,config_real(md_iconf)%nat
       write(531,'(i5,3f30.15)') ia, config_real(md_iconf)%pos_cart(:,ia)
    end do
    close(531, status='keep')

    ! Writting FORCE ....
    end if

     ! Writting COORD ....
     write(532,'(a)') 'ITEM: TIMESTEP'
     write(532,'(a)') '0'
     write(532,'(a)') 'ITEM: NUMBER OF ATOMS'
     write(532,'(i6)') config_real(md_iconf)%nat
     write(532,'(a)') 'ITEM: BOX BOUNDS xy xz yz pp pp pp'
     write(532,'(3e23.15)') 0.d0, config_real(md_iconf)%cell(1,1), 0.d0
     write(532,'(3e23.15)') 0.d0, config_real(md_iconf)%cell(2,2), 0.d0
     write(532,'(3e23.15)') 0.d0, config_real(md_iconf)%cell(3,3), 0.d0
     tmp_val= config_real(md_iconf)%cell(1,2)**2 + config_real(md_iconf)%cell(1,3)**2 + &
     config_real(md_iconf)%cell(2,3)**2 + config_real(md_iconf)%cell(2,1)**2 + &
     config_real(md_iconf)%cell(3,1)**2 + config_real(md_iconf)%cell(3,2)**2
     if (tmp_val >=1.d-15) then
       if (rangml==0) write(6,*) 'This save type is not implemented for triclinic box. Stop in calfo_ml'
       stop 'calfo_ml writting triclinic box error'
     end if
     write(532,'(a)') 'ITEM: ATOMS id fx fy fz '
    do ia=1,config_real(md_iconf)%nat
       write(532,'(i5,3f30.15)') ia, fp_snap(:,ia)
    end do
    close(532, status='keep')

    stop 'TEST FORCES IN ML'

end if
return
end subroutine md_calfo_ml
!<-------------------------------------------------------------
end module 
