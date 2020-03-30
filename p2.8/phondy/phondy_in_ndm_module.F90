module phondy_in_ndm_module
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE T_kind_param_m, ONLY:  double
      use gen_com_m
      use tab_imm_m
      use var_pot
      use jqmod
#if(PARAPH)
      use mpi
      use mod_mpi_phondy
#endif
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      implicit none


      integer                                       :: it_phondy
      integer, parameter:: NN_MAX=1000
      integer :: nmat
      real(double) ::  convert_phondy,rcut_ph, rcut_ph_ang
      real(double),save :: epot0
      real(double),dimension(:),allocatable, save:: w
      real(double),dimension(:,:),allocatable, save:: fp0,matfor,densmodes
      complex(kind(1.d0)), dimension(:,:), allocatable :: dynmat
      integer, allocatable, dimension(:) :: u_local,v_local,imodes
      real(double), allocatable, dimension(:) :: mlocal
      integer :: imax
      real(double)   :: avogadro, electron,two_pi,unit_nu,units_ndm
      integer :: i_start,i_final
     integer ::  isave, iread, nmodes
     logical :: leigenvectors  ! If you want the eigenvalues
     logical :: lmodes
     real(double) :: nu_min,nu_max,width_dos,wmodes, weight=1.d0
     integer :: n_nu_up,nsitedos, nsitedos_loc
     integer, dimension(:), allocatable :: isitedos
     real(double), allocatable, dimension (:) :: dos
     logical :: lldos, l_thermo_atoms
     character*6 :: namefile
     integer :: inamefile
     !LAMMPS interface ...
     integer :: nat,VECSIZE
     character(len=128) :: fnamtin_lammps
     logical ::  firsttime_lammps
     real(8) :: energy_conversion_lammps,  position_conversion_lammps
     real(kind=8) , allocatable, dimension(:)  ::  posa, forca
     real(kind=8) , allocatable, dimension(:)  ::  cm_phondy, cm_phondy_at
     real(8)               :: energy, boxl(3)

     !Qpoints part
     logical :: lq_points
     integer :: iconf_data_ph, no_of_qpoints
     real(kind=8), dimension(:,:), allocatable :: q_point, q_eigenvalues
     real(kind=8), dimension(:), allocatable :: weight_q_point
     integer, parameter :: iconf_ini=1, iconf_big=2
     logical :: debug_ph
     real(kind(0.d0)), dimension(:,:), allocatable :: passage, passage_inv


#if(PARAPH)
     integer :: nratio1,nratio2, nrest, iproc,nb_elements,itempproc


#endif

 contains


subroutine init_phondy()

   implicit none
   it_phondy=0

!  pi=4*datan(1.D0)
   two_pi=2.0d0*4.d0*datan(1.d0)
   avogadro=6.0221367d0
   electron=1.60217733d0
   !hplanck=6.62618
   unit_nu=dsqrt(avogadro*electron/1.D3)
   unit_nu=unit_nu/two_pi
   units_ndm = erg2eV*A2cm*A2cm*umass

! LAMMPS interface ...
   energy_conversion_lammps=1.d0
   position_conversion_lammps=1.d0
   firsttime_lammps=.true.

   call set_limit_for_atoms(rangph,nb_procsph, im,i_start,i_final)
   !rcut_ph=2.d0*rvois!*A2cm
   rcut_ph=rvois!*A2cm
   if (rangph==0) write(*,*) 'Rvois.....(SHOULD BE  AT LEAST 2 x Rue pott .: ', rcut_ph/A2cm

return
end   subroutine init_phondy




subroutine allocate_phondy()
   implicit none
   integer :: i
   nat=im
   nmat=3*nat
  ! isave=0 ! isave=1 we will store de din matrix of HDD for a subequent diagonalization
           ! isave=0 nothing stored.
           ! isave=2 just computing dyn matrix by MPI and diag by threading. nothing stored.
  ! iread=0 ! iread=1 reading the dyn matrix from previous run
           ! iread=0 everything is computed

   VECSIZE=3*nat
   if (allocated(fp0)) deallocate(fp0) ; allocate (fp0(3,imm))

   if (allocated(posa)) deallocate(posa) ; allocate (posa(3*nat))
   if (allocated(forca)) deallocate(forca) ; allocate(forca(3*nat))

# if (LAMMPS_VERSION)
#else
  if (allocated(cm_phondy)) deallocate(cm_phondy) ; allocate (cm_phondy(ntyp))
  if (allocated(cm_phondy_at)) deallocate(cm_phondy_at) ; allocate(cm_phondy_at(im))
  cm_phondy(:) = cm(:)
  do i =1,im
    cm_phondy_at(i)=cm_phondy(ityp(i))
  end do
#endif
return
end   subroutine allocate_phondy

















 subroutine  init_phondy_in_ndm_module
!
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    implicit none
    integer   :: ic_local,iatom


   it_phondy=0



   do ic_local=1,3*im
      iatom = MOD(ic_local,im)
      if (iatom==0) iatom=im
   end do


 return
!
end subroutine init_phondy_in_ndm_module


end module phondy_in_ndm_module



module diago_scalapack_real
          !integer :: N
          !integer :: iproc
          integer, dimension(:), allocatable :: DESCZ
          !real(kind=8) :: W(:)
          real(kind=8), dimension(:,:), allocatable :: Z
end module diago_scalapack_real





subroutine set_limit_for_atoms(rangph,nb_procsph, im,i_start_at,i_final_at)
     implicit none
     integer, intent(in) :: im           ! numbers of atoms to be distributed on procs
     integer, intent(inout)  :: rangph, nb_procsph   ! rang of the proc for MPI
     integer, intent(out) :: i_start_at, i_final_at
#if(PARAPH)
     integer :: nratio1,nratio2, nrest
#endif


#if(PARAPH)
if (im >= nb_procsph) then
  nrest=mod(im,nb_procsph)
  nratio1=im/nb_procsph+1
  nratio2=im/nb_procsph

  if (rangph<nrest) then
     i_start_at=rangph*nratio1+1
     if (rangph /= (nb_procsph-1) ) i_final_at=(rangph+1)*nratio1
  end if

  if (rangph>=nrest) then
    i_start_at = nrest*nratio1 + (rangph-nrest)*nratio2 + 1
    if (rangph /= (nb_procsph-1) ) i_final_at=nrest*nratio1 + (rangph-nrest+1)*nratio2
  endif

  if (rangph == (nb_procsph-1) ) i_final_at= im
else
  if (im == 1) then
     if (rangph==0) then
       i_start_at=1
       i_final_at=1
     else
       i_start_at=0
       i_final_at=0
     end if
  else
    if (rangph+1 <= im ) then
       i_start_at=rangph+1
       i_final_at=rangph+1
    else
       i_start_at=0
       i_final_at=0
    end if
  end if
end if
#else

rangph=0
nb_procsph=1
i_start_at=1
i_final_at=im

#endif



return
end subroutine set_limit_for_atoms
