module thermo_module

  use T_kind_param_m, ONLY:  double
  use phondy_in_ndm_module, ONLY: nsitedos, nsitedos_loc, nat

  implicit none

  real(double),dimension(:),allocatable, save:: Fmin,Fmin_cla,Smin,Smin_cla
  real(double),dimension(:,:),allocatable, save:: Fmin_loc,Fmin_cla_loc,Smin_loc,Smin_cla_loc
  real(double),dimension(:),allocatable, save:: Fmin_at,Fmin_cla_at,Smin_at,Smin_cla_at
  real(double):: thz_to_ev,temperature_to_ev
  real(double):: temp_min,temp_max,dtemp,temperature
  integer, parameter :: ntemp=1001
  contains

subroutine allocate_thermo()

  implicit none

  thz_to_ev=0.004135665538536d0
  temperature_to_ev=8.6173303E-05
  temp_min=1.d0
  temp_max=1000.d0
  dtemp=(temp_max-temp_min)/dble(ntemp-1)
  nsitedos_loc = nat
  if (allocated(Fmin))     deallocate(Fmin)              ; allocate(Fmin(ntemp))
  if (allocated(Fmin_cla)) deallocate(Fmin_cla)          ; allocate(Fmin_cla(ntemp))
  if (allocated(Smin))     deallocate(Smin)              ; allocate(Smin(ntemp))
  if (allocated(Smin_cla)) deallocate(Smin_cla)          ; allocate(Smin_cla(ntemp))

  if (allocated(Fmin_loc))     deallocate(Fmin_loc)           ; allocate(Fmin_loc(nsitedos, ntemp))
  if (allocated(Fmin_cla_loc)) deallocate(Fmin_cla_loc)       ; allocate(Fmin_cla_loc(nsitedos, ntemp))
  if (allocated(Smin_loc))     deallocate(Smin_loc)           ; allocate(Smin_loc(nsitedos, ntemp))
  if (allocated(Smin_cla_loc)) deallocate(Smin_cla_loc)       ; allocate(Smin_cla_loc(nsitedos, ntemp))


  if (allocated(Fmin_at))     deallocate(Fmin_at)         ; allocate(Fmin_at(nat))
  if (allocated(Fmin_cla_at)) deallocate(Fmin_cla_at)     ; allocate(Fmin_cla_at(nat))
  if (allocated(Smin_at))     deallocate(Smin_at)         ; allocate(Smin_at(nat))
  if (allocated(Smin_cla_at)) deallocate(Smin_cla_at)     ; allocate(Smin_cla_at(nat))




return
end subroutine allocate_thermo

end module thermo_module
