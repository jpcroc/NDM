module set_limits


contains

subroutine set_unlimit_for_atoms(rangml,im,i_start_at,i_final_at)
implicit none
integer, intent(in)  :: im, rangml
integer, intent(out) :: i_start_at,i_final_at

i_start_at=1
i_final_at=im

return
end subroutine set_unlimit_for_atoms

subroutine set_unlimit_for_cov (rangml, dim_train, i_final_cov,i_start_cov)
implicit none
integer, intent(in)  :: dim_train, rangml
integer, intent(out) :: i_final_cov, i_start_cov

i_start_cov=1
i_final_cov=dim_train

return
end subroutine set_unlimit_for_cov



subroutine set_limit_for_atoms(rangml,im,i_start_at,i_final_at)
#ifdef PARAML
     use mod_mpi_ml,ONLY:  nb_procsml
     implicit none
#else
     implicit none
#endif
     integer, intent(in) :: im           ! numbers of atoms to be distributed on procs
     integer, intent(inout)  :: rangml   ! rang of the proc for MPI
     integer, intent(out) :: i_start_at, i_final_at
#ifdef PARAML
     integer :: nratio1,nratio2, nrest
#else
     integer :: nb_procsml
#endif


#ifdef PARAML
if (im >= nb_procsml) then
  nrest=mod(im,nb_procsml)
  nratio1=im/nb_procsml+1
  nratio2=im/nb_procsml

  if (rangml<nrest) then
     i_start_at=rangml*nratio1+1
     if (rangml /= (nb_procsml-1) ) i_final_at=(rangml+1)*nratio1
  end if

  if (rangml>=nrest) then
    i_start_at = nrest*nratio1 + (rangml-nrest)*nratio2 + 1
    if (rangml /= (nb_procsml-1) ) i_final_at=nrest*nratio1 + (rangml-nrest+1)*nratio2
  endif

  if (rangml == (nb_procsml-1) ) i_final_at= im
else
  if (im == 1) then
     if (rangml==0) then
       i_start_at=1
       i_final_at=1
     else
       i_start_at=0
       i_final_at=0
     end if
  else
    if (rangml+1 <= im ) then
       i_start_at=rangml+1
       i_final_at=rangml+1
    else
       i_start_at=0
       i_final_at=0
    end if
  end if
end if
#else

rangml=0
nb_procsml=1
i_start_at=1
i_final_at=im

#endif



return
end subroutine set_limit_for_atoms

subroutine find_rang_of_my_atom (iconf)
! find the rang of the proc where my atom is located.
#ifdef PARAML
use mpi
use mod_mpi_ml
#endif
use derived_types, only: config_real
use ml_in_ndm_module, only: rangml, i_start_at, i_final_at
implicit none
integer, intent(in) :: iconf
integer :: i

if (allocated(config_real(iconf)%proc_atom)) deallocate(config_real(iconf)%proc_atom) ; allocate(config_real(iconf)%proc_atom(config_real(iconf)%nat) )

config_real(iconf)%proc_atom(:)=0
if ((i_start_at==0).and.(i_final_at==0)) return

do i=1,config_real(iconf)%nat
! on which proc is my atom
 if ( (i >= i_start_at).and.(i<=i_final_at) ) then
    config_real(iconf)%proc_atom(i) = rangml
 end if
end do

return
end subroutine find_rang_of_my_atom


subroutine  set_limit_for_cov (rangml, dim_train, i_final_cov,i_start_cov)
! set the limit for covariance matrix. all atoms are defined on all procs

#ifdef PARAML
     use mod_mpi_ml,ONLY:  nb_procsml
     implicit none
#else
     implicit none
#endif
     integer, intent(in) :: dim_train ! dimension of the covariance matrix to be distributed
     integer, intent(inout)  :: rangml  ! the rang of the proc for MPI
     integer, intent(out) :: i_final_cov,i_start_cov  !local limits defined on each procs.


#ifdef PARAML
     integer :: nratio1,nratio2, nrest
#else
     integer :: nb_procsml
#endif


#ifdef PARAML

nrest=mod(dim_train,nb_procsml)
nratio1=dim_train/nb_procsml+1
nratio2=dim_train/nb_procsml

 if (rangml<nrest) then
   i_start_cov=rangml*nratio1+1
   if (rangml /= (nb_procsml-1) ) i_final_cov=(rangml+1)*nratio1
 end if

 if (rangml>=nrest) then
 i_start_cov = nrest*nratio1 + (rangml-nrest)*nratio2 + 1
 if (rangml /= (nb_procsml-1) ) i_final_cov=nrest*nratio1 + (rangml-nrest+1)*nratio2
 endif

 if (rangml == (nb_procsml-1) ) i_final_cov= dim_train



!for the future we should pararelize also the computation of fingerprints
!over the atoms paralelization.


#else

rangml=0
nb_procsml=1
i_start_cov=1
i_final_cov=dim_train

#endif

return
end subroutine set_limit_for_cov

end module set_limits
