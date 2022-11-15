module descriptors_interface
    interface compute_descriptors
       subroutine compute_descriptors(xdesc_out,icount)
        implicit none
        double precision,dimension(:,:),allocatable :: xdesc_out
        integer, optional ::  icount
       end subroutine compute_descriptors
    end interface compute_descriptors
end module descriptors_interface



module compute_descriptors_mod
        implicit none
        contains


subroutine compute_descriptors(xdesc_out,icount)

#ifdef PARAML
use mpi
use mod_mpi_ml
#endif
use gen_com_m, ONLY: im,imm
use   ml_in_ndm_module, ONLY : rangml,debug,  &
                               descriptor_type, descriptor_g2, descriptor_g3, descriptor_behler, &
                               descriptor_afs, descriptor_soap, descriptor_pow_so3, descriptor_bispectrum_so3, &
                               descriptor_pow_so4,descriptor_bispectrum_so4, descriptor_g2_bispectrum_so4, descriptor_g2_afs, &
                               descriptor_mtp, &
                               target_energy, target_force,  &
                               i_start_at,i_final_at, &
                               g2,  g2_deriv,  g2_dim, &
                               g3,  g3_deriv,  g3_dim, &
                               afs, afs_deriv, afs_dim, &
                               pow_so3,pow_so3_deriv, pow_so3_dim,   &
                               bispectrum_so3,bispectrum_so3_deriv, bisso3_dim, &
                               pow_so4, pow_so4_deriv, pow_so4_dim, &
                               bispectrum_so4, bispectrum_so4_deriv, &
                               soap, soap_deriv, bisso4_dim, &
                               mtp, mtp_deriv, mtp_dim, soap_dim, &
                               typ,n2_typ,n3_typ, &
                               j_max, &
                               write_desc, imm_neigh, descriptor_g2_pow_so4, desc_forces
use temporary_data_cov, ONLY:  dim_xdesc
use set_limits
use angular_functions
use derived_types, only: config_desc, config_real
use compute_g2_mod
use compute_g3_mod
use compute_afs_mod
use compute_mtp_mod
use compute_bispectrum_so3_mod
use compute_bispectrum_so4_mod
use compute_pow_so4_mod
use compute_pow_so3_mod
use compute_soap_mod




implicit none


double precision, dimension(:,:),allocatable :: xdesc_out
integer, optional:: icount
!local
integer :: dim_reduce_full

!<begin  new version
double precision,dimension(:,:),allocatable :: local_g2
double precision,dimension(:,:,:,:),allocatable :: local_g2_deriv

double precision,dimension(:,:),allocatable :: local_g3
double precision,dimension(:,:,:,:),allocatable :: local_g3_deriv


double precision,dimension(:,:),allocatable :: local_pow_so3
double precision,dimension(:,:,:,:),allocatable :: local_pow_so3_deriv


double precision, dimension(:,:),allocatable :: local_bispectrum_so3
double precision, dimension(:,:,:,:),allocatable :: local_bispectrum_so3_deriv


double precision,dimension(:,:),allocatable   :: local_pow_so4
double precision,dimension(:,:,:,:),allocatable :: local_pow_so4_deriv

double complex,dimension(:,:),allocatable :: local_bispectrum_so4
double complex,dimension(:,:,:,:),allocatable :: local_bispectrum_so4_deriv



double precision,dimension(:,:),allocatable :: local_soap
double precision,dimension(:,:,:,:),allocatable :: local_soap_deriv


double precision,dimension(:,:),allocatable   :: local_mtp
double precision,dimension(:,:,:,:),allocatable :: local_mtp_deriv


double precision,dimension(:,:),allocatable   :: local_afs
double precision,dimension(:,:,:,:),allocatable :: local_afs_deriv

integer, dimension(imm)  :: d_n_neigh, l_d_n_neigh
integer, dimension(imm,imm_neigh) :: l_d_kind_neigh, d_kind_neigh


!<end  new version


 integer :: dim_reduce,dim_reduce1
 integer :: jso4


!interface

!subroutine compute_g2(i_start_at,i_final_at, l_d_n_neigh, l_d_kind_neigh,  local_g2,local_g2_deriv, iconf)
!    use gen_com_m, only: imm
!    use ml_in_ndm_module, only: g2_dim, imm_neigh
!    implicit none
!    integer, intent(in) :: i_start_at, i_final_at
!    integer, dimension(imm), intent(out) :: l_d_n_neigh
!    integer, dimension(imm, imm_neigh), intent(out) :: l_d_kind_neigh
!    double precision,dimension(g2_dim,imm), intent(out) :: local_g2
!    double precision,dimension(g2_dim,imm,0:imm_neigh,3), intent(out) :: local_g2_deriv
!    integer, optional, intent(in) :: iconf
!end subroutine compute_g2
!
!
!subroutine compute_g3(i_start_at,i_final_at, l_d_n_neigh, l_d_kind_neigh,  local_g3,local_g3_deriv, iconf)
!    use gen_com_m, only: imm
!    use ml_in_ndm_module, only: g3_dim, imm_neigh
!    implicit none
!    integer, intent(in) :: i_start_at, i_final_at
!    integer, dimension(imm), intent(out) :: l_d_n_neigh
!    integer, dimension(imm, imm_neigh), intent(out) :: l_d_kind_neigh
!    double precision,dimension(g3_dim,imm), intent(out) :: local_g3
!    double precision,dimension(g3_dim,imm,0:imm_neigh,3), intent(out) :: local_g3_deriv
!    integer, optional, intent(in) :: iconf
!end subroutine compute_g3
!
!subroutine compute_afs(i_start_at,i_final_at,l_d_n_neigh, l_d_kind_neigh, iconf)
!   use gen_com_m, only: imm
!   use ml_in_ndm_module, only: afs_dim, imm_neigh
!   implicit none
!   integer, intent (in) :: i_start_at,i_final_at
!   integer, dimension(imm),intent(out)  :: l_d_n_neigh
!   integer, dimension(imm,imm_neigh), intent(out) :: l_d_kind_neigh
!   !real(kind(0.d0)),dimension(afs_dim,imm),intent(out) :: local_afs_out
!   !real(kind(0.d0)),dimension(afs_dim,imm, 0:imm_neigh, 3),intent(out) :: local_afs_deriv_out
!   integer, optional :: iconf
!
!end subroutine compute_afs
!
!
!subroutine compute_pow_so3(i_start_at,i_final_at,l_d_n_neigh, l_d_kind_neigh, &
!                         local_pow_so3_out,local_pow_so3_deriv_out, iconf)
!   use gen_com_m, only: imm
!   use ml_in_ndm_module, only: pow_so3_dim, imm_neigh
!   implicit none
!   integer, intent (in) :: i_start_at,i_final_at
!   integer, dimension(imm),intent(out)  :: l_d_n_neigh
!   integer, dimension(imm,imm_neigh), intent(out) :: l_d_kind_neigh
!   real(kind(0.d0)),dimension(pow_so3_dim,imm),intent(out) :: local_pow_so3_out
!   real(kind(0.d0)),dimension(pow_so3_dim,imm, 0:imm_neigh, 3),intent(out) :: local_pow_so3_deriv_out
!   integer, optional :: iconf
!
!end subroutine compute_pow_so3
!
!
!subroutine compute_bispectrum_so3(i_start_at,i_final_at,l_d_n_neigh, l_d_kind_neigh, &
!                         local_bispectrum_so3_out,local_bispectrum_so3_deriv_out, iconf)
!   use gen_com_m, only: imm
!   use ml_in_ndm_module, only: bisso3_dim, imm_neigh
!   implicit none
!   integer, intent (in) :: i_start_at,i_final_at
!   integer, dimension(imm),intent(out)  :: l_d_n_neigh
!   integer, dimension(imm,imm_neigh), intent(out) :: l_d_kind_neigh
!   double precision, dimension(bisso3_dim,imm),intent(out) :: local_bispectrum_so3_out
!   double precision, dimension(bisso3_dim,imm, 0:imm_neigh, 3),intent(out) :: local_bispectrum_so3_deriv_out
!   integer, optional :: iconf
!
!end subroutine compute_bispectrum_so3
!
!
!subroutine compute_pow_so4(i_start_at,i_final_at,l_d_n_neigh, l_d_kind_neigh, &
!                           local_pow_so4_out, local_pow_so4_deriv_out, iconf)
!
!   use gen_com_m, only: imm
!   use ml_in_ndm_module, only: jj_max, imm_neigh
!   implicit none
!   integer, intent (in) :: i_start_at,i_final_at
!   integer, dimension(imm),intent(out)  :: l_d_n_neigh
!   integer, dimension(imm,imm_neigh), intent(out) :: l_d_kind_neigh
!   double precision,dimension(0:jj_max,imm),intent(out) :: local_pow_so4_out
!   double precision,dimension(0:jj_max,imm,0:imm_neigh,3), intent(out) :: local_pow_so4_deriv_out
!   integer, optional :: iconf
!
!end subroutine compute_pow_so4
!
!
!subroutine compute_bispectrum_so4(i_start_at,i_final_at,l_d_n_neigh, l_d_kind_neigh,iconf)
!   use gen_com_m, only: imm
!   use ml_in_ndm_module, only: bisso4_dim, imm_neigh
!   implicit none
!   integer, intent (in) :: i_start_at,i_final_at
!   integer, dimension(imm),intent(out)  :: l_d_n_neigh
!   integer, dimension(imm,imm_neigh), intent(out) :: l_d_kind_neigh
!   !double precision,dimension(bisso4_dim,imm),intent(out) :: local_bispectrum_so4_out
!   !double precision,dimension(bisso4_dim,imm, 0:imm_neigh, 3),intent(out) :: local_bispectrum_so4_deriv_out
!   integer, optional :: iconf
!
!end subroutine compute_bispectrum_so4
!
!
!
!subroutine compute_soap(i_start_at,i_final_at,l_d_n_neigh, l_d_kind_neigh, &
!                          iconf)
!   use gen_com_m, only: imm
!   use ml_in_ndm_module, only: soap_dim, imm_neigh
!   implicit none
!   integer, intent (in) :: i_start_at,i_final_at
!   integer, dimension(imm),intent(out)  :: l_d_n_neigh
!   integer, dimension(imm,imm_neigh), intent(out) :: l_d_kind_neigh
!   !double precision, dimension(soap_dim,imm),intent(out) :: local_soap_out
!   !double precision, dimension(soap_dim,imm, 0:imm_neigh, 3),intent(out) :: local_soap_deriv_out
!   integer, optional :: iconf
!
!end subroutine compute_soap
!
!
!
!
!subroutine compute_mtp(i_start_at,i_final_at,l_d_n_neigh, l_d_kind_neigh, &
!                         iconf)
!
!   use gen_com_m, only: imm
!   use ml_in_ndm_module, only: mtp_dim, imm_neigh
!   implicit none
!   integer, intent (in) :: i_start_at,i_final_at
!   integer, dimension(imm),intent(out)  :: l_d_n_neigh
!   integer, dimension(imm,imm_neigh), intent(out) :: l_d_kind_neigh
!   integer, optional :: iconf
!
!end subroutine compute_mtp
!
!end interface


!depending on rang the index of atoms is distributed on procs ...
call set_limit_for_atoms(rangml,im,i_start_at,i_final_at)
!find on which proc are is specific atom ...
call find_rang_of_my_atom(icount)

#ifdef PARAML
    call MPI_ALLREDUCE(MPI_IN_PLACE, config_real(icount)%proc_atom, config_real(icount)%nat,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD,codeml)
#endif
!debug write(6,*) 'rangml', rangml,   config_real(icount)%proc_atom
!without this MPI_BARRIER after ind_rang_of_my_atom the code crash ...
#ifdef PARAML
    call MPI_BARRIER(MPI_COMM_WORLD,codeml)
#endif

!C_DEBUG
call typ(n2_typ,n3_typ)


 if (debug) then
   write(6,'("ML: compute_descriptors in compute_descriptor rangml im i_start_at, i_final_at  nf ",i4, 4i9)')rangml, im, i_start_at, i_final_at, i_final_at-i_start_at+1
#ifdef PARAML
 call MPI_BARRIER(MPI_COMM_WORLD,codeml)
#endif
 end if

 select case(descriptor_type)
  case (descriptor_g2)
    if (allocated(local_g2))       deallocate(local_g2)      ; allocate(local_g2(g2_dim, imm))
    if (allocated(local_g2_deriv)) deallocate(local_g2_deriv); allocate(local_g2_deriv(g2_dim,imm, 0:imm_neigh, 3))

    call compute_g2(i_start_at,i_final_at,l_d_n_neigh, l_d_kind_neigh, local_g2,local_g2_deriv, icount)
#ifdef PARAML
    dim_reduce=g2_dim*imm
    call MPI_ALLREDUCE(local_g2,g2,dim_reduce,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,codeml)

    dim_reduce=g2_dim*imm*(imm_neigh+1)*3
    call MPI_ALLREDUCE(local_g2_deriv,g2_deriv,dim_reduce,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,codeml)
    call MPI_ALLREDUCE(l_d_n_neigh,d_n_neigh,imm,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD,codeml)
    call MPI_ALLREDUCE(l_d_kind_neigh,d_kind_neigh,imm*imm_neigh,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD,codeml)
#else
    g2(:,:) =local_g2(:,:)
    g2_deriv(:,:,:,:)=local_g2_deriv(:,:,:,:)
    d_n_neigh(1:imm) = l_d_n_neigh(1:imm)
    d_kind_neigh(1:imm, 1:imm_neigh) = l_d_kind_neigh(1:imm, 1:imm_neigh)
#endif

    dim_xdesc=g2_dim
    config_desc(icount)%dim_desc=dim_xdesc
    if (allocated(config_desc(icount)%n_neigh))    deallocate(config_desc(icount)%n_neigh)    ; allocate(config_desc(icount)%n_neigh(imm))
    if (allocated(config_desc(icount)%kind_neigh)) deallocate(config_desc(icount)%kind_neigh) ; allocate(config_desc(icount)%kind_neigh(imm,imm_neigh))
    if (allocated(config_desc(icount)%energy))     deallocate(config_desc(icount)%energy)     ; allocate(config_desc(icount)%energy(dim_xdesc,imm))
    if (allocated(config_desc(icount)%force ))     deallocate(config_desc(icount)%force )     ; allocate(config_desc(icount)%force (dim_xdesc,imm,0:imm_neigh,3))
    config_desc(icount)%n_neigh(1:imm)=d_n_neigh(1:imm)
    !write(*,*) 'outter', d_n_neigh(2)
    config_desc(icount)%kind_neigh(1:imm,1:imm_neigh)=d_kind_neigh(1:imm, 1:imm_neigh)
    config_desc(icount)%energy(1:dim_xdesc,1:imm) = g2(1:g2_dim,1:imm)
    config_desc(icount)%force(1:dim_xdesc,1:imm,0:imm_neigh,1:3)=g2_deriv(1:g2_dim,1:imm, 0:imm_neigh,1:3)

#ifdef PARAML
    call MPI_BARRIER(MPI_COMM_WORLD,codeml)
#endif

    if (write_desc) then
      call write_descriptors(icount)
    end if

  case (descriptor_g3)
    if (allocated(local_g3))       deallocate(local_g3)      ; allocate(local_g3(g3_dim, imm))
    if (allocated(local_g3_deriv)) deallocate(local_g3_deriv); allocate(local_g3_deriv(g3_dim, imm, 0:imm_neigh, 3))

   call compute_g3(i_start_at,i_final_at,l_d_n_neigh, l_d_kind_neigh,   local_g3,local_g3_deriv, icount)
#ifdef PARAML
    dim_reduce=g3_dim*imm
    call MPI_ALLREDUCE(local_g3,g3,dim_reduce,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,codeml)
    dim_reduce=g3_dim*imm*(imm_neigh+1)*3
    call MPI_ALLREDUCE(local_g3_deriv,g3_deriv,dim_reduce,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,codeml)

    call MPI_ALLREDUCE(l_d_n_neigh,d_n_neigh,imm,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD,codeml)
    call MPI_ALLREDUCE(l_d_kind_neigh,d_kind_neigh,imm*imm_neigh,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD,codeml)
#else
    g3(:,:) =local_g3(:,:)
    g3_deriv(:,:,:,:)=g3_deriv(:,:,:,:)
    d_n_neigh(1:imm) = l_d_n_neigh(1:imm)
    d_kind_neigh(1:imm, 1:imm_neigh) = l_d_kind_neigh(1:imm, 1:imm_neigh)
#endif
    dim_xdesc=g3_dim
    config_desc(icount)%dim_desc=dim_xdesc
    if (allocated(config_desc(icount)%n_neigh))    deallocate(config_desc(icount)%n_neigh)    ; allocate(config_desc(icount)%n_neigh(imm))
    if (allocated(config_desc(icount)%kind_neigh)) deallocate(config_desc(icount)%kind_neigh) ; allocate(config_desc(icount)%kind_neigh(imm,imm_neigh))
    if (allocated(config_desc(icount)%energy))     deallocate(config_desc(icount)%energy)     ; allocate(config_desc(icount)%energy(dim_xdesc,imm))
    if (allocated(config_desc(icount)%force ))     deallocate(config_desc(icount)%force )     ; allocate(config_desc(icount)%force (dim_xdesc,imm,0:imm_neigh,3))
    config_desc(icount)%n_neigh(1:imm)=d_n_neigh(1:imm)
    config_desc(icount)%kind_neigh(1:imm,1:imm_neigh)=d_kind_neigh(1:imm, 1:imm_neigh)
    config_desc(icount)%energy(1:dim_xdesc,1:imm) = g3(1:g3_dim,1:imm)
    config_desc(icount)%force(1:dim_xdesc,1:imm,0:imm_neigh,1:3)=g3_deriv(1:g3_dim,1:imm, 0:imm_neigh,1:3)
    if (write_desc) then
      call write_descriptors(icount)
    end if


  case (descriptor_behler)
    if (allocated(local_g2))       deallocate(local_g2)      ; allocate(local_g2(g2_dim, imm))
    if (allocated(local_g2_deriv)) deallocate(local_g2_deriv); allocate(local_g2_deriv(g2_dim,imm, 0:imm_neigh, 3))

    if (allocated(local_g3))       deallocate(local_g3)      ; allocate(local_g3(g3_dim, imm))
    if (allocated(local_g3_deriv)) deallocate(local_g3_deriv); allocate(local_g3_deriv(g3_dim, imm, 0:imm_neigh, 3))

    call compute_g3(i_start_at,i_final_at,l_d_n_neigh, l_d_kind_neigh, local_g3,local_g3_deriv, icount)
    call compute_g2(i_start_at,i_final_at,l_d_n_neigh, l_d_kind_neigh, local_g2,local_g2_deriv, icount)
#ifdef PARAML
    dim_reduce=g3_dim*imm
    call MPI_ALLREDUCE(local_g3,g3,dim_reduce,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,codeml)
    dim_reduce=g3_dim*imm*(imm_neigh+1)*3
    call MPI_ALLREDUCE(local_g3_deriv,g3_deriv,dim_reduce,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,codeml)

    dim_reduce1=g2_dim*imm
    call MPI_ALLREDUCE(local_g2,g2,dim_reduce1,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,codeml)
    dim_reduce1=g2_dim*imm*(imm_neigh+1)*3
    call MPI_ALLREDUCE(local_g2_deriv,g2_deriv,dim_reduce1,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,codeml)

    call MPI_ALLREDUCE(l_d_n_neigh,d_n_neigh,imm,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD,codeml)
    call MPI_ALLREDUCE(l_d_kind_neigh,d_kind_neigh,imm*imm_neigh,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD,codeml)

#else
    g2(:,:) =local_g2(:,:)
    g2_deriv(:,:,:,:)=local_g2_deriv(:,:,:,:)
    g3(:,:) =local_g3(:,:)
    g3_deriv(:,:,:,:)=local_g3_deriv(:,:,:,:)
    d_n_neigh(1:imm) = l_d_n_neigh(1:imm)
    d_kind_neigh(1:imm, 1:imm_neigh) = l_d_kind_neigh(1:imm, 1:imm_neigh)
#endif
    dim_xdesc=g2_dim+g3_dim
    config_desc(icount)%dim_desc=dim_xdesc
    if (allocated(config_desc(icount)%n_neigh))    deallocate(config_desc(icount)%n_neigh)    ; allocate(config_desc(icount)%n_neigh(imm))
    if (allocated(config_desc(icount)%kind_neigh)) deallocate(config_desc(icount)%kind_neigh) ; allocate(config_desc(icount)%kind_neigh(imm,imm_neigh))
    if (allocated(config_desc(icount)%energy))     deallocate(config_desc(icount)%energy)     ; allocate(config_desc(icount)%energy(dim_xdesc,imm))
    if (allocated(config_desc(icount)%force ))     deallocate(config_desc(icount)%force )     ; allocate(config_desc(icount)%force (dim_xdesc,imm,0:imm_neigh,3))
    config_desc(icount)%n_neigh(1:imm)=d_n_neigh(1:imm)
    config_desc(icount)%kind_neigh(1:imm,1:imm_neigh)=d_kind_neigh(1:imm, 1:imm_neigh)

    config_desc(icount)%energy(1:g2_dim,1:imm) = g2(1:g2_dim,1:imm)
    config_desc(icount)%force(1:g2_dim,1:imm,0:imm_neigh,1:3)=g2_deriv(1:g2_dim,1:imm, 0:imm_neigh,1:3)
    config_desc(icount)%energy(g2_dim+1:dim_xdesc,1:imm) = g3(1:g3_dim,1:imm)
    config_desc(icount)%force(g2_dim+1:dim_xdesc,1:imm,0:imm_neigh,1:3)=g3_deriv(1:g3_dim,1:imm, 0:imm_neigh,1:3)
    if (write_desc) then
      call write_descriptors(icount)
    end if
!

  case (descriptor_afs)
    !dim_xdesc=afs_dim
    !if (weighted) dim_xdesc=2*afs_dim
    config_desc(icount)%dim_desc=dim_xdesc

    !parallel computing of neighbours ...
    if (allocated(config_desc(icount)%n_neigh))    deallocate(config_desc(icount)%n_neigh)    ; allocate(config_desc(icount)%n_neigh(imm))
    if (allocated(config_desc(icount)%kind_neigh)) deallocate(config_desc(icount)%kind_neigh) ; allocate(config_desc(icount)%kind_neigh(imm,imm_neigh))
    !call pre_compute_neighbours_descriptors(i_start_at, i_final_at,config_desc(icount)%n_neigh, config_desc(icount)%kind_neigh,   icount)

#ifdef PARAML
    !put the list of neighbours in descritors together ...
    !call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(icount)%n_neigh,imm,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD,codeml)
    !call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(icount)%kind_neigh ,imm*imm_neigh,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD,codeml)
    !call MPI_BARRIER(MPI_COMM_WORLD,codeml)
!
#endif

!     if (allocated(config_desc(icount)%n_neigh_ghost)     ) deallocate(config_desc(icount)%n_neigh_ghost )  ; allocate( config_desc(icount)%n_neigh_ghost(imm) )
!     if (allocated(config_desc(icount)%kind_neigh_ghost)) deallocate(config_desc(icount)%kind_neigh_ghost)  ; allocate( config_desc(icount)%kind_neigh_ghost(imm, imm_neigh) )
!     if (allocated(config_desc(icount)%kind_neigh_proc) ) deallocate(config_desc(icount)%kind_neigh_proc )  ; allocate( config_desc(icount)%kind_neigh_proc(imm, imm_neigh) )
!     call snap_pre_pack_force_descriptor(icount)
!     call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(icount)%n_neigh_ghost   ,imm          ,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD,codeml)
!     call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(icount)%kind_neigh_ghost,imm*imm_neigh,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD,codeml)
!     call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(icount)%kind_neigh_proc ,imm*imm_neigh,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD,codeml)
!old-v call compute_afs(i_start_at,i_final_at,config_desc(icount)%n_neigh, config_desc(icount)%kind_neigh, config_desc(icount)%energy,   config_desc(icount)%force,   icount)

    if (allocated(config_desc(icount)%energy))     deallocate(config_desc(icount)%energy)     ; allocate(config_desc(icount)%energy(dim_xdesc,imm))
    if (desc_forces) then
      if (.not.((i_start_at==0).and.(i_final_at==0))) then
        if (allocated(config_desc(icount)%force ))     deallocate(config_desc(icount)%force )     ; allocate(config_desc(icount)%force (dim_xdesc,i_start_at:i_final_at,0:imm_neigh,3))
      end if
    end if

    call compute_afs(i_start_at,i_final_at, config_desc(icount)%n_neigh, config_desc(icount)%kind_neigh,   icount)
    !all the neighbours
    call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(icount)%n_neigh,imm,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD,codeml)
    call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(icount)%kind_neigh ,imm*imm_neigh,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD,codeml)

    !debug call distribute_ghost_descriptors()
#ifdef PARAML
    dim_reduce=dim_xdesc*imm
    call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(icount)%energy ,dim_reduce,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,codeml)

!    dim_reduce=afs_dim*imm*(imm_neigh+1)*3
!    if (desc_forces) call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(icount)%force ,dim_reduce,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,codeml)

    !call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(icount)%n_neigh,imm,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD,codeml)
    !call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(icount)%kind_neigh ,imm*imm_neigh,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD,codeml)
    call MPI_BARRIER(MPI_COMM_WORLD,codeml)
#endif

    if (write_desc) then
      call write_descriptors(icount)
    end if


  case (descriptor_soap)

    dim_xdesc=soap_dim
    config_desc(icount)%dim_desc=dim_xdesc

    if (allocated(config_desc(icount)%n_neigh))    deallocate(config_desc(icount)%n_neigh) ; allocate(config_desc(icount)%n_neigh(imm))
    if (allocated(config_desc(icount)%kind_neigh)) deallocate(config_desc(icount)%kind_neigh) ; allocate(config_desc(icount)%kind_neigh(imm,imm_neigh))
    if (allocated(config_desc(icount)%energy))     deallocate(config_desc(icount)%energy) ; allocate(config_desc(icount)%energy(dim_xdesc,imm))
    if (desc_forces) then
       if (allocated(config_desc(icount)%force ))     deallocate(config_desc(icount)%force ) ; allocate(config_desc(icount)%force (dim_xdesc,imm,0:imm_neigh,3))
    else
       if (allocated(config_desc(icount)%force ))     deallocate(config_desc(icount)%force ) ; allocate(config_desc(icount)%force (dim_xdesc,1,0:imm_neigh,3))
    end if

    !call compute_soap(i_start_at,i_final_at,config_desc(icount)%n_neigh, config_desc(icount)%kind_neigh, config_desc(icount)%energy,   config_desc(icount)%force,   icount)
    call compute_soap(i_start_at,i_final_at,config_desc(icount)%n_neigh, config_desc(icount)%kind_neigh,   icount)

#ifdef PARAML
    dim_reduce=imm*soap_dim
    call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(icount)%energy ,dim_reduce,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,codeml)

    dim_reduce=imm*soap_dim*(imm_neigh+1)*3
    if (desc_forces) call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(icount)%force ,dim_reduce,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,codeml)


    call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(icount)%n_neigh,imm,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD,codeml)
    call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(icount)%kind_neigh ,imm*imm_neigh,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD,codeml)
    call MPI_BARRIER(MPI_COMM_WORLD,codeml)
#endif

    if (write_desc) then
      call write_descriptors(icount)
      !call MPI_FINALIZE(codeml)
      !stop

    endif


  case (descriptor_pow_so3)
    if (allocated(local_pow_so3))       deallocate(local_pow_so3)      ; allocate(local_pow_so3(pow_so3_dim, imm))
    if (allocated(local_pow_so3_deriv)) deallocate(local_pow_so3_deriv); allocate(local_pow_so3_deriv(pow_so3_dim,imm, 0:imm_neigh, 3))

    call compute_pow_so3(i_start_at,i_final_at,l_d_n_neigh, l_d_kind_neigh, local_pow_so3,local_pow_so3_deriv, icount)
#ifdef PARAML
    dim_reduce=(pow_so3_dim)*imm
    call MPI_ALLREDUCE(local_pow_so3,pow_so3,dim_reduce,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,codeml)

    dim_reduce=(pow_so3_dim)*imm*(imm_neigh+1)*3
    if (desc_forces) call MPI_ALLREDUCE(local_pow_so3_deriv,pow_so3_deriv,dim_reduce,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,codeml)
    call MPI_ALLREDUCE(l_d_n_neigh,d_n_neigh,imm,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD,codeml)
    call MPI_ALLREDUCE(l_d_kind_neigh,d_kind_neigh,imm*imm_neigh,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD,codeml)
#else
    pow_so3(:,:) =local_pow_so3(:,:)
    pow_so3_deriv(:,:,:,:)=local_pow_so3_deriv(:,:,:,:)
    d_n_neigh(1:imm) = l_d_n_neigh(1:imm)
    d_kind_neigh(1:imm, 1:imm_neigh) = l_d_kind_neigh(1:imm, 1:imm_neigh)
#endif

    dim_xdesc=pow_so3_dim
    config_desc(icount)%dim_desc=dim_xdesc
    if (allocated(config_desc(icount)%n_neigh))    deallocate(config_desc(icount)%n_neigh)    ; allocate(config_desc(icount)%n_neigh(imm))
    if (allocated(config_desc(icount)%kind_neigh)) deallocate(config_desc(icount)%kind_neigh) ; allocate(config_desc(icount)%kind_neigh(imm,imm_neigh))
    if (allocated(config_desc(icount)%energy))     deallocate(config_desc(icount)%energy)     ; allocate(config_desc(icount)%energy(dim_xdesc,imm))
    if (allocated(config_desc(icount)%force ))     deallocate(config_desc(icount)%force )     ; allocate(config_desc(icount)%force (dim_xdesc,imm,0:imm_neigh,3))
    config_desc(icount)%n_neigh(1:imm)=d_n_neigh(1:imm)
    !write(*,*) 'outter', d_n_neigh(2)
    config_desc(icount)%kind_neigh(1:imm,1:imm_neigh)=d_kind_neigh(1:imm, 1:imm_neigh)
    config_desc(icount)%energy(1:dim_xdesc,1:imm) = pow_so3(1:pow_so3_dim,1:imm)
    config_desc(icount)%force(1:dim_xdesc,1:imm,0:imm_neigh,1:3)=pow_so3_deriv(1:pow_so3_dim,1:imm, 0:imm_neigh,1:3)

#ifdef PARAML
    call MPI_BARRIER(MPI_COMM_WORLD,codeml)
#endif

    if (write_desc) then
      call write_descriptors(icount)
    end if


  case (descriptor_bispectrum_so3)

    dim_xdesc=bisso3_dim
    config_desc(icount)%dim_desc=dim_xdesc
    if (allocated(config_desc(icount)%n_neigh))    deallocate(config_desc(icount)%n_neigh) ; allocate(config_desc(icount)%n_neigh(imm))
    if (allocated(config_desc(icount)%kind_neigh)) deallocate(config_desc(icount)%kind_neigh) ; allocate(config_desc(icount)%kind_neigh(imm,imm_neigh))
    if (allocated(config_desc(icount)%energy))     deallocate(config_desc(icount)%energy) ; allocate(config_desc(icount)%energy(dim_xdesc,imm))
    if (allocated(config_desc(icount)%force ))     deallocate(config_desc(icount)%force ) ; allocate(config_desc(icount)%force (dim_xdesc,imm,0:imm_neigh,3))

    call compute_bispectrum_so3(i_start_at,i_final_at,config_desc(icount)%n_neigh, config_desc(icount)%kind_neigh, config_desc(icount)%energy,   config_desc(icount)%force,   icount)
#ifdef PARAML
    dim_reduce=imm*bisso3_dim
    call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(icount)%energy ,dim_reduce,MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD,codeml)
    dim_reduce=imm*bisso3_dim*(imm_neigh+1)*3
    if (desc_forces) call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(icount)%force ,dim_reduce,MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD,codeml)

    call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(icount)%n_neigh,imm,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD,codeml)
    call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(icount)%kind_neigh ,imm*imm_neigh,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD,codeml)
#endif

    if (write_desc) then
      call write_descriptors(icount)
    endif



  case (descriptor_g2_pow_so4)

    jso4=int(2*j_max)
    dim_xdesc=pow_so4_dim + g2_dim
    config_desc(icount)%dim_desc=dim_xdesc
    if (allocated(config_desc(icount)%energy)) deallocate(config_desc(icount)%energy) ; allocate(config_desc(icount)%energy(dim_xdesc,imm))
    if (allocated(config_desc(icount)%force )) deallocate(config_desc(icount)%force ) ; allocate(config_desc(icount)%force (dim_xdesc,imm,0:imm_neigh,3))
    if (allocated(config_desc(icount)%n_neigh)) deallocate(config_desc(icount)%n_neigh) ; allocate(config_desc(icount)%n_neigh(imm))
    if (allocated(config_desc(icount)%kind_neigh)) deallocate(config_desc(icount)%kind_neigh) ; allocate(config_desc(icount)%kind_neigh(imm,imm_neigh))


    call compute_g2(i_start_at,i_final_at, config_desc(icount)%n_neigh, config_desc(icount)%kind_neigh, g2, g2_deriv, icount)
    call compute_pow_so4(i_start_at,i_final_at, config_desc(icount)%n_neigh, config_desc(icount)%kind_neigh, pow_so4, pow_so4_deriv, icount)
#ifdef PARAML
    dim_reduce=(jso4+1)*imm
    call MPI_ALLREDUCE(MPI_IN_PLACE,pow_so4,dim_reduce,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,codeml)

    dim_reduce_full=(jso4+1)*imm*(imm_neigh+1)*3
    if (desc_forces) call MPI_ALLREDUCE(MPI_IN_PLACE,pow_so4_deriv,dim_reduce_full,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,codeml)

    dim_reduce=g2_dim*imm
    call MPI_ALLREDUCE(MPI_IN_PLACE,g2,dim_reduce,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,codeml)

    dim_reduce=g2_dim*imm*(imm_neigh+1)*3
    if (desc_forces) call MPI_ALLREDUCE(MPI_IN_PLACE,g2_deriv,dim_reduce,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,codeml)

    call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(icount)%n_neigh,imm,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD,codeml)
    call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(icount)%kind_neigh ,imm*imm_neigh,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD,codeml)

#endif
    dim_xdesc=pow_so4_dim + g2_dim
    config_desc(icount)%dim_desc=dim_xdesc

    config_desc(icount)%energy(1:pow_so4_dim,1:imm) = pow_so4(0:jso4,1:imm)
    config_desc(icount)%force(1:pow_so4_dim,1:imm,0:imm_neigh,1:3)=pow_so4_deriv(0:jso4,1:imm, 0:imm_neigh,1:3)

    if (desc_forces) config_desc(icount)%energy(pow_so4_dim+1:dim_xdesc,1:imm) = g2(1:g2_dim,1:imm)
    if (desc_forces) config_desc(icount)%force(pow_so4_dim+1:dim_xdesc,1:imm,0:imm_neigh,1:3)=g2_deriv(1:g2_dim,1:imm, 0:imm_neigh,1:3)

    if (write_desc) then
      call write_descriptors(icount)
    endif


  case (descriptor_pow_so4)
    jso4=int(2*j_max)
    if (allocated(local_pow_so4)) deallocate(local_pow_so4); allocate(local_pow_so4(0:jso4,imm))
    if (allocated(local_pow_so4_deriv)) deallocate(local_pow_so4_deriv); allocate(local_pow_so4_deriv(0:jso4,imm,0:imm_neigh,3))
    call compute_pow_so4(i_start_at,i_final_at,l_d_n_neigh, l_d_kind_neigh, &
                         local_pow_so4, local_pow_so4_deriv, icount)
#ifdef PARAML
    dim_reduce=(jso4+1)*imm
    dim_reduce_full=(jso4+1)*imm*(imm_neigh+1)*3
    call MPI_ALLREDUCE(local_pow_so4,pow_so4,dim_reduce,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,codeml)
    call MPI_ALLREDUCE(local_pow_so4_deriv,pow_so4_deriv,dim_reduce_full,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,codeml)
    call MPI_ALLREDUCE(l_d_n_neigh,d_n_neigh,imm,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD,codeml)
    call MPI_ALLREDUCE(l_d_kind_neigh,d_kind_neigh,imm*imm_neigh,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD,codeml)

#else
    pow_so4(:,:) =local_pow_so4(:,:)
    pow_so4_deriv(:,:,:,:)=local_pow_so4_deriv(:,:,:,:)
    d_n_neigh(1:imm) = l_d_n_neigh(1:imm)
    d_kind_neigh(1:imm, 1:imm_neigh) = l_d_kind_neigh(1:imm, 1:imm_neigh)
#endif
    dim_xdesc=(jso4+1)
    config_desc(icount)%dim_desc=dim_xdesc
    if (allocated(config_desc(icount)%n_neigh)) deallocate(config_desc(icount)%n_neigh) ; allocate(config_desc(icount)%n_neigh(imm))
    if (allocated(config_desc(icount)%kind_neigh)) deallocate(config_desc(icount)%kind_neigh) ; allocate(config_desc(icount)%kind_neigh(imm,imm_neigh))
    if (allocated(config_desc(icount)%energy)) deallocate(config_desc(icount)%energy) ; allocate(config_desc(icount)%energy(dim_xdesc,imm))
    if (allocated(config_desc(icount)%force )) deallocate(config_desc(icount)%force ) ; allocate(config_desc(icount)%force (dim_xdesc,imm,0:imm_neigh,3))
    config_desc(icount)%n_neigh(1:imm)=d_n_neigh(1:imm)
    config_desc(icount)%kind_neigh(1:imm,1:imm_neigh)=d_kind_neigh(1:imm, 1:imm_neigh)
    config_desc(icount)%energy(1:dim_xdesc,1:imm) = pow_so4(0:jso4,1:imm)
    config_desc(icount)%force(1:dim_xdesc,1:imm,0:imm_neigh,1:3)=pow_so4_deriv(0:jso4,1:imm, 0:imm_neigh,1:3)
    if (write_desc) then
      call write_descriptors(icount)
    endif

  case (descriptor_bispectrum_so4)
    !dim_xdesc=bisso4_dim
    config_desc(icount)%dim_desc=dim_xdesc
    if (allocated(config_desc(icount)%n_neigh))    deallocate(config_desc(icount)%n_neigh) ; allocate(config_desc(icount)%n_neigh(imm))
    if (allocated(config_desc(icount)%kind_neigh)) deallocate(config_desc(icount)%kind_neigh) ; allocate(config_desc(icount)%kind_neigh(imm,imm_neigh))

    if (allocated(config_desc(icount)%energy))     deallocate(config_desc(icount)%energy) ; allocate(config_desc(icount)%energy(dim_xdesc,imm))
    if (desc_forces) then
      if (.not.((i_start_at==0).and.(i_final_at==0))) then
        if (allocated(config_desc(icount)%force ))     deallocate(config_desc(icount)%force ) ; allocate(config_desc(icount)%force (dim_xdesc,i_start_at:i_final_at,0:imm_neigh,3))
      end if
    end if

    call compute_bispectrum_so4(i_start_at,i_final_at,config_desc(icount)%n_neigh, config_desc(icount)%kind_neigh, icount)
#ifdef PARAML
    !all the neighbours
    call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(icount)%n_neigh,imm,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD,codeml)
    call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(icount)%kind_neigh ,imm*imm_neigh,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD,codeml)
    !dim_reduce=imm*bisso4_dim*(imm_neigh+1)*3
    !if (desc_forces) call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(icount)%force, dim_reduce,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,codeml)
    dim_reduce=imm*dim_xdesc
    call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(icount)%energy,  dim_reduce,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,codeml)
#endif
    if (desc_forces) then
      if (.not.((i_start_at==0).and.(i_final_at==0))) then
       config_desc(icount)%force(1:dim_xdesc,i_start_at:i_final_at,0:imm_neigh,1:3)=- config_desc(icount)%force(1:dim_xdesc,i_start_at:i_final_at,0:imm_neigh,1:3)
      end if
    end if
    if (write_desc) then
      call write_descriptors(icount)
    endif


 case (descriptor_mtp)

    !C dim_xdesc=mtp_dim
    config_desc(icount)%dim_desc=dim_xdesc
    if (allocated(config_desc(icount)%n_neigh))    deallocate(config_desc(icount)%n_neigh) ; allocate(config_desc(icount)%n_neigh(imm))
    if (allocated(config_desc(icount)%kind_neigh)) deallocate(config_desc(icount)%kind_neigh) ; allocate(config_desc(icount)%kind_neigh(imm,imm_neigh))
    if (allocated(config_desc(icount)%energy))     deallocate(config_desc(icount)%energy) ; allocate(config_desc(icount)%energy(dim_xdesc,imm))


    if (desc_forces) then
      !LM
      if (.not.((i_start_at==0).and.(i_final_at==0))) then
        if (allocated(config_desc(icount)%force ))     deallocate(config_desc(icount)%force ) ; allocate(config_desc(icount)%force (dim_xdesc,imm,0:imm_neigh,3))
     !LM
     end if
    end if

    call compute_mtp(i_start_at,i_final_at,config_desc(icount)%n_neigh, config_desc(icount)%kind_neigh,  icount)
#ifdef PARAML
    !all the neighbours
    call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(icount)%n_neigh,imm,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD,codeml)
    call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(icount)%kind_neigh ,imm*imm_neigh,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD,codeml)
    !dim_reduce=imm*dim_xdesc*(imm_neigh+1)*3
    !if (desc_forces) call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(icount)%force ,dim_reduce,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,codeml)
    dim_reduce=imm*dim_xdesc
    call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(icount)%energy ,dim_reduce,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,codeml)
#endif

    if (write_desc) then
      call write_descriptors(icount)
    endif
#ifdef PARAML
    call MPI_BARRIER(MPI_COMM_WORLD,codeml)
#endif


  case (descriptor_g2_bispectrum_so4)

    dim_xdesc=bisso4_dim + g2_dim
    if (allocated(bispectrum_so4))       deallocate(bispectrum_so4)      ; allocate(bispectrum_so4(bisso4_dim, imm))
    if (allocated(bispectrum_so4_deriv)) deallocate(bispectrum_so4_deriv); allocate(bispectrum_so4_deriv(bisso4_dim, imm, 0:imm_neigh,3))
    if (allocated(g2))       deallocate(g2)       ; allocate(g2(g2_dim, imm))
    if (allocated(g2_deriv)) deallocate(g2_deriv) ; allocate(g2_deriv(g2_dim, imm, 0:imm_neigh, 3))


    if (allocated(config_desc(icount)%n_neigh))    deallocate(config_desc(icount)%n_neigh) ; allocate(config_desc(icount)%n_neigh(imm))
    if (allocated(config_desc(icount)%kind_neigh)) deallocate(config_desc(icount)%kind_neigh) ; allocate(config_desc(icount)%kind_neigh(imm,imm_neigh))

    config_desc(icount)%dim_desc=bisso4_dim
    if (allocated(config_desc(icount)%energy))     deallocate(config_desc(icount)%energy) ; allocate(config_desc(icount)%energy(bisso4_dim,imm))
    if (allocated(config_desc(icount)%force ))     deallocate(config_desc(icount)%force ) ; allocate(config_desc(icount)%force (bisso4_dim,imm,0:imm_neigh,3))
    call compute_bispectrum_so4(i_start_at,i_final_at,config_desc(icount)%n_neigh, config_desc(icount)%kind_neigh, icount)
    bispectrum_so4(:,:)=config_desc(icount)%energy(:,:)
    bispectrum_so4_deriv(:,:,:,:) = config_desc(icount)%force(:,:,:,:)

    call compute_g2(i_start_at,i_final_at, config_desc(icount)%n_neigh, config_desc(icount)%kind_neigh,  g2,   g2_deriv, icount)


#ifdef PARAML
    dim_reduce=imm*bisso4_dim
    call MPI_ALLREDUCE(MPI_IN_PLACE, bispectrum_so4, dim_reduce, MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,codeml)

    dim_reduce=imm*bisso4_dim*(imm_neigh+1)*3
    if (desc_forces) call MPI_ALLREDUCE(MPI_IN_PLACE,bispectrum_so4_deriv,dim_reduce,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,codeml)

    dim_reduce=g2_dim*imm
    call MPI_ALLREDUCE(MPI_IN_PLACE,g2,dim_reduce,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,codeml)

    dim_reduce=g2_dim*imm*(imm_neigh+1)*3
    if (desc_forces) call MPI_ALLREDUCE(MPI_IN_PLACE,g2_deriv,dim_reduce,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,codeml)


    call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(icount)%n_neigh,imm,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD,codeml)
    call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(icount)%kind_neigh ,imm*imm_neigh,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD,codeml)

#endif
    config_desc(icount)%dim_desc=dim_xdesc
    if (allocated(config_desc(icount)%energy))     deallocate(config_desc(icount)%energy) ; allocate(config_desc(icount)%energy(dim_xdesc,imm))
    if (allocated(config_desc(icount)%force ))     deallocate(config_desc(icount)%force ) ; allocate(config_desc(icount)%force (dim_xdesc,imm,0:imm_neigh,3))

    config_desc(icount)%energy(1: bisso4_dim,1:imm) = bispectrum_so4(1:bisso4_dim,1:imm)
    config_desc(icount)%energy(bisso4_dim+1:dim_xdesc,1:imm) = g2(1:g2_dim,1:imm)
    !warning: the derivatives in bispectrum have opposite sign
    if (desc_forces) config_desc(icount)%force(1:bisso4_dim,1:imm,0:imm_neigh,1:3)=-bispectrum_so4_deriv(1:bisso4_dim,1:imm, 0:imm_neigh,1:3)
    if (desc_forces) config_desc(icount)%force(bisso4_dim+1:dim_xdesc,1:imm,0:imm_neigh,1:3)=g2_deriv(1:g2_dim,1:imm, 0:imm_neigh,1:3)

    if (write_desc) then
      call write_descriptors(icount)
    endif

  case (descriptor_g2_afs)

    dim_xdesc=afs_dim + g2_dim
    if (allocated(g2))       deallocate(g2)         ; allocate(g2(g2_dim, imm))
    if (allocated(g2_deriv)) deallocate(g2_deriv)   ; allocate(g2_deriv(g2_dim, imm, 0:imm_neigh, 3))
    if (allocated(afs)) deallocate(afs)             ; allocate(afs(afs_dim,imm))
    if (allocated(afs_deriv)) deallocate(afs_deriv) ; allocate(afs_deriv(afs_dim,imm,0:imm_neigh, 3))

    if (allocated(config_desc(icount)%n_neigh))    deallocate(config_desc(icount)%n_neigh) ; allocate(config_desc(icount)%n_neigh(imm))
    if (allocated(config_desc(icount)%kind_neigh)) deallocate(config_desc(icount)%kind_neigh) ; allocate(config_desc(icount)%kind_neigh(imm,imm_neigh))

    config_desc(icount)%dim_desc=afs_dim
    if (allocated(config_desc(icount)%energy))     deallocate(config_desc(icount)%energy) ; allocate(config_desc(icount)%energy(afs_dim,imm))
    if (allocated(config_desc(icount)%force ))     deallocate(config_desc(icount)%force ) ; allocate(config_desc(icount)%force (afs_dim,imm,0:imm_neigh,3))
    call compute_afs(i_start_at, i_final_at, config_desc(icount)%n_neigh, config_desc(icount)%kind_neigh, icount)
    afs(:,:) = config_desc(icount)%energy(:,:)
    afs_deriv(:,:,:,:) = config_desc(icount)%force(:,:,:,:)



    call compute_g2(i_start_at , i_final_at, config_desc(icount)%n_neigh, config_desc(icount)%kind_neigh,  g2,    g2_deriv, icount)
    !g2(:,:) = config_desc(icount)%energy(:,:)
    !g2_deriv(:,:,:,:) = config_desc(icount)%force(:,:,:,:)


    config_desc(icount)%dim_desc=dim_xdesc
    if (allocated(config_desc(icount)%energy))     deallocate(config_desc(icount)%energy) ; allocate(config_desc(icount)%energy(dim_xdesc,imm))
    if (allocated(config_desc(icount)%force ))     deallocate(config_desc(icount)%force ) ; allocate(config_desc(icount)%force (dim_xdesc,imm,0:imm_neigh,3))
#ifdef PARAML
    dim_reduce=imm*afs_dim
    call MPI_BARRIER(MPI_COMM_WORLD,codeml)
    call MPI_ALLREDUCE(MPI_IN_PLACE,afs,dim_reduce,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,codeml)

    dim_reduce=imm*afs_dim*(imm_neigh+1)*3
    call MPI_ALLREDUCE(MPI_IN_PLACE,afs_deriv,dim_reduce,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,codeml)

    dim_reduce=g2_dim*imm
    call MPI_ALLREDUCE(MPI_IN_PLACE,g2,dim_reduce,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,codeml)

    dim_reduce=g2_dim*imm*(imm_neigh+1)*3
    if (desc_forces) call MPI_ALLREDUCE(MPI_IN_PLACE,g2_deriv,dim_reduce,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,codeml)


    call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(icount)%n_neigh,imm,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD,codeml)
    call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(icount)%kind_neigh ,imm*imm_neigh,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD,codeml)
#endif

    config_desc(icount)%energy(1:afs_dim,1:imm) = afs(1:afs_dim,1:imm)
    config_desc(icount)%energy(afs_dim+1:dim_xdesc,1:imm) = g2(1:g2_dim,1:imm)
    !warning: the derivatives in bispectrum have opposite sign
    if (desc_forces) config_desc(icount)%force(1:afs_dim,1:imm,0:imm_neigh,1:3)=afs_deriv(1:afs_dim,1:imm, 0:imm_neigh,1:3)
    if (desc_forces) config_desc(icount)%force(afs_dim+1:dim_xdesc,1:imm,0:imm_neigh,1:3)=g2_deriv(1:g2_dim,1:imm, 0:imm_neigh,1:3)

    if (write_desc) then
      call write_descriptors(icount)
    endif

  case default
       if (rangml==0) write(6,*) 'ML: No implementation for descriptor_type...', descriptor_type
       stop "fatal in File: descriptor.f90, Subroutine: compute_descriptor"
end select


#ifdef PARAML
if (debug) then
 call MPI_BARRIER(MPI_COMM_WORLD,codeml)
  if (rangml==0) write(6,'("ML: descriptor was computed in compute_descriptors...")')
end if
 call MPI_BARRIER(MPI_COMM_WORLD,codeml)
#endif
!mpi_rangml
if (rangml==0) call renormalize_descritors_iconf(icount)

return
end subroutine compute_descriptors

subroutine renormalize_descritors_iconf(iconf)
use derived_types, only : config_desc
use ml_in_ndm_module, only: val_desc_max
implicit none
integer, intent(in) :: iconf

config_desc(iconf)%energy(:,:) = config_desc(iconf)%energy(:,:)/val_desc_max

return
end subroutine renormalize_descritors_iconf

subroutine val_renormalize_descritors(iconf)
use derived_types, only : config_desc
use ml_in_ndm_module, only: tmp_val_desc_max
implicit none
integer, intent(in) :: iconf

tmp_val_desc_max=MAXVAL(dabs(config_desc(iconf)%energy(:,:)))

return
end subroutine val_renormalize_descritors


!$--------------------------------------------------------------
subroutine  write_descriptors (iconf)
!$--------------------------------------------------------------
#ifdef PARAML
use mpi
use mod_mpi_ml
#endif
use ml_in_ndm_module, only: descriptor_type, descriptor_afs, descriptor_bispectrum_so4, rangml
implicit none
integer, intent(in) :: iconf
integer :: ip

if ((descriptor_type==descriptor_afs).or.(descriptor_type==descriptor_bispectrum_so4)) then
    !do  ip=1,nb_procsml
      !if (rangml==(ip-1)) call para_write_descriptors(ip-1, iconf)
      call para_write_descriptors(0, iconf)
    !end do
  else
    call serial_write_descriptors(iconf)
end if


return
end subroutine write_descriptors
!<-------------------------------------------------------------



subroutine serial_write_descriptors (iconf)
use ml_in_ndm_module, only: db_path, desc_forces, rangml
use derived_types, only : config_desc, config_real
implicit none
integer, intent(in) :: iconf
integer :: ndim
character(len=100) :: efilename, ffilename
character (len=60) :: CHFMT, CHFMTf
integer :: eunit, funit, ia,ja,ix

ndim = config_desc(iconf)%dim_desc
efilename='desc'//trim(adjustl(db_path))//config_real(iconf)%class//'_'//config_real(iconf)%klm//'_'//config_real(iconf)%cnumber//'.eml'
if (desc_forces) ffilename='desc'//trim(adjustl(db_path))//config_real(iconf)%class//'_'//config_real(iconf)%klm//'_'//config_real(iconf)%cnumber//'.fml'

!there is only one procs who writting ...
if (rangml==0) then

eunit=41
funit=42


open(eunit,file=efilename,status='unknown')
if (desc_forces) open(funit,file=ffilename,status='unknown')
write(CHFMT,*)'(i6, 1x, ',int(ndim),'e20.10)'
!write(CHFMTf,*)'(i6, 1x, e20.10,',int(ndim),'e20.10)'
write(CHFMTf,*)'(i6, 1x,i5,1x,i6, ',int(ndim),'e20.10)'
do ia=1,config_real(iconf)%nat
    write(eunit,FMT=CHFMT) ia, real(config_desc(iconf)%energy(:,ia))
    if (desc_forces) then
      write(funit,'(i6)') ia
      do ix=1,3
       write(funit,FMT=CHFMTf) ia, ix, 0,  real(config_desc(iconf)%force(:,ia,0,ix))
      end do
    end if

    if (desc_forces) then
    do ja=1,config_desc(iconf)%n_neigh(ia)
      do ix=1,3
         !write(funit,FMT=CHFMTf) config_desc(iconf)%kind_neigh(ia,ja), config_real(iconf)%u_ij(ia,ja,ix),  config_desc(iconf)%force(:,ia,ja,ix)
         write(funit,FMT=CHFMTf) config_desc(iconf)%kind_neigh(ia,ja), ix, ja,  config_desc(iconf)%force(:,ia,ja,ix)
      end do
    end do
    end if
enddo

close(eunit)
if (desc_forces) close(funit)
end if


return
end subroutine serial_write_descriptors


subroutine para_write_descriptors (iproc, iconf)
use ml_in_ndm_module, only: db_path, desc_forces, rangml, i_start_at, i_final_at, debug
use derived_types, only : config_desc, config_real
implicit none
integer, intent(in) :: iconf, iproc
integer :: ndim
character(len=100) :: efilename, ffilename
character (len=60) :: CHFMT, CHFMTf
integer :: eunit, funit, ia,ja,ix


if (debug) then
  if (rangml==0) write(6,*) 'ML: writting descritors on proc 0 ...'
end if

if ((i_start_at==0).and.(i_final_at==0)) return
ndim = config_desc(iconf)%dim_desc
efilename='desc'//trim(adjustl(db_path))//config_real(iconf)%class//'_'//config_real(iconf)%klm//'_'//config_real(iconf)%cnumber//'.eml'
if (desc_forces) ffilename='desc'//trim(adjustl(db_path))//config_real(iconf)%class//'_'//config_real(iconf)%klm//'_'//config_real(iconf)%cnumber//'.fml'

!there is only one procs who writting ...
if (rangml==iproc) then

    eunit=41
    funit=42

    if (desc_forces) then
        write(6,*) 'desc_forces=.true. and para low memory not implemented. Put desc_forces to .false.'
    end if

    open(eunit,file=efilename,status='unknown')
    if (desc_forces) open(funit,file=ffilename,status='unknown')
    write(CHFMT,*)'(i6, 1x, ',int(ndim),'e20.10)'
    !write(CHFMTf,*)'(i6, 1x, e20.10,',int(ndim),'e20.10)'
    write(CHFMTf,*)'(i6, 1x,i5,1x,i6, ',int(ndim),'e20.10)'
    do ia=1, config_real(iconf)%nat
    !do ia=i_start_at, i_final_at
        write(eunit,FMT=CHFMT) ia, real(config_desc(iconf)%energy(:,ia))
        if (desc_forces) then
            write(funit,'(i6)') ia
            do ix=1,3
                write(funit,FMT=CHFMTf) ia, ix, 0,  real(config_desc(iconf)%force(:,ia,0,ix))
            end do
        end if
        if (desc_forces) then
            do ja=1,config_desc(iconf)%n_neigh(ia)
                do ix=1,3
                    !write(funit,FMT=CHFMTf) config_desc(iconf)%kind_neigh(ia,ja), config_real(iconf)%u_ij(ia,ja,ix),  config_desc(iconf)%force(:,ia,ja,ix)
                    write(funit,FMT=CHFMTf) config_desc(iconf)%kind_neigh(ia,ja), ix, ja,  config_desc(iconf)%force(:,ia,ja,ix)
                end do
            end do
        end if
    enddo

    close(eunit)
    if (desc_forces) close(funit)
end if

return
end subroutine para_write_descriptors

subroutine  init_descriptors

#ifdef PARAML
use mpi
use mod_mpi_ml
#endif
use   ml_in_ndm_module, ONLY : rangml,debug,  &
                               descriptor_type, descriptor_g2, descriptor_g3, descriptor_behler, &
                               descriptor_afs, descriptor_g2_afs, descriptor_soap, descriptor_pow_so3, descriptor_bispectrum_so3, &
                               descriptor_pow_so4,descriptor_bispectrum_so4, descriptor_g2_bispectrum_so4, &
                               descriptor_mtp, &
                               gen_param_behler, j_max, &
                               g2_dim, g3_dim, bisso4_dim, char_desc, descriptor_g2_pow_so4,  &
                               pow_so4_dim, mtp_dim, afs_dim, lbso4_diag, n_rbf, &
                               pow_so3_dim, l_max, lbso3_diag, bisso3_dim, &
                               soap_dim, weighted
use temporary_data_cov, only : dim_xdesc
use compute_pow_so3_mod
use compute_pow_so4_mod
use compute_soap_mod
use compute_mtp_mod
use compute_afs_mod
use compute_bispectrum_so3_mod
use compute_bispectrum_so4_mod
implicit none

if (debug) then
#ifdef PARAML
 call MPI_BARRIER(MPI_COMM_WORLD,codeml)
#endif
  if (rangml==0) write(6,'("ML: descriptor was initialized in init_descriptors...")')
end if


select case(descriptor_type)

  case (descriptor_g2)
     call gen_param_behler
     dim_xdesc = g2_dim
     char_desc = 'bhg2'

  case (descriptor_g3)
     call gen_param_behler
     dim_xdesc = g3_dim
     char_desc = 'bhg3'

  case (descriptor_behler)
     call gen_param_behler
     dim_xdesc = g2_dim+g3_dim
     char_desc = 'bhlr'

  case (descriptor_soap)
     call compute_cg_vector(dble(l_max),1)
     call gen_dimension_for_soap()
     dim_xdesc=soap_dim
     char_desc='soap'

  case (descriptor_pow_so3)
     pow_so3_dim=int((1 + l_max))*n_rbf
     call init_pow_so3_rbf()
     call compute_cg_vector(dble(l_max),1)
     dim_xdesc=int((1 + l_max))*n_rbf
     char_desc='pso3'

  case (descriptor_bispectrum_so3)
     call init_pow_so3_rbf()
     call compute_cg_vector(dble(l_max),1)
     if (lbso3_diag) then
       call gen_dimension_for_bispectrum_so3_diagonal()
     else
       call gen_dimension_for_bispectrum_so3_all()
     end if
     dim_xdesc=bisso3_dim
     char_desc = 'bso3'

  case (descriptor_pow_so4)
     call compute_cg_vector(dble(j_max),2)
     dim_xdesc = int(2*j_max)+1
     char_desc = 'pso4'

  case (descriptor_afs)
     call init_afs_rbf()
     if (weighted) then
       dim_xdesc = 2*afs_dim
     else
       dim_xdesc = afs_dim
     end if
     char_desc = 'afsr'

  case (descriptor_g2_afs)
     call gen_param_behler
     call init_afs_rbf()
     dim_xdesc = g2_dim+afs_dim
     char_desc = 'g2af'

  case (descriptor_g2_pow_so4)
     call gen_param_behler
     call compute_cg_vector(dble(j_max),2)
     dim_xdesc = g2_dim+int(2*j_max)+1
     pow_so4_dim = int(2*j_max)+1
     char_desc = 'g2p4'

  case (descriptor_bispectrum_so4)
     call compute_cg_vector(dble(j_max),2)
     !snap call gen_dimension_for_bispectrum_so4_all()
     ! Gabor version for which are taken only (J J_1 J_1) componenets
     ! Thompsson all componenets ...
     if (lbso4_diag) then
       call gen_dimension_for_bispectrum_so4_diagonal()
     else
       call gen_dimension_for_bispectrum_so4_all()
     end if
     if (weighted) then
       dim_xdesc=2*bisso4_dim
     else
       dim_xdesc=bisso4_dim
     end if
     char_desc = 'bso4'

  case (descriptor_mtp)
     call gen_dimension_for_mtp()
     if (weighted) then
       dim_xdesc=mtp_dim
     else
       dim_xdesc=mtp_dim
     end if
     dim_xdesc=mtp_dim
     char_desc = 'mtp3'

  case (descriptor_g2_bispectrum_so4)
     call gen_param_behler
     call compute_cg_vector(dble(j_max),2)
     !snap call gen_dimension_for_bispectrum_so4_all()
     ! Gabor version for which are taken only (J J_1 J_1) componenets
     if (lbso4_diag) then
       call gen_dimension_for_bispectrum_so4_diagonal()
     else
       call gen_dimension_for_bispectrum_so4_all()
     end if
     dim_xdesc=bisso4_dim+g2_dim
     if (rangml==0)  write(6,'("ML: hybrid descritor in init desc  bso4 + g2:  ", 2i5)') bisso4_dim, g2_dim
     char_desc = 'g2b4'
  case default
       if (rangml==0) write(6,*) 'No implementation for descriptor_type...', descriptor_type
       stop "fatal in File: descriptor.f90, Subroutine: init_descriptor"
end select

  if (rangml==0) write(6,'("ML: descriptor ",a," has the dimension ",i6)') char_desc, dim_xdesc

if (debug) then
#ifdef PARAML
 call MPI_BARRIER(MPI_COMM_WORLD,codeml)
#endif
  if (rangml==0) write(6,'("ML: descriptor was initialized in init_descriptors...")')
end if
#ifdef PARAML
 call MPI_BARRIER(MPI_COMM_WORLD,codeml)
#endif


return
end subroutine  init_descriptors


subroutine build_database_with_function(nd_local_data,dim_local_xdesc,yfunc, xdesc)

use ml_in_ndm_module, ONLY : toy_model,seed
use toy_models
implicit none
integer, intent(in) :: nd_local_data, dim_local_xdesc
!real(kind=kind(1.d0)), dimension(:,:), intent(out) ::  xdesc(dim_local_xdesc,nd_local_data)
real(kind=kind(1.d0)), dimension(:,:), intent(out) ::  xdesc
real(kind=kind(1.d0)), intent(out) ::  yfunc(nd_local_data)
!local
real(kind=kind(1.d0)) ::  rvalue, leng,y,  x(dim_local_xdesc),y_err(nd_local_data)
integer :: i, j


! internal length for data ... to see the units effetct.
leng=3.0

if (toy_model) then

! generate nd_local_data random numbers between (0,L=leng)

call random_seed(seed)

   do i=1,nd_local_data
     do j =1,dim_local_xdesc
      call random_number(rvalue)
      xdesc(j,i)=rvalue*leng
     end do
   end do


 do i=1,nd_local_data
  x(:)=real(xdesc(:,i))
  call toy_nD(x,y,leng,dim_local_xdesc)
  yfunc(i) = y
 end do

 call generate_random_gaussian(y_err,nd_local_data)


 yfunc(:) = yfunc(:) + 0.02d0*y_err(:)

else

  write(*,*) 'ML error: <build_database_with_function> the toy_model=.F. in this subroutine is not yet implemented'
stop

end if


return
end subroutine build_database_with_function


!$-------------------------------------------------------------
subroutine train_deallocate_desc(iconf)
!$-------------------------------------------------------------
use ml_in_ndm_module, only: ml_type_descriptors, ml_type
use derived_types, only:config_desc
implicit none
integer, intent(in):: iconf

if (allocated(config_desc(iconf)%force))       deallocate(config_desc(iconf)%force)
if (allocated(config_desc(iconf)%energy))      deallocate(config_desc(iconf)%energy)
if (allocated(config_desc(iconf)%n_neigh))     deallocate(config_desc(iconf)%n_neigh)
if (allocated(config_desc(iconf)%kind_neigh))  deallocate(config_desc(iconf)%kind_neigh)


if (ml_type==ml_type_descriptors) then
    if (allocated(config_desc(iconf)%type_neigh))    deallocate(config_desc(iconf)%type_neigh)
    if (allocated(config_desc(iconf)%n_neigh_ghost))     deallocate(config_desc(iconf)%n_neigh_ghost)
    if (allocated(config_desc(iconf)%kind_neigh_ghost))  deallocate(config_desc(iconf)%kind_neigh_ghost)
    if (allocated(config_desc(iconf)%kind_neigh_proc))  deallocate(config_desc(iconf)%kind_neigh_proc)
    if (allocated(config_desc(iconf)%pack_force))       deallocate(config_desc(iconf)%pack_force)
    if (allocated(config_desc(iconf)%pack_energy))       deallocate(config_desc(iconf)%pack_energy)
    if (allocated(config_desc(iconf)%pack_stress))       deallocate(config_desc(iconf)%pack_stress)
end if

return
end subroutine train_deallocate_desc
!<-------------------------------------------------------------
end module 
