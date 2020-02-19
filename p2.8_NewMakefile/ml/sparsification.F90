module sparsification_mod
        use svd_mod
        use read_poscar_mod
        use alloc_typ_ml_mod
        implicit none
        contains
subroutine try_sparsification(n_count)

      USE T_kind_param_m, ONLY:  double
      use gen_com_m
      use var_pot
      use tab_imm_m
      use ml_in_ndm_module
      use temporary_data_cov, ONLY:  yfunc, yfunc_nd, xdesc, dim_data, dim_valid, dim_train, dim_xdesc, dim_data_train, dim_data_test, &
                                     i_final_cov, i_start_cov, i_local_cov, &
                                     error_valid, error_test, &
                                     yfunc_valid, yfunc_train, y_test, &
                                     xdesc_train, xdesc_valid,xdesc_test,y_extra,eigen_values, &
                                     dim_yfunc, yfunc_average
     use extrapolation
     use k_cross_validation
     use set_limits
     use def_kernels, ONLY : length_kse
     use opt_marginal_likelihood

#ifdef PARAML
!for curie      use mkl_service
      use mpi
      use mod_mpi_ml
#endif
implicit none
integer :: j, n_count


if ((sparsification_by_entropy) .and.(sparsification_by_acd)) then
  if (rangml==0) then
     write(6,*) 'ML error: <sparsification> The sparsification cannot be simultanously by entropy or acd'
     stop
  end if
end if

if (rangml==0) write(6,*) 'ML: sparsification ...'

if (sparsification_by_entropy) call sub_sparsification_by_entropy()
if (sparsification_by_acd) call sub_sparsification_by_acd()


#ifdef PARAML
       call MPI_BCAST(reject,ns_data,MPI_LOGICAL,0,MPI_COMM_WORLD,codeml)
#endif
       n_count=0
       do j=1,ns_data
          if (reject(j).eqv..false.) n_count=n_count+1
       enddo
       if (allocated(data_im)) deallocate(data_im); allocate(data_im(0:n_count))
       data_im(0)=0

return

end subroutine try_sparsification


subroutine sub_sparsification_by_entropy()


      USE T_kind_param_m, ONLY:  double
      use gen_com_m
      use var_pot
      use tab_imm_m
      use ml_in_ndm_module
      use temporary_data_cov, ONLY:  yfunc, xdesc,    dim_xdesc, &
                                     error_valid, &
                                     dim_yfunc, yfunc_average
     use extrapolation
     use k_cross_validation
     use set_limits
     use opt_marginal_likelihood
     use dynalloccell
     use divid_mod
     use neigcel_mod
     use caltabi_mod
     use caltabt_mod
     use compute_descriptors_mod
     use recips_mod

#ifdef PARAML
!for curie      use mkl_service
      use mpi
      use mod_mpi_ml
#endif
implicit none

!interface compute_descriptors
!   subroutine compute_descriptors(xdesc_out,icount)
!    implicit none
!    double precision,dimension(:,:),allocatable :: xdesc_out
!    integer, optional ::  icount
!   end subroutine compute_descriptors
!end interface compute_descriptors

double precision,dimension(:),allocatable :: diff_entropy
integer,dimension(:),allocatable :: idx
double precision :: y_target,entropy
double precision,dimension(:,:),allocatable :: xdesc_i
integer :: n_temp, nb, j,i, i_count, i2_count


          if (rangml==0) write(6,*) 'ML: sparsification using entropy ...'

          if (allocated(diff_entropy)) deallocate(diff_entropy); allocate(diff_entropy(0:ns_data))
          if (allocated(idx)) deallocate(idx); allocate(idx(ns_data))

          do j=0,ns_data
              reject(:)=.false.
              if (j==0) then
                 n_temp=ns_data
              else
                 n_temp=ns_data-1
                 reject(j)=.true.
              endif

            size_indi2=0
            max_ntyp=0
            i_count=0
            i2_count=0
            select case (target_type)
            case (target_energy)
               if (allocated(xdesc)) deallocate(xdesc)
               do i=i_begin+1,i_begin+ns_data
                  i_count=i_count+1
                  if (reject(i_count)) cycle
                  i2_count=i2_count+1
                  call read_poscar(pref,i,y_target) ! Mise à jour de im,at,bg,xp,ntyp,ityp
                  call convert_A2cm(1,.true.,.true.)
                  im_glob=im                        ! Mise à jour de im_glob pour divid
                  call alloc_typ_ml()
                  call Deallocatecel()
                  natperc=-1                        ! Force le calcul de natperc dans divid
                  call divid(1)
                  call DynamicalAllocationCell()    ! Réallocation du pointeur last(natperc,:noxyz) pour caltabt
                  call neigcel()
                  call caltabt()
                  call caltabi()
                  if (max_ntyp<size(natm)) then
                      max_ntyp=size(natm)
                  endif
                  call allocate_ml()
                  call compute_descriptors(xdesc_i)
                  if (.not.allocated(xdesc)) allocate(xdesc(dim_xdesc,n_temp))
                  xdesc(1:dim_xdesc,i2_count)=xdesc_i(1:dim_xdesc,1)
                  call deallocate_ml()
               enddo
            case (target_force)
               if (allocated(xdesc)) deallocate(xdesc)
               do i=i_begin+1,i_begin+ns_data
                  i_count=i_count+1
                  if (reject(i_count)) cycle
                  call read_poscar(pref,i,y_target) ! Mise à jour de im,at,bg,xp,ntyp,ityp
                  call convert_A2cm(1,.true.,.true.)
                  call recips(at(:,1), at(:,2), at(:,3), bg(:,1), bg(:,2), bg(:,3))
                  im_glob=im                        ! Mise à jour de im_glob pour divid
                  call alloc_typ_ml()
                  call Deallocatecel()
                  natperc=-1                        ! Force le calcul de natperc dans divid
                  call divid(1)
                  call DynamicalAllocationCell()    ! Réallocation du pointeur last(natperc,:noxyz) pour caltabt
                  call neigcel()
                  call caltabt()
                  call caltabi()
                  if (size_indi2<size(indi2)) then
                     size_indi2=size(indi2)
                  endif
                  if (max_ntyp<size(natm)) then
                      max_ntyp=size(natm)
                  endif
                  call allocate_ml()
                  call compute_descriptors(xdesc_i)
                  if (.not.allocated(xdesc)) allocate(xdesc(dim_xdesc,n_temp))
                  xdesc(1:dim_xdesc,nb+1:nb+im)=xdesc_i(1:dim_xdesc,1:im)
                  nb=nb+im
                  call deallocate_ml()
               enddo
            end select

            call svd(entropy,n_temp)
            if (j==0) then
               diff_entropy(j)=entropy
            else
               diff_entropy(j)=diff_entropy(0)-entropy
            endif

          enddo !j

          reject(:)=.true.
          call indexx(ns_data,diff_entropy(1:ns_data),idx)
          do j=0,max_data-1
             reject(idx(ns_data-j))=.false.
          enddo !j

          if (debug.and.(rangml==0)) write(*,*)"diff_entropy",diff_entropy(1:)

          deallocate(diff_entropy,idx)

          if (debug.and.(rangml==0)) write(*,*)"reject",reject(:)

          if (rangml==0) write(6,*) 'ML: end sparsification using entropy ...'

return
end subroutine sub_sparsification_by_entropy


subroutine sub_sparsification_by_acd
USE T_kind_param_m, ONLY:  double
use gen_com_m
use var_pot
use tab_imm_m
use ml_in_ndm_module
use temporary_data_cov, ONLY:  yfunc_nd,     &
                               i_start_cov, &
                               xdesc_train, xdesc_test,&
                               dim_yfunc
use extrapolation
use k_cross_validation
use set_limits
use opt_marginal_likelihood
use dynalloccell
use divid_mod
use neigcel_mod
use caltabi_mod
use caltabt_mod
use compute_acd_mod
use compute_acd_local_mod


#ifdef PARAML
!for curie      use mkl_service
      use mpi
      use mod_mpi_ml
#endif
implicit none
integer ::  i, i_count

double precision :: y_target
          if (rangml==0) write(6,*) 'ML: sparsification using ACD ...'
          reject(:)=.true.

          size_indi2=0
          max_ntyp=0
          if (allocated(data_im)) deallocate(data_im)
          allocate(data_im(0:ns_data))
          data_im(0)=0
          i_count=0
          select case (target_type)
          case (target_energy)
             do i=i_begin+1,i_begin+ns_data
                i_count=i_count+1
                call read_poscar(pref,i,y_target) ! Mise à jour de im,at,bg,xp,ntyp,ityp
                data_im(i_count)=im
                if (max_ntyp<size(natm)) then
                    max_ntyp=size(natm)
                endif
             enddo
          case (target_force)
             do i=i_begin+1,i_begin+ns_data
                i_count=i_count+1
                call read_poscar(pref,i,y_target) ! Mise à jour de im,at,bg,xp,ntyp,ityp
                data_im(i_count)=im
                if (size_indi2<size(indi2)) then
                   size_indi2=size(indi2)
                endif
                if (max_ntyp<size(natm)) then
                    max_ntyp=size(natm)
                endif
             enddo
          end select

          if (allocated(r_acd)) deallocate(r_acd)
          allocate(r_acd(3,maxval(data_im),ns_data))
          if (allocated(data_natm)) deallocate(data_natm)
          allocate(data_natm(max_ntyp,ns_data))
          i_count=0
          select case (target_type)
          case (target_energy)
             do i=i_begin+1,i_begin+ns_data
                i_count=i_count+1
                call read_poscar(pref,i,y_target) ! Mise à jour de im,at,bg,xp,ntyp,ityp
                data_natm(:,i_count)=natm(:)
                r_acd(:,:,i_count)=xp
             enddo
          case (target_force)
             if (allocated(iwmax2_acd)) deallocate(iwmax2_acd)
             if (allocated(indi2_acd)) deallocate(indi2_acd)
             if (allocated(at_acd)) deallocate(at_acd)
             allocate(iwmax2_acd(maxval(data_im),ns_data))
             allocate(indi2_acd(size_indi2,ns_data))
             allocate(at_acd(3,3,ns_data))
             do i=i_begin+1,i_begin+ns_data
                i_count=i_count+1
                call read_poscar(pref,i,y_target) ! Mise à jour de im,at,bg,xp,ntyp,ityp
                data_natm(:,i_count)=natm(:)
                call convert_A2cm(1,.true.,.true.)
                im_glob=im                        ! Mise à jour de im_glob pour divid
                call alloc_typ_ml()
                call Deallocatecel()
                natperc=-1                        ! Force le calcul de natperc dans divid
                call divid(1)
                call DynamicalAllocationCell()    ! Réallocation du pointeur last(natperc,:noxyz) pour caltabt
                call neigcel()
                call caltabt()
                call caltabi()
                call coord_acd_local(i_count)
             enddo
          end select

          if (allocated(k_acd)) deallocate(k_acd)
          allocate(k_acd(ns_data,ns_data))
          if (allocated(distance_acd)) deallocate(distance_acd)
          allocate(distance_acd(ns_data,ns_data))
          select case (target_type)
          case (target_energy)
             call compute_kernel_acd(ns_data,k_acd,distance_acd)
          case (target_force)
             call compute_kernel_acd_local(ns_data,k_acd,distance_acd)
          end select

          deallocate(k_acd,distance_acd,r_acd,data_natm)

          if (debug.and.(rangml==0)) write(*,*)"reject ",reject(:)

          if (rangml==0) write(6,*) 'ML: sparsification using ACD ...'

return
end subroutine sub_sparsification_by_acd
end module
