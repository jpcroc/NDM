module build_kernel_matrix_mod
        implicit none 
        contains
subroutine save_matrix_para(rangml)
 use temporary_data_cov, ONLY:  i_local_cov,&
                                d_mlocal,  mlocal, u_local, v_local
 implicit none
 integer , intent(in) :: rangml

!local
 integer :: inamefile,i
 character*70 :: namefile_m,namefile_u,namefile_v,namefile_d


 inamefile=100000+rangml
 write(namefile_m,'(i6,a4)') inamefile,".mml"
 write(namefile_u,'(i6,a4)') inamefile,".uml"
 write(namefile_v,'(i6,a4)') inamefile,".vml"
 write(namefile_d,'(i6,a4)') inamefile,".dml"
 namefile_m=TRIM(namefile_m)
 namefile_u=TRIM(namefile_u)
 namefile_v=TRIM(namefile_v)
 namefile_d=TRIM(namefile_d)
 open(71,file=namefile_m,status='unknown', access='sequential',form='unformatted')
 open(72,file=namefile_u,status='unknown', access='sequential',form='unformatted')
 open(73,file=namefile_v,status='unknown', access='sequential',form='unformatted')
 open(74,file=namefile_d,status='unknown', access='sequential',form='unformatted')
! write(namefile,'(i6)') inamefile
! namefile=TRIM(namefile)
! open(71,file=namefile//".mml",status='unknown', access='sequential',form='unformatted')
! open(72,file=namefile//".uml",status='unknown', access='sequential',form='unformatted')
! open(73,file=namefile//".vml",status='unknown', access='sequential',form='unformatted')
! open(74,file=namefile//".dml",status='unknown', access='sequential',form='unformatted')
 write(71) i_local_cov
 write(74) i_local_cov
 write(*,*) 'writing on proc  and imax ',rangml, i_local_cov
 !/
 do i=1,i_local_cov
  write(71) mlocal(i)
  write(72) u_local(i)
  write(73) v_local(i)
  write(74) d_mlocal(i)
 end do
 !
 close(71)
 close(72)
 close(73)
 close(74)
return
end subroutine save_matrix_para


subroutine read_matrix_serial()
 use temporary_data_cov, ONLY:  i_local_cov,&
                                d_mlocal,  mlocal, u_local, v_local
 implicit none

!local
 integer :: inamefile,i
 character*70 :: namefile_m,namefile_u,namefile_v,namefile_d


    inamefile=100000
    write(namefile_m,'(i6,a4)') inamefile,".mml"
    write(namefile_u,'(i6,a4)') inamefile,".uml"
    write(namefile_v,'(i6,a4)') inamefile,".vml"
    write(namefile_d,'(i6,a4)') inamefile,".dml"
    namefile_m=TRIM(namefile_m)
    namefile_u=TRIM(namefile_u)
    namefile_v=TRIM(namefile_v)
    namefile_d=TRIM(namefile_d)
    open(71,file=namefile_m,status='unknown', access='sequential',form='unformatted')
    open(72,file=namefile_u,status='unknown', access='sequential',form='unformatted')
    open(73,file=namefile_v,status='unknown', access='sequential',form='unformatted')
    open(74,file=namefile_d,status='unknown', access='sequential',form='unformatted')
!    write(namefile,'(i6)') inamefile
!    namefile=TRIM(namefile)
!    open(71,file=namefile//".mml",status='unknown', access='sequential',form='unformatted')
!    open(72,file=namefile//".uml",status='unknown', access='sequential',form='unformatted')
!    open(73,file=namefile//".vml",status='unknown', access='sequential',form='unformatted')
!    open(74,file=namefile//".dml",status='unknown', access='sequential',form='unformatted')
    read(71) i_local_cov
    read(74) i_local_cov
    write(*,*) 'reading on proc  and imax ', i_local_cov
 !
    if (allocated(mlocal)) deallocate(mlocal)
    if (allocated(d_mlocal)) deallocate(d_mlocal)
    if (allocated(u_local)) deallocate(u_local)
    if (allocated(v_local)) deallocate(v_local)

    allocate(mlocal(i_local_cov),d_mlocal(i_local_cov),u_local(i_local_cov),v_local(i_local_cov))
    do i=1,i_local_cov
     read (71) mlocal(i)
     read (72) u_local(i)
     read (73) v_local(i)
     read (74) d_mlocal(i)
    end do
   close(71)
   close(72)
   close(73)
   close(74)

return
end subroutine read_matrix_serial


subroutine save_matrix_serial()
  use temporary_data_cov, ONLY:  i_local_cov,&
                                 d_mlocal, mlocal, u_local, v_local
 implicit none

!local
 integer :: inamefile,i
 character*70 :: namefile_m,namefile_u,namefile_v,namefile_d


    inamefile=100000
    write(namefile_m,'(i6,a4)') inamefile,".mml"
    write(namefile_u,'(i6,a4)') inamefile,".uml"
    write(namefile_v,'(i6,a4)') inamefile,".vml"
    write(namefile_d,'(i6,a4)') inamefile,".dml"
    namefile_m=TRIM(namefile_m)
    namefile_u=TRIM(namefile_u)
    namefile_v=TRIM(namefile_v)
    namefile_d=TRIM(namefile_d)
    open(71,file=namefile_m,status='unknown', access='sequential',form='unformatted')
    open(72,file=namefile_u,status='unknown', access='sequential',form='unformatted')
    open(73,file=namefile_v,status='unknown', access='sequential',form='unformatted')
    open(74,file=namefile_d,status='unknown', access='sequential',form='unformatted')
!    write(namefile,'(i6)') inamefile
!    namefile=TRIM(namefile)
!    open(71,file=namefile".mml",status='unknown', access='sequential',form='unformatted')
!    open(72,file=namefile".uml",status='unknown', access='sequential',form='unformatted')
!    open(73,file=namefile".vml",status='unknown', access='sequential',form='unformatted')
!    open(74,file=namefile".dml",status='unknown', access='sequential',form='unformatted')
    write(71) i_local_cov
    write(74) i_local_cov
    write(*,*) 'writing on proc  and imax ', i_local_cov
 !
    do i=1,i_local_cov
     write(71) mlocal(i)
     write(72) u_local(i)
     write(73) v_local(i)
     write(74) d_mlocal(i)
    end do
 !
   close(71)
   close(72)
   close(73)
   close(74)

return
end subroutine save_matrix_serial


subroutine set_neighbors_size_cov_ml (i_start_cov,i_final_cov,i_local_cov, dim_train)
 USE T_kind_param_m, ONLY:  double
 use gen_com_m, ONLY: rangml
 use temporary_data_cov, ONLY:  neigh_cov, list_cov
 use ml_in_ndm_module, ONLY: debug

 implicit none
 integer, intent (in) :: i_start_cov,i_final_cov,dim_train
 integer, intent(out) :: i_local_cov

 integer :: i

 if (debug) then
 if (rangml==0) write (6, *) 'ML: building the index of the kernel matrix on each proc .... '
 end if

 ! This situation correspond to no rcut in the finger-print space.
 ! Hence in the following lines each point in data space has dim_data neighbors

 if (allocated(neigh_cov)) deallocate(neigh_cov)
 allocate (neigh_cov(dim_train))

 if (allocated(list_cov)) deallocate(list_cov)
 allocate(list_cov(dim_train))

 do i=i_start_cov,i_final_cov
  neigh_cov(i) = dim_train
 end do

 do i=1,dim_train
  list_cov(i)=i
 end do

! if you want only half of the matrix uncommebnt those lines
!
! do i=i_start_cov,i_final_cov
!  allocate (list_cov(dim_data,i-1))
! end do
!
! do i=i_start_cov,i_final_cov
!  neigh_cov(i) = i
!   do j=1,i-1
!     list_cov(i,j)=j
!  end do
! end do


! size of the matrix which will be stored in the local procs ...
 i_local_cov=(i_final_cov+1-i_start_cov)*dim_train


!this subroutine is usefull if we intyroduce devide and conquer (or rcut in the fingerprint space).
! you should uncomment and addapt the following loop.
!  do j=i_start,i_final
!          do iw=1,dim_data
!               if distance(data_i - data_j) < "distance critical" ilocal=ilocal+1
!          end do
!  end do

 !

end subroutine set_neighbors_size_cov_ml



subroutine partial_covariance_matrix(local_xdesc, dim_tr, dim_xd, i_start_cov,i_final_cov)
USE T_kind_param_m, ONLY:  double
use temporary_data_cov, ONLY: dim_xdesc, dim_train, neigh_cov, list_cov,&
                              mlocal,u_local,v_local
use def_kernels
use ml_in_ndm_module, ONLY: rangml, kernel_type, marginal_likelihood

implicit none
integer, intent(in)      :: i_start_cov,i_final_cov,dim_tr,dim_xd
real(double), intent(in) :: local_xdesc(dim_xdesc, dim_train)
real(double) :: k_ij,d_length_k_ij
integer :: i,j, jj, ilocal,n_j

if (dim_tr/=dim_train) then
  if (rangml==0) write(6,*) 'dim_tr and dim_train have illegal value. ', dim_tr, dim_train
  if (rangml==0) write(6,*) 'both variables should be equals in order to correctly build the covariance matrix'
  stop
end if


if (dim_xd/=dim_xdesc) then
  if (rangml==0) write(6,*) 'dim_xd and dim_xdesc have illegal value. ', dim_tr, dim_train
  if (rangml==0) write(6,*) 'both variables should be equals in order to correctly build the covariance matrix'
  stop
end if



ilocal=0



 do i = i_start_cov,i_final_cov
    n_j=neigh_cov(i)
     do jj=1,n_j
         j=list_cov(jj)
         ilocal=ilocal+1
         select case(kernel_type)
          case (kernel_se)
             k_ij=kernel_square_exp(local_xdesc(:,i),local_xdesc(:,j),dim_xdesc)
             if (marginal_likelihood) d_length_k_ij=d_length_kernel_square_exp(local_xdesc(:,i),local_xdesc(:,j),dim_xdesc)
          case (kernel_ou)
             k_ij=kernel_ornstein_uhlenbeck (local_xdesc(:,i),local_xdesc(:,j),dim_xdesc)
          case (kernel_mc)
             k_ij=kernel_matern_class(local_xdesc(:,i),local_xdesc(:,j),dim_xdesc)
          case (kernel_so)
             k_ij=kernel_soap(local_xdesc(:,i),local_xdesc(:,j),dim_xdesc)
        end select
        mlocal      (ilocal) = k_ij
!       d_mlocal     (ilocal) = d_length_k_ij
        u_local     (ilocal) = i
        v_local     (ilocal) = j

    end do

 end do


return
end subroutine  partial_covariance_matrix




subroutine build_kernel_matrix(local_xdesc,dim_xdesc, dim_train)
! here wdesc is not the general xdesc ... probably i should change that in the future ...
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
!  use tab_imm_m, ONLY : im
  use ml_in_ndm_module
  use temporary_data_cov, ONLY:  i_start_cov, i_final_cov,i_local_cov,&
                                 d_mlocal, mlocal, u_local, v_local
#ifdef PARAML
  use mpi
  use mod_mpi_ml
#endif

  implicit none
  integer, intent(in) :: dim_xdesc, dim_train
  real(double), intent(in) :: local_xdesc(dim_xdesc,dim_train)
  !-----------------------------------------------
  !   L o c a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: count1,count_rate,count_max    ! time counting.


! ML ERROR SHOULB BE CHANGED ....
! ERROR -----
!  integer :: inamefile
!  character*70 :: namefile


!-------END ERROR---------------/

  call system_clock (count1,count_rate,count_max)

  call set_neighbors_size_cov_ml(i_start_cov,i_final_cov,i_local_cov,dim_train)

 if (debug) then
#ifdef PARAML
  write(*,'(" proc i_start_cov i_final_cov natoms",i5,3i7,i10)') rangml,i_start_cov,i_final_cov,i_final_cov-i_start_cov,i_local_cov
#else
  write(*,*) 'i_start_cov i_final_cov i_local_cov ', i_start_cov, i_final_cov, i_local_cov
#endif
  end if

  if (i_final_cov<=i_start_cov) then
   if (rangml==0) write(*,*) 'i_final_cov should be greater than i_start_cov',i_final_cov, ' > ', i_start_cov
   if (rangml==0) write(*,*) 'stop in force_constat.F90. Change no of procs.'
   stop
  end if

  if (allocated(mlocal))   deallocate(mlocal)  ; allocate(mlocal(i_local_cov))
  if (allocated(d_mlocal)) deallocate(d_mlocal); allocate(d_mlocal(i_local_cov))

  if (allocated(u_local)) deallocate(u_local); allocate(u_local(i_local_cov))
  if (allocated(v_local)) deallocate(v_local); allocate(v_local(i_local_cov))

  if ((iread_ml==0).or.(isave_ml==2)) then
     call partial_covariance_matrix(local_xdesc, dim_train, dim_xdesc,  i_start_cov,i_final_cov)
  end if

  if (iread_ml==1) then
    if (rangml==0) then
      write(*,*) 'The KERNEL MATRIX are not computed. Are readed from an previuos run'
      if (isave_ml/=0) then
       write(*,*) 'WARINIIIIIIIIIIIIG isave_ml is set to', isave_ml
       write(*,*) 'We will switch isave_ml to 0'
      end if
     end if
   isave_ml=0
  end if

!This is MPI version ........
! #if ( PARAML && ML ) || ( PARAML && ML && PHONDY && PARAPH )
#if defined PARAML && defined ML 
if (isave_ml==1) then   !Writing on disk in order to be read by ScaLpack
                        !diagonalization program ...
 call save_matrix_para(rangml)
end if


#endif
!This is serial version ...
#ifndef PARAML


  if (iread_ml==1) then     ! serial version reading the matrix
   call read_matrix_serial()
  end if   !iread_ml==1 for serial version

  if (isave_ml==1) then    !serial only store the matrix
   call save_matrix_serial()
  end if !isave_ml==1


  !call system_clock (count2,count_rate,count_max)
  !time=real((count2-count1))/real(count_rate)
  !if (iread_ml /=1) then
  ! if (debug) then
  ! if (rangml==0) write(*,"(' ML: local matrix  were filled in.......:  ',f16.8,' s')") time
  ! end if
  !end if
#endif


if (debug) then
  if (rangml==0) write(*,*) 'ML: ...the kernel matrix was computed'
end if
end subroutine build_kernel_matrix
end module
