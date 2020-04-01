
subroutine get_neighbours_mpi
use mpi
use mod_mpi_phondy
use derived_types_ph, only: config_real
use phondy_in_ndm_module, only: iconf_ini,NN_MAX
implicit none
integer:: dim_reduce

    call MPI_BARRIER(MPI_COMM_WORLD, codeph)
    dim_reduce=config_real(iconf_ini)%nat
    call MPI_ALLREDUCE(MPI_IN_PLACE, config_real(iconf_ini)%n_neigh,dim_reduce,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD,codeph)

    dim_reduce=config_real(iconf_ini)%nat*NN_MAX
    call MPI_ALLREDUCE(MPI_IN_PLACE, config_real(iconf_ini)%kind_neigh ,dim_reduce,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD,codeph)


    dim_reduce=config_real(iconf_ini)%nat*NN_MAX
    call MPI_ALLREDUCE(MPI_IN_PLACE, config_real(iconf_ini)%kind_neigh_big ,dim_reduce,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD,codeph)

    dim_reduce=config_real(iconf_ini)%nat*NN_MAX
    call MPI_ALLREDUCE(MPI_IN_PLACE, config_real(iconf_ini)%type_neigh ,dim_reduce,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD,codeph)

    dim_reduce=config_real(iconf_ini)%nat*NN_MAX
    call MPI_ALLREDUCE(MPI_IN_PLACE,config_real(iconf_ini)%r_ij,dim_reduce,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,codeph)

    dim_reduce=3*config_real(iconf_ini)%nat*NN_MAX
    call MPI_ALLREDUCE(MPI_IN_PLACE,config_real(iconf_ini)%u_ij,dim_reduce,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,codeph)


    dim_reduce=3*config_real(iconf_ini)%nat*NN_MAX
    call MPI_ALLREDUCE(MPI_IN_PLACE,config_real(iconf_ini)%uperiod_ij,dim_reduce,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,codeph)
return
end subroutine get_neighbours_mpi


! ************************************************
!           Sous-programme force_constant
! ************************************************
subroutine force_constant
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use var_pot
!  use tab_imm_m, ONLY : im
  use phondy_in_ndm_module
  use derived_types_ph, only : config_real
  use mpi
  use mod_mpi_phondy

  implicit none
  !-----------------------------------------------
  !   L o c a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: i, j
  real(double)::time, tempi
  integer :: count1,count2,count_rate,count_max
  integer :: nlocal,u,v

  call system_clock (count1,count_rate,count_max)


  call partial_take_size(i_start,i_final,imax)

  if (lq_points) call get_neighbours_mpi


!  if (rangph==0) then
!    write(6,'("positions of atoms 1000 in  ", 3f20.8)') config_real(iconf_ini)%pos_cart(:,config_real(iconf_big)%ia_ini(1000))
!    write(6,'("Rperiod for 1000 in  ", 3f20.8)')        config_real(iconf_big)%Rperiodic(:,1000)
!  end if


  write(*,'("PHONDY:  proc i_start i_final natoms",i5,3i7,i10)') rangph,i_start,i_final,i_final-i_start+1,imax

  nlocal=(i_final-i_start)*3
  if (i_final<i_start) then
   if (rangph==0) write(*,*) 'PHONDY: i_final should be greater or equal to than i_start',i_final, ' >= ', i_start
   if (rangph==0) write(*,*) 'PHONDY: stop in force_constant.F90. Change no of procs (decrease the number of procs).'
   stop
  end if


  if (allocated(mlocal))  deallocate(mlocal) ; allocate (mlocal(imax))
  if (allocated(u_local)) deallocate(u_local); allocate (u_local(imax))
  if (allocated(v_local)) deallocate(v_local); allocate (v_local(imax))
  if (allocated(w))       deallocate(w)      ; allocate (w(nmat))

  if ((iread==0).or.(isave==2)) then
    call partial_hessian (i_start,i_final,imax,mlocal,u_local,v_local,HessianOrder)
  end if

  !ph write(*,*) 'get out of partial index', rangph
  call MPI_BARRIER(MPI_COMM_WORLD,codeph)
  !ph write(*,*) 'after get out of partial index', rangph

  if (iread==1) then
    if (rangph==0) then
      write(*,*) 'PHONDY: The force constants are not computed. we read them from a previuos run.'
      if (isave/=0) then
       write(*,*) 'PHONDY: WARINIIIIIIIIIIIING isave is set to', isave
       write(*,*) 'PHONDY: We will switch isave to 0'
      end if
     end if
   isave=0
  end if

 !Writing on disk in order to be read by ScaLpack diagonalization program ...
  if (isave==1) then
    inamefile=100000+rangph
    write(namefile,'(i6)') inamefile
    namefile=TRIM(namefile)
    open(71,file=namefile//".m",status='unknown', access='sequential',form='unformatted')
    open(72,file=namefile//".u",status='unknown', access='sequential',form='unformatted')
    open(73,file=namefile//".v",status='unknown', access='sequential',form='unformatted')
    write(71) imax
    !ph write(*,*) 'writing on proc  and imax ',rangph, imax
    !
    do i=1,imax
      write(71) mlocal(i)
      write(72) u_local(i)
      write(73) v_local(i)
    end do
    !
    close(71)
    close(72)
    close(73)
  end if !isave==1

  if (isave==2) then
    if (rangph == 0 ) then
      if (lq_points) then
        if (allocated(matfor)) deallocate(matfor) ; allocate(matfor(nmat,3*NN_MAX))
      else
        if (allocated(matfor)) deallocate(matfor) ; allocate(matfor(nmat,nmat))
      end if
      write(6,'("PHONDY: marfor has the size: ",i9," X ", i9)') size(matfor,DIM=1), size(matfor,DIM=2)
      write(6,'("PHONDY: u, v, mlocal has the size outside, force constant : ",i12, i12, i12)') size(u_local,DIM=1), size(v_local,DIM=1),size(mlocal,DIM=1)
      matfor(:,:)=0.d0

        do i=1,imax
          u=u_local(i)
          v=v_local(i)
          matfor(u,v)=mlocal(i)
        end do
    end if

    if (rangph /= 0 ) then
       call MPI_SEND(mlocal,imax,MPI_DOUBLE_PRECISION,0,10000+rangph,MPI_COMM_WORLD,codeph)
       call MPI_SEND(u_local,imax,MPI_INTEGER,0,20000+rangph,MPI_COMM_WORLD,codeph)
       call MPI_SEND(v_local,imax,MPI_INTEGER,0,30000+rangph,MPI_COMM_WORLD,codeph)
    end if

    if (rangph == 0)  then
      iproc=0
      do iproc=1,nb_procsph-1
       call MPI_PROBE(MPI_ANY_SOURCE,MPI_ANY_TAG, MPI_COMM_WORLD,statut,codeph)
       call MPI_GET_COUNT(statut,MPI_DOUBLE_PRECISION,nb_elements,codeph)
       !debug write(*,*) 'on proc  ', iproc, 'imax is', nb_elements
       deallocate(mlocal,u_local,v_local)
       allocate(u_local(nb_elements),v_local(nb_elements), mlocal(nb_elements))
       itempproc=statut(MPI_SOURCE)
       call MPI_RECV(mlocal,nb_elements,MPI_DOUBLE_PRECISION,statut(MPI_SOURCE), &
                     10000+statut(MPI_SOURCE),MPI_COMM_WORLD,statut,codeph)
       call MPI_RECV(u_local,nb_elements,MPI_INTEGER,itempproc,20000+itempproc,MPI_COMM_WORLD,statut,codeph)
       call MPI_RECV(v_local,nb_elements,MPI_INTEGER,itempproc,30000+itempproc,MPI_COMM_WORLD,statut,codeph)
       do i=1,nb_elements
         u=u_local(i)
         v=v_local(i)
         matfor(u,v)=mlocal(i)
       end do
      end do
    end if
    !
  end if ! isave=2
  call MPI_BARRIER(MPI_COMM_WORLD,codeph)
  !call MPI_BARRIER(MPI_COMM_WORLD,codeph)
!<-------here end serial version----->

  call system_clock (count2,count_rate,count_max)
  time=real((count2-count1))/real(count_rate)
  if (rangph==0) write(*,"('PHONDY: MATFOR was filled in.......:  ',f16.8,' s')") time

if (isave==2) then
  ! Here we enter only in the case not Q-points ....
  if (rangph==0) then

    if (.not.lq_points) then
      do i=1,nmat
        do j=i,nmat
          if (i.eq.j) cycle
          matfor(i,j)=0.5d0*(matfor(i,j)+matfor(j,i))
          matfor(j,i)=matfor(i,j)
        end do
      end do


      do i=1,nmat
      !Apply the sum rule for the diagonal elements ...
         tempi=SUM(matfor(i,:))-matfor(i,i)
         !debug write(70,*) i, tempi+matfor(i,i)
         if ( (tempi+matfor(i,i)).gt.1d-10) write(6,'("PHONDY warning: The sum rule is not respected imat, deviance matfor",i9,d20.10)') i, tempi+matfor(i,i)
         matfor(i,i) =-tempi
      end do

    end if !lq_points




  end if  !rangph
end if !isave=2

  call system_clock (count1,count_rate,count_max)
  time=real((count1-count2))/real(count_rate)
  if (rangph==0) write(*,"('PHONDY: MATFOR was symmetrized in...:  ',f16.8,' s')") time

if (isave /= 2 ) then
  if (rangph==0) write(*,'("PHONDY: ...the force constants were  filled")')
    !  call MPI_BARRIER(MPI_COMM_WORLD,codeph)

   if (rangph==0) write(*,*) 'PHONDY: MATFOR', u_local(1), v_local(1), mlocal(1)
   if (rangph==0) write(*,*) 'PHONDY: MATFOR', u_local(2), v_local(2), mlocal(2)
   if (rangph==0) write(*,*) 'PHONDY: MATFOR', u_local(3), v_local(3), mlocal(3)
   if (rangph==0) write(*,*) 'PHONDY: MATFOR', u_local(4), v_local(4), mlocal(4)
end if
!debug   write(*,*) 'test rangph', rangph, imax
end subroutine force_constant
