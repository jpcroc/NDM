subroutine ph_build_box_ml
use derived_types, only: config_real, config_real_copy
use derived_types_ph, only : atom_ph, inv_atom_ph

real(kind(0.d0)) :: Rtemp(3)
real(kind(0.d0)), dimension(:,:), allocatable :: xpnp
integer :: icount, ia, i1,i2,i3, nn1, nn2, nn3

!if (debug) call print_message('Enter in subroutine...','ph_build_box_ml')
nn1=int(config_real(iconf)%nxCell/2)+1
nn2=int(config_real(iconf)%nyCell/2)+1
nn3=int(config_real(iconf)%nzCell/2)+1

imm=(2*nn1+1)*(2*nn2+1)*(2*nn3+1)

if (iconf /= 1) then
  if (rangph==0) write(6,*) 'iconf should be 1 in PH mode'
  stop 'ph_build_box_ml iconf is not 1'
end if

if (allocated(config_real_copy)) deallocate(config_real_copy) ; allocate (config_real_copy(1))

config_real_copy(iconf) = config_real(iconf)
if (allocated(config_real_copy(iconf)%pos_cart)) deallocate(config_real_copy(iconf)%pos_cart) ; allocate (config_real_copy(iconf)%pos_cart(3,imm))
if (allocated(config_real_copy(iconf)%itype)) deallocate(config_real_copy(iconf)%itype) ; allocate (config_real_copy(iconf)%itype(imm))
if (allocated(atom_ph)) deallocate(atom_ph) ; allocate(atom_ph(imm))
if (allocated(inv_atom_ph)) deallocate(inv_atom_ph) ; allocate(inv_atom_ph(config_real(iconf)%nat,-nn1:nn1,-nn2:nn2, -nn3:nn3))

config_real_copy(iconf)%volume = config_real(iconf)%volume*dble(imm)
config_real_copy(iconf)%ntypes=config_real(iconf)%ntypes
config_real(iconf)%nat=imm
im=imm


if (allocated(xpnp)) deallocate(xpnp); allocate (xpnp(3,imm))

icount=0
do  ia=1,config_real(iconf)%nat

       do i1=-nn1,nn1
       do i2=-nn2,nn2
       do i3=-nn3,nn3
         icount = icount + 1
         Rtemp(:) = dble(i1)*config_real(iconf)%cell(:,1) + &
                     dble(i2)*config_real(iconf)%cell(:,2) + &
                     dble(i3)*config_real(iconf)%cell(:,3)
         config_real_copy(iconf)%pos_cart(:,icount)=Rtemp(:) + config_real(iconf)%pos_cart(:,ia)
         config_real_copy(iconf)%itype(icount) = config_real(iconf)%itype(ia)
         atom_ph(icount)%nn1 = i1
         atom_ph(icount)%nn2 = i2
         atom_ph(icount)%nn3 = i3
         atom_ph(icount)%ia = ia
         inv_atom_ph(ia,i1,i2,i3) = icount
       end do
       end do
       end do
end do


!config_real(iconf)%pos_crst(3,1:imm) = xc(3,1:imm)/A2cm
!config_real(iconf)%force(1:3,1:imm) = fp(1:3,1:imm)*(A2cm*erg2ev)
config_real_copy(iconf)%cell(1,:) = config_real(iconf)%cell(1,:)*dble(2*nn1+1)
config_real_copy(iconf)%cell(2,:) = config_real(iconf)%cell(2,:)*dble(2*nn2+1)
config_real_copy(iconf)%cell(3,:) = config_real(iconf)%cell(3,:)*dble(2*nn3+1)
config_real_copy(iconf)%nat = icount
call recips (config_real_copy(iconf)%cell(1,1), config_real_copy(iconf)%cell(1,2), config_real_copy(iconf)%cell(1,3), &
             config_real_copy(iconf)%bg_cell(1,1), config_real_copy(iconf)%bg_cell(1,2), config_real_copy(iconf)%bg_cell(1,3))

config_real_copy(iconf)%mass_per_type(:) = config_real(iconf)%mass_per_type(:)
!future config_real(iconf)%prev_pos_cart(1:3,1:imm)=xpp(1:3,1:imm)/A2cm

!if (debug) call print_message('Exit from subroutine...','ph_build_box_ml')
return
end subroutine ph_build_box_ml


!subroutine compute_force_constants
!use derived_types_ph, only: atom_ph, inv_atom_ph
!
!!Negative displacements
!       do ic=1,icount
!          do ii=1,3
!          do jj
!
!
!       end do
!
!
!        do j=1,config_real(iconf)%nat
!          do ij=1,3
!            do lx=-lxmax,lxmax
!            do ly=-lymax,lymax
!            do lz=-lzmax,lzmax
!              do i=1,config_real(iconf)%nat
!                icount = inv_atom_ph(,)
!                do ii=1,3
!                  xpref (:,:) = config_real_copy(iconf)%pos_cart(:,:)
!                  xp(:,:) = xpref(:,:)
!                  xp(:,icount) = xpref(:,icount) + deltax
!                  call calfo_phondy(j,j)
!                  pp(:,i,lx,ly,lz) = fp(:,i)
!
!                  xpref (:,:) = config_real_copy(iconf)%pos_cart(:,:)
!                  xp(:,:) = xpref(:,:)
!                  xp(:,icount) = xpref(:,icount) - deltax
!                  call calfo_phondy(j,j)
!                  pm(:,i,lx,ly,lz) = fp(:,i)
!
!                enddo
!              end do
!            enddo
!            enddo
!            enddo
!  c Positive displacements
!            do lx=-lxmax,lxmax
!            do ly=-lymax,lymax
!            do lz=-lzmax,lzmax
!              do i=1,natoms
!                read(iunit2,*) (pn(ii,i,lx,ly,lz),ii=1,3)
!              enddo
!            enddo
!            enddo
!            enddo
!  c
!
!
!end subroutine compute_force_constants




!subroutine  dyn_mat_index(unit_vect,tau,Rn,Rn_period,Rmax,NN_MAX,ndir,ndir2, neigh_type,nat_up,ndir_max,itab)
subroutine  dyn_mat_index(iconf)
use derived_types, only: config_real
implicit none
integer, INTENT(IN) :: iconf

!local
integer :: ia, ja, nn1, nn2, nn3, i1, i2, i3, icount
real(kind(0.d0)), parameter  :: zero=1.d-10
real(kind(0.d0))  :: Rtemp(3), norm
integer, dimension(:,:,:,:), allocatable :: itab



!Rmax2=2.0d0*Rmax



nn1=int(config_real(iconf)%nxCell/2)+1
nn2=int(config_real(iconf)%nyCell/2)+1
nn3=int(config_real(iconf)%nzCell/2)+1

if (allocated(itab)) deallocate(itab) ; allocate(itab(config_real(iconf)%nat, -nn1:nn1, -nn2:nn2, -nn3:nn3))


 do  ia=1,config_real(iconf)%nat
   icount=0
   do  ja=1,config_real(iconf)%nat
       do i1=-nn1,nn1
       do i2=-nn2,nn2
       do i3=-nn3,nn3
         Rtemp(:) = dble(i1)*config_real(iconf)%cell(:,1) + &
                     dble(i2)*config_real(iconf)%cell(:,2) + &
                     dble(i3)*config_real(iconf)%cell(:,3)

         Rtemp(:)=Rtemp(:) + config_real(iconf)%pos_cart(:,ja) - config_real(iconf)%pos_cart(:,ia)
         norm=DOT_PRODUCT(Rtemp(:),Rtemp(:))
         ! should be added the RCut ....
         if (norm>zero) then
            icount=icount+1
            !o Rn(ia,icount,:)=Rtemp(:)
	         itab(ja,i1,i2,i3)=icount
            !o neigh_type(ia,icount)=ja
            !o Rn_period(ia,icount,:)=Rn(ia,icount,:)+tau(:,ia)-tau(:,ja)
            config_real(iconf)%type_neigh(ia,icount)=config_real(iconf)%itype(ja)
            config_real(iconf)%kind_neigh(ia,icount)=ja
            config_real(iconf)%r_ij(ia,icount)=dsqrt(norm)
            config_real(iconf)%u_ij(ia,icount,:)=-Rtemp(:)!/config_real(iconf)%r_ij(ia,c)
	      else
	         itab(ja,i1,i2,i3)=0
         endif

       end do
       end do
       end do
   end do
   !o ndir2(ia)=icount
   config_real(iconf)%n_neigh(ia)=icount
end do
return
end subroutine dyn_mat_index
!*********************************************************************!
!*********************************************************************!





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
#if(PARAPH)
  use mpi
  use mod_mpi_phondy
#endif

  implicit none
  !-----------------------------------------------
  !   L o c a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: i, j
  real(double)::time
  integer :: count1,count2,count_rate,count_max
  integer :: nlocal,u,v

  call system_clock (count1,count_rate,count_max)


  call partial_take_size(i_start,i_final,imax)

#if(PARAPH)
  write(*,'("PHONDY:  proc i_start i_final natoms",i5,3i7,i10)') rangph,i_start,i_final,i_final-i_start,imax
#else
  write(*,*) 'PHONDY: i_start i_final imax ', i_start, i_final, imax
#endif

  nlocal=(i_final-i_start)*3
  if (i_final<=i_start) then
   if (rangph==0) write(*,*) 'PHONDY: i_final should be grater than i_start',i_final, ' > ', i_start
   if (rangph==0) write(*,*) 'PHONDY: stop in force_constat.F90. Change no of procs.'
   stop
  end if


  if (allocated(mlocal))  deallocate(mlocal)
  if (allocated(u_local)) deallocate(u_local)
  if (allocated(v_local)) deallocate(v_local)
  if (allocated(w))       deallocate(w)
  allocate (mlocal(imax),u_local(imax), v_local(imax),w(nmat))

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
       write(*,*) 'PHONDY: WARINIIIIIIIIIIIIG isave is set to', isave
       write(*,*) 'PHONDY: We will switch isave to 0'
      end if
     end if
   isave=0
  end if

#if PARAPH && PHONDY
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

   if (allocated(matfor)) deallocate(matfor)
   allocate(matfor(nmat,nmat))

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
  !call MPI_BARRIER(MPI_COMM_WORLD,codeph)
  !call MPI_BARRIER(MPI_COMM_WORLD,codeph)

#else
  if (iread == 0 ) then

    if (allocated(matfor)) deallocate(matfor)
    allocate(matfor(nmat,nmat))

    matfor(:,:)=0.d0
    do i=1,imax
      u=u_local(i)
      v=v_local(i)
      matfor(u,v)=mlocal(i)
    end do
  end if

  if (iread==1) then

    if (allocated(matfor)) deallocate(matfor)
    allocate(matfor(nmat,nmat))

    inamefile=100000
    write(namefile,'(i6)') inamefile
    namefile=TRIM(namefile)
    open(71,file=namefile//".m",status='unknown', access='sequential',form='unformatted')
    open(72,file=namefile//".u",status='unknown', access='sequential',form='unformatted')
    open(73,file=namefile//".v",status='unknown', access='sequential',form='unformatted')
    read(71) imax
    write(*,*) 'reading on proc  and imax ',rangph, imax
 !
    do i=1,imax
     read (72) u_local(i)
     read (73) v_local(i)
     read (71) matfor (u_local(i),v_local(i))
    end do
   close(71)
   close(72)
   close(73)
  end if   !iread==1



  if (isave==1) then

    inamefile=100000
    write(namefile,'(i6)') inamefile
    namefile=TRIM(namefile)
    open(71,file=namefile//".m",status='unknown', access='sequential',form='unformatted')
    open(72,file=namefile//".u",status='unknown', access='sequential',form='unformatted')
    open(73,file=namefile//".v",status='unknown', access='sequential',form='unformatted')
    write(71) imax
    write(*,*) 'writing on proc  and imax ',rangph, imax
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





  call system_clock (count2,count_rate,count_max)
  time=real((count2-count1))/real(count_rate)
  if (iread /=1) then
   if (rangph==0) write(*,"(' PHONDY: MATFOR  was filled in.......:  ',f16.8,' s')") time
   if (rangph==0) write(*,*) 'PHONDY: MATFOR', matfor (1,2) ,matfor(2,1)
   if (rangph==0) write(*,*) 'PHONDY: MATFOR', matfor (1,5) ,matfor(5,1)
!debug  if (rangph==0) write(*,*) 'PHONDY: MATFOR', matfor (1002,503) ,matfor(503,1002)
  end if
#endif

#if(PARAPH)

  call system_clock (count2,count_rate,count_max)
  time=real((count2-count1))/real(count_rate)
  if (rangph==0) write(*,"(' PHONDY: MATFOR  was filled in.......:  ',f16.8,' s')") time

 if (isave==2) then

 if (rangph==0) then
#endif
  do i=1,3*im
    do j=i,3*im

    matfor(i,j)=0.5d0*(matfor(i,j)+matfor(j,i))
    matfor(j,i)=matfor(i,j)
    end do
  end do
#if(PARAPH)
end if

  call system_clock (count1,count_rate,count_max)
  time=real((count1-count2))/real(count_rate)
  if (rangph==0) write(*,"(' PHONDY: MATFOR  was symmetrized in...:  ',f16.8,' s')") time


end if
#else

  call system_clock (count1,count_rate,count_max)
  time=real((count1-count2))/real(count_rate)
  if (rangph==0) write(*,"(' PHONDY: MATFOR  was symmetrized in...:  ',f16.8,' s')") time
#endif

  if (rangph==0) write(*,*) 'PHONDY: ...the force constants were  filled'
    !  call MPI_BARRIER(MPI_COMM_WORLD,codeph)

   if (rangph==0) write(*,*) 'PHONDY: MATFOR', u_local(1), v_local(1), mlocal(1)
   if (rangph==0) write(*,*) 'PHONDY: MATFOR', u_local(2), v_local(2), mlocal(2)
   if (rangph==0) write(*,*) 'PHONDY: MATFOR', u_local(3), v_local(3), mlocal(3)
   if (rangph==0) write(*,*) 'PHONDY: MATFOR', u_local(4), v_local(4), mlocal(4)

!debug   write(*,*) 'test rangph', rangph, imax
end subroutine force_constant
