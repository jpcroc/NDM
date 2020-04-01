subroutine write_dos_total
use phondy_in_ndm_module, only: rangph, isave
use thermo_module, only: allocate_thermo
implicit none
! passing argumets ....
integer :: iproc, i
real(kind=8) :: temps1,temps2


call allocate_thermo

if ((isave==0).or.(isave==2)) then


  if (isave==2) then
    if (rangph==0) call write_dos !() this will be in common for lq_points and gamma
  else
    call write_dos
  end if

end if


return
end subroutine write_dos_total



! ************************************************
!           Sous-programme dos.F90
!           MCM, CL  - 10/2019 update from 08/2012
! ************************************************


subroutine write_dos_ldos
use phondy_in_ndm_module, only: W, nmat, rangph, isave,lq_points,leigenvectors, lmodes, l_thermo_atoms, lldos, unit_nu
use thermo_module, only: allocate_thermo
use diago_scalapack_real, only : Z, DESCZ
use mpi
implicit none
! passing argumets ....
integer :: iproc, i
real(kind=8) :: temps1,temps2

interface
        subroutine write_ldos(iproc,N,W)
          use diago_scalapack_real, only : Z, DESCZ
          implicit none
          integer :: N
          integer :: iproc
          !integer :: DESCZ(:)
          real(kind=8) :: W(:)
          !real(kind=8) :: Z(:,:)
        end subroutine write_ldos

        !
        subroutine write_lmodes(iproc,N,W)
          use diago_scalapack_real, only : Z, DESCZ
          implicit none
          integer :: N
          integer :: iproc
          !integer :: DESCZ(:)
          real(kind=8) :: W(:)
          !real(kind=8) :: Z(:,:)
        end subroutine write_lmodes


        subroutine write_thermo_atoms(iproc,N,W)
          use diago_scalapack_real, only : Z, DESCZ
          implicit none
          integer :: N
          integer :: iproc
          !integer :: DESCZ(:)
          real(kind=8) :: W(:)
          !real(kind=8) :: Z(:,:)
        end subroutine write_thermo_atoms
end interface

call allocate_thermo

if ((isave==0).or.(isave==2)) then
  if (isave==2) then
    if (rangph==0) call write_dos !() this will be in common for lq_points and gamma
  else
    call write_dos
  end if

  return
  if (isave==2) stop 'in dos waiting for further developements ...'
end if

 if (lq_points) then

      if (rangph==0) write(*,'("PHONDY: not yet a implementation for q_points")')

  else

    if (leigenvectors) then
      temps1=MPI_Wtime()
      !This is for LDOS .....
      if (lldos.or.lmodes.or.l_thermo_atoms) call allocate_thermo
      if (lldos) then   !if ldos=.true. compute the Local DOS
        !
        call write_ldos(rangph, nmat,W)  ! all others variables are sended by the module phony_in_ndm_module
        !
      end if !lldos
      temps2=MPI_Wtime()


      if (lmodes) then
        call write_lmodes(rangph,nmat,W)
      end if !lmodes
      temps2=MPI_Wtime()


      if (l_thermo_atoms) then
        call write_thermo_atoms(rangph,nmat,W)
        !call write_thermo_atoms(N,Z,DESCZ)
      end if !lmodes
      temps2=MPI_Wtime()
    end if !leigenvectors

    if (rangph==0) then
     print*,'That ldos time  ',temps2-temps1,' seconds'
    end if
end if

return
end subroutine write_dos_ldos










subroutine write_dos()
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use var_pot
  use phondy_in_ndm_module
  use thermo_module
  implicit none
  integer :: n,jmat
  real(double) :: nu, dnu, beta_dos,fermi,delta_fermi
  real(double) :: fermi_func, fermi_func_deriv
!for future developement (Q-points)
  real(double) :: eigenvalues(nmat)

  real(double) :: kBT, kBT2,eig, xp2
  logical      :: lskip
!  real(double), allocatable, dimension(:) :: Fmin2, Smin2
  integer :: nT, itest_mode, ik
character (len=60) :: CHFMT


!call allocate_thermo()
  if (rangph==0) then


    if (.not.lq_points) then
      no_of_qpoints=1
      if (allocated(q_eigenvalues)) deallocate(q_eigenvalues) ; allocate(q_eigenvalues(nmat, no_of_qpoints))
      if (allocated(weight_q_point)) deallocate(weight_q_point) ; allocate(weight_q_point(no_of_qpoints))
      weight_q_point(:) = 1.d0
      eigenvalues(:)=sign(1.d0,W(:))*dsqrt(dabs(W(:)))*unit_nu*1.D3
      q_eigenvalues(:,1) = eigenvalues(:)
    end if


    open(unit=711,file='dispersion.dat',status='unknown')
    write(CHFMT,*)'(i6, 1x, ',  size(q_eigenvalues(:,:), DIM=1)   ,'f20.10)'
    do n=1,size(q_eigenvalues(:,:), DIM=2)
      write(711, CHFMT) n,  q_eigenvalues(:,n)
    end do
    close(711)

    open(unit=712,file='eigenvalues.dat',status='unknown')
    write(CHFMT,*)'(i6, 1x, ',  size(q_eigenvalues(:,:), DIM=2)   ,'f20.10)'
    do n=1,size(q_eigenvalues(:,:), DIM=1)
      write(712, CHFMT) n,  q_eigenvalues(n,:)
    end do
    close(712)



    beta_dos=1.0D0/width_dos
    dnu=(nu_max-nu_min)/(n_nu_up-1)
    dos(:) = 0.d0
    do ik=1,no_of_qpoints
      eigenvalues(:) = q_eigenvalues(:,ik)
      do n=1,n_nu_up
        nu=nu_min+dble(n-1)*dnu
        do jmat=1,nmat
          fermi=fermi_func(eigenvalues(jmat)-nu,beta_dos)
          delta_fermi=fermi_func_deriv(eigenvalues(jmat)-nu,beta_dos)
          dos(n)=dos(n)+delta_fermi*weight_q_point(ik)
        end do
      end do
    end do


    open(unit=711,file='dos.dat',status='unknown')
    do n=1,n_nu_up
      nu=nu_min+dble(n-1)*dnu
      write(711,*) nu, dos(n)/dble(nmat)
    end do
    close(711)

    ! some thermo is here
    Fmin(:)=0.d0
    Smin(:)=0.d0
    Fmin_cla(:)=0.d0
    Smin_cla(:)=0.d0



    do ik=1,no_of_qpoints
      eigenvalues(:) = q_eigenvalues(:,ik)
      itest_mode=0
      do jmat=1,nmat
        eig=eigenvalues(jmat)

        lskip=.false.
        if ( (eig < 0).and.(dabs(eig)>1.d-4)) then
          write(*,*) 'WARNING negative frequencies: unrelaxed structure or saddle point.'
          lskip=.true.
        end if

        if (dabs(eig) < 1.d-4 ) then
          itest_mode=itest_mode+1
          lskip=.true.
        end if


        do nT=1,ntemp
          temperature=temp_min + dtemp*dble(nT-1)
          kBT=temperature*temperature_to_ev
          kBT2=2.0d0*kBT
          xp2=eig*thz_to_ev/kBT2
          if (.not.(lskip)) then
            Fmin(nT)=Fmin(nT)+kBT*DLOG(2.0d0*DSINH(xp2) )*weight_q_point(ik)
            !Smin(nT)=Smin(nT)+temperature_to_ev*(xp2/DTANH(xp2)-DLOG(2.0d0*DSINH(xp2))) !*weight
            Smin(nT)=Smin(nT)+(xp2/DTANH(xp2)-DLOG(2.0d0*DSINH(xp2)))*weight_q_point(ik)

            Fmin_cla(nT)=Fmin_cla(nT)+kBT*DLOG(2.0d0*xp2 )*weight_q_point(ik)
!           Smin_cla(nT)=Smin_cla(nT)-temperature_to_ev*DLOG(2.d0*xp2)+temperature_to_ev !*weight
            Smin_cla(nT)=Smin_cla(nT)-(DLOG(2.d0*xp2)+1.0d0)*weight_q_point(ik)
          end if
        end do
      end do  !jmat
      if (itest_mode > 3) then
        write(*,*) 'WARNING more that 3 translation modes. Maybe soft modes. Unreliable thermo! ik point', ik
      end if
    end do ! ik

    open (67, file='thermo_global.dat')
    do nT=1,ntemp
      temperature=temp_min + dtemp*dble(nT-1)
      write(67,'(f12.5, 4e25.15)') temperature, Fmin(nT), Fmin_cla(nT), Smin(nT), Smin_cla(nT)
    end do
    close(67, status='keep')



  end if
return
end subroutine write_dos




subroutine write_thermo_atoms(iproc,N,W)
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use var_pot
  use phondy_in_ndm_module, only: nmat, unit_nu
  use thermo_module
  use diago_scalapack_real, only : Z, DESCZ
  implicit none

  ! passing argumets ....
  integer :: N
  integer :: iproc
  !integer :: descz(:)
  !real(kind=8) ::  Z(:,:)
  real(kind=8) ::  W(:)
  ! passing argumets ....

  integer :: jmat, i_loc, jj_mat
  real(double) :: nu, dnu, beta_dos,fermi,delta_fermi
  real(double) :: fermi_func, fermi_func_deriv
  !for future developement (Q-points)
  real(double) :: eigenvalues(nmat), eigenvector1(nmat),eigenvector2(nmat),eigenvector3(nmat)

  real(double) :: kBT, kBT2,eig, xp2, ww
  !debug real(double) :: ww_TOT,stemp
  logical      :: lskip
  integer :: nT, itest_mode


  eigenvalues(:)=sign(1.d0,w(:))*dsqrt(dabs(w(:)))*unit_nu*1.D3
  open (67, file='thermo_atoms.dat')

  Fmin_at(:)=0.d0
  Smin_at(:)=0.d0
  Fmin_cla_at(:)=0.d0
  Smin_cla_at(:)=0.d0

  !debug ww_TOT=0.d0
  do i_loc = 1, nat
    do jj_mat = 1, nmat
      call pdelget('A',' ',eigenvector1(jj_mat),Z,1+(i_loc-1)*3,jj_mat,descz)
      call pdelget('A',' ',eigenvector2(jj_mat),Z,2+(i_loc-1)*3,jj_mat,descz)
      call pdelget('A',' ',eigenvector3(jj_mat),Z,3+(i_loc-1)*3,jj_mat,descz)
    end do

    do jj_mat = 1, nmat
      itest_mode=0
      eig=eigenvalues(jj_mat)

      lskip=.false.
      if ( (eig < 0).and.(dabs(eig)>1.d-4)) then
        write(*,*) 'WARNING negative frequencies: unrelaxed structure or saddle point.'
        lskip=.true.
      end if

      if (dabs(eig) < 1.d-4 ) then
        itest_mode=itest_mode+1
        lskip=.true.
      end if

      !ww =  (DOT_PRODUCT(eigenvector1, eigenvector1) + &
      !      DOT_PRODUCT(eigenvector2, eigenvector2) + &
      !      DOT_PRODUCT(eigenvector3, eigenvector3))/dble(nmat)

      ww =  eigenvector1(jj_mat)**2  + eigenvector2(jj_mat)**2 + eigenvector3(jj_mat)**2
      !debug ww_TOT = ww_TOT + ww

      do nT=ntemp,ntemp
        temperature=temp_min + dtemp*dble(nT-1)
        kBT=temperature*temperature_to_ev
        kBT2=2.0d0*kBT
        xp2=eig*thz_to_ev/kBT2
        if (.not.(lskip)) then
          !if (rangph==0) write(*,*) temperature, i_loc,  ww, DOT_PRODUCT(eigenvector1, eigenvector1), DOT_PRODUCT(eigenvector2, eigenvector2), DOT_PRODUCT(eigenvector3, eigenvector3), Smin_cla_at(i_loc)
          Fmin_at(i_loc)=Fmin_at(i_loc)+kBT*ww*DLOG(2.0d0*DSINH(xp2) ) !*weight
          !Smin(nT)=Smin(nT)+temperature_to_ev*(xp2/DTANH(xp2)-DLOG(2.0d0*DSINH(xp2))) !*weight
          Smin_at(i_loc)=Smin_at(i_loc)+ww*(xp2/DTANH(xp2)-DLOG(2.0d0*DSINH(xp2)))

          Fmin_cla_at(i_loc)=Fmin_cla_at(i_loc)+ww*kBT*DLOG(2.0d0*xp2 ) !*weight
          !Smin_cla(nT)=Smin_cla(nT)-temperature_to_ev*DLOG(2.d0*xp2)+temperature_to_ev !*weight
          Smin_cla_at(i_loc)=Smin_cla_at(i_loc)-ww*(DLOG(2.d0*xp2)-1.0d0) ! *weight)
        end if
      end do !ntemp
    end do  ! jj_mat
      !debug stemp = stemp + Smin_cla_at(i_loc)
      !debug if(rangph==0) write(*,*) i_loc, ww, ww_TOT, stemp
  end do ! i_loc

  do i_loc=1,nat
      !do nT=ntemp,ntemp
      !   temperature=temp_min + dtemp*dble(nT-1)
    write(67,'(i8, 4e25.15)') i_loc, Fmin_at(i_loc), Fmin_cla_at(i_loc), Smin_at(i_loc), Smin_cla_at(i_loc)
      !end do
  end do

  close(67, status='keep')

return
end subroutine write_thermo_atoms

!-----------------------------------------------------------------------
!
  subroutine write_ldos (iproc,N,W)
    use T_kind_param_m, ONLY:  double
    use phondy_in_ndm_module, ONLY : unit_nu,nu_min,nu_max,n_nu_up,nsitedos,width_dos,isitedos,weight
    use thermo_module
    use diago_scalapack_real, only : Z, DESCZ
    implicit none
! passing argumets ....
    integer :: N
    integer :: iproc
    !integer :: descz(:)
    !real(kind=8) ::  Z(:,:)
    real(kind=8) ::  W(:)
    real(kind=8) :: nu
! LDOS and local variables .....
    real(kind=8) :: beta_dos,delta_fermi, fermi_func, fermi, fermi_func_deriv, eigenvalues(N)
    real(kind=8) :: eigenvector(N)
    real(kind=8) :: dnu,ww
    real(kind=8), allocatable, dimension(:,:,:) :: dos_loc
    real(kind=8), allocatable, dimension(:,:) :: dos_loc_tot
    integer :: isite,iatom,inu,ialpha,imat,jmat,inamefile
    character*7 :: namefile
    real(double) :: kBT, kBT2,eig, xp2
    logical      :: lskip
!  real(double), allocatable, dimension(:) :: Fmin2, Smin2
    integer :: nT, itest_mode

    eigenvalues(:)=sign(1.d0,w(:))*dsqrt(dabs(w(:)))*unit_nu*1.D3

    if (allocated(dos_loc)) deallocate(dos_loc)
    if (allocated(dos_loc_tot)) deallocate(dos_loc_tot)
    allocate (dos_loc(3,nsitedos,n_nu_up),dos_loc_tot(nsitedos,n_nu_up))
    dos_loc(:,:,:)=0.d0
    dos_loc_tot(:,:)=0.d0

    beta_dos=1.0D0/width_dos
    dnu=(nu_max-nu_min)/(n_nu_up-1)

    do isite=1,nsitedos
       iatom=isitedos(isite)
      do ialpha=1,3
        imat=ialpha+ (iatom-1)*3
        do jmat = 1, N
         call pdelget('A',' ',eigenvector(jmat),Z,imat,jmat,descz)
        enddo
        do inu=1,n_nu_up
          nu=nu_min+dnu*dble(inu-1)
            do jmat=1,N
              fermi=fermi_func(eigenvalues(jmat)-nu,beta_dos)
              delta_fermi=fermi_func_deriv(eigenvalues(jmat)-nu,beta_dos)
              ! ww=MAT(imat,jmat)*DCONJG(MAT(imat,jmat))
              ww=eigenvector(jmat)*eigenvector(jmat)
              dos_loc(ialpha,isite,inu)=dos_loc(ialpha,isite,inu)+ww*delta_fermi*weight
            end do !
            dos_loc_tot(isite,inu)=dos_loc_tot(isite,inu)+dos_loc(ialpha,isite,inu)
        end do     ! inu
      end do       !ialpha
    end do         !isite


 if (iproc==0) then
    do isite=1,nsitedos
       iatom=isitedos(isite)
       inamefile=1000000+iatom
       write(namefile,'(i7)') inamefile
       namefile=TRIM(namefile)
       open(612,file="ldos_"//namefile//".dat",status='unknown',access='sequential',form='formatted')
        do inu=1,n_nu_up

         nu=nu_min+dnu*(inu-1)
         write(612,'(f12.6,4E15.8)') nu, dos_loc_tot(isite,inu)/N,dos_loc(1,isite,inu)/N,dos_loc(2,isite,inu)/N,dos_loc(3,isite,inu)/N
        end do
       close (612,status='keep')
!! some LOCAL  thermo is here ....untested - 17/10/2017 - !!!!!

!       Fmin_loc(:,:)=0.d0
!       Smin_loc(:,:)=0.d0
!       Fmin_cla_loc(:,:)=0.d0
!       Smin_cla_loc(:,:)=0.d0

!       open (617, file='lthermo_'//namefile//'.dat',status='unknown',access='sequential',form='formatted')
!       do nT=ntemp,ntemp
!          temperature=temp_min + dtemp*dble(nT-1)
!          kBT=temperature*temperature_to_ev
!          kBT2=2.0d0*kBT
!          do inu=2,n_nu_up
!             nu=nu_min+dnu*dble(inu-1)
!             xp2=nu*thz_to_ev/kBT2
!             Fmin_loc(isite,nT)=Fmin_loc(isite,nT)+kBT*DLOG(2.0d0*DSINH(xp2))*3.d0*dos_loc_tot(isite,inu)*dnu !*weight
!             Smin_loc(isite,nT)=Smin_loc(isite,nT)+(xp2/DTANH(xp2)-DLOG(2.0d0*DSINH(xp2)))*dos_loc_tot(isite,inu)*dnu
!            Fmin_cla_loc(isite,nT)=Fmin_cla_loc(isite,nT)+kBT*DLOG(2.0d0*xp2)*3.d0*dos_loc_tot(isite,inu)*dnu !*weight
!             Smin_cla_loc(isite,nT)=Smin_cla_loc(isite,nT)-(DLOG(2.d0*xp2)+1.0d0)*3.d0*dos_loc_tot(isite,inu)*dnu
!          end do
!        write(617,'(f12.5, 4e25.15)') temperature, Fmin(nT), Fmin_cla(nT), Smin(nT), Smin_cla(nT)
!       end do
!       close (617,status='keep')
    end do
end if    !iproc
 return
 end subroutine write_ldos

subroutine write_lmodes(iproc,N,W)
 USE T_kind_param_m, ONLY:  double
 use gen_com_m, ONLY: im,imm,angst,at
 use tab_imm_m, ONLY :xp
 use phondy_in_ndm_module, ONLY : nmodes,imodes,wmodes
use diago_scalapack_real, only : Z, DESCZ
  implicit none
  integer :: N
  integer :: iproc
  !integer :: DESCZ(:)
  real(kind=8) :: W(:)
  !real(kind=8) :: Z(:,:)

  integer :: inu,ialpha,jmat,imat
  real(double) :: ww,eigentemp(3),temp(3)
  integer :: i,isite
  character*7 :: namefile(nmodes)
  integer     :: inamefile, icnt(nmodes)


! if (iam==0) then

!!!!!!!!!!THIS IS NOT CORRECT!!!!!!!!!!!!!!!!
    do inu=1,nmodes
       jmat=imodes(inu)
       inamefile=1000000+jmat
       write(namefile(inu),'(i7)') inamefile
       namefile(inu)=TRIM(namefile(inu))
       if (iproc==0)  open(546,file="lmodes_"//namefile(inu)//".temp",status='unknown',access='sequential',form='formatted')
       icnt(inu)=0
       do isite=1,im
        ww=0.d0
        do ialpha=1,3
          imat=ialpha+ (isite-1)*3
          call pdelget('A',' ',eigentemp(ialpha),Z,imat,jmat,descz)
          ww=ww+eigentemp(ialpha)*eigentemp(ialpha)
        end do
!debug   if (iproc==0) write(778,'(i6,f12.6,3d12.3)') isite,ww,eigentemp(1:3)
        if (ww>=wmodes) then
            icnt(inu)=icnt(inu)+1
            if (iproc==0 )write(546,'(3d18.8,3d17.6)') xp(1:3,isite),eigentemp(1:3)
        end if
       end do   !isite 1,im
     if (iproc==0)  close (546, status='keep')
     end do     !inu   1,imodes



   if (iproc==0 )  then
    write(*,'("PHONDY: ", i7," modes are treated")') nmodes
    write(*,'("PHONDY:       no    mode name_file no_of_atoms")')

    do inu=1,nmodes
     write(*,'("PHONDY:  ",  i7,i7," ",a," ",i7)') inu, imodes(inu),namefile(inu),icnt(inu)
     open(546,file="lmodes_"//namefile(inu)//".temp",status='old',access='sequential',form='formatted')
     open(548,file="modes_"//namefile(inu)//".xyz",status='unknown',access='sequential',form='formatted')

     write(548,'(i7)') icnt(inu)+8
     write(548,'("      ")')

     do i =1,icnt(inu)
      read (546,'(3d18.8,3d17.6)')  temp(1:3),eigentemp(1:3)
      write(548,'("W   ", 3d18.8,3d17.6)') temp(1:3)*angst,eigentemp(1:3)
     end do

     write(548,'("H   ", 6d17.6)') 0.d0,0.d0,0.d0,0.d0,0.d0,0.d0
     write(548,'("H   ", 6d17.6)') at(1:3,1)*angst, 0.d0,0.d0,0.d0
     write(548,'("H   ", 6d17.6)') at(1:3,2)*angst, 0.d0,0.d0,0.d0
     write(548,'("H   ", 6d17.6)') at(1:3,1)*angst+at(1:3,2)*angst, 0.d0,0.d0,0.d0

     write(548,'("H   ", 6d17.6)') at(1:3,3)*angst, 0.d0,0.d0,0.d0
     write(548,'("H   ", 6d17.6)') at(1:3,1)*angst+at(1:3,3)*angst, 0.d0,0.d0,0.d0
     write(548,'("H   ", 6d17.6)') at(1:3,2)*angst+at(1:3,3)*angst, 0.d0,0.d0,0.d0
     write(548,'("H   ", 6d17.6)') at(1:3,1)*angst+at(1:3,2)*angst+at(1:3,3)*angst, 0.d0,0.d0,0.d0
      close(546,status='delete')
      close(548,status='keep')
     end do

    end if
!!!!!!!!!I SHOULD CORRECT THAT!!!!!!!!!!!!!!!

! end if !iam==0




end subroutine write_lmodes


!****|******************************************************************|
        double precision FUNCTION FERMI_FUNC(E,beta)
        USE T_kind_param_m, ONLY:  double
        implicit none
        real(double) :: E,beta,xx

              xx=beta*E
              if(xx.lt.-100.D0) then
                    fermi_func=1.0D0
               elseif(xx.gt.100.D0) then
                   fermi_func=0.0D0
               else
                   fermi_func=1.D0/( 1.D0+EXP(xx) )
              endif
        RETURN
        END

!****|******************************************************************|
        double precision FUNCTION FERMI_FUNC_DERIV(E,beta)
        USE T_kind_param_m, ONLY:  double
        implicit none
	real(double) :: E,beta,xx

              xx=beta*E
              if(xx.lt.-100.D0) then
              fermi_func_deriv=0.0
               elseif(xx.gt.100.D0) then
              fermi_func_deriv=0.0
               else
              fermi_func_deriv=beta*EXP(xx)/( 1.D0+EXP(xx) )**2
              endif
        RETURN
        END
