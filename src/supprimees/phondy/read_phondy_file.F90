subroutine read_phondy_file
use T_kind_param_m, ONLY:  double
use gen_com_m, ONLY: lenfnam,fnam,angst,ev2erg,im
use var_pot, only: ipotentiel

use  phondy_in_ndm_module, ONLY: rangph, nu_min, nu_max, n_nu_up,   &
                                 nsitedos,isitedos, isave,iread,width_dos,lldos,dos,  &
                                 lmodes,leigenvectors,imodes,nmodes, wmodes, fnamtin_lammps, lq_points, l_thermo_atoms, &
                                 no_of_qpoints, q_point,weight_q_point, q_eigenvalues, debug_ph


implicit none

namelist /input_phondy/ nu_min, nu_max, n_nu_up,width_dos,lldos,lmodes,leigenvectors,wmodes,isave,iread, &
                        lq_points, l_thermo_atoms, debug_ph

character(len=128) :: fnamtin,fnamt_ldos,fnamt_lmodes,fnamtin_qpoints
integer :: luphondy
integer :: ii

   isave=0 ! isave=1 we will store dynamical matrix on HDD for a subequent diagonalization
           ! isave=0 nothing stored. Working on the fly.
           ! isave=2 just computing dyn matrix by MPI and diag by threading. Nothing stored.
   iread=0 ! iread=1 reading the dyn matrix from previous run
           ! iread=0 everything is computed. Working on the fly.


    debug_ph=.false.
    isave=0
    nu_min=0.1
    nu_max=20.0
    n_nu_up=1000
    width_dos=0.02
    lldos=.false.     ! if you want the local density of states
    lmodes=.false.    ! if you want nice jmol output
    l_thermo_atoms=.false.
    leigenvectors=.false.
    wmodes=1.d-3  ! the default localization of modes
    lq_points=.false.
 fnamtin = fnam(1:lenfnam)//'.phondy'
 fnamt_ldos = fnam(1:lenfnam)//'.phondy.ldos'
 fnamt_lmodes = fnam(1:lenfnam)//'.phondy.lmodes'
 fnamtin_lammps = fnam(1:lenfnam)//'.in.lammps'
 fnamtin_qpoints = fnam(1:lenfnam)//'.phondy.qpoints'
!debug  write(*,*) fnamtin
 luphondy=621
 open(unit=luphondy, file=fnamtin, status='unknown')
 read (luphondy, nml=input_phondy)
 close (luphondy)


 allocate (dos(n_nu_up))
 if (lldos) then
   open(unit=luphondy, file=fnamt_ldos, status='unknown')
   read(luphondy,*) nsitedos

   if (allocated(isitedos)) deallocate(isitedos)
   allocate (isitedos(nsitedos))

   do ii=1,nsitedos
    !
     read(luphondy,*)  isitedos(ii)
     if (isitedos(ii) > im) then
       if (rangph==0) then
        write(*,*) 'the LDOS cannot be projected on the atom', isitedos(ii)
        write(*,*) 'the value exceeds the number of atoms isitedos, im ',isitedos(ii),im
        write(*,*) 'stop in read_phondy_file.F90'
       end if
        stop
     end if
   !
   end do
 close(luphondy)
end if


if (lq_points) then

   open(unit=luphondy, file=fnamtin_qpoints, status='unknown')
   read(luphondy,*) no_of_qpoints
   if (allocated(q_point)) deallocate(q_point) ; allocate(q_point(3,no_of_qpoints))
   if (allocated(weight_q_point)) deallocate(weight_q_point) ; allocate(weight_q_point(no_of_qpoints))
   do ii=1,no_of_qpoints
      read(luphondy,*) q_point(1:3,ii), weight_q_point(ii)
   end do
   close(luphondy)


   if (.not.((ipotentiel == 20 ).or. (ipotentiel==10))) then
     if (rangph==0) write(6,*) 'PHONDY error: The dispersion relation (lq_points = .true.) are implemented only for EAM (dmtype=10) and ML (ipotential=20)'
     if (rangph==0) write(6,*) 'PHONDY error: lq_points and ipotentiel', lq_points, ipotentiel
     if (rangph==0) write(6,*) 'PHONDY error: For others type of forces please ask someone to imlplement it ...'
     stop  'read_phondy_file: incompatibilty between ipotentiel and lq_points '
   end if

   if (isave/=2) then

     if (rangph==0) write(6,*) 'PHONDY error: In the qpoints version  (lq_points = .true.) the isave should be set to 2'
     if (rangph==0) write(6,*) 'PHONDY error: Now the lq_points and isave ', lq_points, isave
     stop  'read_phondy_file: incompatibilty between  lq_points and isave '
   end if
end if

 if (lmodes) then
   open(unit=luphondy, file=fnamt_lmodes, status='unknown')
   read(luphondy,*) nmodes

   if (allocated(imodes)) deallocate(imodes)
   allocate (imodes(nmodes))

    if (nmodes > 3*im) then
      write(*,*) 'The max no of modes is ...:', 3*im
      write(*,*) 'The current value is .....:', nmodes
      write(*,*) 'stop in read_phondy_file.F90'
      stop
    end if


   do ii=1,nmodes
    !
     read(luphondy,*)  imodes(ii)
     if (imodes(ii) > 3*im) then
       if (rangph==0) then
        write(*,*) 'the densimodes of this mode cannot be projected on the atom', imodes(ii)
        write(*,*) 'the value exceeds the number of modes, 3*im ',imodes(ii),3*im
        write(*,*) 'stop in read_phondy_file.F90'
       end if
        stop
     end if
   !
   end do
 close(luphondy)
end if



#if (PARAPH)
if (rangph==0) then
        select case (isave)
           case (0)
                  write(6,*) 'PHONDY:  The force constant matrix is not stored.'
           case (1)
                  write(6,*) 'PHONDY:  The force constant matrix is stored  on  HDD for a subequent diagonalization'
           case (2)
                  write(6,*) 'PHONDY:  Just computed dyn matrix by MPI and then diag by threading. Nothing stored '
        end select


        select case (iread)
           case (0)
                  write(6,*) 'PHONDY:  Nothing to read. '
           case (1)
                  write(6,*) 'PHONDY:  Reading the dyn matrix from the previous run. '
        end select

end if

#endif

if (rangph==0) then
    write(*,'("PHONDY: isave .......................................:",I6)') isave
    write(*,'("PHONDY: iread .......................................:",I6)') iread

   if ((isave==1).and.(iread==1)) then
    write(*,'("this case is not LOGICAL. You cannot save and read the DM in the same time")')
    write(*,'("isave=1 and iread=0 - means saving")')
    write(*,'("isave=0 and iread=1 - means reading and diagonalizing previous saving")')
    write(*,'("isave=0 and iread=0 - means computing on the fly and diagonalizing ")')
    write(*,'("Please choose one of this options ")')
    write(*,'("STOP in read_phondy_file ")')
    stop
   end if

end if





end subroutine  read_phondy_file





subroutine print_phondy(rangloc)
implicit none
integer, intent(in) :: rangloc
character(len=1) :: quote,dquote


 quote=char(39)
dquote=char(34)

if (rangloc==0) then


write(6,'("/---------------------------------------------\")')
write(6,'("                                               ")')
write(6,'("            _                     _            ")')
write(6,'("           | |                   | |           ")')
write(6,'("      _ __ | |__   ___  _ __   __| |_   _      ")')
write(6,'("     | ",a,"_ \| ",a,"_ \ / _ \| ",a,"_ \ / _",a," | | | |     ")')quote,quote,quote,quote
write(6,'("     | |_) | | | | (_) | | | | (_| | |_| |     ")')
write(6,'("     | .__/|_| |_|\___/|_| |_|\__,_|\__, |     ")')
write(6,'("     | |                             __/ |     ")')
write(6,'("     |_|                            |___/      ")')
write(6,'("                                               ")')
write(6,'("copyright mihai-cosmin.marinica@cea.fr         ")')
write(6,'("contributions:                                 ")')
write(6,'("C. Lapointe, MCM                               ")')
write(6,'("\---------------------------------------------/")')

end if

end subroutine print_phondy
