module storage
contains
  !> ART store
  !!    This subroutine stores the configurations at minima
  !!    and activated points.
  !!    By definition, it uses pos, box and scala (actually x, y , z -> pos
  subroutine store( fname )

    use defs
    use atomconfig,only: atom_config
    use cellconfig,only:cell_config
    use boxconfig,only:box_config
    use ndm2art2ndm,only:art2ndm,atcfart,boxart
    use gen_com_m,only:lwgin
    USE rasmolT_mod,only: rasmolT

  
    implicit none

    !Arguments
    character(len=*), intent(in)  :: fname

    !Local variables
    integer :: i, ierror
    real(kind=8),dimension(3) :: boxl
    character(len=*), parameter :: extension = ".xyz"
    character(len=60) :: fnamexyz,namemol
    character(len=4), dimension(natoms) :: frzchain
    type(atom_config)::ataux
    type(box_config)::boxaux
    type(cell_config)::celaux


    
    ! Update the box size
    boxl = box * scala


    ! Added by Fedwa El-Mellouhi July 2002,
    ! writes the configuration in .xyz format.
    ! Modified by E. Machado-charry for v_sim and BigDFT.
    if ( write_xyz ) then

       ! If there is a constraint over a given atom,
       ! it is written in the geometry file.
       do i = 1, NATOMS
          if      ( constr(i)== 0) then
             frzchain(i)='    '
          else if (constr(i) == 1) then
             frzchain(i)='   f'
          else if (constr(i) == 2) then
             frzchain(i)='  fy'
          else if (constr(i) == 3) then
             frzchain(i)=' fxz'
          end if
       end do

       fnamexyz = trim(fname) // extension
       write(unit6P,*) ' Writing to file : ', fnamexyz
       open(unit=XYZ, file=fnamexyz, status='unknown', &
            action='write', iostat=ierror)
       write(XYZ,*) NATOMS
       write(unit6P,*) 'boxl and scala are: ', (boxl(i),i=1,3), scalaref
       if (boundary == 'P') then
          write(XYZ,'(a,3(1x,1p,e24.17,0p))')'Periodic',  (boxl(i),i=1,3)
          ! write(XYZ,'(A,9(1X,F8.4),1X,A,1X,A)') 'Lattice="',boxl(1),0.0,0.0 &
          !                                               &,0.0,boxl(2),0.0 &
          !                                               &,0.0,0.0,boxl(3), '"'
       else if (boundary == 'S') then
          write(XYZ,'(a,3(1x,1p,e24.17,0p))')'Surface',   (boxl(i),i=1,3)
       else if (boundary == 'T') then
          ! write(XYZ,'(a,3(1x,1p,e24.17,0p))') 'Triclinic', cell(:,1)
          ! write(XYZ,'(3(e24.17))') cell(:,2)
          ! write(XYZ,'(3(e24.17))') cell(:,3)
          write(XYZ,'(a,9(1x,1p,e24.17,0p),a)')'Lattice="', cell(1,:), cell(2,:), cell(3,:),'"'
       else
          write(XYZ,*)'Free'
       end if

       do i=1, NATOMS
          if (constr(i)== 0) then
             write(XYZ,'(i2,3(1x,es23.16))') typat(i), x(i), y(i), z(i)
          else
             write(XYZ,'(i2,3(1x,es23.16),1x,a4)') typat(i), &
                  x(i), y(i), z(i), frzchain(i)
          end if
       end do

       close(XYZ)

    else

       write(unit6P,*) ' Writing to file : ', fname

       open(unit=FCONF, file=fname, status='unknown', &
            action='write', iostat=ierror)
       write(FCONF,*) 'run_id: ', mincounter
       write(FCONF,*) 'total_energy: ', total_energy

       if ( boundary == 'T' ) then
          write(FCONF,*) boundary, cell(:,1)
          write(FCONF,*) ' ', cell(:,2)
          write(FCONF,*) ' ', cell(:,3)
       else
          write(FCONF,*) boundary, boxl
       end if

       do i=1, NATOMS
          write(FCONF,'(i2,3(1x,es23.16))') typat(i), x(i), y(i), z(i)
       end do

       close(FCONF)

    end if

    if (ivisuart.gt.0) then
       call art2ndm(ataux,boxaux,celaux,linit=.true.)
!!$       if (index(fname, 'min').ne.0) then
!!$          namemol='CONF'//trim(fname(4:))//'2'
!!$       end if
!!$       if (index(fname, 'sad').ne.0) then
!!$          namemol='CONF'//trim(fname(4:))//'1'
!!$       end if
!!$       write(unit6P,*)'name MOL',namemol
!!$       call rasmolT(ataux,boxaux,namefr=namemol,latcomp=.true.,ivisumol=ivisuart)
       namemol=trim(fname)
       call rasmolT(atcfart,boxart,namefr=namemol,latcomp=.true.,ivisumol=ivisuart)
       if (lwgin)       call rasmolT(atcfart,boxart,namefr=namemol,latcomp=.true.,ivisumol=5)
    end if
       
    
  END SUBROUTINE store

  !> ART newunit
  !! This function returns the lowest available unit to open a file
  !! http://www.fortran90.org/src/best-practices.html#file-input-output
  !! http://fortranwiki.org/fortran/show/newunit
  subroutine newunit(nunit)
    integer, intent(out) :: nunit
    ! local
    integer, parameter :: LUN_MIN=500, LUN_MAX=1000
    logical :: opened
    integer :: lun
    ! begin
    nunit=-1
    do lun=LUN_MIN,LUN_MAX
       inquire(unit=lun,opened=opened)
       if (.not. opened) then
          nunit=lun
          exit
       end if
    end do
  end subroutine newunit

end module storage

