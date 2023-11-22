module newunit_mod
contains
  subroutine newunit(nunit)
    integer, intent(out) :: nunit
    ! local
    integer, parameter :: LUN_MIN=500, LUN_MAX=999
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

end module newunit_mod


  
