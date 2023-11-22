module init_all_atoms_mod
  use update_invcell_mod,only: update_invcell
contains
  subroutine init_all_atoms( nat, typa, posa, const_, boxl, boxtype, nproc_, me_, inputfile )

  use defs, only : FREFCONFIG, cell, invcell,unit6P

  implicit none
  character(len=20) :: dummy
  integer           :: ierror, nproc
  logical           :: flag
  real(kind=8), dimension(:), pointer :: xa, ya, za

  !Arguments
  integer,      intent(in)               :: nat
  integer, intent(out)                                :: typa(nat)
  real(kind=8), target, intent(inout)                    :: posa(3*nat)
  integer, intent(out)                    :: const_(nat)
  real(kind=8), dimension(3), intent(out) :: boxl
  character(len=1), intent(out)           :: boxtype
  integer,      intent(in)                :: nproc_
  integer,      intent(in)                :: me_
  character(len=*), intent(in)            :: inputfile

  !Local variables
  integer                                 :: i
  character(len=2)                        :: symbol
  character(len=10)                       :: name

  !_______________________

  nproc = nproc_
  cell = 0.0d0
  invcell = 0.0d0

  ! Read the atomic positions
  inquire(file = inputfile, exist = flag )
  if (.not.flag) then
     write(unit6P,*) "you have not given an initial configuration"
     write(unit6P,*) "the program will stop"
     stop
  else
     write(unit6P,*) "initializing atomic positions with: ", inputfile
  endif

  xa => posa(1:nat)
  ya => posa(nat+1:2*nat)
  za => posa(2*nat+1:3*nat)

  open(unit=FREFCONFIG,file=inputfile,status='old',action='read',iostat=ierror)
  read(FREFCONFIG, '(A10,i5)' ) dummy
  if (inputfile(len(trim(inputfile))-3:) .ne. '.xyz') &
       read(FREFCONFIG,*) dummy  !! M-A Malouin
  read(FREFCONFIG,*) boxtype, boxl(1), boxl(2), boxl(3)
  if (boxtype=='T') then
     cell(:,1) = boxl(:)
     read(FREFCONFIG,*) cell(:,2)
     read(FREFCONFIG,*) cell(:,3)
     call update_invcell( )
     write(unit6P,*)  "Triclinic box info " 
     write(unit6P,*)  cell(:,1)
     write(unit6P,*)  cell(:,2)
     write(unit6P,*)  cell(:,3)
  end if

  do i = 1, nat
     read(FREFCONFIG,*) typa(i),xa(i),ya(i),za(i)
  end do
  close(FREFCONFIG)

  const_ = 0
END SUBROUTINE init_all_atoms


end module init_all_atoms_mod
