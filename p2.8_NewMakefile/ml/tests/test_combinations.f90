subroutine read_database
use ml_in_ndm_module, only : rangml, db_file
implicit none
integer :: n_db_class
logical :: ok
inquire (file=db_file, exist=ok)
if (ok) then
  open (file=db_file, unit=50, action='read')
  !reading how many lines are in db_file
  call read_n_db_class(50,n_db_class)
  !allocate and fill the db_model object
  call read_db_file(50,n_db_class)
else
  if (rangml==0) write(6,*) 'ML: Fatal Error the database file is not there. The present name: ', db_file
  stop "read_database"
endif
close(50)

!read the files of database
!call read_files_from_db

return
end subroutine read_database

subroutine read_n_db_class(inp,nlines)
implicit none
integer, intent(in) :: inp
integer, intent(out) :: nlines
character(len=80) :: ctmp
integer :: io

nlines = 0
DO
  read(inp,*,iostat=io) ctmp
  !debug write(*,*) ctmp(1:1)
  if (ctmp(1:1)=='#') cycle
  if (io/=0) exit
  nlines = nlines + 1
END DO
rewind(inp)
return
end subroutine read_n_db_class

subroutine read_db_file(inp,n_db_class)
use ml_in_ndm_module, only : rangml, debug
use derived_types, only: db_model
implicit none
integer, intent(in) :: inp, n_db_class
character(len=80) :: ctmp
integer :: io,i

if (allocated(db_model)) deallocate(db_model) ; allocate(db_model(n_db_class))
i=0
DO
  read(inp,'(a)',iostat=io) ctmp
  write(*,*) ctmp(1:1), ctmp
  if (ctmp(1:1)=='#') cycle
  if (io/=0) exit
  backspace(inp)
  i = i + 1
  read(inp,*) db_model(i)%class, db_model(i)%klm, db_model(i)%no_total, db_model(i)%no_selec
END DO
rewind(inp)


if (i .ne. n_db_class) then
  if (rangml==0) write(6,*) 'Problems in readind db_file',inp
  if (rangml==0) write(6,*) 'These values should be equal', n_db_class, i
  stop "read_db_file"
end if

if (debug) then
   if (rangml==0) write(6,*) 'ML: Reading DB file class, klm , no_total, no_selec'
   do i=1,n_db_class
      if (rangml==0) write(*,'(i3, "  ",(a)," ",(a),i6,i6)') i, db_model(i)%class, db_model(i)%klm, db_model(i)%no_total, db_model(i)%no_selec
   end do
end if
return
end subroutine read_db_file










subroutine buildsubdata(nd_data)

use ml_in_ndm_module, only : rangml,selection_type,pref,ns_data,seed,kelem

implicit none

integer,intent(in) :: nd_data

integer :: i,j
integer,dimension(ns_data,kelem) :: sample
double precision :: y_target

namelist /input_ml/ selection_type,pref,ns_data,kelem,seed

select case(selection_type)
  case(1)
     if (rangml==0) write(6,*)"ML: Selects first",ns_data,"elements of database"
     do i=1,ns_data
        call read_poscar(pref,i,y_target)
     enddo
  case(2)
     if (rangml==0) write(6,*)"ML: Selects last",ns_data,"elements of database"
     do i=nd_data-(ns_data-1),nd_data
        call read_poscar(pref,i,y_target)
     enddo
  case (3)
     if (rangml==0) write(6,*)"ML: Random selection of",ns_data,"subset of",kelem,"elements of database"
     if (ns_data .gt. kelem) then
          if (rangml==0) write(6,*) "ML error: ns_data=", ns_data,"should be lower than kelem=",kelem,"and obviously is not the case"
          stop
     end if


     call combinations(nd_data,kelem,ns_data,seed,sample)
     do i=1,ns_data
        do j=1,kelem
           call read_poscar(pref,sample(i,j),y_target)
        enddo
     enddo
   case default
      if (rangml==0) write(*,*) 'ML error: Only 3 possible choices for selection_type 1, 2 or 3. Read the manual'
      stop

end select

return
end subroutine buildsubdata


subroutine get_number_of_files(nfiles)

use ml_in_ndm_module, only : path,lpath,rangml

implicit none

integer,intent(out) :: nfiles
character(len=120) :: command
lpath=len_trim(path)

if (rangml==0) then
 write(command,*)"ls $(echo '",path,"*.poscar' | sed 's/ *//g') | wc -l > nfiles.txt"
 call system(command)
 open(31,file='nfiles.txt',action="read")
 read(31,*)nfiles
 close(31)
! call system('rm nfiles.txt')
end if


return
end subroutine get_number_of_files


subroutine combinations(n,k,n_sample,seedin,b)

implicit none

integer(kind=4),intent(in) :: n,k,seedin,n_sample
integer,dimension(n_sample,k),intent(out) :: b

integer,dimension(k) :: a
integer(kind=4) :: i,j,c,seed
logical :: test

seed=seedin
call rks2 (n,k,seed,a)
b(1,:)=a(:)
i=2
c=0
do while (i-1.ne.n_sample)
   c=c+1
   seed=seed+c
   call rks2 (n,k,seed,a)

   test = .true.
   do j = 1,i-1
      if (all(a(:) == b(j,:))) then
         test = .false.
         exit
      endif
   end do

   if (test) then
      b(i,:)=a(:)
      i=i+1
   endif
enddo

return
end subroutine combinations


subroutine rks2 ( n, k, seed, a )

use ml_in_ndm_module, only : rangml

!*****************************************************************************80
!
!! KSUB_RANDOM2 selects a random subset of size K from a set of size N.
!
!  Discussion:
!
!    This algorithm is designated Algorithm RKS2 in the reference.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 May 2015
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Second Edition,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the set.
!
!    Input, integer ( kind = 4 ) K, the size of the subset.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, integer ( kind = 4 ) A(K), the indices of the selected elements.
!
implicit none

integer (kind = 4) :: k,c1,c2,i,k0,n,seed
integer (kind = 4),dimension(k) :: a
real (kind = 8) :: r,r8_uniform_01

if ( k < 0 .or. n < k ) then
  if (rangml==0) write ( *, '(a)' ) ''
  if (rangml==0) write ( *, '(a)' ) 'KSUB_RANDOM2 - Fatal error!'
  if (rangml==0) write ( *, '(a,i8)' ) '  N = ', n
  if (rangml==0) write ( *, '(a,i8)' ) '  K = ', k
  if (rangml==0) write ( *, '(a)' ) '  but 0 <= K <= N is required!'
  stop 1
end if

if ( k == 0 ) then
  return
end if

c1 = k
c2 = n
k0 = 0
i = 0

do i = 1,n

  r = r8_uniform_01 ( seed )

  if ( real ( c2, kind = 8 ) * r <= real ( c1, kind = 8 ) ) then

    c1 = c1 - 1
    k0 = k0 + 1
    a(k0) = i

    if ( c1 <= 0 ) then
      exit
    end if

  end if

  c2 = c2 - 1

end do
return
end subroutine rks2

function r8_uniform_01 ( seed )

use ml_in_ndm_module, only : rangml

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2^1 - 1 )
!      r8_uniform_01 = seed / ( 2^31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Pierre LEcuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer (kind = 4) :: i4_huge,k,seed
  real (kind = 8) :: r8_uniform_01

  if ( seed == 0 ) then
    if (rangml==0) write ( *, '(a)' ) ''
    if (rangml==0) write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
    if (rangml==0) write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop 1
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge ( )
  end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end function r8_uniform_01

function i4_huge ( )

!*****************************************************************************80
!
!! I4_HUGE returns a "huge" I4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) I4_HUGE, a "huge" integer.
!
  implicit none

  integer ( kind = 4 ) i4_huge

  i4_huge = 2147483647

  return
end function i4_huge
