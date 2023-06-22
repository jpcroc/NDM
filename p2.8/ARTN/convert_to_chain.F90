module convert_to_chain_mod
contains
  subroutine convert_to_chain( init_number, length, chain )

  implicit none

  !Arguments
  integer,          intent(in)  :: init_number
  integer,          intent(in)  :: length      !can be up to 4
  character(len=4), intent(out) :: chain

  !Local variables
  character(len=10) :: digit = '0123456789'
  integer :: i, decades, divider, remainder, number, lm

  number = init_number
  if ( number == 0 ) then
     chain =''
     do i = 1, length
        chain = trim(chain)//"0"
     end do
     return
  else
     decades = log10( 1.0d0 * number) + 1
  end if

  if ( decades > length ) then   ! WARNING: chain will be is zero.
     chain =''
     do i = 1, length
        chain = trim(chain)//"0"
     end do
     return
  end if

  divider = 1
  do i = 2, decades
     divider =  divider * 10
  end do

  lm = length - decades - 1
  chain = ''
  do i = 1, decades
     remainder = number / divider  + 1
     chain = trim(chain)// digit(remainder:remainder)
     remainder = remainder -1
     number = number - remainder * divider
     divider = divider / 10
  end do

  do i= 1, length
     if ( len(trim(chain)) == length ) exit
     chain = "0"//trim(chain)
  end do

END SUBROUTINE convert_to_chain
end module convert_to_chain_mod
