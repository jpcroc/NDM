subroutine prtplz(xp,ityp)


  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  implicit none

  integer  :: ityp(imm)
  real(double)  :: xp(3,imm)


  real(double), pointer,save :: natz(:)
  integer :: i,j,k,l,itr,ii,itempo
  integer,save :: ntrl,icall
  real(double):: ltr
  character :: fnamtampon*10, fnamfilm2it*20, extension*10
  ltr=0.1d-8

  if (it==0) then
     ntrl=Int(zl(3)*1.1/ltr)
     write(6,*) 'ntrl =',ntrl
     allocate (natz(0:ntrl+1))
     natz=0.
!    if (itmax.ge.10000000000) then
     if (itmax.ge.100000000) then
        write(6,*) 'trop diterations'
        stop
     end if
     return

  else
     itempo=iteplz-mod(it,iteplz)

     if ((iteplz.le.nplz).or.(itempo.lt.nplz)) then
        call cryst_to_cart (imm, xp, bg, -1)    !cart vers cryst
        do i = 1, imd
           itr=Int(ntrl*xp(3,i))
           natz(itr)=natz(itr)+1.
        end do
        call cryst_to_cart (imm, xp, at, 1)     !cryst vers cart
     end if
     if (mod(it,iteplz)==0) then
        natz=natz/min(nplz,iteplz)
        open(unit=17, file='tampon', form='formatted', status='unknown')
        if (it<=9) write (17, '(I1)') it
        if (it<=99.and.it>9) write (17, '(I2)') it
        if (it<=999.and.it>99) write (17, '(I3)') it
        if (it<=9999.and.it>999) write (17, '(I4)') it
        if (it<=99999.and.it>9999) write (17, '(I5)') it
        if (it<=999999.and.it>99999) write (17, '(I6)') it
        if (it<=9999999.and.it>99999) write (17, '(I7)') it
        if (it<=99999999.and.it>999999) write (17, '(I8)') it
        if (it<=999999999.and.it>9999999) write (17, '(I9)') it
!       if (it<=9999999999.and.it>99999999) write (17, '(I10)') it
!       if (it>=9999999999) then
        if (it>=999999999) then
           write (6, *) 'probleme de format dans prtplz.f'
           stop
        endif
        !        stop
        rewind(17)
        !        close(17)
        !        open(unit=17, file='tampon', form='formatted', status='unknown')
        read (17, '(A10)') extension
        !        write(6,*)'extension',extension
        !        stop
        close(17)
        fnamfilm2it = 'prtplz.'//extension
        open(unit=941, file=fnamfilm2it, form='formatted')
        do itr=1,ntrl
           write(941,*)itr,natz(itr)
        end do
        close(941)
        natz=0.
     end if


  end if

  return
end subroutine prtplz
