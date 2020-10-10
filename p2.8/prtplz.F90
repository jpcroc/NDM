module prtplz_mod
  USE cryst_to_cart_mod,only: cryst_to_cart
  USE gen_com_m, ONLY:at,bg,itmax,iteplz,nplz,zl,imm,it,im
  implicit none
contains

  subroutine prtplz(xp,ityp)


    USE T_kind_param_m, ONLY:  double

    implicit none

    integer  :: ityp(imm)
    real(double)  :: xp(3,imm)


    real(double), allocatable,save :: natz(:)
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
          do i = 1, im
             itr=Int(ntrl*xp(3,i))
             natz(itr)=natz(itr)+1.
          end do
          call cryst_to_cart (imm, xp, at, 1)     !cryst vers cart
       end if
       if (mod(it,iteplz)==0) then
          natz=natz/min(nplz,iteplz)

          write(extension,'(i10.10)') it

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
end module prtplz_mod
