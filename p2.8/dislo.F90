! ************************************************
!           Sous-programme at_bord
! alloue et remplit le tableau latdebord
! ************************************************
module dislo_mod
  USE cryst_to_cart_mod,only: cryst_to_cart
  USE temp_com, ONLY:imm,im,at,bg,zl
  USE gen_com_m, ONLY:fdislo,epcoudis,latdebord
  implicit none
contains

  subroutine at_bord(xp)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double

    implicit none
    !-----------------------------------------------
    !   G l o b a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    real(double)  :: xp(3,imm)

    !-----------------------------------------------
    !   L o c a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: i,is

    allocate (latdebord(imm))
    latdebord(:)=0
    epcoudis=epcoudis/zl(3)
    call cryst_to_cart (imm, xp, bg, -1)
    do i=1,im
       if ((xp(3,i)< epcoudis)) latdebord=-1
       if ((xp(3,i)>1.- epcoudis)) latdebord=1
    end do
    call cryst_to_cart (imm, xp, at, 1)
    epcoudis=epcoudis*zl(3)

    return  
  end subroutine at_bord


  !           Sous-programme calfodislo
  ! alloue et remplit le tableau latdebord
  ! ************************************************

  subroutine calfodislo(fp)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double
    USE gen_com_m, ONLY:
    implicit none
    !-----------------------------------------------
    !   G l o b a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    real(double)  :: fp(3,imm)

    !-----------------------------------------------
    !   L o c a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer ::i

    do i=1,im
       if(latdebord(i).ne.0) then
          fp(3,i)=fp(3,i)+dsign(fdislo,dfloat(latdebord(i)))
       end if
    end do

    return  
  end subroutine calfodislo
end module dislo_mod
