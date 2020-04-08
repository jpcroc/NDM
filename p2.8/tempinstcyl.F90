module tempinstcyl_mod
  USE gen_com_m, ONLY:imm,bk,ncyl,cyl
  implicit none
contains
  !c******************************************************************
  function tempinstcyl(vp,ityp)     !calcul de la T instant.
    !c******************************************************************

    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m
    USE var_pot, ONLY:cm
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    real(double),dimension (3,imm) ::  vp
    integer :: ityp(imm)
    real(double)::  tempinstcyl
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------

    real(double) ::  mv2,v2
    integer :: i
    mv2=0.0

    do i = 1,imm
       if (cyl(i).EQV..true.) then
          v2= vp(1,i)**2+ vp(2,i)**2+ vp(3,i)**2
          mv2= mv2 + cm(ityp(i))*v2
       end if
    enddo

    tempinstcyl=mv2/(3.d0*float(ncyl)*bk)

    return
  end function tempinstcyl
end module tempinstcyl_mod
