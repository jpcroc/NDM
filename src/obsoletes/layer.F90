module layer_mod
!  USE cryst_to_cart_mod,only: cryst_to_cart
!  USE gen_com_m, ONLY:at,bg,im,imd,imfree,imgi,imgs,imm,it,rang,rulayer,nzl
  implicit none
contains

  subroutine layer
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double

    USE tab_imm_m,only:

    write(6,*)'LAYER DETRUIT'
    stop
    

  end subroutine layer
end module layer_mod
