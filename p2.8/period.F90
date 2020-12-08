module period_mod
  USE cryst_to_cart_mod,only: cryst_to_cart
  USE temp_com,only:at,bg ! A EFFACER
  implicit none 
contains
  ! *****************************************************************
  subroutine period (imm,xp,xpp,ax)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double
    USE gen_com_m, ONLY:DECAL_bc,ldecal_bc,low_limit,lperiod,zero

    !       version du 09 decembre 2003

    ! *****************************************************************
    ! Cette routine applique les conditions périodiques par décalage
    !  ou +/-1 (triclinique).
    !Le décalage est fait sur xp xp ET ax. Ce décalage sur ax 
    !permet de mesurer correctement le déplacements des atomes
    ! depuis leur position de départ

    implicit none
    !-----------------------------------------------
    !   G l o b a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer,intent(in)::imm
    real(double),intent(inout)::xp(3,imm)
    real(double),optional,intent(inout)::xpp(3,imm),ax(3,imm)
    !-----------------------------------------------
    !   L o c a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: i, ic,icp
    real(double)::dz,trav,ecav,ecap
    real(double):: cpp,xpici,cppzl,ctest
    !      integer,save  :: iperiod
    !  if (rang==0) write(6,*)'PARA-T entree period'
    if(.not.lperiod) return

    !      iperiod=iperiod+1

    !      write(*,*) 'PBC PBC PBC capitala tarii e ....             period',iperiod

    IF (ldecal_bc.EQV..FALSE.) THEN

       call cryst_to_cart (imm, xp,  bg,  -1) !cart vers cryst
       if(present(xpp))      call cryst_to_cart (imm, xpp, bg,  -1)
       if(present(ax))      call cryst_to_cart (imm, ax,  bg,  -1)
       do i=1,imm
          do ic=1,3
             xpici=xp(ic,i)
             if ( (xpici < 0.d0 ).OR.( xpici >= 1.d0 ) ) then
                if ( (xpici > -low_limit).and.(xpici<0.d0) ) then
                   xp(ic,i)=zero
                else
                   cpp  = Dble(Floor(xp(ic,i)))
                   if(present(xpp))                 xpp(ic,i) = xpp(ic,i) - cpp
                   if(present(ax))                 ax (ic,i) = ax(ic,i)  - cpp
                   xp (ic,i) = xpici     - cpp
                end if
             end if
          end do
       end do
       call cryst_to_cart (imm, xp , at,  1)  !cryst vers cart
       if(present(xpp))      call cryst_to_cart (imm, xpp, at,  1)
       if(present(ax))      call cryst_to_cart (imm, ax , at,  1)

    ELSE IF (ldecal_bc.EQV..TRUE.) THEN

       call cryst_to_cart (imm, xp,  bg,  -1) !cart vers cryst
       if(present(xpp))         call cryst_to_cart (imm, xpp, bg,  -1)
       do i=1,imm
          XP(3,i)  = XP(3,i)  - DECAL_bc*int(XP(1,i))

          if(present(xpp))		XPP(1,i) = XPP(1,i) - int(XP(1,i))
          XP(1,i) = XP(1,i) - int(XP(1,i))

          if(present(xpp))		XPP(2,i) = XPP(2,i) - int(XP(2,i))
          XP(2,i) = XP(2,i) - int(XP(2,i))

          if(present(xpp))		XPP(3,i) = XPP(3,i) - int(XP(3,i))
          XP(3,i) = XP(3,i) - int(XP(3,i))
       end do
       call cryst_to_cart (imm, xp , at,  1)  !cryst vers cart
       if(present(xpp))       call cryst_to_cart (imm, xpp, at,  1)
    END IF


    !  if (rang==0) write(6,*)'PARA-T sortie period'


  end subroutine period
end module
