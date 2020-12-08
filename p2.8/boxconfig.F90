module boxconfig
  USE T_kind_param_m
  use recips_mod,only:recips,calcvol
  use atomconfig,only:atom_config,atom_config_d,atom_config_e
  implicit none
  type box_config
     real(double):: at(3,3),h0(3,3)
     real(double):: bg(3,3)
     real(double):: zl(3),zls2(3),nzl(3),volu,normat(3),normbg(3)
     integer(long)::icaltabt 

   contains
     procedure, pass::print=>boxprint
  end type box_config
contains
  subroutine initbox(boxnew,at,zl)
    type(box_config),intent(out)::boxnew
    real(double),intent(in),optional::at(3,3)
    real(double),optional,intent(in)::zl(3)
    integer::i,ic
    if((present(zl).eqv..false.).and.(present(at).eqv..false.)) then
       write(6,*)'box init at ET zl indéfinis : STOP'
       stop
    end if
    if(present(zl).and.(present(at))) then
       write(6,*)'box init at ET zl définis : STOP'
       stop
    end if
    if (present(at)) then
       boxnew%at=at
       !                    write(6,*)at
    else
       boxnew%at=0
       do ic = 1, 3
          boxnew%at(ic,ic)=boxnew%zl(ic)
       end do
    end if

    call recips (at(1,1), at(1,2), at(1,3), boxnew%bg(1,1), boxnew%bg(1,2), boxnew%bg(1,3))
    do ic = 1, 3
       boxnew%normat(ic) = 0
       boxnew%normat(ic) = boxnew%normat(ic)+sum(at(:,ic)**2)
       boxnew%normat(ic) = sqrt(boxnew%normat(ic))
       boxnew%zl(ic) = boxnew%normat(ic)
       boxnew%normbg(ic)=sqrt(sum(boxnew%bg(:,ic)**2))
        boxnew%nzl(ic)=1.0/ boxnew%normbg(ic)
    end do
    boxnew%zls2 = boxnew%zl/2.0
    boxnew%volu=calcvol(at(1:3,1),at(1:3,2),at(1:3,3))
    boxnew%h0=at
    return
  end subroutine initbox


  subroutine ndm2boxconfig(at,bg,zl,zls2,nzl,volu,normat,boxnew)
    real(double):: at(3,3)
    real(double):: bg(3,3)
    real(double):: zl(3),zls2(3),nzl(3),volu,normat(3)
    integer(long)::icaltabt 
    type(box_config)::boxnew
    boxnew%at=at
    boxnew%bg=bg
    boxnew%zl=zl
    boxnew%zls2=zls2
    boxnew%nzl=nzl
    boxnew%volu=volu
    boxnew%normat=normat
    boxnew%icaltabt=0
    return
  end subroutine ndm2boxconfig
  subroutine boxconfig2ndm(at,bg,zl,zls2,nzl,volu,normat,boxnew)
    real(double):: at(3,3)
    real(double):: bg(3,3)
    real(double):: zl(3),zls2(3),nzl(3),volu,normat(3)
    integer(long)::icaltabt 
    type(box_config)::boxnew
    at=boxnew%at
    bg=boxnew%bg
    zl=boxnew%zl
    zls2=boxnew%zls2
    nzl=boxnew%nzl
    volu=boxnew%volu
    normat=boxnew%normat
    return
  end subroutine boxconfig2ndm

  subroutine boxprint(boxprt)
    class(box_config)::boxprt
     write(6,*)'boxprt at',boxprt%at(:,:)
     write(6,*)'boxprt bg',boxprt%bg(:,:)
     write(6,*)'boxprt volu',boxprt%volu
     write(6,*)'boxprt icaltabt',boxprt%icaltabt
   end subroutine boxprint

   subroutine periodbox(box,atcf)
     !-----------------------------------------------
     !   M o d u l e s
     !-----------------------------------------------
     USE cryst_to_cart_mod,only: cryst_to_cart
     USE T_kind_param_m, ONLY:  double
     USE gen_com_m, ONLY:lperiod,zero

    USE gen_com_m, ONLY:low_limit,zero

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
     type(box_config),intent(in)::box
     class(atom_config)::atcf

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


     call cryst_to_cart (atcf%imm, atcf%xp,  box%bg,  -1) !cart vers cryst
     select type (atcf)
     type is (atom_config_d)
        call cryst_to_cart (atcf%imm, atcf%xpp, box%bg,  -1)
     type is (atom_config_e)
        call cryst_to_cart (atcf%imm, atcf%xpp, box%bg,  -1)
        if(atcf%lax)      call cryst_to_cart (atcf%imm, atcf%ax,  box%bg,  -1)
     end select
     do i=1,atcf%imm
        do ic=1,3
           xpici=atcf%xp(ic,i)
           if ( (xpici < 0.d0 ).OR.( xpici >= 1.d0 ) ) then
              if ( (xpici > -low_limit).and.(xpici<0.d0) ) then
                 atcf%xp(ic,i)=zero
              else
                 cpp  = Dble(Floor(atcf%xp(ic,i)))
                 select type (atcf)
                 type is (atom_config_d)
                    atcf%xpp(ic,i) = atcf%xpp(ic,i) - cpp
                 type is (atom_config_e)
                    atcf%xpp(ic,i) = atcf%xpp(ic,i) - cpp
                    if(atcf%lax)   atcf%ax (ic,i) = atcf%ax(ic,i)  - cpp
                 end select
                 atcf%xp (ic,i) = xpici     - cpp
              end if
           end if
        end do
     end do
     call cryst_to_cart (atcf%imm, atcf%xp , box%at,  1)  !cryst vers cart
     select type (atcf)
     type is (atom_config_d)
        call cryst_to_cart (atcf%imm, atcf%xpp , box%at,  1)  !cryst vers cart
     type is (atom_config_e)
        call cryst_to_cart (atcf%imm, atcf%xpp , box%at,  1)  !cryst vers cart
        if(atcf%lax) call cryst_to_cart (atcf%imm, atcf%ax , box%at,  1)  !cryst vers cart
     end select
   end subroutine periodbox


 end module boxconfig


