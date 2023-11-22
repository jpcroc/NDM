module tccontr
  USE T_kind_param_m
  use cryst_to_cart_mod,only:cryst_to_cart
  USE atomconfig,only:atom_config_d
  use cellconfig,only: cell_config,caltabtC
  use boxconfig,only:box_config,periodbox
  use gen_com_m,only:lperiod,tfcou
  USE calctemp_mod,only: calctemp
  use notperiod_mod,only:notperiod

  implicit none

contains

  subroutine contrTcou(atcf,celcf,boxcf,epcou)

    ! **************************************************************
    implicit none
    class(atom_config_d),intent(inout)::atcf
    class(cell_config),intent(inout)::celcf
    class(box_config)::boxcf
    real(double)::epcou
    integer :: i, ic,j,k,l,m,n,im,imm
    type(atom_config_d)::atcou
    type(cell_config):: celcou
    real(double),allocatable::xpnp(:,:)
    real(double)::epc1,epc2,epc3,kinecou,tempcou
    imm =atcf%imm; im=atcf%im
    ALLOCATE(xpnp(3,imm))
    call notperiod(im,atcf%xp,xpnp,boxcf%at,boxcf%bg,lperiod)       

    call cryst_to_cart (im, xpnp, boxcf%bg, -1) ! cart vers cryst
    epc1 = epcou/boxcf%at(1,1)
    epc2 = epcou/boxcf%at(2,2)
    epc3 = epcou/boxcf%at(3,3)
    tempcou = 0.d0


    atcf%lgul=.false.
    celcou=celcf

    do i = 1, im
       if (.not.(xpnp(1,i)<epc1.or.xpnp(1,i)>1.0-epc1.or.xpnp(2,i)<&
            epc2.or.xpnp(2,i)>1.0-epc2.or.xpnp(3,i)<epc3.or.xpnp(3,i)>&
            1.0-epc3)) cycle
       atcf%lgul(i)=.true.
    end do

    call atcf%fab(atcou,lback=.false.)
    call caltabtC(celcou,atcou,lperiod,boxcf)
    call calctemp(tempcou,kinecou,atcou,celcou)
    do i=1,im
       if (atcf%lgul(i)) then
          atcf%xpp(:,i)=atcf%xp(:,i)-(atcf%xp(:,i)-atcf%xpp(:,i))*sqrt(tfcou/tempcou)
          atcf%vp(:,i)=atcf%vp(:,i)*sqrt(tfcou/tempcou)
       end if
    end do

  end subroutine contrTcou

end module tccontr
