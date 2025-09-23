module tccontr
  USE T_kind_param_m
  use cryst_to_cart_mod,only:cryst_to_cart
  USE atomconfig,only:atom_config_d,atom_config_e
  use cellconfig,only: cell_config,caltabtC
  use boxconfig,only:box_config,periodbox
  use gen_com_m,only:lperiod,tfcou,epcou,couxyz
  USE calctemp_mod,only: calctemp
  use notperiod_mod,only:notperiod

  implicit none

contains

  subroutine contrTcou(atcf,celcf,boxcf)

    ! **************************************************************
    implicit none
    class(atom_config_d),intent(inout)::atcf
    class(cell_config),intent(inout)::celcf
    class(box_config)::boxcf
    integer :: i, im,imm,ic
    type(atom_config_e)::atcou
    type(cell_config):: celcou
    real(double),allocatable::xpnp(:,:)
    real(double)::epc(3),kinecou,tempcou
    logical:: lgs(atcf%imm)


    imm =atcf%imm; im=atcf%im
    ALLOCATE(xpnp(3,imm))
    call notperiod(im,atcf%xp,xpnp,boxcf%at,boxcf%bg,lperiod)       

    call cryst_to_cart (im, xpnp, boxcf%bg, -1) ! cart vers cryst
    do ic=1,3
       epc(ic) = epcou/boxcf%at(ic,ic)
    end do
    tempcou = 0.d0

    lgs=atcf%lgul
    
    atcf%lgul=.false.
    celcou=celcf

    loopi: do i = 1, im
       do ic=1,3
          if (couxyz(ic)==1) then 
             if (xpnp(ic,i)<epc(ic).or.xpnp(ic,i)>1.0-epc(ic)) then
                atcf%lgul(i)=.true.
                cycle loopi
             end if
          end if
       end do
    end do loopi
    call atcf%fab(atcou,lback=.true.)
    call caltabtC(celcou,atcou,lperiod,boxcf,lchktrav=.false.)
    call calctemp(tempcou,kinecou,atcou,celcou)
    
          if (atcou%lxpp)          atcou%xpp(:,:)=atcou%xp(:,:)-(atcou%xp(:,:)-atcou%xpp(:,:))*sqrt(tfcou/tempcou)
    
          atcou%vp(:,:)=atcou%vp(:,:)*sqrt(tfcou/tempcou)
    call atcou%backto(atcf)

    atcf%lgul=lgs
  end subroutine contrTcou

end module tccontr
