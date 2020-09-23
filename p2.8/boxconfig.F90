module boxconfig
  USE T_kind_param_m
  use recips_mod,only:recips,calcvol
  implicit none
  type box_config
     real(double):: at(3,3)
     real(double):: bg(3,3)
     real(double):: zl(3),zls2(3),nzl(3),volu,normat(3),normbg(3)
     integer(long)::icaltabt 
   contains

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


end module boxconfig


