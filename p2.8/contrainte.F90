module contrainte

  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  implicit none

  real(double),pointer,save :: dg0(:,:)
  real(double),save::dg2,contrval
  integer ::ic
contains


  ! **************************************************************
  subroutine initcontr(xp, xpp, vp, ax,ityp )
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double
    use gen_com_m
    implicit none
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    real(double)  :: xp(3,imm)
    real(double)  :: xpp(3,imm)
    real(double)  :: vp(3,imm)
    real(double)  :: ax(3,imm)
    integer:: ityp(im)
    integer :: i
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------

    real(double)::masstot,sigg0(3)

    integer::is,ist(1)
    real(double)::z1
    !-----------------------------------------------

    write(6,*)'entree initcontr'
    allocate (dg0(3,imm))
    open (98,file='contrainte')
    do i=1,im
       read(98,*)dg0(1:3,i)
    end do

    !  call system_clock (is) 
    !  
    ! ist(1)=12
    ! call    random_seed (put=ist)

    !  do i=1,im
    !     do ic =1,3
    !        call random_number(z1)
    !        dg0(ic,i)=z1*3.0d-8
    !     end do
    !  end do

    !dg0(:,:)=0.
    !dg0(1,1)=1.

    do i=1,im
       do ic =1,3
          sigg0(ic)=sigg0(ic)+dg0(ic,i)*cm(ityp(i))
       end do
       masstot=masstot+cm(ityp(i))
    end do
    sigg0(:)=sigg0(:)/masstot

    do i=1,im
       do ic=1,3
          dg0(ic,i)=dg0(ic,i)- sigg0(ic)
          dg2=dg2+dg0(ic,i)**2
       end do
    end do
    contrval=0.
    do i=1,im
       do ic=1,3
          contrval=contrval+xp(ic,i)*dg0(ic,i)
       end do
    end do
    write(6,*)'contrval=',contrval

    vp=0.
    xpp=xp

    !  close(98)

    return
  end subroutine initcontr

  ! **************************************************************
  subroutine contr(xp, vp, fp,ityp )
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double
    use gen_com_m
    implicit none
    !-----------------------------------------------
    !   G l o b a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    real(double)  :: xp(3,imm)
    real(double)  :: vp(3,imm)
    real(double)  :: fp(3,imm)
    integer :: ityp(im)
    !-----------------------------------------------
    !   L o c a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    real(double) :: lamda
    integer :: i
    !-----------------------------------------------

    contrval=0.
    do i=1,im
       do ic=1,3
          contrval=contrval+xp(ic,i)*dg0(ic,i)
       end do
    end do
    write(6,*)'contrval=',contrval
    !return
    lamda=0.
    do i=1,im
       do ic=1,3
          lamda=lamda-(1./dg2)*(fp(ic,i)*dg0(ic,i))
       end do
    end do

    do i=1,im
       do ic=1,3
          fp(ic,i)=fp(ic,i)+lamda*dg0(ic,i)
          !        write(6,*)'fp',fp(ic,i),lamda*dg0(ic,i)
       end do
    end do

    !write(6,*)'fp_p', fp(1,1),fp(1,2)

    return
  end subroutine contr
end module contrainte
