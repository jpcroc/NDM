module heat_mod
        USE notperiod_mod,only: notperiod
        USE cryst_to_cart_mod,only: cryst_to_cart
        USE gen_com_m, ONLY:bk,eheat,lperiod,rheat,theat
        USE temp_com,only:at,bg
        implicit none
        contains
! *************************************************************
subroutine heat(im,xp,vp,ityp)
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  USE var_pot, ONLY:cm
  USE jqmod

  implicit none
    integer,intent(in)::im
    real(double),allocatable,intent(in)::xp(:,:)
    real(double),allocatable::vp(:,:)
    integer,allocatable,intent(in)::ityp(:)
  
  
  real(double), dimension(:,:), allocatable :: xpnp
  integer ::i,natheat,j
  real(double):: dxp(3),r2,tsph,xtr(3)
  integer, allocatable,dimension (:):: iatheat


  write(6,*)'entree heat'
  ALLOCATE(xpnp(3,im))
  allocate (iatheat(im))

  if (lperiod) then
   xpnp(:,:)=xp(:,:)
  else
   call notperiod(im,xp,xpnp,at,bg)
  end if

  call cryst_to_cart (im, xpnp, bg, -1)    !cart vers cryst
  natheat=0
  tsph=0.

  do i=1,im
     
     dxp(1:3) = xpnp(1:3,i) - 0.5
     WHERE ( (dxp.GT.0.5d0).OR.(dxp.LT.-0.5d0) )
        dxp(1:3) = dxp(1:3) - Dble(Nint(dxp(1:3)))
     END WHERE
     ! Transformation des coordonnÃ©es rÃ©duites en cartÃ©siennes

     
     dxp = MatMul(at,dxp)
     r2 = Sum( dxp(1:3)**2 )
!     write(6,*)i,r2,rheat**2
     if (r2.le.rheat**2) then
        natheat=natheat+1
        iatheat(natheat)=i
        tsph = tsph+(vp(1,i)**2+vp(2,i)**2+vp(3,i)**2)*cm(ityp(i))/(3.0*bk)        
     end if
  end do
  tsph=tsph/float(natheat)
  if (Eheat.ne.0) Theat=2*Eheat/(3*bk*natheat)
  write(6,*)'Tsphere Theat natheat', tsph,theat,natheat

  do j=1,natheat
     i=iatheat(j)
!     xtr(:) = xp(:,i)
!     xpp(:,i)=xp(:,i)-(xp(1,:)-xpp(:,i))*sqrt(theat/tsph)
     vp(:,i) = vp(:,i)*sqrt((theat+tsph)/tsph)
  end do
  deallocate(xpnp)
  deallocate(iatheat)

  return
end subroutine heat

end module
