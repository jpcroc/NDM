module d_at_at_mod
  USE T_kind_param_m, ONLY:  double
  USE atomconfig,only : atom_config
  USE cellconfig, only:cell_config 
  USE boxconfig,only:box_config
  USE gen_com_m, ONLY: lperiod
USE cryst_to_cart_mod,only: cryst_to_cart
 USE notperiod_mod,only: notperiod
  implicit none 
contains
  ! boucle de DM pour velocity Verlet
  ! ************************************************

  subroutine d_at_at(atdml,celndm,boxndm)

   
    implicit none
    type(box_config)::boxndm
    class(atom_config)::atdml
    type(cell_config):: celndm

    REAL(double), dimension(1:3) :: cp, dxp
    REAL(double), dimension(1,1:3) :: cv

    integer::i,j,imin,imax,jmin,jmax
    real(double),allocatable :: xpnp(:,:)
    real(double)::dmin,dmax,dist
    allocate(xpnp(3,atdml%imm))
    ! Vecteurs de la boîte et grandeurs associées à l'instant initial
!    call cryst_to_cart (atdml%im, atdml%xp, boxndm%bg, -1)
    if (lperiod) then
       xpnp(:,:)=atdml%xp(:,:)
    else 
       call notperiod(atdml%imm,atdml%xp,xpnp,boxndm%at,boxndm%bg)
    end if
    dmin=1e10
    dmax=-100
    do i=1,atdml%im
       do j=i+1,atdml%im
          cv(1,1:3) = xpnp(1:3,i)-xpnp(1:3,j)
          call cryst_to_cart (1, cv, boxndm%bg, -1) !cart vers cryst sur cv
          WHERE ( (cv.GT.0.5d0).OR.(cv.LT.-0.5d0) )
             cv(:,1:3) = cv(:,1:3) - Dble(Nint(cv(:,1:3)))
          END WHERE
          call cryst_to_cart (1, cv, boxndm%at, 1) !cryst vers cart sur cv
          dist =dsqrt(cv(1,1)**2+cv(1,2)**2+cv(1,3)**2)
          if (dist.le.dmin)then
             dmin=dist
             imin=i; jmin=j
             write(6,*)dmin,i,j
          end if
          if (dist.gt.dmax)then
             dmax=dist
             imax=i;jmax=j
          end if
!          write(16,*) i,j, dist
!          write(17,*)dist
       end do
    end do
    write(6,*)'dmin',dmin,imin,jmin
    write(6,*)'dmax',dmax,imax,jmax
  end subroutine d_at_at
end module d_at_at_mod
