module d_at_at_mod
  USE T_kind_param_m, ONLY:  double
  USE atomconfig,only : atom_config
  USE cellconfig, only:cell_config,caltabtC 
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

    integer::i,j,imin,imax,jmin,jmax,ko1,koo,i2,i1,ncelvois
    real(double),allocatable :: xpnp(:,:)
    real(double)::dmin,dmax,dist,c1,c2,c3,c1p,c2p,c3p
    call caltabtC(celndm,atdml,lperiod,boxndm)
    allocate(xpnp(3,atdml%imm))
    ! Vecteurs de la boîte et grandeurs associées à l'instant initial
!    call cryst_to_cart (atdml%im, atdml%xp, boxndm%bg, -1)
    if (lperiod) then
       xpnp(:,:)=atdml%xp(:,:)
    else 
       call notperiod(atdml%imm,atdml%xp,xpnp,boxndm%at,boxndm%bg)

    end if
    ncelvois = min(celndm%noxyz,27)-1

    dmin=1e10
    dmax=-100
!    call celndm%print
    do i=1,atdml%im
       koo = atdml%ielat(i)                          ! Numero de la cellule

       do i1 = 0, ncelvois

          ko1 = celndm%ncel(koo,i1)
          c1p = xpnp(1,i)+sum(boxndm%at(1,:)*celndm%deltadist(:,i1,koo))
          c2p = xpnp(2,i)+sum(boxndm%at(2,:)*celndm%deltadist(:,i1,koo))
          c3p = xpnp(3,i)+sum(boxndm%at(3,:)*celndm%deltadist(:,i1,koo))

          ! pour chaque atome ds la cel. voisine
          do i2 = 1, celndm%nato(ko1)
             j = celndm%atincel(i2,ko1)
                if (atdml%num_at_glob(i).ge.atdml%num_at_glob(j)) cycle !terme déja calculé
             
!#else
!             if (num_at_glob(i).ge.num_at_glob(j)) cycle !terme déja calculé
!#endif

             if (celndm%noxyz.ne.1) then
                c1 = c1p-xpnp(1,j)
                c2 = c2p-xpnp(2,j)
                c3 = c3p-xpnp(3,j)
                cv(1,1) = c1
                cv(1,2) = c2
                cv(1,3) = c3

             else
                c1 = c1p-xpnp(1,j)
                c2 = c2p-xpnp(2,j)
                c3 = c3p-xpnp(3,j)
                cv(1,1) = c1
                cv(1,2) = c2
                cv(1,3) = c3
                call cryst_to_cart (1, cv, boxndm%bg, -1) !cart vers cryst sur cv
                WHERE ( (cv.GT.0.5d0).OR.(cv.LT.-0.5d0) )
                   cv(:,1:3) = cv(:,1:3) - Dble(Nint(cv(:,1:3)))
                END WHERE
                call cryst_to_cart (1, cv, boxndm%at, 1) !cryst vers cart sur cv
                c1=cv(1,1)
                c2=cv(1,2)
                c3=cv(1,3)

             end if
             dist =1d8*dsqrt(cv(1,1)**2+cv(1,2)**2+cv(1,3)**2)
             

       
          if (dist.lt.dmin)then
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
 end do
    write(6,*)'dmin',dmin,imin,jmin
    write(6,*)'dmax',dmax,imax,jmax
  end subroutine d_at_at
end module d_at_at_mod
