module d_at_at_mod
  USE T_kind_param_m, ONLY:  double
  USE atomconfig,only : atom_config
  USE cellconfig, only:cell_config,caltabtC 
  USE boxconfig,only:box_config
  USE gen_com_m, ONLY: lperiod
  use vect_dist_mod,only:vect_dist
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

    integer::i,j,imin,imax,jmin,jmax,ko1,koo,i2,i1
    real(double)::dmin,dmax,dist,c1,c2,c3,c1p,c2p,c3p
    
    call caltabtC(celndm,atdml,lperiod,boxndm)
    ! Vecteurs de la boîte et grandeurs associées à l'instant initial

    dmin=1e10
    dmax=-100
!    call celndm%print
    do i=1,atdml%im
       koo = atdml%ielat(i)                          ! Numero de la cellule

       do i1 = 0, celndm%ncelvois(koo)

          ko1 = celndm%ncel(koo,i1)

          ! pour chaque atome ds la cel. voisine
          do i2 = 1, celndm%nato(ko1)
             j = celndm%atincel(i2,ko1)
                if (atdml%num_at_glob(i).ge.atdml%num_at_glob(j)) cycle !terme déja calculé
                
                call vect_dist(atdml,celndm,boxndm,i,j,dist=dist,lperiod=lperiod)

             dist=dist*1d8

       
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
