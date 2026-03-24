module d_at_at_mod
  USE T_kind_param_m, ONLY:  double
  USE atomconfig,only : atom_config
  USE cellconfig, only:cell_config,caltabtC 
  USE boxconfig,only:box_config
  USE gen_com_m, only:uwrt,lwrt, lperiod
  use vect_dist_mod,only:vect_dist
  use montecarlo_mod,only:bublcenter
  use vect_dist_mod,only:closest_at
    USE cryst_to_cart_mod, ONLY: cryst_to_cart
  implicit none 
contains

  subroutine best_at_pos(atdml,celndm,boxndm)
    implicit none
    type(box_config)::boxndm
    class(atom_config)::atdml
    type(cell_config):: celndm

    integer::i,ic,j,k,nstp(3),idec(3),imin,icl
    real(double),parameter::step=0.1d-8
    real(double)::decxyz(3),post(3,1),posopt(3,1),distmax,distcl
    nstp(:)=boxndm%zl(:)/step
    write(uwrt,*)'NSTP',nstp
    distmax=-1000.0
    do i=1,nstp(1)
       do j=1,nstp(2)
          do k=1,nstp(3)
             idec(1)=i;idec(2)=j;idec(3)=k
             decxyz(:)=float(idec(:))/nstp(:)
             post(:,1)=decxyz(1)*boxndm%at(:,1)+decxyz(2)*boxndm%at(:,2)+decxyz(3)*boxndm%at(:,3)
             call closest_at(post(:,1),atdml,celndm,boxndm,.true.,dist=distcl)
!             write(uwrt,*)'distmax= ',distmax,i,j,k,post
             if (distmax.lt.distcl) then
                distmax=distcl
                posopt(:,1)=post(:,1)
                write(uwrt,*)'distmax= ',distmax,i,j,k,post
             end if
          end do
       end do
    end do
    call cryst_to_cart(1,posopt,boxndm%bg,-1)
    
    write(uwrt,*)'POSopt', posopt
    write(uwrt,*)'distmax', distmax*1d8
  end subroutine best_at_pos
  
  subroutine d_at_at(atdml,celndm,boxndm)

   
    implicit none
    type(box_config)::boxndm
    class(atom_config)::atdml
    type(cell_config):: celndm




    integer::i,j,imin,imax,jmin,jmax,ko1,koo,i2,i1
    real(double)::dmin,dmax,dist
    
    call caltabtC(celndm,atdml,lperiod,boxndm,lchktrav=.false.)
    ! Vecteurs de la boîte et grandeurs associées à l'instant initial

    dmin=1e10
    dmax=-100
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
             write(uwrt,*)dmin,i,j
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
    write(uwrt,*)'dmin',dmin,imin,jmin
    write(uwrt,*)'dmax',dmax,imax,jmax
  end subroutine d_at_at
end module d_at_at_mod
