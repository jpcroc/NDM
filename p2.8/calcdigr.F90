module calcdigr_mod
  USE notperiod_mod,only: notperiod
  USE cryst_to_cart_mod,only: cryst_to_cart
  USE gen_com_m, ONLY:noxyz,ncel,atincel,lperiod,rang,atincel,rcrdf,at,nato,&
       &deltadist,celsize,bg

  implicit none
contains
  subroutine calcdigr(im,xp,ityp,ielat)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double
    USE var_pot, ONLY:nkmax,ntyp,nad,digr,gdertot



    implicit none

    integer,intent(in)::im
    real(double),intent(in),allocatable::xp(:,:)
    integer,allocatable,intent(in)::ityp(:),ielat(:)

    integer :: i, iti, itj, i1, i2, icell, kx, ky, kz, koo, ko1, j, &
         ic, k, m,m1,n,iti1,iti2
    real(double) :: rij, c1, c2, c3, x1, x2, x3, rmax2,rmax, incre
    real(double) :: rspace2,invincre
    real(double) :: aaa, bbb, ccc,cv(1,3)
    real(double) :: xpnp(3,im)


    rmax=rcrdf*1.0d-8
    if (rmax.gt.minval(celsize)) then
       write(6,*)'diminuer nox, noy, noz'
       stop
    end if
    !      write(6,*)'rmax ',rmax
    rmax2 = rmax**2
    incre = rmax/nkmax
    invincre = 1/incre
    if (rang==0) write(6,*) 'nkmax incre',nkmax,incre

    if (lperiod) then
       xpnp(:,:)=xp(:,:)
    else 
       call notperiod(im,xp,xpnp)
    end if

    do i = 1, im
       koo = ielat(i)
       do i1 = 0, 26
          ko1 = ncel(koo,i1)
          do i2 = 1, nato(ko1)
             j = atincel(i2,ko1)
             if(j==i) cycle 

             c1 = xpnp(1,i)-xpnp(1,j) 
             c2 = xpnp(2,i)-xpnp(2,j) 
             c3 = xpnp(3,i)-xpnp(3,j) 

             c1 = c1+sum(at(1,:)*deltadist(:,i1,koo))
             c2 = c2+sum(at(2,:)*deltadist(:,i1,koo))
             c3 = c3+sum(at(3,:)*deltadist(:,i1,koo))
             if (noxyz.ne.1) then
                if (abs(c1)>rmax) cycle
                if (abs(c2)>rmax) cycle
                if (abs(c3)>rmax) cycle 
             else
                cv(1,1) = c1
                cv(1,2) = c2
                cv(1,3) = c3
                call cryst_to_cart (1, cv, bg, -1) !cryst vers cart sur cv
                WHERE ( (cv.GT.0.5d0).OR.(cv.LT.-0.5d0) )
                   cv(:,1:3) = cv(:,1:3) - Dble(Nint(cv(:,1:3)))
                END WHERE
                call cryst_to_cart (1, cv, at, 1) !cryst vers cart sur cv
                c1=cv(1,1)
                c2=cv(1,2)
                c3=cv(1,3)
             end if
             if ((c1**2+c2**2+c3**2)> rmax2) cycle



             rij = sqrt(c1*c1+c2*c2+c3*c3)   


             k= int(rij*invincre)
             write(16,'(2I4,G15.7,2I4)') i,j,rij*1d8,ityp(i),ityp(j)
             m=k+1
             digr(ityp(i),ityp(j),m) = digr(ityp(i),ityp(j),m)+1.
             !                 write(6,*)' digr ', digr(ityp(i),ityp(j),m)
          end do
       end do
    end do
    !------------------------------------------------
    ! Calcul de la RDF total
    !------------------------------------------------

    do iti1=1,ntyp
       if(nad(iti1)==0) cycle
       do iti2=1,ntyp
          if(nad(iti2)==0) cycle
          do m=1,nkmax
             gdertot(m)=gdertot(m)+digr(iti1,iti2,m)/(nad(iti1)*nad(iti2))
          enddo

       enddo
    enddo

    !  if (rang==0) write(6,*) 'PARA-T sortie calcdigr'

    return
  end subroutine calcdigr

end module calcdigr_mod
