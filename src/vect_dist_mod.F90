module vect_dist_mod
  USE arret_ndm_mod,only:arret_ndm
  USE T_kind_param_m
  USE notperiod_mod,only: notperiod
  USE cryst_to_cart_mod,only: cryst_to_cart
  USE atomconfig,only:atom_config
  use cellconfig,only: cell_config
  use boxconfig,only:box_config
  implicit none
contains
  subroutine vect_dist(atcf,celcf,boxcf,i,j,VJI,indcv,lperiod,rum,dist,linter)
    class(atom_config),intent(in)::atcf
    class(cell_config),intent(in)::celcf
    class(box_config)::boxcf
    integer,intent(in)::i,j !i j indices des atomes,indcv= rang dans le voisinage de la cellule de i !! (de 1 à 27)
    integer,intent(in),optional ::indcv !indcv= rang dans le voisinage de la cellule de i !! (de 1 à 27)
    logical::lperiod
    real(double),intent(out),optional::VJI(3)
    real(double)::XJI(3)
    logical,intent(out),optional::linter
    real(double),optional,intent(in)::rum
    real(double),optional,intent(out)::dist

    real(double)::xpnp(3,2),xp(3,2),cv(1,3),dx2(3)
    integer::ic,ic2
    real(double)::distance

    if (((present(rum)).and.(.not.(present(linter)))).or.((present(linter)).and.(.not.(present(rum))))) then
       write(6,*)'incohérence dans appel a vect_dist'
       call arret_ndm
    end if
    xp(:,1)=atcf%xp(:,i)
    xp(:,2)=atcf%xp(:,j)
    call notperiod(2,xp,xpnp,boxcf%at,boxcf%bg,lperiod)
    XJI(:)= xpnp(:,1)-xpnp(:,2)
    if (present(indcv))then
       do ic2=1,3
          if (celcf%nox(ic2).ge.3) then 
             do ic=1,3
                XJI(ic)=XJI(ic)+boxcf%at(ic,ic2)*celcf%deltadist(ic2,indcv,atcf%ielat(i))
             end do
          else
             cv(1,:) = XJI(:)
             call cryst_to_cart (1, cv, boxcf%bg, -1) !cart vers cryst cryst vers cart sur cv
             if (boxcf%ipbc(ic2)==1) then
                if ( (cv(1,ic2).GT.0.5d0).OR.(cv(1,ic2).LT.-0.5d0) )then
                   cv(1,ic2) = cv(1,ic2) - Dble(Nint(cv(1,ic2)))
                end if
             end if
             call cryst_to_cart (1, cv, boxcf%at, 1) !cryst vers cart sur cv
             XJI(:)=cv(1,:)
          end if
       end do
    else
       cv(1,:) = XJI(:)
       call cryst_to_cart (1, cv, boxcf%bg, -1) !cart vers cryst cryst vers cart sur cv
       do ic=1,3
          if (boxcf%ipbc(ic)==1) then
             if ( (cv(1,ic).GT.0.5d0).OR.(cv(1,ic).LT.-0.5d0) )then
                cv(1,ic) = cv(1,ic) - Dble(Nint(cv(1,ic)))
             end if
          end if
       end do
       call cryst_to_cart (1, cv, boxcf%at, 1) !cryst vers cart sur cv
       XJI(:)=cv(1,:)
    end if
    if (present(VJI))VJI=XJI
    if (present(rum)) then
       if (any(abs(XJI)>rum)) then
          linter=.false.
          return
       end if
       dx2(:)=XJI(:)*XJI(:)
       distance=sqrt(sum(dx2(:)))
       if (distance.gt.rum) then
          linter=.false.
          if (present(dist)) dist=distance
       else
          linter=.true.
          if (present(dist)) dist=distance
       end if
    else
       if (present(dist)) then
          dx2(:)=XJI(:)*XJI(:)
          dist=sqrt(sum(dx2(:)))
       end if
    end if
    return   
  end subroutine vect_dist

  subroutine closest_at(xPtest,atcf,celcf,boxcf,lperiod,iclose,rumin,dist,lclose,itypt)
    class(atom_config),intent(in)::atcf
    class(cell_config),intent(in)::celcf
    class(box_config)::boxcf
    logical::lperiod
    real(double),intent(in)::xptest(3)
    logical,intent(out),optional::lclose ! true si distmin < rumin false sinon
    real(double),optional,intent(in)::rumin !
    real(double),optional,intent(out)::dist ! distance minimel effective
    integer,intent(out),optional ::iclose !i= indice du plus proche
    integer,intent(in),optional ::itypt !type des atomes à tester

    real(double)::xpnp(3,2),xp(3,2),dx2(3),XJI(3)
    integer::ic,i
    real(double)::distance,distance0

    if (((present(rumin)).and.(.not.(present(lclose)))).or.((present(lclose)).and.(.not.(present(rumin))))) then
       write(6,*)'incohérence dans appel a closest_at'
       call arret_ndm
    end if

    if (present(lclose))lclose=.false.
    distance0=1d10
    do i=1,atcf%im
       
       if (present(itypt).and.(atcf%ityp(i).ne.itypt)) cycle
       xp(:,1)=xptest(:)
       xp(:,2)=atcf%xp(:,i)
       call notperiod(2,xp,xpnp,boxcf%at,boxcf%bg,lperiod)
       XJI(:)= xpnp(:,1)-xpnp(:,2)
       !    if ((celcf%noxyz.ne.1).and.(i1.ge.1).and.(i1.le.27)) then
!       cv(1,:) = XJI(:)
       call cryst_to_cart (1, XJI, boxcf%bg, -1) !cart vers cryst 
       do ic=1,3
          if (boxcf%ipbc(ic)==1) then
             if ( (XJI(ic).GT.0.5d0).OR.(XJI(ic).LT.-0.5d0) )then
                XJI(ic) = XJI(ic) - Dble(Nint(XJI(ic)))
             end if
          end if
       end do
       call cryst_to_cart (1, XJI, boxcf%at, 1) !cryst vers cart sur cv


       dx2(:)=XJI(:)*XJI(:)
       distance=sqrt(sum(dx2(:)))

       if (present(rumin)) then   
          if (distance.le.rumin) then
             if (present(lclose))             lclose=.true.
             if (present(iclose)) iclose=i
             if (present(dist)) dist=distance
             return
          end if
       end if
       if (distance.lt.distance0) then
          distance0=distance
          if (present(iclose)) iclose=i
       end if
    end do
    if (present(dist)) dist=distance0
             
    
    return
    
  end subroutine closest_at

  subroutine distat(xi,x0,box,dist,x0red,rum,linter)
    class(box_config),intent(in)::box
    real(double), dimension(3),intent(in)::xi
    real(double),optional,intent(in)::x0red(3),x0(3)
    real(double),optional,intent(in)::rum !
    real(double),optional,intent(out)::dist ! distance minimel effective
    logical,optional::linter
    real(double),dimension(3)::dx
    real(double),dimension (3,2)::xat
    real(double)::distc
    integer::ns=2,ic

    if (((present(rum)).and.(.not.(present(linter)))).or.((present(linter)).and.(.not.(present(rum))))) then
       write(6,*)'incohérence dans appel a vect_dist'
       call arret_ndm
    end if

    if (present(x0).and.present(x0red)) then
       write(6,*)'x0 and x0red in distat STOP'
       call arret_ndm
    end if
    if ((.not.(present(x0))).and.(.not.(present(x0red)))) then
       write(6,*)'no x0 nor  x0red in distat STOP'
       call arret_ndm
    end if
    xat(:,1)=xi(:)
    if (present(x0)) then
       xat(:,2)=x0(:)
       call cryst_to_cart (ns,xat,box%bg,-1)
       dx(1:3)=xat(1:3,1)-xat(1:3,2)
    else
       xat(:,2)=0
       call cryst_to_cart (ns,xat,box%bg,-1)
       dx(1:3)=xat(1:3,1)-x0red(1:3)
    end if
    do ic=1,3
       if (box%ipbc(ic)==1) then
          if ( (dx(ic).GT.0.5d0).OR.(dx(ic).LT.-0.5d0) )then
             dx(ic) = dx(ic) - Dble(Nint(dx(ic)))
          end if
       end if
    end do

!!$    WHERE ( (dx.GT.0.5d0).OR.(dx.LT.-0.5d0) )
!!$       dx(1:3) = dx(1:3) - Dble(Nint(dx(1:3)))
!!$    END WHERE
    dx = MatMul(box%at,dx)
    distc = sqrt(Sum( dx(1:3)**2 ))
    if (present(dist))dist=distc
    if (present(rum)) then
       if (distc.gt.rum) then
          linter=.false.
       else
          linter=.true.
       end if
    end if
    

       
    return
  end subroutine distat


end module vect_dist_mod
