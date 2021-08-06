module vect_dist_mod
  USE T_kind_param_m
  USE notperiod_mod,only: notperiod
  USE cryst_to_cart_mod,only: cryst_to_cart
  USE atomconfig,only:atom_config
  use cellconfig,only: cell_config
    use boxconfig,only:box_config
contains
  subroutine vect_dist(atcf,celcf,boxcf,i,j,XJI,i1,lperiod,rum,dist,linter)
    class(atom_config),intent(in)::atcf
    class(cell_config),intent(in)::celcf
    type(box_config)::boxcf
    integer,intent(in)::i,j,i1 !i j indices des atomes,i1= rang dans le voisinage de la cellule de i !! (de 1 à 27)
    logical::lperiod
    real(double),intent(out)::XJI(3)

    logical,intent(out),optional::linter
    real(double),optional,intent(in)::rum
    real(double),optional,intent(out)::dist

    real(double)::xpnp(3,2),xp(3,2),cv(1,3),dx2(3)
    integer::ic
    real(double)::distance

    if (((present(rum)).and.(.not.(present(linter)))).or.((present(linter)).and.(.not.(present(linter))))) then
       write(6,*)'incohérence dans appel a vect_dist'
       stop
    end if
    xp(:,1)=atcf%xp(:,i)
    xp(:,2)=atcf%xp(:,j)
    if (lperiod) then
       xpnp(:,1)=atcf%xp(:,i)
       xpnp(:,2)=atcf%xp(:,j)
    else 
       call notperiod(2,xp,xpnp,boxcf%at,boxcf%bg)
    end if
    XJI(:)= xpnp(:,1)-xpnp(:,2)
    !    if ((celcf%noxyz.ne.1).and.(i1.ge.1).and.(i1.le.27)) then
    do ic=1,3
       XJI(ic)=XJI(ic)+sum(boxcf%at(ic,:)*celcf%deltadist(:,i1,atcf%ielat(i)))
       !       c1 = c1+sum(boxcf%at(1,:)*celcf%deltadist(:,i1,koo))
       !       c2 = c2+sum(boxcf%at(2,:)*celcf%deltadist(:,i1,koo))
       !       c3 = c3+sum(boxcf%at(3,:)*celcf%deltadist(:,i1,koo))
    end do
    !    else
    cv(1,:) = XJI(:)
    call cryst_to_cart (1, cv, boxcf%bg, -1) !cryst vers cart sur cv
    WHERE ( (cv.GT.0.5d0).OR.(cv.LT.-0.5d0) )
       cv(:,1:3) = cv(:,1:3) - Dble(Nint(cv(:,1:3)))
    END WHERE
    call cryst_to_cart (1, cv, boxcf%at, 1) !cryst vers cart sur cv
    XJI(:)=cv(1,:)
    if (present(rum)) then
       if (any(abs(XJI)>rum)) then
          linter=.false.
          return
       end if
       dx2(:)=XJI(:)*XJI(:)
       distance=sqrt(sum(dx2(:)))
       if (distance.gt.rum) then
          linter=.false.
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
end module vect_dist_mod
