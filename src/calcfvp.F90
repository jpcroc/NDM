module calcfvp_mod
  USE T_kind_param_m,only:double
  use var_pot,only : cm
  real(double),allocatable,dimension(:)::aspl,bspl,vmin,vmax,S0spl,xs

contains
  subroutine calcfvp(fvp,iti,vp,mov)
    real(double),intent(out)::fvp
    real(double),intent(in)::vp(3)
    integer,intent(in)::iti,mov
    real(double)::m,vn,xpar
    if (mov==2) then
       fvp=1
    else
       
       m=cm(iti)
       vn=norm2(vp(:))
       xpar=vn*m-vmin(iti)*m
       !                write(uwrt,*)'xpar',xpar,xs(iti)
       if (xpar.le.0)  then
          fvp=0
       else if ((xpar.lt.xs(iti)).and.(xpar.gt.0)) then
          fvp=pol(xpar,aspl(iti),bspl(iti))/vn
       else
          fvp=1
          !                   write(uwrt,*)'POOOOO'
          !                   call arret_ndm
       end if
    end if
    return
  end subroutine calcfvp
  function ppol(x,a,b)
    real*8::x,ppol,a,b
    ppol=a*x**4/4+b*x**3/3

  end function ppol

  function pol(x,a,b)
    real*8::x,pol,a,b
    pol=a*x**3+b*x**2

  end function pol

  function dpol(x,a,b)
    real*8::x,dpol,a,b
    dpol=3*a*x**2+2*b*x
  end function dpol

  function uf(x,a,b,s0,xs,pmin,m)
    real*8::uf,x,a,b,s0,xs,pmin,m
    integer::mov
    if (x.le.0)  then
       uf =S0
       mov=0
    else if ((x.lt.xs).and.(x.gt.0)) then
       uf=ppol(x,a,b)+S0
       mov=1
    else
       uf=0.5*(x+pmin)**2/m
       mov=2
    end if
  end function uf

end module calcfvp_mod
