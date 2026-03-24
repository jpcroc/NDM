!  **********************************************************
!      Calcul du potentiel de Ziegler-Biersack-Littmark
!           
!  **********************************************************
module zieg3_mod
  USE arret_ndm_mod,only:arret_ndm
  USE Mat_utils_mod,only: reslin33
  use var_pot,only:auxe,lprtpot
  USE T_kind_param_m, ONLY:  double
  USE gen_com_m, only:uwrt,lwrt,erg2eV
  implicit none 
contains

  subroutine zieg3(pot, csive,ngrid, ntyp,npair, catom, roff1, roff2,lu_roff_pair,ipotentiel,typ_pot_pair,ipo)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double

    implicit none
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer, dimension(:,:), allocatable  :: ipo                      ! indice des paires d'atomes
    integer, allocatable:: typ_pot_pair(:) ! donne le type d'interaction de la paire
    integer , intent(in) :: ngrid,ipotentiel
    integer  :: ntyp
    integer  :: npair
    real(double) , intent(in) :: csive
    real(double) , intent(inout) :: pot(4,npair,0:ngrid+1)
    real(double)  :: catom(ntyp)
    real(double)  :: roff1(npair)
    real(double)  :: roff2(npair)
    logical :: lu_roff_pair(npair)


    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: j1, l, i1, i2,k,iroff1,iroff2,lw
    real(double) :: aux1, aux2, aux3, rr, r3, a0, b1, b2, b3, b4, som, r4, r5&
         , rbohr, r2, c1, c2, c3, c4, som1,aa,bb,bet,cc,fp0,dd,f0,gx,rr1,rr2,dzr1,&
         &alph,gpx,pr2,prp2,rrp1,rrp2,uu,vv,x,zr1,zrp1,r,gam,delt,zie,zieg1,&
         &ABC(3),BM(3),AM(3,3),ee,ff,g2px,rrm1,rrm2,f2p0,zrm1

    !-----------------------------------------------


    data rbohr/ 0.529D-8/
!          write(uwrt,*)'***********************************ZIEG3 l '
    c1 = 0.1818
    c2 = 0.5099
    c3 = 0.2802
    c4 = 0.02817
    if (lprtpot) then 
       do i1 = 1, ntyp
          do i2 = i1, ntyp
             l = ipo(i1,i2)
             do j1=1,ngrid
                r = j1*csive
                call calc_ziegp(zieg1,catom,ntyp,i1,i2,r)
                lw=370+l
                write(lw,'(2I6,5D15.6)')l,j1,r*1d8,zieg1,zieg1*erg2eV
             end do
          end do
       end do
    end if
    l = 0
    do i1 = 1, ntyp
       do i2 = i1, ntyp
          l = ipo(i1,i2)
          if(typ_pot_pair(l).ne.ipotentiel)cycle
          if(lu_roff_pair(l).EQV..false.) cycle
          
          iroff1=int(roff1(l)/csive)
          iroff2=int(roff2(l)/csive)
!          write(uwrt,*)'IROFF',roff1(l),iroff1,roff2(l),iroff2
          rr1=iroff1*csive; rrp1=(iroff1+1)*csive ; rrm1=(iroff1-1)*csive
          rr2=iroff2*csive; rrp2=(iroff2+1)*csive; rrm2=(iroff2-1)*csive
          call calc_ziegp(zr1,catom,ntyp,i1,i2,rr1) ! valeur de Zieg à roff1
          call calc_ziegp(zrp1,catom,ntyp,i1,i2,rrp1) ! valeur de ZiegPRIME à roff1
          call calc_ziegp(zrm1,catom,ntyp,i1,i2,rrm1) ! valeur de ZiegPRIME à roff1
!          dzr1=(zrp1-zr1)*csive
          f0=zr1
          fp0=(zrp1-zr1)/csive
          f2p0=(zrp1+zrm1-2*zr1)/(csive**2)

          gx=pot(1,l,iroff2)
          gpx=(pot(1,l,iroff2)-pot(1,l,iroff2-1))/csive
          g2px=(-2*pot(1,l,iroff2)+pot(1,l,iroff2+1)+pot(1,l,iroff2-1))/(csive**2)
          
!          write(uwrt,*)'POT Z', zr1, pr2,prp2
          x=rr2-rr1
          dd=f2p0/2
          ee=fp0
          ff=f0
!          write(uwrt,*)'dd ee ff',dd,ee,ff
          AM(1,1)=x**5
          AM(1,2)=x**4
          AM(1,3)=x**3
          AM(2,1)=5*x**4
          AM(2,2)=4*x**3
          AM(2,3)=3*x**2
          AM(3,1)=20*x**3
          AM(3,2)=12*x**2
          AM(3,3)=6*x

          BM(1)=gx-dd*x**2-ee*x-ff
          BM(2)=gpx-2*dd*x-ee
          BM(3)=g2px-2*dd

          call reslin33(ABC,BM,AM)
!          write(uwrt,*)'ABC',ABC
!          stop
          aa=ABC(1);bb=ABC(2);cc=ABC(3)
          
!!$          alph=(gx-dd-x*cc)/(x**3)
!!$          bet=1/x
!!$          gam=2/(3*x)
!!$          delt=(gpx-cc)/(3*x**2)
!!$          bb=(alph-delt)/(bet-gam)
!!$          aa=alph-bet*bb
!          write(uwrt,*)'***********************************ZIEG3 l '
!          write(uwrt,*)'ZIEG3 l ',l,aa,bb,cc,dd

          do j1 = 1, ngrid
             !               k = k+1
             r = j1*csive
             if (r<=roff1(l))then
                call calc_ziegp(pot(1,l,j1),catom,ntyp,i1,i2,r)
                
             else                              !Potentiel polynomiale de raccord entre ROFF1 et ROFF2
                if (r<roff2(l)) then
                   rr=r-rr1
                   r2=rr*rr
                   r3 = r2*rr
                   r4=r2*r2
                   r5=r2*r3
                   pot(1,l,j1) =aa*r5+bb*r4+cc*r3+dd*r2+ee*rr+ff
                   pot(2,l,j1) = 5*aa*r4+4*bb*r3+3*cc*r2+2*dd*rr+ee
                endif
             endif
          end do
          
       end do
    end do
  end subroutine zieg3



  subroutine calc_ziegp(zie,catom,ntyp,i1,i2,r)
    real(double),intent(out)::zie
    integer,intent(in)::i1,i2
    integer,intent(in)::ntyp
    real(double),intent(in)::catom(ntyp)
    real(double),intent(in)::r
    real(double) :: aux1, aux2, aux3, r3, a0, b1, b2, b3, b4, som, r4, r5&
         , rbohr, r2, c1, c2, c3, c4, som1
    data rbohr/ 0.529D-8/    
    c1 = 0.1818
    c2 = 0.5099
    c3 = 0.2802
    c4 = 0.02817
    r2 = r*r
    r3 = r2*r
    a0 = 0.88534*rbohr/(catom(i1)**0.23+catom(i2)**0.23)
    b1 = 3.2/a0
    b2 = 0.9423/a0
    b3 = 0.4029/a0
    b4 = 0.2016/a0
    aux2 = auxe*catom(i1)*catom(i2)
    som = c1*exp((-b1*r))+c2*exp((-b2*r))+c3*exp((-b3*r))+c4*&
         exp((-b4*r))
    som1 = -c1*b1*exp((-b1*r))-c2*b2*exp((-b2*r))-c3*b3*exp((-b3*r))+c4*&
         b4*((-b4*r))
    zie = aux2/r*som
!    zie(2) = -aux2/r2*som + aux2/r*som1
    return
  end subroutine calc_ziegp

end module zieg3_mod
