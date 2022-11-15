!  **********************************************************
!      Calcul du potentiel de Ziegler-Biersack-Littmark
!           
!  **********************************************************


subroutine zieg2(pot, pot_d, csive,ngrid, catom, roff1, roff2,potzt,potzt2)
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------


  implicit none
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  integer , intent(in) :: ngrid
  real*8 , intent(in) :: csive
  real*8  :: auxe= 23.06134575D-20 
  real*8 , intent(inout) :: pot(4,0:ngrid+1),pot_d(4,0:ngrid+1),potzt(0:ngrid+1),potzt2(0:ngrid+1)
  real*8, intent(in)  :: catom(2)
  real*8, intent(in)  :: roff1
  real*8, intent(in)  :: roff2



  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: j1, l, i1, i2,k
  real*8 :: aux1, aux2, aux3, r, r3, a0, b1, b2, b3, b4, som, r4, r5&
       , rbohr, r2, c1, c2, c3, c4, som1
  real*8, dimension(0:5) :: zie
  real*8 :: decal
  !-----------------------------------------------

  interface
subroutine zieg(zie,decal,catom,roff1,roff2,auxe,  pot,ngrid,csive)
   
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------

  implicit none
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------

  real*8  :: auxe
  real*8 , intent(inout) :: zie(0:5)
  real*8 , intent(inout) :: decal
  real*8 , intent(in) :: catom(2)
  real*8 , intent(in) :: roff1
  real*8 , intent(in) :: roff2
  integer :: ngrid
  real*8 , intent(inout) :: pot(4,0:ngrid+1)


  real*8::csive
end subroutine zieg
end interface


  data rbohr/ 0.529D-8/
  !      stop

  !      do l=1,npair
  !         roff1(l)=csive*Int(roff1(l)/csive)
  !         roff2(l)=csive*Int(roff2(l)/csive)
  !      end do
  ! calcul du polynome de reccordement
!csive=csive*1e-8
call zieg(zie,decal,catom,roff1,roff2,auxe,  pot,ngrid,csive)
!  call zieg (zie,decal,catom,roff1,roff2, auxe,ntyp,npair,pot,ngrid,csive,lu_roff_pair,ipotentiel,typ_pot_pair,ipo)
!write(6,*)'zie',zie

  c1 = 0.1818
  c2 = 0.5099
  c3 = 0.2802
  c4 = 0.02817
  l = 0
!  do i1 = 1, ntyp
!     do i2 = i1, ntyp
!        l = ipo(i1,i2)
!        if(typ_pot_pair(l).ne.ipotentiel)cycle
!        if(lu_roff_pair(l).EQV..false.) cycle

        !            k = 0.D0
!  write(6,*)'csive', csive,ngrid
        do j1 = 1, ngrid
           !               k = k+1
           r = j1*csive
           r2 = r*r
           r3 = r2*r
           if (r<=roff1) then              !Potentiel de Ziegler
              !      L=0
              !         L=L+1
              a0 = 0.88534*rbohr/(catom(1)**0.23+catom(2)**0.23)
              b1 = 3.2/a0
              b2 = 0.9423/a0
              b3 = 0.4029/a0
              b4 = 0.2016/a0
              aux2 = auxe*catom(1)*catom(2)

              som = c1*exp((-b1*r))+c2*exp((-b2*r))+c3*exp((-b3*r))+c4*&
                   exp((-b4*r))
              som1 = -c1*b1*exp((-b1*r))-c2*b2*exp((-b2*r))-c3*b3*exp((-b3*r))+c4*&
                    b4*((-b4*r))
!              write(6,*)'a', a0,som,som1
!              write(6,*)'aux2', aux2,decal,r ,som
              pot(1,j1+1) = aux2/r*som+decal
              potzt(j1+1) = aux2/r*som+decal
              potzt2(j1+1) = aux2/r*som
              pot_d(1,j1+1) = -aux2/r2*som + aux2/r*som1
           else                              !Potentiel polynomiale de raccord entre ROFF1 et ROFF2
              if (r<roff2) then
                 a0 = 0.88534*rbohr/(catom(1)**0.23+catom(2)**0.23)
                 b1 = 3.2/a0
                 b2 = 0.9423/a0
                 b3 = 0.4029/a0
                 b4 = 0.2016/a0
                 aux2 = auxe*catom(1)*catom(2)
                 
                 som = c1*exp((-b1*r))+c2*exp((-b2*r))+c3*exp((-b3*r))+c4*&
                      exp((-b4*r))
                 som1 = -c1*b1*exp((-b1*r))-c2*b2*exp((-b2*r))-c3*b3*exp((-b3*r))+c4*&
                      b4*((-b4*r))
                 potzt(j1+1) = aux2/r*som+decal
                 potzt2(j1+1) = aux2/r*som



                 r3 = r2*r
                 r4 = r3*r
                 r5 = r4*r

!                 write(6,*)'rrr',r,r2,r3,r4,r5
!                 write(6,*)'z1',zie(0)
!                 write(6,*)'z1',zie(1),r
!                 write(6,*)'z2',zie(2),r2
!                 write(6,*)'z3',zie(3),r3

                 pot(1,j1) = (zie(5)*r5+zie(4)*r4+zie(3)*r3+zie(2)*r2+zie(1)*r+zie(0))
                 pot_d(1,j1) = 5.d0*zie(5)*r3*r+4.d0*zie(4)*r3+3.d0*zie(3)*r2+2.d0*zie(2)*r+zie(1)
!                 write(6,*)'pot',pot(1,j1)
              else
                 a0 = 0.88534*rbohr/(catom(1)**0.23+catom(2)**0.23)
                 b1 = 3.2/a0
                 b2 = 0.9423/a0
                 b3 = 0.4029/a0
                 b4 = 0.2016/a0
                 aux2 = auxe*catom(1)*catom(2)
                 
                 som = c1*exp((-b1*r))+c2*exp((-b2*r))+c3*exp((-b3*r))+c4*&
                      exp((-b4*r))
                 som1 = -c1*b1*exp((-b1*r))-c2*b2*exp((-b2*r))-c3*b3*exp((-b3*r))+c4*&
                      b4*((-b4*r))
                 potzt(j1+1) = aux2/r*som+decal
                 potzt2(j1+1) = aux2/r*som
              end if
           endif
        end do

!     end do
!  end do

!rite(6,*) potzt
!rite(6,*) potzt2




  return
end subroutine zieg2


! --- Potentiel de Ziegler calcul du polynome de raccordement ---
subroutine zieg(zie,decal,catom,roff1,roff2,auxe,pot,ngrid,csive)
   
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------

  implicit none
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------

  real*8  :: auxe
  real*8 , intent(out) :: zie(0:5)
  real*8 , intent(out) :: decal
  real*8 , intent(in) :: catom(2)
  real*8 , intent(in) :: roff1
  real*8 , intent(in) :: roff2
  integer :: ngrid
  real*8 , intent(inout) :: pot(4,0:ngrid+1)


  real*8::csive


  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: l, i1, i2, l1, l2, ll,k
  real*8 :: sk
  real*8 :: aux3, pi, rbohr, roff12, roff13, roff14, roff15, roff22, &
       roff23, roff25, rdif, rsom, rmid, roff24, eta, auxpi, c1, c2, c3, c4, &
       alpi, a0, b1, b2, b3, b4, abmh, aux1, aux2, r, r2, r3, som, v1&
       , vd1, vdd1, v1mid, vd1mid, damp, ar, v2, vd2, v2mid, x1, xd1, xdd1, &
       x2, xd2, xdd2, vala1, vala2, vala3, valc1, valc2, valc3, valc4, vdd2, &
       r4,r5,r6, r7, r8,Nv2,Nvd2,Nvdd2,Nv2mid,drk
  !-----------------------------------------------
  !   E x t e r n a l   F u n c t i o n s
  !-----------------------------------------------
  real*8 , external :: fac
  !-----------------------------------------------
  data pi/ 3.141592654D0/
  data rbohr/ 0.529D-8/
  ! --- Choix des donnees ---
!  write(6,*)'roff',roff1,roff2,catom,csive
  eta = 0.2
  auxpi = 2.0/sqrt(pi)
  c1 = 0.1818
  c2 = 0.5099
  c3 = 0.2802
  c4 = 0.02817
  alpi = -2.0/sqrt(pi)
  ! --- Fin du choix ---
  ! --- Initialisation des valeurs a trouver ---

!  do i1 = 1, ntyp
!     do i2 = i1, ntyp
!        l = ipo(i1,i2)
!        if(typ_pot_pair(l).ne.ipotentiel)cycle
!        if(lu_roff_pair(l).EQV..false.) cycle
        roff12 = roff1**2
        roff13 = roff1**3
        roff14 = roff1**4
        roff15 = roff1**5
        roff22 = roff2**2
        roff23 = roff2**3
        roff24 = roff2**4
        roff25 = roff2**5
        rdif = roff1-roff2
        rsom = roff1+roff2
        rmid = rsom/2.0
        a0 = 0.88534*rbohr/(catom(1)**0.23+catom(2)**0.23)
        b1 = 3.2/a0
        b2 = 0.9423/a0
        b3 = 0.4029/a0
        b4 = 0.2016/a0

        aux2 = 9.D18*1.602D-19**2*catom(1)*catom(2)
        !potentiel de Ziegler à roff1
        r = roff1
        r2 = roff12
        r3 = roff13
        som = c1*exp((-b1*r))+c2*exp((-b2*r))+c3*exp((-b3*r))+c4*exp((-&
             b4*r))
        v1 = aux2/r*som                      !Potentiel de Ziegler
!        write(6,*)'aB', a0,som
!        write(6,*)'aux2B', aux2,r 

        aux3 = aux2/r*((-c1*exp((-b1*r))*b1)-c2*exp((-b2*r))*b2-c3*exp((&
             -b3*r))*b3-c4*exp((-b4*r))*b4)
        vd1 = (-aux2/r2*som)+aux3            !Derivee de Ziegler
        vdd1 = 2.*aux2/r3*som-2.*aux3/r+aux2/r*(c1*b1*b1*exp((-b1*r))+c2*&
             b2*b2*exp((-b2*r))+c3*b3*b3*exp((-b3*r))+c4*b4*b4*exp((-b4*r)&
             ))                                !Derivee seconde de Ziegler
        !potentiel de Ziegler à mi-chemin de roff1 et roff2
        r = rmid
        som = c1*exp((-b1*r))+c2*exp((-b2*r))+c3*exp((-b3*r))+c4*exp((-&
             b4*r))
        v1mid = aux2/r*som                   !Potentiel de Ziegler
        !                                                !Derivee de Ziegler
        vd1mid = (-aux2/r/r*som)+aux2/r*((-c1*exp((-b1*r))*b1)-c2*exp((-&
             b2*r))*b2-c3*exp((-b3*r))*b3-c4*exp((-b4*r))*b4)

        !potentiel de répulsion standard à roff2
        r = roff2
        sk = r/(csive)
        k=int(sk)
        drk=r-k*csive
        Nv2=pot(1,k)+pot(2,k)*drk+pot(3,k)*drk**2+pot(4,k)*drk**3
!        write(6,*)'POT A ROFF2', r*1d8, NV2
        Nvd2=pot(2,k)+2.0*pot(3,k)*drk+3.0*pot(4,k)*drk**2
        Nvdd2=2.0*pot(3,k)+6.0*pot(4,k)*drk
        v2=Nv2 ; vd2=Nvd2 ; vdd2=Nvdd2
        !potentiel de répulsion standard à mi-chemin de roff1 et roff2
        r  = rmid

        k=Int(r/csive)
        drk=r-k*csive
!                    write(6,*)'k r drk ',k,r,drk
!                    write(6,*)pot(1,k),pot(2,k),pot(3,k),pot(4,k)
        Nv2mid=pot(1,k)+pot(2,k)*drk+pot(3,k)*drk**2+pot(4,k)*drk**3
!        write(6,*)'POT A ROFFint', r*1d8, NV2mid
        v2mid=Nv2mid

!                    write(6,*)'l r1,r2 ',roff1,roff2
!                    write(6,*)' v2, vd2, vdd2 ',v2,vd2,vdd2
!                    write(6,*)'Nv2,Nvd2,Nvdd2 ',Nv2,Nvd2,Nvdd2
        ! 
!                   write(6,*)'v2mid Nv2mid',v2mid,Nv2mid
!                    write(6,*)


        ! --- Decalage du potentiel de Ziegler
!        decal = (-v1mid)+v2mid+0.1D-10
        decal=-(vd1+vd2)*(roff2-roff1)/2.0 +v2-v1
        v1 = v1+decal
        x1 = v1
        xd1 = vd1
        xdd1 = vdd1
        x2 = v2
        xd2 = vd2
        xdd2 = vdd2
        ! --- Calcul du polynome
        vala1 = (xdd1-xdd2)/rdif/2.
        vala2 = (xd1-xd2-roff1*xdd1+roff2*xdd2)/rdif
        vala3 = (x1-x2-roff1*xd1+roff2*xd2+xdd1*roff12/2.-xdd2*roff22&
             /2.)/rdif
        valc1 = rsom
        valc2 = roff12+roff1*roff2+roff22
        valc3 = roff13+roff12*roff2+roff1*roff22+roff23
        valc4 = roff14+roff13*roff2+roff12*roff22+roff1*roff23+roff24
        aux1 = (valc2*vala1-3.*vala3)-(valc1*vala1+vala2)*(6.*valc1*valc2-&
             9.*valc3)/(6.*valc1**2-8.*valc2)
        aux2 = (10.*valc2**2-18.*valc4)-(10.*valc1*valc2-15.*valc3)*(6.*&
             valc1*valc2-9.*valc3)/(6.*valc1**2-8.*valc2)
!        write(6,*)vala1,vala2,vala3
        zie(5) = aux1/aux2
        zie(4) = (valc1*vala1+vala2-zie(5)*(10.*valc1*valc2-15.*valc3))&
             /(6.*valc1**2-8.*valc2)
        zie(3) = (vala1-10.*zie(5)*valc2-6.*zie(4)*valc1)/3.
        zie(2) = (xdd1-20.*zie(5)*roff13-12.*zie(4)*roff12-6.*zie(3&
             )*roff1)/2.
        zie(1) = xd1-5.*zie(5)*roff14-4.*zie(4)*roff13-3.*zie(3)*&
             roff12-2.*zie(2)*roff1
        zie(0) = x1-zie(5)*roff15-zie(4)*roff14-zie(3)*roff13-zie(&
             2)*roff12-zie(1)*roff1
!     end do
!  end do
  return
end subroutine zieg

