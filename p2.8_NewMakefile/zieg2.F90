!  **********************************************************
!      Calcul du potentiel de Ziegler-Biersack-Littmark
!           
!  **********************************************************
module zieg2_mod
        implicit none 
        contains

subroutine zieg2(pot, pot_d, csive,ngrid, ntyp,npair, catom, roff1, roff2,lu_roff_pair,ipotentiel,typ_pot_pair,ipo)
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double

  implicit none
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  integer, dimension(:,:), pointer  :: ipo                      ! indice des paires d'atomes
  integer, pointer:: typ_pot_pair(:) ! donne le type d'interaction de la paire
  integer , intent(in) :: ngrid,ipotentiel
  integer  :: ntyp
  integer  :: npair
  real(double) , intent(in) :: csive
  real(double)  :: auxe= 23.06134575D-20 
  real(double) , intent(inout) :: pot(4,npair,0:ngrid+1),pot_d(4,npair,0:ngrid+1)
  real(double)  :: catom(ntyp)
  real(double)  :: roff1(npair)
  real(double)  :: roff2(npair)
  logical :: lu_roff_pair(npair)


  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: j1, l, i1, i2,k
  real(double) :: aux1, aux2, aux3, r, r3, a0, b1, b2, b3, b4, som, r4, r5&
       , rbohr, r2, c1, c2, c3, c4, som1
  real(double), dimension(npair,0:5) :: zie
  real(double), dimension(npair) :: decal
  !-----------------------------------------------

!  interface
!subroutine zieg(zie,decal,catom,roff1,roff2,auxe, ntyp, npair, pot,ngrid,csive,lu_roff_pair,ipotentiel,typ_pot_pair,ipo)
!   
!  !-----------------------------------------------
!  !   M o d u l e s
!  !-----------------------------------------------
!  USE T_kind_param_m, ONLY:  double
!  implicit none
!  !-----------------------------------------------
!  !   D u m m y   A r g u m e n t s
!  !-----------------------------------------------
!  integer , intent(in) :: ntyp,ipotentiel
!  integer, pointer:: typ_pot_pair(:) ! donne le type d'interaction de la paire
!  integer, dimension(:,:), pointer  :: ipo                      ! indice des paires d'atomes
!
!  integer , intent(in) :: npair
!  real(double)  :: auxe
!  real(double) , intent(inout) :: zie(npair,0:5)
!  real(double) , intent(inout) :: decal(npair)
!  real(double) , intent(in) :: catom(ntyp)
!  real(double) , intent(in) :: roff1(npair)
!  real(double) , intent(in) :: roff2(npair)
!  integer :: ngrid
!  real(double) , intent(inout) :: pot(4,npair,0:ngrid+1)
!  logical :: lu_roff_pair(npair)
!
!  real(double)::csive
!end subroutine zieg
!end interface


  data rbohr/ 0.529D-8/
  !      stop

  !      do l=1,npair
  !         roff1(l)=csive*Int(roff1(l)/csive)
  !         roff2(l)=csive*Int(roff2(l)/csive)
  !      end do
  ! calcul du polynome de reccordement
  call zieg (zie,decal,catom,roff1,roff2, auxe,ntyp,npair,pot,ngrid,csive,lu_roff_pair,ipotentiel,typ_pot_pair,ipo)


  c1 = 0.1818
  c2 = 0.5099
  c3 = 0.2802
  c4 = 0.02817
  l = 0
  do i1 = 1, ntyp
     do i2 = i1, ntyp
        l = ipo(i1,i2)
        if(typ_pot_pair(l).ne.ipotentiel)cycle
        if(lu_roff_pair(l).EQV..false.) cycle

        !            k = 0.D0
        do j1 = 1, ngrid
           !               k = k+1
           r = j1*csive
           r2 = r*r
           r3 = r2*r
           if (r<=roff1(l)) then              !Potentiel de Ziegler
              !      L=0
              !         L=L+1
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
              pot(1,l,j1+1) = aux2/r*som+decal(l)
              pot_d(1,l,j1+1) = -aux2/r2*som + aux2/r*som1
           else                              !Potentiel polynomiale de raccord entre ROFF1 et ROFF2
              if (r<roff2(l)) then
                 r3 = r2*r
                 r4 = r3*r
                 r5 = r4*r
                 pot(1,l,j1) = (zie(l,5)*r5+zie(l,4)*r4+zie(l,3)*r3+zie&
                      (l,2)*r2+zie(l,1)*r+zie(l,0))
                 pot_d(1,l,j1) = 5.d0*zie(l,5)*r3*r+4.d0*zie(l,4)*r3+3.d0*zie(l,3)*r2+2.d0*zie(l,2)*r+zie(l,1)
              endif
           endif
        end do

     end do
  end do



  return
end subroutine zieg2


! --- Potentiel de Ziegler calcul du polynome de raccordement ---
subroutine zieg(zie,decal,catom,roff1,roff2,auxe, ntyp, npair, pot,ngrid,csive,lu_roff_pair,ipotentiel,typ_pot_pair,ipo)
   
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  implicit none
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  integer , intent(in) :: ntyp,ipotentiel
  integer, pointer:: typ_pot_pair(:) ! donne le type d'interaction de la paire
  integer, dimension(:,:), pointer  :: ipo                      ! indice des paires d'atomes

  integer , intent(in) :: npair
  real(double)  :: auxe
  real(double) , intent(inout) :: zie(npair,0:5)
  real(double) , intent(inout) :: decal(npair)
  real(double) , intent(in) :: catom(ntyp)
  real(double) , intent(in) :: roff1(npair)
  real(double) , intent(in) :: roff2(npair)
  integer :: ngrid
  real(double) , intent(inout) :: pot(4,npair,0:ngrid+1)
  logical :: lu_roff_pair(npair)

  real(double)::csive
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: l, i1, i2, l1, l2, ll,k
  real(double) :: sk
  real(double) :: aux3, pi, rbohr, roff12, roff13, roff14, roff15, roff22, &
       roff23, roff25, rdif, rsom, rmid, roff24, eta, auxpi, c1, c2, c3, c4, &
       alpi, a0, b1, b2, b3, b4, abmh, aux1, aux2, r, r2, r3, som, v1&
       , vd1, vdd1, v1mid, vd1mid, damp, ar, v2, vd2, v2mid, x1, xd1, xdd1, &
       x2, xd2, xdd2, vala1, vala2, vala3, valc1, valc2, valc3, valc4, vdd2, &
       r4,r5,r6, r7, r8,Nv2,Nvd2,Nvdd2,Nv2mid,drk
  !-----------------------------------------------
  !   E x t e r n a l   F u n c t i o n s
  !-----------------------------------------------
  real(double) , external :: fac
  !-----------------------------------------------
  data pi/ 3.141592654D0/
  data rbohr/ 0.529D-8/
  ! --- Choix des donnees ---
  eta = 0.2
  auxpi = 2.0/sqrt(pi)
  c1 = 0.1818
  c2 = 0.5099
  c3 = 0.2802
  c4 = 0.02817
  alpi = -2.0/sqrt(pi)
  ! --- Fin du choix ---
  ! --- Initialisation des valeurs a trouver ---
  l = 0
  do i1 = 1, ntyp
     do i2 = i1, ntyp
        l = ipo(i1,i2)
        if(typ_pot_pair(l).ne.ipotentiel)cycle
        if(lu_roff_pair(l).EQV..false.) cycle
        roff12 = roff1(l)**2
        roff13 = roff1(l)**3
        roff14 = roff1(l)**4
        roff15 = roff1(l)**5
        roff22 = roff2(l)**2
        roff23 = roff2(l)**3
        roff24 = roff2(l)**4
        roff25 = roff2(l)**5
        rdif = roff1(l)-roff2(l)
        rsom = roff1(l)+roff2(l)
        rmid = rsom/2.0
        a0 = 0.88534*rbohr/(catom(i1)**0.23+catom(i2)**0.23)
        b1 = 3.2/a0
        b2 = 0.9423/a0
        b3 = 0.4029/a0
        b4 = 0.2016/a0

        aux2 = 9.D18*1.602D-19**2*catom(i1)*catom(i2)
        !potentiel de Ziegler a roff1
        r = roff1(l)
        r2 = roff12
        r3 = roff13
        som = c1*exp((-b1*r))+c2*exp((-b2*r))+c3*exp((-b3*r))+c4*exp((-&
             b4*r))
        v1 = aux2/r*som                      !Potentiel de Ziegler
        aux3 = aux2/r*((-c1*exp((-b1*r))*b1)-c2*exp((-b2*r))*b2-c3*exp((&
             -b3*r))*b3-c4*exp((-b4*r))*b4)
        vd1 = (-aux2/r2*som)+aux3            !Derivee de Ziegler
        vdd1 = 2.*aux2/r3*som-2.*aux3/r+aux2/r*(c1*b1*b1*exp((-b1*r))+c2*&
             b2*b2*exp((-b2*r))+c3*b3*b3*exp((-b3*r))+c4*b4*b4*exp((-b4*r)&
             ))                                !Derivee seconde de Ziegler
        !potentiel de Ziegler a mi-chemin de roff1 et roff2
        r = rmid
        som = c1*exp((-b1*r))+c2*exp((-b2*r))+c3*exp((-b3*r))+c4*exp((-&
             b4*r))
        v1mid = aux2/r*som                   !Potentiel de Ziegler
        !                                                !Derivee de Ziegler
        vd1mid = (-aux2/r/r*som)+aux2/r*((-c1*exp((-b1*r))*b1)-c2*exp((-&
             b2*r))*b2-c3*exp((-b3*r))*b3-c4*exp((-b4*r))*b4)

        !potentiel de répulsion standarda roff2
        r = roff2(l)
        sk = r/csive
        k=int(sk)
        drk=r-k*csive
        Nv2=pot(1,l,k)+pot(2,l,k)*drk+pot(3,l,k)*drk**2+pot(4,l,k)*drk**3
        Nvd2=pot(2,l,k)+2.0*pot(3,l,k)*drk+3.0*pot(4,l,k)*drk**2
        Nvdd2=2.0*pot(3,l,k)+6.0*pot(4,l,k)*drk
        v2=Nv2 ; vd2=Nvd2 ; vdd2=Nvdd2
        !potentiel de répulsion standarda mi-chemin de roff1 et roff2
        r  = rmid

        k=Int(r/csive)
        drk=r-k*csive
        !            write(6,*)'k r drk ',k,r,drk
        !            write(6,*)pot(1,l,k),pot(2,l,k),pot(3,l,k),pot(4,l,k)
        Nv2mid=pot(1,l,k)+pot(2,l,k)*drk+pot(3,l,k)*drk**2+pot(4,l,k)*drk**3
        v2mid=Nv2mid

        !            write(6,*)'l r1,r2 ',l,roff1(l),roff2(l)
        !            write(6,*)' v2, vd2, vdd2 ',v2,vd2,vdd2
        !            write(6,*)'Nv2,Nvd2,Nvdd2 ',Nv2,Nvd2,Nvdd2
        ! 
        !           write(6,*)'v2mid Nv2mid',v2mid,Nv2mid
        !            write(6,*)


        ! --- Decalage du potentiel de Ziegler
        decal(l) = (-v1mid)+v2mid+0.1D-10
        decal(l)=-(vd1+vd2)*(roff2(l)-roff1(l))/2.0 +v2-v1
        v1 = v1+decal(l)
        x1 = v1
        xd1 = vd1
        xdd1 = vdd1
        x2 = v2
        xd2 = vd2
        xdd2 = vdd2
        ! --- Calcul du polynome
        vala1 = (xdd1-xdd2)/rdif/2.
        vala2 = (xd1-xd2-roff1(l)*xdd1+roff2(l)*xdd2)/rdif
        vala3 = (x1-x2-roff1(l)*xd1+roff2(l)*xd2+xdd1*roff12/2.-xdd2*roff22&
             /2.)/rdif
        valc1 = rsom
        valc2 = roff12+roff1(l)*roff2(l)+roff22
        valc3 = roff13+roff12*roff2(l)+roff1(l)*roff22+roff23
        valc4 = roff14+roff13*roff2(l)+roff12*roff22+roff1(l)*roff23+roff24
        aux1 = (valc2*vala1-3.*vala3)-(valc1*vala1+vala2)*(6.*valc1*valc2-&
             9.*valc3)/(6.*valc1**2-8.*valc2)
        aux2 = (10.*valc2**2-18.*valc4)-(10.*valc1*valc2-15.*valc3)*(6.*&
             valc1*valc2-9.*valc3)/(6.*valc1**2-8.*valc2)
        zie(l,5) = aux1/aux2
        zie(l,4) = (valc1*vala1+vala2-zie(l,5)*(10.*valc1*valc2-15.*valc3))&
             /(6.*valc1**2-8.*valc2)
        zie(l,3) = (vala1-10.*zie(l,5)*valc2-6.*zie(l,4)*valc1)/3.
        zie(l,2) = (xdd1-20.*zie(l,5)*roff13-12.*zie(l,4)*roff12-6.*zie(l,3&
             )*roff1(l))/2.
        zie(l,1) = xd1-5.*zie(l,5)*roff14-4.*zie(l,4)*roff13-3.*zie(l,3)*&
             roff12-2.*zie(l,2)*roff1(l)
        zie(l,0) = x1-zie(l,5)*roff15-zie(l,4)*roff14-zie(l,3)*roff13-zie(l&
             ,2)*roff12-zie(l,1)*roff1(l)
     end do
  end do
  return
end subroutine zieg
end module
