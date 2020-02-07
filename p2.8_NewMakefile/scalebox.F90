module scalebox_mod
        use dynalloccell
        use neigcel_mod
        use caltabt_mod
        use period_mod
        use recips_mod 
        use caltabi_mod
        implicit none
        contains
! ******************************************************************
subroutine scalebox(xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use var_pot
  use gen_com_m
  !           (version du 09 juin 2000)
  ! ******************************************************************

  implicit none
  !-----------------------------------------------
  !   G l o b a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  integer  :: ielat(imm)
  integer  :: iwmax(imm)
  integer  :: ityp(imm)
  real(double)  :: xp(3,imm)
  real(double)  :: xpp(3,imm)
  real(double)  :: vp(3,imm)
  real(double)  :: ax(3,imm)
  real(double)  :: fp(3,imm)
  !-----------------------------------------------
  !   L o c a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: i, nb1, nb2, nb3, i1, l,noxn,noyn,nozn
  real(double) :: zlx, zly, zlz, ux, uy, uz,  pi2, fact, fact1&
       , fact2, hk2, ex, ex1, ex2
  !real(double), external :: calcvol
  !-----------------------------------------------
  !
  !
  !
  !
  ! reciprocal space ewald summation



  !debug      write (*,*) 'sub scalebox',it,xp(1,1)
  ! -------------------------------------------------------------
  ! Rescaling des positions
  ! -------------------------------------------------------------
  if (lperiod)    call period

  !debug       write (*,*) 'sub scalebox',it,xp(1,1)


  ! ---------------------------------------------------------------
  ! Recalcul des quantites dependantes de la dimension
  ! ---------------------------------------------------------------
  if (lpr) then
     zl(1) = Sqrt( Sum(at(1:3,1)**2 ) )
     zl(2) = Sqrt( Sum(at(1:3,2)**2 ) )
     zl(3) = Sqrt( Sum(at(1:3,3)**2 ) )
  endif

  volu=calcvol(at(1:3,1),at(1:3,2),at(1:3,3))
  zls2(1:3) = 0.5d0*zl(1:3)

  noxn = int(zl(1)/rumax)
  noyn = int(zl(2)/rumax)
  nozn = int(zl(3)/rumax)
  
#if(ML)
  if (noxn==0) noxn=1
  if (noyn==0) noyn=1
  if (nozn==0) nozn=1

  if (nox==0) nox=1
  if (noy==0) noy=1
  if (noz==0) noz=1
#endif

  if ((noxn==2).or.(noyn==2).or.(nozn==2))then
     noxn=1 ;noyn=1; nozn=1
  end if

  if ((nox.ne.noxn).or.(noy.ne.noyn).or.(noz.ne.nozn).or.((dmtype.eq.9).and.(it==1)))then
     call Deallocatecel
     nox=noxn; noy=noyn; noz=nozn

     if (dmtype.ne.9) then
        if (rang==0) write (6, *) 'IT =',IT,'chgt nox noy noz  = '&
             , nox, noy, noz
     end if

     celsize(1) = zl(1)/float(nox)
     celsize(2) = zl(2)/float(noy)
     celsize(3) = zl(3)/float(noz)
     noxy = nox*noy
     noxyz = nox*noy*noz

     !write(*,*) 'inside scalebox1', nox, noy, noz
     !write(*,*) 'inside scalebox2', noxn, noyn, nozn
     !write(*,*) 'inside scalebox3', noxyz, zl(1), rumax,  im
     natperc= INT(im/noxyz)

     natperc=max(3*natperc,10)
     nvat=10*natperc

     if (dmtype.ne.9) then
        if (rang==0)       write(6,*) ' natperc ', natperc
     end if
     call DynamicalAllocationCell
     call neigcel  

  end if


  call caltabt 

  if (ltabvois.and.(dmtype==9).and.((it==1).or.(mod(it,itetabvois)==0))) then
     call caltabi 
  end if


  if (iewald>0) then

     ! --- Tableaux des troisiemes termes de la sommation d'Ewald ---
     auxe = 23.06134575D-20                  ! en erg.cm (charge electron^2/4*pi*permitivite vide)
     pi2 = pi*pi
     !volu = zl(1)*zl(2)*zl(3)
     volu=calcvol(at(1:3,1),at(1:3,2),at(1:3,3))
     fact = pi2/alpha**2
     fact1 = auxe/2./pi/volu
     fact2 = auxe*2./volu
     do nb1 = -ncoucx, ncoucx
        do nb2 = -ncoucy, ncoucy
           do nb3 = -ncoucz, ncoucz
              if (nb1==0.and.nb2==0.and.nb3==0) cycle
              hk2 = nb1*nb1/zl(1)**2+nb2*nb2/zl(2)**2+nb3*nb3/zl(3)**2
              ex = exp((-hk2*fact))/hk2
              ex1 = ex*fact1
              ex2 = ex*fact2
              tabv3(nb1,nb2,nb3) = ex1
              tabf3(:,nb1,nb2,nb3) = ex2*q(:)
           end do
        end do
     end do

  endif



  return
end subroutine scalebox
end module
