! *********************************************************************
subroutine Hcyl 
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use var_pot
  !     Version  du 12 avril 2007
  ! *********************************************************************
  use tab_imm_m
#if(PARA)
  use mod_mpi
#endif

  implicit none
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  real(double) :: xd, yd, zd ! coordonnes du point projete sur la droite generatrice
  integer :: i,j,k,l,m
  real(double):: k1, d2, f2, tinitc, tempsauv, vv

  !-----------------------------------------------
  !   E x t e r n a l    F u n c t i o n s
  !-----------------------------------------------
  real(double) :: tempinstcyl
  if (rang==0) write(6,*)'generation de pointe thermique !!!!!!!!!!'
  open (111,file="cyl.mol")
  ! open (112,file="vcyl.mol")
  allocate (cyl(imm))

  do k=1, im
     cyl(k)=.false.
  end do

  ncyl=0

  do i=1, im
     ! construction du cylindre
     k1=(vdc(1)*(xp(1,i)- pc(1))+vdc(2)*(xp(2,i)- pc(2))+vdc(3)*(xp(3,i)- pc(3)))&
          &/(vdc(1)**2+vdc(2)**2+vdc(3)**2)
     xd=vdc(1)*k1+pc(1)
     yd=vdc(2)*k1+pc(2)
     zd=vdc(3)*k1+pc(3)

     ! calcul de la distance au carre 
     d2=(xp(1,i)-xd)**2+(xp(2,i)-yd)**2+(xp(3,i)-zd)**2

     ! selection des atomes a l'interieur du cylindre
     f2=((xd-pc(1))*vdc(1)+(yd-pc(2))*vdc(2)+(zd-pc(3))*vdc(3))&
          &/sqrt(vdc(1)**2+vdc(2)**2+vdc(3)**2)

     ! write(6,*)i,d
     ! write(6,*)xp(1,i),xp(2,i),xp(3,i)
     ! write(6,*)xd,yd,zd

     if ((d2<=rayonc**2) .and. (f2<=lgc).and.(f2>=0)) then
        ncyl=ncyl+1
        cyl(i)=.true.

     end if
  end do


  write (111,*) ncyl
  at=at*1.d8
  write (111,'(9F12.6)')at(1,1),at(2,1),at(3,1),at(1,2),at(2,2),at(3,2),&
       at(1,3),at(2,3),at(3,3)
  at=at/1.d8 

  ! write (unit=111,fmt='(9F7.2)') 
  do j=1, im
     if (cyl(j).EQV..true.) then
       write (111,*) ty(ityp(j)), xp (1,j)*1d8, xp (2,j)*1d8, xp (3,j)*1d8

     end if
  end do

  ! modifications des vitesses
  tinitc=Ecyl*ev2erg/(3.d0*bk*ncyl)
  tempsauv=tempinstcyl(vp, ityp)
  vv=sqrt(1+tinitc/tempsauv)
  !	write (6,*)'VVV',vv,ncyl,tinitc,tempsauv
  ! write(6,*)'xpP',xpp(1,1),xpp(1,2),xpp(1,3),xpp(1,4)

  do l=1, im
     if (cyl(l).EQV..true.) then
        !		write(6,*) l
        !		xpp(:,l)=xp(:,l)-(xp(:,l)-xpp(:,l))*vv
        vp(:,l)=vp(:,l)*vv
        xpp(:,l) = xp(:,l)-vp(:,l)*tstep
        !        xpp(2,:im) = xp(2,:im)-vp(2,:im)*tstep
        !        xpp(3,:im) = xp(3,:im)-vp(3,:im)*tstep

     end if
  end do
  call period
  !	write(6,*)'xpP',xpp(1,1),xpp(1,2),xpp(1,3),xpp(1,4)

  !do m=1,im
  !	if (cyl(m)==.false.) then
  !	write(112,*) vp(:,m)
  !	else
  !	write(6,*) xp(:,m)
  !	write(112,*) vp(:,m), 'vit modifiee', tinitc, tempsauv
  !	end if
  !end do

end subroutine Hcyl
