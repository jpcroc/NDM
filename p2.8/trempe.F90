subroutine trempe(xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use var_pot

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
  integer :: i, ic, iti
  real(double), dimension(ntyp) :: aux
  real(double) :: forctot, xprov
  !-----------------------------------------------
  !

  aux(:ntyp) = tstep**2/cm(:ntyp)
  !if (rang==0) write(6,*) 'entree trempe and the mass', cm(:ntyp), tstep
          do i = 1, im
             do ic = 1, 3
                if (vp(ic,i)*fp(ic,i)>0) then
                   xprov = xp(ic,i)-xpp(ic,i)+xp(ic,i)+aux(ityp(i))*fp(ic,i)
                else
                   xprov = xp(ic,i)+fp(ic,i)*aux(ityp(i))
                endif
                vp(ic,i) = (xprov-xpp(ic,i))*usdh
                xpp(ic,i) = xp(ic,i)
                xp(ic,i) = xprov
             end do
          end do

  IF (lperiod) call period


  forctot = 0.0
  do i = 1, im
     forctot = forctot+sum(fp(:,i)**2)
  end do
  forctot = sqrt(forctot)
  !debug if(rang==0) write (6, *) 'ITTRP ', it, potist, forctot, usdh, aux(:ntyp)
  !if (rang==0) write(6,*) 'PARA-T sortie trempe'
  return
end subroutine trempe
