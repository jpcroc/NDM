module trempe_mod
  USE gen_com_m, ONLY:lperiod,tstep,usdh

  implicit none
contains
  subroutine trempe(xp, xpp, vp, fp, ityp,im)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double

    USE var_pot, ONLY:cm,ntyp
!    USE calctemp_mod,only: calctemp

    implicit none
    !-----------------------------------------------
    !   G l o b a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer,allocatable   :: ityp(:)
    real(double),allocatable   :: xp(:,:)
    real(double) ,allocatable  :: xpp(:,:)
    real(double) ,allocatable  :: vp(:,:)
    real(double),allocatable   :: fp(:,:)
    real(double), dimension(ntyp) :: temptyp
    integer::im
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
    !  if (rang==0) write(6,*) 'entree trempe and the mass', cm(:ntyp), tstep
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




    forctot = 0.0
    do i = 1, im
       forctot = forctot+sum(fp(:,i)**2)
    end do
    forctot = sqrt(forctot)
    !  write(6,*)'vp',vp
    !  call calctemp (temptyp) ; write (6,*) 'temp',temptyp

    !  if(rang==0) write (6, *) 'ITTRP ', it, pqotist, forctot, usdh, aux(:ntyp),vp(:,1)
    !if (rang==0) write(6,*) 'PARA-T sortie trempe'
    return
  end subroutine trempe
end module trempe_mod
