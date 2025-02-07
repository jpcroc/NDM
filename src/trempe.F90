module trempe_mod
  USE gen_com_m, ONLY:lperiod,tstep,usdh,rang
  use atomconfig,only:atom_config_d,atom_config_e
  implicit none
contains
  subroutine trempe(atdml)
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
    class(atom_config_d)::atdml
    integer::im
    !-----------------------------------------------
    !   L o c a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: i, ic
    real(double), dimension(ntyp) :: aux
    real(double) :: forctot, xprov
    !-----------------------------------------------
    !
    im=atdml%im
    aux(:ntyp) = tstep**2/cm(:ntyp)
    select type (atdml)
    class is (atom_config_e)
       if (atdml%lxpp) then 
          do i = 1, im
             do ic = 1, 3
                if (atdml%vp(ic,i)*atdml%fp(ic,i)>0) then
                   xprov = atdml%xp(ic,i)-atdml%xpp(ic,i)+atdml%xp(ic,i)+aux(atdml%ityp(i))*atdml%fp(ic,i)
                else
                   xprov = atdml%xp(ic,i)+atdml%fp(ic,i)*aux(atdml%ityp(i))
                endif
                atdml%vp(ic,i) = (xprov-atdml%xpp(ic,i))*usdh
                atdml%xpp(ic,i) = atdml%xp(ic,i)
                atdml%xp(ic,i) = xprov
             end do
             

          end do
       end if
    end select



    forctot = 0.0
    do i = 1, im
       forctot = forctot+sum(atdml%fp(:,i)**2)
    end do
    forctot = sqrt(forctot)
    !  write(6,*)'vp',vp
    !  call calctemp (temptyp) ; write (6,*) 'temp',temptyp

    !  if(rang==0) write (6, *) 'ITTRP ', it, pqotist, forctot, usdh, aux(:ntyp),vp(:,1)
    !if (rang==0) write(6,*) 'PARA-T sortie trempe'
    return
  end subroutine trempe
end module trempe_mod
