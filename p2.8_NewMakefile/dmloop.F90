module dmloop_mod
        use calfo_mod
        use dyn_mod
        use analyse_mod
        use controle_mod
        use trempe_mod
        use sauvegarde_mod
        use correl_mod
        use sauveforce_mod
        use sauveposition_mod
        implicit none
        contains
! ************************************************
!           Sous-programme dmloop.f
!          Version MPI du 21 fevrier 2001
! ************************************************

subroutine dmloop
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use tab_imm_m
  USE FireModule
#ifdef PARA
  use mod_para
#endif

  implicit none
  !-----------------------------------------------
  !   G l o b a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: i, iti
  REAL(double) :: fire_dt, fire_alph
  INTEGER :: fire_nstep
  !-----------------------------------------------
  !

  ! MPI
  if (rang==0) write (6, *) '***** PREMIERE ITERATION  ****'

  ! Initialization
  IF ((dmtype.EQ.2).AND.lFire) THEN
          CALL init_trempe_fire(fire_dt, fire_nstep, fire_alph)
  END IF

  !      write(6,*)'im',im
1 continue
  it = it+1

  !      write(6,*)
  !      write(6,*)'***** ITERATION  ****', it

  ! appel de la routine generale des forces

  call calfo
  select case (dmtype)

  case(1)

     call dyn 
     if (lcorrelvp) call correlvp(xp,xpp,vp,ax,fp,ityp)

  case (2) 
          IF (lFire) THEN
                  call trempe_fire (xp, xpp, vp, ax, fp, ielat, iwmax, ityp, &
                        fire_dt, fire_nstep, fire_alph)
          ELSE
                  call trempe (xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
          END IF

  case default
     write (6, *) 'ne sait pas quoi faire stop'
     stop
  end select


  call analyse 
  ! MPI
  if (rang==0) then
!          write(6,*)'analyse -> sauvegarde'
     if (itesauv.GT.0) then
        if (mod(it,itesauv)==0) call sauvegarde 
     endif

!          write(6,*)'analyse -> sauveposition'
     if (itesauvposition.GT.0) then
        if (mod(it,itesauvposition)==0) call sauveposition ( it)
     endif
     if (itesauvforce.GT.0) then
        if (mod(it,itesauvforce)==0) call sauveforce ( it)
     endif
!          write(6,*)'sauvposition -> control'
  endif                                   ! fin rang=0

  call controle 


  !     write(6,*)' controle ->'


  go to 1


  return

end subroutine dmloop
end module
