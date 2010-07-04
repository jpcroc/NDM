! ************************************************
!           Sous-programme dmloop.f
!          Version MPI du 21 fevrier 2001
! ************************************************

subroutine art_loop_min ()
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use tab_imm_m




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
  !-----------------------------------------------
  !

  ! MPI
  !if (rang==0) write (6, *) '***** PREMIERE ITERATION  ****'

  !      write(6,*)'im',im

1 continue
  it_min = it_min+1

  !      write(6,*)
  !      write(6,*)'***** ITERATION  ****', it

  ! appel de la routine generale des forces

  call calfo
  select case (min_art)

  case(1)

     call trempe (xp, xpp, vp, ax, fp, ielat, iwmax, ityp)

  case (2) 
     
     call tr_fire(vp,fp,xpp,xp)

  case default
     write (6, *) 'ne sait pas quoi faire stop'
     stop
  end select
!-----------------? all the rest ...........

  call analyse 
  ! MPI
  if (rang==0) then
     !     write(6,*)'analyse -> sauvegarde'
     if (itesauv/=0) then
        if (mod(it,itesauv)==0) call sauvegarde 
     endif

     !     write(6,*)'analyse -> sauveposition'
     if (itesauvposition/=0) then
        if (mod(it,itesauvposition)==0) call sauveposition ( it)
     endif
     !     write(6,*)'sauvposition -> control'
  endif                                   ! fin rang=0

  call controle 


  !     write(6,*)' controle ->'


  go to 1


  return

end subroutine dmloop
