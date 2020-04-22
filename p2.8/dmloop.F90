module dmloop_mod
        USE calfo_mod
        USE dyn_mod
        USE analyse_mod
        USE controle_mod
        USE trempe_mod
        USE sauvegarde_mod
        USE correl_mod
        USE sauveforce_mod
        USE sauveposition_mod
        USE gen_com_m, ONLY:itesauvforce,itesauvposition,lcorrelvp,lfire
        USE atomconfig

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

  USE tab_imm_m
  USE FireModule
#ifdef PARA
  USE mod_para
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
    type(atom_config_d)::atdml
    integer, allocatable ::iwmaxCF(:),indiCF(:)
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
  call ndm2config(atdml,im,imm,xp,fp,vp,xpp,ityp,ielat,num_at_glob,ltabvois,iwmax,indi)
  CALL CalFo(atdml)
!  write(6,*)'dml potist ',potist,atdml%potist
    call config2ndm(atdml,im,imm,xp,fp,vp,xpp,ityp,num_at_glob,ielat,ltabvois,iwmax,indi)
!    iwmax=iwmaxCF
!    indi=indiCF

!  call calfo
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
