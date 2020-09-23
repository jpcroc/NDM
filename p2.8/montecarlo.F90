module montecarlo_mod
  USE gen_com_m,only:  at, bg, lperiod, timel,& 
                      &tstep, timel,tstep, sig, potist,&
             &at,bg,zl,zls2,nzl,volu,normat
  USE atomconfig,only:atom_config,atom_config_d
  USE period_mod,only: period 
  USE cellconfig, only:cell_config, caltabtC
  USE var_pot,only:ntyp,cm
  USE calfo_mod,only: calfo
  USE T_kind_param_m, ONLY:  double
  USE cryst_to_cart_mod, ONLY: cryst_to_cart
  USE caltabi_mod,only: caltabi
  USE boxconfig,only:box_config,boxconfig2ndm,ndm2boxconfig

  implicit none
  type(box_config)::boxndm

  type(atom_config_d)::config_atom_n !type derive atom_config du systeme a n atomes
  type(atom_config_d)::config_atom_nplus1 !type derive atom_config du systeme a n+1 atomes

  type(cell_config):: cells_n !type derive cell_config du systeme a n atomes
  type(cell_config):: cells_nplus1 !type derive cell_config du systeme a n+1 atomes
integer::  pas_lambda_mc
contains 
  subroutine montecarlo

    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    implicit none
    !-----------------------------------------------
    !   G l o b a l   P a r a m e t e r s
    !-----------------------------------------------
    !type(atom_config_d)::config_atom_n !type derive atom_config du systeme a n atomes
    !type(atom_config_d)::config_atom_nplus1 !type derive atom_config du systeme a n+1 atomes

    !type(cell_config):: cells_n !type derive cell_config du systeme a n atomes
    !type(cell_config):: cells_nplus1 !type derive cell_config du systeme a n+1 atomes
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    real(double), dimension(ntyp) :: aux  !pour les calculs d'acceleration
    real(double) :: x_nplus1, y_nplus1, z_nplus1 !position initiale aleatoire de la N+1eme particule
    !real(double), dimension(3) :: vec_nplus1
    real(double), dimension(3,1) :: cart_vec_nplus1
    real(double) :: lambda_mc !lambda compris entre 0 et 1
    real(double), dimension(3,3) :: sig_n, sig_nplus1
    real(double) :: potist_n, potist_nplus1

    integer :: i
    logical :: lextend   

!initialisation variables
lambda_mc = 0.0
lextend = .true.

lperiod = .true.
call ndm2boxconfig(at,bg,zl,zls2,nzl,volu,normat,boxndm)

!defintion de la boite du syst a N atomes
    call caltabtC(cells_n,config_atom_n,lperiod,bg)
   !call cells_n%print

!definir le systeme a N+1 en tirant une position aleatoire pour le N+1eme atome
   call random_number(x_nplus1)
   call random_number(y_nplus1)
   call random_number(z_nplus1)

   cart_vec_nplus1(1,1) = x_nplus1
   cart_vec_nplus1(2,1) = y_nplus1
   cart_vec_nplus1(3,1) = z_nplus1

   call cryst_to_cart(1,cart_vec_nplus1,at,1) !at vecteur de base de la boite en cm, defini dans gen_com_m

   !copie du syst n dans n+1 et addition de la n+1eme particule
   call config_atom_nplus1%init(config_atom_n%im+1,config_atom_n%imm,config_atom_n%ltabvois)

   config_atom_nplus1%ltabvois=config_atom_n%ltabvois
   call config_atom_n%copy_config(config_atom_nplus1,lrescl=.false.)

   config_atom_nplus1%xp(1:3,config_atom_nplus1%im) = cart_vec_nplus1(1:3,1)
   config_atom_nplus1%fp(1:3,config_atom_nplus1%im) = 0
   config_atom_nplus1%vp(1:3,config_atom_nplus1%im) = 0
   config_atom_nplus1%xpp(1:3,config_atom_nplus1%im) = 0
   config_atom_nplus1%ityp(config_atom_nplus1%im) = 1

!copie de cell puis caltabtC pour redecouper avec la n+1eme particule
  

   !call cells_n%print
   call cells_nplus1%init(cells_n%nox,cells_n%noy,cells_n%noz, cells_n%natperc)
   !write(6,*) 'postinit' 
   !cells_nplus1=cells_n
   call cells_n%copy_cell(cells_nplus1)
   
   !call cells_nplus1%print
   call caltabtC(cells_nplus1,config_atom_nplus1,lperiod,bg)
   !call cells_nplus1%print

!lorsque ltabvois = true, attention, il faut la recalculer pour le syst n+1
!   if (config_atom_nplus1%ltabvois)&
!   &call caltabi(config_atom_nplus1%atom_config,cells_nplus1)


!boucle sur lambda
DO WHILE (lambda_mc < 1)

!calcul des forces des systemes N et N+1
CALL CalFo(sig_n,potist_n,config_atom_n,cells_n,boxndm)
CALL CalFo(sig_nplus1,potist_nplus1,config_atom_nplus1,cells_nplus1,boxndm)


! melange des deux systemes
    DO i=1,config_atom_n%im
     config_atom_nplus1%fp(:,i) = (1-lambda_mc)*config_atom_n%fp(:,i) + lambda_mc*config_atom_nplus1%fp(:,i)
    END DO

! faire le pas de velocity verlet pour determiner les nouvelles forces et positions

    ! step 1 First half-step velocities update, v(t) -> v(t+dt/2)
    timel = timel+tstep
    aux(:ntyp) = tstep/cm(:ntyp)/2.d0
    DO i=1, config_atom_nplus1%im
      config_atom_nplus1%vp(1:3,i) = config_atom_nplus1%vp(1:3,i) + aux(config_atom_nplus1%ityp(i))*config_atom_nplus1%fp(1:3,i)
    END DO
    
    ! step 2  Coordinate update, x(t)-> x(t+dt)
    DO i=1, config_atom_nplus1%im
       config_atom_nplus1%xpp(1:3,i)=config_atom_nplus1%xp(1:3,i)
       config_atom_nplus1%xp(1:3,i) = config_atom_nplus1%xp(1:3,i) + tstep*config_atom_nplus1%vp(1:3,i)
    END DO

    !conditions periodiques ?
    write(6,*) 'config_atom_nplus1%im =' ,config_atom_nplus1%im
    if (lperiod)    call period(config_atom_nplus1%im,config_atom_nplus1%xp,config_atom_nplus1%xpp)
    
    ! repartition des atomes dans la nouvelle boite
    call caltabtC(cells_nplus1,config_atom_nplus1,lperiod,bg)

    ! Force calculation
    CALL CalFo(sig_nplus1,potist_nplus1,config_atom_nplus1,cells_nplus1,boxndm)
    
    ! Second half-step velocities update, v(t+1/2dt) -> v(t+dt)
    DO i=1, config_atom_nplus1%im
      config_atom_nplus1%vp(1:3,i) = config_atom_nplus1%vp(1:3,i) + aux(config_atom_nplus1%ityp(i))*config_atom_nplus1%fp(1:3,i)
    END DO

!repartir les nouvelles positions et forces dans les syst N et N+1
   !Systeme a N

    DO i=1, config_atom_n%im
       config_atom_n%xpp(1:3,i) = config_atom_nplus1%xpp(1:3,i)
       config_atom_n%xp(1:3,i) = config_atom_nplus1%xp(1:3,i)
       config_atom_n%vp(1:3,i) = config_atom_nplus1%vp(1:3,i) 
    END DO

    call caltabtC(cells_n,config_atom_n,lperiod,bg)

!incrementation de lambda
lambda_mc = lambda_mc + 1/pas_lambda_mc

END DO 


  end subroutine montecarlo
end module montecarlo_mod
