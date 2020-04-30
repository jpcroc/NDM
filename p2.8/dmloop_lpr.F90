module dmloop_lpr_mod
  USE analyse_mod,only: analyse
  USE controle_mod,only: controle 
  USE sauvegarde_mod,only: sauvegarde
  USE sauveposition_mod,only: sauveposition
  USE sauveforce_mod,only: sauveforce
  USE gen_com_m, ONLY: itesauvforce, itesauvposition,itesauv,ltnose

  USE atomconfig

#ifdef PARA
  USE recips_mod,only: recips
#endif
  implicit none 
contains
  ! boucle de DM pour velocity Verlet
  ! ************************************************

  subroutine dmloop_lpr
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double
    USE Parrinello_Rahman
    USE Parrinello_Rahman_Nose
    USE tab_imm_m
#ifdef PARA
    USE mod_para
#endif
    implicit none

    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    ! MPI

#ifdef PARA
    real(double)::wbox_tot
    real(double) sigkine_tot(3,3)
    integer :: nb1, nb2, nb3, i1, l,noxn,noyn,nozn
    real(double) :: zlx, zly, zlz, ux, uy, uz,  pi2, fact, fact1&
         , fact2, hk2, ex, ex1, ex2
    !real(double), external :: calcvol

#endif
    type(atom_config_d)::atdml



    if (rang==0) write (6, *) '***** PREMIERE ITERATION  ****'


    ! Initialization -------------------------------------------------------
    IF (lTNose) THEN ! Parrinello-Rahman with Nose thermostat
       call initlprNose(xp,xpp,vp,ityp)
    ELSE ! Parinello-Rahman with Nose-Hoover thermostat or constant energy
       call initlpr(xp, xpp, vp, ax, fp, ielat, iwmax, ityp,num_at_glob)
    END IF

    call analyse 

    ! MD loop -------------------------------------------------------------
1   continue 
    it = it+1
    IF (lTNose) THEN ! Parrinello-Rahman with Nose thermostat
!       call calfo
    call ndm2config(atdml,im,imm,xp,fp,vp,xpp,ityp,ielat,num_at_glob,ltabvois,iwmax,indi)
    CALL CalFo(sig,potist,atdml)
    call config2ndm(atdml,im,imm,xp,fp,vp,xpp,ityp,num_at_glob,ielat,ltabvois,iwmax,indi)

       call prNose(xp,xpp,vp,fp,ityp)

#ifdef PARA
       zl(1) = Sqrt( Sum(at(1:3,1)**2 ) )
       zl(2) = Sqrt( Sum(at(1:3,2)**2 ) )
       zl(3) = Sqrt( Sum(at(1:3,3)**2 ) )
       volu=calcvol(at(1:3,1),at(1:3,2),at(1:3,3))
       zls2(1:3) = 0.5d0*zl(1:3)


    call caltabt(im,xp,ielat) 

       temps_debpara=MPI_Wtime()
       ! Mise a jour des atomes (locaux/frontieres/fantomes) sur tous les processeurs
       call maj_atomes_frt_ftm
       temps_para=temps_para+MPI_Wtime()-temps_debpara


       if (iewald>0) then

          ! --- Tableaux des troisiemes termes de la sommation d'Ewald ---
          auxe = 23.06134575D-20                  ! en erg.cm (charge electron^2/4*pi*permitivite vide)
          pi2 = pi*pi
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

#else

       CALL ScaleBox(xp, xpp, vp, ax, fp, ielat, iwmax, ityp)

#endif
       timel=timel+fNose*tstep
    ELSE ! Parinello-Rahman with Nose-Hoover thermostat or constant energy
       call pr(xp, xpp, vp, ax, fp, ielat, iwmax, ityp,num_at_glob)
       timel=timel+tstep
    END IF

    call analyse 


    if (itesauv.GT.0) then
       if (mod(it,itesauv)==0) call sauvegarde
    endif
    if (itesauvposition.GT.0) then
       if (mod(it,itesauvposition)==0) call sauveposition (it)
    endif
    if (itesauvforce.GT.0) then
       if (mod(it,itesauvforce)==0) call sauveforce (it)
    endif


    call controle

    go to 1

    return
  end subroutine dmloop_lpr

end module dmloop_lpr_mod
