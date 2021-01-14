module neb_controle_mod
  USE endrun_mod,only: endrun
  USE dynalloccell,only:
  USE sauveposition_mod,only: sauveposition
  USE recips_mod,only: recips
  USE deftimestep_mod,only: deftimestep
  USE gen_com_m, ONLY:fpstop,fsumstop,tempstop,nebtype,temp,rang,potist,leev,itmax,itetimestep,itetemp,&
       &angst,erg2ev,it
  USE dynalloccell,only:deallocateall
  implicit none
contains


  ! ***********************************************************
  !           Neb control
  !           MCM for JPC 08/02/2007
  ! ***********************************************************

  subroutine neb_controle (ii,xp,fp,im) 
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double

    USE neb_module,only: forctot,formax,nebtest,dragtest

    implicit none
    !-----------------------------------------------
    !   G l o b a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer,intent(in)  :: ii,im
    real(double),allocatable  :: xp(:,:)
    real(double),allocatable  :: fp(:,:)
    !-----------------------------------------------
    !   L o c a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: nacou, i,  iti
    real(double) :: vv, a1, a2, a3, c1, c2, c3
    real(double), dimension(1,3) :: g1,aux
    real(double) :: ltc, ctime, tdev, tcool, epc1, epc2, epc3,massa,tclt
    real(double), dimension(1,3) :: xtr, cv
    real(double) :: fpmax,fpn

    save ltc
    !-----------------------------------------------
    !
    !
    !


    ! last iteration ?
    if (it>=itmax) then
       if (rang==0) write (6, *) '*******Derniere iteration PAS CONVERGE !!**** '
       dragtest=1
       nebtest(ii)=1
       return
    endif

    select case (nebtype)

    case(1)

       if (lEev.EQV..true.) then 
          forctot = sqrt(SUM(fp(1:3,1:im)**2))*erg2eV/angst
          formax  = sqrt(MAXVAL(fp(1,:im)**2+fp(2,1:im)**2+fp(3,1:im)**2))*erg2eV/angst
!                      write(*,'("GC: ",i6,3E20.10)') it,forctot, formax, potist*erg2eV
          if (fpstop>0) then   
             if (formax.le.fpstop) then
                write(6,*)'force par atome  max  ev/Ang ', formax
                write (6, *) 'energie ', potist*erg2eV
                dragtest=1
             end if
          end if
          if (fsumstop>0) then   
             if (forctot.le.fsumstop) then
                write(6,*)'  sqrt ( sum_f F_i^2 ):   ev/Ang ', forctot
                write (6, *) 'energie ', potist*erg2eV
                dragtest=1
             end if
          end if

       else
          forctot = sqrt(SUM(fp(1:3,1:im)**2))
          formax  = sqrt(MAXVAL(fp(1,1:im)**2+fp(2,1:im)**2+fp(3,1:im)**2))
          !            write(*,'("GC: ",i6,3E20.10)') it,forctot, formax, potist
          if (fpstop>0) then   
             if (formax.le.fpstop) then
                write(6,*)'force par atome  max cgs ',formax
                write (6, *) 'energie ', potist

                dragtest=1                  
             end if
          end if

          if (fsumstop>0) then   
             if (forctot.le.fsumstop) then
                write(6,*)'  sqrt ( sum_f F_i^2 ):  cgs  ', forctot
                write (6, *) 'energie ', potist
                dragtest=1  
             end if
          end if

       end if

    case(2)


       nebtest(ii)=0 

       if (lEev.EQV..true.) then 
          forctot = sqrt(SUM(fp(1:3,1:im)**2))*erg2eV/angst
          formax  = sqrt(MAXVAL(fp(1,1:im)**2+fp(2,1:im)**2+fp(3,1:im)**2))*erg2eV/angst
          
!          formaxperp  = sqrt(MAXVAL(fp_par(1,:)**2+fp_par(2,:)**2+     &
!               fp_par(3,:)**2))*erg2eV/angst
!          formaxparl  = sqrt(MAXVAL((fp(1,:)-fp_par(1,:))**2+(fp(2,:)-  &
!               fp_par(2,:))**2+(fp(3,:)-fp_par(3,:))**2))*erg2eV/angst

          if (fpstop>0) then   
             if (formax.le.fpstop) then
                !                  write(6,*)'force par atome  max  ev/Ang ', formax
                !                  write (6, *) 'energie ', potist*erg2eV
                !               write(*,*) 'NEB : image  and force .....:', formax 
                nebtest(ii)=1
             end if
          end if
          if (fsumstop>0) then   
             if (forctot.le.fsumstop) then
                write(6,*)'  sqrt ( sum_f F_i^2 ):   ev/Ang ', forctot
                write (6, *) 'energie ', potist*erg2eV
                nebtest(ii)=1
             end if
          end if

       else
          forctot = sqrt(SUM(fp(1:3,1:im)**2))
          formax  = sqrt(MAXVAL(fp(1,1:im)**2+fp(2,1:im)**2+fp(3,1:im)**2))
          !            write(*,'("GC: ",i6,3E15.5)') it,forctot, formax, potist
          if (fpstop>0) then   
             if (formax.le.fpstop) then
                write(6,*)'force par atome  max cgs ',formax
                write (6, *) 'energie ', potist

                nebtest(ii)=1
             end if
          end if

          if (fsumstop>0) then   
             if (forctot.le.fsumstop) then
                write(6,*)'  sqrt ( sum_f F_i^2 ):  cgs  ', forctot
                write (6, *) 'energie ', potist
                nebtest(ii)=1
             end if
          end if

       end if


    end select






    return
  end subroutine neb_controle
end module neb_controle_mod
