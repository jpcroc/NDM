module force_tersoff_facteurs
  USE gen_com_m, only:uwrt,lwrt
  USE arret_ndm_mod,only:arret_ndm
  USE T_kind_param_m
  USE var_pot, ONLY:  ipotentiel
  implicit none

  real(double),dimension (:), allocatable :: lambda1,lambda2,lambda3,Ater,Bter,psi
  real(double),dimension (:), allocatable :: Rter,Ster, beta,nter,cter, dter, hter,deltater
  !real(double) , external :: fr, fa !fc
  real(double) :: pi=3.141592654D0
contains

  !++++++++++++++++++++++ fonction fc(r) ++++++++++++++++++++++++++++++++++  
  subroutine facteur_amortissement (r, ptyp, fc, dfc)

    implicit none

    real(double), intent(in) :: r
    real(double), intent(out), optional :: fc, dfc
    integer, intent(in) :: ptyp
    real(double)::x,rms,rps,fact
    ! Nouveau fc(r) de Lisa Porter & Ju Li de 89
!    write(uwrt,*)'entree',ipotentiel
    select case (ipotentiel)
    case(13)
       if (ster(ptyp).lt.rter(ptyp))then
          write(uwrt,*)'contradiction entre tersoff.potin et ipotentiel'
          write(uwrt,*)'ipotentiel=',ipotentiel,'ster<rter'
          call arret_ndm
       end if

       IF (present(fc)) then
          if (r<=Rter(ptyp)) fc = 1.
          if (Rter(ptyp)<r .and. r<Ster(ptyp)) &
               fc = 0.5+0.5*cos(pi*(r-Rter(ptyp))/(Ster(ptyp)-Rter(ptyp)))
          if (r>=Ster(ptyp)) fc = 0.
       END IF

       IF (present(dfc)) then
          if (Rter(ptyp)<r .and. r<Ster(ptyp)) then
             dfc = -0.5*pi/(Ster(ptyp)-Rter(ptyp))*sin(pi*(r-Rter(ptyp))/(Ster(ptyp)-Rter(ptyp)))
          else
             dfc = 0.
          end if
       END IF

    case (14)

       if (ster(ptyp).gt.rter(ptyp))then
          write(uwrt,*)'contradiction entre tersoff.potin et ipotentiel'
          write(uwrt,*)'ipotentiel=',ipotentiel,'ster>rter'
          call arret_ndm
       end if
!       fact=50.
       fact=10.
       rms=rter(ptyp)-ster(ptyp)
       rps=rter(ptyp)+ster(ptyp)
       x=fact*(r-rter(ptyp))/ster(ptyp)
       IF (present(fc)) then
!          if (r<=rms)then
!             fc = 1.
!          else if (r>rps) then
!             fc = 0.
!          else
             fc=1./(1.+exp(x))
!write(uwrt,*)r,x,fc
!          END IF	
!	  write(uwrt,*)r,fc
       end IF
       IF (present(dfc)) then
          if (Rms<r .and. r<rps) then
             dfc = -(fact*exp(x))/(ster(ptyp)*(1+exp(x))**2)
          else
             dfc = 0.
          end if
       END IF

       
       
       
    case(15)
       ! ANCIEN FONCTION fc(r) De Tersoff 1988
       IF (present(fc)) then
          if (r<=Rter(ptyp)-Ster(ptyp)) fc = 1.
          if (Rter(ptyp)-Ster(ptyp)<r .and. r<Rter(ptyp)+Ster(ptyp)) &
               fc = 0.5-0.5*sin(0.5*pi*(r-Rter(ptyp))/Ster(ptyp))
          if (r>=Rter(ptyp)+Ster(ptyp)) fc = 0.
!          write(uwrt,*)r,fc
       end IF
       IF (present(dfc)) then
          if (Rter(ptyp)-Ster(ptyp)<r .and. r<Rter(ptyp)+Ster(ptyp)) then
             dfc = -0.25*pi/Ster(ptyp)*cos(0.5*pi*(r-Rter(ptyp))/Ster(ptyp))
          else
             dfc = 0.
          end if
       END IF
    case default
       write(uwrt,*)'quel tersoff ?',ipotentiel
       call arret_ndm
    end select
    
    RETURN

  end subroutine facteur_amortissement

  !+++++++++++++++++++ fonction g(cos theta) ++++++++++++++++++++++++++++
  subroutine facteur_angulaire (costheta, ptyp, g, dg)

    implicit none

    real(double), intent(in) :: costheta
    real(double), intent(out), optional :: g, dg
    integer, intent(in) :: ptyp
    !    real(double) :: d=16.217, h=-0.59825, c=1.0039E+5 

    IF (present(g)) g = 1.+cter(ptyp)**2/dter(ptyp)**2-cter(ptyp)**2/ &
         (dter(ptyp)**2+(hter(ptyp)-costheta)**2)

    IF (present(dg)) dg = -2.*cter(ptyp)**2*(hter(ptyp)-costheta)/ &
         (dter(ptyp)**2+(hter(ptyp)-costheta)**2)**2


    RETURN

  end subroutine facteur_angulaire

  !++++++++++++ fonction fr(r) ++++++++++++++++++
  real(double) function fr (r,A,lamb1)
    USE T_kind_param_m, ONLY:  double
    implicit none
    real(double) , intent(in) :: r, A, lamb1

    fr = A*exp(-lamb1*r)

  end function fr

  !++++++++++++ fonction fa(r) ++++++++++++++++++
  real(double) function fa (r,B,lamb2)
    USE T_kind_param_m, ONLY:  double
    implicit none
    real(double) , intent(in) :: r, B, lamb2

    fa = -B*exp(-lamb2*r)

  end function fa
  !=================== FIN ========================================

end module force_tersoff_facteurs
