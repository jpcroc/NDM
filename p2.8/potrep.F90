
!  **********************************************************
!      Calcul du potentiel repulsif a courte distance
!           de type : Vrep=V0+(rij-rc)**n
!  **********************************************************


subroutine potrep(csive,r0rep,V0rep,ngrid,ntyp,npair)
  !----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  USE gen_com_m, ONLY : ipo, pot,lue_paire

  implicit none
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  integer , intent(in) :: ngrid
  integer :: ntyp
  integer :: npair
  real(double) , intent(in) :: csive
  real(double), intent(in) :: V0rep(npair)
  real(double) , intent(in) :: r0rep(npair)
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: j1,l, i1, i2,k
  real(double) :: r
  real(double), dimension(npair) :: err
  real(double) :: A
  !-----------------------------------------------


  A=1e25

  err(:npair)=0

  l = 0

  do i1 = 1, ntyp
     do i2 = i1, ntyp
        l = ipo(i1,i2)
        if(lue_paire(l)) then
           k = 0
           if(err(l)==20) cycle
           if(r0rep(l)==0.D0) cycle
           ! Calcul du pot. repulsif a courte distance
           do j1 = 1, ngrid
              k = k+1
              r = k*csive
              if (r<=r0rep(l)) then
                 pot(1,l,k) = V0rep(l) + A*(r-r0rep(l))**4
              else
                 pot(1,l,k) = pot(1,l,k)
              endif
           enddo
           err(l)=20
        endif
     enddo
  enddo
  return
end subroutine potrep
