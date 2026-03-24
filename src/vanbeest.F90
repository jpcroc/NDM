module dervbeest_mod
  USE gen_com_m, only:uwrt,lwrt
  USE arret_ndm_mod,only:arret_ndm
  USE calerf_mod,only: calerf 
  implicit none 
  contains  
subroutine deriVBEEST(rrep,itdp,l,dp,auxe,alpha,ngrid, &
     ntyp,npair,pau,dip,ro,zz,csive)
  !------------------------------------------------------------
  !   M o d u l e s
  !------------------------------------------------------------
  USE T_kind_param_m, ONLY:  double

  implicit none
  !------------------------------------------------------------
  !   D u m m y   A r g u m e n t s
  !------------------------------------------------------------
  integer :: ngrid
  integer :: ntyp
  integer :: npair
  real(double) :: auxe
  real(double) :: alpha
  real(double) :: pau(npair)
  real(double) :: dip(npair)
  real(double) :: ro(npair)
  real(double) :: zz(npair)
  real(double), intent(in) :: csive

  real(double) :: rrep
  integer , intent(in) :: l
  integer , intent(in) :: itdp
  real(double) :: dp(3)

  !------------------------------------------------------------
  !   L o c a l   V a r i a b l e s
  !------------------------------------------------------------
  integer :: k
  real(double) :: r
  real(double) :: r2,r3,r7
  real(double) :: dampr, ar

  !real(double), external :: derfc
  !------------------------------------------------------------

  do k = -1,1
     r=rrep+k*csive*1.D0/2.D0**itdp
     r2 = r*r
     r3 = r2*r
     r7 = r3*r3*r
     ar = alpha*r
     dampr = derfc(ar)
     dp(k+2)=  -pau(l)/ro(l)*exp((-r)/ro(l)) + 6*dip(l)/r7 - &
          auxe*zz(l)*dampr/r2
  enddo
  return
end subroutine deriVBEEST

! ***********************************************************
! ***********************************************************



!      Calcul du maximum local du pot de Van BEEST
!
!  **********************************************************

subroutine maxVBEEST(rrep,csive,l,auxe,alpha,ngrid,ntyp, &
     npair,pau,dip,ro,zz,convrep)
  !-----------------------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------------------
  USE T_kind_param_m, ONLY:  double, extended

  implicit none
  !-----------------------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------------------
  real(double) , intent(in) :: csive

  integer , intent(in) :: ngrid
  integer :: ntyp
  integer :: npair
  real(double) :: auxe
  real(double) :: alpha
  real(double) :: pau(npair)
  real(double) :: dip(npair)
  real(double) :: ro(npair)
  real(double) :: zz(npair)

  real(double), dimension(3):: dp
  real(double) :: rrep
  integer , intent(in) :: l
  integer :: convrep
  !-----------------------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------------------
  integer :: j1,m,k
  integer :: itdp
  integer :: err
  integer :: istop
  integer :: ierr
  real(double) :: t1,t2

  real(double), parameter :: delta=1.D-13
  integer, parameter :: maxconv=30000

  !-----------------------------------------------------------


  ! initialisation des parametres
  ierr=0
  istop=0
  m=1
  itdp=0
  rrep=2*csive
  k=0
  err=0

  do j1=1,maxconv
     if(err==20) cycle
     call deriVBEEST(rrep,itdp,l,dp,auxe,alpha,ngrid, &
          ntyp,npair,pau,dip,ro,zz,csive)

     t1=abs((dp(1)-dp(2)))
     t2=abs((dp(2)-dp(3)))
     if(((t2<=delta.and.istop==1).or.(t1<=delta.and.istop==1)).and. &
          ierr==1) then

        err=20
        convrep=20
        write(uwrt,*) 'l=',l,'rrep=',rrep
        return

     else
        if(k==ngrid-1.and.istop==0) then
           !               write(uwrt,*) 'pas de maximum local pour l=',l
           rrep=0.D0
           err=20
           convrep=20
           return

        else

           if(j1==maxconv.and.istop==1) then

              write(uwrt,*)'perdu max loc pour l=',l
              write(uwrt,*)'problemes de convergence dans la routine max2VBEEST'
              write(uwrt,*)'Verifiez les parametres'
              call arret_ndm
           endif
        endif
     endif

     if(sign(1.0D0,dp(1)).ne.sign(1.0D0,dp(2))) then
        if(sign(1.D0,dp(2))==-1.D0) then
           itdp=itdp+1
           rrep=rrep - csive*1/2**itdp
           istop=1
           ierr=1
        endif
     else
        if(sign(1.0D0,dp(2)).ne.sign(1.0D0,dp(3))) then
           if(sign(1.D0,dp(3))==-1.D0) then
              itdp=itdp+1
              rrep=rrep + csive*1.D0/2.D0**itdp
              istop=1
              ierr=1
           endif
        else
           k=k+1
           if(istop==1) then
              rrep=rrep-csive*1.D0/(2.D0**(itdp+1))
              ierr=0
           else
              rrep=rrep+csive*1.D0/(2.D0**itdp)
           endif
        endif
     endif
  enddo
  return
end subroutine maxVBEEST

!      Calcul du pot de Van BEEST
!  **********************************************************

subroutine potVBEEST(Vpot,r,l,auxe,alpha,ngrid, &
     ntyp,npair,pau,dip,ro,zz)
  !------------------------------------------------------------
  !   M o d u l e s
  !------------------------------------------------------------
  USE T_kind_param_m, ONLY:  double
  implicit none
  !------------------------------------------------------------
  !   D u m m y   A r g u m e n t s
  !------------------------------------------------------------
  integer , intent(in) :: ngrid
  integer :: ntyp
  integer :: npair
  real(double) :: auxe
  real(double) :: alpha
  real(double) :: pau(npair)
  real(double) :: dip(npair)
  real(double) :: ro(npair)
  real(double) :: zz(npair)

  real(double) , intent(in) :: r
  integer , intent(in) :: l
  real(double) :: Vpot
  !------------------------------------------------------------
  !   L o c a l   V a r i a b l e s
  !------------------------------------------------------------
  real(double) :: r2,r3,r6
  real(double) :: dampr, ar

  !real(double), external :: derfc
  !------------------------------------------------------------

  if(r==0.D0) then
     Vpot=0.D0
  else
     r2 = r*r
     r3 = r2*r
     r6 = r3*r3
     ar = alpha*r
     dampr = derfc(ar)
     Vpot = pau(l)*exp((-r)/ro(l)) - dip(l)/r6 + &
          auxe*zz(l)*dampr/r
  end if
  return

end subroutine potVBEEST
end module dervbeest_mod

