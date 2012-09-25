module teledyn_module
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use jqmod
  use random_art
  use lanczos_defs
  !-----------------------------------------------
  !   _c -> courant
  !   _s -> selection
  !   _d -> depart

      integer,     dimension(:),allocatable,save  :: ielat_c   ,ielat_s   ,ielat_d
      integer,     dimension(:),allocatable,save  :: iwmax_c   ,iwmax_s   ,iwmax_d
      integer,     dimension(:),allocatable,save  :: ityp_c	 ,ityp_s    ,ityp_d
      integer,     dimension(:),allocatable,save  :: ipovois
      real(double),dimension(:,:),allocatable,save  :: xp_c	 ,xp_s	    ,xp_d
      real(double),dimension(:,:),allocatable,save  :: xpp_c	 ,xpp_s     ,xpp_d
      real(double),dimension(:,:),allocatable,save  :: vp_c	 ,vp_s	    ,vp_d
      real(double),dimension(:,:),allocatable,save  :: ax_c	 ,ax_s	    ,ax_d
      real(double),dimension(:,:),allocatable,save  :: fp_c	 ,fp_s	    ,fp_d

      real(double),dimension(:),allocatable,save    :: rec_work , rec_para, rec_correction
      real(double),dimension(:),allocatable,save    :: rec_q6
      real(double),dimension(:,:),allocatable       :: rec_energy
      real(double),dimension(:,:),allocatable,save  :: m_i
      real(double),dimension(:,:),allocatable,save  :: gau
      real(double),dimension(:,:),allocatable,save  :: sig_i
      real(double),dimension(:,:),allocatable,save  :: rga_i
      real(double),dimension(:,:),allocatable,save  :: fadd
      real(double),dimension(:),allocatable,save    :: xbarini
      real(double),dimension(:),allocatable,save    :: d2vois
      real(double),dimension(:,:),allocatable,save  :: xpvois
      real(double),dimension(:,:),allocatable,save  :: xpvoisini
      real(double),dimension(:),allocatable,save    :: xtransla
!      real(double),dimension(:),allocatable,save    :: wq4

      integer                                       :: it_teledyn
      integer                                       :: lupout 
      integer                                       :: luwout
      integer                                       :: lufout 
      integer                                       :: lusout
      integer                                       :: luvout
      integer                                       :: luhout
      
      integer                                       :: lucout
      integer                                       :: ncomptinter
      integer                                       :: idistance
      integer :: iteration
      integer :: iloop
!      integer :: nchemin
      integer :: nfrequence,im1

      real(double),dimension(:),allocatable,save    :: tmass_teledyn
      real(double) ::  tstep_teledyn,usdh_teledyn, convert_teledyn
      real(double) :: sigma
      real(double) :: xdist0
      real(double) :: xdist1
      real(double) :: kappas2 , potistadd
      real(double) :: betaq,barbetaq
      real(double) :: barbeta,barbetaE
      real(double) :: betaweff
      real(double) :: Ecinetique
      real(double) :: Ecinetique0
      real(double) :: dlambda
      real(double) :: lambdax
      real(double) :: lambday
      real(double) :: lambdaz
      real(double) :: deltawork
      real(double) :: work
      real(double) :: deltaU
      real(double) :: xbar0
      real(double) :: m_tot
      real(double) :: T,Tav
      real(double) :: potist0
      real(double) :: dlambdai,dlambdaf
      real(double) :: rdist,rdist2
      real(double) :: x111
      real(double) :: dis2protect
      character*80 :: fnampout
      character*80 :: fnamwout
      character*80 :: fnamfout
      character*80 :: fnamsout
      character*80 :: fnamvout
      character*80 :: fnamhout
      character*80 :: fnamcin
      logical      :: ldeter
      logical      :: ldistance
      logical      :: lta
      logical      :: lvac
      logical      :: lljc
      logical      :: lq6
!     parametre du lennard-jones
      character*80 :: fnamfin
      character*80 :: fnamhin

      integer N
      !PARAMETER(N=256)
      !PARAMETER(N=864)
      PARAMETER(N=38)
      integer nfenetre,lufin,luhin
!      integer 
      parameter(nfenetre=100)
      real(double) :: rsig_lj      
      real(double) :: sig_lj    
      real(double) :: eps_4lj
      real(double)    xij(N,N)
      real(double)    yij(N,N)
      real(double)    zij(N,N)
      real(double)    dij(N,N)
      real(double)    FQ4X(N),FQ4Y(N),FQ4Z(N)
      real(double)    FQ6X(N),FQ6Y(N),FQ6Z(N)
      real(double) :: xpq4,xpq6
      real(double) :: xmaxq6
      real(double) :: pot_auxi(-10:nfenetre+10)
      real(double) :: pot_auxiliary
      real(double) :: contour_auxiliary
      real(double) :: kappaE,kappaEs2
      real(double) :: potistaddE
      real(double) :: alphadd(2)
      real(double) :: massadd(2)
      real(double) :: pot_contour(-10:nfenetre+10,0:nfenetre)
      real(double) :: energie_min,energie_max
      real(double) :: bargamma,seed
      real(double) :: xtempmin

real (double),dimension(N) :: fxx
real(double),dimension(N) :: fyy
real(double),dimension(N) :: fzz
real(double),dimension(N) :: vsecx
real(double),dimension(N) :: vsecy
real(double),dimension(N) :: vsecz



contains

  subroutine allocate_teledyn
  implicit none
    allocate ( ielat_c(imm),iwmax_c(imm), ityp_c(imm),xp_c(3,imm),  xpp_c(3,imm), &
               vp_c(3,imm), ax_c(3,imm), fp_c(3,imm),                             &
               ielat_s(imm),iwmax_s(imm), ityp_s(imm),xp_s(3,imm),  xpp_s(3,imm), &
               vp_s(3,imm), ax_s(3,imm), fp_s(3,imm),                             &
  	       ielat_d(imm),iwmax_d(imm), ityp_d(imm),xp_d(3,imm),  xpp_d(3,imm), &
               vp_d(3,imm), ax_d(3,imm), fp_d(3,imm),                             &
               m_i(3,imm),gau(6,imm),sig_i(3,imm),rga_i(3,imm),fadd(3,imm),       &
               xbarini(3))
 
  allocate ( rec_work(0:niteration),rec_para(0:niteration),tmass_teledyn(imm),rec_q6(0:niteration) )
  allocate ( rec_energy(0:niteration,2),rec_correction(0:niteration))
  allocate (ipovois(15), d2vois(15),xpvois(3,15),xpvoisini(3,15),xtransla(3))

  end subroutine allocate_teledyn
  
  subroutine ndm_into_depart (xp, xpp, vp, ax, fp, ielat, iwmax, ityp)

    implicit none
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer  :: ielat(imm)
    integer  :: iwmax(imm)
    integer  :: ityp(imm)
    real(double)  :: xp(3,imm)
    real(double)  :: xpp(3,imm)
    real(double)  :: vp(3,imm)
    real(double)  :: ax(3,imm)
    real(double)  :: fp(3,imm)

    ielat_d   (:) = ielat   (:)
    iwmax_d   (:) = iwmax   (:)
    ityp_d    (:) = ityp    (:)
    xp_d (:,:)    = xp (:,:)
    xpp_d(:,:)    = xpp(:,:)
    vp_d (:,:)    = vp (:,:)
    ax_d (:,:)    = ax (:,:)
    fp_d (:,:)    = fp (:,:)


  end subroutine ndm_into_depart


 subroutine control_output(xp)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    implicit none
    !-----------------------------------------------
!    real(double)  :: denominant,numerat
    !-----------------------------------------------
    real(double)  :: xp(3,imm)
!        call angular_velocities(xp,vp,ityp)
!        Tav=Tav+T

   betaweff= potistadd/(bk*text_teledyn) + ( potist-potist0+Ecinetique - Ecinetique0)/(bk*text_teledyn)-betaq-barbetaq
   
   IF ((MOD(ITERATION,nfrequence).EQ.0).and.(MOD(iloop,2).EQ.1)) THEN

      write(luwout,'(e18.7,e18.7,e18.7,e18.7,e18.7,e18.7)') betaweff,xp(1:3,1),xp(1,im+1),work
      write(luvout,'(e18.7,e18.7,e18.7,e18.7,e18.7,e18.7,e18.7,e18.7,e18.7,e18.7,e18.7,e18.7,e18.7,e18.7,e18.7)') d2vois(1:15)
   ENDIF
   
   IF (MOD(ITERATION,nfrequence*10).EQ.0) THEN
      write(*,*) ' betaweff',betaweff , potist/(bk*text_teledyn),betaq,barbetaq
      write(*,*) ' Ecinetique/(bk*text_teledyn)' , Ecinetique/(bk*text_teledyn)
      write(*,*) ' T/(bk*text_teledyn) Tav/(bk*text_teledyn)  = ',T/(bk*text_teledyn),Tav/(bk*text_teledyn)/real(iteration)
      write(*,*) 'potist0 ',potist0/(bk*text_teledyn),potistadd/(bk*text_teledyn)
   ENDIF
   
   IF (ITERATION.EQ.niteration) THEN
      write(6,*) ' potist  potistadd  ',potist/(bk*text_teledyn) , potistadd/(bk*text_teledyn)
      write(6,*) ' betaq   barbetaq   ',betaq,barbetaq   , Ecinetique/(bk*text_teledyn)
      write(6,*) ' betaweff  ', betaweff ,one/(bk*text_teledyn)
!      write(6,*) ' ncomptinter  ',ncomptinter
      write(6,*) ' taux d''acceptation  ',ncomptinter/real(ITERATION)
      write(6,*) ' taux de dépassement  ',idistance/real(ITERATION)
   ENDIF
   
   return
 end subroutine control_output

  subroutine genere_bruit 
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
   implicit none
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    !-----------------------------------------------
   real(double) :: u1,u2,b1
   real(double):: deriv(6)
   integer   :: ic,iatom,i
    !-----------------------------------------------
   do i=1,im+1
      do ic=1,6
         call random_number(u1)
         call random_number(u2)
         b1=sqrt(-2.*log(u1))*cos(2.*pi*u2)
         gau(ic,i) = b1
      enddo
   enddo

   gau(1:3,1:im+1) = gau(1:3,1:im+1)*sig_i(1:3,1:im+1)
   gau(4:6,1:im+1) = gau(4:6,1:im+1)*sig_i(1:3,1:im+1)

   do ic=1,6
      deriv(ic)    = sum(gau(ic,1:im))/real(im)
      gau(ic,1:im) = gau(ic,1:im)-deriv(ic)
   enddo

  end subroutine genere_bruit


subroutine calfoljc(xp,fp) !calcul des positions reciproques,des forces et qppelle q4 et q6 à partir des positions 
  !  USE T_kind_param_m, ONLY:  double
  !  use gen_com_m
  implicit none
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
 
  real(double)  :: xp(3,imm)
  real(double)  :: fp(3,imm)
  real(double)  :: c(3)
  real(double)  :: r2,r2sq,potemp,phu
  real(double)  :: frap(1:3)
  real(double)  :: fraptot(1:3),xbar(1:3)
  real(double)  :: f(1:3)
  real(double)  :: r6,r8,r12,r14
  real(double)  :: potistlj
  real(double)  :: ra(1:3)
  real(double)  :: q4,q6
  integer  :: i,j,icri,ic
   real(double)    xij(im,im)
   real(double)    yij(im,im)
   real(double)    zij(im,im)
   real(double)    dij(im,im)
  potistlj = zero
  fp=zero
  fraptot=zero
  icri = 0

  do ic=1,3
     xbar(ic)=sum(xp(ic,1:im))/real(im)
  enddo
  
  do 10 i=1,im
     
     !     anti-évaporation 
     c(1:3)=xp(1:3,i)-xbar(1:3)
     r2=sum(c(1:3)**2)
     r2sq=sqrt(r2)
     !      write(*,*) 'rsig_lj ',rsig_lj
     if (r2sq.gt.rsig_lj) then 
        potemp = (r2sq-rsig_lj)**thr
        potistlj = potistlj + potemp
        phu    = -thr*(r2sq-rsig_lj)**two/r2sq
        frap(1:3)    = phu*c(1:3)
        fraptot(1:3) = phu*c(1:3) + fraptot(1:3)
        f(1:3)=fp(1:3,i)+frap(1:3)
        !         write(*,*) ' r2sq    ', r2sq
        icri=2
     else
        f(1:3)=fp(1:3,i)
     endif

! --------------------------------------------------------
!                  ouverture boucle sur j
! --------------------------------------------------------
      do 30 j=i+1,im
      c(1:3)=xp(1:3,i)-xp(1:3,j)
      r2=sum(c(1:3)**2)
      r2sq=sqrt(r2)
      r6=(sig_lj/r2sq)**6
      r8=r6/r2
      r12=(r6)**2
      r14=r12/r2
      potemp = eps_4lj*(r12-r6) 
      potistlj = potistlj + potemp

      phu = eps_4lj*(12.0*r14 - 6.0*r8)  

! --------------------------------------------------------
!            calcul de la force vectorielle
! --------------------------------------------------------
      ra(1:3)=phu*c(1:3)
      f(1:3) = f(1:3) + ra(1:3)
      fp(1:3,j)=fp(1:3,j)-ra(1:3)
      pist=pist+phu*r2 ! deuxieme partie du viriel

!---------------------------------------------------------
!            tabulation pour Q4 et grad_Q4  
! --------------------------------------------------------
         xij(i,j)=-c(1)
         yij(i,j)=-c(2)
         zij(i,j)=-c(3)
         dij(i,j)=sqrt(xij(i,j)**2+yij(i,j)**2+zij(i,j)**2)
         xij(j,i)=-xij(i,j)
         yij(j,i)=-yij(i,j)
         zij(j,i)=-zij(i,j)
         dij(j,i)= dij(i,j)

   30 continue
   20 continue
      fp(1:3,i)=f(1:3)
   10 continue

! ----------------------- fin du calcul -----------------------

      potist=potistlj

   !xpq4 = q4(im,xij,yij,zij,dij)
   !call fq4(im,xij,yij,zij,dij,fq4x,fq4y,fq4z)

   !if (lq6) then
    !  xpq6 = q6(im,xij,yij,zij,dij)
!   call fq6(im,xij,yij,zij,dij,fq6x,fq6y,fq6z)
   !endif


   do ic=1,3
      fp(ic,1:im) = fp(ic,1:im)- fraptot(ic)/real(im)
   enddo

!   if (icri.gt.1) then 
!       write(*,*) ' calfo sum fp   ',sum(fp(1,1:im)),sum(fp(2,1:im)),sum(fp(3,1:im))
!       write(*,*) ' calfo fp en 1  ',sum(fp(1,1:1)),sum(fp(2,1:1)),sum(fp(3,1:1))
!    endif
!   write(*,*) 'xpq4 ',xpq4
   !call calpot_auxi(xpq4)

end subroutine calfoljc

subroutine calfoljq4(xp,fp)
   integer ic,i
   real(double)  :: xp(3,imm)
   real(double)  :: fp(3,imm)
    !real(double)    FQ4X(im),FQ4Y(im),FQ4Z(im)
   !   real(double)    FQ6X(N),FQ6Y(N),FQ6Z(N)
   real(double)  :: xrel(1:2)
   real(double)  :: xadd(3)
!    real(double)  :: potistadd

   xadd(1:2)    =  xp(1:2,im+1)

   xrel(1)      = -kappa*(xpq4-xadd(1))    ! moins le gradient = la force
   potistadd    =  kappas2*(xpq4-xadd(1))**2

   fadd(1,1:im) =    xrel(1)*fq4x(1:im) ! fq4 est le gradient donc égal à  moins  une force
   fadd(2,1:im) =    xrel(1)*fq4y(1:im)
   fadd(3,1:im) =    xrel(1)*fq4z(1:im)

!   fadd(1,im1)  =  xrel(1)

   if (lta) then 
   xrel(2)        =  -kappaE*(potist-xadd(2))
   potistaddE     =  kappaEs2*(potist-xadd(2))**2
   
   fadd(1:3,1:im) = fadd(1:3,1:im) - xrel(2)*fp(1:3,1:im)! fp est moins un gradient 
   fadd(1:2,im1)  = - xrel(1:2) 
   endif

end subroutine calfoljq4

subroutine calpot_auxi (xpq4) ! calc potenziale ausiliario. In realtà l'unica csa che serve é la definiwzione di ifenetre per l'istogramma
   integer icfen,icener,i
   real(double)  :: xp(3,imm)
   real(double)  :: xrel(3)
   real(double)  :: xadd(3)
   real(double)  :: xbar(3)
   real(double)  :: rdist,rdist2
   real(double)  :: frsr
   real(double)  :: lambdaq4,lambdaE
   real (double), dimension (niteration)::w4
   integer :: ifenetre,ienergy
   real (double) ::xpq4

  

   ifenetre = int((xpq4/x111)*real(nfenetre))
   lambdaq4 = (xpq4/x111)*real(nfenetre)-real(ifenetre)

   if (ifenetre.gt.nfenetre+10) ifenetre = nfenetre+10
   icfen = ifenetre+1 
   if (icfen.gt.nfenetre+10) icfen = nfenetre+10
   if (ifenetre.lt.-10) ifenetre = -10

   pot_auxiliary = pot_auxi(ifenetre)*(one-lambdaq4) + pot_auxi(icfen)*lambdaq4

!   write(*,*) 'pot_auxiliary ',pot_auxiliary,pot_auxi(ifenetre),pot_auxi(ifenetre+1)

!   ifenetre = nint((xpq4/x111)*real(nfenetre))
!   if (ifenetre.gt.nfenetre+10) ifenetre = nfenetre+10
!   if (ifenetre.lt.-10) ifenetre = -10

   ienergy  = int((potist-energie_min)/(energie_max-energie_min)*real(nfenetre))
   lambdaE = ((potist-energie_min)/(energie_max-energie_min))*real(nfenetre)-real(ienergy)

   if (ienergy.gt.nfenetre) ienergy = nfenetre
   if (ienergy.lt.0) ienergy = 0
   icener = ienergy+1 
   if (icener.gt.nfenetre) icener = nfenetre

   contour_auxiliary = (pot_contour(ifenetre,ienergy)*(one-lambdaq4) + pot_contour(icfen,ienergy)*lambdaq4)*(one-lambdaE) +(pot_contour(ifenetre,icener)*(one-lambdaq4) + pot_contour(icfen,icener)*lambdaq4)*lambdaE


end subroutine calpot_auxi

Subroutine langevin_ljc (xp, vp, fp,dt)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    implicit none
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------

    integer ic,i
    real(double)  :: xp(3,N)
    real(double)  :: xpp(3,N)
    real(double)  :: vp(3,N)
    real(double)  :: fp(3,N)
    real(double)  :: pp(3,N)
    real(double) :: xbar(3)
    real(double) :: Ecin0
    real(double) :: Ecin1
    real(double) :: Ecin3
    !real(double) :: Ecin4
    real(double) :: erreur
    real(double) :: xalea
    real(double) :: xprob, dt
    real(double) :: Ecinetique_d

   !-----------------------------------------------
    integer   :: ic_local,iatom
   !-----------------------------------------------

    im=N
   ! call calfoljq4(xp,fp)
    deltawork = potistadd
    !xp(1,im1) =dlambdai+(dlambdaf-dlambdai)*real(iteration)/real(niteration)
    !call calfoljq4(xp,fp)
    !deltawork = potistadd -deltawork
    !work       = work + deltawork/(bk*text_teledyn)
    call genere_bruit

    !Ecin0 = zero
    !Ecin1 = zero
    !Ecin3 = zero
   ! Ecin4 = zero

    pp(1:3,1:im) = vp(1:3,1:im)!*m_i(1:3,1:im)

   ! do ic =1,3
    !   Ecin0 = Ecin0   + DOT_PRODUCT(pp(ic,1:im),vp(ic,1:im))/two
   ! enddo

   ! xp_d(1:3,1:im)=xp(1:3,1:im)
    !vp_d(1:3,1:im)=vp(1:3,1:im)
    !Ecinetique_d = Ecin0

    pp(1:3,1:im) = pp(1:3,1:im)*rga_i(1:3,1:im) + gau(1:3,1:im)
    vp(1:3,1:im) = pp(1:3,1:im)!/m_i(1:3,1:im)
   ! do ic =1,3
   !    Ecin1 = Ecin1   + DOT_PRODUCT(pp(ic,1:im),vp(ic,1:im))/two
    !enddo
    !pp(1:3,1:im) =  pp(1:3,1:im) + (fp(1:3,1:im)+fadd(1:3,1:im))*tstep/two
    xp(1:3,1:im)=  xp(1:3,1:im) + pp(1:3,1:im)*dt!/m_i(1:3,1:im)
!    write(*,*) ' dxp ' ,pp(1:3,1:im),tstep/m_i(1:3,1:im)
    do ic=1,3
       xbar(ic)    = sum(xp(ic,1:im))/real(im)-xbarini(ic) ! déplacement du barycentre
       xp(ic,1:im) = xp(ic,1:im) - xbar(ic)  ! on recentre tout le systeme 
    enddo

    deltaU = potist + potistadd

    call calfoljc(xp,fp)
    !call calfoljq4(xp,fp)
    deltaU = -deltaU + potist + potistadd

    pp(1:3,1:im) = pp(1:3,1:im) + (fp(1:3,1:im)+fadd(1:3,1:im))*dt/two
    vp(1:3,1:im) = pp(1:3,1:im) !/ m_i(1:3,1:im)
   ! do ic=1,3
    !   Ecin3 = Ecin3 + DOT_PRODUCT(pp(ic,1:im),vp(ic,1:im))/two
    !enddo

    pp(1:3,1:im) = pp(1:3,1:im)*rga_i(1:3,1:im) + gau(4:6,1:im)
    vp(1:3,1:im) = pp(1:3,1:im)!/m_i(1:3,1:im)

   ! do ic=1,3
    !  Ecin4 = Ecin4 + DOT_PRODUCT(pp(ic,1:im),vp(ic,1:im))/two
    !enddo

    !erreur = ( deltaU + Ecin3 - Ecin1  )/(bk*text_teledyn)



   ! xprob = dexp(-erreur)
   ! call random_number(xalea)

!    write(*,*) 'erreur xprob xalea ', erreur ,xprob, xalea,bk*text_teledyn
!    write(*,*) ' deltaU ',deltaU,potistadd,potist
    !if ((xalea.lt.xprob).and.(ldistance)) then  
!       print *,'acceptation'
    !   ncomptinter = ncomptinter + 1
     !  betaq      = betaq     + erreur ! (Ecin4-Ecin3+Ecin1-Ecin0)/(bk*text_teledyn)
       !Ecinetique = Ecin4
    !else
     !  xp(1:3,1:im)=  xp_d(1:3,1:im)
      ! vp(1:3,1:im)=- vp_d(1:3,1:im)
      ! Ecinetique = Ecinetique_d
      ! call calfoljc(xp,fp)
!       print *,'refus'
      ! if (.not.ldistance) idistance = idistance +1
      ! ldistance = .true.
    !endif
    return

End subroutine langevin_ljc


Subroutine langevin(xp, vp, fp)!on s'en fiche pour l'instant
    !   M o d u l e s
    !-----------------------------------------------
    implicit none
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer ic,i

    real(double)  :: xp(3,imm)
    real(double)  :: vp(3,imm)
    real(double)  :: fp(3,imm)
    real(double)  :: pp(3,imm)

    real(double) :: xbar(3)

    real(double) Ecin0
    real(double) Ecin1
    real(double) Ecin3
    real(double) Ecin4

    real(double) Ebar0
    real(double) Ebar1
    real(double) Ebar3
    real(double) Ebar4

   !-----------------------------------------------
    integer   :: ic_local,iatom
   !-----------------------------------------------

!    write(6,*) ' vp = ',vp(1,1),vp(1,1)**2*cm(1)/(bk*text_teledyn),fp(1,1),tstep/two

    call genere_bruit

    Ecin0 = zero
    Ecin1 = zero
    Ecin3 = zero
    Ecin4 = zero

    Ebar0 = zero
    Ebar1 = zero
    Ebar3 = zero
    Ebar4 = zero

    pp(1:3,1:im1) = vp(1:3,1:im1)*m_i(1:3,1:im1)

!    write(*,*) ' potist  module =',potist,potistadd
!   write(*,*) 'fadd avant langevin', fadd(1:3,1),xpq4,xp(1,im1)
!   write(*,*) 'pp    pp(1:3,im1) ', pp(1:3,im1)
!    write(*,*) ' potistadd potistaddE', potistadd,potistaddE
!    deltaU =  potist/(bk*text_teledyn) + potistadd*barbeta +  potistaddE *barbeta 
    deltaU =  (potist + potistadd +  potistaddE)/(bk*text_teledyn)

    do ic =1,3
       Ecin0 = Ecin0   + DOT_PRODUCT(pp(ic,1:im),vp(ic,1:im))/two
    enddo

    Ebar0 =   barbeta*pp(1,im1)*vp(1,im1)/two + barbetaE*pp(2,im1)*vp(2,im1)/two 

    pp(1:3,1:im1) = pp(1:3,1:im1)*rga_i(1:3,1:im1) + gau(1:3,1:im1)
    vp(1:3,1:im1) = pp(1:3,1:im1)/m_i(1:3,1:im1)

    do ic =1,3
       Ecin1 = Ecin1   + DOT_PRODUCT(pp(ic,1:im),vp(ic,1:im))/two
    enddo

    Ebar1 =   barbeta*pp(1,im1)*vp(1,im1)/two +  barbetaE*pp(2,im1)*vp(2,im1)/two 
!   call calfo  ! car déjà appelé
    pp(1:3,1:im) =  pp(1:3,1:im) + (fp(1:3,1:im)+fadd(1:3,1:im))*tstep/two
    pp(1:2,im1)  =  pp(1:2,im1)  + fadd(1:2,im1)*tstep/two*alphadd(1:2)
    do ic=1,3
       xbar(ic)    = sum(xp(ic,1:im))/real(im) ! barycentre sur les particules
    enddo
    xp(1:3,1:im)=  xp(1:3,1:im) + pp(1:3,1:im)*tstep/m_i(1:3,1:im)

    xp(1,im1)=  xp(1,im1) + pp(1,im1)*tstep/m_i(1,im1)
    xp(2,im1)=  xp(2,im1) + pp(2,im1)*tstep/m_i(2,im1)

    do ic=1,3
       xbar(ic)    = sum(xp(ic,1:im))/real(im)-xbarini(ic) ! déplacement du barycentre
       xp(ic,1:im) = xp(ic,1:im) - xbar(ic)  ! on recentre tout le systeme 
    enddo
!    write(*,*) ' fp en 1  ',
    call calfoljc(xp,fp)
    call calfoljq4(xp,fp) ! a besoin de fp pour calculer fadd  

    pp(1:3,1:im) = pp(1:3,1:im) + (fp(1:3,1:im)+fadd(1:3,1:im))*tstep/two
    vp(1:3,1:im) = pp(1:3,1:im) / m_i(1:3,1:im)

    pp(1:2,im1) = pp(1:2,im1) +  fadd(1:2,im1)*tstep/two*alphadd(1:2)
    vp(1:2,im1) = pp(1:2,im1) / m_i(1:2,im1)

!    write(*,*) ' langevin  pp(1:2,im1) ',pp(1:2,im1)
!    write(*,*) ' langevin alphadd(1:2) ', fadd(1:2,im1)*alphadd(1:2)
!    write(*,*) '                       ',rga_i(1:3,im1)
!    stop

    do ic=1,3
       Ecin3 = Ecin3 + DOT_PRODUCT(pp(ic,1:im),vp(ic,1:im))/two
    enddo
    Ebar3 =   barbeta*pp(1,im1)*vp(1,im1)/two + barbetaE*pp(2,im1)*vp(2,im1)/two 

    pp(1:3,1:im1) = pp(1:3,1:im1)*rga_i(1:3,1:im1) + gau(4:6,1:im1)
    vp(1:3,1:im1) = pp(1:3,1:im1)/m_i(1:3,1:im1)

    do ic=1,3
      Ecin4 = Ecin4 + DOT_PRODUCT(pp(ic,1:im),vp(ic,1:im))/two
    enddo

    Ebar4 =  barbeta*pp(1,im1)*vp(1,im1)/two + barbetaE*pp(2,im1)*vp(2,im1)/two 

   ! betaq      = betaq     + (Ecin4-Ecin3+Ecin1-Ecin0)/(bk*text_teledyn)
    barbetaq   = barbetaq  + (Ebar4-Ebar3+Ebar1-Ebar0)
    Ecinetique = Ecin4     +  Ebar4*barbeta*(bk*text_teledyn)

!    write(*,*) ' xp en im1   ', xp(1,im1),xp(2,im1)
!    write(*,*) ' pp en im1   ', pp(1,im1),pp(2,im1)
!    write(*,*) ' pp^2 en im1 ', pp(1,im1)**2*barbeta,pp(2,im1)**2*barbeta

!    deltaU =  (Ecin3-Ecin1)/(bk*text_teledyn) + (Ebar3-Ebar1)*barbeta +   (potist )/(bk*text_teledyn) + potistadd*barbeta +  potistaddE *barbeta -deltaU    
    deltaU =  (Ecin3-Ecin1)/(bk*text_teledyn) + (Ebar3-Ebar1) +   (potist+potistadd+potistaddE )/(bk*text_teledyn) -deltaU

!    write(*,*) ' deltaU potistadd potistaddE work ',deltaU,potistadd,potistaddE,work
!    write(*,*) 'fadd après langevin', fadd(1:3,1),xpq4,xp(1,im1)
!    write(*,*) ' potistadd potistaddE', potistadd,potistaddE
    work = work + deltaU

!   write(6,*) 'barbetaq ',barbetaq, potistadd/(bk*text_teledyn),xp(1,im1),im1
!   write(6,*) 'Ebar ',Ebar0/(bk*text_teledyn),Ebar1/(bk*text_teledyn),Ebar3/(bk*text_teledyn),Ebar4/(bk*text_teledyn)

!      pp_1 = pp+pp
!      Ecini = pp*pp/2.0
!      pp          = pp + force*dtos2
!      xp          = xp   + pp*dto
!      call calfo()
!      pp          = pp + force*dtos2
!      Ecinf = pp*pp/2.0
!      fdto_bruit  = gau(2,itera)
!      pp_2        = pp
!      pp          = pp*rgammas2 + fdto_bruit
!      pp_2 = pp_2 + pp

    return

End subroutine langevin


Subroutine init_ljc(xp,vp,fp)
   implicit none
   integer NEMAX,NLSTEP,nmax,i,j,ic_local,iatom,ic,itera
   real(double) r6,Q4init,q4,q6,vtot,gamma,xbar(3)
   PARAMETER(NEMAX=100)
   PARAMETER(NLSTEP=10,nmaX=180)
   
   real(double) X0(N),Y0(N),Z0(N)
   real(double) PX0(N),PY0(N),PZ0(N)
   
   real(double) Xref(N),Yref(N),Zref(N)
   real(double) VXref(N),VYref(N),VZref(N)
   
   real(double) X(N),Y(N),Z(N)
   real(double) VX(N),VY(N),VZ(N)
   
   real(double) Xt(N),Yt(N),Zt(N)
   real(double) VXt(N),VYt(N),VZt(N)
   real(double) FX(N),FY(N),FZ(N)
   real(double) AX(N),AY(N),AZ(N)
   real(double) RIJ(N,N),RCIJ(N,N)
   real(double) H,PX(N),PY(N),PZ(N)
   real(double) XG,YG,ZG,ENTOTS,H2
   real(double) xtoto
!   real(double) 

   CHARACTER*25 filename
   real(double) xij(im,im)
   real(double) yij(im,im)
   real(double) zij(im,imm)
   real(double) dij(im,im)
   !common/pairdist/xij,yij,zij,dij
   real(double) q4coef(3,0:4)
   real(double) q6coef(4,0:6)

   common/q4coef/q4coef
   common/q6coef/q6coef

   character*2 type
   integer nsize

   real(double)  :: xp(3,imm)
   real(double)  :: vp(3,imm)
   real(double)  :: fp(3,imm)

  integer ifenetre,iselect,icourant,ienergy
  real(double) :: tab_lfe(-10:nfenetre+10,2)
  real(double) :: tab_cumul(-10:nfenetre+10,4)
  real(double) :: tab_energy(-10:nfenetre+10,2)
  real(double) :: sum_tot
  real(double) :: xalea,xexp_work,sumexp_work
  real(double) :: xtemp0,xtemp1,xtemp2,xtemp3,xtemp4,dt
  integer      :: ivois
! 
   call shinit
!   open(10,file='xyz.in')
!   i=0
!0001 i=i+1
!   read(10,*,end=0009) xp(1:3,i)
!   goto 0001
!0009 continue
   nsize=im
   print *,'Found ',nsize,' atoms.'
!   close(10)
   
   vtot=0.D0
   do i=1,nsize
     ! write(*,*) ' pos atomes ',xp(1:3,i)
      do j=i+1,nsize
         xij(i,j)=xp(1,j)-xp(1,i)
         yij(i,j)=xp(2,j)-xp(2,i)
         zij(i,j)=xp(3,j)-xp(3,i)
         dij(i,j)=sqrt(xij(i,j)**2+yij(i,j)**2+zij(i,j)**2)
         xij(j,i)=-xij(i,j)
         yij(j,i)=-yij(i,j)
         zij(j,i)=-zij(i,j)
         dij(j,i)= dij(i,j)
         
         r6=1./dij(i,j)**6
         
         vtot=vtot+r6*(r6-1.D0)
      enddo
   enddo
   
   vtot=vtot*4.
   print *,'Energy of configuration:',vtot
   
   Q4init=q4(nsize,xij,yij,zij,dij)
   
   print *,'Q4=',q4init
   Q4init=q6(nsize,xij,yij,zij,dij)
   print *,'Q6=',q4init
!   stop
   !       GRADIENT OF Q4
   
   call fq4(nsize,xij,yij,zij,dij,fq4x,fq4y,fq4z)

  ! print *,'Gradient:'

   do i=1,nsize
   !!!!!   print *,fq4x(i),fq4y(i),fq4z(i)
   enddo

   !       GRADIENT OF Q6
   
   call fq6(nsize,xij,yij,zij,dij,fq6x,fq6y,fq6z)

  !!!! print *,'Gradient:'

   do i=1,nsize
   !!!   print *,fq6x(i),fq6y(i),fq6z(i)
   enddo

   lq6 = .true. 
!   stop

   !write(*,*) 'centre de masse',sum(xp(1,1:nsize)),sum(xp(2,1:nsize)),sum(xp(3,1:nsize))

   sig_lj  = one
   rsig_lj = 2.25*sig_lj ! valeur dans la publie de florent      
   eps_4lj = two*two ! *cb

   im=nsize

   dlambdai = 0.2000000000
   dlambdaf = 0.0000000000
   x111   = 0.2
   xmaxq6 = 0.6
   energie_min = -174.0
   energie_max = -144.0
   

   call allocate_teledyn

   it_teledyn=0
!   tstep_teledyn=tstep  / utemps 
!   tstep = tstep/utemps  ! retour a des unites LJ 
!   tstep = tstep/10.0   ! loop3 et loop4
!   tstep = tstep/20.0  trop lent  
   write(*,*) ' tstep =  ',tstep 
   write(*,*) 'utemps =',utemps
!   stop
   im1=im+1
   sigma=1.d-2 
   write(6,*) ' kappa ',kappa 
   kappas2 = kappa/two

   ! initialisation des parametres du langevin

   write(*,*) 'text_teledyn ' ,text_teledyn

   text_teledyn = text_teledyn/bk

   write(*,*) 'text_teledyn ' ,text_teledyn
   gamma=10.0/(tstep)!*1d2)

   write(*,*) 'gamma*tstep/two ',gamma*tstep/two
!   stop
   m_i(1:3,1:im) = one
   rga_i(:,:)    = exp(-gamma*tstep/two)
   sig_i(:,:)    = sqrt(m_i(:,:)*text_teledyn*bk*(one-exp(-gamma*tstep)))

   write(6,*) ' rga_i  ', rga_i(:,1:1)
   write(6,*) ' sig_i ', sig_i(:,1:1)
   write(6,*) ' bk text_teledyn',bk,text_teledyn,one-exp(-gamma*tstep/two)
!   write(6,*) ' cm ',m_i(1:3,1:im)
!   stop

!   parametres additionnels

!  parametres sur l'énergie
   xp(1,im1)  = 0.02  !0.185

   call calfoljc(xp,fp)
   xp(1,im1)  = xpq4  !0.185

   print*,'xpq4 ',xpq4
   xp(2,im+1) = -163.0  ! xp(1,1)+xp(2,1)+xp(3,1)
   xp(3,im+1) = zero

   kappas2 = kappa/2.0
   kappaEs2 = kappaE/2.0
   barbeta        = one/(text_teledyn*bk)
   barbetaE       = one/(text_teledyn*bk)
   m_i(1:2,im1) = m_i(1,1)*real(im)*massadd(1:2)


   vp(1,im1) = -sqrt(one/barbeta/m_i(1,im1)) ! LJCOMB 5
   vp(2,im1) = sqrt(one/barbetaE/m_i(2,im1))! LJCOMB 5

   print *,m_i(1:3,im+1),m_i(1,1)

   bargamma=one/(1d6)

   write(*,*) ' bargamma ', bargamma
!   stop
   rga_i(1:2,im1) = exp(-bargamma/two)

   sig_i(1,im1)   = sqrt(m_i(1,im+1)*(one-exp(-bargamma))/barbeta)
   sig_i(2,im1)   = sqrt(m_i(2,im+1)*(one-exp(-bargamma))/barbetaE)

   rga_i(3,im1) = zero
   sig_i(3,im1)   = zero
   Ecinetique = zero
   betaq      = zero
   barbetaq   = zero

   m_tot = sum(m_i(1,1:im))

   fnampout = fnam(1:lenfnam)//'.pout'
   fnamwout = fnam(1:lenfnam)//'.wout'
   fnamfout = fnam(1:lenfnam)//'.fout'
   fnamsout = fnam(1:lenfnam)//'.sout'
   fnamvout = fnam(1:lenfnam)//'.vout'
   fnamhout = fnam(1:lenfnam)//'.hout'
   fnamfin  = fnam(1:lenfnam)//'.fin'
   fnamcout = fnam(1:lenfnam)//'.cout'
   fnamhin  = 'contour.hin'
   write(6,*) ' fnampout',fnampout

   lupout  = 17
   luwout  = 18
   lufout  = 19
   lusout  = 20
   luvout  = 21
   luhout  = 22
   
   lufin   = 23
   luhin   = 24
   

   open(unit=luvout, file=fnamvout, status='unknown')
   open(unit=lusout, file=fnamsout, status='unknown')
   open(unit=lupout, file=fnampout, status='unknown')
   open(unit=luwout, file=fnamwout, status='unknown')
   open(unit=lucout, file=fnamcout, status='unknown')

   if (nchemin_teledyn.eq.0) then 
      print*, 'Traitemement des données uniquement'
      call traitement
      stop
   endif

!   xtempmin = 3.61e-12
!   xtempmin = 1.0e-9
!   xtempmin = 1.0e-4
!   xtempmin = 5.0e-2 
!   xtempmin = 1.0e-10
!   xtempmin = one


   write(*,*) ' xtempmin = ', xtempmin

!   stop
   open(unit=luhin, file=fnamhin, status='old')

   do ifenetre=-10,nfenetre+10
      do ienergy=0,nfenetre
         read(luhin,'(e18.7,e18.7,e18.7,e18.7)') xtemp0,xtemp1,xtemp2,xtemp3
         pot_contour(ifenetre,ienergy) = log(max(xtemp2/xtemp3,xtempmin))
!         if ((xtemp2.gt.zero).and.(xtemp2.lt.xtempmin)) xtempmin = xtemp2
      enddo
      read(luhin,*) 
   enddo
!   print*, xtempmin
!   stop
   do ifenetre=-10,nfenetre+10
      do ienergy=0,nfenetre
         xtemp2 = pot_contour(ifenetre,ienergy)
         !            if (xtemp2.eq.zero) pot_contour(ifenetre,ienergy) = xtempmin
!         print*,pot_contour(ifenetre,ienergy) 
      enddo
   enddo




   betaq   =zero
   barbetaq=zero
   dlambda =zero
   
   do ic=1,3
      xbarini(ic)  = sum(xp(ic,1:im)*m_i(ic,1:im))/m_tot
   enddo
   
   write(6,*) ' xbarini =',xbarini
   
   ! stop
   xbar = xbarini


   if (.not.lta) then 
      xp(1,im+1)= dlambdai
   endif
   

! test
!   xpq4=0.04
!   potist = -158.0
!   call calpot_auxi
!   print * , contour_auxiliary,exp(contour_auxiliary)
!   stop
! fin test


   call calfoljc(xp,fp)
  !!! write(*,*) ' potist ',potist
!   write(*,*) 'fp',fp(1:3,1:im)
!      call fq4(nsize,xij,yij,zij,dij,fq4x,fq4y,fq4z)
     !!!!!!  write(*,*) 'q4 ',q4(nsize,xij,yij,zij,dij),xpq4,xp(1,im1)

   iteration=0
   do ic=1,100
      do itera=1, 10
         if (lta) then
            call langevin(xp,vp,fp)
         else
            call langevin_ljc(xp,vp,fp,dt)
         endif
      enddo
!         write(*,*) ' fadd fp en 1 ',fadd(1,1),fp(1,1),xpq4 
 !        write(*,*) ' potistadd  xp(1,im1) ' , potistadd,xp(1,im1),xpq4

 !     write(*,*)
      
!!!!!!!!!!!!!!     !write(*,*) ' equilibrage potist ',potist

      call calfoljq4(xp,fp)
!      call fq4(nsize,xij,yij,zij,dij,fq4x,fq4y,fq4z)
      !!!!!! write(*,*) 'q4 ',q4(nsize,xij,yij,zij,dij),xpq4,xp(1,im1)

   enddo
   
   potist0 = potist
!!!!!!!!!!!!!!!!!!  ! write(*,*) ' potistadd ',potistadd
!   stop
   xp_s(1:3,1:im+1)=xp(1:3,1:im+1)
   vp_s(1:3,1:im+1)=vp(1:3,1:im+1)
   potist0 = potist+potistadd
   Ecinetique0 = Ecinetique
   
!   nfrequence =  niteration/10
   if (nfrequence.eq.0) nfrequence=1
   dlambda=(dlambdaf-dlambdai)/real(niteration)

   rang=1   
   write(*,*) ' nchemin niteration kappa text_teledyn ',nchemin_teledyn , niteration , kappa , text_teledyn

!   write(*,*) 'fp',fp(1:3,1:im)
!   stop  'Normal end.'

   call calfoljq4(xp,fp)

 !!  write(*,*) 'fadd(1,1:im) ',fadd(1,1:im)
 !  write(*,*) 
 !  write(*,*) ' sum fadd ',sum(fadd(1,1:im)),sum(fadd(2,1:im)),sum(fadd(3,1:im))
 !  write(*,*) 
  ! write(*,*) ' fq4x ',fq4x(1:im),fq4y(1:im),fq4z(1:im)
 !  write(*,*) 
 !   write(*,*) ' sum fq4x ',sum(fq4x(1:im)),sum(fq4y(1:im)),sum(fq4z(1:im))
!   stop

  

 end subroutine init_ljc
 
Subroutine loop_ljc(xp, vp, fp)

   implicit none
   integer ic,ifenetre,iselect,icourant,ifen_q6
   
   real(double)  :: xp(3,imm)
   real(double)  :: vp(3,imm)
   real(double)  :: fp(3,imm)
   
   
   real(double) :: tab_lfe(-10:nfenetre+10,2)
   real(double) :: tab_cumul(-10:nfenetre+10,4)
   real(double) :: tab_energy(-10:nfenetre+10,2)
   real(double) :: tab_contour(-10:nfenetre+10,0:nfenetre)
   real(double) :: tab_contour_q6(-10:nfenetre+10,0:nfenetre)
   real(double) :: cumul_contour(-10:nfenetre+10,0:nfenetre)
   real(double) :: cumul_contour_q6(-10:nfenetre+10,0:nfenetre)
   real(double) :: sum_tot,tab_correction
   real(double) :: xalea,xexp_work,sumexp_work
   integer      :: ivois,ienergy
   real(double) :: xtemp0,xtemp1,xtemp2,xtemp3,xtemp4,dt
   real(double) :: pot_auxi_select,delta_auxi
   real(double) :: contour_auxi_select

   iselect = 0
   tab_cumul(-10:nfenetre+10,1:3)=zero
   cumul_contour(-10:nfenetre+10,0:nfenetre)=zero
   cumul_contour_q6(-10:nfenetre+10,0:nfenetre)=zero

   write(*,*) 'nchemin ',nchemin_teledyn
      write(*,*) 'lta = ' ,lta

   do iloop=1,nchemin_teledyn
      betaq         = zero
      barbetaq      = zero
      work          = zero
      delta_auxi    = zero

      call tirage_vitesse(vp_s)

      xp(1:3,1:im1)  = xp_s(1:3,1:im1)
      vp(1:3,1:im1)  = vp_s(1:3,1:im1)
      if (.not.lta) then 
      xp(1,im+1)    = dlambdai+(dlambdaf-dlambdai)*real(iselect)/real(niteration)
      print *, xp(1:2,im+1)
      stop
      endif
      call calfoljc(xp,fp)
      call calfoljq4(xp,fp)

      pot_auxi_select = pot_auxiliary
      contour_auxi_select = contour_auxiliary

      write(6,*) 'xp en im+1 ' , xp(1:3,im+1)
      write(*,*) 'potist0 ',potist0/(bk*text_teledyn),potistadd/(bk*text_teledyn)
      ncomptinter = 0
      sumexp_work = one
     ! write(*,*) ' sum fadd ',sum(fadd(1,1:im)),sum(fadd(2,1:im)),sum(fadd(3,1:im))
      IF (MOD(ITERATION,nfrequence).EQ.0)    write(luwout,'(e18.7,e18.7,e18.7,e18.7,e18.7,e18.7,e18.7)') work,xpq4,potist,xp(1:2,im+1),xpq6,delta_auxi

      rec_work(iselect)   = work
      rec_para(iselect)   = xpq4
      rec_q6(iselect)     = xpq6
      rec_energy(iselect,1) = potist!/(bk*text_teledyn)
      rec_correction(iselect) = dexp(contour_auxiliary)
      
      write(*,*) ' potist work ' ,potist,work,contour_auxiliary
      
      do iteration=iselect+1,niteration

!         call calfoljq4(xp)
!         deltawork = potistadd
!         xp(1,im1) =dlambdai+(dlambdaf-dlambdai)*real(iteration)/real(niteration)
!         call calfoljq4(xp)
!         deltawork = potistadd -deltawork
!         work       = work + deltawork/(bk*text_teledyn)
         if (iteration.lt.iselect+4) then
            write(*,*) ' iteration = ',iteration
            write(*,*) ' work, potist , potistadd,potistaddE ', work ,potist , potistadd,potistaddE
         endif

         if (lta) then 
            call langevin(xp,vp,fp) 
          write(*,*) ' xpq4 potist Edd',xpq4,potist,xp(1:2,im1)
         else
            call langevin_ljc (xp, vp, fp,dt)
         endif
         !        call control_output(xp)
         IF (MOD(ITERATION,nfrequence).EQ.0)    then 
            write(luwout,'(e18.7,e18.7,e18.7,e18.7,e18.7,e18.7,e18.7)') work,xpq4,potist,xp(1:2,im+1),xpq6,delta_auxi
!            write(*,*) ' sum fadd ',sum(fadd(1,1:im)),sum(fadd(2,1:im)),sum(fadd(3,1:im))
!            write(*,*) ' sum fp   ',sum(fp(1,1:im)),sum(fp(2,1:im)),sum(fp(3,1:im))
!            write(*,*) ' fp en 1  ',sum(fp(1,1:1)),sum(fp(2,1:1)),sum(fp(3,1:1))
         endif

!         delta_auxi = pot_auxiliary - pot_auxi_select
         delta_auxi = contour_auxiliary - contour_auxi_select

         rec_work(iteration)       = work + delta_auxi
         rec_para(iteration)       = xpq4
         rec_q6(iteration)       = xpq6
         rec_energy(iteration,1)   = potist!/(bk*text_teledyn)
         rec_correction(iteration) = exp(contour_auxiliary)
!         write(*,*) ' delta_auxi,contour_auxiliary ', delta_auxi,contour_auxiliary
!         write(*,*) 'rec_correction(iteration)' ,rec_correction(iteration),xpq4,potist
         xexp_work   = dexp(-(work+delta_auxi)) 
         sumexp_work = sumexp_work + xexp_work
         xexp_work   = xexp_work  / sumexp_work
         call random_number(xalea)
         
         if (xalea.lt.xexp_work) then  
            icourant = iteration
            xp_c = xp
            vp_c = vp
            write(*,*) ' xexp_work xalea ',  xexp_work,xalea
            print *,'icourant=',icourant
         endif
         
!!         IF (MOD(ITERATION,nfrequence).EQ.0) THEN
!!           print *,'icourant=',icourant
!!            print *,'iteration=',iteration
!!            print *,'work = ',work,xpq4,potist
!!            write(*,*) ' xexp_work xalea ',  xexp_work,xalea
!!         ENDIF
         !        stop
      enddo
      
      betaq    = zero
      barbetaq = zero
      work     = zero
      xp(1:3,1:im1)=xp_s(1:3,1:im1)
      vp(1:3,1:im1)=-vp_s(1:3,1:im1)  ! signe moins car on remonte le temps
      if (.not.lta) then
         xp(1,im+1) = dlambdai+(dlambdaf-dlambdai)*real(iselect)/real(niteration)
         print*, ' failure '
         stop
      endif
      call calfoljc(xp,fp)
      call calfoljq4(xp,fp)
      ! write(6,*) 'xp en im+1 ' , xp(1:3,im+1)
      write(*,*) 'potist0 ',potist0/(bk*text_teledyn),potistadd/(bk*text_teledyn)
      write(*,*) ' xp(1:2,im1) xpq4,potist ',xp(1:2,im1), xpq4,potist,iselect
!      if (mod(iloop,2).eq.1) then
!      if (iloop.gt.1) then
         write(lupout,*)
         write(luwout,*)
!      endif
      
      ! sumexp_work = one
      
      do iteration=iselect-1,0,-1
         
!         call calfoljq4(xp)
!         deltawork = potistadd
!         xp(1,im1) =dlambdai+(dlambdaf-dlambdai)*real(iteration)/real(niteration)
!         call calfoljq4(xp)
!         deltawork = potistadd -deltawork
!         work       = work + deltawork/(bk*text_teledyn)

         if (lta) then
            call langevin(xp, vp, fp)
         else
            call langevin_ljc(xp, vp, fp,dt)
         endif
         !        call control_output(xp)
         IF (MOD(ITERATION,nfrequence).EQ.0)    then 
            write(luwout,'(e18.7,e18.7,e18.7,e18.7,e18.7,e18.7,e18.7)') work,xpq4,potist,xp(1:2,im+1),xpq6,delta_auxi
!            write(*,*) ' sum fadd ',sum(fadd(1,1:im)),sum(fadd(2,1:im)),sum(fadd(3,1:im))
         endif
         delta_auxi = pot_auxiliary -pot_auxi_select
         delta_auxi = contour_auxiliary -contour_auxi_select
         rec_work(iteration) = work+delta_auxi
         rec_para(iteration) = xpq4
         rec_q6(iteration) = xpq6
         rec_energy(iteration,1) = potist!/(bk*text_teledyn)
         rec_correction(iteration) = dexp(contour_auxiliary)
         xexp_work   = dexp(-(work+delta_auxi))
         sumexp_work = sumexp_work + xexp_work
         xexp_work   = xexp_work  / sumexp_work
         call random_number(xalea)

         if (xalea.lt.xexp_work) then  
            icourant = iteration
            write(*,*) ' xexp_work xalea ',  xexp_work,xalea
            print *,'icourant=',icourant
         endif

!!           IF (MOD(ITERATION,nfrequence).EQ.0) THEN
!!              print *,'icourant  =',icourant
!!              print *,'iteration =',iteration
!!              print *,'work = ',work
!!              write(*,*) ' xexp_work xalea ',  xexp_work,xalea
!!           ENDIF

        if (xalea.lt.xexp_work) then
           icourant = iteration
           xp_c = xp
           vp_c = -vp  ! car on remonte le temps
           write(*,*) ' xexp_work xalea ',  xexp_work,xalea
           print *,'icourant=',icourant,'inferieur a ',iselect
        endif
     enddo

     if (iselect.gt.0) then
!        if (mod(iloop,10).eq.1) then 
           write(lupout,*)
           write(luwout,*)
!        endif
     endif
     
     xp_s = xp_c
     vp_s = vp_c
     
     sum_tot = sum(dexp(-rec_work(0:niteration)))
     tab_correction = sum(dexp(-rec_work(0:niteration))*rec_correction(0:niteration))
     tab_correction = tab_correction /sum_tot
     write(*,*) ' sum_tot sumexp_work tab_correction = ',sum_tot , sumexp_work , tab_correction
     write(*,*)
     write(*,*)
!     write(*,*) rec_work(0:niteration)
!     write(*,*) rec_correction(0:niteration)
!     write(*,*) ' rec_work iselect et +1 ',rec_work(iselect),rec_work(iselect+1)

!     iselect = niteration-icourant
     iselect = icourant

     if (iloop.ge.1) then 

     tab_lfe(-10:nfenetre+10,1:2)   =zero
     tab_energy(-10:nfenetre+10,1:2)=zero
     tab_contour(-10:nfenetre+10,0:nfenetre)=zero
     tab_contour_q6(-10:nfenetre+10,0:nfenetre)=zero
     
     do iteration=1,niteration
        ifenetre = nint((rec_para(iteration)/x111)*real(nfenetre))
        ienergy  = nint((rec_energy(iteration,1)-energie_min)/(energie_max-energie_min)*real(nfenetre))

        if ((ifenetre.le.nfenetre+10).and.(ifenetre.ge.-10)) then
           tab_lfe(ifenetre,1)=tab_lfe(ifenetre,1)+dexp(-rec_work(iteration))*rec_correction(iteration)
           tab_lfe(ifenetre,2)=tab_lfe(ifenetre,2)+one*rec_correction(iteration)
           tab_energy(ifenetre,1)=tab_energy(ifenetre,1)+dexp(-rec_work(iteration))*rec_energy(iteration,1)*rec_correction(iteration)
           if ((ienergy.le.nfenetre).and.(ifenetre.ge.0)) then
              tab_contour(ifenetre,ienergy) = tab_contour(ifenetre,ienergy) + dexp(-rec_work(iteration))*rec_correction(iteration)
              tab_contour_q6(ifenetre,ienergy) = tab_contour_q6(ifenetre,ienergy) +dexp(-rec_work(iteration))*rec_correction(iteration)

           endif
        endif

        ifen_q6 = nint((rec_q6(iteration)/xmaxq6)*real(nfenetre))
        if ((ifen_q6.le.nfenetre+10).and.(ifen_q6.ge.-10)) then
           if ((ienergy.le.nfenetre).and.(ifenetre.ge.0)) then
              tab_contour_q6(ifen_q6,ienergy) = tab_contour_q6(ifen_q6,ienergy) + dexp(-rec_work(iteration))*rec_correction(iteration)

           endif
        endif
        
     enddo
     
     open(unit=lufout, file=fnamfout, status='unknown')
     open(unit=luhout, file=fnamhout, status='unknown')
 
    do ifenetre=-10,nfenetre+10

        tab_lfe(ifenetre,1)   = tab_lfe(ifenetre,1)/sum_tot
        tab_lfe(ifenetre,2)   = tab_lfe(ifenetre,2)/sum_tot

        tab_energy(ifenetre,1)= tab_energy(ifenetre,1)  /sum_tot

        tab_cumul(ifenetre,1) = tab_cumul(ifenetre,1)+tab_lfe(ifenetre,1)
        tab_cumul(ifenetre,2) = tab_cumul(ifenetre,2)-log(tab_lfe(ifenetre,1))
        tab_cumul(ifenetre,3) = tab_cumul(ifenetre,3)+tab_energy(ifenetre,1)
        tab_cumul(ifenetre,4) = tab_cumul(ifenetre,4)+tab_correction
     enddo
     
     do ifenetre=-10,nfenetre+10
!        if (MOD(iloop,10).EQ.1) then
           write(lupout,'(e18.7,e18.7,e18.7,e18.7,e18.7)') real(ifenetre)*x111/100.0, tab_lfe(ifenetre,1),tab_lfe(ifenetre,2)
!        endif

write(lufout,'(e18.7,e18.7,e18.7,e18.7,e18.7)') real(ifenetre)*x111/100.0,-log(tab_cumul(ifenetre,1)/tab_cumul(ifenetre,4)), &
 tab_cumul(ifenetre,2)/tab_cumul(ifenetre,4),tab_cumul(ifenetre,3)/tab_cumul(ifenetre,1),tab_cumul(ifenetre,4)
     enddo

     do ifenetre=-10,nfenetre+10
        do ienergy=0,nfenetre
           cumul_contour(ifenetre,ienergy)=cumul_contour(ifenetre,ienergy)+tab_contour(ifenetre,ienergy)/sum_tot
           cumul_contour_q6(ifenetre,ienergy)=cumul_contour_q6(ifenetre,ienergy)+tab_contour_q6(ifenetre,ienergy)/sum_tot

           xtemp0 = real(ifenetre)*x111/100.0 ! pas sur l'axe des Q
           xtemp1 = real(ienergy)/real(nfenetre)*(energie_max-energie_min)+energie_min !pas sur l'axe de l'energie
           xtemp2 = cumul_contour(ifenetre,ienergy)!/tab_cumul(ifenetre,4) ! somme sur les nombres de clones
           xtemp3 = tab_cumul(ifenetre,4)
           xtemp4 = cumul_contour_q6(ifenetre,ienergy)! pour q6

           write(luhout,'(e18.7,e18.7,e18.7,e18.7,e18.7)') xtemp0,xtemp1,xtemp2,xtemp3,xtemp4
!real(ifenetre)*x111/100.0 ,real(ienergy)/real(nfenetre)*(174.0-144.0)-174.0,cumul_contour(ifenetre,ienergy)/real(iloop)
        enddo
         write(luhout,*) 
     enddo

     endif

     xp(1:3,1:im) = xp_s(1:3,1:im)
     vp(1:3,1:im) = vp_s(1:3,1:im)
     call calfoljc(xp,fp)
     call calfoljq4(xp,fp)

     write(lusout,'(a29,i7,e18.7,e18.7,e18.7,e18.7)') 'select q4 energie sum q6: ',icourant,xpq4,potist,sum_tot , xpq6
!     write(lusout,*) ' sum_tot tab_correction :',sum_tot , tab_correction
     if (.not.lta) then
     write(lusout,*) ' taux d''acceptation    :',ncomptinter/real(ITERATION)
     write(lusout,*) ' taux de dépassement    :',idistance/real(ITERATION)
!     write(lusout,*) ' sum_tot tab_correction :',sum_tot , tab_correction
     endif
     close(lufout)
     close(luhout)

  enddo

end subroutine loop_ljc

Subroutine tirage_vitesse(vp)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
   implicit none
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    !-----------------------------------------------
   real(double) :: vp(3,imm)
   real(double) :: u1,u2,b1
   integer      :: ic
    !-----------------------------------------------
         call random_number(u1)
         call random_number(u2)
         b1=sqrt(-2.*log(u1))*cos(2.*pi*u2)
!         vp(ic,im1) = b1*sqrt(bk*text_teledyn/m_i(ic,im1))
         vp(1,im1) = b1*sqrt(one/m_i(1,im1)/barbeta)

         call random_number(u1)
         call random_number(u2)
         b1=sqrt(-2.*log(u1))*cos(2.*pi*u2)
         vp(2,im1) = b1*sqrt(one/m_i(2,im1)/barbetaE)

end subroutine tirage_vitesse

Subroutine traitement
  implicit none
  real(double) xtemp0,xtemp1,xtemp2,xtemp3
  real(double) t_energie(0:nfenetre),t_q4(-10:nfenetre+10)
  real(double) beta0,beta1
  real(double) prob_contour(-10:nfenetre+10,0:nfenetre)
  real(double) prob_q4(-10:nfenetre+10)
  integer ienergy,ifenetre,ite
  real(double) partfunction0,partfunction1
  real(double) tempmin,tempmax

  open(unit=luhout, file=fnamhout, status='old')
  
  do ifenetre=-10,nfenetre+10
     do ienergy=0,nfenetre
        read(luhout,'(e18.7,e18.7,e18.7,e18.7)') xtemp0,xtemp1,xtemp2,xtemp3
         prob_contour(ifenetre,ienergy) = xtemp2/xtemp3
         t_energie(ienergy) = xtemp1-energie_min
      enddo
      t_q4(ifenetre) = xtemp0
      read(luhout,*)
   enddo

  ! changement de temperature 
   print*,' text_teledyn =',text_teledyn,energie_min

   beta0=one/(text_teledyn*bk)
   beta1=one/0.05
   print*,' beta0 beta1 ',beta0,beta1

   ! calcul de l'énergie libre à beta0 et beta1

   partfunction0 = zero
   do ifenetre=-10,nfenetre+10
      do ienergy=0,nfenetre
         partfunction0 = partfunction0 + prob_contour(ifenetre,ienergy)*exp(-zero*t_energie(ienergy))
         ! print*,'partfunction0 ',partfunction0,prob_contour(ifenetre,ienergy),exp(-beta0*t_energie(ienergy)),beta0*t_energie(ienergy)
      enddo
   enddo

   partfunction1 = zero
   do ifenetre=-10,nfenetre+10
      do ienergy=0,nfenetre
         partfunction1 = partfunction1 + prob_contour(ifenetre,ienergy)*exp(-(beta1-beta0)*t_energie(ienergy))
      enddo
   enddo

   print*,'partfunction0 partfunction1',partfunction0,partfunction1

   open(unit=29, file='contour.05')
   
   do ifenetre=-10,nfenetre+10
      do ienergy=0,nfenetre
!         write(29,'(d19.7,d19.7,d19.7,d19.7)') 
         write(29,*) t_q4(ifenetre),t_energie(ienergy)+energie_min, prob_contour(ifenetre,ienergy)*exp(-(beta1-beta0)*t_energie(ienergy))*partfunction0/partfunction1
      enddo
      write(29,*) 
   enddo

   tempmax = one/beta0
   tempmin = 0.02
   open(unit=30, file='landscape')
   do ite=0,nfenetre
      beta1=one/(tempmax + (tempmin-tempmax)*real(ite)/real(nfenetre))
      prob_q4=zero
      partfunction1 = zero
      do ifenetre=-10,nfenetre+10
         do ienergy=0,nfenetre
            partfunction1     = partfunction1 + prob_contour(ifenetre,ienergy)*exp(-(beta1-beta0)*t_energie(ienergy))
            prob_q4(ifenetre) = prob_q4(ifenetre) + prob_contour(ifenetre,ienergy)*exp(-(beta1-beta0)*t_energie(ienergy))
         enddo
      enddo
      do ifenetre=-10,nfenetre+10
         prob_q4(ifenetre) = prob_q4(ifenetre)/partfunction1 
         write(30,*) t_q4(ifenetre), one/beta1 , -log(prob_q4(ifenetre))/beta1
      enddo
      write(30,*)
   enddo
  
   stop

end subroutine traitement


 Subroutine neufdoublewell (xp)

implicit none

real (double), dimension (:), allocatable :: Weight
integer, dimension (:),allocatable ::ToCopy
integer, dimension (:),allocatable ::ToDelete
integer, dimension (:),allocatable ::Tofollow
integer, dimension (:),allocatable ::Nb
integer, dimension (:),allocatable ::Number
!integer, parameter ::N=38 ! numero atomes du cluster
 
real (double) ::r     !distance interatomique
real (double) ::eta    !paramètre d'énergie du potentiel
real (double) ::sigma ! dist. éq. pour le potentiel LJ 
real (double),dimension(0:N-1) :: fxx
real (double),dimension(0:N-1) :: fyy
real(double),dimension(0:N-1) :: fzz
real(double),dimension(0:N-1) :: vsecx
real (double),dimension(0:N-1) :: vsecy
real (double),dimension(0:N-1) :: vsecz

real(double),dimension (N,N)::   xij
real (double),dimension (N,N)::   yij
real(double),dimension (N,N)::   zij
real (double),dimension (N,N)::   dij

real (double),dimension (:),allocatable::w4,w6,wp4,wp6, xpq4,xpq6,tv,tvmoy
integer, dimension (:),allocatable ::Nvite
integer, dimension(:), allocatable ::multip

real (double),dimension (N):: fq4x,fq4y,fq4z,fq6x, fq6y,fq6z

real (double), dimension (3,N):: xp
real (double), dimension (3,N):: vp
real(double),dimension (3,N):: fp
real (double),dimension (3,N)::fppro
real (double), dimension (3,N):: xppro
real (double), dimension(1:N):: delt
real (double), dimension(1:N):: gQ4x
real (double), dimension(1:N):: gQ4y
real (double), dimension(1:N):: gq4z
real (double), dimension(1:N):: gQ6x
real (double), dimension(1:N):: gQ6y
real (double), dimension(1:N):: gq6z
real (double), dimension(1:N):: gEx
real (double), dimension(1:N):: gEy
real (double), dimension(1:N):: gEz

real (double) :: q4,q6, q4moy, q6moy,w4moy, w6moy, lambdaq4,lambdaE, lambdaq6, logt, logt2! linear displacement
type Particule
    real (double), dimension (0:N-1):: qx  
    real (double), dimension (0:N-1):: qy 
    real (double), dimension (0:N-1):: qz  ! vector of position
    
    real (double), dimension (0:N-1):: px  
    real (double), dimension (0:N-1):: py 
    real (double), dimension (0:N-1):: pz  ! vector of impulsion
    
    real (double), dimension (0:N-1):: uqx 
    real (double), dimension (0:N-1):: uqy
    real (double), dimension (0:N-1):: uqz ! linear displacement
    
    real (double), dimension (0:N-1):: upx 
    real (double), dimension (0:N-1):: upy 
    real (double), dimension (0:N-1):: upz! linear displacement

    !RealNumber Lyap
 endtype

type (Particule), dimension (:), allocatable :: P

real (double) ::seed, utot, tot, x, U,x4,x6,tconf, tconftot,tconfmoy,kinetotdt, upro, q4ref
 !     Declaration of the variables
real(double)::fftot, fftotdt,laputot,laputotdt,ff,lapu, pi2, uproject
real (double) ::inter,AverageTime,NbTimeStep,dt,TotalTime,Time,w,Norm,TimeStore,c1,c2,sqrtgTh,eta1,eta2,gamma,eps, Temperature ! Lyapunov=0,LyapunovAv=0,mu=0,
integer :: NbClones,NMax,j,Ndelete,Ncopy,l,k,ran,m, a,i,b, iw4,iw6,lo,cont,ien,copiati,nfollow,nq4,teq,n0, n4

character(len=128) :: sortie
character(len=128) :: histotot
character(len=128) :: histopart
character(len=128) ::posfinal
character(len=128) ::snap
character(len=128) ::stock
character(len=128) ::stq4t
character(len=128) ::stq6t
character(len=128) ::q4q6

real (double)::energy,rga,beq,rien,q,kine,kinetot,pif, trasm1,trasm2, normgradq4, normgradq6, q6ref
real (double) ::xtemp0,xtemp1,xtemp2,xtemp3,xtemp4,xtemp5
integer, parameter::nfenetre=100
integer :: ienergy, ifenetre, ifen_q6, iter,kaq4,kaq6,niter,SYS

real (double):: stat(0:nfenetre, 0:nfenetre)
real (double):: statneg(0:nfenetre, 0:nfenetre)
real (double):: statpos(0:nfenetre, 0:nfenetre)
real (double):: stat6(0:nfenetre, 0:nfenetre)
real(double)::tvmoytot,xn,enprmoy,jtrans, jtot,trasm,trastot

real(double) :: tab_contour(-10:nfenetre+10,0:nfenetre)
real(double) :: tab_contour_q6(-10:nfenetre+10,0:nfenetre)
real(double) :: tab_cont_q4q6(-10:nfenetre+10,-10:nfenetre+10)
real(double) :: cumul_contour(-10:nfenetre+10,0:nfenetre)
real(double) :: cumul_contour_q6(-10:nfenetre+10,0:nfenetre)
real(double) :: cumul_contour_q4q6(-10:nfenetre+10,-10:nfenetre+10)
real(double) :: tab_energy(-10:nfenetre+10,2),toto

real (double)::x111, xmaxq6,energie_min,energie_max,ecinetique,tequilib,delta,beqtot,dhtot,beqtotdt,dhtotdt,intbeq, intdh,e
integer::icfen,icener,icfen_q6,a1, a2, a3, a4, a5

 real(double), dimension(:), allocatable::betaq,dh,ener,enerpro

write(6,*)'usage %s:\n'

!    Initialisation of the variables
write(6,*)'TotalTime '
read (*,*) TotalTime      ! Total duration of the simulation 
write(*,*) TotalTime 
write(6,*)'dt' 
read (*,*)  dt          ! Time step
write(*,*) dt
write(6,*)'NbClones'  
read (*,*) NbClones  
write(*,*) NbClones    ! Number of clones
write(6,*)'temperature'  
read (*,*) temperature  ! temperature (kT)
write (*,*) temperature
write(6,*)'gamma*dt'  
read (*,*) gamma  ! friction*dt
write (*,*) gamma
write(6,*)'inter'  
read (*,*) inter  
write(*,*) inter     ! interval between 2 data records
write(6,*)'sortie'  
read (*,*) sortie   
write (*,*) sortie! Name of the file where data are stored
write(6,*)'interval-annealing'  
read (*,*) inter 
write (*,*) inter ! annealing
write(6,*)'lo'
read (*,*) lo  
write (*,*) lo        ! cloner (lo=1) ou pas (lo=0)
write(6,*)'cont'
read (*,*) cont
write (*,*) cont
write(6,*)'niter'
read (*,*) niter
write (*,*) niter
write(6,*)'t_equilib'
read (*,*) teq
write (*,*) teq


open (unit=11, FILE=sortie, action='write', status='replace')
open (unit=31, action='write', status='replace')

histotot  =sortie(1:lenfnam)//'.thout'
histopart =sortie(1:lenfnam)// '.phout'
posfinal = sortie(1:lenfnam)//'.cin'
snap = sortie(1:lenfnam)//'.film'
stock = sortie(1:lenfnam)//'.stock'
stq4t = sortie(1:lenfnam)//'.stq4t'
stq6t= sortie(1:lenfnam)//'.stq6t'
q4q6= sortie(1:lenfnam)//'.qhout'

open(unit=14, file=histopart, action='write', status='replace')
open(unit=26, file=histotot,  action='write', status='replace')
open(unit=28, file=snap,  action='write', status='replace')
open(unit=29, file=stock,  action='write', status='replace')
open(unit=33, file=stq4t,  action='write', status='replace')
open(unit=34, file=stq6t,  action='write', status='replace')
open(unit=31, file=q4q6,  action='write', status='replace')

NMax= 100*NbClones      ! Max number of clones
!TimeStore = TotalTime-AverageTime ! Next time at which data will be stored 
                                           
! Array where the weight of every clones is stored
allocate (Weight(0:nmax-1)) 
allocate (w4(0:nmax))
allocate (w6(0:nmax))
allocate (wp4(0:nmax))
allocate (wp6(0:nmax))

allocate (tv(0:nmax-1))
allocate (tvmoy(0:nmax-1))
allocate (nvite(0:nmax-1))

allocate (xpq4(0:Nmax))
allocate (xpq6(0:Nmax))
allocate(multip(0:nmax))

allocate (betaq(0:nmax-1))
allocate (dh(0:nmax-1))
! Arrays where the clone to copy or to kill are stored
allocate(ToCopy(0:NMax-1))
allocate(ToDelete(0:NMax-1))
allocate(ener(0:nmax-1))
allocate(Tofollow(0:NMax-1))
allocate(enerpro(0:nmax-1))
 
! Array which contains the numbers from 1 to Nmax
allocate(Nb(0:NMax-1)) 

ToCopy(0:NMax-1)=0
ToDelete(0:NMax-1)=0
ener(:)=0
sigma=1
eta=1
enerpro(:)=0
xpq4(:)=0
xpq6(:)=0
gq4x(:)=0
gq4y(:)=0
gq4z(:)=0
gq6x(:)=0
gq6y(:)=0
gq6z(:)=0
gEx(:)=0
gEy(:)=0
gEz(:)=0
n0=0

SYS=1
gamma=gamma/(dt)

write(*,*)'gamma' , gamma, gamma*dt
!temperature=0.19
stat(:,:)=0
stat6(:,:)=0
!statneg(:,:)=0
!statpos(:,:)=0
kine=0
laputot=0
fftot=0
laputotdt=0
fftotdt=0
kinetotdt=0
kinetot=0
nfollow=0
enprmoy=0
ienergy=0
nq4=0
dhtotdt=0
e=0.000001
beqtotdt=0
x111   = 0.2
xmaxq6 = 0.6
energie_min = -174.0
energie_max = -144.0


     tab_contour(-10:nfenetre+10,0:nfenetre)=zero
     tab_contour_q6(-10:nfenetre+10,0:nfenetre)=zero
     tab_cont_q4q6(-10:nfenetre+10,-10:nfenetre+10)=zero
     cumul_contour_q4q6(-10:nfenetre+10,-10:nfenetre+10)=0
     cumul_contour(-10:nfenetre+10,0:nfenetre)=0
     cumul_contour_q6(-10:nfenetre+10,0:nfenetre)=0

 do j=0,NMax-1
    Nb(j)=j   ! Number(j) is the number of descendants of clone j at the cloning step
 enddo

allocate(Number(0:NbClones-1)) 

! Some constants used in the quasisymplectic stochastic integrator (cfr Mannella)
! c1 = 1 - gamma * dt / 4
! c2 =  1 / ( 1 + gamma * dt / 4 )
! sqrtgTh = sqrt(2*gamma*temperature*dt)

beq=0
!ecinetique=0

!!!!!!!!!!!!T EQUILIBRAGE
tequilib=(dt)*teq
write(*,*) 'tquilib',tequilib
!!!!!!!!!!!! T EQUILIBRAGE
  
write(*,*) 'sortie affichée'


allocate(P(0:NMax-1)) ! cas  sans lyapunov
!close(27)
!open(unit=27, file=fnamcin, status='old')

n4=4*nbclones
!!!EQULIBRAGE
xppro(:,:)=0
Time=0
iter=0

rien=0
Q=0
logt=0
logt2=0
pif=0
pi2=0
trastot=1
uproject=0
tv(:)=0
tvmoy(:)=0
nvite(:)=0

if ((gamma*dt).le.100000) then
 
   rga = exp(-gamma*dt/two)

else 

  rga=0

endif


write(*,*) 'rga',gamma, rga,two

write (*,*) 'initialisation' 

do  m = 0,NbClones-1 
     
       do j=0,N-1
       P(m)%uqx(j)=genrand()!-0.5
       P(m)%uqy(j)=genrand()!-0.5
       P(m)%uqz(j)=genrand()!-0.5

       P(m)%upx(j)=0+genrand()!-0.5
       P(m)%upy(j)=0+genrand()!-0.5
       P(m)%upz(j)=0+genrand()!-0.5
       enddo
 
       call Normalizza(P(m)%uqx,P(m)%uqy,P(m)%uqz,P(m)%upx,P(m)%upy,P(m)%upz, N)

       call CalculNorme(P(m)%uqx,P(m)%uqy, P(m)%uqz,P(m)%upx, P(m)%upy,P(m)%upz,N,Norm)

enddo



if(cont.eq.1) then

open(unit=27, file=posfinal, status='old')

do m=0,NbClones-1
   read (27,*) P(m)%qx(:)
   read (27,*) P(m)%qy(:)
   read (27,*) P(m)%qz(:)
   read (27,*) 
   read (27,*) P(m)%px(:)
   read (27,*) P(m)%py(:)
   read (27,*) P(m)%pz(:)      
   read (27,*)
!enddo 


!do m=0,NbClones-1
  read (27,*) P(m)%uqx(:)
  read (27,*) P(m)%uqy(:)
  read (27,*) P(m)%uqz(:)
  read (27,*) 
  read (27,*) P(m)%upx(:)
  read (27,*) P(m)%upy(:)
  read (27,*) P(m)%upz(:)      
  read (27,*)
enddo

close(27)

!do l=0,NbClones-1   
!do while ((Time.lt.Tequilib).and.(time.lt.totaltime))
!   do j=0,NbClones-1
              
!  call propagateA(P(j)%qx,P(j)%qy,P(j)%qz,P(j)%px,P(j)%py,P(j)%pz,P(j)%uqx,P(j)%uqy,P(j)%uqz,P(j)%upx,P(j)%upy,P(j)%upz,&
 !                 rga,temperature,dt,N,betaq(j),dH(j),ecinetique,ff,lapu)
!enddo
!enddo

!    do i=1,N-1
!      do j=i+1,N
!        xij(i,j)=(P(l)%qx(i-1)-P(l)%qx(j-1))
!        yij(i,j)=(P(l)%qy(i-1)-P(l)%qy(j-1))
!        zij(i,j)=(P(l)%qz(i-1)-P(l)%qz(j-1))
!        dij(i,j)=sqrt(xij(i,j)**2+yij(i,j)**2+zij(i,j)**2)
!        xij(j,i)=-xij(i,j)
!        yij(j,i)=-yij(i,j)
!        zij(j,i)=-zij(i,j)
!        dij(j,i)= dij(i,j)
!      enddo
!    enddo

!    xpq4(L)=q4(N,xij,yij,zij,dij)

!call Force(P(l)%qx,P(l)%qy,P(l)%qz,fp,Utot)

!write(*,*)'energie',Time,xpq4(L),Utot

!enddo

!endif

else  if (cont.eq.0) then

do while ((Time.lt.Tequilib).and.(time.lt.totaltime))

!write (*,*)'atomo', xp(1,1),vp(1,1)
lapu=0
ff=0

  call calfoljc2(xp, fp,utot)

  call langevin2(xp,vp,fp,dt,rga, temperature,Q,rien,ecinetique)

ff=sum(fp(1,:)**2+fp(2,:)**2+fp(3,:)**2)

!write (*,*) 'ff',utot

fftot=fftot+ff

call Hessien(xp(1,:), xp(2,:),xp(3,:), P(0)%uqx, P(0)%uqy, P(0)%uqz, Vsecx, vsecy,vsecz,lapU)

!write (*,*) 'lapu', time, lapu

laputot=laputot+lapu

do i=1,N-1
      do j=i+1,N
        xij(i,j)=xp(1,i)-xp(1,j)
        yij(i,j)=xp(2,i)-xp(2,j)
        zij(i,j)=xp(3,i)-xp(3,j)
        dij(i,j)=sqrt(xij(i,j)**2+yij(i,j)**2+zij(i,j)**2)
        xij(j,i)=-xij(i,j)
        yij(j,i)=-yij(i,j)
        zij(j,i)=-zij(i,j)
        dij(j,i)=dij(i,j)
      enddo
enddo

x4=q4(N,xij,yij,zij,dij)
x6=q6(N,xij,yij,zij,dij)



!if ((x6.le.0.08).and.(utot.lt.-164.0)) then
!time=tequilib
!endif
!write(*,*) 'energie' ,time, x4, x6

!write (*,*) 'dissipation', rien, q
!if (time.gt.99) temperature=0.05


time=time+dt

enddo


write(*,*)'tconf', time, fftot/laputot

write(*,*)'tcin', time, ecinetique/(1.5d0*real(N))


!write(*,*) 'ctrl', xp(1,:)
           
     do  m = 0,NbClones-1  
       P(m)%qx=xp(1,:)
       P(m)%qy=xp(2,:)
       P(m)%qz=xp(3,:)
     enddo

   do  m = 0,NbClones-1
       P(m)%px=vp(1,1:N)
       P(m)%py=vp(2,1:N)
       P(m)%pz=vp(3,1:N)
    enddo



endif
 
  
  do m = NbClones,NMax-1
    do j=0,N-1
    P(m)%qx(j) =0
    P(m)%px(j) =0
    P(m)%uqx(j)=0
    P(m)%upx(j)=0
 
    P(m)%qy(j) =0
    P(m)%py(j) =0
    P(m)%uqy(j)=0
    P(m)%upy(j)=0

    P(m)%qz(j) =0
    P(m)%pz(j) =0
    P(m)%uqz(j)=0
    P(m)%upz(j)=0
   enddo
   ! P(m).Lyap = 0
  enddo


  
write(*,*) 'initialisation faite' 


fftotdt=fftot
laputotdt=laputot
betaq(:)=q
dh(:)=0
ecinetique=0
kinetot=ecinetique/(1.5d0*real(N))
ecinetique=0
fftot=0
laputot=0
e=1.e-6

!!!!!!!!trsm2
trasm2=1

temperature=0.05

do while (Time.lt.TotalTime)

    Ncopy=0
    Ndelete=0
    copiati=0
    n0=0
    w=0
    u=0
    q4moy=0
    q6moy=0
    tconf=0
    beqtot=0
    dhtot=0
    kine=0
    jtrans=0
    jtot=0
    trasm1=0
    !trasm2=0   

    iter=int(time/dt)

!!!!! PARTICLES PROPAGATION

 do j=0,NbClones-1
              
  call propagateA(P(j)%qx,P(j)%qy,P(j)%qz,P(j)%px,P(j)%py,P(j)%pz,P(j)%uqx,P(j)%uqy,P(j)%uqz,P(j)%upx,P(j)%upy,P(j)%upz,&
                  rga,temperature,dt,N,betaq(j),dH(j),ecinetique,ff,lapu)

    !  call hessien()

     kine=kine+ecinetique

   ! write(*,*) 'cin',kine, lapu

     fftot=fftot+ff
     laputot=laputot+lapu
     dhtot=dhtot+dh(j)/nbclones
     beqtot=beqtot+betaq(j)/nbclones
     
    !tconftot=tconftot+tconf   
    ! write (*,*) 'dissipation', dH(j), betaq(j)

     !call propagateC(P(j)%qx,P(j)%qy,P(j)%qz,P(j)%uqx,P(j)%uqy,P(j)%uqz,temperature,dt,N) !!!! senza inerzia
 enddo

fftotdt=fftotdt+fftot/real(Nbclones)
laputotdt=laputotdt+laputot/real(nbclones)

!write(*,*) 'cin', time, (kine/real(nbclones))/(1.5d0*real(N))


!write(*,*) (kine/real(nbclones))/((3/2)*N)
kine=(kine/real(nbclones))/(1.5d0*real(N))


kinetot=kinetot+kine
!write(*,*) 'ctrl', time,  kine!, kinetot



beqtotdt=beqtotdt+beqtot
dhtotdt=dhtotdt+dhtot
intbeq=(beqtotdt)/((time-tequilib)/dt+1.0)
intdh=(dhtotdt)/((time-tequilib)/dt+1.0)

   ! write (*,*) 'dissipation', intdh,intbeq

!write(*,*) 'ctrl', kine, kinetot


 
tconfmoy=fftotdt/laputotdt

!if (tconfmoy.ge.0.19) then
!   write (*,*) 'equilibre', time
!  stop
!endif
!if(mod(iter,10).eq.0) then
!write(*,*) 'tcin1',kine
!endif

!
!if(time.eq.tequilib)then
!kinetotdt=kine
!else if (time.gt.tequilib) then 
kinetotdt=(kinetot)/((time-tequilib)/dt+1.0)
!endif


if(mod(iter,5000).eq.0) then
write(*,'(a, e11.4,e11.4, e11.4)') 'tcin',time, kinetotdt,tconfmoy
endif

open(unit=lucout, file=fnamcout, status='unknown')
     
!!!!CALCUL Q4,Q6,E DELLA SERIE DI CLONI DOPO PROPAGAZIONE ED ISTOGRAMMA

do l=0,NbClones-1   

    do i=1,N-1
      do j=i+1,N
        xij(i,j)=(P(l)%qx(i-1)-P(l)%qx(j-1))
        yij(i,j)=(P(l)%qy(i-1)-P(l)%qy(j-1))
        zij(i,j)=(P(l)%qz(i-1)-P(l)%qz(j-1))
        dij(i,j)=sqrt(xij(i,j)**2+yij(i,j)**2+zij(i,j)**2)
        xij(j,i)=-xij(i,j)
        yij(j,i)=-yij(i,j)
        zij(j,i)=-zij(i,j)
        dij(j,i)= dij(i,j)
      enddo
    enddo

    xpq4(L)=q4(N,xij,yij,zij,dij) 
!    call fq4(N,xij,yij,zij,dij,fq4x,fq4y,fq4z)
  
    xpq6(L)=q6(N,xij,yij,zij,dij)
!   call fq6(N,xij,yij,zij,dij,fq6x,fq6y,fq6z)

    call Force(P(l)%qx,P(l)%qy,P(l)%qz,fp,Utot)

    ener(l)=utot

!    energy=Utot  

    q4moy=q4moy+xpq4(l)
    q6moy=q6moy+xpq6(l)
    u=u+utot


!     w4(l)=0
!     w6(l)=0

 enddo

!if(mod(iter,100).eq.0) then
!write(*,'(a, e18.7, e18.7, e18.7, e18.7)')'param', time, q4moy ,q6moy, U!, logt, pif
!write (*,'(a, e11.4, e11.4, e11.4, e12.4, e12.4, e11.4)') 'param', time, q4moy ,q6moy,U, logt, pif!, logt!,tvmoytot
!endif


!!!!!!!!!!!!!!MD!!!!!!!!!!!!
!if ((time.gt.2000).or.(q4moy.ge.(0.07))) then
!if(mod(iter,niter).eq.0) then
!lo=1
!else
!lo=0
!endif
!!!!!!!!!!!!MD!!!!!!!!!!!!!!!!!!!!!


!if (time.gt.500) then

!write(*,*) 'inizio clonaggio', time

!lo=1

!else

!lo=0
!endif

!!!!!! STATISTIQUE DE CLONAGE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (lo.eq.1) then
 !write(*,*) 'sto clonando'
!!!!!! NORME DU CLONE PROPAGÉ POUR EN TIRER LE POIDS DE PROBABILITÉ      

do j=0,NbClones-1

   ! New norm of the tangent vector
   call CalculNorme(P(j)%uqx,P(j)%uqy,P(j)%uqz,P(j)%upx,P(j)%upy,P(j)%upz,N,norm)
  ! write (*,*) 'renorm', Norm  
  
     ! Lyap computation
     ! P(j).Lyap  += log(Norm)
      
   ! Compute the weight for the cloning step
   Weight(j) = Norm

   w=w+ Weight(j)
 
enddo

!jtot=w  



w=w/real(NbClones)
!write (*,*) 'renorm2', w, jtot
   !write (*,*) 'renorm2', w
!xn=genrand()
   
do j=0,NbClones-1

      Weight(j)=Weight(j)/w

!!!!!!!!!!!!!!!!! SYSTEMATIC SAMPLING!!!!!!!!!!!!!!!

!multip(j)=(floor(nbclones*sum(weight(0:j))+xn)-floor(nbclones*sum(weight(0:j-1))+xn))!/real(nbclones)

!if (multip(j).lt.1) then

!	     ToDelete(Ndelete)=j 
!            Ndelete=Ndelete+one
!           ien  = int((ener(j)-energie_min)/(energie_max-energie_min)*real(nfenetre))
!           kaq4 = int((xpq4(j)/x111)*nfenetre)
!           kaq6 = int((xpq6(j)/xmaxq6)*nfenetre)
!            stat(kaq4,ienergy)=stat(kaq4,ienergy)-1
!           stat6(kaq6,ien)=stat6(kaq6,ien)-1
!          if(Ndelete.gt.NbClones) then
!	         write(6,*)'All the clones have been killed'
!!	         stop
!            endif
!     endif

 !do while (multip(j).gt.1) 
             	
!       Ncopy=Ncopy+one
!      ToCopy(Ncopy-1)=j
!        kaq4 = int((xpq4(j)/x111)*nfenetre)
!       kaq6 = int((xpq6(j)/xmaxq6)*nfenetre)
!        ien  = int((ener(j)-energie_min)/(energie_max-energie_min)*real(nfenetre))
!       stat(kaq4,ien)=stat(kaq4,ien)+1
!   !  statpos(kaq4,ien)=statpos(kaq4,ien)+1
!      stat6(kaq6,ien)=stat6(kaq6,ien)+1
!      !   write (*,*) 'clonage clone # avec norme',j, number(j), ncopy, tocopy(ncopy-1)
!	multip(j)=multip(j)-1
!	if((Ncopy.gt.NMax).or.(Ncopy.eq.NMax)) then
!	   write(6,*)'Too many clones'
!	    stop
!      endif
!enddo

!write(*,*) 'systr',ncopy, ndelete
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 
!if (sys.eq.0) then
  !How many offsprings for the clone j ?
      !Number(j)=floor(Weight(j)-0.5+genrand())
      Number(j)=floor(Weight(j)+genrand())
      !  write (*,*) 'norme intiere de j' ,j, number(j)
      
 !if zero -> delete it
      if (Number(j).lt.1) then

	    ToDelete(Ndelete)=j 
            Ndelete=Ndelete+one
            ien  = int((ener(j)-energie_min)/(energie_max-energie_min)*real(nfenetre))
            kaq4 = int((xpq4(j)/x111)*nfenetre)
            kaq6 = int((xpq6(j)/xmaxq6)*nfenetre)
            stat(kaq4,ienergy)=stat(kaq4,ienergy)-1
            !statneg(kaq4,ienergy)=statneg(kaq4,ienergy)+1
            stat6(kaq6,ien)=stat6(kaq6,ien)-1
            nvite(j)=nvite(j)+1
             tvmoy(j)=((nvite(j)-1)*tvmoy(j)+tv(j))/real(nvite(j))
            tv(j)=0

             !write (*,*)'detruit clone # avec norme', j, number(j)!, Ndelete,ToDelete(Ndelete-1) 
             if(Ndelete.gt.NbClones) then
	      !   write(6,*)'All the clones have been killed'
	         stop
             endif
          
       else if (Number(j).gt.1) then
         copiati=copiati+1

        !  write (*,*) 'clonage clone # avec norme',j, weight(j),number(j)
       !  if (time.gt.1800) then
         !!!!!!!!!trsm
         !jtrans=jtrans+weight(j)

         !!!!!!!!trsm2
       !  jtrans=jtrans+(weight(j)-1.d0)
       !  endif
         !write(*,*) 'jtrans', jtrans

         tv(j)=tv(j)+dt
           if (nvite(j).eq.0) then 
               tvmoy(j)=tv(j)
           endif
       else if (Number(j).eq.1) then
         n0=n0+1
         tv(j)=tv(j)+dt
         if (nvite(j).eq.0) then 
            tvmoy(j)=tv(j)
         endif
       endif

 ! if strictly more than one offspring -> cloning
      do while (Number(j).gt.1) 
        Ncopy=Ncopy+one
        ToCopy(Ncopy-1)=j
         kaq4 = int((xpq4(j)/x111)*nfenetre)
         kaq6 = int((xpq6(j)/xmaxq6)*nfenetre)
         ien  = int((ener(j)-energie_min)/(energie_max-energie_min)*real(nfenetre))
         stat(kaq4,ien)=stat(kaq4,ien)+1
         ! statpos(kaq4,ien)=statpos(kaq4,ien)+1
         stat6(kaq6,ien)=stat6(kaq6,ien)+1
        !   write (*,*) 'clonage clone # avec norme',j, number(j), ncopy, tocopy(ncopy-1)
	Number(j)=Number(j)-1
	if((Ncopy.gt.NMax).or.(Ncopy.eq.NMax)) then
	!   write(6,*)'Too many clones'
	    stop
        endif
      enddo

 

!!!!!!!SYSTR!!!!!!!!
!ENDIF
!!!!!!!!!!!!!
enddo

!if(time.gt.1800) then

!!!!!!!trsm
!jtrans=jtrans+real(n0*w)

!trasm=real(jtrans)/real(jtot)

!!!!!!trsm
!trasm1=real(jtrans)/real(N)
!trasm2=real(jtrans)/real(ndelete)

!!!!!!!!trsm2
!trasm1=real(jtrans)/real(N)
!trasm2=trasm2*trasm1

!endif

!write(*,*) 'trasmette', trasm

!trastot=trastot*trasm
!pi2=real(Ncopy)/real(nbclones)
!pi2=real(n0+ncopy)/real(nbclones)
pif=real(Nbclones-ndelete+ncopy)/real(nbclones)


logt=logt+log(pif)

!logt2=logt2+log(pi2)


!if(mod(iter,300).eq.0) then
!  write(*,*) 'statistica', time, Ndelete, Ncopy, n0, nbclones
!endif
     ! At this stage, we have a list of N1 clones to delete in
      !ToDelete and a list of N2 clones to copy in ToCopy. If N1=N2,
      !then the clones to be copied are copied where those who have
      !to be deleted are stored.
     
 !     If N1.gt.N2, there are more clones to delete than to
 !     clone. First, we copy the N2 clones to be copied on N2 clones
 !     to be deleted. We are left with N1-N2 clones still to be
 !     deleted. The idea is to pull at random N1-N2 clones _among the
  !    rest_ and to put them in place of the the N1-N2 still to be deleted.
    
 !     If N2.gt.N1, there are more clones to be copied that to be
 !     killed. We first copy N1 clones to be copied on the N1 clones
 !     to be deleted. Then, one chooses N2-N1 clones at random _among
 !     the NbClones+N2-N1 which we have and delete them to go back to
 !     NbClones clones.
  
  !    we copy what's possible to copy
!if (lo.eq.1) then
   do j=0,MIN(Ndelete-1,Ncopy-1)
      ! write (*,*) 'copie',j, ndelete, todelete(j), tocopy(j)
      do l=0,N-1
        P(ToDelete(j))%qx(l) =P(ToCopy(j))%qx(l) 
	P(ToDelete(j))%px(l) =P(ToCopy(j))%px(l) 
	P(ToDelete(j))%uqx(l)=P(ToCopy(j))%uqx(l)
	P(ToDelete(j))%upx(l)=P(ToCopy(j))%upx(l)
  
        P(ToDelete(j))%qy(l) =P(ToCopy(j))%qy(l) 
	P(ToDelete(j))%py(l) =P(ToCopy(j))%py(l) 
	P(ToDelete(j))%uqy(l)=P(ToCopy(j))%uqy(l)
	P(ToDelete(j))%upy(l)=P(ToCopy(j))%upy(l)

        P(ToDelete(j))%qz(l) =P(ToCopy(j))%qz(l) 
	P(ToDelete(j))%pz(l) =P(ToCopy(j))%pz(l) 
	P(ToDelete(j))%uqz(l)=P(ToCopy(j))%uqz(l)
	P(ToDelete(j))%upz(l)=P(ToCopy(j))%upz(l)

        enddo
        betaq(ToDelete(j))=betaq(ToCopy(j))
	dh(ToDelete(j))=dh(ToCopy(j))
     enddo
     !  P(ToDelete(j)).Lyap =P(ToCopy(j)).Lyap 
 
!!!!!!!!!!!!SYSTR!!!!!!!!!!!!!!
!if(sys.eq.0) then   
!!!!!!!!!!SYSTR!!!!!!!!!!!!

    
    if(Ncopy.gt.Ndelete) then
 ! We put the clones still to be copied at the end of the
 !  array of clones.
      do j=Ndelete,Ncopy-1
     !  write(*,*) ndelete, ncopy, j, NbClones+j-Ndelete
	do l=0,N-1
	  P(NbClones+j-Ndelete)%qx(l) =P(ToCopy(j))%qx(l) 
	  P(NbClones+j-Ndelete)%px(l) =P(ToCopy(j))%px(l) 
	  P(NbClones+j-Ndelete)%uqx(l)=P(ToCopy(j))%uqx(l)
	  P(NbClones+j-Ndelete)%upx(l)=P(ToCopy(j))%upx(l)

          P(NbClones+j-Ndelete)%qy(l) =P(ToCopy(j))%qy(l) 
	  P(NbClones+j-Ndelete)%py(l) =P(ToCopy(j))%py(l) 
	  P(NbClones+j-Ndelete)%uqy(l)=P(ToCopy(j))%uqy(l)
	  P(NbClones+j-Ndelete)%upy(l)=P(ToCopy(j))%upy(l)

          P(NbClones+j-Ndelete)%qz(l) =P(ToCopy(j))%qz(l) 
	  P(NbClones+j-Ndelete)%pz(l) =P(ToCopy(j))%pz(l) 
	  P(NbClones+j-Ndelete)%uqz(l)=P(ToCopy(j))%uqz(l)
	  P(NbClones+j-Ndelete)%upz(l)=P(ToCopy(j))%upz(l)

         enddo
         betaq(NbClones+j-Ndelete)=betaq(ToCopy(j))
	 dh(NbClones+j-Ndelete)=dh(ToCopy(j))

	!P(NbClones+j-Ndelete).Lyap =P(ToCopy(j)).Lyap 
      enddo
   
      
      ! We now have to delete Ncopy - Ndelete clones to keep the
      ! population constant
      
      Ndelete=Ncopy-Ndelete
     

!write(*,*) 'copie clones'
      
   !   We pull at random Ndelete among the clones
   !   and store them in ToDelete
!write (*,*)'combinaison',Ndelete,NbClones+Ndelete
     call Combinaison(Ndelete,NbClones+Ndelete,ToDelete,Nb)
   ! In 'todelete' we now have Ndelete different rank chosen at random.

!   we reorder Nb() for a future use

      do j=0,Ndelete-1
       !write (*,*) j, ndelete,ToDelete(j)
	Nb(ToDelete(j))=ToDelete(j)
      enddo    
  
    !   Sort the array ToDelete in ascending numerical order
      call shell(Ndelete,ToDelete-1)
      
     !  We delete the clones, starting by the last one.
  !    To do so, we replace them by the last clone of the array
  !    of particles
      do j=Ndelete-1,0,-1
         
	 do l=0,N-1
	  P(ToDelete(j))%qx(l) =P(NbClones+Ndelete-1)%qx(l) 
	  P(ToDelete(j))%px(l) =P(NbClones+Ndelete-1)%px(l) 
	  P(ToDelete(j))%uqx(l)=P(NbClones+Ndelete-1)%uqx(l)
	  P(ToDelete(j))%upx(l)=P(NbClones+Ndelete-1)%upx(l)

          P(ToDelete(j))%qy(l) =P(NbClones+Ndelete-1)%qy(l) 
	  P(ToDelete(j))%py(l) =P(NbClones+Ndelete-1)%py(l) 
	  P(ToDelete(j))%uqy(l)=P(NbClones+Ndelete-1)%uqy(l)
	  P(ToDelete(j))%upy(l)=P(NbClones+Ndelete-1)%upy(l)

          P(ToDelete(j))%qz(l) =P(NbClones+Ndelete-1)%qz(l) 
	  P(ToDelete(j))%pz(l) =P(NbClones+Ndelete-1)%pz(l) 
	  P(ToDelete(j))%uqz(l)=P(NbClones+Ndelete-1)%uqz(l)
	  P(ToDelete(j))%upz(l)=P(NbClones+Ndelete-1)%upz(l)
	enddo
        betaq(ToDelete(j))=betaq(NbClones+Ndelete-1)
	 dh(ToDelete(j))=dh(NbClones+Ndelete-1)
        
         	!P(ToDelete(j)).Lyap=P(NbClones+Ndelete-1).Lyap
	Ndelete=ndelete-1
      enddo

    endif

    if(Ncopy.lt.Ndelete) then
    !MOD  if((Ncopy.lt.Ndelete).and.(Ncopy.ne.0))then
  ! write(*,*)'nd gt nc', ndelete, ncopy
     !MOD  ran = floor(genrand()*(ncopy))
      
    !M do m=ndelete-1,Ncopy,-1
   !write(*,*)  m, ToDelete(m),Tocopy(ran)
   !M	do k=0,N-1
	!M  P(todelete(m))%qx(k) =P(Tocopy(ran))%qx(k)  
	!m  P(todelete(m))%px(k) =P(tocopy(ran))%px(k)  
         !m P(todelete(m))%uqx(k)=P(tocopy(ran))%uqx(k)  
	!m  P(todelete(m))%upx(k)=P(tocopy(ran))%upx(k)  

        !m  P(todelete(m))%qy(k) =P(tocopy(ran))%qy(k)  
	!m  P(todelete(m))%py(k) =P(tocopy(ran))%py(k)  
	!m  P(todelete(m))%uqy(k)=P(tocopy(ran))%uqy(k)  
	!  P(todelete(m))%upy(k)=P(tocopy(ran))%upy(k)

        !  P(todelete(m))%qz(k) =P(tocopy(ran))%qz(k)  
	!  P(todelete(m))%pz(k) =P(tocopy(ran))%pz(k)  
	!  P(todelete(m))%uqz(k)=P(tocopy(ran))%uqz(k)  
	!  P(todelete(m))%upz(k)=P(tocopy(ran))%upz(k)
	!enddo
!        betaq(todelete(m))=betaq(tocopy(ran))
	! dh(todelete(m))=dh(tocopy(ran))
	!P(NbClones-l).Lyap=P(ran).Lyap  
    ! enddo


  ! else if((Ncopy.lt.Ndelete).and.(Ncopy.eq.0))then


! write(*,*) 'copie-delete'   
    
 !   if there are more clone to delete than to copy
 !   if(Ncopy.lt.Ndelete) then
!write(*,*)'nd gt nc=0', ndelete, ncopy
      l=1
  !  We put the dead guys at the end

      do j=Ndelete-1,Ncopy,-1
       ! write(*,*)  j, ToDelete(j),NbClones-l
        do k=0,N-1
	  P(ToDelete(j))%qx(k) =P(NbClones-l)%qx(k)  
	  P(ToDelete(j))%px(k) =P(NbClones-l)%px(k)  
	  P(ToDelete(j))%uqx(k)=P(NbClones-l)%uqx(k)
	  P(ToDelete(j))%upx(k)=P(NbClones-l)%upx(k)

          P(ToDelete(j))%qy(k) =P(NbClones-l)%qy(k)  
	  P(ToDelete(j))%py(k) =P(NbClones-l)%py(k)  
	  P(ToDelete(j))%uqy(k)=P(NbClones-l)%uqy(k)
	  P(ToDelete(j))%upy(k)=P(NbClones-l)%upy(k)
 
          P(ToDelete(j))%qz(k) =P(NbClones-l)%qz(k)  
	  P(ToDelete(j))%pz(k) =P(NbClones-l)%pz(k)  
	  P(ToDelete(j))%uqz(k)=P(NbClones-l)%uqz(k)
	  P(ToDelete(j))%upz(k)=P(NbClones-l)%upz(k)
        enddo
         tv(ToDelete(j))=tv(NbClones-l)
         tvmoy(ToDelete(j))=tv(NbClones-l)
        betaq(ToDelete(j))=betaq(NbClones-l)
	 dh(ToDelete(j))=dh(NbClones-l)
	l=l+1
      enddo

      
   !   We have to clone Ndelete-Ncopy clone to keep the population constant
      
      Ncopy=Ndelete-Ncopy
      l=1
      do m=0,Ncopy-1
	ran = floor(genrand()*(NbClones-Ncopy))
!write(*,*) NbClones-l,ran
	do k=0,N-1
	  P(NbClones-l)%qx(k) =P(ran)%qx(k)  
	  P(NbClones-l)%px(k) =P(ran)%px(k)  
	  P(NbClones-l)%uqx(k)=P(ran)%uqx(k)  
	  P(NbClones-l)%upx(k)=P(ran)%upx(k)  

          P(NbClones-l)%qy(k) =P(ran)%qy(k)  
	  P(NbClones-l)%py(k) =P(ran)%py(k)  
	  P(NbClones-l)%uqy(k)=P(ran)%uqy(k)  
	  P(NbClones-l)%upy(k)=P(ran)%upy(k)

          P(NbClones-l)%qz(k) =P(ran)%qz(k)  
	  P(NbClones-l)%pz(k) =P(ran)%pz(k)  
	  P(NbClones-l)%uqz(k)=P(ran)%uqz(k)  
	  P(NbClones-l)%upz(k)=P(ran)%upz(k)
	enddo
         ! tv(NbClones-l)=tv(ran)
        ! tvmoy(NbClones-l)=tv(ran)
        betaq(NbClones-l)=betaq(ran)
	 dh(NbClones-l)=dh(ran)
	!P(NbClones-l).Lyap=P(ran).Lyap  
	l=l+1
      enddo
     endif

 
do m=0,NbClones-1!
   call Normalizza(P(m)%uqx,P(m)%uqy,P(m)%uqz,P(m)%upx,P(m)%upy,P(m)%upz, N)
   !call CalculNorme(P(m)%uqx,P(m)%uqy,P(m)%uqz,P(m)%upx,P(m)%upy,P(m)%upz,N,norm)
  ! write(*,*) 'renorm' , norm
  enddo






   
endif !!!!!!!!!!!!! FINE CLONING
!!!!!!!!!!!!!!!!!!!SYSTR!!!!!!
!ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!q4moy=0
!q6moy=0
!u=0
!xpq4(:)=0
!xpq6(:)=0
!ener(:)=0


!write (*,*) 'cloning is over'

do l=0,NbClones-1   !! cloni


!    do i=1,N-1
!      do j=i+1,N
!        xij(i,j)=(P(l)%qx(i-1)-P(l)%qx(j-1))
!        yij(i,j)=(P(l)%qy(i-1)-P(l)%qy(j-1))
!        zij(i,j)=(P(l)%qz(i-1)-P(l)%qz(j-1))
!        dij(i,j)=sqrt(xij(i,j)**2+yij(i,j)**2+zij(i,j)**2)
!        xij(j,i)=-xij(i,j)
!        yij(j,i)=-yij(i,j)
 !       zij(j,i)=-zij(i,j)
!        dij(j,i)= dij(i,j)
!      enddo
!    enddo

!    xpq4(L)=q4(N,xij,yij,zij,dij)

!write(*,*)'www', xpq4(l)

 !  call fq4(N,xij,yij,zij,dij,fq4x,fq4y,fq4z) 
!    xpq6(L)=q6(N,xij,yij,zij,dij)

 !   call fq6(N,xij,yij,zij,dij,fq6x,fq6y,fq6z)
!    call Force(P(l)%qx,P(l)%qy,P(l)%qz,fp,Utot)

!    ener(l)=utot
!    q4moy=q4moy+xpq4(l)
!    q6moy=q6moy+xpq6(l)
!    u=u+ener(L)


     w4(l)=0
     w6(l)=0


   tvmoytot=sum(tv(:))/nbclones

!!!!!!!!!!ISTOGRAMMA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ifenetre = int((xpq4(l)/x111)*real(nfenetre))
   ifen_q6 = int((xpq6(l)/xmaxq6)*real(nfenetre))
   ienergy  = int((ener(l)-energie_min)/(energie_max-energie_min)*real(nfenetre))
   tab_contour(-10:nfenetre+10,0:nfenetre)=zero
   tab_contour_q6(-10:nfenetre+10,0:nfenetre)=zero
   tab_cont_q4q6(-10:nfenetre+10,-10:nfenetre+10)=zero

   
   !!!!!CALPOT_AUXI
   lambdaq4 = (xpq4(L)/x111)*real(nfenetre)-real(ifenetre)

   if (ifenetre.gt.nfenetre+10) ifenetre = nfenetre+10
   icfen = ifenetre+1 
   if (icfen.gt.nfenetre+10) icfen = nfenetre+10
   if (ifenetre.lt.-10) ifenetre = -10

   lambdaE = ((potist-energie_min)/(energie_max-energie_min))*real(nfenetre)-real(ienergy)

   if (ienergy.gt.nfenetre) ienergy = nfenetre
   if (ienergy.lt.0) ienergy = 0
   icener = ienergy+1 
   if (icener.gt.nfenetre) icener = nfenetre

   lambdaq6 = (xpq6(L)/xmaxq6)*real(nfenetre)-real(ifen_q6)

   if (ifen_q6.gt.nfenetre+10) ifen_q6 = nfenetre+10
   icfen_q6 = ifen_q6+1 
   if (icfen.gt.nfenetre+10) icfen_q6 = nfenetre+10
   if (ifen_q6.lt.-10) ifen_q6 = -10

 !!!!!!!FINE CALPOT AUXI

  !!!!histogramme q4 - energie

        if ((ifenetre.le.nfenetre+10).and.(ifenetre.ge.-10)) then
           if ((ienergy.le.nfenetre).and.(ifenetre.ge.0)) then
             tab_contour(ifenetre,ienergy) = tab_contour(ifenetre,ienergy) +one
             cumul_contour(ifenetre,ienergy)=1+cumul_contour(ifenetre,ienergy)
           endif
        endif

   !histogramme q6-energie
        
        if ((ifen_q6.le.nfenetre+10).and.(ifen_q6.ge.-10)) then
            if ((ienergy.le.nfenetre).and.(ienergy.ge.0)) then
               tab_contour_q6(ifen_q6,ienergy)=1+tab_contour_q6(ifen_q6,ienergy)
               cumul_contour_q6(ifen_q6,ienergy)=1+cumul_contour_q6(ifen_q6,ienergy)
            endif
        endif


   !histogramme q4-q6 (fait par moi)
          if ((ifenetre.le.nfenetre+10).and.(ifenetre.ge.-10)) then
            if ((ifen_q6.le.nfenetre+10).and.(ifen_q6.ge.-10)) then
                tab_cont_q4q6(ifenetre,ifen_q6)=1+tab_cont_q4q6(ifenetre,ifen_q6)
                cumul_contour_q4q6(ifenetre,ifen_q6)=1+cumul_contour_q4q6(ifenetre,ifen_q6)
                !write(*,*) 'rr' ,  cumul_contour_q4q6(ifenetre,ifen_q6)
          endif 
         endif
!endif

!!!!!!!!!!!!!!! FINE ISTOGRAMMA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
        do i=1,N-1
          do j=i+1,N
            xij(i,j)=(P(l)%qx(i-1)+P(l)%uqx(i-1)-P(l)%qx(j-1)-P(l)%uqx(j-1))
            yij(i,j)=(P(l)%qy(i-1)+P(l)%uqy(i-1)-P(l)%qy(j-1)-P(l)%qy(j-1))
            zij(i,j)=(P(l)%qz(i-1)+P(l)%uqz(i-1)-P(l)%qz(j-1)-P(l)%qz(j-1))
            dij(i,j)=sqrt(xij(i,j)**2+yij(i,j)**2+zij(i,j)**2)
            xij(j,i)=-xij(i,j)
            yij(j,i)=-yij(i,j)
            zij(j,i)=-zij(i,j)
            dij(j,i)= dij(i,j)
          enddo
       enddo
       w4(l)=q4(N,xij,yij,zij,dij)
       w6(l)=q6(N,xij,yij,zij,dij)
     


   xppro(1,1:N)=P(l)%qx+P(l)%uqx
   xppro(2,1:N)=P(l)%qy+P(l)%uqy
   xppro(3,1:N)=P(l)%qz+P(l)%uqz

   call calfoljc2(xppro,fppro,uproject)

    enerpro(l)=uproject


 
if (lo.eq.2) then 
!!!!!!!!!!!!!!!!!!!METODO 1!!!!!!!!!!!!!!!!
!POUR CALCULER LA LONGUEUR DES FLECHES DES VECTEURS on met le q+dq  dans l'espace de phases (Q4,Q6) mais PAS les p+dp car ils sont pas dans Q4,Q6  
! on calcule le Q pour la configuration courante

  q4ref=xpq4(l)
  q6ref=xpq6(l)
   
! On sait que Q=Q(xij,yin,zij), c'est-à-dire des distances réciproques entre les atomes. Les éléments de grad Q sont donc donnés par 
!dQ/dx_i=(dQ/dx_ij)*(dx_ij/dx_i)=(dQ/dx_ij) car la deuxième dérivée est égale à + ou - 1. On calcule (dQ/dx_ij)={Q(eps_ij+x_ij, …)-Q(…)}/eps
! Première partie du gradient: boucle sur x
!write(*,*)'x', k, delt(k)

   do k=0,N-1
     P(l)%qx(k)=P(l)%qx(k)+e
        do i=1,N-1
          do j=i+1,N
            xij(i,j)=(P(l)%qx(i-1)-P(l)%qx(j-1))
            yij(i,j)=(P(l)%qy(i-1)-P(l)%qy(j-1))
            zij(i,j)=(P(l)%qz(i-1)-P(l)%qz(j-1))
            dij(i,j)=sqrt(xij(i,j)**2+yij(i,j)**2+zij(i,j)**2)
            xij(j,i)=-xij(i,j)
            yij(j,i)=-yij(i,j)
            zij(j,i)=-zij(i,j)
            dij(j,i)= dij(i,j)
          enddo
       enddo
       gQ4x(k+1)=((q4(N,xij,yij,zij,dij)-q4ref)/real(e))
       gQ6x(k+1)=((q6(N,xij,yij,zij,dij)-q6ref)/real(e))
       p(l)%qx(k)=p(l)%qx(k)-e
    enddo

    do k=0,N-1
      P(l)%qy(k)=P(l)%qy(k)+e
      do i=1,N-1
        do j=i+1,N
          xij(i,j)=(P(l)%qx(i-1)-P(l)%qx(j-1))
          yij(i,j)=(P(l)%qy(i-1)-P(l)%qy(j-1))
          zij(i,j)=(P(l)%qz(i-1)-P(l)%qz(j-1))
          dij(i,j)=sqrt(xij(i,j)**2+yij(i,j)**2+zij(i,j)**2)
          xij(j,i)=-xij(i,j)
          yij(j,i)=-yij(i,j)
          zij(j,i)=-zij(i,j)
          dij(j,i)= dij(i,j)
        enddo
      enddo
      gQ4y(k+1)=((q4(N,xij,yij,zij,dij)-q4ref)/real(e))
     gQ6y(k+1)=((q6(N,xij,yij,zij,dij)-q6ref)/real(e))
     p(l)%qy(k)=p(l)%qy(k)-e
   enddo

    do k=0,N-1
       P(l)%qz(k)=P(l)%qz(k)+e
        do i=1,N-1
          do j=i+1,N
            xij(i,j)=(P(l)%qx(i-1)-P(l)%qx(j-1))
            yij(i,j)=(P(l)%qy(i-1)-P(l)%qy(j-1))
            zij(i,j)=(P(l)%qz(i-1)-P(l)%qz(j-1))
            dij(i,j)=sqrt(xij(i,j)**2+yij(i,j)**2+zij(i,j)**2)
            xij(j,i)=-xij(i,j)
            yij(j,i)=-yij(i,j)
            zij(j,i)=-zij(i,j)
           dij(j,i)= dij(i,j)
         enddo
        enddo
     gQ4z(k+1)=((q4(N,xij,yij,zij,dij)-q4ref)/real(e))
    gQ6z(k+1)=((q6(N,xij,yij,zij,dij)-q6ref)/real(e))
   p(l)%qz(k)=p(l)%qz(k)-e
  enddo

   normgradq4=sum(gQ4x(:)**2+gQ4y(:)**2+gQ4z(:)**2)
   normgradq6=sum(gQ6x(:)**2+gQ6y(:)**2+gQ6z(:)**2)

   w4(l)=(dot_product(P(l)%uqx(:),gQ4x(:))+dot_product(P(l)%uqy(:),gQ4y(:))+dot_product(P(l)%uqz(:),gQ4z(:)))/normgradq4
   w6(l)=(dot_product(P(l)%uqx(:),gQ6x(:))+dot_product(P(l)%uqy(:),gQ6y(:))+dot_product(P(l)%uqz(:),gQ6z(:)))/normgradq6


   xppro(1,1:N)=P(l)%qx
   xppro(2,1:N)=P(l)%qy
   xppro(3,1:N)=P(l)%qz

   do k=0,N-1
    xppro(1,k+1)=P(l)%qx(k)+e
    call calfoljc2(xppro,fppro,uproject)
    gEx(k+1)=(uproject-ener(l))/real(e)
    xppro(1,k+1)=xppro(1,k+1)-e
  enddo

  uproject=0

   do k=0,N-1
    xppro(2,k+1)=P(l)%qy(k)+e
     call calfoljc2(xppro,fppro,uproject)
     gEy(k+1)=(uproject-ener(l))/real(e)
    xppro(2,k+1)=xppro(2,k+1)-e
     enddo

       uproject=0


    do k=0,N-1
     xppro(3,k+1)=P(l)%qz(k)+e
     call calfoljc2(xppro,fppro,uproject)
     gEz(k+1)=(uproject-ener(l))/real(e)
     xppro(3,k+1)=xppro(3,k+1)-e
   enddo

     uproject=0


     enerpro(l)=(dot_product(P(l)%uqx(:),gEx(:))+dot_product(P(l)%uqy(:),gEy(:))+dot_product(P(l)%uqz(:),gEz(:)))

       xppro(1,1:N)=0
       xppro(2,1:N)=0
       xppro(3,1:N)=0

endif

enddo

q4moy=q4moy/NbClones
q6moy=q6moy/NbClones
u=u/NbClones

w4(:)=0
w6(:)=0
enerpro(:)=0

nq4=0

tot=(NbClones)*(TotalTime/dt) !

!if (lo.eq.3) then
if(mod(iter,100).eq.0) then
write(*,'(a, e18.7, e18.7, e18.7, e18.7)')'param', time, q4moy ,q6moy, U!, logt, pif
!write (*,'(a e11.4, e11.4, e11.4, e12.4, e12.4, e11.4)') 'param', time, q4moy ,q6moy,U, logt, pif!, logt!,tvmoytot
endif

if(iter.eq.0) then
do i=0, Nbclones-1
write(28,'(e11.4, e11.4, e11.4, e12.4, e12.4, e12.4)') time ,xpq4(i),xpq6(i),ener(i), weight(i),w
enddo
endif

if(mod(iter,500).eq.0) then
do i=0, Nbclones-1
write(28,'(e11.4, e11.4, e11.4, e12.4, e12.4, e12.4)') time ,xpq4(i),xpq6(i),ener(i), weight(i),w
enddo
endif
!if (time.gt.2000) then
!if(mod(iter,500).eq.0) then
!do i=0, Nbclones-1
!write(28,'(e11.4, e11.4, e11.4, e12.4, e12.4, e12.4, e12.4, e12.4, e12.4)') time ,xpq4(i),xpq6(i),ener(i),w4(i),w6(i),enerpro(i), weight(i),w
!enddo
!endif
!write(28,*)
!else if(mod(iter,200).eq.10) then
!do i=0, Nbclones-1
!write(28,'(e11.4, e11.4, e11.4, e12.4, e12.4, e12.4, e12.4, e12.4, e12.4)') time ,xpq4(i),xpq6(i),ener(i),w4(i),w6(i),enerpro(i), weight(i),w
!enddo
!write(28,*)
!else if(mod(iter,200).eq.20) then
!do i=0, Nbclones-1
!write(28,'(e11.4, e11.4, e11.4, e12.4, e12.4, e12.4, e12.4, e12.4, e12.4)') time ,xpq4(i),xpq6(i),ener(i),w4(i),w6(i),enerpro(i), weight(i),w
!enddo
!write(28,*)
!else if(mod(iter,200).eq.30) then
!do i=0, Nbclones-1
!write(28,'(e11.4, e11.4, e11.4, e12.4, e12.4, e12.4, e12.4, e12.4, e12.4)') time ,xpq4(i),xpq6(i),ener(i),w4(i),w6(i),enerpro(i), weight(i),w
!enddo
!write(28,*)
!else if(mod(iter,200).eq.40) then
!do i=0, Nbclones-1
!write(28,'(e11.4, e11.4, e11.4, e12.4, e12.4, e12.4, e12.4, e12.4, e12.4)') time ,xpq4(i),xpq6(i),ener(i),w4(i), w6(i),enerpro(i), weight(i),w
!enddo
!write(28,*)
!endif
!endif

if(mod(iter,1000).eq.0) then
write(29,*) time
do m=0,NbClones-1
!write(29,*) m
   write (29,*) P(m)%qx(:)
   write (29,*) P(m)%qy(:)
   write (29,*) P(m)%qz(:)
   write (29,*) 
   write (29,*) P(m)%px(:)
   write (29,*) P(m)%py(:)
   write (29,*) P(m)%pz(:)      
   write (29,*)
enddo 

!do m=0,NbClones-1
!   write (29,*) P(m)%uqx(:)
!   write (29,*) P(m)%uqy(:)
!   write (29,*) P(m)%uqz(:)
!   write (29,*) 
 !  write (29,*) P(m)%upx(:)
!   write (29,*) P(m)%upy(:)
!   write (29,*) P(m)%upz(:)      
!   write (29,*)
!enddo
!do m=0, Nbclones-1
!   write(29,*) '##clone', m
!   write(29,*)
!   write (29,'(e11.4)') real(P(m)%qx(:))
!   write (29,'(e11.4)') real(P(m)%qy(:))
!   write (29,'(e11.4)') real(P(m)%qz(:))
!   write (29,'(e11.4)') 
!   write (29,'(e11.4)') real(P(m)%px(:))
!   write (29,'(e11.4)') real(P(m)%py(:))
!   write (29,'(e11.4)') real(P(m)%pz(:))      
!   write (29,'(e11.4)')
!   write (29,'(e11.4)') real(P(m)%uqx(:))
!   write (29,'(e11.4)') real(P(m)%uqy(:))
!   write (29,'(e11.4)') real(P(m)%uqz(:))
!   write (29,'(e11.4)') 
!   write (29,'(e11.4)') real(P(m)%upx(:))
!   write (29,'(e11.4)') real(P(m)%upy(:))
!   write (29,'(e11.4)') real(P(m)%upz(:))      
!   write (29,'(e11.4)')
!enddo
!write (29,*) '#end ', time
endif


if(mod(iter,5000).eq.0) then
write(33,*)
write(33,*) '#debut',time
do i=0,nfenetre
    do j=0,nfenetre
      write (33,'(e12.4, e12.4, e12.4 )') real(i)/500,real(j)*3/10-174,stat(i,j)
     enddo
   write (33,*)
enddo
write (33,*) '#end ', time
write (33,*)
write (34,*)
write(34,*) '#debut',time
do i=0,nfenetre
    do j=0,nfenetre
      write (34,'(e12.4, e12.4, e12.4 )') real(i)*3/500,real(j)*3/10-174,stat6(i,j)
     enddo
   write (34,*)
enddo
write (34,*) '#end ', time
write (34,*)
endif

!endif
!if (time.lt.8000) then
!if(mod(iter,100000).eq.0) then
!do i=0, Nbclones-1
!write(*,'(a, e18.7, e18.7, e18.7, e18.7, e18.7)') 'snap',xpq4(i),xpq6(i),ener(i),weight(i),w
!enddo
!write(*,*)!'tv', time, iter
!endif
!endif

tvmoytot=0
 

!!!!!ADESSO RINORMALIZZO IL VETTORE USCENTE DAL CLONAGGIO, CHE SARÀ PRONTO PER LA PROPAGAZIONE AL TEMPO t+dt

!if(mod(iter,niter).eq.0) then
!do m=0,NbClones-1!
!   call Normalizza(P(m)%uqx,P(m)%uqy,P(m)%uqz,P(m)%upx,P(m)%upy,P(m)%upz, N)
   !call CalculNorme(P(m)%uqx,P(m)%uqy,P(m)%uqz,P(m)%upx,P(m)%upy,P(m)%upz,N,norm)
  ! write(*,*) 'renorm' , norm
!  enddo

!endif
!write(*,*) 'temp',  time, temperature,lo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! annealing pour JT
!if ((time.gt.10).and.(time.lt.800)) then
!temperature=0.16
!endif

!if ((time.gt.300).and.(time.lt.3200)) then
!lo=1
!endif!
!!!!!!!!!!!!!!!!!!!!!! ANNEALING
if (temperature.gt.0.05) then
!if ((time.ge.800).and.(time.lt.3200)) then
    temperature=temperature-inter
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!if (time.ge.3200) then
!  temperature=0.05  !lo=1
!gamma=100.d0
!write (*,*) 'temp',time, temperature, kinetotdt
!  endif
!!!!!!!!!!!!!!!!!!!!!! quenching
  ! if (time.gt.500) then
  !  temperature=0.03  !lo=1
!gamma=100.d0
!if(mod(iter,1000).eq.0) then
!write (*,*) 'temp',time, temperature, kinetotdt
! endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!! CLONING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!if ((time.gt.28000).and.(temperature.le.0.04)) then
!if (time.gt.50000) then
! lo=1
!endif!

!if (time.gt.100) then
!temperature=0.19
!endif






!else if (lo.eq.0) then

!if(mod(iter,1000).eq.0) then
! write (*,'(a, e18.7, e18.7, e18.7, e18.7, e18.7,e18.7)')'param',  time, q4moy ,q6moy,U, logt, pif
 !write (*,'(a, e18.7, e18.7, e18.7, e18.7, e18.7,e18.7, e18.7)') 'metodo', time, q4moy,q6moy,w4moy,w6moy,U,enprmoy
 !write (*,*) 'noclon'
!endif 

!endif


!if (time.gt.4000)then

!lo=0

!endif



!if(mod(iter,500).eq.0) then
! write (*,'(a, e18.7, e18.7, e18.7, e18.7)')'param',  time, q4moy ,q6moy,U!, logt, pi!,tvmoytot
   !write (*,'(a, e18.7, e18.7, e18.7, e18.7, e18.7,e18.7, e18.7)') 'metodo', time, q4moy,q6moy,w4moy,w6moy,U,enprmoy
 !write (*,*) 'tvmoytot',  q4moy, 
!endif


!do m=0,NbClones-1 !!!!!!!!!!!!!!
!call CalculNorme(P(m)%uqx,P(m)%uqy,P(m)%uqz,P(m)%upx,P(m)%upy,P(m)%upz,N,norm)
   !write(*,*) 'renorm0' , norm
 !!!!!!!call Normalizza(P(m)%uqx,P(m)%uqy,P(m)%uqz,P(m)%upx,P(m)%upy,P(m)%upz, N)!!!!!!!!!!!!!
!enddo !!!!!!!!!!!!!!!!

!write(11,*)'non clono'!!!!!!!!!!!!!!!!!!

!endif


!if(mod(iter,100).eq.0) then
!write(*,*) 'gna', nbclones
!endif
!  mu += log(weight/(RealNumber) NbClones)

!!!TEMPS FINI.RIPARTE IL CICLO
Time=Time+dt

enddo

 !'tvita medio finale',tvmoytot 

stat(:,:)=stat(:,:)/(((totaltime-tequilib)/(niter*dt))*Nbclones)
stat6(:,:)=stat6(:,:)/(((totaltime-tequilib)/(niter*dt))*Nbclones)

!statneg(:,:)=statneg(:,:)/(((totaltime-tequilib)/(niter*dt))*Nbclones)
!statpos(:,:)=statpos(:,:)/(((totaltime-tequilib)/(niter*dt))*Nbclones)

do i=0,nfenetre
    do j=0,nfenetre
      write (11,*) 'stat',real(i)/500,real(j)*3/10-174,stat(i,j)
     enddo
   write (11,*)
enddo

do i=0,nfenetre
    do j=0,nfenetre
      write (11,*) 'stat6',real(i)*3/500,real(j)*3/10-174,stat6(i,j)
     enddo
   write (11,*)
enddo

!do i=0,nfenetre
!    do j=0,nfenetre
     ! write (11,*) 'statneg',real(i)/500,real(j)*3/10-174,statneg(i,j)
!     enddo
  ! write (11,*)
!enddo

!do i=0,nfenetre
!    do j=0,nfenetre
   !   write (11,*) 'statpos',real(i)/500,real(j)*3/10-174,statpos(i,j)
!     enddo
  ! write (11,*)
!enddo

!kine=kine/(((time-tequilib)/dt)*Nbclones)

!write(*,*) 'ecintot',kine/(3/2*38)

open(unit=27, file=posfinal, status='replace')

!write (27,*) NbClones
!write (27,*) N
do m=0,NbClones-1
   write (27,*) P(m)%qx(:)
   write (27,*) P(m)%qy(:)
   write (27,*) P(m)%qz(:)
   write (27,*) 
   write (27,*) P(m)%px(:)
   write (27,*) P(m)%py(:)
   write (27,*) P(m)%pz(:)      
   write (27,*)
enddo 

!do m=0,NbClones-1
!   write (27,*) P(m)%uqx(:)
!   write (27,*) P(m)%uqy(:)
!   write (27,*) P(m)%uqz(:)
!   write (27,*) 
!   write (27,*) P(m)%upx(:)
!   write (27,*) P(m)%upy(:)
!   write (27,*) P(m)%upz(:)      
!   write (27,*)
!enddo

close(27)



do l=0,NbClones-1   

    do i=1,N-1
      do j=i+1,N
        xij(i,j)=(P(l)%qx(i-1)-P(l)%qx(j-1))
        yij(i,j)=(P(l)%qy(i-1)-P(l)%qy(j-1))
        zij(i,j)=(P(l)%qz(i-1)-P(l)%qz(j-1))
        dij(i,j)=sqrt(xij(i,j)**2+yij(i,j)**2+zij(i,j)**2)
        xij(j,i)=-xij(i,j)
        yij(j,i)=-yij(i,j)
        zij(j,i)=-zij(i,j)
        dij(j,i)= dij(i,j)
      enddo
    enddo

    xpq4(L)=q4(N,xij,yij,zij,dij)

call Force(P(l)%qx,P(l)%qy,P(l)%qz,fp,Utot)

!write(*,*)'energie',Time,xpq4(L),Utot

enddo



!write (*,*) 'atomi', P(0)%qx(:),P(0)%qy(:),P(0)%qz(:)
!write(*,*)'q4',xpq4(0)
!write (*,*) 'velocità', P(0)%px(:),P(0)%py(:),P(0)%pz(:)

!write (*,*)'enfin bref', Ecin4+utot+173.4, betaq*temperature

!!!!ALLA FINE FACCIO STATISTICA COMPLESSIVA DELLO SPAZIO DELLE FASI ESPLORATO DAI CLONI

!do i=-10,nfenetre+10
 !               do j=-10,nfenetre+10
  !                     x=cumul_contour_q4q6(i,j)
                        !if (x.gt.0) then  
                          !write(*,*) 'final q4 q6',i,j,x !cumul_contour_q4q6(i,j)
   !                       write(26,*)'final q4 q6',i,j,x,tot
                        !endif
    !         enddo
     !        write (26,*)
     !        enddo   

 do ifenetre=-10,nfenetre+10
                do ienergy=0,nfenetre
                        cumul_contour(ifenetre,ienergy)=cumul_contour(ifenetre,ienergy)
                        !if (cumul_contour(ifenetre,ienergy).gt.0) then  
                          !write(*,*) 'final q4 E',ifenetre,ienergy,cumul_contour(ifenetre,ienergy)
                         write(26,*) real(ifenetre)/500,real(ienergy)*3/10-174,cumul_contour(ifenetre,ienergy), tot
                        !endif
             enddo
              write (26,*)
            enddo   



 do ifen_q6=-10,nfenetre+10
                do ienergy=0,nfenetre
                      cumul_contour_q6(ifen_q6,ienergy)=cumul_contour_q6(ifen_q6,ienergy)
                        !if (cumul_contour_q6(ifen_q6,ienergy).gt.0) then  
                         !write(*,*) 'final q6 E',ifen_q6,ienergy,cumul_contour_q6(ifen_q6,ienergy)
                         write(14,*)real(ifen_q6)*3/500,real(ienergy)*3/10-174,cumul_contour_q6(ifen_q6,ienergy),tot
                      !endif
            enddo
               write (14,*)
            enddo

 do ifenetre=-10,nfenetre+10
                do ifen_q6=-10,nfenetre+10
                        cumul_contour_q4q6(ifenetre,ifen_q6)=cumul_contour_q4q6(ifenetre,ifen_q6)
                        !if (cumul_contour(ifenetre,ienergy).gt.0) then  
                          !write(*,*) 'final q4 E',ifenetre,ienergy,cumul_contour_q4q6(ifenetre,ienergy)
                         write(31,*) real(ifenetre)/500,real(ifen_q6)*3/500,cumul_contour_q4q6(ifenetre,ifen_q6), tot
                        !endif
             enddo
              write (31,*)
            enddo  
  
!write (*,*) 'usc', nbclones  
!do while(Time.gt.TimeStore)
    !  Lyapunov=0
     ! do l=0,NbClones-1
	!do j=0,N-1
	 ! write(11,*) Time, P(l)%qx(j),P(l)%qy(j), P(l)%qz(j), P(l)%px(j),P(l)%py(j), P(l)%pz(j), P(l)%uqx(j),P(l)%uqy(j),P(l)%uqz(j),P(l)%upx(j), P(l)%upy(j),P(l)%upz(j)
	!Lyapunov+=P(l).Lyap/(RealNumber) NbClones
     !   enddo
     ! enddo
      !write(sortie,'#Lyap\t%lg\t%lg\n',Time,P(4).Lyap/Time)
     ! write(sortie,*)'#Average\t%lg\t%lg\t%lg\n',Time,Lyapunov/Time,mu/Time
     
     ! TimeStore=TimeStore+inter
 !enddo
  
  !Lyapunov=0
  !for(l=0l.lt.NbClonesl++)
   ! Lyapunov+=P(l).Lyap/(RealNumber) NbClones
!  write(sortie,*)'#Final\t%lg\t%lg\t%lg\n'Time !,Lyapunov/Time,mu/Time)


stop

end subroutine neufdoublewell

Subroutine propagateA(qx,qy, qz, px,py,pz, uqx, uqy, uqz, upx,upy, upz,rga,temperature,dt, N,beq,dh,ecinetique,ff,lapU)

 implicit none
 integer :: i,alfa
  integer::N
  real (double), dimension (0:N-1) :: qx   ! vector of position
  real (double),dimension (0:N-1) :: px   ! vector of impulsion
  real (double),dimension (0:N-1) :: uqx  ! linear displacement
  real (double),dimension (0:N-1) :: upx  ! linear displacement
  real (double), dimension (0:N-1) :: qy   ! vector of position
  real (double),dimension (0:N-1) :: py   ! vector of impulsion
  real (double),dimension (0:N-1) :: uqy  ! linear displacement
  real (double),dimension (0:N-1) :: upy  ! linear displacement
  real (double), dimension (0:N-1) ::qz   ! vector of position
  real (double),dimension (0:N-1) :: pz   ! vector of impulsion
  real (double),dimension (0:N-1) :: uqz  ! linear displacement
  real (double),dimension (0:N-1) :: upz  ! linear displacement
  real (double) ::c1,a1,a2,b1,b2
  real (double) ::c2,Utot,lapU,tconf,ff
  real (double) ::dt
  real (double) ::gamma,rga,beq,fftot
  real (double) ::temperature,dh,ecinetique
  real (double) ::eta2 
  real (double),dimension (0:N-1) ::q1x 
  real (double),dimension (0:N-1) ::p1x
  real (double),dimension (0:N-1) ::q2x
  real (double),dimension (0:N-1) ::uq1x
  real (double),dimension (0:N-1) ::up1x
  real (double),dimension (0:N-1) ::uq2x
  real (double),dimension (0:N-1) ::q1y 
  real (double),dimension (0:N-1) ::p1y
  real (double),dimension (0:N-1) ::q2y
  real (double),dimension (0:N-1) ::uq1y
  real (double),dimension (0:N-1) ::up1y
  real (double),dimension (0:N-1) ::uq2y
  real (double),dimension (0:N-1) ::q1z 
  real (double),dimension (0:N-1) ::p1z
  real (double),dimension (0:N-1) ::q2z
  real (double),dimension (0:N-1) ::uq1z
  real (double),dimension (0:N-1) ::up1z
  real (double),dimension (0:N-1) ::uq2z
  real (double),dimension(0:N-1):: fxx
  real (double),dimension(0:N-1):: fzz
  real (double),dimension(0:N-1):: fyy
  real (double),dimension (3,N)::xp
  real (double),dimension (3,N)::x2p
  real (double),dimension (3,N)::vp
  real (double),dimension (3,N)::fp
  real (double),dimension (3,N)::f2p
  real(double), dimension (n)::vsecx,vsecy,vsecz,vsecx1,vsecy1, vsecz1

!write(*,*)'gamma', gamma  

!if ((gamma*dt).le.1000) then
 
!   rga = exp(-gamma*dt/two)

!  else 

!  rga=0

! endif
lapu=0
!write (*,*) 'rga', rga
call Hessien(qx, qy,qz, uqx, uqy, uqz, Vsecx, vsecy,vsecz,lapU)
!do i=1,n
!write (*,*) 'vsecx', vsecx(i), vsecy(i), vsecz(i)
!enddo
!write (*,*)'vsec', vsecx(3)

do i=0,N-1
  up1x(i)= upx(i)*rga-dt/2d0* Vsecx(i+1) 
  up1y(i)= upy(i)*rga-dt/2d0* Vsecy(i+1) 
  up1z(i)= upz(i)*rga-dt/2d0* Vsecz(i+1) 
enddo



do i=0,N-1
  uq1x(i)= uqx(i) + dt* up1x(i) 
  uq1y(i)= uqy(i) + dt* up1y(i)  
  uq1z(i)= uqz(i) + dt* up1z(i) 
enddo

!write(*,*) uq1x(3)

!call Force(qx,qy,qz,fp,Utot)

 xp(1,1:N)=qx(0:N-1)
 xp(2,1:N)=qy(0:N-1)
 xp(3,1:N)=qz(0:N-1)

 vp(1,1:N)=px
 vp(2,1:N)=py
 vp(3,1:N)=pz

!write(*,*) xp(1,4)


call calfoljc2(xp,fp,utot)

ff=sum(fp(1,:)**2+fp(2,:)**2+fp(3,:)**2)

!write (*,*) 'force',ff


!tconf=ff/lapU

!write(*,*) 'ctrl energ', utot

call langevin2(xp,vp,fp,dt,rga, temperature,beq,dh,ecinetique)
 

!write (*,*)'ecin', ecinetique 

qx(0:N-1)=xp(1,1:N)
qy(0:N-1)=xp(2,1:N)
qz(0:N-1)=xp(3,1:N)

px(0:N-1)=vp(1,1:N)
py(0:N-1)=vp(2,1:N)
pz(0:N-1)=vp(3,1:N) 

!write(*,*) qx(3)



!if (utot.lt.(-169.0)) then 
   alfa=-1
!else 
!   alfa=1
!endif

call Hessien(qx, qy,qz, uq1x, uq1y, uq1z, Vsecx1, vsecy1,vsecz1,lapU)

!write (*,*)  'rga',rga
!stop
!write (*,*)'vsec', vsecx1(3)
!
do i=0,N-1
   upx(i) = (up1x(i)-dt/2.0*Vsecx1(i+1))*rga 
   upy(i) = (up1y(i)-dt/2.0*Vsecy1(i+1))*rga
   upz(i) = (up1z(i)-dt/2.0*Vsecz1(i+1))*rga
enddo

do i=0,N-1
   uqx(i) = uq1x(i)
   uqy(i) = uq1y(i) 
   uqz(i) = uq1z(i) 
enddo

end subroutine propagateA


 Subroutine largedeviations (xp)

implicit none

real (double), dimension (:), allocatable :: Weight
integer, dimension (:),allocatable ::ToCopy
integer, dimension (:),allocatable ::ToDelete
integer, dimension (:),allocatable ::Tofollow
integer, dimension (:),allocatable ::Nb
integer, dimension (:),allocatable ::Number
!real(double), dimension (:),allocatable ::Number
!integer, parameter::N= 864
!integer, parameter ::N=38 ! numero atomes du cluster
!integer, parameter::N=256
real (double) ::r     !distance interatomique
real (double) ::eta    !paramètre d'énergie du potentiel
real (double) ::sigma ! dist. éq. pour le potentiel LJ 
real (double),dimension(0:N-1) :: fxx
real (double),dimension(0:N-1) :: fyy
real(double),dimension(0:N-1) :: fzz
real(double),dimension(0:N-1) :: vsecx
real (double),dimension(0:N-1) :: vsecy
real (double),dimension(0:N-1) :: vsecz

real(double),dimension (N,N)::   xij
real (double),dimension (N,N)::   yij
real(double),dimension (N,N)::   zij
real (double),dimension (N,N)::   dij

real (double),dimension (:),allocatable::w4,w6,wp4,wp6, xpq4,xpq6,tv,tvmoy
integer, dimension (:),allocatable ::Nvite
integer, dimension(:), allocatable ::multip

real (double),dimension (N):: fq4x,fq4y,fq4z,fq6x, fq6y,fq6z

real (double), dimension (3,N):: xp
real (double), dimension (3,N):: vp
real(double),dimension (3,N):: fp
real (double),dimension (3,N)::fppro
real (double), dimension (3,N):: xppro


real (double) :: q4,q6, q4moy, q6moy,w4moy, w6moy, lambdaq4,lambdaE, lambdaq6, logt, logt2,boh! linear displacement
type Particule
    real (double), dimension (0:N-1):: qx  
    real (double), dimension (0:N-1):: qy 
    real (double), dimension (0:N-1):: qz  ! vector of position
    
    real (double), dimension (0:N-1):: px  
    real (double), dimension (0:N-1):: py 
    real (double), dimension (0:N-1):: pz  ! vector of impulsion
    !RealNumber Lyap
 endtype

type (Particule), dimension (:), allocatable :: P

real (double) ::seed, utot, tot, x, U,x4,x6,tconf, tconftot,tconfmoy,kinetotdt, upro
 !     Declaration of the variables
real(double)::fftot, fftotdt,laputot,laputotdt,ff,lapu, pi2, potist,box_lenght
real (double) ::inter,AverageTime,NbTimeStep,dt,TotalTime,Time,w,Norm,TimeStore,c1,c2,sqrtgTh,eta1,eta2,gamma,eps, Temperature, lambda ! Lyapunov=0,LyapunovAv=0,mu=0,
integer :: NbClones,NMax,j,Ndelete,Ncopy,l,k,ran,m, a,i,b, iw4,iw6,lo,cont,ien,copiati,nfollow,nq4,teq,n0, n4

character(len=128) :: sortie
character(len=128) :: histotot
character(len=128) :: histopart
character(len=128) ::posfinal
character(len=128) :: snap
character(len=128) :: stock

real (double)::energy,rga,beq,rien,q,kine,kinetot,pif, trasm1,trasm2
real (double) ::xtemp0,xtemp1,xtemp2,xtemp3,xtemp4,xtemp5, taille
integer, parameter::nfenetre=100
integer :: ienergy, ifenetre, ifen_q6, iter,kaq4,kaq6,niter,SYS, liq, ico, nplace, num

real (double):: stat(0:nfenetre, 0:nfenetre)
real (double):: statneg(0:nfenetre, 0:nfenetre)
real (double):: statpos(0:nfenetre, 0:nfenetre)
real (double):: stat6(0:nfenetre, 0:nfenetre)
real(double)::tvmoytot,xn,enprmoy,jtrans, jtot,trasm,trastot

real(double) :: tab_contour(-10:nfenetre+10,0:nfenetre)
real(double) :: tab_contour_q6(-10:nfenetre+10,0:nfenetre)
real(double) :: tab_cont_q4q6(-10:nfenetre+10,-10:nfenetre+10)
real(double) :: cumul_contour(-10:nfenetre+10,0:nfenetre)
real(double) :: cumul_contour_q6(-10:nfenetre+10,0:nfenetre)
real(double) :: cumul_contour_q4q6(-10:nfenetre+10,-10:nfenetre+10)
real(double) :: tab_energy(-10:nfenetre+10,2),toto

real (double)::x111, xmaxq6,energie_min,energie_max,ecinetique,tequilib,delta,beqtot,dhtot,beqtotdt,dhtotdt,intbeq, intdh
integer::icfen,icener,icfen_q6,a1, a2, a3, a4, a5

real(double), dimension(:), allocatable::betaq,dh,ener,enerpro

write(6,*)'usage %s:\n'

!    Initialisation of the variables
read (*,*)
write(6,*)'TotalTime '
read (*,*) TotalTime      ! Total duration of the simulation 
write(*,*) TotalTime 
read (*,*)
write(6,*)'dt' 
read (*,*) dt          ! Time step
write(*,*) dt
read (*,*)
write(6,*)'NbClones'  
read (*,*) NbClones  
write(*,*) NbClones    ! Number of clones
read (*,*)
write(6,*)'temperature'  
read (*,*) temperature  ! temperature (kT)
write (*,*) temperature
read (*,*)
write(6,*)'gamma*dt'  
read (*,*) gamma  ! friction*dt
write (*,*) gamma       
read (*,*)
write(6,*)'sortie'  
read (*,*) sortie   
write (*,*) sortie! Name of the file where data are stored
read (*,*)
write(6,*)'lambda'
read (*,*) lambda
write (*,*) lambda
read (*,*)
write(6,*)'deltaT'
read (*,*) inter
write (*,*) inter
read (*,*)
write(6,*)'lo'
read (*,*) lo  
write (*,*) lo
read (*,*)        ! cloner (lo=1) ou pas (lo=0)
write(6,*)'cont'
read (*,*) cont
write (*,*) cont
read (*,*)
write(6,*)'niter'
read (*,*) niter
write (*,*) niter
read (*,*)
write(6,*)'t_equilib'
read (*,*) teq
write (*,*) teq
read (*,*)

open (unit=11, FILE=sortie, action='write', status='replace')
open (unit=31, action='write', status='replace')


histotot  =sortie(1:lenfnam)//'.thout'
histopart =sortie(1:lenfnam)// '.phout'
posfinal = sortie(1:lenfnam)//'.cin'
snap = sortie(1:lenfnam)//'.film'
stock = sortie(1:lenfnam)//'.stock'

open(unit=28, file=snap,  action='write', status='replace')
open(unit=14, file=histopart, action='write', status='replace')
open(unit=26, file=histotot,  action='write', status='replace')
open(unit=29, file=stock,  action='write', status='replace')

NMax= 100*NbClones      ! Max number of clones
!TimeStore = TotalTime-AverageTime ! Next time at which data will be stored 
                                           
! Array where the weight of every clones is stored
allocate (Weight(0:nmax-1)) 
allocate (w4(0:nmax))
allocate (w6(0:nmax))
allocate (wp4(0:nmax))
allocate (wp6(0:nmax))

allocate (tv(0:nmax-1))
allocate (tvmoy(0:nmax-1))
allocate (nvite(0:nmax-1))

allocate (xpq4(0:Nmax))
allocate (xpq6(0:Nmax))
allocate(multip(0:nmax))

allocate (betaq(0:nmax-1))
allocate (dh(0:nmax-1))
! Arrays where the clone to copy or to kill are stored
allocate(ToCopy(0:NMax-1))
allocate(ToDelete(0:NMax-1))
allocate(ener(0:nmax-1))
allocate(Tofollow(0:NMax-1))
allocate(enerpro(0:nmax-1))
 
! Array which contains the numbers from 1 to Nmax
allocate(Nb(0:NMax-1)) 

ToCopy(0:NMax-1)=0
ToDelete(0:NMax-1)=0
ener(:)=0
sigma=1
eta=1
enerpro(:)=0
xpq4(:)=0
xpq6(:)=0

n0=0

SYS=1
gamma=gamma/(dt)

write(*,*)'gamma' , gamma, gamma*dt
!temperature=0.19
stat(:,:)=0
stat6(:,:)=0
!statneg(:,:)=0
!statpos(:,:)=0
kine=0
laputot=0
fftot=0
laputotdt=0
fftotdt=0
kinetotdt=0
kinetot=0
nfollow=0
enprmoy=0
nq4=0
x6=0
x111   = 0.2
xmaxq6 = 0.6
energie_min = -174.0
energie_max = -144.0
box_lenght=6.35

     tab_contour(-10:nfenetre+10,0:nfenetre)=0
     tab_contour_q6(-10:nfenetre+10,0:nfenetre)=0
     tab_cont_q4q6(-10:nfenetre+10,-10:nfenetre+10)=0
     cumul_contour_q4q6(-10:nfenetre+10,-10:nfenetre+10)=0
     cumul_contour(-10:nfenetre+10,0:nfenetre)=0
     cumul_contour_q6(-10:nfenetre+10,0:nfenetre)=0

 do j=0,NMax-1
    Nb(j)=j   ! Number(j) is the number of descendants of clone j at the cloning step
 enddo

allocate(Number(0:NbClones-1)) 

! Some constants used in the quasisymplectic stochastic integrator (cfr Mannella)
! c1 = 1 - gamma * dt / 4
! c2 =  1 / ( 1 + gamma * dt / 4 )
! sqrtgTh = sqrt(2*gamma*temperature*dt)

beq=0
!ecinetique=0

!!!!!!!!!!!!T EQUILIBRAGE
tequilib=(dt)*teq
write(*,*) 'tquilib',tequilib
!!!!!!!!!!!! T EQUILIBRAGE

 write(*,*) 'sortie affichée'


allocate(P(0:NMax-1)) ! cas  sans lyapunov
close(27)
!open(unit=27, file=fnamcin, status='old')

n4=4*nbclones
!!!EQULIBRAGE
xppro(:,:)=0
Time=0
iter=0

rien=0
Q=0
logt=0
logt2=0
pif=0
pi2=0
trastot=1

tv(:)=0
tvmoy(:)=0
nvite(:)=0

if ((gamma*dt).le.100000) then
 
   rga = exp(-gamma*dt/2.d0)

else 

  rga=0

endif

write(*,*) 'largedeviations!'

write(*,*) 'rga',gamma, rga,2.d0

write (*,*) 'initialisation' 

if(cont.eq.1) then  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   inizio da vecchio file!!!!!!!!!!!!!!!!!!!!

open(unit=27, file=posfinal, status='old')

do m=0,NbClones-1
   read (27,*) P(m)%qx(:)
   read (27,*) P(m)%qy(:)
   read (27,*) P(m)%qz(:)
   read (27,*) 
   read (27,*) P(m)%px(:)
   read (27,*) P(m)%py(:)
   read (27,*) P(m)%pz(:)      
   read (27,*)
enddo 


close(27)

do l=0,NbClones-1   

    do i=1,N-1
      do j=i+1,N
        xij(i,j)=(P(l)%qx(i-1)-P(l)%qx(j-1))
        yij(i,j)=(P(l)%qy(i-1)-P(l)%qy(j-1))
        zij(i,j)=(P(l)%qz(i-1)-P(l)%qz(j-1))
        dij(i,j)=sqrt(xij(i,j)**2+yij(i,j)**2+zij(i,j)**2)
        xij(j,i)=-xij(i,j)
        yij(j,i)=-yij(i,j)
        zij(j,i)=-zij(i,j)
        dij(j,i)= dij(i,j)
      enddo
    enddo

    xpq4(L)=q4(N,xij,yij,zij,dij)

call Force(P(l)%qx,P(l)%qy,P(l)%qz,fp,Utot)

!write(*,*)'energie',Time,xpq4(L),Utot

enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! inizio da nuovo file
else  if (cont.eq.0) then  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! EQUILIBRAGGIO  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do while ((Time.lt.Tequilib).and.(time.lt.totaltime))
!write (*,*)'atomo', xp(1,1),vp(1,1)
 lapu=0
 ff=0
!call force_lj_MD(xp,fp,box_lenght,utot)

!call langevinlj_md(xp,vp,fp,dt,rga, temperature,Q,rien,ecinetique)
 ! call calfoljc2(xp, fp,utot)

 ! call langevin2(xp,vp,fp,dt,rga, temperature,Q,rien,ecinetique)

  ff=sum(fp(1,:)**2+fp(2,:)**2+fp(3,:)**2)
  fftot=fftot+ff

!call Hessien(xp(1,:), xp(2,:),xp(3,:), P(1)%uqx, P(1)%uqy, P(1)%uqz, Vsecx, vsecy,vsecz,lapU)

!write (*,*) 'lapu', lapu

!laputot=laputot+lapu


    do i=1,N-1
      do j=i+1,N
        xij(i,j)=xp(1,i)-xp(1,j)
        yij(i,j)=xp(2,i)-xp(2,j)
        zij(i,j)=xp(3,i)-xp(3,j)
        dij(i,j)=sqrt(xij(i,j)**2+yij(i,j)**2+zij(i,j)**2)
        xij(j,i)=-xij(i,j)
        yij(j,i)=-yij(i,j)
        zij(j,i)=-zij(i,j)
        dij(j,i)= dij(i,j)
      enddo
   enddo

x4=q4(N,xij,yij,zij,dij)

!call fq6(N,xij,yij,zij,dij,fq6x,fq6y,fq6z)
x6=q6(N,xij,yij,zij,dij)

write(*,'(e18.7, e18.7, e18.7, e18.7)') time, x4, x6,utot
!write (*,*) 'dissipation', rien, q

time=time+dt

enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FINE EQUILIBRAGGIO !!!!!!!!!!!!!!!!!!!!!!!!!

write(*,*)'tcin', time, ecinetique/(1.5d0*real(N)), temperature

     do  m = 0,NbClones-1  
       P(m)%qx=xp(1,:)
       P(m)%qy=xp(2,:)
       P(m)%qz=xp(3,:)
     enddo

   do  m = 0,NbClones-1
       P(m)%px=vp(1,1:N)
       P(m)%py=vp(2,1:N)
       P(m)%pz=vp(3,1:N)
    enddo

else  if (cont.eq.2) then 

    Num = 4!int((real(N)/4.0)**(1.0d0/3.0d0)+1.5)
      Nplace = 1!0
 
      taille =real(box_lenght)/real(Num)
 
      !Place = 0.2d0*Size
 !!!!!!!!!!!!!!!!!!! FCC CRYSTAL STRUCTURE !!!!!!!!!!!!!!! 
  !    Do I = 1,Num
  !       Do J = 1,Num
  !          Do K = 1,Num
  !              If (Nplace.Le.N) Then
   !               xp(1,Nplace) = (real(I) + 0.1+ 0.01d0*(genrand()-0.5d0))*taille
   !!               xp(2,Nplace) = (real(J) + 0.1 + 0.01d0*(genrand()-0.5d0))*taille
    !              xp(3,Nplace) = (real(K) + 0.1 +0.01d0*(genrand()-0.5d0))*taille
    !              Nplace = Nplace + 1
    !              xp(1,Nplace) = (real(I) + 0.1+ 0.01d0*(genrand()-0.5d0))*taille
    !              xp(2,Nplace) = (real(J) + 0.6 + 0.01d0*(genrand()-0.5d0))*taille
   !               xp(3,Nplace) = (real(K) + 0.6 +0.01d0*(genrand()-0.5d0))*taille
   !               Nplace = Nplace + 1
   !               xp(1,Nplace) = (real(I) + 0.6+ 0.01d0*(genrand()-0.5d0))*taille
    !!              xp(2,Nplace) = (real(J) + 0.6 + 0.01d0*(genrand()-0.5d0))*taille
    !              xp(3,Nplace) = (real(K) + 0.1 +0.01d0*(genrand()-0.5d0))*taille
    !              Nplace = Nplace + 1
     !             xp(1,Nplace) = (real(I) + 0.6+ 0.01d0*(genrand()-0.5d0))*taille
     !!             xp(2,Nplace) = (real(J) + 0.1 + 0.01d0*(genrand()-0.5d0))*taille
     !             xp(3,Nplace) = (real(K) + 0.6 +0.01d0*(genrand()-0.5d0))*taille
      !            Nplace = Nplace + 1
   !            Endif
   !         Enddo
   !      Enddo
    !  Enddo
write(*,*) N, num, nplace, num**3, taille
write(*,*) 1, real(1)

     Do I = 0,Num-1
         Do J = 0,Num-1
            Do K = 0,Num-1
                If (Nplace.Le.N) Then
                  xp(1,Nplace) = (real(I) + 0.1d0)*taille
                  xp(2,Nplace) = (real(J) + 0.1d0)*taille
                  xp(3,Nplace) = (real(K) + 0.1d0)*taille
                  !Nplace = Nplace + 1
                  xp(1,Nplace+1) = (real(I) + 0.1d0)*taille
                  xp(2,Nplace+1) = (real(J) + 0.6d0)*taille
                  xp(3,Nplace+1) = (real(K) + 0.6d0)*taille
                  !Nplace = Nplace + 1
                  xp(1,Nplace+2) = (real(I) + 0.6d0)*taille
                  xp(2,Nplace+2) = (real(J) + 0.6d0)*taille
                  xp(3,Nplace+2) = (real(K) + 0.1d0)*taille
                 ! Nplace = Nplace + 1
                  xp(1,Nplace+3) = (real(I) + 0.6d0)*taille
                  xp(2,Nplace+3) = (real(J) + 0.1d0)*taille
                  xp(3,Nplace+3) = (real(K) + 0.6d0)*taille
                  Nplace = Nplace + 4!1
               Endif
            Enddo
         Enddo
      Enddo


write(*,*)'??x',nplace !xp(1,:)
write(*,*) xp(1,N)
vp(:,:)=0


!do i=1,3
! do j=1,N
!  xp(i,j)=(genrand()-0.5)*box_lenght
!  vp(i,j)=0.0!genrand()!
! enddo 
!enddo
do while ((Time.lt.Tequilib).and.(time.lt.totaltime))

!write(*,*)'??v', vp(1,:)

!call force_lj_MD(xp,fp,box_lenght,utot)
!write(*,*)'??f', fp(1,:)


!call langevinlj_md(xp,vp,fp,dt,rga, temperature,Q,rien,ecinetique)
!do j=0,n
! write(*,*) xp(1,j), xp(2,j), xp(3,j)
!enddo
!write(*,*)'??v', vp(1,250)

 do i=1,N-1
      do j=i+1,N
        xij(i,j)=xp(1,i)-xp(1,j)
        yij(i,j)=xp(2,i)-xp(2,j)
        zij(i,j)=xp(3,i)-xp(3,j)
        dij(i,j)=sqrt(xij(i,j)**2+yij(i,j)**2+zij(i,j)**2)
        xij(j,i)=-xij(i,j)
        yij(j,i)=-yij(i,j)
        zij(j,i)=-zij(i,j)
        dij(j,i)= dij(i,j)
      enddo
   enddo

   !   write (*,*) 'xij', xij(:,:)
x4=q4(N,xij,yij,zij,dij)

!stop
!call fq6(N,xij,yij,zij,dij,fq6x,fq6y,fq6z)
x6=q6(N,xij,yij,zij,dij)
!
!call energy_lj_MD(xp,fp,box_lenght,utot)
call force_lj_MD(xp,fp,box_lenght,utot)
!write(*,*)'??f', fp(1,:)
call langevinlj_md(xp,vp,fp,dt,rga, temperature,Q,rien,ecinetique)

!call force_lj_MD(xp,fp,box_lenght,utot)
!write(*,'(e18.7, e18.7, e18.7, e18.7)') time, x4, x6, utot/real(n)
!write (*,*) temperature
time=time+dt

enddo

! write(*,*)'??xv',xp(1,:), vp(1,1:N)


     do  m = 0,NbClones-1  
       P(m)%qx=xp(1,:)
       P(m)%qy=xp(2,:)
       P(m)%qz=xp(3,:)
     enddo

   do  m = 0,NbClones-1
       P(m)%px=vp(1,1:N)
       P(m)%py=vp(2,1:N)
       P(m)%pz=vp(3,1:N)
    enddo


endif
 
write(*,*) 'initialisation faite' 
!do  m = 0,NbClones-1  
!       write(*,*)'??',P(m)%qx
!      write(*,*)'??',P(m)%px
!     enddo
!stop


!fftotdt=fftot
!laputotdt=laputot
betaq(:)=q
dh(:)=0
kinetot=ecinetique/(1.5d0*real(N))
ecinetique=0
fftot=0
laputot=0
fftotdt=0!fftot
laputotdt=0!laputot
box_lenght=6.5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  INIZIO ALGORITMO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do while (Time.lt.TotalTime)

    Ncopy=0
    Ndelete=0
    copiati=0
    n0=0
    w=0
    u=0
    q4moy=0
    q6moy=0
    tconf=0
    beqtot=0
    dhtot=0
    kine=0
    jtrans=0
    jtot=0
    trasm1=0
    ico=0
    liq=0
    beqtotdt=0
dhtotdt=0
    !trasm2=0   
!write(*,*) 'te', temperature
    iter=int(time/dt)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  PROPAGAZIONE DELEL PARTICELLE E CALCOLO TEMPERATURA  !!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do j=0,NbClones-1
     
       call propagLD(P(j)%qx,P(j)%qy,P(j)%qz,P(j)%px,P(j)%py,P(j)%pz,rga,temperature,dt,N,betaq(j),dH(j),ecinetique,ff,lapu,utot,cont,box_lenght)

       kine=kine+ecinetique
       fftot=fftot+ff
      ! write (*,*) 'ff', ff
       !laputot=laputot+lapu
       dhtot=dhtot+dh(j)/nbclones
       beqtot=beqtot+betaq(j)/nbclones
       !tconftot=tconftot+tconf   
       !write (*,*) 'dissipation', ecinetique!, dH(j), betaq(j)
       !call propagateC(P(j)%qx,P(j)%qy,P(j)%qz,P(j)%uqx,P(j)%uqy,P(j)%uqz,temperature,dt,N) !!!! senza inerzia
    enddo


fftotdt=fftotdt+fftot/real(Nbclones)
!laputotdt=laputotdt+laputot/real(Nbclones)
!write(*,*) 'cin', time, (kine/real(nbclones))/(1.5d0*real(N))
!write(*,*) (kine/real(nbclones))/((3/2)*N)
kine=(kine/real(nbclones))/(1.5d0*real(N))


   kinetot=kinetot+kine

!write(*,*) 'ctrl',ff! time,  kine, kinetot

beqtotdt=beqtotdt+beqtot
dhtotdt=dhtotdt+dhtot
intbeq=(beqtotdt)/((time-tequilib)/dt+1.0)
intdh=(dhtotdt)/((time-tequilib)/dt+1.0)

   ! write (*,*) 'dissipation', intdh,intbeq
!write(*,*) 'ctrl', kine, kinetot

!tconfmoy=fftotdt/laputotdt
 
kinetotdt=(kinetot)/((time-tequilib)/dt+1.0)
!write(*,*) 'ctrl', kine, kinetot
if(mod(iter,10000).eq.0) then
write(*, '(a, e18.7, e18.7, e18.7, e18.7)')'tcin',time, kinetotdt, temperature!0.5*tconfmoy
endif


!!!!!!!!!!!!!!!!! ANALISI REPLICHE DEL SISTEMA (CON O SENZA CLONING) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do l=0,NbClones-1   

    do i=1,N-1
      do j=i+1,N
        xij(i,j)=(P(l)%qx(i-1)-P(l)%qx(j-1))
        yij(i,j)=(P(l)%qy(i-1)-P(l)%qy(j-1))
        zij(i,j)=(P(l)%qz(i-1)-P(l)%qz(j-1))
        dij(i,j)=sqrt(xij(i,j)**2+yij(i,j)**2+zij(i,j)**2)
        xij(j,i)=-xij(i,j)
        yij(j,i)=-yij(i,j)
        zij(j,i)=-zij(i,j)
        dij(j,i)= dij(i,j)
      enddo
    enddo

    xpq4(L)=q4(N,xij,yij,zij,dij)

  
    !call fq4(N,xij,yij,zij,dij,fq4x,fq4y,fq4z)
  
    xpq6(L)=q6(N,xij,yij,zij,dij)

   if (xpq6(l).lt.0.1) then

      liq=liq+1

   else if (xpq6(l).gt.0.1) then

      ico=ico+1

   endif 

       !call fq6(N,xij,yij,zij,dij,fq6x,fq6y,fq6z)

       call energy_lj_MD(P(l)%qx,P(l)%qy,P(l)%qz,fp,Utot,box_lenght)

       ener(l)=utot

!write (*,*) '??', xpq4(l), xpq6(l), ener(l)

    q4moy=q4moy+xpq4(l)
    q6moy=q6moy+xpq6(l)
    u=u+utot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ISTOGRAMMA
if (lo.eq.5) then
   ifenetre = int((xpq4(l)/x111)*real(nfenetre))
   ifen_q6 = int((xpq6(l)/xmaxq6)*real(nfenetre))
   ienergy  = int((ener(L)-energie_min)/(energie_max-energie_min)*real(nfenetre))
   tab_contour(-10:nfenetre+10,0:nfenetre)=0
   tab_contour_q6(-10:nfenetre+10,0:nfenetre)=0
   tab_cont_q4q6(-10:nfenetre+10,-10:nfenetre+10)=0
  
   !!!!!CALPOT_AUXI
   lambdaq4 = (xpq4(L)/x111)*real(nfenetre)-real(ifenetre)

   if (ifenetre.gt.nfenetre+10) ifenetre = nfenetre+10
   icfen = ifenetre+1 
   if (icfen.gt.nfenetre+10) icfen = nfenetre+10
   if (ifenetre.lt.-10) ifenetre = -10

   lambdaE = ((potist-energie_min)/(energie_max-energie_min))*real(nfenetre)-real(ienergy)

   if (ienergy.gt.nfenetre) ienergy = nfenetre
   if (ienergy.lt.0) ienergy = 0
   icener = ienergy+1 
   if (icener.gt.nfenetre) icener = nfenetre

   lambdaq6 = (xpq6(L)/xmaxq6)*real(nfenetre)-real(ifen_q6)

   if (ifen_q6.gt.nfenetre+10) ifen_q6 = nfenetre+10
   icfen_q6 = ifen_q6+1 
   if (icfen.gt.nfenetre+10) icfen_q6 = nfenetre+10
   if (ifen_q6.lt.-10) ifen_q6 = -10

 !!!!!!!!!!!!!!!!!!!!!!!!!!!FINE calpot auxi (???)

  !!!!histogramme q4 - energie

        if ((ifenetre.le.nfenetre+10).and.(ifenetre.ge.-10)) then
           if ((ienergy.le.nfenetre).and.(ifenetre.ge.0)) then
             tab_contour(ifenetre,ienergy) = tab_contour(ifenetre,ienergy) +1.d0
             cumul_contour(ifenetre,ienergy)=1+cumul_contour(ifenetre,ienergy)
           endif
        endif

   !histogramme q6-energie
        
        if ((ifen_q6.le.nfenetre+10).and.(ifen_q6.ge.-10)) then
            if ((ienergy.le.nfenetre).and.(ienergy.ge.0)) then
               tab_contour_q6(ifen_q6,ienergy)=1+tab_contour_q6(ifen_q6,ienergy)
               cumul_contour_q6(ifen_q6,ienergy)=1+cumul_contour_q6(ifen_q6,ienergy)
            endif
        endif


   !histogramme q4-q6 (fait par moi)
          if ((ifenetre.le.nfenetre+10).and.(ifenetre.ge.-10)) then
            if ((ifen_q6.le.nfenetre+10).and.(ifen_q6.ge.-10)) then
                tab_cont_q4q6(ifenetre,ifen_q6)=tab_cont_q4q6(ifenetre,ifen_q6)+1
                cumul_contour_q4q6(ifenetre,ifen_q6)=1+cumul_contour_q4q6(ifenetre,ifen_q6)
          endif 
         endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  fine istogramma
endif
!!!!!!!!!!!!!!!!!!!!!!!!!
enddo

!!!!!!!!!!!!!!!!!!!!!! fine analisi repliche

!lo=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    INIZIO CLONAGGIO SE l=1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (lo.eq.1) then
 
   do j=0,NbClones-1
     ! Lyap computation
     ! P(j).Lyap  += log(Norm)
  ! Compute the weight for the cloning step
       Weight(j) = exp(-lambda*(xpq6(j)))
!  w=w+ Weight(j)
   enddo

!w=w/real(NbClones)
   
do j=0,NbClones-1

 !     Weight(j)=Weight(j)/w
      !How many offsprings for the clone j ?
      !Number(j)=floor(Weight(j)-0.5+genrand())
      Number(j)=floor(Weight(j)+genrand())
        !write (*,*) 'norme intiere de j' ,j, number(j)
      
      !if 0 -> delete it
      if (Number(j).lt.1) then

	    ToDelete(Ndelete)=j 
            Ndelete=Ndelete+1
            ien  = int((ener(j)-energie_min)/(energie_max-energie_min)*real(nfenetre))
            kaq4 = int((xpq4(j)/x111)*nfenetre)
            kaq6 = int((xpq6(j)/xmaxq6)*nfenetre)

            ! write (*,*)'detruit clone # avec norme', j, number(j)!, Ndelete,ToDelete(Ndelete-1) 
             if(Ndelete.gt.NbClones) then
	         write(6,*)'All the clones have been killed'
	         stop
             endif
          
           endif

      ! if strictly more than one offspring -> cloning
      
     do while (Number(j).gt.1) 
             	
        Ncopy=Ncopy+1
        ToCopy(Ncopy-1)=j
         kaq4 = int((xpq4(j)/x111)*nfenetre)
         kaq6 = int((xpq6(j)/xmaxq6)*nfenetre)
         ien  = int((ener(j)-energie_min)/(energie_max-energie_min)*real(nfenetre))
        
         ! write (*,*) 'clonage clone # avec norme',j, number(j), ncopy, tocopy(ncopy-1)
	Number(j)=Number(j)-1
	if((Ncopy.gt.NMax).or.(Ncopy.eq.NMax)) then
	   write(6,*)'Too many clones'
	    stop
        endif
     enddo

 
 enddo

!  write(*,*) 'statistica', time, Ndelete, Ncopy, n0, nbclones

     ! At this stage, we have a list of N1 clones to delete in
      !ToDelete and a list of N2 clones to copy in ToCopy. If N1=N2,
      !then the clones to be copied are copied where those who have
      !to be deleted are stored.
     
 !     If N1.gt.N2, there are more clones to delete than to
 !     clone. First, we copy the N2 clones to be copied on N2 clones
 !     to be deleted. We are left with N1-N2 clones still to be
 !     deleted. The idea is to pull at random N1-N2 clones _among the
  !    rest_ and to put them in place of the the N1-N2 still to be deleted.
    
 !     If N2.gt.N1, there are more clones to be copied that to be
 !     killed. We first copy N1 clones to be copied on the N1 clones
 !     to be deleted. Then, one chooses N2-N1 clones at random _among
 !     the NbClones+N2-N1 which we have and delete them to go back to
 !     NbClones clones.
  
  !    we copy what's possible to copy
!if (lo.eq.1) then
do j=0,MIN(Ndelete-1,Ncopy-1)
      ! write (*,*) 'copie',j, ndelete, todelete(j), tocopy(j)
      do l=0,N-1
        P(ToDelete(j))%qx(l) =P(ToCopy(j))%qx(l) 
	P(ToDelete(j))%px(l) =P(ToCopy(j))%px(l) 
	  
        P(ToDelete(j))%qy(l) =P(ToCopy(j))%qy(l) 
	P(ToDelete(j))%py(l) =P(ToCopy(j))%py(l) 
	
        P(ToDelete(j))%qz(l) =P(ToCopy(j))%qz(l) 
	P(ToDelete(j))%pz(l) =P(ToCopy(j))%pz(l) 
	
        enddo
        betaq(ToDelete(j))=betaq(ToCopy(j))
	dh(ToDelete(j))=dh(ToCopy(j))
     enddo
     !  P(ToDelete(j)).Lyap =P(ToCopy(j)).Lyap 
 
!!!!!!!!!!!!SYSTR!!!!!!!!!!!!!!
!if(sys.eq.0) then   
!!!!!!!!!!SYSTR!!!!!!!!!!!!

    
    if(Ncopy.gt.Ndelete) then
 ! We put the clones still to be copied at the end of the
 !  array of clones.
      do j=Ndelete,Ncopy-1
     !  write(*,*) ndelete, ncopy, j, NbClones+j-Ndelete
	do l=0,N-1
	  P(NbClones+j-Ndelete)%qx(l) =P(ToCopy(j))%qx(l) 
	  P(NbClones+j-Ndelete)%px(l) =P(ToCopy(j))%px(l) 
	  
          P(NbClones+j-Ndelete)%qy(l) =P(ToCopy(j))%qy(l) 
	  P(NbClones+j-Ndelete)%py(l) =P(ToCopy(j))%py(l) 
	  
          P(NbClones+j-Ndelete)%qz(l) =P(ToCopy(j))%qz(l) 
	  P(NbClones+j-Ndelete)%pz(l) =P(ToCopy(j))%pz(l) 
	 
         enddo
         betaq(NbClones+j-Ndelete)=betaq(ToCopy(j))
	 dh(NbClones+j-Ndelete)=dh(ToCopy(j))

	!P(NbClones+j-Ndelete).Lyap =P(ToCopy(j)).Lyap 
      enddo
     
      ! We now have to delete Ncopy - Ndelete clones to keep the
      ! population constant
      
      Ndelete=Ncopy-Ndelete
      
   !   We pull at random Ndelete among the clones
   !   and store them in ToDelete
!write (*,*)'combinaison',Ndelete,NbClones+Ndelete
     call Combinaison(Ndelete,NbClones+Ndelete,ToDelete,Nb)
   ! In 'todelete' we now have Ndelete different rank chosen at random.

!   we reorder Nb() for a future use

      do j=0,Ndelete-1
       !write (*,*) j, ndelete,ToDelete(j)
	Nb(ToDelete(j))=ToDelete(j)
      enddo    
  
    !   Sort the array ToDelete in ascending numerical order
      call shell(Ndelete,ToDelete-1)
      
     !  We delete the clones, starting by the last one.
  !    To do so, we replace them by the last clone of the array
  !    of particles
      do j=Ndelete-1,0,-1
         
	 do l=0,N-1
	  P(ToDelete(j))%qx(l) =P(NbClones+Ndelete-1)%qx(l) 
	  P(ToDelete(j))%px(l) =P(NbClones+Ndelete-1)%px(l) 
	
          P(ToDelete(j))%qy(l) =P(NbClones+Ndelete-1)%qy(l) 
	  P(ToDelete(j))%py(l) =P(NbClones+Ndelete-1)%py(l) 
	  
          P(ToDelete(j))%qz(l) =P(NbClones+Ndelete-1)%qz(l) 
	  P(ToDelete(j))%pz(l) =P(NbClones+Ndelete-1)%pz(l) 
	  
	enddo
        betaq(ToDelete(j))=betaq(NbClones+Ndelete-1)
	 dh(ToDelete(j))=dh(NbClones+Ndelete-1)
        
         	!P(ToDelete(j)).Lyap=P(NbClones+Ndelete-1).Lyap
	Ndelete=ndelete-1
      enddo

 endif

    if(Ncopy.lt.Ndelete) then

      l=1

  !  We put the dead guys at the end

      do j=Ndelete-1,Ncopy,-1
       ! write(*,*)  j, ToDelete(j),NbClones-l
        do k=0,N-1
	  P(ToDelete(j))%qx(k) =P(NbClones-l)%qx(k)  
	  P(ToDelete(j))%px(k) =P(NbClones-l)%px(k)  
	
          P(ToDelete(j))%qy(k) =P(NbClones-l)%qy(k)  
	  P(ToDelete(j))%py(k) =P(NbClones-l)%py(k)  
	 
          P(ToDelete(j))%qz(k) =P(NbClones-l)%qz(k)  
	  P(ToDelete(j))%pz(k) =P(NbClones-l)%pz(k)  

        enddo
         tv(ToDelete(j))=tv(NbClones-l)
         tvmoy(ToDelete(j))=tv(NbClones-l)
        betaq(ToDelete(j))=betaq(NbClones-l)
	 dh(ToDelete(j))=dh(NbClones-l)
	l=l+1
      enddo
      
   !   We have to clone Ndelete-Ncopy clone to keep the population constant
      
      Ncopy=Ndelete-Ncopy
      l=1
      do m=0,Ncopy-1
	ran = floor(genrand()*(NbClones-Ncopy))
!write(*,*) NbClones-l,ran
	do k=0,N-1
	  P(NbClones-l)%qx(k) =P(ran)%qx(k)  
	  P(NbClones-l)%px(k) =P(ran)%px(k)  
	  
          P(NbClones-l)%qy(k) =P(ran)%qy(k)  
	  P(NbClones-l)%py(k) =P(ran)%py(k)  
	
          P(NbClones-l)%qz(k) =P(ran)%qz(k)  
	  P(NbClones-l)%pz(k) =P(ran)%pz(k)  
	
	enddo
         ! tv(NbClones-l)=tv(ran)
        ! tvmoy(NbClones-l)=tv(ran)
        betaq(NbClones-l)=betaq(ran)
	 dh(NbClones-l)=dh(ran)
	!P(NbClones-l).Lyap=P(ran).Lyap  
	l=l+1
      enddo
   endif


ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  FINE CLONAGGIO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


q4moy=q4moy/NbClones
q6moy=q6moy/NbClones
u=u/real(NbClones)
nq4=0

tot=(NbClones)*(TotalTime/dt) !renormalisation histog q4 q6


if(mod(iter,10000).eq.0) then
 write (*,'(a, e18.7, e18.7, e18.7, e18.7, e18.7, e18.7)')'param', time, q4moy ,q6moy,U/real(n), real(ico), real(liq)!,temperature
 !11write (*,*) 'trasm',  ico, liq 
endif

if(mod(iter,1000).eq.0) then
do i=0, Nbclones-1
write(28,'(e11.4, e11.4, e11.4, e12.4, e12.4)') time ,xpq4(i),xpq6(i),ener(i)/real(n), temperature!,w4(i),w6(i),enerpro(i), weight(i),w
enddo
endif




!boh=real(temperature)-real(inter)
!!!!!!!!!!!!!!!!!!!!!! ANNEALING!!!!!!!!!!!!!!!!!!!!!!
if ((time.gt.1500).and.(temperature.gt.0.5)) then
    temperature=temperature-inter!temperature= real(temperature)-real(inter)
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

!!!!!!!!!!!!!!!!!!!!!! quenching!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !if (time.gt.500) then
!if(mod(iter,10000).eq.0) then
  !  temperature=0.03  !lo=1
!gamma=100.d0
! write (*,*) 'temp',time, temperature, kinetotdt!, !temperature-inter
! endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!! CLONING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!if ((time.gt.28000).and.(temperature.le.0.04)) then
if (time.gt.1500) then
 !lo=1
endif
!  mu += log(weight/(RealNumber) NbClones)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


if(mod(iter,5000).eq.0) then
write(29,*) time
do m=0,NbClones-1
!write(29,*) m
   write (29,*) P(m)%qx(:)
   write (29,*) P(m)%qy(:)
   write (29,*) P(m)%qz(:)
   write (29,*) 
   write (29,*) P(m)%px(:)
   write (29,*) P(m)%py(:)
   write (29,*) P(m)%pz(:)      
   write (29,*)
enddo 
endif


!!!TEMPS FINI.RIPARTE IL CICLO
Time=Time+dt

enddo



write(29,*) time
do m=0,NbClones-1
!write(29,*) m
   write (29,*) P(m)%qx(:)
   write (29,*) P(m)%qy(:)
   write (29,*) P(m)%qz(:)
   write (29,*) 
   write (29,*) P(m)%px(:)
   write (29,*) P(m)%py(:)
   write (29,*) P(m)%pz(:)      
   write (29,*)
enddo
close (29)

 write (*,*) 'finito'


!!!!ALLA FINE FACCIO STATISTICA COMPLESSIVA DELLO SPAZIO DELLE FASI ESPLORATO DAI CLONI

!do i=-10,nfenetre+10
 !               do j=-10,nfenetre+10
  !                     x=cumul_contour_q4q6(i,j)
                        !if (x.gt.0) then  
                          !write(*,*) 'final q4 q6',i,j,x !cumul_contour_q4q6(i,j)
   !                       write(26,*)'final q4 q6',i,j,x,tot
                        !endif
    !         enddo
     !        write (26,*)
     !        enddo   
if (lo.eq.5) then
 do ifenetre=-10,nfenetre+10
                do ienergy=0,nfenetre
                        cumul_contour(ifenetre,ienergy)=cumul_contour(ifenetre,ienergy)
                        !if (cumul_contour(ifenetre,ienergy).gt.0) then  
                          !write(*,*) 'final q4 E',ifenetre,ienergy,cumul_contour(ifenetre,ienergy)
                         write(26,*) real(ifenetre)/500,real(ienergy)*3/10-174,cumul_contour(ifenetre,ienergy), tot
                        !endif
             enddo
              write (26,*)
            enddo   



 do ifen_q6=-10,nfenetre+10
                do ienergy=0,nfenetre
                      cumul_contour_q6(ifen_q6,ienergy)=cumul_contour_q6(ifen_q6,ienergy)/tot
                        !if (cumul_contour_q6(ifen_q6,ienergy).gt.0) then  
                         !write(*,*) 'final q6 E',ifen_q6,ienergy,cumul_contour_q6(ifen_q6,ienergy)
                         write(14,*)real(ifen_q6)*3/500,real(ienergy)*3/10-174,cumul_contour_q6(ifen_q6,ienergy),tot
                      !endif
            enddo
               write (14,*)
            enddo
endif  
!write (*,*) 'usc', nbclones  
!do while(Time.gt.TimeStore)
    !  Lyapunov=0
     ! do l=0,NbClones-1
	!do j=0,N-1
	 ! write(11,*) Time, P(l)%qx(j),P(l)%qy(j), P(l)%qz(j), P(l)%px(j),P(l)%py(j), P(l)%pz(j), P(l)%uqx(j),P(l)%uqy(j),P(l)%uqz(j),P(l)%upx(j), P(l)%upy(j),P(l)%upz(j)
	!Lyapunov+=P(l).Lyap/(RealNumber) NbClones
     !   enddo
     ! enddo
      !write(sortie,'#Lyap\t%lg\t%lg\n',Time,P(4).Lyap/Time)
     ! write(sortie,*)'#Average\t%lg\t%lg\t%lg\n',Time,Lyapunov/Time,mu/Time
     
     ! TimeStore=TimeStore+inter
 !enddo
  
  !Lyapunov=0
  !for(l=0l.lt.NbClonesl++)
   ! Lyapunov+=P(l).Lyap/(RealNumber) NbClones
!  write(sortie,*)'#Final\t%lg\t%lg\t%lg\n'Time !,Lyapunov/Time,mu/Time)


stop

end subroutine largedeviations


Subroutine propagLD(qx,qy, qz, px,py,pz, rga,temperature,dt, N,beq,dh,ecinetique,ff,lapU,utot,cont,box)

 implicit none
 integer :: i,alfa,cont
  integer::N
  real (double), dimension (0:N-1) :: qx   ! vector of position
  real (double),dimension (0:N-1) :: px   ! vector of impulsion
  
  real (double), dimension (0:N-1) :: qy   ! vector of position
  real (double),dimension (0:N-1) :: py   ! vector of impulsion
  
  real (double), dimension (0:N-1) ::qz   ! vector of position
  real (double),dimension (0:N-1) :: pz   ! vector of impulsion
  
  real (double) ::c1,a1,a2,b1,b2
  real (double) ::c2,Utot,lapU,tconf,ff
  real (double) ::dt, box
  real (double) ::gamma,rga,beq,fftot
  real (double) ::temperature,dh,ecinetique
  real (double) ::eta2 
  real (double),dimension (0:N-1) ::q1x 
  real (double),dimension (0:N-1) ::p1x
  real (double),dimension (0:N-1) ::q2x
 

  real (double),dimension (0:N-1) ::q1y 
  real (double),dimension (0:N-1) ::p1y
  real (double),dimension (0:N-1) ::q2y

  real (double),dimension (0:N-1) ::q1z 
  real (double),dimension (0:N-1) ::p1z
  real (double),dimension (0:N-1) ::q2z

  real (double),dimension (0:N-1) ::uqx
  real (double),dimension (0:N-1) ::uqy
  real (double),dimension (0:N-1) ::uqz

  real (double),dimension(0:N-1):: fxx
  real (double),dimension(0:N-1):: fzz
  real (double),dimension(0:N-1):: fyy
  real (double),dimension (3,N)::xp
  real (double),dimension (3,N)::x2p
  real (double),dimension (3,N)::vp
  real (double),dimension (3,N)::fp
  real (double),dimension (3,N)::f2p

uqx(:)=0
uqy(:)=0
uqz(:)=0

lapu=0
!box=10.0
!call Hessien(qx, qy,qz, uqx, uqy, uqz, Vsecx, vsecy,vsecz,lapU)

 xp(1,1:N)=qx(0:N-1)
 xp(2,1:N)=qy(0:N-1)
 xp(3,1:N)=qz(0:N-1)

 vp(1,1:N)=px
 vp(2,1:N)=py
 vp(3,1:N)=pz

!if (cont.EQ.2) then
  call force_lj_MD(xp,fp,box,utot)
!write(*,*) 'ctrl', utot/real(256)
  call langevinlj_md(xp,vp,fp,dt,rga, temperature,beq,dh,ecinetique)
!else
!  call calfoljc2(xp,fp,utot)
!  call langevin2(xp,vp,fp,dt,rga, temperature,beq,dh,ecinetique)
!endif

ff=sum(fp(1,:)**2+fp(2,:)**2+fp(3,:)**2)

if (ff.gt.10000) then
!write(*,*) 'cazz', ff! write (*,*) fp(:,:)
!stop
endif!write(*,*) 'f', ff

ff=0.0d0 

qx(0:N-1)=xp(1,1:N)
qy(0:N-1)=xp(2,1:N)
qz(0:N-1)=xp(3,1:N)

px(0:N-1)=vp(1,1:N)
py(0:N-1)=vp(2,1:N)
pz(0:N-1)=vp(3,1:N) 


   alfa=-1

end subroutine propagLD



 Subroutine LWDTPS (xp)

implicit none

real (double), dimension (:), allocatable :: Weight
integer, dimension (:),allocatable ::ToCopy
integer, dimension (:),allocatable ::ToDelete
integer, dimension (:),allocatable ::Tofollow
integer, dimension (:),allocatable ::Nb
integer, dimension (:),allocatable ::Number
!integer, parameter ::N=38 ! numero atomes du cluster
 
real (double) ::r     !distance interatomique
real (double) ::eta    !paramètre d'énergie du potentiel
real (double) ::sigma ! dist. éq. pour le potentiel LJ 
real (double),dimension(0:N-1) :: fxx
real (double),dimension(0:N-1) :: fyy
real(double),dimension(0:N-1) :: fzz
real(double),dimension(0:N-1) :: vsecx
real (double),dimension(0:N-1) :: vsecy
real (double),dimension(0:N-1) :: vsecz

real(double),dimension (N,N)::   xij
real (double),dimension (N,N)::   yij
real(double),dimension (N,N)::   zij
real (double),dimension (N,N)::   dij

real (double),dimension (:),allocatable::w4,w6,wp4,wp6
real (double),dimension (:,:),allocatable::xpq4,xpq6,ener
integer, dimension (:),allocatable ::Nvite
integer, dimension(:), allocatable ::multip

real (double),dimension (N):: fq4x,fq4y,fq4z,fq6x, fq6y,fq6z

real (double), dimension (3,N):: xp
real (double), dimension (3,N):: vp
real (double), dimension (3,N):: fp
real (double), dimension (3,N)::fppro
real (double), dimension (3,N):: xppro
real (double), dimension(1:N):: delt
real (double), dimension(1:N):: gQ4x
real (double), dimension(1:N):: gQ4y
real (double), dimension(1:N):: gq4z
real (double), dimension(1:N):: gQ6x
real (double), dimension(1:N):: gQ6y
real (double), dimension(1:N):: gq6z
real (double), dimension(1:N):: gEx
real (double), dimension(1:N):: gEy
real (double), dimension(1:N):: gEz

real (double) :: q4,q6, q4moy, q6moy,w4moy, w6moy! linear displacement

type Trajectoire
    real (double), dimension (0:N-1):: qx  
    real (double), dimension (0:N-1):: qy 
    real (double), dimension (0:N-1):: qz  ! vector of position
    
    real (double), dimension (0:N-1):: px  
    real (double), dimension (0:N-1):: py 
    real (double), dimension (0:N-1):: pz  ! vector of impulsion
   
    real (double), dimension (0:N-1):: uqx 
    real (double), dimension (0:N-1):: uqy
    real (double), dimension (0:N-1):: uqz ! linear displacement
    
    real (double), dimension (0:N-1):: upx 
    real (double), dimension (0:N-1):: upy 
    real (double), dimension (0:N-1):: upz! linear displacement!

    real (double)  Lyap
 endtype

type (Trajectoire), dimension (:,:), allocatable :: P
type (Trajectoire), dimension (:,:), allocatable :: Pfw
type (Trajectoire), dimension (:,:), allocatable :: Pbw

real (double) ::mcconf,ranf, timebw, timefw,dp,l0, alpha,tempo,teq, kapa
real (double) ::seed, utot, tot, x, U,x4,tconf, tconftot,tconfmoy,kinetotdt, upro, q4ref
real(double)::fftot, fftotdt,laputot,laputotdt,ff,lapu, pi2, uproject, dqx, dqy, dqz, dpx, dpy, dpz
real (double) ::inter,AverageTime,NbTimeStep,dt,TotalTime,Time,w,Norm,TimeStore,c1,c2,sqrtgTh,eta1,eta2,gamma,eps, Temperature ! Lyapunov=0,LyapunovAv=0,mu=0,
integer :: nbClones, j,l,k,ran,m, a,i,b, cont,ien,nq4,n0, n4, ltot,mcmoves,totalmcmoves

character(len=128) :: sortie
character(len=128) :: histotot
character(len=128) :: histopart
character(len=128) ::posfinal
character(len=128) ::snap
character(len=128) ::stock
character(len=128) ::stq4t
character(len=128) ::stq6t
character(len=128) ::q4q6

real (double)::energy,rga,beq,rien,q,kine,kinetot,pif, trasm1,trasm2, normgradq4, normgradq6, q6ref,phi, theta, pi,rlimed
real (double) ::xtemp0,xtemp1,xtemp2,xtemp3,xtemp4,xtemp5,ss,delta_x, a_sto
integer, parameter::nfenetre=100
integer :: ienergy, ifenetre, ifen_q6, iter,kaq4,kaq6,niter,kappa,ix, scrivi

real (double):: stat(0:nfenetre, 0:nfenetre)
real (double):: statneg(0:nfenetre, 0:nfenetre)
real (double):: statpos(0:nfenetre, 0:nfenetre)
real (double):: stat6(0:nfenetre, 0:nfenetre)
real(double)::tvmoytot,xn,enprmoy,jtrans, jtot,trasm,trastot,tau

real(double) :: tab_contour(-10:nfenetre+10,0:nfenetre)
real(double) :: tab_contour_q6(-10:nfenetre+10,0:nfenetre)
real(double) :: tab_cont_q4q6(-10:nfenetre+10,-10:nfenetre+10)
real(double) :: cumul_contour(-10:nfenetre+10,0:nfenetre)
real(double) :: cumul_contour_q6(-10:nfenetre+10,0:nfenetre)
real(double) :: cumul_contour_q4q6(-10:nfenetre+10,-10:nfenetre+10)
real(double) :: tab_energy(-10:nfenetre+10,2),toto
real (double), dimension (:), allocatable::norm0
real (double), dimension (:,:), allocatable::dex
real (double)::x111, xmaxq6,energie_min,energie_max,ecinetique,tequilib,delta,beqtot,dhtot,beqtotdt,dhtotdt,intbeq, intdh,e
integer::icfen,icener,icfen_q6,a1, a2, a3, a4, a5,iterfw, iterbw,nmax,totiter, acc


real(double), dimension(:), allocatable::betaq,enerpro
real (double) ,dimension(:), allocatable:: ener0 
real (double) ,dimension(:), allocatable:: ener0bw
real (double) ,dimension(:), allocatable::RLInew
real (double) ,dimension(:), allocatable::rliold 
real (double) ,dimension(:), allocatable::rapporto
real (double) ,dimension(:,:), allocatable:: dh
real (double) ,dimension(:,:), allocatable:: dhx
write(6,*)'usage %s:\n'

!    Initialisation of the variables
write(6,*)'TotalMcmoves '
read (*,*)
read (*,*) Totalmcmoves      ! Total duration of the simulation 
write(*,*) Totalmcmoves 
read (*,*)
write(6,*)'TotalTime '
read (*,*) TotalTime      ! Total duration of the simulation 
write(*,*) TotalTime 
read (*,*)
write(6,*)'dt' 
read (*,*) dt          ! Time step
write(*,*) dt
read (*,*)
write(6,*)'NbClones'  
read (*,*) NbClones  
write(*,*) NbClones    ! Number of clones
read (*,*)
write(6,*)'temperature'  
read (*,*) temperature  ! temperature (kT)
write (*,*) temperature
read (*,*)
write(6,*)'gamma*dt'  
read (*,*) gamma  ! friction*dt
write (*,*) gamma       ! interval between 2 data records
read (*,*)
write(6,*)'sortie'  
read (*,*) sortie   
write (*,*) sortie! Name of the file where data are stored
read (*,*)
write(6,*)'alpha'
read (*,*) alpha 
write (*,*) alpha
read (*,*)
write(6,*)'sigmap'
read (*,*) ss
write (*,*) ss
read (*,*)
write(6,*)'niter'
read (*,*) niter
write (*,*) niter
read (*,*)
write(6,*)'t_equilib'
read (*,*) teq
write (*,*) teq
read (*,*)
write(6,*)'delta_X'
read (*,*) delta_x
write (*,*) delta_x
read (*,*)
write(6,*)'a_sto'
read (*,*) a_sto
write (*,*) a_sto
read (*,*)
write(6,*)'k ressort'
read (*,*) kapa
write (*,*) kapa


open (unit=11, FILE=sortie, action='write', status='replace')
open (unit=31, action='write', status='replace')

histotot  =sortie(1:lenfnam)//'.thout'
histopart =sortie(1:lenfnam)// '.phout'
posfinal = sortie(1:lenfnam)//'.cin'
snap = sortie(1:lenfnam)//'.film'
stock = sortie(1:lenfnam)//'.stock'
stq4t = sortie(1:lenfnam)//'.stq4t'
stq6t= sortie(1:lenfnam)//'.stq6t'
q4q6= sortie(1:lenfnam)//'.qhout'

open(unit=14, file=histopart, action='write', status='replace')
open(unit=26, file=histotot,  action='write', status='replace')
open(unit=28, file=snap,  action='write', status='replace')
open(unit=29, file=stock,  action='write', status='replace')
open(unit=33, file=stq4t,  action='write', status='replace')
open(unit=34, file=stq6t,  action='write', status='replace')
open(unit=31, file=q4q6,  action='write', status='replace')

NMax=NbClones*10      ! Max number of clones
!TimeStore = TotalTime-AverageTime ! Next time at which data will be stored 
                                           
! Array where the weight of every clones is stored
!allocate (w4(0:nmax))
!allocate (w6(0:nmax))
!allocate (wp4(0:nmax))
!allocate (wp6(0:nmax))
allocate (betaq(0:nmax-1))
allocate (ener0(0:Nmax))
allocate (ener0bw(0:Nmax))
allocate(enerpro(0:nmax-1))
! Array which contains the numbers from 1 to Nmax
allocate(Nb(0:NMax-1)) 
allocate(norm0(0:NMax-1)) 
allocate(dex(3,0:N-1))
kappa=1
scrivi=0
ltot=int(totaltime/real(kappa*dt))
totiter=nint(totaltime/dt)
write(*,*)'totiter', totiter
gamma=gamma/(dt)
write(*,*)'gamma' , gamma, gamma*dt,ltot
allocate (xpq4(0:Nmax,0:totiter+10))
allocate (xpq6(0:Nmax,0:totiter+10))
allocate(ener(0:nmax-1,0:totiter+10))

ener(:,:)=0
stat(:,:)=0
stat6(:,:)=0
sigma=1
eta=1
enerpro(:)=0
xpq4(:,:)=0
xpq6(:,:)=0
!statneg(:,:)=0
!statpos(:,:)=0
kine=0
laputot=0
fftot=0
laputotdt=0
fftotdt=0
kinetotdt=0
kinetot=0
mcmoves=0
tau=0
enprmoy=0
nq4=0
e=0.000001
!ss=0.05
l0=0.0!0.12
acc=0

x111   = 0.2
xmaxq6 = 0.6
energie_min = -174.0
energie_max = -144.0

     tab_contour(-10:nfenetre+10,0:nfenetre)=zero
     tab_contour_q6(-10:nfenetre+10,0:nfenetre)=zero
     tab_cont_q4q6(-10:nfenetre+10,-10:nfenetre+10)=zero
     cumul_contour_q4q6(-10:nfenetre+10,-10:nfenetre+10)=0
     cumul_contour(-10:nfenetre+10,0:nfenetre)=0
     cumul_contour_q6(-10:nfenetre+10,0:nfenetre)=0

allocate(Number(0:NbClones-1)) 

!!!!!!!!!!!!T EQUILIBRAGE
tequilib=(dt)*teq
write(*,*) 'tquilib',tequilib
!!!!!!!!!!!! T EQUILIBRAGE
  
write(*,*) 'sortie affichée'

allocate(P(0:NMax-1,0:totiter+10)) 
allocate(Pfw(0:NMax-1,0:totiter+10)) 
allocate(Pbw(0:NMax-1,0:totiter+10)) 
allocate (dh(0:nmax-1,0:totiter+10))
allocate(rliold(0:nbclones-1))
allocate(rlinew(0:nbclones-1))
allocate(rapporto(0:nbclones-1))
allocate (dhx(0:nmax-1,0:totiter+10))

!!EQUiLIBRAGE
xppro(:,:)=0
rapporto(:)=0
iter=0
rliold(:)=0
rlinew(:)=0
rlimed=0
rien=0
Q=0
uproject=0
kappa=1
pi=3.14159
if ((gamma*dt).le.100000) then
 
   rga = exp(-gamma*dt/two)

else 

  rga=0

endif


write(*,*) 'rga',gamma, rga,two,pi

tempo=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! INITIALISATION !!!!!!!!!!!!!!!!!!!!!!
write(*,*) 'conditions initiales pour traj de reference avc distrib stoch à temperature T=',temperature
write(*,*) 'temps equilibrage', tequilib
!!!!!!!!!!!!!!!!!!!!!!!! EQUILIBRAGE initial (STOCH DYN)!!!!!!!!!!!!!!!!!!
     do while ((tempo.lt.Tequilib))!.and.(time.lt.totaltime))
            call calfoljc2(xp,fp,utot)
            call langevin2(xp,vp,fp,dt,rga, temperature,Q,rien,ecinetique)
           ! ff=sum(fp(1,:)**2+fp(2,:)**2+fp(3,:)**2)
            !do i=1,N-1
            !     do j=i+1,N
            !         xij(i,j)=xp(1,i)-xp(1,j)
            !         yij(i,j)=xp(2,i)-xp(2,j)
            !         zij(i,j)=xp(3,i)-xp(3,j)
            !         dij(i,j)=sqrt(xij(i,j)**2+yij(i,j)**2+zij(i,j)**2)
            !         xij(j,i)=-xij(i,j)
            !        yij(j,i)=-yij(i,j)
            !        zij(j,i)=-zij(i,j)
            !        dij(j,i)= dij(i,j)
            !      enddo
            !  enddo
              !x4=q4(N,xij,yij,zij,dij)
           tempo=tempo+dt
       enddo
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

write(*,*) 'point de départ traj de reference deterministe' 
     do  m = 0,NbClones-1  
       P(m,0)%qx=xp(1,:)
       P(m,0)%qy=xp(2,:)
       P(m,0)%qz=xp(3,:)
     enddo

    do m = 0,NbClones-1
       P(m,0)%px=vp(1,1:N)
       P(m,0)%py=vp(2,1:N)
       P(m,0)%pz=vp(3,1:N)
    enddo

 phi=2*pi*genrand()
       theta=pi*genrand()

write(*,*) 'point de départ traj de reference deterministe compagne, pour RLI' 
     do  m = 0,NbClones-1  
       if(mod(m,2).eq.1) then

      
       !write (*,*) 'angoli', theta,phi
       P(m,0)%qx=P(m-1,0)%qx+delta_x*(1.0/sqrt(2*real(N)))*sin(phi)*cos(theta)!genrand()
       P(m,0)%qy=P(m-1,0)%qy+delta_x*(1.0/sqrt(2*real(N)))*sin(phi)*sin(theta)!*genrand()
       P(m,0)%qz=P(m-1,0)%qz+delta_x*(1.0/sqrt(2*real(N)))*cos(phi)!*genrand()
       !phi=2*pi*genrand()
       !theta=2*pi*genrand()
       P(m,0)%px=P(m-1,0)%px+delta_x*(1.0/sqrt(2*real(N)))*sin(phi)*cos(theta)!*genrand()
       P(m,0)%py=P(m-1,0)%py+delta_x*(1.0/sqrt(2*real(N)))*sin(phi)*sin(theta)!*genrand()
       P(m,0)%pz=P(m-1,0)%pz+delta_x*(1.0/sqrt(2*real(N)))*cos(phi)!*genrand()
      
      !call CalculNorme(P(m,0)%qx-P(m-1,0)%qx,P(m,0)%qy-P(m-1,0)%qy,P(m,0)%qz-P(m-1,0)%qz,&
      !                P(m,0)%px-P(m-1,0)%px,P(m,0)%py-P(m-1,0)%py,P(m,0)%pz-P(m-1,0)%pz,N,norm)
      ! write (*,*) 'normdelta_x', m, norm
       endif
    enddo

!write(*,*) 'ctrel' , delta_x!P(0,0)%qx, P(1,0)%QX
      !  do j=0, Nbclones-1
      !         k=0
      !         do i=1,N-1
      !            do a=i+1,N
      !              xij(i,a)=(P(j,k)%qx(i-1)-P(j,k)%qx(a-1))
      !              yij(i,a)=(P(j,k)%qy(i-1)-P(j,k)%qy(a-1))
      !              zij(i,a)=(P(j,k)%qz(i-1)-P(j,k)%qz(a-1))
      !              dij(i,a)=sqrt(xij(i,a)**2+yij(i,a)**2+zij(i,a)**2)
      !              xij(a,i)=-xij(i,a)
      !              yij(a,i)=-yij(i,a)
      !              zij(a,i)=-zij(i,a)
      !             dij(a,i)= dij(i,a)
      !            enddo
      !          enddo

       !     xpq4(j,k)=q4(N,xij,yij,zij,dij)
       !     xpq6(j,k)=q6(N,xij,yij,zij,dij)
       !     call Force(P(j,k)%qx,P(j,k)%qy,P(j,k)%qz,fp,Utot)
       !     ener(j,k)=utot
           ! write(*,'(a,i, i, e11.4, e12.4, e12.4)') 'trajinit0',j, k,xpq4(j,k),xpq6(j,k),ener(j,k)
      ! enddo

do j=0,N-1
       phi=2*pi*genrand()
       theta=pi*(genrand())
        
       dex(1,j)=0.0001*(1.0/sqrt(2*real(N)))*sin(phi)*cos(theta)
       dex(2,j)=0.0001*(1.0/sqrt(2*real(N)))*sin(phi)*sin(theta)    
       dex(3,j)=0.0001*(1.0/sqrt(2*real(N)))*cos(phi)

enddo
!dqx=0.5*genrand()
!dqy=0.5*genrand()
!dqz=0.5*genrand()
!dpx=0.5*genrand()
!dpy=0.5*genrand()
!dpz=0.5*genrand()
!dqx=genrand()
!dqy=genrand()
!dqz=genrand()
!dpx=genrand()
!dpy=genrand()
!dpz=genrand()

 !Norm=dqx*dqx+dqy*dqy+dqz*dqz+dpz*dpz+dpx*dpx+dpy*dpy
 !Norm=sqrt(Norm)
! write(*,*) 'dqx', dqx, dpx  

       phi=2*pi*genrand()
       theta=pi*(genrand())

write (*,*) 'initialisation DES UQ, UP pour calcul des lambda' 
   do  m = 0,NbClones-1 
        
     do j=0,N-1

       
       P(m,0)%uqx(j)=P(m,0)%qx(j)+dex(1,j)
       P(m,0)%uqy(j)=P(m,0)%qy(j)+dex(2,j)
       P(m,0)%uqz(j)=dex(3,j)+P(m,0)%qz(j)
       
       P(m,0)%upx(j)=dex(1,j)+P(m,0)%px(j)
       P(m,0)%upy(j)=dex(2,j)+P(m,0)%py(j)
       P(m,0)%upz(j)=dex(3,j)+P(m,0)%pz(j)
      enddo
 
       !call Normalizza(P(m,0)%uqx,P(m,0)%uqy,P(m,0)%uqz,P(m,0)%upx,P(m,0)%upy,P(m,0)%upz, N)
!call CalculNorme(0.05*genrand(),0.05*genrand(),0.05*genrand(),0.05*genrand(),&
 !                       0.05*genrand(),0.05*genrand(),N,Norm)
       call CalculNorme(P(m,0)%uqx-P(m,0)%qx,P(m,0)%uqy-P(m,0)%qy,P(m,0)%uqz-P(m,0)%qz,P(m,0)%upx-P(m,0)%px,&
                        P(m,0)%upy-P(m,0)%py,P(m,0)%upz-P(m,0)%pz,N,Norm)
       norm0(m)=norm
      ! write (*,*) 'dx2', norm0(m)
!write(*,*) 'in', P(m,0)%qx(1), P(m,0)%uqx(1)
    enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! metto a zero il resto avanzato
 P(:,:)%Lyap = 0

 do l=1,totiter 
  do m = NbClones,NMax-1
    do j=0,N-1
    P(m,l)%qx(j) =0
    P(m,l)%px(j) =0
    P(m,l)%uqx(j)=0
    P(m,l)%upx(j)=0
 
    P(m,l)%qy(j) =0
    P(m,l)%py(j) =0
    P(m,l)%uqy(j)=0
    P(m,l)%upy(j)=0

    P(m,l)%qz(j) =0
    P(m,l)%pz(j) =0
    P(m,l)%uqz(j)=0
    P(m,l)%upz(j)=0
   enddo
   ! P(m).Lyap = 0
  enddo
enddo


ecinetique=0
dh(:,:)=0
write(*,*) 'initialisation faite: go with RLI' 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  c'est parti !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do j=0,NbClones-1 
      time=0
 !write(*,*) 'in1', P(j,0)%qx(1), P(j,0)%uqx(1) 

     do while (Time.lt.Totaltime) !!!!!!!!!!!!!!!!!!!!!!  propagazione prima traiettoria   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           iter=nint(real(time)/real(dt))
           !write(*,*)'!!', iter, P(j,iter)%qx(1), P(j,iter)%uqx(1)
           ! write(*,*) 'in1', P(j,0)%qx(1), P(j,0)%uqx(1) 
           call propagatfwdet(P(j,iter)%qx,P(j,iter)%qy,P(j,iter)%qz,P(j,iter)%px,P(j,iter)%py,P(j,iter)%pz,dt,N,dh(j,iter))
          !write(*,*) '1', P(j,iter)%qx(1)
           call propagatfwdet(P(j,iter)%Uqx,P(j,iter)%uqy,P(j,iter)%uqz,P(j,iter)%upx,P(j,iter)%upy,P(j,iter)%upz,dt,N,dhx(j,iter))
           !write(*,*) '1', P(j,iter)%uqx(1)
           call CalculNorme(P(j,iter)%uqx-P(j,iter)%qx,P(j,iter)%uqy-P(j,iter)%qy,P(j,iter)%uqz-P(j,iter)%qz,&
                            P(j,iter)%upx-P(j,iter)%px,P(j,iter)%upy-P(j,iter)%py,P(j,iter)%upz-P(j,iter)%upz,N,norm)
           !!!!!!!!!!!!!!!!!!!!!! STOCH DYN
           ! call propagateA (P(j,iter)%qx,P(j,iter)%qy,P(j,iter)%qz,P(j,iter)%px,P(j,iter)%py,P(j,iter)%pz,&
           !                  P(j,iter)%uqx,P(j,iter)%uqy,P(j,iter)%uqz,P(j,iter)%upx,P(j,iter)%upy,&
           !                  P(j,iter)%upz,rga,temperature,dt,N,beq,dh(j,iter),ecinetique,ff,lapU)
           ! call CalculNorme(P(j,iter)%uqx,P(j,iter)%uqy,P(j,iter)%uqz,P(j,iter)%upx,P(j,iter)%upy,P(j,iter)%upz,N,norm)
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !write (*,*) 'norm',j, norm!if (mod(iter,kappa).eq.0) then
            !if (iter.eq.0) then 
            P(j,iter)%Lyap=log(Norm/norm0(j))
            !write (*,*) 'lyap',j, iter, P(j,iter)%Lyap!/real(iter)
            !else
            !P(j,iter)%Lyap =log(Norm/norm0(j))+ P(j,iter-1)%Lyap
            !write (*,*) 'lyap', j, iter, P(j,iter)%Lyap!/real(iter)
            !endif
            !write (*,*) 'dh', mcmoves, iter, j, dh(j,iter)
            !endif
            P(j,iter+1)%qx=P(j,iter)%qx
            P(j,iter+1)%qy=P(j,iter)%qy
            P(j,iter+1)%qz=P(j,iter)%qz

            P(j,iter+1)%px=P(j,iter)%px
            P(j,iter+1)%py=P(j,iter)%py
            P(j,iter+1)%pz=P(j,iter)%pz
   
            P(j,iter+1)%uqx=P(j,iter)%uqx
            P(j,iter+1)%uqy=P(j,iter)%uqy
            P(j,iter+1)%uqz=P(j,iter)%uqz

            P(j,iter+1)%upx=P(j,iter)%upx
            P(j,iter+1)%upy=P(j,iter)%upy
            P(j,iter+1)%upz=P(j,iter)%upz
         
       time=time+dt
   enddo  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !write(*,*) P(1,totiter-1)%qx, P(1,totiter-1)%uqx
!write(*,*) 'verolyap', P(j,totiter-1)%Lyap/real(totiter-1)
!stop
     call Force(P(j,0)%qx,P(j,0)%qy,P(j,0)%qz,fp,Utot)
       !write (*,*) 'yy', utot
        ener0(j)=dh(j,totiter-1)
        !write (*,*) 'yy', ener0(j)
           do k=0,totiter-1
               do i=1,N-1
                  do a=i+1,N
                    xij(i,a)=(P(j,k)%qx(i-1)-P(j,k)%qx(a-1))
                    yij(i,a)=(P(j,k)%qy(i-1)-P(j,k)%qy(a-1))
                    zij(i,a)=(P(j,k)%qz(i-1)-P(j,k)%qz(a-1))
                    dij(i,a)=sqrt(xij(i,a)**2+yij(i,a)**2+zij(i,a)**2)
                    xij(a,i)=-xij(i,a)
                    yij(a,i)=-yij(i,a)
                    zij(a,i)=-zij(i,a)
                    dij(a,i)= dij(i,a)
                  enddo
                enddo

            xpq4(j,k)=q4(N,xij,yij,zij,dij)
            xpq6(j,k)=q6(N,xij,yij,zij,dij)
               if ((xpq4(j,k).ge.0.1).and.(xpq4(j,k).ge.0.3)) then
               tau=tau+1
               end if
            call Force(P(j,k)%qx,P(j,k)%qy,P(j,k)%qz,fp,Utot)
            ener(j,k)=utot
            if(mod(k,100).eq.0) then
            write(*,'(a, i, i,e11.4, e12.4, e12.4, e12.4)') 'traj0', j, k,xpq4(j,k),xpq6(j,k),ener(j,k),dh(j,k)
            endif
 enddo
          
xpq4(:,:)=0
xpq6(:,:)=0
xij(:,:)=0
yij(:,:)=0
zij(:,:)=0
dij(:,:)=0

          !do k=0,totiter-1
          !     do i=1,N-1
          !        do a=i+1,N
          !          xij(i,a)=(P(j,k)%uqx(i-1)-P(j,k)%uqx(a-1))
          !          yij(i,a)=(P(j,k)%uqy(i-1)-P(j,k)%uqy(a-1))
          !          zij(i,a)=(P(j,k)%uqz(i-1)-P(j,k)%uqz(a-1))
          !          dij(i,a)=sqrt(xij(i,a)**2+yij(i,a)**2+zij(i,a)**2)
          !          xij(a,i)=-xij(i,a)
          !          yij(a,i)=-yij(i,a)
           !         zij(a,i)=-zij(i,a)
          !          dij(a,i)= dij(i,a)
          !        enddo
          !      enddo

           ! xpq4(j,k)=q4(N,xij,yij,zij,dij)
           ! xpq6(j,k)=q6(N,xij,yij,zij,dij)
           ! call Force(P(j,k)%uqx,P(j,k)%uqy,P(j,k)%uqz,fp,Utot)
           ! ener(j,k)=utot
           ! write(*,'(a, i, e11.4, e12.4, e12.4, e12.4)') 'trajdelta', k,xpq4(j,k),xpq6(j,k),ener(j,k),dh(j,k)


           !!!!!!!!!! traiettorie reattive !!!!!!!!!!!!!!!!
           !if ((xpq4(j,k).ge.0.1).and.(xpq4(j,k).ge.0.3)) then
           !tau=tau+1
           !endif
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !enddo
 
enddo
              !  do k=0,totiter-1
                  !    if (mod(k,100).eq.0) then
                  !    write(*,'(a, i, e11.4, e12.4, e12.4, e12.4)') 'trajinit', k,xpq4(0,k),xpq6(0,k),ener(0,k),dh(0,k)
                  !    endif
                  !  enddo

do j=0,nbclones-1
  do k=1,totiter-1
      if(mod(j,2).eq.0)then
      RLIold(j)=Rliold(j)+abs(P(j+1,k)%Lyap-P(j,k)%Lyap)/real(k)
      else 
      RLIold(j)=Rliold(j)+abs(P(j-1,k)%Lyap-P(j,k)%Lyap)/real(k)
      endif
      !rlimed=rlimed+RLIold(j)
  enddo
      if(mod(j,2).eq.0)then
      RLIold(j)=Rliold(j)+abs(P(j+1,0)%Lyap-P(j,0)%Lyap)
      else 
      RLIold(j)=Rliold(j)+abs(P(j-1,0)%Lyap-P(j,0)%Lyap)
      endif

!write(*,*) 'rliold' , j, rliold(j)/real(totiter-1)
enddo

!write(*,*)'rlimedio' , rlimed/real(nbclones)/real(totiter-1)

RLIold(:)=Rliold(:)/real(totiter-1)

! rliold(:)=0

norm0(:)=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do while (mcmoves.lt.Totalmcmoves) !!!!!!!!!!!!! MO' inizio montecarlo su traiettorie

!!!!!!!!!!!!  fACCIAMO SHOOTING:
  ix=nint(Ltot*kappa*genrand()) ! PUNTO DI SHOOTING ltot numero di punti stockati (segmenti di traiettoria)
  timefw=ix*dt!*kappa
  timebw=ix*dt!*kappa
  !write(*,*) 'shooting time', ix, timefw, timebw

!!!! shooting deterministico 
  !call gauss(ss,l0,dp)
  
  !write (*,*) 'dp', dp
!!!per lo shooting stocastico perturbo anche le posizioni ....

do j=0,nbclones-1
 !!!!! inizializzo backward 
 if(mod(j,2).eq.0) then
     Pbw(j,ix)%qx=P(j,ix)%qx
     Pbw(j,ix)%qy=P(j,ix)%qy
     Pbw(j,ix)%qz=P(j,ix)%qz

!!!! shooting classique
    
do m=0, N-1

  call gauss(ss,l0,dpx)
  call gauss(ss,l0,dpy)
  call gauss(ss,l0,dpz)

    Pbw(j,ix)%px(m)=P(j,ix)%px(m)+dpx
    Pbw(j,ix)%py(m)=P(j,ix)%py(m)+dpy
    Pbw(j,ix)%pz(m)=P(j,ix)%pz(m)+dpz
    enddo

!!!! shooting 'Stoltz'
    !Pbw(j,ix)%px(m)=a_sto*P(j,ix)%px(m)+sqrt(1.0-a_sto*a_sto)*dpx 
    !Pbw(j,ix)%py(m)=a_sto*P(j,ix)%py(m)+sqrt(1.0-a_sto*a_sto)*dpy
    !Pbw(j,ix)%pz(m)=a_sto*P(j,ix)%pz(m)+sqrt(1.0-a_sto*a_sto)*dpz
    !enddo

 else 
     Pbw(j,ix)%qx=P(j-1,ix)%qx
     Pbw(j,ix)%qy=P(j-1,ix)%qy
     Pbw(j,ix)%qz=P(j-1,ix)%qz

     do m=0, N-1

     call gauss(ss,l0,dpx)
     call gauss(ss,l0,dpy)
     call gauss(ss,l0,dpz)

!!!!! shooting classique
     Pbw(j,ix)%px(m)=P(j-1,ix)%px(m)+dpx
     Pbw(j,ix)%py(m)=P(j-1,ix)%py(m)+dpy
     Pbw(j,ix)%pz(m)=P(j-1,ix)%pz(m)+dpz
     enddo

!!!! shooting 'Stoltz'
     !Pbw(j,ix)%px(m)=a_sto*P(j-1,ix)%px(m)+sqrt(1.0-a_sto*a_sto)*dpx 
     !Pbw(j,ix)%py(m)=a_sto*P(j-1,ix)%py(m)+sqrt(1.0-a_sto*a_sto)*dpy
     !Pbw(j,ix)%pz(m)=a_sto*P(j-1,ix)%pz(m)+sqrt(1.0-a_sto*a_sto)*dpz
     !enddo

 endif

!write(*,*) 'emmo', j, Pbw(j,ix)%qx

timebw=ix*dt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! propagazione BACKWARD a partire da punto di shoting (senza lyapunov)
  do while (timebw.ge.zero)
        iterbw=nint(timebw/dt)
        call propagatbckw(Pbw(j,iterbw)%qx,Pbw(j,iterbw)%qy,Pbw(j,iterbw)%qz,Pbw(j,iterbw)%px,Pbw(j,iterbw)%py,&
                          Pbw(j,iterbw)%pz,rga,temperature,dt,N,dh(j,iterbw))
    
      Pbw(j,iterbw-1)%qx=Pbw(j,iterbw)%qx
      Pbw(j,iterbw-1)%qy=Pbw(j,iterbw)%qy
      Pbw(j,iterbw-1)%qz=Pbw(j,iterbw)%qz

      Pbw(j,iterbw-1)%Px=Pbw(j,iterbw)%Px
      Pbw(j,iterbw-1)%Py=Pbw(j,iterbw)%Py
      Pbw(j,iterbw-1)%Pz=Pbw(j,iterbw)%Pz

    timebw=timebw-dt
 enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!! calcolo energia di config iniziale traiettoria BW
 ! call Force(Pbw(j,0)%qx,Pbw(j,0)%qy,Pbw(j,0)%qz,fp,Utot)
  ener0bw(j)=dh(j,ix)!utot
enddo

phi=2*pi*genrand()
theta=pi*(genrand())


do j=0,nbclones-1
 
!!!!!!!!!!!!!!!!!!!!! la traiettoria bw é arrivata a tempo t=0
!!!!!!!!!metto delta_x per rli su clone 1
if (mod(j,2).eq.1) then

       !phi=2*pi*genrand()
       !theta=2*pi*genrand()
       !write (*,*) 'angoli', theta,phi
       
     Pbw(j,0)%qx=Pbw(j,0)%qx+delta_x*(1.0/sqrt(2*real(N)))*sin(phi)*cos(theta)
     Pbw(j,0)%qy=Pbw(j,0)%qy+delta_x*(1.0/sqrt(2*real(N)))*sin(phi)*sin(theta)
     Pbw(j,0)%qz=Pbw(j,0)%qz+delta_x*(1.0/sqrt(2*real(N)))*cos(phi)

     Pbw(j,0)%px=Pbw(j,0)%px+delta_x*(1.0/sqrt(2*real(N)))*sin(phi)*cos(theta)
     Pbw(j,0)%py=Pbw(j,0)%py+delta_x*(1.0/sqrt(2*real(N)))*sin(phi)*sin(theta)
     Pbw(j,0)%pz=Pbw(j,0)%pz+delta_x*(1.0/sqrt(2*real(N)))*cos(phi)

    ! call CalculNorme(Pbw(j,0)%qx-Pbw(0,0)%qx,Pbw(j,0)%qy-Pbw(0,0)%qy,Pbw(j,0)%qz-Pbw(0,0)%qz,&
    !                  Pbw(j,0)%px-Pbw(0,0)%px,Pbw(j,0)%py-Pbw(0,0)%py,Pbw(j,0)%pz-Pbw(0,0)%pz,N,norm)
     !write (*,*) 'normdelta_xbw', j,  norm!*sqrt(2.0)
endif

enddo
 
   !phi=2.0*pi*genrand()
   !theta=pi*(genrand())

do j=0,nbclones-1
!!!!!!!!!!!!!!!!!!!!!!!!!!! Lyapunov DETERMINISTICO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!! a t=0 re-inizializzo i vettorini per ricalcolare il lyapunov
       
      do m=0,N-1
       
        Pbw(j,0)%uqx(m)=dex(1,m)+Pbw(j,0)%qx(m)
        Pbw(j,0)%uqy(m)=dex(2,m)+Pbw(j,0)%qy(m)
        Pbw(j,0)%uqz(m)=dex(3,m)+Pbw(j,0)%qz(m)

        Pbw(j,0)%upx(m)=dex(1,m)+Pbw(j,0)%px(m)
        Pbw(j,0)%upy(m)=dex(2,m)+Pbw(j,0)%py(m)
        Pbw(j,0)%upz(m)=dex(3,m)+Pbw(j,0)%pz(m)
       enddo

       call CalculNorme(Pbw(j,0)%uqx-Pbw(j,0)%qx,Pbw(j,0)%uqy-Pbw(j,0)%qy,Pbw(j,0)%uqz-Pbw(j,0)%qz,&
                        Pbw(j,0)%upx-Pbw(j,0)%px,Pbw(j,0)%upy-Pbw(j,0)%py,Pbw(j,0)%upz-Pbw(j,0)%pz,N,norm)
       norm0(j)=norm
!call Force(Pbw(j,0)%uqx,Pbw(j,0)%uqy,Pbw(j,0)%uqz,fp,Utot)
 !write (*,*) 'norbw0',norm0(j)
!enddo
   timebw=0
!do j=0,nbclones-1
  do while (timebw.le.ix*dt) !!!!!!!!!!!!!!!!!!!!!!! adesso ripercorro traiettoria BW calcolando pero` i lyapunov
        iterbw=nint(timebw/dt)
           ! if (iterbw.eq.ix) then
           !         exit
           ! else 
                call propagatfwdet(Pbw(j,iterbw)%uqx,Pbw(j,iterbw)%uqy,Pbw(j,iterbw)%uqz,Pbw(j,iterbw)%upx,Pbw(j,iterbw)%upy,Pbw(j,iterbw)%upz,dt,N,dhx(j,iterBw))
                call CalculNorme(Pbw(j,iterbw)%uqx-Pbw(j,iterbw)%qx,Pbw(j,iterbw)%uqy-Pbw(j,iterbw)%qy,Pbw(j,iterbw)%uqz-Pbw(j,iterbw)%qz,&
                                 Pbw(j,iterbw)%upx-Pbw(j,iterbw)%px,Pbw(j,iterbw)%upy-Pbw(j,iterbw)%py,Pbw(j,iterbw)%upz-Pbw(j,iterbw)%upz,N,norm)
                !if (iterbw.eq.0) then 
                Pbw(j,iterbw)%Lyap=log(Norm/norm0(j))
                !write (*,*) 'lyapbw',norm
                !else
                !Pbw(j,iterbw)%Lyap=Pbw(j,iterbw-1)%Lyap+log(Norm/norm0(j))
                !write (*,*) 'lyapbw',j, iterbw, Pbw(j,iterbw)%Lyap!/real(iterbw),norm
                !endif
                !write (*,*) 'dhbw', mcmoves, iterbw, j, dh(j,iterbw)
                !endif
                Pbw(j,iterbw+1)%uqx=Pbw(j,iterbw)%uqx
                Pbw(j,iterbw+1)%uqy=Pbw(j,iterbw)%uqy
                Pbw(j,iterbw+1)%uqz=Pbw(j,iterbw)%uqz

                Pbw(j,iterbw+1)%upx=Pbw(j,iterbw)%upx
                Pbw(j,iterbw+1)%upy=Pbw(j,iterbw)%upy
                Pbw(j,iterbw+1)%upz=Pbw(j,iterbw)%upz
         !endif
         timebw=timebw+dt
   enddo
enddo



do j=0,nbclones-1
  do k=1,ix
   !write (*,*) 'deliap', k, abs(Pbw(1,k)%Lyap-Pbw(0,k)%Lyap)/real(k)
      if(mod(j,2).eq.0)then
         RLInew(j)=Rlinew(j)+(abs(Pbw(j+1,k)%Lyap-Pbw(j,k)%Lyap)/real(k))
      else 
         RLInew(j)=Rlinew(j)+(abs(Pbw(j-1,k)%Lyap-Pbw(j,k)%Lyap)/real(k))
      endif
     !write(*,*) 'rli2', k, RLinew(j)
  enddo
   if(mod(j,2).eq.0)then
         RLInew(j)=Rlinew(j)+abs(Pbw(j+1,0)%Lyap-Pbw(j,0)%Lyap)
      else 
         RLInew(j)=Rlinew(j)+abs(Pbw(j-1,0)%Lyap-Pbw(j,0)%Lyap)
      endif
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! STOCASTICO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Do m=0,N-1
!        Pbw(j,0)%uqx(m)=genrand()
!        Pbw(j,0)%uqy(m)=genrand()
!        Pbw(j,0)%uqz(m)=genrand()
!
!        Pbw(j,0)%upx(m)=genrand()
!        Pbw(j,0)%upy(m)=genrand()
!        Pbw(j,0)%upz(m)=genrand() 
 !      enddo

!       call Normalizza(Pbw(j,0)%uqx,Pbw(j,0)%uqy,Pbw(j,0)%uqz,Pbw(j,0)%upx,Pbw(j,0)%upy,Pbw(j,0)%upz,N)
!       call CalculNorme(Pbw(j,0)%uqx,Pbw(j,0)%uqy,Pbw(j,0)%uqz,Pbw(j,0)%upx, Pbw(j,0)%upy,Pbw(j,0)%upz,N,Norm)

!   do while (timebw.lt.ix*dt) !!!!!!!!!!!!!!!!!!!!!!! adesso ripercorro traiettoria BW calcolando pero` i lyapunov
!        iterbw=nint(timebw/dt)
!write (*,*) 'lyapbw', iterbw
!            if (iterbw.eq.ix) then
!                    exit
!            else 
!                call reconstructLy(Pbw(j,iterbw)%qx,Pbw(j,iterbw)%qy,Pbw(j,iterbw)%qz,Pbw(j,iterbw+1)%qx,Pbw(j,iterbw+1)%qy,&
!                             Pbw(j,iterbw+1)%qz,Pbw(j,iterbw)%uqx,Pbw(j,iterbw)%uqy,Pbw(j,iterbw)%uqz,Pbw(j,iterbw)%upx,&
!                             Pbw(j,iterbw)%upy,Pbw(j,iterbw)%upz,rga,temperature,dt,N)
                   !if (mod(iterbw,kappa)==0) then

!                call CalculNorme(Pbw(j,iterbw)%uqx,Pbw(j,iterbw)%uqy,Pbw(j,iterbw)%uqz,Pbw(j,iterbw)%upx,Pbw(j,iterbw)%upy,Pbw(j,iterbw)%upz,N,norm)
!                if (iterbw.eq.0) then 
!                Pbw(j,iterbw)%Lyap=log(Norm)
!                write (*,*) 'lyapbw',norm
!                else
!                Pbw(j,iterbw)%Lyap=Pbw(j,iterbw-1)%Lyap+log(Norm)
!                write (*,*) 'lyapbw', iterbw, Pbw(j,iterbw)%Lyap/real(iterbw),norm
!                endif
       
          ! write (*,*) 'dhbw', mcmoves, iterbw, j, dh(j,iterbw)
       !endif
!                Pbw(j,iterbw+1)%uqx=Pbw(j,iterbw)%uqx
!                Pbw(j,iterbw+1)%uqy=Pbw(j,iterbw)%uqy
!                Pbw(j,iterbw+1)%uqz=Pbw(j,iterbw)%uqz

!                Pbw(j,iterbw+1)%upx=Pbw(j,iterbw)%upx
!                Pbw(j,iterbw+1)%upy=Pbw(j,iterbw)%upy
!                Pbw(j,iterbw+1)%upz=Pbw(j,iterbw)%upz
!         endif

!         timebw=timebw+dt
!   enddo
!!!!!!!!!!!!!!!!!!!! PUNTO DI CONGIUNZIONE BW-FW
  
do j=0,nbclones-1

     Pfw(j,ix)%qx=Pbw(j,ix)%qx
     Pfw(j,ix)%qy=Pbw(j,ix)%qy
     Pfw(j,ix)%qz=Pbw(j,ix)%qz

     Pfw(j,ix)%Px=Pbw(j,ix)%Px
     Pfw(j,ix)%Py=Pbw(j,ix)%Py
     Pfw(j,ix)%Pz=Pbw(j,ix)%Pz

     Pfw(j,ix)%uqx=Pbw(j,ix)%uqx
     Pfw(j,ix)%uqy=Pbw(j,ix)%uqy
     Pfw(j,ix)%uqz=Pbw(j,ix)%uqz

     Pfw(j,ix)%upx=Pbw(j,ix)%upx
     Pfw(j,ix)%upy=Pbw(j,ix)%upy
     Pfw(j,ix)%upz=Pbw(j,ix)%upz
 
     !Pfw(j,ix-1)%Lyap=Pbw(j,ix-1)%Lyap

 !write (*,*) 'lup', Pbw(j,ix-1)%upx 
! write (*,*) 'lupx',Pbw(j,ix)%upx
!write (*,*) 'lyapfcong', ix, Pfw(j,ix-1)%Lyap,Pbw(j,ix-1)%Lyap 
timefw=ix*dt
!!!!!!!!!!!! trajectoire forward !!!!!!!!!!!!!!!!!!!!
  do while (timefw.lt.totaltime-dt)
          iterfw=nint(timefw/dt)

                call propagatfwDET(Pfw(j,iterfw)%qx,Pfw(j,iterfw)%qy,Pfw(j,iterfw)%qz,Pfw(j,iterfw)%px,Pfw(j,iterfw)%py,&
                                   Pfw(j,iterfw)%pz,dt,N,dh(j,iterfw))
                call propagatfwDET(Pfw(j,iterfw)%Uqx,Pfw(j,iterfw)%uqy,Pfw(j,iterfw)%uqz,Pfw(j,iterfw)%upx,Pfw(j,iterfw)%upy,&
                                   Pfw(j,iterfw)%upz,dt,N,dhx(j,iterfw))
                !if (mod(iter,kappa)==0) then
                call CalculNorme(Pfw(j,iterfw)%uqx-Pfw(j,iterfw)%qx,Pfw(j,iterfw)%uqy-Pfw(j,iterfw)%qy,Pfw(j,iterfw)%uqz-Pfw(j,iterfw)%qz,&
                                 Pfw(j,iterfw)%upx-Pfw(j,iterfw)%px,Pfw(j,iterfw)%upy-Pfw(j,iterfw)%py,Pfw(j,iterfw)%upz-Pfw(j,iterfw)%pz,N,norm)
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  STOCASTICO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !call propagatfw(Pfw(j,iterfw)%qx,Pfw(j,iterfw)%qy,Pfw(j,iterfw)%qz,Pfw(j,iterfw)%px,Pfw(j,iterfw)%py,&
                !                Pfw(j,iterfw)%pz,Pfw(j,iterfw)%uqx,Pfw(j,iterfw)%uqy,Pfw(j,iterfw)%uqz,&
                !                Pfw(j,iterfw)%upx,Pfw(j,iterfw)%upy,Pfw(j,iterfw)%upz,rga,temperature,dt,N,dh(j,iterfw))
                !if (mod(iter,kappa)==0) then
                !call CalculNorme(Pfw(j,iterfw)%uqx,Pfw(j,iterfw)%uqy,Pfw(j,iterfw)%uqz,Pfw(j,iterfw)%upx,Pfw(j,iterfw)%upy,Pfw(j,iterfw)%upz,N,norm)
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !if (iterfw.gt.ix-1) then 
                !Pfw(j,iterfw)%Lyap=Pfw(j,iterfw-1)%Lyap+log(Norm/norm0(j))
                Pfw(j,iterfw+1)%Lyap=log(Norm/norm0(j)) 
               ! write (*,*) 'lyapfw',j, iterfw+1, Pfw(j,iterfw+1)%Lyap!/real(iterfw),norm
                !endif
                !write (*,*) 'dhfw', mcmoves, iterfw, j, dh(j,iterfw)
                !endif

                Pfw(j,iterfw+1)%qx=Pfw(j,iterfw)%qx
                Pfw(j,iterfw+1)%qy=Pfw(j,iterfw)%qy
                Pfw(j,iterfw+1)%qz=Pfw(j,iterfw)%qz

                Pfw(j,iterfw+1)%Px=Pfw(j,iterfw)%Px
                Pfw(j,iterfw+1)%Py=Pfw(j,iterfw)%Py
                Pfw(j,iterfw+1)%Pz=Pfw(j,iterfw)%Pz
   
                Pfw(j,iterfw+1)%uqx=Pfw(j,iterfw)%uqx
                Pfw(j,iterfw+1)%uqy=Pfw(j,iterfw)%uqy
                Pfw(j,iterfw+1)%uqz=Pfw(j,iterfw)%uqz

                Pfw(j,iterfw+1)%uPx=Pfw(j,iterfw)%UPx
                Pfw(j,iterfw+1)%uPy=Pfw(j,iterfw)%UPy
                Pfw(j,iterfw+1)%uPz=Pfw(j,iterfw)%UPz
           
           timefw=timefw+dt
    enddo

             !  do i=1,N-1
             !     do a=i+1,N
             !       xij(i,a)=(Pfw(j,totiter-1)%qx(i-1)-Pfw(j,totiter-1)%qx(a-1))
             !       yij(i,a)=(Pfw(j,totiter-1)%qy(i-1)-Pfw(j,totiter-1)%qy(a-1))
             !       zij(i,a)=(Pfw(j,totiter-1)%qz(i-1)-Pfw(j,totiter-1)%qz(a-1))
             !       dij(i,a)=sqrt(xij(i,a)**2+yij(i,a)**2+zij(i,a)**2)
             !       xij(a,i)=-xij(i,a)
             !       yij(a,i)=-yij(i,a)
             !       zij(a,i)=-zij(i,a)
             !       dij(a,i)= dij(i,a)
              !    enddo
              !  enddo

!write(*,*) 'q4q6', q4(N,xij,yij,zij,dij), q6(N,xij,yij,zij,dij)
enddo

do j=0,nbclones-1
  do k=ix+1,totiter-1
   !    write (*,*) 'deliap', k, abs(Pfw(1,k)%Lyap-Pfw(0,k)%Lyap)/real(k)
      if(mod(j,2).eq.0)then
      RLInew(j)=Rlinew(j)+abs(Pfw(j+1,k)%Lyap-Pfw(j,k)%Lyap)/real(k)
      else 
      RLInew(j)=Rlinew(j)+abs(Pfw(j-1,k)%Lyap-Pfw(j,k)%Lyap)/real(k)
      endif
  enddo
!write(*,*) 'rlinew' , j, rlinew(j)/real(totiter-1)
enddo




RLinew(:)= RLinew(:)/real(totiter-1)

!write(*,*) 'rlinew', mcmoves, RLinew(:)
!stop
!write(*,*) 'deltarli1', RLinew(:)-RLIold(:)

do j=0,nbclones-1
     !!!!!!!!!!!! shooting classique
     ! rapporto(j)=exp(alpha*(RLInew(j)-RLIold(j)))*exp(-(ener0bw(j)-ener0(j))/temperature)
      rapporto(j)=exp(alpha*(RLInew(j)-RLIold(j)))*exp(-(ener0bw(j)-ener0(j))/temperature)
     
   ! write(*,'(a, e12.4, e12.4, e12.4)') 'deltarli', RLinew(j)-RLIold(j), ener0bw(j)-ener0(j),rapporto(j)
!write(*,*) 'uno', j,mcmoves,exp(alpha*(RLInew(j)-RLIold(j)))
!write(*,'(a, e12.4, e12.4, e12.4, e12.4)')'deltarli', RLinew(j)-RLIold(j),exp(alpha*(RLInew(j)-RLIold(j))), exp(-(ener0bw(j)-ener0(j))/temperature), rapporto(j)
 
 ranf=genrand()
  mcconf=min(1.0,rapporto(j))
!write(*,*) 'sel',rapporto(j),mcconf,ranf
     
     if (ranf.lt.mcconf) then
       if (j.eq.0) acc=acc+1
     !write (*,*) 'cambio traiettoria'
      
       do i=0,ix-1
           P(j,i)%qx=Pbw(j,i)%qx
           P(j,i)%qy=Pbw(j,i)%qy
           P(j,i)%qz=Pbw(j,i)%qz

           P(j,i)%px=Pbw(j,i)%px
           P(j,i)%py=Pbw(j,i)%py
           P(j,i)%pz=Pbw(j,i)%pz
        enddo

        do i=ix,totiter
           P(j,i)%qx=Pfw(j,i)%qx
           P(j,i)%qy=Pfw(j,i)%qy
           P(j,i)%qz=Pfw(j,i)%qz

           P(j,i)%px=Pfw(j,i)%px
           P(j,i)%py=Pfw(j,i)%py
           P(j,i)%pz=Pfw(j,i)%pz
         enddo

        Rliold(j)=rlinew(j)
        ener0(j)=ener0bw(j)
      
     ! write (*,*) 'tiro altro num aleatorio'
     endif
       
           do k=0,totiter-1
               do i=1,N-1
                  do a=i+1,N
                    xij(i,a)=(P(j,k)%qx(i-1)-P(j,k)%qx(a-1))
                    yij(i,a)=(P(j,k)%qy(i-1)-P(j,k)%qy(a-1))
                    zij(i,a)=(P(j,k)%qz(i-1)-P(j,k)%qz(a-1))
                    dij(i,a)=sqrt(xij(i,a)**2+yij(i,a)**2+zij(i,a)**2)
                    xij(a,i)=-xij(i,a)
                    yij(a,i)=-yij(i,a)
                    zij(a,i)=-zij(i,a)
                    dij(a,i)= dij(i,a)
                  enddo
                enddo

            xpq4(j,k)=q4(N,xij,yij,zij,dij)
            xpq6(j,k)=q6(N,xij,yij,zij,dij)
            call Force(P(j,k)%qx,P(j,k)%qy,P(j,k)%qz,fp,Utot)
            ener(j,k)=utot
            
            if((xpq4(j,k).le.0.15).or.(xpq6(j,k).le.0.4)) then
            tau=tau+1
            scrivi=1
            endif
        enddo

 enddo !! clones

           if ((mod(mcmoves,50).eq.0).or.(scrivi.eq.1)) then 
              do k=0,totiter-1   
              if (mod(k,100).eq.0) then
              write(*,'(a,i, i, e11.4, e12.4, e12.4, e12.4)') 'traj',mcmoves, k, xpq4(0,k), xpq6(0,k),ener(0,k), dh(0,k) 
              endif
              enddo
           endif 


!write(*,*) 'rliold', mcmoves, rliold(:)
rlinew(:)=0
scrivi=0
!!!!!!!!!!!!! CALCOLO Q4, Q6, E PER LE TRAIETTORIa finale
norm0(:)=0

   mcmoves=mcmoves+1
 enddo

   !    do k=0,totiter
   !         do l=0,NbClones-1   
   !             do i=1,N-1
   !               do j=i+1,N
   !                 xij(i,j)=(P(l,k)%qx(i-1)-P(l,k)%qx(j-1))
   !                 yij(i,j)=(P(l,k)%qy(i-1)-P(l,k)%qy(j-1))
   !                 zij(i,j)=(P(l,k)%qz(i-1)-P(l,k)%qz(j-1))
   !                 dij(i,j)=sqrt(xij(i,j)**2+yij(i,j)**2+zij(i,j)**2)
   !                 xij(j,i)=-xij(i,j)
    !                yij(j,i)=-yij(i,j)
    !!                zij(j,i)=-zij(i,j)
     !               dij(j,i)= dij(i,j)
     !             enddo
   !             enddo

     !       xpq4(L,k)=q4(N,xij,yij,zij,dij)
     !       xpq6(L,k)=q6(N,xij,yij,zij,dij)
     !       call Force(P(l,k)%qx,P(l,k)%qy,P(l,k)%qz,fp,Utot)
     !       ener(l,k)=utot
     !      enddo
     !   enddo
   
                do k=0,totiter-1
                   !do i=0, Nbclones-1
                   if (mod(k,100).eq.0) then
                     write(*,'(a,i, i, e11.4, e12.4, e12.4, e12.4)') 'trajclonefin',i, k,xpq4(0,k),xpq6(0,k),ener(0,k),dh(0,k)
                    endif
                   !enddo
               enddo

 write(*,*) 'tau', tau, real(acc)/real(totalmcmoves)
 !tot=(NbClones)*(TotalTime/dt) !

!if (time.gt.2000) then

!endif

!if(mod(iter,1000).eq.0) then
!write(29,*) time
!do m=0,NbClones-1
!write(29,*) m
!   write (29,*) P(m)%qx(:)
!   write (29,*) P(m)%qy(:)
!   write (29,*) P(m)%qz(:)
!   write (29,*) 
!   write (29,*) P(m)%px(:)
!   write (29,*) P(m)%py(:)
!   write (29,*) P(m)%pz(:)      
!   write (29,*)
!enddo 

!endif


!!!!!ADESSO RINORMALIZZO IL VETTORE USCENTE DAL CLONAGGIO, CHE SARÀ PRONTO PER LA PROPAGAZIONE AL TEMPO t+dt

!do m=0,NbClones-1
!   call Normalizza(P(m)%uqx,P(m)%uqy,P(m)%uqz,P(m)%upx,P(m)%upy,P(m)%upz, N)
   !call CalculNorme(P(m)%uqx,P(m)%uqy,P(m)%uqz,P(m)%upx,P(m)%upy,P(m)%upz,N,norm)
  ! write(*,*) 'renorm' , norm
!  enddo

!  mu += log(weight/(RealNumber) NbClones)

!!!TEMPS FINI.RIPARTE IL CICLO



!!!!ALLA FINE FACCIO STATISTICA COMPLESSIVA DELLO SPAZIO DELLE FASI ESPLORATO DAI CLONI

  
!write (*,*) 'usc', nbclones  
!do while(Time.gt.TimeStore)
    !  Lyapunov=0
     ! do l=0,NbClones-1
	!do j=0,N-1
	 ! write(11,*) Time, P(l)%qx(j),P(l)%qy(j), P(l)%qz(j), P(l)%px(j),P(l)%py(j), P(l)%pz(j), P(l)%uqx(j),P(l)%uqy(j),P(l)%uqz(j),P(l)%upx(j), P(l)%upy(j),P(l)%upz(j)
	!Lyapunov+=P(l).Lyap/(RealNumber) NbClones
     !   enddo
     ! enddo
      !write(sortie,'#Lyap\t%lg\t%lg\n',Time,P(4).Lyap/Time)
     ! write(sortie,*)'#Average\t%lg\t%lg\t%lg\n',Time,Lyapunov/Time,mu/Time
     
     ! TimeStore=TimeStore+inter
 !enddo
  
  !Lyapunov=0
  !for(l=0l.lt.NbClonesl++)
   ! Lyapunov+=P(l).Lyap/(RealNumber) NbClones
!  write(sortie,*)'#Final\t%lg\t%lg\t%lg\n'Time !,Lyapunov/Time,mu/Time)


stop

end subroutine LWDTPS


 Subroutine LyapLanczos (xp)

implicit none

real (double), dimension (:), allocatable :: Weight
integer, dimension (:),allocatable ::ToCopy
integer, dimension (:),allocatable ::ToDelete
integer, dimension (:),allocatable ::Tofollow
integer, dimension (:),allocatable ::Nb
integer, dimension (:),allocatable ::Number
integer, parameter ::N=38 ! numero atomes du cluster
 
real (double) ::r     !distance interatomique
real (double) ::eta    !paramètre d'énergie du potentiel
real (double) ::sigma ! dist. éq. pour le potentiel LJ 
real (double),dimension(0:N-1) :: fxx
real (double),dimension(0:N-1) :: fyy
real(double),dimension(0:N-1) :: fzz
real(double),dimension(0:N-1) :: vsecx
real (double),dimension(0:N-1) :: vsecy
real (double),dimension(0:N-1) :: vsecz

real(double),dimension (N,N)::   xij
real (double),dimension (N,N)::   yij
real(double),dimension (N,N)::   zij
real (double),dimension (N,N)::   dij

real (double),dimension (:),allocatable::w4,w6,wp4,wp6
real (double),dimension (:,:),allocatable::xpq4,xpq6,ener
real (double),dimension (:,:),allocatable::correl_2
integer, dimension (:),allocatable ::Nvite
integer, dimension(:), allocatable ::multip

real (double),dimension (N):: fq4x,fq4y,fq4z,fq6x, fq6y,fq6z

real (double), dimension (3,N):: xp
real (double), dimension (3,N):: vp
real (double), dimension (3,N):: fp
real (double), dimension (3,N)::fppro
real (double), dimension (3,N):: xppro
real (double), dimension(1:N):: delt
real (double), dimension(1:N):: xref
real (double), dimension(1:N):: yref
real (double), dimension(1:N):: zref
real (double), dimension(1:N):: gQ6x
real (double), dimension(1:N):: gQ6y
real (double), dimension(1:N):: gq6z
real (double), dimension(1:N):: gEx
real (double), dimension(1:N):: gEy
real (double), dimension(1:N):: gEz

real (double) :: q4,q6, q4moy, q6moy,w4moy, w6moy! linear displacement

type Trajectoire
    real (double), dimension (0:N-1):: qx  
    real (double), dimension (0:N-1):: qy 
    real (double), dimension (0:N-1):: qz  ! vector of position
    
    real (double), dimension (0:N-1):: px  
    real (double), dimension (0:N-1):: py 
    real (double), dimension (0:N-1):: pz  ! vector of impulsion
   
    real(double), dimension(3*N) :: project
    real (double)  Lyap
 endtype

type (Trajectoire), dimension (:,:), allocatable :: P
type (Trajectoire), dimension (:,:), allocatable :: Pfw
type (Trajectoire), dimension (:,:), allocatable :: Pbw
type (Trajectoire), dimension (:,:), allocatable :: Pshift

real (double) ::mcconf,ranf, timebw, timefw,dp,l0, alpha,tempo,teq,z, a_sto
real (double) ::seed, utot, tot, x, U,x4,tconf, tconftot,tconfmoy,kinetotdt, upro, q4ref
real(double)::fftot, fftotdt,laputot,laputotdt,ff,lapu, pi2, uproject, dqx, dqy, dqz, dpx, dpy, dpz,n_react
real (double) ::inter,AverageTime,NbTimeStep,dt,TotalTime,Time,w,Norm,TimeStore,c1,c2,sqrtgTh,eta1,eta2,gamma,eps, Temperature ! Lyapunov=0,LyapunovAv=0,mu=0,
integer :: nbClones, j,l,k,ran,m, a,i,b, cont,ien,nq4,n0, n4, ltot,mcmoves,totalmcmoves, newtraj,tsecond
logical ::  new_projection
logical ::  waste_recycling
logical ::  FCC
logical ::  defaut
character(len=128) :: sortie
!character(len=128) :: histotot
!character(len=128) :: histopart
character(len=128) ::posfinal
!character(len=128) ::snap
!character(len=128) ::stock
!character(len=128) ::stq4t
!character(len=128) ::stq6t
!character(len=128) ::q4q6
character(len=128) ::data_mbar
character(len=128) ::moyennes_mbar
character(len=128) ::moyennes_mbar_denom

character(len=128) ::kappaF
character(len=128) ::kappaFd
character(len=128) ::kappadI
character(len=128) ::kappaFI
character(len=128) ::kappadF

real (double)::energy,rga,beq,rien,q,kine,kinetot,pif, trasm1,trasm2, normgradq4, normgradq6, q6ref,phi, theta, pi,rlimed,fhi
real (double) ::xtemp0,xtemp1,xtemp2,xtemp3,xtemp4,xtemp5,ss,delta_x, kapa
integer, parameter::nfenetre=100
integer :: ienergy, ifenetre, ifen_q6, iter,kaq4,kaq6,niter,kappa,ix, scrivi, continue

real (double):: stat(0:nfenetre, 0:nfenetre)
real (double):: statneg(0:nfenetre, 0:nfenetre)
real (double):: statpos(0:nfenetre, 0:nfenetre)
real (double):: stat6(0:nfenetre, 0:nfenetre)
real(double)::tvmoytot,xn,enprmoy,jtrans, jtot,trasm,trastot

real(double) :: tab_contour(-10:nfenetre+10,0:nfenetre)
real(double) :: tab_contour_q6(-10:nfenetre+10,0:nfenetre)
real(double) :: tab_cont_q4q6(-10:nfenetre+10,-10:nfenetre+10)
real(double) :: cumul_contour(-10:nfenetre+10,0:nfenetre)
real(double) :: cumul_contour_q6(-10:nfenetre+10,0:nfenetre)
real(double) :: cumul_contour_q4q6(-10:nfenetre+10,-10:nfenetre+10)
real(double) :: tab_energy(-10:nfenetre+10,2),toto
real (double), dimension (:), allocatable::norm0

real (double)::x111, xmaxq6,energie_min,energie_max,ecinetique,tequilib,delta,beqtot,dhtot,beqtotdt,dhtotdt,intbeq, intdh,e,pix, timebsh, timefsh, alpha_current
integer::icfen,icener,icfen_q6,a1, a2, a3, a4, a5,iterfw, iterbw,nmax,totiter,  lo, tprim, nbclones_mbar, depart_boucle_nbclones, tprimo

real(double), dimension(:), allocatable::betaq,enerpro,poids,alpha_bias, q4bw, q6bw
real (double) ,dimension(:), allocatable:: ener0, xpq4ref, xpq6ref
real (double) ,dimension(:), allocatable:: ener0bw
real (double) ,dimension(:), allocatable::newLyap
real (double) ,dimension(:), allocatable::oldlyap 
real (double) ,dimension(:), allocatable::rapporto
real (double) ,dimension(:), allocatable::react
real (double) ,dimension(:), allocatable::reactivity
real (double) ,dimension(:), allocatable::h_A
real (double) ,dimension(:), allocatable::h_A_average
real (double) ,dimension(:), allocatable::tau
real (double) ,dimension(:), allocatable:: hamilt
real (double) ,dimension(:,:), allocatable:: dh
real (double) ,dimension(:,:), allocatable:: dhx
real (double) ,dimension(:,:), allocatable:: S
real (double) ,dimension(:,:,:), allocatable:: U_KNL
integer ,dimension(:), allocatable::acc

real (double) ,dimension(:), allocatable::h_F
real (double) ,dimension(:), allocatable::react_F
real (double) ,dimension(:), allocatable::h_Fd
real (double) ,dimension(:), allocatable::react_Fd
real (double) ,dimension(:), allocatable::h_dI
real (double) ,dimension(:), allocatable::react_dI
real (double) ,dimension(:), allocatable::h_dF
real (double) ,dimension(:), allocatable::react_dF
real (double) ,dimension(:), allocatable::h_FI
real (double) ,dimension(:), allocatable::react_FI


write(6,*)'usage %s:\n'

!    Initialisation of the variables
write(6,*)'TotalMcmoves '
read (*,*)
read (*,*) Totalmcmoves      ! Total duration of the simulation 
write(*,*) Totalmcmoves 
read (*,*)
write(6,*)'TotalTime '
read (*,*) TotalTime      ! Total duration of the simulation 
write(*,*) TotalTime 
read (*,*)
write(6,*)'dt' 
read (*,*) dt          ! Time step
write(*,*) dt
read (*,*)
write(6,*)'NbClones or current clone'  
read (*,*) NbClones  
write(*,*) NbClones    ! Number of clones
read (*,*)
write(6,*)'temperature'  
read (*,*) temperature  ! temperature (kT)
write (*,*) temperature
read (*,*)
write(6,*)'gamma*dt'  
read (*,*) gamma  ! friction*dt
write (*,*) gamma       ! interval between 2 data records
read (*,*)
write(6,*)'sortie data_mbar'  
read (*,*) sortie   
write (*,*) sortie! Name of the file where data are stored
read (*,*)
write(6,*)'alpha max'
read (*,*) alpha 
write (*,*) alpha
read (*,*)
write(6,*)'sigmap'
read (*,*) ss
write (*,*) ss
read (*,*)
write(6,*)'niter'
read (*,*) niter
write (*,*) niter
read (*,*)
write(6,*)'t_equilib'
read (*,*) teq
write (*,*) teq
read (*,*)
write(6,*)'delta_X'
read (*,*) delta_x
write (*,*) delta_x
read (*,*)
write(6,*)'a_sto'
read (*,*) a_sto
write (*,*) a_sto
read (*,*)
write(6,*)'k ressort'
read (*,*) kapa
write (*,*) kapa
read (*,*)
write(6,*)'continue'
read (*,*) continue
write (*,*) continue
read (*,*)
write(6,*)'nbclones_mbar'
read (*,*) nbclones_mbar
write (*,*) nbclones_mbar
read (*,*)
write(6,*)'cluster? 0 no, nbclones yes'
read (*,*) depart_boucle_nbclones
write (*,*) depart_boucle_nbclones
read (*,*)
write(6,*)'waste recycling?'
read (*,*) waste_recycling
write (*,*) waste_recycling
read (*,*)
write(6,*)'ressort on FCC?'
read (*,*) fcc
write (*,*) fcc
read (*,*)
write(6,*)'ressort on defaut?'
read (*,*) defaut
write (*,*) defaut
read (*,*)
write(6,*)'tprimo (int) valutazione k_react def-ico e def-fcc'
read (*,*) tprimo
write (*,*) tprimo

open (unit=11, FILE=sortie, action='write', status='replace')
!open (unit=31, action='write', status='replace')
moyennes_mbar  =sortie(1:lenfnam)//'.data_moy'
moyennes_mbar_denom  =sortie(1:lenfnam)//'.data_moy2'
data_mbar  =sortie(1:lenfnam)//'.data'
!histotot  =sortie(1:lenfnam)//'.thout'
!histopart =sortie(1:lenfnam)// '.phout'
posfinal = sortie(1:lenfnam)//'.cin'

kappaF = sortie(1:lenfnam)//'.kappaF'
kappaFd = sortie(1:lenfnam)//'.kappaFd'
kappadI = sortie(1:lenfnam)//'.kappadI'
kappaFI= sortie(1:lenfnam)//'.kappaFI'
kappadF= sortie(1:lenfnam)//'.kappadF'

!snap = sortie(1:lenfnam)//'.film'
!stock = sortie(1:lenfnam)//'.stock'
!stq4t = sortie(1:lenfnam)//'.stq4t'
!stq6t= sortie(1:lenfnam)//'.stq6t'
!q4q6= sortie(1:lenfnam)//'.qhout'

!open(unit=14, file=histopart, action='write', status='replace')
!open(unit=26, file=histotot,  action='write', status='replace')
!open(unit=28, file=snap,  action='write', status='replace')
!open(unit=29, file=stock,  action='write', status='replace')

!!!!!!!!!! output per calcolo correl function reactivity
open(unit=30, file=kappaF,  action='write', status='replace')
open(unit=31, file=kappaFd,  action='write', status='replace')
open(unit=32, file=kappadI,  action='write', status='replace')
open(unit=33, file=kappaFI,  action='write', status='replace')
open(unit=34, file=kappadF,  action='write', status='replace')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
open(unit=24, file=data_mbar, action='write', status='replace')
open(unit=25, file=moyennes_mbar, action='write', status='replace')
open(unit=23, file=moyennes_mbar_denom, action='write', status='replace')
NMax=NbClones*3.0+1    ! Max number of clones
!TimeStore = TotalTime-AverageTime ! Next time at which data will be stored 
                                           
! Array where the weight of every clones is stored
allocate (acc(0:nmax))
allocate (tau(0:nmax))
allocate (h_A_average(0:nmax))
allocate (react(0:nmax))
allocate (betaq(0:nmax-1))
allocate (ener0(0:Nmax))
allocate (ener0bw(0:Nmax))
allocate(enerpro(0:nmax-1))
allocate(hamilt(0:nmax-1))
! Array which contains the numbers from 1 to Nmax
allocate(Nb(0:NMax-1)) 
allocate(norm0(0:NMax-1)) 
allocate (q4bw(0:nmax))
allocate (q6bw(0:nmax))
kappa=1
scrivi=0
ltot=int(totaltime/real(kappa*dt))
totiter=nint(totaltime/dt)
write(*,*)'totiter', totiter
gamma=gamma/(dt)
write(*,*)'gamma' , gamma, gamma*dt,ltot
!waste_recycling=.true.
!waste_recycling=.false.
write(*,*) 'waste_recycling?',  waste_recycling
allocate (h_A(0:2*totiter+10))
allocate (reactivity(0:2*totiter+10))
allocate (correl_2(0:2*totiter+10,0:2*totiter+1))
allocate (xpq4(0:Nmax,0:2*totiter+10))
allocate (xpq6(0:Nmax,0:2*totiter+10))
allocate (xpq4ref(0:Nmax))
allocate (xpq6ref(0:Nmax))
allocate(ener(0:nmax-1,0:totiter+10))
allocate(S(0:nmax-1,0:2*totiter+10))
allocate(u_knl(0:nmax-1,0:nmax-1,0:Totalmcmoves))
allocate(old_projection(3*N))
!allocate(projection(3*N))
allocate(first_projection(3*N))
allocate(old_before_sc_projection(3*N))
!allocate(project(0:totiter+10,3*N))

allocate (h_F(0:2*totiter+10))
allocate (react_F(0:nmax))
allocate (h_Fd(0:2*totiter+10))
allocate (react_Fd(0:nmax))
allocate (h_dF(0:2*totiter+10))
allocate (react_dF(0:nmax))
allocate (h_dI(0:2*totiter+10))
allocate (react_dI(0:nmax))
allocate (h_FI(0:2*totiter+10))
allocate (react_FI(0:nmax))

ener(:,:)=0
stat(:,:)=0
stat6(:,:)=0
sigma=1
eta=1
enerpro(:)=0
xpq4(:,:)=0
xpq6(:,:)=0
!statneg(:,:)=0
!statpos(:,:)=0
hamilt(:)=0
kine=0
laputot=0
fftot=0
laputotdt=0
fftotdt=0
kinetotdt=0
kinetot=0
mcmoves=0
tau(:)=0
enprmoy=0
nq4=0
e=0.000001
eigenvalue=0
!ss=0.05
l0=0.0
acc(:)=0
n_react=0.0
x111   = 0.2
xmaxq6 = 0.6
energie_min = -174.0
energie_max = -144.0

     tab_contour(-10:nfenetre+10,0:nfenetre)=zero
     tab_contour_q6(-10:nfenetre+10,0:nfenetre)=zero
     tab_cont_q4q6(-10:nfenetre+10,-10:nfenetre+10)=zero
     cumul_contour_q4q6(-10:nfenetre+10,-10:nfenetre+10)=0
     cumul_contour(-10:nfenetre+10,0:nfenetre)=0
     cumul_contour_q6(-10:nfenetre+10,0:nfenetre)=0

allocate(Number(0:NbClones-1)) 

!!!!!!!!!!!!T EQUILIBRAGE
tequilib=(dt)*teq
write(*,*) 'tquilib',tequilib
!!!!!!!!!!!! T EQUILIBRAGE
  
write(*,*) 'sortie affichée'

allocate(P(0:NMax-1,0:totiter+10)) 
allocate(Pfw(0:NMax-1,0:totiter+10)) 
allocate(Pbw(0:NMax-1,0:totiter+10))
allocate(Pshift(0:NMax-1,0:2*totiter+10))

allocate (dh(0:nmax-1,0:2*totiter+10))
allocate(oldLyap(0:nmax))
allocate(newLyap(0:nmax))
allocate(rapporto(0:nmax))
allocate (dhx(0:nmax-1,0:2*totiter+10))
allocate(poids(1:totiter))
allocate (alpha_bias(0:nbclones_mbar))!!EQUILIBRAGE
xppro(:,:)=0
rapporto(:)=0
iter=0

rien=0
Q=0
uproject=0
kappa=1
pi=3.14159
if ((gamma*dt).le.100000) then
 
   rga = exp(-gamma*dt/two)
else 
  rga=0
endif

write(*,*) 'rga',gamma, rga,two,pi

if(continue.eq.1) then
  open(unit=27, file=posfinal, status='old')
  do j=1,N
  read (27,*) xp(1,j), xp(2,j), xp(3,j)
  enddo
  close(27)
endif

xref(:)=xp(1,:)
yref(:)=xp(2,:)
zref(:)=xp(3,:)

tempo=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! INITIALISATION !!!!!!!!!!!!!!!!!!!!!!
write(*,*) 'conditions initiales pour traj de reference avc distrib stoch à temperature T=',temperature
write(*,*) 'temps equilibrage', tequilib
!!!!!!!!!!!!!!!!!!!!!!!! EQUILIBRAGE initial (STOCH DYN)!!!!!!!!!!!!!!!!!!
     do while ((tempo.lt.Tequilib))!.and.(time.lt.totaltime))
            call calfoljc2(xp,fp,utot)
            call langevin2(xp,vp,fp,dt,rga, temperature,Q,rien,ecinetique)
           ! ff=sum(fp(1,:)**2+fp(2,:)**2+fp(3,:)**2)
            !do i=1,N-1
            !     do j=i+1,N
            !         xij(i,j)=xp(1,i)-xp(1,j)
            !         yij(i,j)=xp(2,i)-xp(2,j)
            !         zij(i,j)=xp(3,i)-xp(3,j)
            !         dij(i,j)=sqrt(xij(i,j)**2+yij(i,j)**2+zij(i,j)**2)
            !         xij(j,i)=-xij(i,j)
            !        yij(j,i)=-yij(i,j)
            !        zij(j,i)=-zij(i,j)
            !        dij(j,i)= dij(i,j)
            !      enddo
            !  enddo
              !x4=q4(N,xij,yij,zij,dij)
           tempo=tempo+dt
       enddo
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

write(*,*) 'point de départ traj de reference deterministe' 
     do  m = 0,NbClones 
       P(m,0)%qx=xp(1,:)
       P(m,0)%qy=xp(2,:)
       P(m,0)%qz=xp(3,:)
     enddo

    do m = 0,NbClones
       P(m,0)%px=vp(1,1:N)
       P(m,0)%py=vp(2,1:N)
       P(m,0)%pz=vp(3,1:N)
    enddo

!write(*,*) 'ctrel' , P(0,0)%qx, P(1,0)%QX
     !   do j=0, Nbclones-1
     !          k=0
     !          do i=1,N-1
     !             do a=i+1,N
     !               xij(i,a)=(P(j,k)%qx(i-1)-P(j,k)%qx(a-1))
     !               yij(i,a)=(P(j,k)%qy(i-1)-P(j,k)%qy(a-1))
     !               zij(i,a)=(P(j,k)%qz(i-1)-P(j,k)%qz(a-1))
     !               dij(i,a)=sqrt(xij(i,a)**2+yij(i,a)**2+zij(i,a)**2)
     !               xij(a,i)=-xij(i,a)
     !               yij(a,i)=-yij(i,a)
     !               zij(a,i)=-zij(i,a)
      !              dij(a,i)= dij(i,a)
      !            enddo
      !!          enddo!

      !      xpq4(j,k)=q4(N,xij,yij,zij,dij)
      !      xpq6(j,k)=q6(N,xij,yij,zij,dij)
       !     call Force(P(j,k)%qx,P(j,k)%qy,P(j,k)%qz,fp,Utot)
      !      ener(j,k)=utot
            !if (mod(k,100).eq.0) then
       !     write(*,'(a,i, i, e11.4, e12.4, e12.4)') 'trajinit0',j, k,xpq4(j,k),xpq6(j,k),ener(j,k)
           ! endif      
       !  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! metto a zero il resto avanzato
 P(:,:)%Lyap = 0

 do l=1,totiter 
  do m = NbClones+1,NMax-1
    do j=0,N-1
    P(m,l)%qx(j) =0
    P(m,l)%px(j) =0
  
    P(m,l)%qy(j) =0
    P(m,l)%py(j) =0
    
    P(m,l)%qz(j) =0
    P(m,l)%pz(j) =0
    enddo
   ! P(m).Lyap = 0
  enddo
enddo

write(*,*) 'initialisation faite: go with Lanczos'
ecinetique=0
dh(:,:)=0
!write(*,*) 'initialisation faite'

!! initailisation paramètres de bias alfa pour reconstruction
do j=0, nbclones_mbar
alpha_bias(j)=real(j*(alpha/real(nbclones_mbar)))
enddo 
write(*,*)'bias', alpha_bias(:)

!do while (mcmoves.lt.Totalmcmoves)
!do j=0,nbclones
!write (*,*)'boh',  alpha_bias(j)
!enddo
!mcmoves=mcmoves+1
!enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  c'est parti !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!write(*,*) 'initialisation faite'
j=0
!j=depart_boucle_nbclones


new_projection=.true.

      time=0
     do while (Time.lt.Totaltime) !!!!!!!!!!!!!!!!!!!!!!  propagazione prima traiettoria   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            iter=nint(real(time)/real(dt))
            call propagatfwdet(P(j,iter)%qx,P(j,iter)%qy,P(j,iter)%qz,P(j,iter)%px,P(j,iter)%py,P(j,iter)%pz,dt,N,dh(j,iter))
            call lanczos(N,P(j,iter)%qx,P(j,iter)%qy,P(j,iter)%qz,new_projection,eigenvalue,P(j,iter)%project)!!!posizioni date dopo la propagazione
               !write(*,*)'proj', iter, P(j,iter)%project(1)
              !if (iter.eq.0) then 
                 if (eigenvalue.lt.0.0) then 
                   P(j,iter)%Lyap=log(real(1.0+dt*sqrt(abs(eigenvalue))))
                 else 
                   P(j,iter)%Lyap=0.0  !! autoval max del Lanczos
                 endif
           !!!write (*,*) 'lyap', j,iter,  P(j,iter)%Lyap
              !else
              !   if (eigenvalue.le.0.0) then 
              !     P(j,iter)%Lyap =P(j,iter-1)%Lyap+log(real(1.0+dt*sqrt(abs(eigenvalue))))
              !   else 
              !    P(j,iter)%Lyap=P(j,iter-1)%Lyap+0.0  !! autoval max del Lanczos
              !   endif   ! autoval max del Lanczos
           !write (*,*) 'lyap', j, iter, P(j,iter)%Lyap/real(iter)
            !endif
          !  write (*,*) 'lyap', j, iter, P(j,iter)%Lyap, eigenvalue!/real(iter)
          !  write (*,*) 'dh', mcmoves, dh(j,iter)
            !endif
            P(j,iter+1)%qx=P(j,iter)%qx
            P(j,iter+1)%qy=P(j,iter)%qy
            P(j,iter+1)%qz=P(j,iter)%qz

            P(j,iter+1)%px=P(j,iter)%px
            P(j,iter+1)%py=P(j,iter)%py
            P(j,iter+1)%pz=P(j,iter)%pz
            P(j,iter+1)%project=P(j,iter)%project
       time=time+dt
       new_projection=.false.
   enddo  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !write(*,*) 'lyap',P(j,totiter-1)%Lyap/real(totiter-1)

!oldLyap(j)=sum(P(j,:)%Lyap)/real(totiter)
!write(*,*) 'lyap', oldLyap(j)
!stop
        call Force(P(j,0)%qx,P(j,0)%qy,P(j,0)%qz,fp,Utot)
       !write (*,*) 'yy', utot
        ener0(j)=dh(j,totiter-1)
        !write (*,*) 'yy', ener0(j)
           do k=0,totiter-1
               do i=1,N-1
                  do a=i+1,N
                    xij(i,a)=(P(j,k)%qx(i-1)-P(j,k)%qx(a-1))
                    yij(i,a)=(P(j,k)%qy(i-1)-P(j,k)%qy(a-1))
                    zij(i,a)=(P(j,k)%qz(i-1)-P(j,k)%qz(a-1))
                    dij(i,a)=sqrt(xij(i,a)**2+yij(i,a)**2+zij(i,a)**2)
                    xij(a,i)=-xij(i,a)
                    yij(a,i)=-yij(i,a)
                    zij(a,i)=-zij(i,a)
                    dij(a,i)= dij(i,a)
                  enddo
                enddo

            xpq4(j,k)=q4(N,xij,yij,zij,dij)
            xpq6(j,k)=q6(N,xij,yij,zij,dij)
            call Force(P(j,k)%qx,P(j,k)%qy,P(j,k)%qz,fp,Utot)
            ener(j,k)=utot
           !!!!!!!!!! traiettorie reattive !!!!!!!!!!!!!!!!
          ! if ((xpq4(j,k).ge.0.1).and.(xpq4(j,k).ge.0.3)) then
          ! tau=tau+1
          ! endif
             
              if (mod(k,100).eq.0) then
              write(*,'(a,i, i, i, e11.4, e12.4, e12.4, e12.4)') 'traj0',mcmoves, j, k, xpq4(j,k), xpq6(j,k),ener(j,k),dh(j,k)/real(3*38)
              endif
          enddo 

          xpq4ref(j)=xpq4(j,0)
          xpq6ref(j)=xpq6(j,0)

!write(*,*) 'qref', xpq4ref(j), xpq6ref(j)
  

do j=depart_boucle_nbclones,NbClones

P(j,:)=P(0,:)
ener0(:)=ener0(0)
oldLyap(j)=sum(P(0,:)%Lyap)/real(totiter)

xpq4ref(j)=xpq4ref(0)
xpq6ref(j)=xpq6ref(0)

enddo




!stop

do while (mcmoves.lt.Totalmcmoves) !!!!!!!!!!!!! inizio montecarlo su traiettorie
!!!!!!!!!!!!  fACCIAMO SHOOTING:
  ix=nint((totiter-1)*genrand()) ! PUNTO DI SHOOTING ltot numero di punti stockati (segmenti di traiettoria)
  timefw=ix*dt!*kappa
  timebw=ix*dt!*kappa
 !write(*,*) 'shooting time', ix, timefw, timebw

!!!! shooting deterministico 
  !call gauss(ss,l0,dp)
  dh(:,:)=0
  !write (*,*) 'dp', dp
!!!per lo shooting stocastico perturbo anche le posizioni ....

do j=depart_boucle_nbclones,nbclones
!write(*,*) 'qref2', xpq4ref(j), xpq6ref(j)
 !!!!! inizializzo backward 
 !if(mod(j,2).eq.0) then
     Pbw(j,ix)%qx=P(j,ix)%qx
     Pbw(j,ix)%qy=P(j,ix)%qy
     Pbw(j,ix)%qz=P(j,ix)%qz
     Pbw(j,ix)%project=P(j,ix)%project

!!!! shooting classique
    
   do m=0, N-1

   call gauss(ss,l0,dpx)
   call gauss(ss,l0,dpy)
   call gauss(ss,l0,dpz)
!!!!! shooting classique
    !Pbw(j,ix)%px(m)=P(j,ix)%px(m)+dpx
    !Pbw(j,ix)%py(m)=P(j,ix)%py(m)+dpy
    !Pbw(j,ix)%pz(m)=P(j,ix)%pz(m)+dpz
    !enddo

!!!! shooting 'Stoltz'
    Pbw(j,ix)%px(m)=a_sto*P(j,ix)%px(m)+sqrt(1.0-a_sto*a_sto)*dpx 
    Pbw(j,ix)%py(m)=a_sto*P(j,ix)%py(m)+sqrt(1.0-a_sto*a_sto)*dpy
    Pbw(j,ix)%pz(m)=a_sto*P(j,ix)%pz(m)+sqrt(1.0-a_sto*a_sto)*dpz
    enddo

 !else 
 !    Pbw(j,ix)%qx=P(j-1,ix)%qx
 !    Pbw(j,ix)%qy=P(j-1,ix)%qy
 !    Pbw(j,ix)%qz=P(j-1,ix)%qz
 !    Pbw(j,ix)%project=P(j-1,ix)%project

   !  do m=0, N-1!

   !  call gauss(ss,l0,dpx)
   !  call gauss(ss,l0,dpy)
   !  call gauss(ss,l0,dpz)

!!!!! shooting classique
!     Pbw(j,ix)%px(m)=P(j-1,ix)%px(m)+dpx
!     Pbw(j,ix)%py(m)=P(j-1,ix)%py(m)+dpy
!     Pbw(j,ix)%pz(m)=P(j-1,ix)%pz(m)+dpz
!     enddo

!!!! shooting 'Stoltz'
    ! Pbw(j,ix)%px(m)=a_sto*P(j-1,ix)%px(m)+sqrt(1.0-a_sto*a_sto)*dpx 
    ! Pbw(j,ix)%py(m)=a_sto*P(j-1,ix)%py(m)+sqrt(1.0-a_sto*a_sto)*dpy
    ! Pbw(j,ix)%pz(m)=a_sto*P(j-1,ix)%pz(m)+sqrt(1.0-a_sto*a_sto)*dpz
    ! enddo

 !endif

timebw=ix*dt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! propagazione BACKWARD a partire da punto di shooting (senza lyapunov)

new_projection=.false.
!projection=project(ix,:)
iterbw=nint(timebw/dt)
  do while (iterbw.gt.0)
      ! iterbw=nint(timebw/dt)
       call propagatbckw(Pbw(j,iterbw)%qx,Pbw(j,iterbw)%qy,Pbw(j,iterbw)%qz,Pbw(j,iterbw)%px,Pbw(j,iterbw)%py,&
                         Pbw(j,iterbw)%pz,rga,temperature,dt,N,dh(j,iterbw))
       call lanczos(N,Pbw(j,iterbw)%qx,Pbw(j,iterbw)%qy,Pbw(j,iterbw)%qz,new_projection,eigenvalue, Pbw(j,iterbw)%project)
           ! write(*,*) 'putt', timebw,iterbw-1
             ! write(*,*) 'dh' , mcmoves, dh(j,iterbw)
          !if (iterbw.eq.ix) then 
                if (eigenvalue.lt.0.0) then 
                  Pbw(j,iterbw-1)%Lyap=log(real(1.0+dt*sqrt(abs(eigenvalue))))
                else 
                  Pbw(j,iterbw-1)%Lyap=0.0  !! autoval max del Lanczos
                endif
             !write (*,*) 'bw',Pbw(j,iterbw)%qx(1)
          !else
          !      if (eigenvalue.le.0.0) then 
          !        Pbw(j,iterbw)%Lyap =Pbw(j,iterbw+1)%Lyap+log(real(1.0+sqrt(abs(eigenvalue))))
          !      else 
          !        Pbw(j,iterbw)%Lyap=Pbw(j,iterbw+1)%Lyap+0.0  !! autoval max del Lanczos
          !      endif   ! autoval max del Lanczos
    ! write (*,*) 'lyapbw',j, iterbw-1, Pbw(j,iterbw-1)%Lyap,eigenvalue!/real(iterbw)
          !endif
      Pbw(j,iterbw-1)%qx=Pbw(j,iterbw)%qx
      Pbw(j,iterbw-1)%qy=Pbw(j,iterbw)%qy
      Pbw(j,iterbw-1)%qz=Pbw(j,iterbw)%qz

      Pbw(j,iterbw-1)%Px=Pbw(j,iterbw)%Px
      Pbw(j,iterbw-1)%Py=Pbw(j,iterbw)%Py
      Pbw(j,iterbw-1)%Pz=Pbw(j,iterbw)%Pz
      Pbw(j,iterbw-1)%Project=Pbw(j,iterbw)%Project
      iterbw=iterbw-1
   ! timebw=timebw-dt
     new_projection=.false.
 enddo

         
               do i=1,N-1
                  do a=i+1,N
                    xij(i,a)=(Pbw(j,iterbw)%qx(i-1)-Pbw(j,iterbw)%qx(a-1))
                    yij(i,a)=(Pbw(j,iterbw)%qy(i-1)-Pbw(j,iterbw)%qy(a-1))
                    zij(i,a)=(Pbw(j,iterbw)%qz(i-1)-Pbw(j,iterbw)%qz(a-1))
                    dij(i,a)=sqrt(xij(i,a)**2+yij(i,a)**2+zij(i,a)**2)
                    xij(a,i)=-xij(i,a)
                    yij(a,i)=-yij(i,a)
                    zij(a,i)=-zij(i,a)
                    dij(a,i)= dij(i,a)
                  enddo
                enddo
            q4bw(j)=q4(N,xij,yij,zij,dij)
            q6bw(j)=q6(N,xij,yij,zij,dij)
     
!write (*,*) 'gna' , iterbw, q4bw(j), q6bw(j)

!write(*,*) 'qref21', xpq4ref(j), xpq6ref(j)
!stop
  !!!!!!!!!!!!!!!!!!!!!!!!!!! calcolo energia di config iniziale traiettoria BW
 ! call Force(Pbw(j,0)%qx,Pbw(j,0)%qy,Pbw(j,0)%qz,fp,Utot)
  ener0bw(j)=dh(j,ix)!utot
 
!!!!!!!!!!!!!!!!!!!!! la traiettoria bw é arrivata a tempo t=0
  
!do j=0,nbclones-1

     Pfw(j,ix)%qx=Pbw(j,ix)%qx
     Pfw(j,ix)%qy=Pbw(j,ix)%qy
     Pfw(j,ix)%qz=Pbw(j,ix)%qz
     Pfw(j,ix)%project=Pbw(j,ix)%project

     Pfw(j,ix)%Px=Pbw(j,ix)%Px
     Pfw(j,ix)%Py=Pbw(j,ix)%Py
     Pfw(j,ix)%Pz=Pbw(j,ix)%Pz

    Pfw(j,ix)%Lyap=Pbw(j,ix)%Lyap

 !write (*,*) 'lup', Pbw(j,ix-1)%upx 
! write (*,*) 'lupx',Pbw(j,ix)%upx
!write (*,*) 'lyapfcong', ix, Pfw(j,ix-1)%Lyap,Pbw(j,1)%Lyap 
timefw=(ix)*dt

 new_projection=.false.
iterfw=nint(timefw/dt)
!!!!!!!!!!!! trajectoire forward !!!!!!!!!!!!!!!!!!!!
  do while (iterfw.lt.totiter-1)
!do while (timefw.lt.totaltime-dt)
  !        iterfw=nint(timefw/dt)
     ! write(*,*) iterfw+1
                call propagatfwDET(Pfw(j,iterfw)%qx,Pfw(j,iterfw)%qy,Pfw(j,iterfw)%qz,Pfw(j,iterfw)%px,Pfw(j,iterfw)%py,&
                                   Pfw(j,iterfw)%pz,dt,N,dh(j,iterfw))
                call lanczos(N,Pfw(j,iterfw)%qx,Pfw(j,iterfw)%qy,Pfw(j,iterfw)%qz,new_projection,eigenvalue,Pfw(j,iterfw)%project)!!!!!
              !  write(*,*) 'dh' , mcmoves, dh(j,iterfw)!write (*,*) 'fw',iterfw, Pfw(j,iterfw)%project(1)
! write(*,*) 'ctrlfw' , eigenvalue!if (mod(iter,kappa)==0) then
               ! if (iterfw.gt.ix-1) then 
                  if (eigenvalue.le.0.0) then 
                    Pfw(j,iterfw+1)%Lyap =log(real(1.0+dt*sqrt(abs(eigenvalue))))!Pfw(j,iterfw-1)%Lyap+
                  else 
                    Pfw(j,iterfw+1)%Lyap=0.0!Pfw(j,iterfw-1)%Lyap+0.0  !! autoval max del Lanczos
                  endif  
               
               !write (*,*) 'lyapfw',j, iterfw+1, Pfw(j,iterfw+1)%Lyap,eigenvalue!/real(iterfw)
               ! endif
                !write (*,*) 'dhfw', mcmoves, iterfw, j, dh(j,iterfw)
                !endif
                Pfw(j,iterfw+1)%qx=Pfw(j,iterfw)%qx
                Pfw(j,iterfw+1)%qy=Pfw(j,iterfw)%qy
                Pfw(j,iterfw+1)%qz=Pfw(j,iterfw)%qz

                Pfw(j,iterfw+1)%Px=Pfw(j,iterfw)%Px
                Pfw(j,iterfw+1)%Py=Pfw(j,iterfw)%Py
                Pfw(j,iterfw+1)%Pz=Pfw(j,iterfw)%Pz
                Pfw(j,iterfw+1)%Project=Pfw(j,iterfw)%Project
              iterfw=iterfw+1
            !timefw=timefw+dt
       new_projection=.false.
    enddo
             !  do i=1,N-1
             !     do a=i+1,N
             !       xij(i,a)=(Pfw(j,totiter-1)%qx(i-1)-Pfw(j,totiter-1)%qx(a-1))
             !       yij(i,a)=(Pfw(j,totiter-1)%qy(i-1)-Pfw(j,totiter-1)%qy(a-1))
             !       zij(i,a)=(Pfw(j,totiter-1)%qz(i-1)-Pfw(j,totiter-1)%qz(a-1))
             !       dij(i,a)=sqrt(xij(i,a)**2+yij(i,a)**2+zij(i,a)**2)
             !       xij(a,i)=-xij(i,a)
             !       yij(a,i)=-yij(i,a)
             !       zij(a,i)=-zij(i,a)
             !       dij(a,i)= dij(i,a)
              !    enddo
              !  enddo
!write(*,*) 'fw', Pfw(j,totiter-1)%qx(1)
!stop

do k=0,ix-1
newlyap(j)=newlyap(j)+Pbw(j,k)%Lyap
enddo

do k=ix, totiter-1
newlyap(j)=newlyap(j)+Pfw(j,k)%Lyap
enddo

newlyap(j)=newlyap(j)/real(totiter)

!write(*,*) 'qref2a', xpq4ref(j), xpq6ref(j)

!newLyap(j)=(sum(Pfw(j,:)%Lyap)+sum(Pbw(j,:)%Lyap)+P(j,ix)%Lyap)/real(totiter-1)

!write(*,*) 'new' ,newLyap(j)-oldlyap(j)

!write(*,*) 'q4q6', q4(N,xij,yij,zij,dij), q6(N,xij,yij,zij,dij)
!stop

if (defaut) then
!!!!!!!!!!!!!!!!!!!!!!!!! ressort      
!   rapporto(j)=exp(alpha_bias(j)*(newLyap(j)-oldLyap(j)))*exp(-kapa*10.0*(q4bw(j)-0.12)**2-kapa*(q6bw(j)-0.47)**2)*&
!                  exp(kapa*10.0*(xpq4ref(j)-0.12)**2+kapa*(xpq6ref(j)-0.47)**2) !*exp(-(ener0bw(j)-ener0(j))/temperature)

!!!!!!!!!!!!!!!! echelon
  if ((xpq4ref(j).gt.0.1).and.(xpq4ref(j).lt.0.13).and.(q4bw(j).gt.0.1).and.(q4bw(j).lt.0.13)) then
  rapporto(j)= exp(alpha_bias(j)*(newLyap(j)-oldLyap(j)))
  else
  rapporto(j)=0.0
  endif

endif

if (fcc) then
if (alpha_bias(j).ne.0.0) then
   if (fcc) then
      rapporto(j)=exp(alpha_bias(j)*(newLyap(j)-oldLyap(j)))*exp(-kapa*10.0*(q4bw(j)-0.18)**2-kapa*(q6bw(j)-0.57)**2)*&
                  exp(kapa*10.0*(xpq4ref(j)-0.18)**2+kapa*(xpq6ref(j)-0.57)**2)!*exp(-(ener0bw(j)-ener0(j))/temperature)
   else if (defaut) then
      rapporto(j)=exp(alpha_bias(j)*(newLyap(j)-oldLyap(j)))*exp(-kapa*10.0*(q4bw(j)-0.12)**2-kapa*(q6bw(j)-0.47)**2)*&
                  exp(kapa*10.0*(xpq4ref(j)-0.12)**2+kapa*(xpq6ref(j)-0.47)**2) !*exp(-(ener0bw(j)-ener0(j))/temperature)
   endif
else
    rapporto(j)=1.0!*exp(-(ener0bw(j)-ener0(j))/temperature)
endif
endif

!write(*,*) 'rapp' , rapporto(j)!, exp(-kapa*10.0*(q4bw(j)-0.18)**2-kapa*(q6bw(j)-0.57)**2)
!write(*,*) 'qref2', xpq4ref(j), xpq6ref(j)
!q4bw(:)=0.0
!q6bw(:)=0.0
!*exp(-(ener0bw(j)-ener0(j))/temperature)

! write(*,'(a, e12.4, e12.4, e12.4)') 'deltarli', newLyap(j)-oldLyap(j),ener0bw(j)/real(3*38),rapporto(j)!, dh(j,totiter-1)

!write(*,*) 'uno', j,mcmoves!exp(alpha*(RLInew(j)-RLIold(j)))-ener0(j)
!write(*,'(a, i, i, e12.4, e12.4, e12.4, e12.4)')'deltarli', mcmoves, j, newLyap(j)-oldlyap(j),exp(alpha*(newLyap(j)-oldLyap(j))), , rapporto(j)
 
 ranf=genrand()
  mcconf=min(1.0,rapporto(j))
!write(*,*) 'sel',rapporto(j),mcconf,ranf
     
     if (ranf.lt.mcconf) then
       !if (j.eq.0) 
       acc(j)=acc(j)+1
    ! write (*,*) 'cambio traiettoria'
      
       do i=0,ix-1
           P(j,i)%qx=Pbw(j,i)%qx
           P(j,i)%qy=Pbw(j,i)%qy
           P(j,i)%qz=Pbw(j,i)%qz

           P(j,i)%px=Pbw(j,i)%px
           P(j,i)%py=Pbw(j,i)%py
           P(j,i)%pz=Pbw(j,i)%pz
           P(j,i)%project=Pbw(j,i)%project
           P(j,i)%Lyap=Pbw(j,i)%Lyap
        enddo

        do i=ix,totiter-1
           P(j,i)%qx=Pfw(j,i)%qx
           P(j,i)%qy=Pfw(j,i)%qy
           P(j,i)%qz=Pfw(j,i)%qz

           P(j,i)%px=Pfw(j,i)%px
           P(j,i)%py=Pfw(j,i)%py
           P(j,i)%pz=Pfw(j,i)%pz
           P(j,i)%project=Pfw(j,i)%project
           P(j,i)%Lyap=Pfw(j,i)%Lyap
         enddo

      ener0(j)=ener0bw(j)
      xpq4ref(j)=q4bw(j)
      xpq6ref(j)=q6bw(j)
     
    
     endif

!write(*,*) 'qref', xpq4ref(j), xpq6ref(j)
q4bw(:)=0.0
q6bw(:)=0.0


enddo


!!!!!!!! mo' siccome siamo stronzi facciamo pure lo shifting: !!!!!!!!!!
!tiro numerello alla cazzo
lo=4

if (lo.EQ.4) then

do j=depart_boucle_nbclones,Nbclones

pix=nint(genrand()*totiter)
!write(*,*) 'pix', pix

do k=0,totiter-1
Pshift(j,pix+k+1)=P(j,k)
enddo

!write(*,*) 'eh!', Pshift(j,pix+totiter)%qx(1)
!write(*,*) 'eh!', P(j,totiter-1)%qx(1)
!stop
timebsh=dt*(pix+1)
 new_projection=.false.
 iterbw=pix+1
 !iterbw=nint(timebsh/dt)
         do while (iterbw.gt.1)
!do while (timebsh.gt.dt)
  !  iterbw=nint(timebsh/dt)
       call propagatbckw(Pshift(j,iterbw)%qx,Pshift(j,iterbw)%qy,Pshift(j,iterbw)%qz,Pshift(j,iterbw)%px,Pshift(j,iterbw)%py,&
                          Pshift(j,iterbw)%pz,rga,temperature,dt,N,dh(j,iterbw))
          call lanczos(N,Pshift(j,iterbw)%qx,Pshift(j,iterbw)%qy,Pshift(j,iterbw)%qz,new_projection,eigenvalue,Pshift(j,iterbw)%project)
          ! write(*,*) 'ctrlbw' , iterbw,Pshift(j,iterbw)%project(1)!eigenvalue, log(real(1.0+sqrt(abs(eigenvalue))))
          !if (iterbw.eq.pix) then 
                if (eigenvalue.lt.0.0) then 
                  Pshift(j,iterbw-1)%Lyap=log(real(1.0+dt*sqrt(abs(eigenvalue))))
                else 
                  Pshift(j,iterbw-1)%Lyap=0.0  !! autoval max del Lanczos
                endif
            !  write (*,*) 'lyaPshiftbw',iterbw-1, Pshift(j,iterbw-1)%Lyap
           !else
           !     if (eigenvalue.le.0.0) then 
           !       Pshift(j,iterbw)%Lyap = log(real(1.0+sqrt(abs(eigenvalue))))!+ Pshift(j,iterbw+1)%Lyap
           !      else 
           !       Pshift(j,iterbw)%Lyap=0.0!Pshift(j,iterbw+1)%Lyap  !! autoval max del Lanczos
           !     endif   ! autoval max del Lanczos
             
           ! endif
             !write (*,*) 'dhbw', mcmoves, iterbw, dh(j,iterbw)
      Pshift(j,iterbw-1)%qx=Pshift(j,iterbw)%qx
      Pshift(j,iterbw-1)%qy=Pshift(j,iterbw)%qy
      Pshift(j,iterbw-1)%qz=Pshift(j,iterbw)%qz

      Pshift(j,iterbw-1)%Px=Pshift(j,iterbw)%Px
      Pshift(j,iterbw-1)%Py=Pshift(j,iterbw)%Py
      Pshift(j,iterbw-1)%Pz=Pshift(j,iterbw)%Pz
      Pshift(j,iterbw-1)%Project=Pshift(j,iterbw)%Project
   iterbw=iterbw-1
         !timebsh=timebsh-dt  
   new_projection=.false.
enddo

do k=0,totiter-1
Pshift(j,pix+k+1)=P(j,k)
enddo


  Pshift(j,totiter+pix)%qx=P(j,totiter-1)%qx
  Pshift(j,totiter+pix)%qy=P(j,totiter-1)%qy
  Pshift(j,totiter+pix)%qz=P(j,totiter-1)%qz

  Pshift(j,totiter+pix)%Px=P(j,totiter-1)%Px
  Pshift(j,totiter+pix)%Py=P(j,totiter-1)%Py
  Pshift(j,totiter+pix)%Pz=P(j,totiter-1)%Pz

  Pshift(j,totiter+pix)%Lyap=P(j,totiter-1)%Lyap
  !Pshift(j,totiter+pix)%project=P(j,totiter-1)%project
timefsh=dt*(totiter+pix)
!write (*,*)'ma perché?' , P(j,totiter-1)%qx(1)
!write (*,*) timefsh, totiter+pix,real(2*totiter)*dt

 new_projection=.true.
 iterfw=nint(timefsh/dt)
    do while (iterfw.lt.2*totiter) 
 !do while ((timefsh.lt.(real(2*totiter)*dt)).and.(timefsh.gt.totaltime))
          ! iterfw=nint(timefsh/dt)
               ! write(*,*) 'fw' ,iterfw
                call propagatfwDET(Pshift(j,iterfw)%qx,Pshift(j,iterfw)%qy,Pshift(j,iterfw)%qz,Pshift(j,iterfw)%px,Pshift(j,iterfw)%py,&
                                   Pshift(j,iterfw)%pz,dt,N,dh(j,iterfw))
                call lanczos(N,Pshift(j,iterfw)%qx,Pshift(j,iterfw)%qy,Pshift(j,iterfw)%qz,new_projection,eigenvalue,Pshift(j,iterfw)%project)!!!!!
               !write(*,*) 'ctrlfw' ,iterfw+1, Pshift(j,iterfw)%project, eigenvalue!if (mod(iter,kappa)==0) then
                !if (iterfw.gt.totiter+pix) then 
                  if (eigenvalue.le.0.0) then 
                    Pshift(j,iterfw+1)%Lyap=log(real(1.0+dt*sqrt(abs(eigenvalue))))
                       !Pshift(j,iterfw)%Lyap=Pshift(j,iterfw-1)%Lyap+log(real(1.0+sqrt(abs(eigenvalue))))
                  else 
                    Pshift(j,iterfw+1)%Lyap=0.0!Pshift(j,iterfw-1)%Lyap !! autoval max del Lanczos
                  endif  
               !write (*,*) 'lyaPshiftfw',j, iterfw+1, Pshift(j,iterfw+1)%Lyap!/real(iterfw)
               ! endif
                 !write (*,*) 'dhfw', mcmoves, iterfw, dh(j,iterfw)
                !write (*,*) 'shfw', mcmoves, iterfw, j, dh(j,iterfw)
                !endif
                Pshift(j,iterfw+1)%qx=Pshift(j,iterfw)%qx
                Pshift(j,iterfw+1)%qy=Pshift(j,iterfw)%qy
                Pshift(j,iterfw+1)%qz=Pshift(j,iterfw)%qz

                Pshift(j,iterfw+1)%Px=Pshift(j,iterfw)%Px
                Pshift(j,iterfw+1)%Py=Pshift(j,iterfw)%Py
                Pshift(j,iterfw+1)%Pz=Pshift(j,iterfw)%Pz
                Pshift(j,iterfw+1)%Project=Pshift(j,iterfw)%Project
        iterfw=iterfw+1
         !timefsh=timefsh+dt
       new_projection=.false.
enddo

hamilt(j)=ener0(j)
!write(*,*) 'fh', hamilt(j)

S(j,:)=0

!!!!!!!!!!!! fare i pesi secondo le azioni 
do l=1,totiter
     do a=l-1,l+totiter-2
        S(j,l)=S(j,l)+Pshift(j,a)%Lyap
     enddo
enddo

S(j,:)=S(j,:)/(real(totiter))

!write(*,*)'s', s(j,:)
!!! ora L per la traiettoria data é fissato. Si devono calcolare i termini per stati di biais incrociati, ovvero  moltiplicare L (S) per i diversi valori di alpha

z=0.0
poids(:)=0.0

!Ora z rappresenta il exp(-"L_m) per la traiettoria da shiftare. Per il MBAR, si deve calcolare gli z diversi e scriverli su file:
! serve array del tipo Z(1,K) tale per cui

!!!!!!!!!!!!!! per ressort su Q4,Q6
         do k=1,2*totiter
               do i=1,N-1
                  do a=i+1,N
                    xij(i,a)=(Pshift(j,k)%qx(i-1)-Pshift(j,k)%qx(a-1))
                    yij(i,a)=(Pshift(j,k)%qy(i-1)-Pshift(j,k)%qy(a-1))
                    zij(i,a)=(Pshift(j,k)%qz(i-1)-Pshift(j,k)%qz(a-1))
                    dij(i,a)=sqrt(xij(i,a)**2+yij(i,a)**2+zij(i,a)**2)
                    xij(a,i)=-xij(i,a)
                    yij(a,i)=-yij(i,a)
                    zij(a,i)=-zij(i,a)
                    dij(a,i)= dij(i,a)
                  enddo
                enddo
            xpq4(j,k)=q4(N,xij,yij,zij,dij)
            xpq6(j,k)=q6(N,xij,yij,zij,dij)
        enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ESSAY PER MBAR: calcolo dei diversi pesi u_knl con waste recycling
if (waste_recycling) then
do l=0,nbclones
do i=1,totiter

if (defaut) then
!!!!!!!!!!!!!!!!!!!! ressort
!poids(i)=exp(alpha_bias(l)*S(j,i))*exp(-kapa*(10.0*(xpq4(j,i)-0.12d0)**2+(xpq6(j,i)-0.47d0)**2))*&
!         exp(-alpha_bias(l)*S(j,pix+1))*exp(kapa*(10.0*(xpq4(j,pix+1)-0.12d0)**2+(xpq6(j,pix+1)-0.47d0)**2))  !*exp(-hamilt(j)/temperature)

!!!!!!!!!!!!!!!!!!!!!!!!!! echelon
if((xpq4(j,i).gt.0.1).and.(xpq4(j,i).lt.0.13)) then
poids(i)=exp(alpha_bias(l)*S(j,i))
else
poids(i)=0.0
endif
!!!!!!!!!!!!! fine echelon

endif

if (fcc) then
if (alpha_bias(l).gt.0.0) then
   if (fcc) then
      poids(i)=exp(alpha_bias(l)*S(j,i))*exp(-kapa*(10.0*(xpq4(j,i)-0.18d0)**2+(xpq6(j,i)-0.57d0)**2))!*&
        ! exp(-alpha_bias(l)*S(j,pix+1))*exp(kapa*(10.0*(xpq4(j,pix+1)-0.18d0)**2+(xpq6(j,pix+1)-0.57d0)**2))!*exp(-hamilt(j)/temperature) !!! ressort sur Q4, Q6
   else if (defaut) then
     poids(i)=exp(alpha_bias(l)*S(j,i))*exp(-kapa*(10.0*(xpq4(j,i)-0.12d0)**2+(xpq6(j,i)-0.47d0)**2))!*&
         !exp(-alpha_bias(l)*S(j,pix+1))*exp(kapa*(10.0*(xpq4(j,pix+1)-0.12d0)**2+(xpq6(j,pix+1)-0.47d0)**2))  !*exp(-hamilt(j)/temperature)
   endif
else 
poids(i)=1.0
endif
endif


!poids(i)=exp(alpha_bias(l)*S(j,i))*exp(-kapa*(alpha_bias(l)/1.0d3)*(sum((Pshift(j,i)%qx(:)-xref(:))**2)+sum((Pshift(j,i)%qy(:)-yref(:))**2)+&
!         sum((Pshift(j,i)%qy(:)-yref(:))**2)))*exp(-hamilt(j)/temperature)
z=z+poids(i)
enddo
u_knl(j,l,mcmoves)=-log(z)
z=0.0
poids(:)=0.0
enddo
poids(:)=0.0
z=0.0
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
poids(:)=0.0
z=0.0




!!!!!!!!!!!! PER IL BIAS CORRENTE: mettere sempre !!!!!!!!!!!!!!!!!!!!
!!calcolo della cuvette (il secchiello) per tenere l'inizio della traiettoria nel bacino --> serve per h_A(x0) che ha forma gaussiana ATTENZIONE A RESSORT
do k=1,totiter

if (defaut) then
!!!!!!!!!!!!!!!!! using ressort 
!        poids(k)=exp(alpha_bias(j)*S(j,k))*exp(-kapa*(10.0*(xpq4(j,k)-0.12d0)**2+(xpq6(j,k)-0.47d0)**2))*&
!         exp(-alpha_bias(j)*S(j,pix+1))*exp(kapa*(10.0*(xpq4(j,pix+1)-0.12d0)**2+(xpq6(j,pix+1)-0.47d0)**2))!*exp(-hamilt(j)/temperature)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!using echelon
if((xpq4(j,k).gt.0.1).and.(xpq4(j,k).lt.0.13)) then
poids(k)=exp(alpha_bias(j)*S(j,k))*exp(-alpha_bias(j)*S(j,pix+1))
else
poids(k)=0.0
endif
!!!!!!!!!!! end echelon
endif


if (fcc) then
if (alpha_bias(j).ne.0.0) then
   if (fcc) then
      poids(k)=exp(alpha_bias(j)*S(j,k))*exp(-kapa*(10.0*(xpq4(j,k)-0.18d0)**2+(xpq6(j,k)-0.57d0)**2))*&
              exp(-alpha_bias(j)*S(j,pix+1))*exp(kapa*(10.0*(xpq4(j,pix+1)-0.18d0)**2+(xpq6(j,pix+1)-0.57d0)**2))!*exp(-hamilt(j)/temperature) !!! ressort sur Q4, Q6
   else if (defaut) then
       poids(k)=exp(alpha_bias(j)*S(j,k))*exp(-kapa*(10.0*(xpq4(j,k)-0.12d0)**2+(xpq6(j,k)-0.47d0)**2))*&
         exp(-alpha_bias(j)*S(j,pix+1))*exp(kapa*(10.0*(xpq4(j,pix+1)-0.12d0)**2+(xpq6(j,pix+1)-0.47d0)**2))!*exp(-hamilt(j)/temperature)
   endif
else 
poids(k)=1.0
endif
endif
!poids(k)=exp(alpha_bias(j)*S(j,k))*exp(-kapa*(alpha_bias(j)/1.0d3)*(9.0*(xpq4(j,k)-0.18d0)**2+(xpq6(j,k)-0.57d0)**2))*exp(-hamilt(j)/temperature)
!poids(k)=exp(alpha_bias(j)*S(j,k))*exp(-kapa*(alpha_bias(j)/real(1000))*(sum((Pshift(j,k)%qx(:)-xref(:))**2)+sum((Pshift(j,k)%qy(:)-yref(:))**2)+&
!                                       sum((Pshift(j,k)%qy(:)-yref(:))**2)))*exp(-hamilt(j)/real(temperature))
z=z+poids(k)
enddo

!write(*,*) 'pesi1',  z, poids(1:totiter)
poids(:)=poids(:)/z
!write(*,*) 'pesi2',  z, poids(1:totiter)
!write(*,*) 'kappa',exp(-kapa*(10.0*(xpq4(j,1:totiter)-0.12d0)**2+(xpq6(j,1:totiter)-0.47d0)**2))
!write(*,*)'alfa', exp(alpha_bias(j)*S(j,1:totiter))*exp(-alpha_bias(j)*S(j,pix+1))
!write(*,*) exp(-hamilt(j)/temperature)
!write(*,*) 'pesix',  z, poids(1:totiter)

!!!!!!!!!! (1)  calcolo della reattività mediata sulle traiettorie per MBAR con waste recycling: parto da defaut
if (waste_recycling) then

h_F(:)=0
h_Fd(:)=0
h_FI(:)=0
h_dF(:)=0
h_dI(:)=0

react_F(j)=0.0
react_Fd(j)=0.0
react_dI(j)=0.0
react_FI(j)=0.0
react_dF(j)=0.0

if (defaut) then !! parto da defaut 

!!!!!!!!!!!!!!!!!!!!!!!inizio old defaut

!!!!!!!!!!!! calcolo <h_d(0)*h_I(t)>
!!do k=1,totiter
!  if ((xpq4(j,k).lt.0.13).and.(xpq4(j,k).gt.0.1)) then
!     h_A(k)=1.0
!     if (xpq4(j,k+totiter).lt.0.04) then  !!! defaut -> ICO FUNNEL 
!      h_dI(k)=1.0
!     else
!      h_dI(k)=0.0
!     endif
!  else
 !  h_A(k)=0.0
!  endif
!enddo

!do k=1,totiter
!react_dI(j)=react_dI(j)+h_dI(k)*poids(k)
!enddo
!write (32,*) j, mcmoves, react_dI(j)

!!!!!!!!!!!! calcolo <h_d(0)*h_F(t)>
!do k=1,totiter
!  if ((xpq4(j,k).lt.0.13).and.(xpq4(j,k).gt.0.1)) then
!    if (xpq4(j,k+totiter).gt.0.13) then  !!! defaut -> fcc FUNNEL 
!      h_dF(k)=1.0
!    else
!      h_dF(k)=0.0
!    endif
!  endif
!enddo

!do k=1,totiter
!react_dF(j)=react_dF(j)+h_dF(k)*poids(k)
!enddo
!write (34,*) j, mcmoves, react_df(j)

!!!!!!!!!!!! calcolo <h_F(0)*h_d(t)>
!do k=1,totiter
!   if ((xpq4(j,k+totiter).lt.0.13).and.(xpq4(j,k+totiter).gt.0.1)) then  !!! defaut -> fcc FUNNEL 
!       if (xpq4(j,k).gt.0.13) then    !
!          h_Fd(k)=1.0
!       else
!          h_Fd(k)=0.0
!       endif
!   endif
!enddo

!do k=1,totiter
!react_Fd(j)=react_Fd(j)+h_Fd(k)*poids(k)
!enddo
!write (31,*) j, mcmoves, react_fd(j)

!!!!!!!!!!!!!!!! denominatore correl_function popolazioni <h_d(0)>
!do k=1,totiter
!h_A_average(j)=h_A_average(j)+h_A(k)*poids(k)
!enddo
!write (23,*) j, mcmoves,h_A_average(j)


!!!!!!!!!!!!!!!!!!! fine old defaut

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! inizio  prova nuovo deafult (per correl funct)

do tprim=1,100!,700,100
  do k=1,totiter
  ! passaggio F-d ed F-I
    if ((xpq4(j,k).lt.0.13).and.(xpq4(j,k).gt.0.1)) then! if (xpq4(j,k).gt.0.13) then !!FCC -> DEFAULT FCC
      h_A(k)=1.0
      if (xpq4(j,k+tprim*tprimo).gt.0.13) then !!!FCC -> DEFAULT FCC
          h_dF(k)=1.0
          !n_react=n_react+1.0
          !if (xpq4(j,k+totiter).gt.0.13) then
          !  h_dF(k)=1.0
          !endif
      else if ((xpq4(j,k+tprim*tprimo).lt.0.04).and.(xpq6(j,k+tprim*tprimo).gt.0.1)) then
        h_dI(k)=1.0
      endif
    endif
 ! passaggio d-I:
    ! if (((xpq4(j,k).gt.0.1).and.(xpq4(j,k).lt.0.13)).and.((xpq4(j,k+tprim*tprimo).lt.0.04).and.(xpq6(j,k+tprim*tprimo).gt.0.1))) then
    !  h_dI(k)=1.0
    !else  if (((xpq4(j,k+tprim*tprimo).gt.0.1).and.(xpq4(j,k+tprim*tprimo).lt.0.13)).and.((xpq4(j,k+2*tprim*tprimo).lt.0.04).and.(xpq6(j,k+2*tprim*tprimo).gt.0.1))) then
    !  h_dI(k+tprim*tprimo)=1.0
    !endif
   enddo ! k

do k=1,totiter
react_F(j)=react_F(j)+h_F(k)*poids(k)
react_Fd(j)=react_Fd(j)+h_Fd(k)*poids(k)
react_dI(j)=react_dI(j)+h_dI(k)*poids(k)
react_FI(j)=react_FI(j)+h_FI(k)*poids(k)
react_dF(j)=react_dF(j)+h_dF(k)*poids(k)
enddo

write (30,'(i,i,i,e12.4)') tprim, j, mcmoves, react_F(j)
write (31,'(i,i,i,e12.4)') tprim, j, mcmoves, react_Fd(j)
write (32,'(i,i,i,e12.4)') tprim, j, mcmoves, react_dI(j)
write (33,'(i,i,i,e12.4)') tprim, j, mcmoves, react_FI(j)
write (34,'(i,i,i,e12.4)') tprim, j, mcmoves, react_dF(j)

react_F(j)= 0.0
react_Fd(j)=0.0
react_dI(j)=0.0
react_FI(j)=0.0
react_dF(j)=0.0


enddo ! tprimo

!endif ! tprimo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! fine prova

endif !!!! fine cost partendo da defaut
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!(2) calcolo fonction de correlation pour constantes de reaction

if (fcc) then

h_F(:)=0
h_Fd(:)=0
h_FI(:)=0
h_dF(:)=0
h_dI(:)=0
correl_2(:,:)=0.0

if (tprimo.ne.0) then

do tprim=1,65!50!,700,100
 !write(*,*)'t1', tprim
 do k=1,totiter
  ! passaggio F-d ed F-I
    if (xpq4(j,k).gt.0.15) then !!FCC -> DEFAULT FCC
      h_F(k)=1.0
      if ((xpq4(j,k+tprim*tprimo).gt.0.1).and.(xpq4(j,k+tprim*tprimo).lt.0.13)) then !!!FCC -> DEFAULT FCC
          h_Fd(k)=1.0
          n_react=n_react+1.0
          if (xpq4(j,k+totiter).gt.0.13) then
            h_dF(k)=1.0
          endif
      else if ((xpq4(j,k+tprim*tprimo).lt.0.04).and.(xpq6(j,k+tprim*tprimo).gt.0.1)) then
        h_FI(k)=1.0
      endif
    endif
 ! passaggio d-I:
       tsecond=tprim+1
      do tsecond=tprim+1,nint(real(totiter)/tprimo)-tprim
 !write(*,*) tsecond, tsecond-tprim
   if (((xpq4(j,k+tprim*tprimo).gt.0.1).and.(xpq4(j,k+tprim*tprimo).lt.0.13)).and.((xpq4(j,k+tsecond*tprimo).lt.0.04).and.(xpq6(j,k+tsecond*tprimo).gt.0.1))) then
    !if (((xpq4(j,k+tprim*tprimo).gt.0.185).and.(xpq4(j,k+tprim*tprimo).lt.0.19)).and.((xpq4(j,k+tsecond*tprimo).lt.0.185).and.(xpq6(j,k+tsecond*tprimo).gt.0.1))) then
 ! h_dI(k)=1.0
     correl_2(tsecond-tprim,k)=correl_2(tsecond-tprim,k)+1.0/real(int(real(totiter)/real(tsecond*tprimo-tprim*tprimo)))   
   ! write(*,*)tsecond, tsecond-tprim,correl_2(tsecond-tprim,k),1.0/int(real(totiter)/real(tsecond*tprimo-tprim*tprimo)),real(tsecond*tprimo-tprim*tprimo)
!else  if (((xpq4(j,k+tprim*tprimo).gt.0.1).and.(xpq4(j,k+tprim*tprimo).lt.0.13)).and.((xpq4(j,k+2*tprim*tprimo).lt.0.04).and.(xpq6(j,k+2*tprim*tprimo).gt.0.1))) then
    !  h_dI(k+tprim*tprimo)=1.0
    endif
    enddo! tsecond
    !tsecond=tprim+1
   enddo ! k
!stop

do k=1,totiter
react_F(j)=react_F(j)+h_F(k)*poids(k)
react_Fd(j)=react_Fd(j)+h_Fd(k)*poids(k)
!react_dI(j)=react_dI(j)+h_dI(k)*poids(k)
react_FI(j)=react_FI(j)+h_FI(k)*poids(k)
react_dF(j)=react_dF(j)+h_dF(k)*poids(k)
enddo

write (30,'(i,i,i,e12.4)') tprim, j, mcmoves, react_F(j)
write (31,'(i,i,i,e12.4)') tprim, j, mcmoves, react_Fd(j)
!write (32,'(i,i,i,e12.4)') tprim, j, mcmoves, react_dI(j)
write (33,'(i,i,i,e12.4)') tprim, j, mcmoves, react_FI(j)
write (34,'(i,i,i,e12.4)') tprim, j, mcmoves, react_dF(j)

react_F(j)= 0.0
react_Fd(j)=0.0
react_dI(j)=0.0
react_FI(j)=0.0
react_dF(j)=0.0

enddo! tprimo


!do tprim=1,100
!do tsecond=tprim,nint(totiter/tprimo)-tprim
do i=1,nint(real(totiter)/real(tprimo))-1
do k=1,totiter
react_dI(j)=react_dI(j)+correl_2(i,k)*poids(k)
enddo
!write(*,'(i,i,i,e12.4)') i,j,mcmoves,react_dI(j)
write(32,'(i,i,i,e12.4)') i,j,mcmoves,react_dI(j)
react_dI(j)=0.0
enddo
!enddo
!if (((xpq4(j,tprim*tprimo).gt.0.1).and.(xpq4(j,tprim*tprimo).lt.0.13)).and.((xpq4(j,tprim*tprimo+tsecond*tprimo).lt.0.04).and.(xpq6(j,tprim*tprimo+tsecond*tprimo).gt.0.1))) !then
!      h_dI(tprim*tprimo)=1.0
!    endif
!   react_dI(j)=react_dI(j)+h_dI(tprim*tprimo)*poids(tprim*tprimo)
!enndo ! tsecondo

react_dI(j)=0.0
!enddo ! tprimo
!correl_2(i,k)
endif ! tprimo

endif ! fine partenza da fcc

n_react=n_react/(real(totiter)*5.0)
!!!!!!!!!!!!! fine correl funct

endif !!! fine waste recycling
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

x=poids(1)
 newtraj=1
do l=2,totiter
 if(poids(l).gt.x) then 
    newtraj=l 
  endif
enddo

!write(*,*) 'newtr', newtraj
oldlyap(j)=0.0

!!!!!!!! mo copio traiettoria shiftata
do k=0,totiter-1
P(j,k)=Pshift(j,k+newtraj)
oldLyap(j)=oldLyap(j)+P(j,k)%Lyap
enddo

oldLyap(j)=oldLyap(j)/real(totiter)

!!!!!!!!!!!!!!!!!! calcolo dei pesi per MBAR  post-shifting senza waste recycling CON RESSORT SU POSIZIONI X,Y,Z
!if(.not.waste_recycling) then
!do l=0,nbclones
!z=exp(alpha_bias(l)*S(j,newtraj))*exp(-kapa*(alpha_bias(l)/1.0d2)*(sum((Pshift(j,newtraj)%qx(:)-xref(:))**2)+sum((Pshift(j,newtraj)%qy(:)-yref(:))**2)+&
!         sum((Pshift(j,newtraj)%qy(:)-yref(:))**2)))*exp(-hamilt(j)/real(temperature))
!u_knl(j,l,mcmoves)=-log(z)
!enddo
!z=0.0
!endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
           do k=0,totiter-1
               do i=1,N-1
                  do a=i+1,N
                    xij(i,a)=(P(j,k)%qx(i-1)-P(j,k)%qx(a-1))
                    yij(i,a)=(P(j,k)%qy(i-1)-P(j,k)%qy(a-1))
                    zij(i,a)=(P(j,k)%qz(i-1)-P(j,k)%qz(a-1))
                    dij(i,a)=sqrt(xij(i,a)**2+yij(i,a)**2+zij(i,a)**2)
                    xij(a,i)=-xij(i,a)
                    yij(a,i)=-yij(i,a)
                    zij(a,i)=-zij(i,a)
                    dij(a,i)= dij(i,a)
                  enddo
                enddo

            xpq4(j,k)=q4(N,xij,yij,zij,dij)
            xpq6(j,k)=q6(N,xij,yij,zij,dij)
            call Force(P(j,k)%qx,P(j,k)%qy,P(j,k)%qz,fp,Utot)
            ener(j,k)=utot
            enddo !! su k


              do k=0,totiter-1   
              if (mod(k,100).eq.0) then
              write(*,'(a,i,i,i, e11.4, e12.4, e12.4)') 'traj',mcmoves, j, k, xpq4(j,k), xpq6(j,k), ener(j,k)!, dh(j,k) 
              endif
              enddo
              write(*,'(a,i,i,i, e11.4, e12.4, e12.4)') 'traj',mcmoves, j, totiter-1, xpq4(j,totiter-1), xpq6(j,totiter-1),ener(j,totiter-1)!, dh(j,totiter-1)
              write(*,*) 'traj'

!!!!!!!!!!!!!!!!!! calcolo dei pesi per MBAR  post-shifting senza waste recycling CON RESSORT SU Q4,Q6
if(.not.waste_recycling) then

do l=0,nbclones
!z=exp(alpha_bias(l)*S(j,newtraj))*exp(-kapa*(alpha_bias(l)/1.0d2)*(sum((P(j,0)%qx(:)-xref(:))**2)+sum((P(j,0)%qy(:)-yref(:))**2)+&
!         sum((P(j,0)%qy(:)-yref(:))**2)))*exp(-hamilt(j)/real(temperature))
if (alpha_bias(l).gt.0.0) then  
 if (fcc) then
    z=exp(alpha_bias(l)*S(j,newtraj))*exp(-kapa*(10.0*(xpq4(j,0)-0.18d0)**2+(xpq6(j,0)-0.57d0)**2))!*exp(-hamilt(j)/real(temperature))
  else if(defaut) then
    z=exp(alpha_bias(l)*S(j,newtraj))*exp(-kapa*(10.0*(xpq4(j,0)-0.12d0)**2+(xpq6(j,0)-0.47d0)**2))!*exp(-hamilt(j)/real(temperature))
  endif
else 
   z=1.0
endif
  u_knl(j,l,mcmoves)=-log(z)
enddo
z=0.0
       

!!!!!!!!!!!! reattività <h_A*h_B> senza waste recycling!!!!!!!!!!!!

h_F(:)=0
h_Fd(:)=0
h_FI(:)=0
h_dF(:)=0
h_dI(:)=0

!!!!!!!!!!!! reaction rate  partendo da defaut

if (defaut) then

!!!!!!!!!!!!!! calcolo <h_d*h_I>
 if ((xpq4(j,0).lt.0.13).and.(xpq4(j,0).gt.0.1)) then
    h_A(j)=1.0 
      if(xpq4(j,totiter-1).lt.0.04) then !! def->ICO
      h_dI(j)=1.0
      else
      h_dI(j)=0.0
      endif
  endif
write (32,*) j, mcmoves, h_dI(j)

!!!!!!!!!!!!! calcolo <h_d*h_F>
 if ((xpq4(j,0).lt.0.13).and.(xpq4(j,0).gt.0.1)) then
     if(xpq4(j,totiter-1).gt.0.13) then !! def->fcc
        h_dF(j)=1.0
      else
        h_dF(j)=0.0
      endif
  endif
write (34,*) j, mcmoves, h_dF(j)

!!!!!!!!!!! calcolo <h_f*h_d>
 
  if (xpq4(j,0).gt.0.13) then!! fcc->def
   if((xpq4(j,totiter-1).lt.0.13).and.(xpq4(j,totiter-1).gt.0.1)) then  
    h_Fd(j)=1.0
     else
      h_Fd(j)=0.0
     endif
  endif
write (31,*) j, mcmoves, h_Fd(j)

!!!!!!!! calcolo <h_d>
write (23,*) j, mcmoves,h_A(j)

endif !! defaut


!!!!!!!!!!fonction de correlation pour constantes de reaction: calcolo (2) se parto da fcc
if (fcc) then

h_F(:)=0
h_Fd(:)=0
h_FI(:)=0
h_dF(:)=0
h_dI(:)=0

if (tprimo.ne.0) then 

do tprim=1,50!,700,100
     if (xpq4(j,0).gt.0.13) then !!FCC -> DEFAULT FCC
        h_F(0)=1.0
        if ((xpq4(j,tprim*tprimo).gt.0.1).and.(xpq4(j,tprim*tprimo).lt.0.13)) then !!!FCC -> DEFAULT FCC
           h_Fd(0)=1.0
           n_react=n_react+1
           if (xpq4(j,totiter).gt.0.13) then
           h_dF(0)=1.0
           endif
        else if ((xpq4(j,totiter).lt.0.04).and.(xpq6(j,totiter).gt.0.1)) then
          h_FI(0)=1.0
        endif
      else if (((xpq4(j,totiter).lt.0.04).and.(xpq6(j,totiter).gt.0.1)).and.((xpq4(j,tprim*tprimo).gt.0.1).and.(xpq4(j,tprim*tprimo).lt.0.13))) then  !!! FCC FUNNEL -> ICO FUNNEL 
        h_dI(0)=1.0
      endif
   

write (30,'(i,i,i,e12.4)') tprim*tprimo, j, mcmoves, h_F(0)
write (31,'(i,i,i,e12.4)') tprim*tprimo, j, mcmoves, h_Fd(0)
write (32,'(i,i,i,e12.4)') tprim*tprimo, j, mcmoves, h_dI(0)
write (33,'(i,i,i,e12.4)') tprim*tprimo, j, mcmoves, h_FI(0)
write (34,'(i,i,i,e12.4)') tprim*tprimo, j, mcmoves, h_dF(0)

enddo

endif !!! fine tprimo

endif ! fine fcc

!!!!!!!!!!!!! fine correl funct
endif !! fine no waste recycling

enddo !! clones

endif !! shifting

reactivity(:)=0.0
h_A(:)=0.0
react(:)=0.0
h_A_average(:)=0.0
tau(:)=0.0
!!!!!!!!!!!!!! re-initialize correl function
h_F(:)=0
h_dF(:)=0
h_Fd(:)=0
h_FI(:)=0
h_dI(:)=0
react_F(:)=0
react_Fd(:)=0
react_dI(:)=0
react_FI(:)=0
react_dF(:)=0
!!!!!!!!!!!!!!!!!!!!!!!!!

   mcmoves=mcmoves+1
 enddo

write(24,'(i,i,i)') nbclones_mbar+1, nbclones_mbar+1, mcmoves
do k=depart_boucle_nbclones,nbclones
  do l=0,nbclones_mbar
     do i=0,mcmoves-1
      write(24,'(i,i,i,e15.4e3)') k,l,i, u_knl(k,l,i)
     enddo
   enddo
enddo

!write(*,*)'tau',  real(acc)/real(mcmoves)

       do k=0,totiter-1
            do l=depart_boucle_nbclones,NbClones  
                do i=1,N-1
                  do j=i+1,N
                    xij(i,j)=(P(l,k)%qx(i-1)-P(l,k)%qx(j-1))
                    yij(i,j)=(P(l,k)%qy(i-1)-P(l,k)%qy(j-1))
                    zij(i,j)=(P(l,k)%qz(i-1)-P(l,k)%qz(j-1))
                    dij(i,j)=sqrt(xij(i,j)**2+yij(i,j)**2+zij(i,j)**2)
                    xij(j,i)=-xij(i,j)
                    yij(j,i)=-yij(i,j)
                    zij(j,i)=-zij(i,j)
                    dij(j,i)= dij(i,j)
                  enddo
               enddo

            xpq4(L,k)=q4(N,xij,yij,zij,dij)
            xpq6(L,k)=q6(N,xij,yij,zij,dij)
            call Force(P(l,k)%qx,P(l,k)%qy,P(l,k)%qz,fp,Utot)
            ener(l,k)=utot
           enddo
        enddo
                  do j=depart_boucle_nbclones, Nbclones
                      do k=0,totiter-1
                      if (mod(k,100).eq.0) then
                       write(*,'(a,i, i,i, e11.4, e12.4, e12.4)') 'trajclonefin', mcmoves, j, k, xpq4(j,k),xpq6(j,k),ener(j,k)
                      endif
                      enddo
                   write(*,'(a,i,i,i, e11.4, e12.4, e12.4)') 'trajclonefin', mcmoves, j,totiter-1, xpq4(j,totiter-1),xpq6(j,totiter-1),ener(j,totiter-1)!,dh(j,totiter-1)
                enddo

do j=depart_boucle_nbclones, nbclones
write(*,*) 'tau',j,  alpha_bias(j) , tau(j)
write(*,*) 'acceptance ratio', j, alpha_bias(j), real(acc(j))/real(totalmcmoves)
write(*,*) '# traiettorie reattive', n_react/real(totalmcmoves)
enddo


stop

end subroutine LyapLanczos


subroutine lanczos(N,qx,qy,qz,new_projection,eigenvalue,projection)
  !use defs
  use lanczos_defs
  !use random_art
  !use art_in_ndm_module
  implicit none
  integer::N
  logical ::  new_projection
  logical :: lanczos_failed=.false.
  integer, dimension(5*N) :: iscratch
  real(8), dimension( 2 * 5*N -1 ) :: scratcha
  real(8), dimension(5*N) :: diag
  real(8), dimension(5*N-1) :: offdiag
  real(8), dimension(5*N, 5*N) :: vector
  real(8), dimension(3*N, 5*N), target :: lanc, proj
  real(8):: sum_forcenew, sum_force,total_energy

  ! Vectors used to build the matrix for Lanzcos algorithm 
  real(8), dimension(:), pointer :: z0, z1, z2
  ! Projection direction based on lanczos computations of lowest eigenvalues
  
  integer :: i,j,k, i_err, scratcha_size,ivec, nl_iter,nl_failed, it_art,evalf_number
  real(8) :: a1,a0,b2,b1,increment, eigenvalue
  real(8) :: excited_energy,c1,norm
  real(8) :: xsum, ysum, zsum, sum2, invsum
  real(8), dimension(3*N) :: pos, newpos,newforce,ref_force
  real(8), dimension(3*N) :: newforce1,newforce2,projection
  real (double), dimension (0:N-1) :: qx   ! vector of position
  real (double), dimension (0:N-1) :: qy   ! vector of position
  real (double), dimension (0:N-1) ::qz   ! vector of position




  !real(8) :: ran3
  do i =1,N
  pos(i)=qx(i-1)
  pos(i+N)=qy(i-1)
  pos(i+2*N)=qz(i-1)
  enddo

 lanczos_step = 1.d-6
  !boxl(:) = box(:) * scala
  increment = lanczos_step  ! Increment, convert in box units
  overlap=0.0
  evalf_number=0
  nl_iter=0
  nl_failed=0
!write(*,*)'boh1',evalf_number

  if(.not. new_projection ) then
    old_before_sc_projection = projection ! Vectorial operation
  end if

!write(*,*)'boh1',projection 
  ! We now take the current position as the reference point and will make 
  ! a displacement in a random direction or using the previous direction as
  ! the starting point.
  34 continue
  it_art=0
  call calcforce(N,pos,ref_force,total_energy)

  evalf_number = evalf_number + 1
  z0 => lanc(:,1)

  if(((.not. new_projection ).or. self_consistent).and.(.not. lanczos_failed)) then
    z0 = projection             ! Vectorial operation
    old_projection = projection ! Vectorial operation
    old_eigenvalue = eigenvalue 
  else
    do i=1, 3*N
      z0(i) = 0.5d0 - ran3()
    end do

    z1 => lanc(1,:)

    xsum = 0.0d0
    ysum = 0.0d0
    zsum = 0.0d0
    do i=1, N
      xsum = xsum + z0(i)
      ysum = ysum + z0(i+N)
      zsum = zsum + z0(i + 2*N)
    end do
    xsum = xsum / real(N)
    ysum = ysum / real(N)
    zsum = zsum / real(N)
    do i=1, N
      z0(i)  = z0(i)  - xsum
      z0(i+N )   = z0(i+N)   - ysum
      z0(i+2*N)  = z0(i+2*N) - zsum 
    end do

    if(first_time) then
      old_projection = z0   ! Vectorial operation
      first_time = .false.
    else
      old_projection = projection
    endif 
  endif
  ! We normalize the displacement to 1 total
  sum2 = 0.0d0
  do i=1, 3*N
    sum2 = sum2 + z0(i) * z0(i)
  end do
  invsum = 1.0/sqrt(sum2)
  z0 = z0 * invsum

    newpos = pos + z0 * increment
    call calcforce(N,newpos,newforce1,excited_energy)
    evalf_number = evalf_number + 1
    newpos = pos - z0 * increment
    call calcforce(N,newpos,newforce2,excited_energy)
    evalf_number = evalf_number + 1
    newforce = (newforce1 - newforce2)/2.d0  
   !write(*,*)'boh2',evalf_number
  ! We extract lanczos(1)

  ! We get a0
  a0 = 0.0d0
  do i=1, 3*N
    a0 = a0 + z0(i) * newforce(i)
  end do
  diag(1) = a0

  z1 => lanc(:,2)
  z1 = newforce - a0 * z0    ! Vectorial operation

  b1=0.0d0
  do i=1, 3*N
    b1 = b1 + z1(i) * z1(i)
  end do
  offdiag(1) = sqrt(b1)

  invsum = 1.0d0 / sqrt ( b1 )
  z1 = z1 * invsum           ! Vectorial operation
  
  ! We can now repeat this game for the next vectors
  do ivec = 2, 5*N-1
    z1 => lanc(:,ivec)
    newpos = pos + z1 * increment
    call calcforce(N,newpos,newforce1,excited_energy)
    evalf_number = evalf_number + 1
    newpos = pos - z1 * increment
    call calcforce(N,newpos,newforce2,excited_energy)
    evalf_number = evalf_number + 1
    newforce = (newforce1 - newforce2)/2.d0  

    a1 = 0.0d0
    do i=1, 3*N
      a1 = a1 + z1(i) * newforce(i)
    end do
    diag(ivec) = a1

    b1 = offdiag(ivec-1)
    z0 => lanc(:,ivec-1)
    z2 => lanc(:,ivec+1)
    z2 = newforce - a1*z1 -b1*z0

    b2=0.0d0
    do i=1, 3*N
      b2 = b2 + z2(i)*z2(i)
    end do

    offdiag(ivec) = sqrt(b2)
    
    invsum = 1.0/sqrt(b2)
    z2 = z2 * invsum
  end do

  ! We now consider the last line of our matrix
  ivec = 5*N
  z1 => lanc(:,5*N)
    newpos = pos + z1 * increment
    call calcforce(N,newpos,newforce1,excited_energy)
    evalf_number = evalf_number + 1
    newpos = pos - z1 * increment
    call calcforce(N,newpos,newforce2,excited_energy)
    evalf_number = evalf_number + 1
    newforce = (newforce1 - newforce2)/2.d0  
    sum_force = 0.0 
    sum_forcenew = 0.0

    do i = 1 , 3*N
       sum_force = sum_force + ref_force(i)
       sum_forcenew =  sum_forcenew + newforce(i) 
   end do
!   write(*,*) 'the sum of the forces before the move', sum_force

!   write(*,*) 'the sum of the forces After the move', sum_forcenew
! newforce = newforce - ref_force
  sum_forcenew = 0.0
  a1 = 0.0d0
  do i=1, 3*N
    a1 = a1 + z1(i) * newforce(i)
    sum_forcenew =  sum_forcenew + newforce(i)
  end do
!  write(*,*) 'the difference between the forces ' , sum_forcenew
!  write(*,*)
  diag(5*N) = a1

  ! We now have everything we need in order to diagonalise and find the
  ! eigenvectors.

  diag = -1.0d0 * diag
  offdiag = -1.0d0 * offdiag

  ! We now need the routines from Lapack. We define a few values
  i_err = 0

  ! We call the routine for diagonalizing a tridiagonal  matrix
  call dstev('V',5*N,diag,offdiag,vector,5*N,scratcha,i_err)


  ! We now reconstruct the eigenvectors in the real space
  ! Of course, we need only the first 5*N elements of vec

  projection = 0.0d0    ! Vectorial operation
  do k=1, 5*N
    z1 => lanc(:,k)
    a1 = vector(k,1)
    projection = projection + a1 * z1   ! Vectorial operation
  end do 
  c1=0.0d0
  do i=1, 3*N
     c1 = c1 + projection(i) * projection(i)
  end do

     norm = 1.0/sqrt(c1)
     projection = projection *norm 

  ! The following lines are probably not needed.
  !newpos = pos + projection * increment   ! Vectorial operation
  !call calcforce(N,type,newpos,boxl,newforce,excited_energy)
  !evalf_number = evalf_number + 1
  !newforce = newforce - ref_force

  eigenvalue=diag(1)/increment
  do i=1, 4
    eigenvals(i) = diag(i) / increment
  end do
 !write(*,*) 'lambda', eigenvalue

  a1=0.0d0
  b1=0.0d0
  do i=1, 3*N
    a1 = a1 + old_projection(i) * projection(i)
    b1 = b1 + projection(i) * projection(i)
    c1 = c1 + old_projection(i) * old_projection(i)

  end do
!ooo  
!ooo
!ooo!   The condition on the scalar product , we reject the point where we loose the eigen value and try to reduce the step size
!ooo   if(abs(a1)<=0.2) reject = .true. 
!ooo
 if(a1<0.0d0) then
    projection = -1.0d0 * projection
    overlap=-overlap
 end if    
!ooo
!ooo  
  call center(projection,3*N)
  
 !write(*,*)  eigenvalue!, 1.d-2, overlap 
 if ( dabs((old_eigenvalue-eigenvalue)) .gt. 5.d-5) then
   self_consistent=.true.
   lanczos_failed=.false.
   nl_iter=nl_iter+1
   nl_failed=nl_failed+1
   if (nl_failed.gt.30) then
      nl_failed=0
      lanczos_failed=.true.
      write(*,*) 'WARNING: LANCZOS FAILED ... we start with new random Krylov space'
   go to 35
   end if
   go to 34
 end if
35 continue 
 self_consistent=.false.
 
  a1=0.0d0
  do i=1, 3*N
    a1 = a1 + old_before_sc_projection(i) * projection(i)
  end do
  
  overlap=a1 
  lanczos_iter=nl_iter*5*N

!   The condition on the scalar product , we reject the point where we loose the eigen value and try to reduce the step size
   if(abs(a1)<=0.2) reject = .true. 

 overlap=a1

 if(a1<0.0d0) then
    projection = -1.0d0 * projection
    !overlap=-overlap
 end if    

 call center(projection,3*N)

!write(*,*) 'nliter', nl_iter,overlap
!write(*,*)'newproj',projection(1) 

end subroutine lanczos



subroutine center(vector,VECSIZE)
  integer, intent(IN) :: VECSIZE
  real(8), dimension(VECSIZE),intent(inout), target :: vector

  integer :: i, natoms
  real(8), dimension(:), pointer :: x, y, z     ! Pointers for coordinates
  real(8) :: xtotal, ytotal, ztotal

  natoms = VECSIZE / 3

  ! We first set-up pointers for the x, y, z components 
  x => vector(1:natoms)
  y => vector(natoms+1:2*natoms)
  z => vector(2*natoms+1:3*natoms)

  xtotal = 0.0d0
  ytotal = 0.0d0
  ztotal = 0.0d0

  do i = 1, natoms
    xtotal = xtotal + x(i)
    ytotal = ytotal + y(i)
    ztotal = ztotal + z(i)
  enddo 

  xtotal = xtotal / natoms
  ytotal = ytotal / natoms
  ztotal = ztotal / natoms

  do i = 1, natoms
    x(i) = x(i) - xtotal
    y(i) = y(i) - ytotal
    z(i) = z(i) - ztotal
  end do
end subroutine


Subroutine propagatfwDET(qx,qy, qz, px,py,pz,dt, N,dh)

 implicit none
 integer :: i
  integer::N
  real (double), dimension (0:N-1) :: qx   ! vector of position
  real (double),dimension (0:N-1) :: px   ! vector of impulsion
  real (double), dimension (0:N-1) :: qy   ! vector of position
  real (double),dimension (0:N-1) :: py   ! vector of impulsion
     real (double), dimension (0:N-1) ::qz   ! vector of position
  real (double),dimension (0:N-1) :: pz   ! vector of impulsion
  real (double) ::Utot,lapU,tconf,ff
  real (double) ::dt
  real (double) ::beq,fftot
  real (double) ::dh,ecinetique 
  real (double),dimension(0:N-1):: fxx
  real (double),dimension(0:N-1):: fzz
  real (double),dimension(0:N-1):: fyy
  real (double),dimension (3,N)::xp
  real (double),dimension (3,N)::x2p
  real (double),dimension (3,N)::vp
  real (double),dimension (3,N)::fp
  real (double),dimension (3,N)::f2p
  


 xp(1,1:N)=qx(0:N-1)
 xp(2,1:N)=qy(0:N-1)
 xp(3,1:N)=qz(0:N-1)

!write(*,*) 'eh?', xp(1,:)

 vp(1,1:N)=px
 vp(2,1:N)=py
 vp(3,1:N)=pz

call calfoljc2(xp,fp,utot)

ff=sum(fp(1,:)**2+fp(2,:)**2+fp(3,:)**2)

call deterministe2(xp,vp,fp,dt,beq,dh,ecinetique)

qx(0:N-1)=xp(1,1:N)
qy(0:N-1)=xp(2,1:N)
qz(0:N-1)=xp(3,1:N)

px(0:N-1)=vp(1,1:N)
py(0:N-1)=vp(2,1:N)
pz(0:N-1)=vp(3,1:N) 
!write(*,*) 'eh?', qx

end subroutine propagatfwDET


Subroutine propagatfw(qx,qy, qz, px,py,pz, uqx, uqy, uqz, upx,upy, upz,rga,temperature,dt, N,dh)

 implicit none
 integer :: i,alfa
  integer::N
  real (double), dimension (0:N-1) :: qx   ! vector of position
  real (double),dimension (0:N-1) :: px   ! vector of impulsion
  real (double),dimension (0:N-1) :: uqx  ! linear displacement
  real (double),dimension (0:N-1) :: upx  ! linear displacement
  real (double), dimension (0:N-1) :: qy   ! vector of position
  real (double),dimension (0:N-1) :: py   ! vector of impulsion
  real (double),dimension (0:N-1) :: uqy  ! linear displacement
  real (double),dimension (0:N-1) :: upy  ! linear displacement
  real (double), dimension (0:N-1) ::qz   ! vector of position
  real (double),dimension (0:N-1) :: pz   ! vector of impulsion
  real (double),dimension (0:N-1) :: uqz  ! linear displacement
  real (double),dimension (0:N-1) :: upz  ! linear displacement
  real (double) ::c1,a1,a2,b1,b2
  real (double) ::c2,Utot,lapU,tconf,ff
  real (double) ::dt
  real (double) ::gamma,rga,beq,fftot
  real (double) ::temperature,dh,ecinetique
  real (double) ::eta2 
  real (double),dimension (0:N-1) ::q1x 
  real (double),dimension (0:N-1) ::p1x
  real (double),dimension (0:N-1) ::q2x
  real (double),dimension (0:N-1) ::uq1x
  real (double),dimension (0:N-1) ::up1x
  real (double),dimension (0:N-1) ::uq2x
  real (double),dimension (0:N-1) ::q1y 
  real (double),dimension (0:N-1) ::p1y
  real (double),dimension (0:N-1) ::q2y
  real (double),dimension (0:N-1) ::uq1y
  real (double),dimension (0:N-1) ::up1y
  real (double),dimension (0:N-1) ::uq2y
  real (double),dimension (0:N-1) ::q1z 
  real (double),dimension (0:N-1) ::p1z
  real (double),dimension (0:N-1) ::q2z
  real (double),dimension (0:N-1) ::uq1z
  real (double),dimension (0:N-1) ::up1z
  real (double),dimension (0:N-1) ::uq2z
  real (double),dimension(0:N-1):: fxx
  real (double),dimension(0:N-1):: fzz
  real (double),dimension(0:N-1):: fyy
  real (double),dimension (3,N)::xp
  real (double),dimension (3,N)::x2p
  real (double),dimension (3,N)::vp
  real (double),dimension (3,N)::fp
  real (double),dimension (3,N)::f2p
  real(double), dimension (n)::vsecx,vsecy,vsecz,vsecx1,vsecy1, vsecz1

lapu=0

call Hessien (qx,qy,qz, uqx, uqy, uqz, Vsecx, vsecy,vsecz,lapU)

do i=0,N-1
  up1x(i)= upx(i)-dt/2d0*Vsecx(i+1)
  up1y(i)= upy(i)-dt/2d0*Vsecy(i+1)
  up1z(i)= upz(i)-dt/2d0*Vsecz(i+1)
enddo

do i=0,N-1
  uq1x(i)= uqx(i) + dt* up1x(i) 
  uq1y(i)= uqy(i) + dt* up1y(i)  
  uq1z(i)= uqz(i) + dt* up1z(i) 
enddo

 xp(1,1:N)=qx(0:N-1)
 xp(2,1:N)=qy(0:N-1)
 xp(3,1:N)=qz(0:N-1)

 vp(1,1:N)=px
 vp(2,1:N)=py
 vp(3,1:N)=pz

call calfoljc2(xp,fp,utot)

ff=sum(fp(1,:)**2+fp(2,:)**2+fp(3,:)**2)

call deterministe2(xp,vp,fp,dt,beq,dh,ecinetique)
 
qx(0:N-1)=xp(1,1:N)
qy(0:N-1)=xp(2,1:N)
qz(0:N-1)=xp(3,1:N)

px(0:N-1)=vp(1,1:N)
py(0:N-1)=vp(2,1:N)
pz(0:N-1)=vp(3,1:N) 

   alfa=-1

call Hessien(qx, qy,qz, uq1x, uq1y, uq1z, Vsecx1, vsecy1,vsecz1,lapU)
!
do i=0,N-1
   upx(i) = (up1x(i)-dt/2.0*Vsecx1(i+1))!*rga 
   upy(i) = (up1y(i)-dt/2.0*Vsecy1(i+1))!*rga
   upz(i) = (up1z(i)-dt/2.0*Vsecz1(i+1))!*rga
enddo

do i=0,N-1
   uqx(i) = uq1x(i)
   uqy(i) = uq1y(i) 
   uqz(i) = uq1z(i) 
enddo

end subroutine propagatfw


Subroutine propagatbckw(qx,qy, qz, px,py,pz, rga,temperature,dt, N,dh)

 implicit none
 integer :: i,alfa
  integer::N
  real (double), dimension (0:N-1) :: qx   ! vector of position
  real (double),dimension (0:N-1) :: px   ! vector of impulsion
  real (double), dimension (0:N-1) :: qy   ! vector of position
  real (double),dimension (0:N-1) :: py   ! vector of impulsion
  real (double), dimension (0:N-1) ::qz   ! vector of position
  real (double),dimension (0:N-1) :: pz   ! vector of impulsion
  ! linear displacement
   real (double) ::c2,Utot,lapU,tconf,ff
  real (double) ::dt
  real (double) ::gamma,rga,beq,fftot
  real (double) ::temperature,dh,ecinetique
  real (double) ::eta2 
  real (double),dimension (0:N-1) ::q1x 
  real (double),dimension (0:N-1) ::p1x
  real (double),dimension (0:N-1) ::q2x
  real (double),dimension (0:N-1) ::q1y 
  real (double),dimension (0:N-1) ::p1y
  real (double),dimension (0:N-1) ::q2y
  real (double),dimension (0:N-1) ::q1z 
  real (double),dimension (0:N-1) ::p1z
  real (double),dimension (0:N-1) ::q2z
  real (double),dimension(0:N-1):: fxx
  real (double),dimension(0:N-1):: fzz
  real (double),dimension(0:N-1):: fyy
  real (double),dimension (3,N)::xp
  real (double),dimension (3,N)::x2p
  real (double),dimension (3,N)::vp
  real (double),dimension (3,N)::fp
  real (double),dimension (3,N)::f2p


 xp(1,1:N)=qx(0:N-1)
 xp(2,1:N)=qy(0:N-1)
 xp(3,1:N)=qz(0:N-1)

 vp(1,1:N)=px
 vp(2,1:N)=py
 vp(3,1:N)=pz

call calfoljc2(xp,fp,utot)

ff=sum(fp(1,:)**2+fp(2,:)**2+fp(3,:)**2)

!call deterministe2(xp,vp,fp,dt,beq,dh,ecinetique)
call deterministebckw(xp,vp,fp,dt,rga, temperature,beq,dh,ecinetique)
 
qx(0:N-1)=xp(1,1:N)
qy(0:N-1)=xp(2,1:N)
qz(0:N-1)=xp(3,1:N)

px(0:N-1)=vp(1,1:N)
py(0:N-1)=vp(2,1:N)
pz(0:N-1)=vp(3,1:N) 

end subroutine propagatbckw



Subroutine reconstructLy (qx,qy, qz, q1x,q1y,q1z,uqx, uqy, uqz, upx,upy, upz,rga,temperature,dt, N)

 implicit none
 integer :: i,alfa
  integer::N
  real (double), dimension (0:N-1) :: qx   ! vector of position
  real (double),dimension (0:N-1) :: px   ! vector of impulsion
  real (double),dimension (0:N-1) :: uqx  ! linear displacement
  real (double),dimension (0:N-1) :: upx  ! linear displacement
  real (double), dimension (0:N-1) :: qy   ! vector of position
  real (double),dimension (0:N-1) :: py   ! vector of impulsion
  real (double),dimension (0:N-1) :: uqy  ! linear displacement
  real (double),dimension (0:N-1) :: upy  ! linear displacement
  real (double), dimension (0:N-1) ::qz   ! vector of position
  real (double),dimension (0:N-1) :: pz   ! vector of impulsion
  real (double),dimension (0:N-1) :: uqz  ! linear displacement
  real (double),dimension (0:N-1) :: upz  ! linear displacement
  real (double) ::c2,Utot,lapU,tconf,ff
  real (double) ::dt
  real (double) ::gamma,rga,beq,fftot
  real (double) ::temperature,dh,ecinetique
  real (double) ::eta2 
  real (double),dimension (0:N-1) ::q1x 
  real (double),dimension (0:N-1) ::p1x
  real (double),dimension (0:N-1) ::q2x
  real (double),dimension (0:N-1) ::uq1x
  real (double),dimension (0:N-1) ::up1x
  real (double),dimension (0:N-1) ::uq2x
  real (double),dimension (0:N-1) ::q1y 
  real (double),dimension (0:N-1) ::p1y
  real (double),dimension (0:N-1) ::q2y
  real (double),dimension (0:N-1) ::uq1y
  real (double),dimension (0:N-1) ::up1y
  real (double),dimension (0:N-1) ::uq2y
  real (double),dimension (0:N-1) ::q1z 
  real (double),dimension (0:N-1) ::p1z
  real (double),dimension (0:N-1) ::q2z
  real (double),dimension (0:N-1) ::uq1z
  real (double),dimension (0:N-1) ::up1z
  real (double),dimension (0:N-1) ::uq2z
  real (double),dimension(0:N-1):: fxx
  real (double),dimension(0:N-1):: fzz
  real (double),dimension(0:N-1):: fyy
  real (double),dimension (3,N)::xp
  real (double),dimension (3,N)::x2p
  real (double),dimension (3,N)::vp
  real (double),dimension (3,N)::fp
  real (double),dimension (3,N)::f2p
  real(double), dimension (n)::vsecx,vsecy,vsecz,vsecx1,vsecy1, vsecz1

lapU=0
call Hessien(qx, qy,qz, uqx, uqy, uqz, Vsecx, vsecy,vsecz,lapU)

do i=0,N-1
  up1x(i)= upx(i)-dt/2.d0* Vsecx(i+1) 
  up1y(i)= upy(i)-dt/2.d0* Vsecy(i+1) 
  up1z(i)= upz(i)-dt/2.d0* Vsecz(i+1) 
enddo

do i=0,N-1
  uq1x(i)= uqx(i) - dt* up1x(i) 
  uq1y(i)= uqy(i) - dt* up1y(i)  
  uq1z(i)= uqz(i) - dt* up1z(i) 
enddo

call Hessien(q1x, q1y,q1z, uq1x, uq1y, uq1z, Vsecx1, vsecy1,vsecz1,lapU)

do i=0,N-1
   upx(i) = up1x(i)-dt/2.d0*Vsecx1(i+1)!*rga 
   upy(i) = up1y(i)-dt/2.d0*Vsecy1(i+1)!*rga
   upz(i) = up1z(i)-dt/2.d0*Vsecz1(i+1)!*rga
enddo

do i=0,N-1
   uqx(i) = uq1x(i)
   uqy(i) = uq1y(i) 
   uqz(i) = uq1z(i) 
enddo

end subroutine reconstructLy

Subroutine propagateB(qx,qy, qz, px,py,pz, uqx, uqy, uqz, upx,upy, upz, c1, c2, dt, sqrtgTh, N)


 implicit none
 integer :: i
  integer::N
  real (double), dimension (0:N-1) :: qx   ! vector of position
  real (double),dimension (0:N-1) :: px   ! vector of impulsion
  real (double),dimension (0:N-1) :: uqx  ! linear displacement
  real (double),dimension (0:N-1) :: upx  ! linear displacement
  real (double), dimension (0:N-1) :: qy   ! vector of position
  real (double),dimension (0:N-1) :: py   ! vector of impulsion
  real (double),dimension (0:N-1) :: uqy  ! linear displacement
  real (double),dimension (0:N-1) :: upy  ! linear displacement
  real (double), dimension (0:N-1) ::qz   ! vector of position
  real (double),dimension (0:N-1) :: pz   ! vector of impulsion
  real (double),dimension (0:N-1) :: uqz  ! linear displacement
  real (double),dimension (0:N-1) :: upz  ! linear displacement
  real (double) ::c1,a1,a2,b1,b2
  real (double) ::c2,Utot,lapu
  real (double) ::dt
  real (double) ::sqrtgTh
  real (double) ::eta1
  real (double) ::eta2 
  real (double),dimension (0:N-1) ::q1x 
  real (double),dimension (0:N-1) ::p1x
  real (double),dimension (0:N-1) ::q2x
  real (double),dimension (0:N-1) ::uq1x
  real (double),dimension (0:N-1) ::up1x
  real (double),dimension (0:N-1) ::uq2x
  real (double),dimension (0:N-1) ::q1y 
  real (double),dimension (0:N-1) ::p1y
  real (double),dimension (0:N-1) ::q2y
  real (double),dimension (0:N-1) ::uq1y
  real (double),dimension (0:N-1) ::up1y
  real (double),dimension (0:N-1) ::uq2y
  real (double),dimension (0:N-1) ::q1z 
  real (double),dimension (0:N-1) ::p1z
  real (double),dimension (0:N-1) ::q2z
  real (double),dimension (0:N-1) ::uq1z
  real (double),dimension (0:N-1) ::up1z
  real (double),dimension (0:N-1) ::uq2z
  real (double),dimension(0:N-1):: fxx
  real (double),dimension(0:N-1):: fzz
  real (double),dimension(0:N-1):: fyy
  real (double),dimension (3,N)::xp
  real (double),dimension (3,N)::x2p
  real (double),dimension (3,N)::vp
  real (double),dimension (3,N)::fp
  real (double),dimension (3,N)::f2p
  real (double) ::sigma
real (double) ::l0

!call Force(qx,qy,qz,fp,Utot)

     ! xp(1,1:N)=qx
    ! xp(2,1:N)=qy
    !  xp(3,1:N)=qz

     ! vp(1,1:N)=px
     !vp(2,1:N)=py
     !vp(3,1:N)=pz

!write (*,*)'prima', xp

!call langevin_ljc(xp,vp,fp,dt)
  
!write (*,*)'dopo', xp

 !evalutation of the new speed and position*/

a1 = -1.0691860043307065
a2 = -0.1533230407019893
b1 = 0.3044913128854065
b2 = -1.0363164126095790

sigma=1
l0=0
!!!!!!!!METTRE BRUIT GAUSSIEN   
 call gauss(sigma,l0,eta1) !??? mettre le bon gasdev!!!
  
 call gauss(sigma,l0,eta2)!call random_number(eta2) !???ibidem
  
 do i=0,N-1 
   q1x(i)= qx(i) + dt / 4 * px(i)
   q1y(i)= qy(i) + dt / 4 * py(i)
   q1z(i)= qz(i) + dt / 4 * pz(i) 
!write (*,*) 'boh', q1z(i)
enddo

do i=0,N-1
  uq1x(i)= uqx(i) + dt / 4 * upx(i) 
  uq1y(i)= uqy(i) + dt / 4 * upy(i) 
  uq1z(i)= uqz(i) + dt / 4 * upz(i)
enddo

xp(1,1:N)=q1x(0:N-1)
xp(2,1:N)=q1y(0:N-1)
xp(3,1:N)=q1z(0:N-1)
!write (*,*) 'forze prima', fp(1,:)

call calfoljc2(xp,fp,utot)

!write (*,*) 'forze dopo', fp(1,:)
 
fxx(0:N-1)=fp(1,1:N)
fyy(0:N-1)=fp(2,1:N)
fzz(0:N-1)=fp(3,1:N)
  
do i=0,N-1
   p1x(i) = c2 * ( c1 * px(i)  + sqrtgTh * (a1 * eta1 + a2 * eta2)) + c2*dt/2 * fxx(i)
   p1y(i) = c2 * ( c1 * py(i)  + sqrtgTh * (a1 * eta1 + a2 * eta2)) + c2*dt/2 * fyy(i)  
   p1z(i) = c2 * ( c1 * pz(i)  + sqrtgTh * (a1 * eta1 + a2 * eta2)) + c2*dt/2 * fzz(i)
enddo
!

!do i=0,N-1
!qx(i)=xp(1,i+1)
!qy(i)=xp(2,i+1)
!qz(i)=xp(3,i+1)
!enddo
!write (*,*) 'boh3', qx

call Hessien(q1x, q1y,q1z, uq1x, uq1y, uq1z, Vsecx, vsecy,vsecz,lapu)

do i=0,N-1 
   up1x(i)= c2 *  c1 * upx(i) - c2*dt/2 * Vsecx(i+1)
   up1y(i)= c2 *  c1 * upy(i) - c2*dt/2 * Vsecy(i+1)
   up1z(i)= c2 *  c1 * upz(i) - c2*dt/2 * Vsecz(i+1)
!write (*,*) 'esselè upx', upx(i)
enddo

do i=0,N-1 
   q2x(i)= q1x(i) + dt / 2 * p1x(i)
   q2y(i)= q1y(i) + dt / 2 * p1y(i)
   q2z(i)= q1z(i) + dt / 2 * p1z(i)  
   uq2x(i) = uq1x(i) + dt/2 * up1x(i)    
   uq2y(i) = uq1y(i) + dt/2 * up1y(i)
   uq2z(i) = uq1z(i) + dt/2 * up1z(i)
!write (*,*) 'esselè2 uqx', q2x(i)
enddo

x2p(1,1:N)=q2x(0:N-1)
x2p(2,1:N)=q2y(0:N-1)
x2p(3,1:N)=q2z(0:N-1)

call calfoljc2(x2p,f2p,utot)

!write (*,*) 'forze2', f2p(1,:)

fxx(0:N-1)=f2p(1,1:N)
fyy(0:N-1)=f2p(2,1:N)
fzz(0:N-1)=f2p(3,1:N)

do i=0,N-1
  px(i) = c2 * ( c1 * p1x(i)  + sqrtgTh * ( b1 * eta1 + b2 * eta2 ) )+ c2*dt/2 * fxx(i)
  py(i) = c2 * ( c1 * p1y(i)  + sqrtgTh * ( b1 * eta1 + b2 * eta2 ) )+ c2*dt/2 * fyy(i)
  pz(i) = c2 * ( c1 * p1z(i)  + sqrtgTh * ( b1 * eta1 + b2 * eta2 ) )+ c2*dt/2 * fzz(i)
enddo

  call Hessien(q2x,q2y,q2z,uq2x, uq2y, uq2z,Vsecx, vsecy,vsecz,lapu)

do i=0,N-1
   upx(i) = c2 * c1 * up1x(i) -c2*dt/2 * Vsecx(i+1)  
   upy(i) = c2 * c1 * up1y(i) - c2*dt/2 * Vsecy(i+1)   
   upz(i) = c2 * c1 * up1z(i) - c2*dt/2 * Vsecz(i+1)   

   uqx(i) = uq2x(i) + dt/4* upx(i)
   uqy(i) = uq2y(i) + dt/4* upy(i)
   uqz(i) = uq2z(i) + dt/4* upz(i)
 ! Evaluate the tangent vector 
! write (*,*) 'uqx dopo', uqx(i) 
! write (*,*) 'upx dopo', upx(i) 
enddo

do i=0,N-1
   qx(i) = q2x(i) + dt/4 * px(i)
   qy(i) = q2y(i) + dt/4 * py(i)
   qz(i) = q2z(i) + dt/4 * pz(i)
enddo



end subroutine propagateB

Subroutine propagateC(qx,qy, qz, uqx, uqy, uqz, temperature,dt, N)

 implicit none
 integer :: i
  integer::N
  real (double), dimension (0:N-1) :: qx   ! vector of position
  real (double),dimension (0:N-1) :: uqx  ! linear displacement
  real (double), dimension (0:N-1) :: qy   ! vector of position
  real (double),dimension (0:N-1) :: uqy  ! linear displacement
  real (double), dimension (0:N-1) ::qz   ! vector of position
  real (double),dimension (0:N-1) :: uqz  ! linear displacement
  real (double) ::c1,a1,a2,b1,b2
  real (double) ::c2,Utot,lapu
  real (double) ::dt
  !real (double) ::gamma,rga
  real (double) ::temperature
 
  real (double),dimension (3,N)::fp
  real(double), dimension(6,N+1)::gau
  real(double), dimension(3,N) ::ss

!rga = exp(-gamma*dt/two)

call Hessien(qx, qy,qz, uqx, uqy, uqz, Vsecx, vsecy,vsecz,lapu)

!write (*,*)'dopo', up1x

do i=0,N-1
  uqx(i)= uqx(i) - (dt**2)/2* Vsecx(i+1)  
  uqy(i)= uqy(i) - (dt**2)/2* Vsecy(i+1)  
  uqz(i)= uqz(i) - (dt**2)/2* Vsecz(i+1)
enddo



!call calfoljc2(xp,fp,utot)

ss(:,:)=sqrt((dt**2)*temperature)

call genere_bruit2(ss,gau)

do i=0,N-1
 ! qx(i)= qx(i) + (dt**2)/2*fp(1,i+1) +gau(1,i+1)
 ! qy(i)= qy(i) + (dt**2)/2*fp(2,i+1) +gau(2,i+1)
!  qz(i)= qz(i) + (dt**2)/2*fp(3,i+1) +gau(3,i+1)
enddo

end subroutine propagateC


Subroutine gauss(ss,l0,l)
real (double) ::ss
real (double) ::l0
real (double) ::l
real (double) ::r,v1,v2
real (double)::x1,x2

r=2
do while (r.ge.1)
   call random_number(x1)
   call random_number(x2)
   v1=2.0*x1-1.0
   v2=2.0*x2-1.0
   r=v1**2+v2**2
enddo

l=v1*sqrt(-2.0*log(r)/(r))
l=l0+ss*l
return
end subroutine


Subroutine langevin_ljc2 (xp, vp, fp,dt,gamma,temperature) !!chiama genere_bruit2
    
    implicit none
    integer ic,i
    real(double)  :: xp(3,N)
    real(double)  :: xpp(3,N)
    real(double)  :: vp(3,N)
    real(double)  :: fp(3,N)
    real(double)  :: pp(3,N)
    real(double) :: xbar(3)
    real(double) :: gamma,temperature
    real(double) :: xalea,rga
    real(double) :: xprob, dt
    real(double) :: Ecinetique_d,potist
    integer   :: ic_local,iatom
    real(double), dimension(6,N+1)::gau
    real(double), dimension (3,N)::sig
    real(double) :: Ecin0
    real(double) :: Ecin1
    real(double) :: Ecin3
    real(double) :: Ecin4
    real(double) :: erreur
    
    
    im=N
    
    call calfoljq4(xp,fp)
     !write (*,*) 'gamma', gamma ,dt 
    deltawork = potistadd
    work=work+deltawork/temperature
    m_i(:,:) = one
    rga=0
    
    
    sig(:,:) = sqrt(m_i(:,:)*temperature*(one-exp(-gamma*dt)))

    call genere_bruit2(sig, gau)

     Ecin0 = zero
     Ecin1 = zero
     Ecin3 = zero
     Ecin4 = zero

    !write (*,*) 'rumore', gau(1:3,7)
   ! write (*,*) 'velocità', vp
!!!!langevin inerziale salto di rana
    
    pp(1:3,1:im) = vp(1:3,1:im)*m_i

     do ic =1,3
       Ecin0 = Ecin0   + DOT_PRODUCT(pp(ic,1:im),vp(ic,1:im))/two
    enddo

    rga = exp(-gamma*dt/two)
    
    xp_d(1:3,1:im)=xp(1:3,1:im)
    vp_d(1:3,1:im)=vp(1:3,1:im)
    Ecinetique_d=Ecin0
    !  write (*,*) 'rga', rga

    pp(1:3,1:im) = pp(1:3,1:im)*rga+gau(1:3,1:im)
    vp(1:3,1:im) = pp(1:3,1:im)/m_i(1:3,1:im)

     do ic =1,3
       Ecin1 = Ecin1   + DOT_PRODUCT(pp(ic,1:im),vp(ic,1:im))/two
     enddo
   
  ! write (*,*) 'velocità2', pp
   !write (*,*) 'forza', fp

    pp(1:3,1:im) =  pp(1:3,1:im) + (fp(1:3,1:im)+fadd(1:3,1:im))*dt/two 
    xp(1:3,1:im)= xp(1:3,1:im) + pp(1:3,1:im)*dt/m_i(1:3,1:im)

    do ic=1,3
       xbar(ic)    = sum(xp(ic,1:im))/real(im)-xbarini(ic) ! déplacement du barycentre
       xp(ic,1:im) = xp(ic,1:im) - xbar(ic)  ! on recentre tout le systeme 
    enddo
!write (*,*) 'pos dopo2', xp

    deltaU = potist + potistadd
    
    !call calfoljc2(xp,fp)
    call calfoljq4(xp,fp)
    deltaU = -deltaU + potist+ potistadd
    !write (*,*) 'forze', fp

    pp(1:3,1:im) = pp(1:3,1:im) + (fp(1:3,1:im)+fadd(1:3,1:im))*dt/two   
    vp(1:3,1:im) = pp(1:3,1:im) / m_i(1:3,1:im)
    
    do ic=1,3
       Ecin3 = Ecin3 + DOT_PRODUCT(pp(ic,1:im),vp(ic,1:im))/two
    enddo
   
    !write (*,*) 'velocità3', vp
    pp(1:3,1:im) = pp(1:3,1:im)*rga + gau(4:6,1:im)
    vp(1:3,1:im) = pp(1:3,1:im)/m_i(1:3,1:im)

    do ic=1,3
      Ecin4 = Ecin4 + DOT_PRODUCT(pp(ic,1:im),vp(ic,1:im))/two
    enddo

    erreur = ( deltaU + Ecin3 - Ecin1  )/(temperature)

    xprob = dexp(-erreur)
    call random_number(xalea)

!    write(*,*) 'erreur xprob xalea ', erreur ,xprob, xalea,bk*text_teledyn
!    write(*,*) ' deltaU ',deltaU,potistadd,potist
    if ((xalea.lt.xprob).and.(ldistance)) then  
!       print *,'acceptation'
       ncomptinter = ncomptinter + 1
       !betaq      = betaq     + erreur ! (Ecin4-Ecin3+Ecin1-Ecin0)/(bk*text_teledyn)
      Ecinetique = Ecin4
    else
       xp(1:3,1:im)=  xp_d(1:3,1:im)
       vp(1:3,1:im)=- vp_d(1:3,1:im)
       Ecinetique = Ecinetique_d
 
       !call calfoljc2(xp,fp)
!       print *,'refus'
       if (.not.ldistance) idistance = idistance +1
       ldistance = .true.
    endif

    return

End subroutine langevin_ljc2


Subroutine deterministe2(xd, vd, fd,dt,beq,dh,ecinetique)

    implicit none
    integer ic,i,itest
    
    real(double)  :: xd(3,N)
    real(double)  :: vd(3,N)
    real(double)  :: fd(3,N)
    real(double)  :: pp(3,N)

    real(double) :: xbar(3)

    real(double) Ecin0
    real(double) Ecin1
    real(double) Ecin3,ecinetique
    real(double) Ecin4,dh

    real(double) Ebar0
    real(double) Ebar1
    real(double) Ebar3
    real(double) Ebar4

    real(double), dimension (3,N)::sig
    real(double) ::temperature,dt,rga, tstep,utot,beq,gna,rum
    integer   :: ic_local,iatom
   
    im=N
    tstep=dt
     m_i(:,:) = one 

    Ecin0 = zero
    Ecin1 = zero
    Ecin3 = zero
    Ecin4 = zero
    ecinetique=0

    pp(1:3,1:im) = vd(1:3,1:im)!*m_i(1:3,1:im)
    Ecin0 =(DOT_PRODUCT(pp(1,1:im),vd(1,1:im))+DOT_PRODUCT(pp(2,1:im),vd(2,1:im))+DOT_PRODUCT(pp(3,1:im),vd(3,1:im)))/two
    Ecin1 = (DOT_PRODUCT(pp(1,1:im),vd(1,1:im))+DOT_PRODUCT(pp(2,1:im),vd(2,1:im))+DOT_PRODUCT(pp(3,1:im),vd(3,1:im)))/two
    pp(1:3,1:im) =  pp(1:3,1:im) + (fd(1:3,1:im))*tstep/two

    do ic=1,3
       xbar(ic)    = sum(xd(ic,1:im))/real(im) ! barycentre sur les particules
    enddo

    xd(1:3,1:im)=  xd(1:3,1:im) + pp(1:3,1:im)*tstep!/m_i(1:3,1:im)
!write(*,*) 'crist', xd(1,:)
    do ic=1,3
       xbar(ic)    = sum(xd(ic,1:im))/real(im)-xbarini(ic) ! déplacement du barycentre
       xd(ic,1:im) = xd(ic,1:im) - xbar(ic)  ! on recentre tout le systeme 
    enddo
!write(*,*) 'crist2',  xbar(1),xbarini(1) ,xd(1,:)
    call calfoljc2(xd,fd,utot)
!write(*,*) 'crist3', xd(1,:)
    pp(1:3,1:im) = pp(1:3,1:im) + (fd(1:3,1:im))*tstep/two
    vd(1:3,1:im) = pp(1:3,1:im)! / m_i(1:3,1:im)

       Ecin3 = (DOT_PRODUCT(pp(1,1:im),vd(1,1:im))+DOT_PRODUCT(pp(2,1:im),vd(2,1:im))+DOT_PRODUCT(pp(3,1:im),vd(3,1:im)))/two
 
    pp(1:3,1:im) = pp(1:3,1:im)!*rga 
    vd(1:3,1:im) = pp(1:3,1:im)!/m_i(1:3,1:im)

      Ecin4 = (DOT_PRODUCT(pp(1,1:im),vd(1,1:im))+DOT_PRODUCT(pp(2,1:im),vd(2,1:im))+DOT_PRODUCT(pp(3,1:im),vd(3,1:im)))/two


    beq=beq+(Ecin4-Ecin3+Ecin1-Ecin0)

    !dh=utot+173.252378415650
    !dh=ecin4/real(57)+utot+173.928
    dh=ecin4+utot+173.252378415650

     Ecinetique = Ecin4! +  ecinetique!          Ebar4*barbeta*(temperature)

    work = work + deltaU

    return

End subroutine deterministe2


Subroutine deterministebckw(xp, vp, fp,dt,rga,temperature,beq,dh,ecinetique)

    implicit none
    integer ic,i,itest
    
    real(double)  :: xp(3,imm)
    real(double)  :: vp(3,imm)
    real(double)  :: fp(3,imm)
    real(double)  :: pp(3,imm)

    real(double) :: xbar(3)

    real(double) Ecin0
    real(double) Ecin1
    real(double) Ecin3,ecinetique
    real(double) Ecin4,dh

    real(double) Ebar0
    real(double) Ebar1
    real(double) Ebar3
    real(double) Ebar4

    real(double), dimension (3,N)::sig
    real(double) ::temperature,dt,rga, tstep,utot,beq,gna,rum
    integer   :: ic_local,iatom
   
    im=N
    tstep=dt
     m_i(:,:) = one 

    Ecin0 = zero
    Ecin1 = zero
    Ecin3 = zero
    Ecin4 = zero
    ecinetique=0

    pp(1:3,1:im) = vp(1:3,1:im)!*m_i(1:3,1:im)

       Ecin0 =(DOT_PRODUCT(pp(1,1:im),vp(1,1:im))+DOT_PRODUCT(pp(2,1:im),vp(2,1:im))+DOT_PRODUCT(pp(3,1:im),vp(3,1:im)))/two

    pp(1:3,1:im) = -pp(1:3,1:im)!*rga 
    vp(1:3,1:im) = pp(1:3,1:im)!/m_i(1:3,1:im)

       Ecin1 = (DOT_PRODUCT(pp(1,1:im),vp(1,1:im))+DOT_PRODUCT(pp(2,1:im),vp(2,1:im))+DOT_PRODUCT(pp(3,1:im),vp(3,1:im)))/two
 
   pp(1:3,1:im) =  -pp(1:3,1:im)- (fp(1:3,1:im))*tstep/two
    
    do ic=1,3
       xbar(ic)    = sum(xp(ic,1:im))/real(im) ! barycentre sur les particules
    enddo

    xp(1:3,1:im)=  xp(1:3,1:im) - pp(1:3,1:im)*tstep!/m_i(1:3,1:im)

    do ic=1,3
       xbar(ic)    = sum(xp(ic,1:im))/real(im)-xbarini(ic) ! déplacement du barycentre
       xp(ic,1:im) = xp(ic,1:im) - xbar(ic)  ! on recentre tout le systeme 
    enddo
 
    call calfoljc2(xp,fp,utot)

    pp(1:3,1:im) = -pp(1:3,1:im) + (fp(1:3,1:im))*tstep/two
    vp(1:3,1:im) = pp(1:3,1:im)! / m_i(1:3,1:im)

       Ecin3 = (DOT_PRODUCT(pp(1,1:im),vp(1,1:im))+DOT_PRODUCT(pp(2,1:im),vp(2,1:im))+DOT_PRODUCT(pp(3,1:im),vp(3,1:im)))/two
 
    pp(1:3,1:im) = -pp(1:3,1:im)!*rga 
    vp(1:3,1:im) = pp(1:3,1:im)!/m_i(1:3,1:im)

  
      Ecin4 = (DOT_PRODUCT(pp(1,1:im),vp(1,1:im))+DOT_PRODUCT(pp(2,1:im),vp(2,1:im))+DOT_PRODUCT(pp(3,1:im),vp(3,1:im)))/two
   
    beq=beq+(Ecin4-Ecin3+Ecin1-Ecin0)

    dh=ecin4+utot+173.252378415650

    Ecinetique = Ecin4! +  ecinetique!         
    work = work + deltaU
 
    return

End subroutine deterministebckw

Subroutine langevinLJ_MD(xp, vp, fp,dt,rga,temperature,beq,dh,ecinetique)

    !   M o d u l e s
    !-----------------------------------------------
    implicit none
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer ic,i,itest,im
    
    real(double)  :: xp(3,N)
    real(double)  :: vp(3,N)
    real(double)  :: fp(3,N)
    real(double)  :: pp(3,N)

    real(double) :: xbar(3)
    real(double) :: xbarini(3)
    real(double) Ecin0
    real(double) Ecin1
    real(double) Ecin3,ecinetique
    real(double) Ecin4,dh

    real(double) Ebar0
    real(double) Ebar1
    real(double) Ebar3
    real(double) Ebar4

   real(double), dimension(6,N+1)::gau
    real(double), dimension (3,N)::sig, m_i
    real(double) ::temperature,dt,rga, tstep,utot,beq,gna,rum,box

   !-----------------------------------------------
    integer   :: ic_local,iatom
   !-----------------------------------------------

    im=N
    tstep=dt
     m_i(:,:) = one 
    box=10.0

    sig(:,:) = sqrt(temperature*(one-rga**2))

    call genere_bruit2(sig, gau)

    Ecin0 = zero
    Ecin1 = zero
    Ecin3 = zero
    Ecin4 = zero
    ecinetique=0


    !write (*,*) 'fp', fp(:,:)

    pp(1:3,1:im) = vp(1:3,1:im)!*m_i(1:3,1:im)

       Ecin0 =(DOT_PRODUCT(pp(1,1:im),vp(1,1:im))+DOT_PRODUCT(pp(2,1:im),vp(2,1:im))+DOT_PRODUCT(pp(3,1:im),vp(3,1:im)))/two

    pp(1:3,1:im) = pp(1:3,1:im)*rga + gau(1:3,1:im)
    vp(1:3,1:im) = pp(1:3,1:im)!/m_i(1:3,1:im)

       Ecin1 = (DOT_PRODUCT(pp(1,1:im),vp(1,1:im))+DOT_PRODUCT(pp(2,1:im),vp(2,1:im))+DOT_PRODUCT(pp(3,1:im),vp(3,1:im)))/two
  

   pp(1:3,1:im) =  pp(1:3,1:im) + (fp(1:3,1:im))*tstep/two
     !write (*,*) 'p', pp(1,45)
    do ic=1,3
       xbar(ic)    = sum(xp(ic,1:im))/real(im) ! barycentre sur les particules
    enddo

    xp(1:3,1:im)=  xp(1:3,1:im) + pp(1:3,1:im)*tstep!/m_i(1:3,1:im)

    do ic=1,3
       xbar(ic)    = sum(xp(ic,1:im))/real(im)-xbarini(ic) ! déplacement du barycentre
       xp(ic,1:im) = xp(ic,1:im) - xbar(ic)  ! on recentre tout le systeme 
    enddo

    call force_lj_md(xp,fp,box,utot)

    pp(1:3,1:im) = pp(1:3,1:im) + (fp(1:3,1:im))*tstep/two
    vp(1:3,1:im) = pp(1:3,1:im)! / m_i(1:3,1:im)

     Ecin3 = (DOT_PRODUCT(pp(1,1:im),vp(1,1:im))+DOT_PRODUCT(pp(2,1:im),vp(2,1:im))+DOT_PRODUCT(pp(3,1:im),vp(3,1:im)))/two

    pp(1:3,1:im) = pp(1:3,1:im)*rga + gau(4:6,1:im)
    vp(1:3,1:im) = pp(1:3,1:im)!/m_i(1:3,1:im)

      Ecin4 = (DOT_PRODUCT(pp(1,1:im),vp(1,1:im))+DOT_PRODUCT(pp(2,1:im),vp(2,1:im))+DOT_PRODUCT(pp(3,1:im),vp(3,1:im)))/(2.0)
   
   ! beq=beq+(Ecin4-Ecin3+Ecin1-Ecin0)

    !dh=ecin4+utot+173.252378415650

    Ecinetique = Ecin4! +  ecinetique!          Ebar4*barbeta*(temperature)

   ! write(*,*) 'ecin', ecinetique
   
    work = work + deltaU
    return

End subroutine langevinlj_md


Subroutine langevin2(xp, vp, fp,dt,rga,temperature,beq,dh,ecinetique)

    !   M o d u l e s
    !-----------------------------------------------
    implicit none
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer ic,i,itest
    
    real(double)  :: xp(3,imm)
    real(double)  :: vp(3,imm)
    real(double)  :: fp(3,imm)
    real(double)  :: pp(3,imm)

    real(double) :: xbar(3)
    real(double) :: xbarini(3)
    real(double) Ecin0
    real(double) Ecin1
    real(double) Ecin3,ecinetique
    real(double) Ecin4,dh

    real(double) Ebar0
    real(double) Ebar1
    real(double) Ebar3
    real(double) Ebar4

   real(double), dimension(6,N+1)::gau
    real(double), dimension (3,N)::sig, m_i
    real(double) ::temperature,dt,rga, tstep,utot,beq,gna,rum

   !-----------------------------------------------
    integer   :: ic_local,iatom
   !-----------------------------------------------

    im=N
    tstep=dt
     m_i(:,:) = one 

 !  write (*,*) 'ctrl', xp(1,:)
 !write (*,*) 'ctrlvel' ,rga,temperature

    
  !  rum=-log(rga**2)

    sig(:,:) = sqrt(temperature*(one-rga**2))

  !  sig(:,:) = sqrt(temperature*(rum))
    

!write(*,*) 'sig', rum, rga, sig(1,1) , sqrt(temperature*(one-rga**2))

!stop



    call genere_bruit2(sig, gau)

    Ecin0 = zero
    Ecin1 = zero
    Ecin3 = zero
    Ecin4 = zero
    ecinetique=0

  

   ! Ebar0 = zero
   ! Ebar1 = zero
   ! Ebar3 = zero
    !Ebar4 = zero

    pp(1:3,1:im) = vp(1:3,1:im)!*m_i(1:3,1:im)

!    write(*,*) ' potist  module =',potist,potistadd
!   write(*,*) 'fadd avant langevin', fadd(1:3,1),xpq4,xp(1,im1)
!   write(*,*) 'pp    pp(1:3,im1) ', pp(1:3,im1)
!    write(*,*) ' potistadd potistaddE', potistadd,potistaddE
!    deltaU =  potist/(bk*text_teledyn) + potistadd*barbeta +  potistaddE *barbeta 
   ! deltaU =  (potist + potistadd +  potistaddE)/(temperature)

   ! do ic =1,3
       Ecin0 =(DOT_PRODUCT(pp(1,1:im),vp(1,1:im))+DOT_PRODUCT(pp(2,1:im),vp(2,1:im))+DOT_PRODUCT(pp(3,1:im),vp(3,1:im)))/two
   ! enddo

  !write (*,*) 'ecin0', ecin0

    !Ebar0 =   barbeta*pp(1,im1)*vp(1,im1)/two + barbetaE*pp(2,im1)*vp(2,im1)/two 

    pp(1:3,1:im) = pp(1:3,1:im)*rga + gau(1:3,1:im)
    vp(1:3,1:im) = pp(1:3,1:im)!/m_i(1:3,1:im)

    !do ic =1,3
       Ecin1 = (DOT_PRODUCT(pp(1,1:im),vp(1,1:im))+DOT_PRODUCT(pp(2,1:im),vp(2,1:im))+DOT_PRODUCT(pp(3,1:im),vp(3,1:im)))/two
    !enddo

   ! write (*,*) 'ecin', ecin1

    !Ebar1 =   barbeta*pp(1,im1)*vp(1,im1)/two +  barbetaE*pp(2,im1)*vp(2,im1)/two 
!   call calfo  ! car déjà appelé

   pp(1:3,1:im) =  pp(1:3,1:im) + (fp(1:3,1:im))*tstep/two
    
    do ic=1,3
       xbar(ic)    = sum(xp(ic,1:im))/real(im) ! barycentre sur les particules
    enddo

    !write(*,*) ' baryc ', xbar

    xp(1:3,1:im)=  xp(1:3,1:im) + pp(1:3,1:im)*tstep!/m_i(1:3,1:im)

    !xp(1,im1)=  xp(1,im1) + pp(1,im1)*tstep/m_i(1,im1)
    !xp(2,im1)=  xp(2,im1) + pp(2,im1)*tstep/m_i(2,im1)

    do ic=1,3
       xbar(ic)    = sum(xp(ic,1:im))/real(im)-xbarini(ic) ! déplacement du barycentre
       xp(ic,1:im) = xp(ic,1:im) - xbar(ic)  ! on recentre tout le systeme 
    enddo
 
    !call force(xp(1,:),xp(2,:),xp(3,:),fp,utot)

    !write(*,*) ' fp en 1  ', fp(1,1)

    call calfoljc2(xp,fp,utot)

   !write(*,*) 'ctrl energ',utot

   !  write(11,*) ' fp en 1  ', SUM(fp(1,1:iM))

    !call calfoljq4(xp,fp) ! a besoin de fp pour calculer fadd  

    pp(1:3,1:im) = pp(1:3,1:im) + (fp(1:3,1:im))*tstep/two
    vp(1:3,1:im) = pp(1:3,1:im)! / m_i(1:3,1:im)

   ! write (*,*) 'velocità2', pp(1,1)
    !write (*,*) 'forza dopo', fp(1,1)

    !pp(1:2,im1) = pp(1:2,im1)
    !vp(1:2,im1) = pp(1:2,im1) / m_i(1:2,im1)

!    write(*,*) ' langevin  pp(1:2,im1) ',pp(1:2,im1)
!    write(*,*) ' langevin alphadd(1:2) ', fadd(1:2,im1)*alphadd(1:2)
!    write(*,*) '                       ',rga_i(1:3,im1)
!    stop

    !do ic=1,3
     Ecin3 = (DOT_PRODUCT(pp(1,1:im),vp(1,1:im))+DOT_PRODUCT(pp(2,1:im),vp(2,1:im))+DOT_PRODUCT(pp(3,1:im),vp(3,1:im)))/two
    !enddo

   ! write (*,*) 'ecin', ecin3
   ! Ebar3 =   barbeta*pp(1,im1)*vp(1,im1)/two + barbetaE*pp(2,im1)*vp(2,im1)/two 

    pp(1:3,1:im) = pp(1:3,1:im)*rga + gau(4:6,1:im)
    vp(1:3,1:im) = pp(1:3,1:im)!/m_i(1:3,1:im)

    !do ic=1,3
      Ecin4 = (DOT_PRODUCT(pp(1,1:im),vp(1,1:im))+DOT_PRODUCT(pp(2,1:im),vp(2,1:im))+DOT_PRODUCT(pp(3,1:im),vp(3,1:im)))/(2.0)
    !enddo

    !write (*,*) 'ecin', ecin4

    !Ebar4 =  barbeta*pp(1,im1)*vp(1,im1)/two + barbetaE*pp(2,im1)*vp(2,im1)/two 

   ! beq=beq+(Ecin4-Ecin3+Ecin1-Ecin0) !/(temperature)

    beq=beq+(Ecin4-Ecin3+Ecin1-Ecin0)

    dh=ecin4+utot+173.252378415650

    !write(*,*) 'betaq',Ecin4-Ecin3+Ecin1-Ecin0
   ! barbetaq   = barbetaq  + (Ebar4-Ebar3+Ebar1-Ebar0)
    Ecinetique = Ecin4! +  ecinetique!          Ebar4*barbeta*(temperature)

!    deltaU =  (Ecin3-Ecin1)/(bk*text_teledyn) + (Ebar3-Ebar1)*barbeta +   (potist )/(bk*text_teledyn) + potistadd*barbeta +  potistaddE *barbeta -deltaU    
  !  deltaU =  (Ecin3-Ecin1)/(bk*text_teledyn) +   (potist+potistadd+potistaddE )/(temperature) -deltaU

    work = work + deltaU
   
  !  write (*,*) 'dissipation', ecin4+utot+173.252378415650 , beq
   !write (*,*) xp(1,4)

    return

End subroutine langevin2

Subroutine langevinbckw(xp, vp, fp,dt,rga,temperature,beq,dh,ecinetique)

    !   M o d u l e s
    !-----------------------------------------------
    implicit none
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer ic,i,itest
    
    real(double)  :: xp(3,imm)
    real(double)  :: vp(3,imm)
    real(double)  :: fp(3,imm)
    real(double)  :: pp(3,imm)

    real(double) :: xbar(3)

    real(double) Ecin0
    real(double) Ecin1
    real(double) Ecin3,ecinetique
    real(double) Ecin4,dh

    real(double) Ebar0
    real(double) Ebar1
    real(double) Ebar3
    real(double) Ebar4

   real(double), dimension(6,N+1)::gau
    real(double), dimension (3,N)::sig
    real(double) ::temperature,dt,rga, tstep,utot,beq,gna,rum

   !-----------------------------------------------
    integer   :: ic_local,iatom
   !-----------------------------------------------

    im=N
    tstep=dt
     m_i(:,:) = one 

    sig(:,:) = sqrt(temperature*(one-rga**2))

    call genere_bruit2(sig, gau)

    Ecin0 = zero
    Ecin1 = zero
    Ecin3 = zero
    Ecin4 = zero
    ecinetique=0

    pp(1:3,1:im) = vp(1:3,1:im)!*m_i(1:3,1:im)

       Ecin0 =(DOT_PRODUCT(pp(1,1:im),vp(1,1:im))+DOT_PRODUCT(pp(2,1:im),vp(2,1:im))+DOT_PRODUCT(pp(3,1:im),vp(3,1:im)))/two

    pp(1:3,1:im) = -pp(1:3,1:im)*rga + gau(1:3,1:im)
    vp(1:3,1:im) = pp(1:3,1:im)!/m_i(1:3,1:im)

       Ecin1 = (DOT_PRODUCT(pp(1,1:im),vp(1,1:im))+DOT_PRODUCT(pp(2,1:im),vp(2,1:im))+DOT_PRODUCT(pp(3,1:im),vp(3,1:im)))/two
 
   pp(1:3,1:im) =  -pp(1:3,1:im)- (fp(1:3,1:im))*tstep/two
    
    do ic=1,3
       xbar(ic)    = sum(xp(ic,1:im))/real(im) ! barycentre sur les particules
    enddo

    xp(1:3,1:im)=  xp(1:3,1:im) - pp(1:3,1:im)*tstep!/m_i(1:3,1:im)

    do ic=1,3
       xbar(ic)    = sum(xp(ic,1:im))/real(im)-xbarini(ic) ! déplacement du barycentre
       xp(ic,1:im) = xp(ic,1:im) - xbar(ic)  ! on recentre tout le systeme 
    enddo
 
    call calfoljc2(xp,fp,utot)

    pp(1:3,1:im) = -pp(1:3,1:im) + (fp(1:3,1:im))*tstep/two
    vp(1:3,1:im) = pp(1:3,1:im)! / m_i(1:3,1:im)

       Ecin3 = (DOT_PRODUCT(pp(1,1:im),vp(1,1:im))+DOT_PRODUCT(pp(2,1:im),vp(2,1:im))+DOT_PRODUCT(pp(3,1:im),vp(3,1:im)))/two
 
    pp(1:3,1:im) = -pp(1:3,1:im)*rga - gau(4:6,1:im)
    vp(1:3,1:im) = pp(1:3,1:im)!/m_i(1:3,1:im)

  
      Ecin4 = (DOT_PRODUCT(pp(1,1:im),vp(1,1:im))+DOT_PRODUCT(pp(2,1:im),vp(2,1:im))+DOT_PRODUCT(pp(3,1:im),vp(3,1:im)))/two
   
    beq=beq+(Ecin4-Ecin3+Ecin1-Ecin0)

    dh=ecin4+utot+173.252378415650

    Ecinetique = Ecin4! +  ecinetique!         
    work = work + deltaU
 
    return

End subroutine langevinbckw


subroutine genere_bruit2 (sig,gau)
   
   implicit none

   real(double), dimension(6,N+1)::gau
   real(double),dimension(3,N)::sig
   real(double) :: u1,u2,b1
   real(double):: deriv(6)
   integer   :: ic,iatom,i
 
   im=N
   gau(:,:)=0   

!write(*,*) 'sig' ,sig(:,1)

   do i=1,im+1
      do ic=1,6
         call random_number(u1)
         call random_number(u2)
         b1=sqrt(-2.*log(u1))*cos(2.*pi*u2)
         gau(ic,i) = b1
      enddo
   enddo

 
!write(*,*) 'gau', gau(:,1)


   gau(1:3,1:im) = gau(1:3,1:im)*sig(1:3,1:im)
   gau(4:6,1:im) = gau(4:6,1:im)*sig(1:3,1:im)
!write(*,*) 'gau2', gau(:,1)

   do ic=1,6
      deriv(ic)    = sum(gau(ic,1:im))/real(im)
      gau(ic,1:im) = gau(ic,1:im)-deriv(ic)
   enddo

  end subroutine genere_bruit2


subroutine genere_bruit3 (sig,gau)
   
   implicit none

   real(double), dimension(6,N+1)::gau
   real(double),dimension(3,N)::sig
   real(double) :: u1,u2,b1,r
   real(double):: deriv(6)
   integer   :: ic,iatom,i
 
   im=N
   gau(:,:)=0   

write(*,*) 'sig' ,sig(:,1)
   
   do i=1,im+1
      do ic=1,6
         r=2
         do while(r.ge.1)
         call random_number(u1)
         u1=2*u1-1
         call random_number(u2)
         u2=2*u2-1
         r=u1*u1+u2*u2
         enddo
         b1=u1*sqrt(-2.*log(r)/r)
         gau(ic,i) = b1
      enddo
   enddo

write(*,*) 'gau', gau(:,1)

   gau(1:3,1:im) = gau(1:3,1:im)*sig(1:3,1:im)
   gau(4:6,1:im) = gau(4:6,1:im)*sig(1:3,1:im)
write(*,*) 'gau2', gau(:,1)

   do ic=1,6
      deriv(ic)    = sum(gau(ic,1:im))/real(im)
      gau(ic,1:im) = gau(ic,1:im)-deriv(ic)
   enddo

  end subroutine genere_bruit3


subroutine calfoljc2(xp,fp,potist) !calcul des positions reciproques,des forces et appelle q4 et q6 à partir des positions 
  !  USE T_kind_param_m, ONLY:  double
  !  use gen_com_m
  implicit none
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
 
  real(double)  :: xp(3,N)
  real(double)  :: fp(3,N)
  real(double)  :: c(3)
  real(double)  :: r2,r2sq,potemp,phu,phu2
  real(double)  :: frap(1:3)
  real(double)  :: fraptot(1:3),xbar(1:3)
  real(double)  :: f(1:3)
  real(double)  :: r6,r8,r12,r14
  real(double)  :: potistlj,potist
  real(double)  :: ra(1:3)
  real(double)  :: q4,q6
  integer  :: i,j,icri,ic
  real(double) :: uu(3,N)
!potist=0
potistlj=0 
fp(:,:)=0
  fraptot(:)=zero
  icri = 0
uu(:,:)=0
  
    sig_lj  = one
  rsig_lj = 2.25*sig_lj 
!rsig_lj = 1.75*sig_lj
!write(*,*) 'prova', xp(1,:)

!write(*,*) potist


  do ic=1,3
     xbar(ic)=sum(xp(ic,1:N))/real(N)
  enddo
  
  do  i=1,N
     !     anti-évaporation 
     c(1:3)=xp(1:3,i)-xbar(1:3)
     r2=sum(c(1:3)**2)
     r2sq=sqrt(r2)
          
     if (r2sq.gt.rsig_lj) then 
        potemp = (r2sq-rsig_lj)**3
        potistlj = potistlj + potemp
        phu    = -3*((r2sq-rsig_lj)**2)/r2sq
        frap(1:3)    = phu*c(1:3)
        fraptot(1:3) = phu*c(1:3) + fraptot(1:3)
        f(1:3)=fp(1:3,i)+frap(1:3)
                
        !icri=2
     else
        f(1:3)=fp(1:3,i)
     endif

         ! ouverture boucle sur j
        do j=i+1,N
         c(1:3)=xp(1:3,i)-xp(1:3,j)
         r2=sum(c(1:3)**2)
         r2sq=sqrt(r2)
         r6=(sig_lj/r2sq)**6
         r8=r6/r2
         r12=(r6)**2
         r14=r12/r2
         potemp = eps_4lj*(r12-r6) 
         potistlj = potistlj + potemp
         phu = eps_4lj*(12.0*r14 - 6.0*r8) 
         !phu2= eps_4lj*(156.0*r14-42.0*r8)
         !calcul de la force vectorielle
         ra(1:3)=phu*c(1:3)
         f(1:3) = f(1:3) + ra(1:3)
         fp(1:3,j)=fp(1:3,j)-ra(1:3)
         pist=pist+phu*r2 ! deuxieme partie du viriel

        ! Uu(1:3,i)=uu(1:3,i)-phu2*c(1:3)/r2-phu*(r2**2-c(1:3)**2)/(r2**3)

        ! tabulation pour Q4 et grad_Q4  
       !  xij(i,j)=-c(1)
       !  yij(i,j)=-c(2)
       !  zij(i,j)=-c(3)
       !  dij(i,j)=sqrt(xij(i,j)**2+yij(i,j)**2+zij(i,j)**2)
       !  xij(j,i)=-xij(i,j)
       !  yij(j,i)=-yij(i,j)
       !  zij(j,i)=-zij(i,j)
       !  dij(j,i)= dij(i,j)

      enddo
        fp(1:3,i)=f(1:3)
   enddo
!write(*,*) 'prova2', xp(1,:)
      potist=potistlj

    !write (*,*)'uu', uu(:,:)

   !xpq4 = q4(im,xij,yij,zij,dij)
   !call fq4(im,xij,yij,zij,dij,fq4x,fq4y,fq4z)

   !if (lq6) then
    !  xpq6 = q6(im,xij,yij,zij,dij)
!   call fq6(im,xij,yij,zij,dij,fq6x,fq6y,fq6z)
   !endif

   do ic=1,3
      fp(ic,1:N) = fp(ic,1:N)-fraptot(ic)/real(N)
   enddo
!write (*,*)'uu', fp(:,:)

end subroutine calfoljc2


subroutine calcforce(N, pos,nforce,potist) !calcul des forces POUR LANCZOS
  !  USE T_kind_param_m, ONLY:  double
  !  use gen_com_m
  implicit none
  integer::N
  real(double)  :: xp(3,N)
  real(double)  :: pos(3*N)
  real(double)  :: fp(3,N)
  real(double)  :: nforce(3*N)
  real(double)  :: c(3)
  real(double)  :: r2,r2sq,potemp,phu,phu2
  real(double)  :: frap(1:3)
  real(double)  :: fraptot(1:3),xbar(1:3)
  real(double)  :: f(1:3)
  real(double)  :: r6,r8,r12,r14
  real(double)  :: potistlj,potist
  real(double)  :: ra(1:3)
  !real(double)  :: q4,q6
  integer  :: i,j,icri,ic
  real(double) :: uu(3,N)
   
do i=1,N
  xp(1,i)=pos(i)
  xp(2,i)=pos(i+N)
  xp(3,i)=pos(i+2*N)
enddo

do i=1,N
  fp(1,i)=nforce(i)
  fp(2,i)=nforce(i+N)
  fp(3,i)=nforce(i+2*N)
enddo

 
   potistlj=0 
   fp(:,:)=0
   fraptot(:)=zero
   icri = 0
   uu(:,:)=0
  
   sig_lj  = one
   rsig_lj = 2.25*sig_lj 

  do ic=1,3
     xbar(ic)=sum(xp(ic,1:N))/real(N)
  enddo
  
  do  i=1,N
     !     anti-évaporation 
     c(1:3)=xp(1:3,i)-xbar(1:3)
     r2=sum(c(1:3)**2)
     r2sq=sqrt(r2)
     if (r2sq.gt.rsig_lj) then 
        potemp = (r2sq-rsig_lj)**3
        potistlj = potistlj + potemp
        phu    = -3*((r2sq-rsig_lj)**2)/r2sq
        frap(1:3)    = phu*c(1:3)
        fraptot(1:3) = phu*c(1:3) + fraptot(1:3)
        f(1:3)=fp(1:3,i)+frap(1:3)
     else
        f(1:3)=fp(1:3,i)
     endif
        do j=i+1,N
         c(1:3)=xp(1:3,i)-xp(1:3,j)
         r2=sum(c(1:3)**2)
         r2sq=sqrt(r2)
         r6=(sig_lj/r2sq)**6
         r8=r6/r2
         r12=(r6)**2
         r14=r12/r2
         potemp = eps_4lj*(r12-r6) 
         potistlj = potistlj + potemp
         phu = eps_4lj*(12.0*r14 - 6.0*r8) 
         !phu2= eps_4lj*(156.0*r14-42.0*r8)
         !calcul de la force vectorielle
         ra(1:3)=phu*c(1:3)
         f(1:3) = f(1:3) + ra(1:3)
         fp(1:3,j)=fp(1:3,j)-ra(1:3)
         pist=pist+phu*r2 ! deuxieme partie du viriel
      enddo
        fp(1:3,i)=f(1:3)
   enddo
      potist=potistlj

   do ic=1,3
      fp(ic,1:N) = fp(ic,1:N)-fraptot(ic)/real(N)
   enddo

do i=1,N
  pos(i)=xp(1,i)
  pos(i+N)=xp(2,i)
  pos(i+2*N)=xp(3,i)
enddo

do i=1,N
  nforce(i)=fp(1,i)
  nforce(i+N)=fp(2,i)
  nforce(i+2*N)=fp(3,i)
enddo


end subroutine calcforce

recursive Subroutine Combinaison(C,N,Cb,Nb)

 integer ::C
 integer ::N
integer, dimension (0:C-1)::Cb 
integer, dimension (0:N-1) ::Nb
 integer ::ran
  ! we pick up C guys among N.
  ! Nb is an array of numbers.
 ! Cb is the array of the guys which have been chosen
!------------- 
  integer :: i

do i=0,N-1 
  Nb(i)=i
 enddo
!------------

  if(C.eq.0) then
    return;
  elseif ((C.lt.0).and.(C.gt.N)) then
         write (6,*) 'Bug Combinaison\n'
	  stop
	
      else
	  ran = floor(genrand()*N)
	  Cb(0)=Nb(ran)
	  Nb(ran) = Nb(0)
	  call Combinaison(C-1,N-1,Cb+1,Nb+1)
      endif
end subroutine



Subroutine Force_lj_MD (xp,fp,box,upot)
      Implicit None
 
!     Calculate The Forces And Potential Energy
       Integer:: I,J
      real(double):: Dx,Dy,Dz,Ff,R2i,R6i,upot,Rcutsq
      real (double), dimension(1:n)::fxx,fyy,fzz,rxx,ryy,rzz
      real(double):: Box,hbox
      real (double), dimension(3,n)::xp,fp
!     Set Forces, Potential Energy And Pressure To Zero
 
      Do I = 1,N
         Fxx(I) = 0.0d0
         Fyy(I) = 0.0d0
         Fzz(I) = 0.0d0
      Enddo
 
      Upot = 0.0d0
     ! Press = 0.0d0
     ! box=11.0
      hbox=box*0.5d0
      Rcutsq=3.17d0

       rxx(:)=xp(1,1:N)
       ryy(:)=xp(2,1:N)
       rzz(:)=xp(3,1:N)




 !C     Loop Over All Particle Pairs
      Do I = 1,N - 1
         Do J = I + 1,N
 
!C     Calculate Distance And Perform Periodic
!C     Boundary Conditions
 
            Dx = Rxx(I) - Rxx(J)
            Dy = Ryy(I) - Ryy(J)
            Dz = Rzz(I) - Rzz(J)


           ! dx=dx-box*nint(dx/box)
           ! dy=dy-box*nint(dy/box)
           ! dx=dz-box*nint(dz/box)


            If (Dx.Gt.hbox) Then
               Dx = Dx - Box
            Elseif (Dx.Lt.-hbox) Then
               Dx = Dx + Box
            Endif
 
            If (Dy.Gt.hbox) Then
               Dy = Dy - Box
            Elseif (Dy.Lt. -hbox) Then
               Dy = Dy + Box
            Endif
 
            If (Dz.Gt.hbox) Then
               Dz = Dz - Box
            Elseif (Dz.Lt. -hbox) Then
               Dz = Dz + Box
            Endif
 
            R2i = Dx*Dx + Dy*Dy + Dz*Dz
            !write (*,*) 'r2i', r2i
           !if (r2i.lt.1e-1) stop

!C     Check If The Distance Is Within The Cutoff Radius
          If (R2i.Lt.Rcutsq) Then  
               R2i =1.0d0/R2i
               R6i = R2i*R2i*R2i
               Upot  = Upot+4.0d0*R6i*(R6i - 1.0d0)-0.121627043968990
               Ff    = 48.0d0*R6i*(R6i - 0.5d0)
               !Press = Press + Ff
               Ff    = Ff*R2i
 
               Fxx(I) = Fxx(I) + Ff*Dx
               Fyy(I) = Fyy(I) + Ff*Dy
               Fzz(I) = Fzz(I) + Ff*Dz
 
               Fxx(J) = Fxx(J) - Ff*Dx
               Fyy(J) = Fyy(J) - Ff*Dy
               Fzz(J) = Fzz(J) - Ff*Dz
 
            Endif
         Enddo
      Enddo

do i=1,N
  fp(1,i)=fxx(i)
  fp(2,i)=fyy(i)
  fp(3,i)=fzz(i)
enddo

! write(*,*) 'utot', upot/real(n)
!C     Scale The Pressure
 
    !  Press = Press/(3.0d0*Box*Box*Box)
 
      Return
    end subroutine


Subroutine energy_lj_MD (rxx,ryy,rzz,fp,Upot,box)
      Implicit None
 
!     Calculate The Forces And Potential Energy
       Integer:: I,J
      real(double):: Dx,Dy,Dz,Ff,R2i,R6i,upot,rcutsq
      real (double), dimension(n)::fxx,fyy,fzz,rxx,ryy,rzz
      real(double):: Box,hbox
  real (double), dimension(3,n)::xp,fp
!     Set Forces, Potential Energy And Pressure To Zero
 
      Do I = 1,N
         Fxx(I) = 0.0d0
         Fyy(I) = 0.0d0
         Fzz(I) = 0.0d0
      Enddo
 
      Upot = 0.0d0
     ! Press = 0.0d0
     ! box=
      hbox=0.5*box
      rcutsq=3.17d0

 !C     Loop Over All Particle Pairs
      Do I = 1,N - 1
         Do J = I + 1,N
 
!C     Calculate Distance And Perform Periodic
!C     Boundary Conditions
 
            Dx = Rxx(I) - Rxx(J)
            Dy = Ryy(I) - Ryy(J)
            Dz = Rzz(I) - Rzz(J)
 

             If (Dx.Gt.hbox) Then
               Dx = Dx - Box
            Elseif (Dx.Lt.-hbox) Then
               Dx = Dx + Box
            Endif
 
            If (Dy.Gt.hbox) Then
               Dy = Dy - Box
            Elseif (Dy.Lt. -hbox) Then
               Dy = Dy + Box
            Endif
 
            If (Dz.Gt.hbox) Then
               Dz = Dz - Box
            Elseif (Dz.Lt. -hbox) Then
               Dz = Dz + Box
            Endif 
            
            R2i = Dx*Dx + Dy*Dy + Dz*Dz

!C     Check If The Distance Is Within The Cutoff Radius
             If (R2i.Lt.Rcutsq) Then  
            !If ((R2i.Lt.Rcutsq).and.(r2i.gt.0.5d0)) Then
               R2i = 1.0d0/R2i
               R6i = R2i*R2i*R2i
               Upot  = Upot + 4.0d0*R6i*(R6i - 1.0d0)-0.121627035471831
               !Ff    = 48.0d0*R6i*(R6i - 0.5d0)
               !Press = Press + Ff
               !Ff    = Ff*R2i
               !Fxx(I) = Fxx(I) + Ff*Dx
               !Fyy(I) = Fyy(I) + Ff*Dy
               !Fzz(I) = Fzz(I) - Ff*Dz
 
               !Fxx(J) = Fxx(J) - Ff*Dx
               !Fyy(J) = Fyy(J) - Ff*Dy
               !Fzz(J) = Fzz(J) - Ff*Dz
 
            Endif
         Enddo
      Enddo

!write(*,*) 't', upot
!do i=1,N
!  fp(1,i)=fxx(i)
!  fp(2,i)=fyy(i)
!  fp(3,i)=fzz(i)
!enddo

 
!C     Scale The Pressure
 
    !  Press = Press/(3.0d0*Box*Box*Box)
 
      Return
    end subroutine



Subroutine CalculNorme (uqx,uqy,uqz,upx,upy,upz, N,norm)
  integer ::N
  real (double), dimension(0:N-1) ::uqx
  real (double), dimension (0:N-1) ::upx
  real (double), dimension(0:N-1) ::uqy
  real (double), dimension (0:N-1) ::upy
  real (double), dimension(0:N-1) ::uqz
  real (double), dimension (0:N-1) ::upz

  real (double) ::Norm
  integer ::j
  
  Norm = zero
   do j=0,N-1
     Norm=Norm+uqx(j)*uqx(j)+uqy(j)*uqy(j)+uqz(j)*uqz(j)+upz(j)*upz(j)+upx(j)*upx(j)+upy(j)*upy(j)
   enddo
    
   Norm=sqrt(Norm)
    
  end subroutine


Subroutine Normalizza(uqx,uqy,uqz,upx,upy,upz,N)
  integer ::N
  real (double), dimension (0:N-1):: qx   ! vector of position
  real (double), dimension (0:N-1):: px   ! vector of impulsion
  real (double), dimension (0:N-1):: uqx  ! linear displacement
  real (double), dimension (0:N-1):: upx  ! linear displacement
  real (double), dimension (0:N-1):: qy   ! vector of position
  real (double), dimension (0:N-1):: py   ! vector of impulsion
  real (double), dimension (0:N-1):: uqy  ! linear displacement
  real (double), dimension (0:N-1):: upy  ! linear displacement
  real (double), dimension (0:N-1):: qz   ! vector of position
  real (double), dimension (0:N-1):: pz   ! vector of impulsion
  real (double), dimension (0:N-1):: uqz  ! linear displacement
  real (double), dimension (0:N-1):: upz  ! linear displacement
  
  integer ::j
  real (double) ::Norm

    Norm=0
  do j=0,N-1
      Norm=Norm+uqx(j)*uqx(j)+upx(j)*upx(j)+uqy(j)*uqy(j)+upy(j)*upy(j)+uqz(j)*uqz(j)+upz(j)*upz(j)
  enddo

  do j=0,N-1
      uqx(j) = uqx(j)*1.0/sqrt(Norm)
      upx(j) = upx(j)*1.0/sqrt(Norm)
      uqy(j) = uqy(j)*1.0/sqrt(Norm)
      upy(j) = upy(j)*1.0/sqrt(Norm)
      uqz(j) = uqz(j)*1.0/sqrt(Norm)
      upz(j) = upz(j)*1.0/sqrt(Norm)
  enddo
end subroutine


Subroutine Hessien(qx,qy,qz,uqx,uqy,uqz,Vsecx,Vsecy,Vsecz,lapU)

  implicit none
  
  integer,parameter :: N=38
  real (double) , dimension (0:N-1)::qx,qy,qz
  real (double) , dimension (N)::Vsecx,vsecy,vsecz
  !real (double), dimension (3*N) ::x_vect
  integer:: m, i,a,b,j,ai, aj, bj,k
  real (double), dimension (N,3) ::pos_atom !, pos_atom2
  
  !real (double), dimension (N_pair)::V
  real (double)::U_tot, u, sigma, eps,r, h
  real (double), dimension (3*N,3*N)::Usec !,test
  !real (double), dimension (N)::hessx,hessy, hessz, Fxx,Fyy,Fzz, F2xx,F2yy,F2zz,&
  !rx,ry,rz, r2x,r2y,r2z
  real (double),dimension(N,N) ::rsq,dU,d2U
  real (double) :: x,lapU
  real (double), dimension(0:N-1) ::uqx
  real (double), dimension (0:N-1) ::upx
  real (double), dimension(0:N-1) ::uqy
  real (double), dimension (0:N-1) ::upy
  real (double), dimension(0:N-1) ::uqz
  real (double), dimension (0:N-1) ::upz
 real(double), dimension(3*n)::usec2

  eps=1
  sigma=1


  do m=0,N-1
      pos_atom(m+1,1)=qx(m)
      pos_atom(m+1,2)=qy(m)
      pos_atom(m+1,3)=qz(m)      
    enddo

  !calculate élément matrice hessienne
rsq(:,:)=0  
du(:,:)=0
d2u(:,:)=0

  !1)calcul distance interatomica
  do a=1,N-1
  do b=a+1,N
    do i=1,3
    rsq(a,b)=rsq(a,b)+(pos_atom(a,i)-pos_atom(b,i))**2
      enddo
    rsq(a,b)=sqrt(rsq(a,b))
    rsq(b,a)=rsq(a,b)
  enddo
  enddo

do a=1,N-1
  do b=a+1,N
    dU(a,b)=-4.d0*(12.0/(rsq(a,b)**13)-6.0/(rsq(a,b)**7))
    d2u(a,b)=4.d0*(156.0/(rsq(a,b)**14)-42.0/(rsq(a,b)**8))
    du(b,a)=du(a,b)
    d2u(b,a)=d2u(a,b)
  enddo
enddo

 Usec(:,:)=0
 
  !j'ai matrice des distances: chaque rsq(a,b) indique la distance pour pair d'atomes. je veux maintenant calcolare les elements della matrice hessiana.
  !case sur bloc diagonal: alfa égale beta (derivée par rapport à la variation position du meme atome)
do a=1,N
  do i=1,3
    ai=3*(a-1)+i
    do j=1,3
        aj=3*(a-1)+j

       do b=1,N
       if(b.ne.a) then
!       if ((b.ne.a).and.(j.ne.i)) then
       Usec(ai,aj)=Usec(ai,aj)+0.5*(((pos_atom(a,i)-pos_atom(b,i))*(pos_atom(a,j)-pos_atom(b,j))/(rsq(a,b)**2)*d2u(a,b)-&
                   du(a,b)*(pos_atom(a,i)-pos_atom(b,i))*(pos_atom(a,j)-pos_atom(b,j))/(rsq(a,b)**3)))
       
        if (j.eq.i) then     
                Usec(ai,aj)=Usec(ai,aj)+0.5*du(a,b)/(rsq(a,b))!+((pos_atom(a,i)-pos_atom(b,i))*(pos_atom(a,j)-pos_atom(b,j))/(rsq(a,b)**2)*d2u(a,b))-du(a,b)*(pos_atom(a,i)-pos_atom(b,i))*(pos_atom(a,j)-pos_atom(b,j))/(rsq(a,b)**3)
        endif 
       endif
       enddo 

       usec(aj,ai)=usec(ai,aj) 
       enddo !j
  enddo !i
 enddo !a

do i=1,3*n
lapU=lapU+usec(i,i)
enddo

!write (*,*) 'laplaciano', lapu
!stop



  !casella fouri diagonale: un solo termine
  do a=1,N-1
    do i=1,3
    ai=3*(a-1)+i
      do b=a+1,N
       do j=1,3
       bj=3*(b-1)+j
	 !if ((b.ne.a).and.(i.ne.j)) then
	   Usec(ai,bj)=-((pos_atom(a,i)-pos_atom(b,i))*(pos_atom(a,j)-pos_atom(b,j)))/(rsq(a,b)**2)*d2u(a,b)+&
                       (pos_atom(a,i)-pos_atom(b,i))*(pos_atom(a,j)-pos_atom(b,j))/(rsq(a,b)**3)*du(a,b)
            
           if (i.eq.j) then
	   Usec(ai,bj)=Usec(ai,bj)-du(a,b)/(rsq(a,b))
           endif
             !matrice simmetrica
           Usec(bj,ai)=Usec(ai,bj)
        enddo
      enddo
    enddo !i
  enddo !a

!do j=1,3*n
!usec2(j)=sum(usec(1:3*n,j))
!enddo

!x=0

!x=x+usec(8,3)+usec(8,6)

!do b=4,N
!x=x+usec(8,3*(b-1)+3)
!enddo

!do b=1,2
!x=x+usec(8,3*(b-1)+2)
!enddo

!do i=1,3*n
!write (*,'(a,e18.7, e18.7, e18.7)') 'cfr', usec(14,14), usec(7,26)
!enddo



!write(*,*) 'prova', usec(8,9),x

!do a=1,N
!do i = 1, 3
!  ai = 3*(a-1)+i
!      write(*,'(e18.7, e18.7, e18.7, e18.7, e18.7, e18.7, e18.7, e18.7, e18.7)') Usec(ai,:)
! end do
!enddo

!stop
!write (*,*) 'fatto'

 Vsecx(:)=0
 Vsecy(:)=0
 Vsecz(:)=0

!write (*,*) 'cambio'

!write (*,'(e18.7, e18.7, e18.7, e18.7, e18.7, e18.7)')  uqx(0),uqy(0),uqz(0),uqx(1),uqy(1),uqz(1),uqx(2),uqy(2),uqz(2)


!!!!ATTENZIONE ALL'ORDIEN DEGLI ARGOMENTI DI USEC
do a=1,3!N
 do b=1,3
     bj=3*(b-1)
     Vsecx(a)=Vsecx(a)+(Usec(bj+1,3*(a-1)+1)*uqx(b-1)+Usec(bj+2,3*(a-1)+1)*uqy(b-1)+Usec(bj+3,3*(a-1)+1)*uqz(b-1))
     !Vsecx(a)=Vsecx(a)+(Usec(3*(a-1)+1,bj+1)*uqx(b-1)+Usec(3*(a-1)+1,bj+2)*uqy(b-1)+Usec(3*(a-1)+1,bj+3)*uqz(b-1))

     Vsecy(a)=Vsecy(a)+(Usec(bj+1,3*(a-1)+2)*uqx(b-1)+Usec(bj+2,3*(a-1)+2)*uqy(b-1)+Usec(bj+3,3*(a-1)+2)*uqz(b-1))
     !Vsecy(a)=Vsecy(a)+(Usec(3*(a-1)+2,bj+1)*uqx(b-1)+Usec(3*(a-1)+2,bj+2)*uqy(b-1)+Usec(3*(a-1)+2,bj+3)*uqz(b-1))

     Vsecz(a)=Vsecz(a)+(Usec(bj+1,3*(a-1)+3)*uqx(b-1)+Usec(bj+2,3*(a-1)+3)*uqy(b-1)+Usec(bj+3,3*(a-1)+3)*uqz(b-1))
     !Vsecz(a)=Vsecz(a)+(Usec(3*(a-1)+3,bj+1)*uqx(b-1)+Usec(3*(a-1)+3,bj+2)*uqy(b-1)+Usec(3*(a-1)+3,bj+3)*uqz(b-1))
  enddo
   Vsecx(a)=Vsecx(a)
   Vsecy(a)=Vsecy(a)
   Vsecz(a)=Vsecz(a)
enddo

!write (*,*) 'vsecx', vsecx(:)!, vsecy(:), vsecz(:)
 


end subroutine

real(8) Function ran3()
      use random_art
      implicit none

      integer, parameter :: mbig=1000000000
      integer, parameter :: mseed=161803398
      integer, parameter :: mz=0
      real(8), parameter :: fac=1./mbig

      integer :: i,mj, mk, ii, k

      ! Any large mbig, and any smaller (but still large) mseed can be
      !  substituted for the above values.

      if(idum.lt.0.or.iff.eq.0)then
           iff=1
           mj=mseed-iabs(idum)
           mj=mod(mj,mbig)
           ma(55)=mj
           mk=1
           do i=1,54
             ii=mod(21*i,55)
             ma(ii)=mk
             mk=mj-mk
             if(mk.lt.mz)mk=mk+mbig
             mj=ma(ii)
           enddo
           do k=1,4
             do i=1,55
               ma(i)=ma(i)-ma(1+mod(i+30,55))
               if(ma(i).lt.mz)ma(i)=ma(i)+mbig
             enddo
           enddo
           inext=0
           inextp=31
           idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.mz)mj=mj+mbig
      ma(inext)=mj
      ran3=mj*fac
   end function ran3



subroutine testhessien(qx,qy,qz)


!integer::N
real(double), dimension(0:N-1):: qx,qy,qz,q1x,q1y,q1z,uqx,uqy,uqz
real(double), dimension (0:N-1)::hx,hy,hz,gx,gy,gz
real(double)::alfa,utot,x1,x2,x3,norm,lapu
real(double), dimension(3,N)::xp,xp1,fp,fp1,uprova,ff
real(double), dimension(3*n)::usec2
integer::i

alfa=1e-2
call random_number(x1)
hx(0:N-1)=x1
call random_number(x2)
hy(0:N-1)=x2
call random_number(x3)
hz(0:N-1)=x3 
gx(0:N-1)=0
gy(0:N-1)=0
gz(0:N-1)=0
call normalizza(hx,hy,hz,gx,gy,gz,N)
call calculnorme(hx,hy,hz,gx,gy,gz,n,Norm)
!write (*,*) 'norma h' , norm
call force(qx,qy,qz,ff,utot)
!write(*,*)'forzamia', ff(1,:)
!xp(1,1:N)=qx(0:N-1)
!xp(2,1:N)=qy(0:N-1)
!xp(3,1:N)=qz(0:N-1)
!call calfoljc(xp,fp)
!write(*,*)'forza1', fp(1,:)
q1x(:)=qx(:)+alfa*hx(:)
q1y(:)=qy(:)+alfa*hy(:)
q1z(:)=qz(:)+alfa*hz(:)

!xp(1,1:N)=q1x(0:N-1)
!xp(2,1:N)=q1y(0:N-1)
!xp(3,1:N)=q1z(0:N-1)

call force(q1x,q1y,q1z,fp1,utot)
!call calfoljc(xp,fp1)
!write(*,*) 'forza 2',fp1(1,:) 
do i=0,N-1
 uprova(1,i+1)=(fp1(1,i+1)-fp(1,i+1))/alfa
 uprova(2,i+1)=(fp1(2,i+1)-fp(2,i+1))/alfa
 uprova(3,i+1)=(fp1(3,i+1)-fp(3,i+1))/alfa
enddo

!write(*,*) 'uprova ', uprova (1,:)

call hessien(qx,qy,qz,hx,hy,hz,Vsecx,Vsecy,Vsecz,lapu)



end subroutine testhessien

Function genrand()
 real (double) :: genrand 
  real (double) :: x
  call random_number(x)
  genrand=x

endfunction

Subroutine Force(rxx,ryy,rzz,fp,Utot)
        Implicit None
  !     Calculate The Forces And Potential Energy
    !     Include 'system.inc'
        Integer:: I,J
       !integer, parameter::N=38
        real (double):: Dx,Dy,Dz,F,R2i,R6i,eta,sigma, Utot
        !real (double)::rx,ry,rz,
  real (double),dimension(0:N-1):: Rxx,Ryy,Rzz,fyy
  real (double),dimension(0:N-1):: fxx
  real (double),dimension(0:N-1):: fzz
  real (double), dimension(3,N)::fp
sigma=1
 eta=1
utot=0     

! Do I = 1,N
         Fxx(:) = 0.0d0
         Fyy(:) = 0.0d0
         Fzz(:) = 0.0d0

!write(*,*) 'xx',rxx
      ! Enddo
  !    Loop Over All Particle Pairs
        Do I = 0,N - 2
	 Do J =i+1,N-1
	     Dx = Rxx(i) - Rxx(J)
	     Dy = Ryy(i) - Ryy(J)
	     Dz = Rzz(i) - Rzz(J)
	     R2i = Dx*Dx + Dy*Dy + Dz*Dz

  ! write (*,*) 'dx', dx
	     
	       R2i = 1.0/R2i
	       R6i = R2i*R2i*R2i
	       F = 48*eta*R2i*(R6i*(sigma**12) - 0.5*(sigma**6))
	       F = F*R6i
               Utot =Utot+4*eta*(R6i*R6i*(sigma**12) - R6i*(sigma**6))
  !write (*,*) 'energ', F
   
	       Fxx(i) = Fxx(i) + F*Dx
	       Fyy(i) = Fyy(i) + F*Dy
	       Fzz(i) = Fzz(i) + F*Dz
  
	       Fxx(J) = Fxx(J) - F*Dx
	       Fyy(J) = Fyy(J) - F*Dy
	       Fzz(J) = Fzz(J) - F*Dz
  
	 Enddo

        Enddo

!write (*,*)'forse', Fxx(:), Fyy(:), Fzz(:)

do i=0,N-1
  fp(1,i+1)=fxx(i)
  fp(2,i+1)=fyy(i)
  fp(3,i+1)=fzz(i)
enddo



        Return
End subroutine

SUBROUTINE shell(n,a)


INTEGER :: n
integer, dimension(n) ::a
!Sorts an array a(1:n) into ascending numerical order by Shell’s method (diminishing increment sort). n is input; a is replaced on output by its sorted rearrangement.
INTEGER i,j,inc
REAL v
inc=1 !Determine the starting increment.
1 inc=3*inc+1
if(inc.le.n)goto 1
2 continue !Loop over the partial sorts.
    inc=inc/3
    do i=inc+1,n !Outer loop of straight insertion.
       v=a(i)
        j=i
     3 if(a(j-inc).gt.v)then !Inner loop of straight insertion.
          a(j)=a(j-inc)
          j=j-inc
          if(j.le.inc)goto 4
     goto 3
     endif
   4 a(j)=v
   enddo 
if(inc.gt.1)goto 2
return
end subroutine shell

end module teledyn_module  




Subroutine fq4(N,xij,yij,zij,dij,fq4x,fq4y,fq4z)
  USE T_kind_param_m, ONLY:  double

        real(double)  OUT,gmaxmax,gminmin,rcut,delta_r,apol,bpol,cpol
        integer NEMAX,nwmax,N,nsize,i,k,j
	parameter(OUT=1D40,NEMAX=2500,nwmax=50)
	parameter(gmaxmax=1D10,gminmin=1.D-50)
        parameter(rcut=1.392,delta_r=0.1,apol=-3.D0/16.D0/delta_r**5)
        parameter(bpol=5.D0/8.D0/delta_r**3,cpol=-15.D0/16.D0/delta_r)

	real(double) x(N),y(N),z(N)
	real(double)  fwx(N),fwy(N),fwz(N)
	real(double)  fQ0x(N),fQ0y(N),fQ0z(N)
	real(double) fQRx(4,N),fQRy(4,N),fQRz(4,N)
	real(double) fQIx(4,N),fQIy(4,N),fQIz(4,N)
	real(double) fq4x(N),fq4y(N),fq4z(N)
	real(double) COSPHI(4), SINPHI(4)
	real(double) q4coef(3,0:4)
	real(double) qr(4),qi(4)
	common/q4coef/q4coef
	real(double)  xij(N,N)
	real(double)  yij(N,N)
	real(double)  zij(N,N)
	real(double)  dij(N,N)
        real(double) Q0,DX,DY,DZ,WDB,dist,weight,dweight,w2
        real(double) costh,costh2,costh4,sinth,sinth2,scth,twocph
        real(double) dctdx,dctdy,dctdz,dn
        real(double) dcpdx,dcpdy,dcpdz
        real(double) dspdx,dspdy,dspdz
        real(double) wq0,dwq0,wq0_1,wq1,dwq1,wq1_1,wq2,dwq2,wq2_1,wq3,dwq3,wq3_1
        real(double) dfrx,dfry,dfrz
        real(double) dc2pdx,dc2pdy,dc2pdz
        real(double) dfix,dfiy,dfiz
        real(double) ds2pdx,ds2pdy,ds2pdz
        real(double) dcp,dcp2
        real(double) dc3pdx,dc3pdy,dc3pdz
        real(double) dc4pdx,dc4pdy,dc4pdz
        real(double) ds4pdx,ds4pdy,ds4pdz
        real(double) ds3pdx,ds3pdy,ds3pdz
        real(double) WQ4,WQ4_1,dWQ4,WNB,Q4

!	common/pairdist/xij,yij,zij,dij

	do i=1,N
	   fq4x(i)=0.D0
	   fq4y(i)=0.D0
	   fq4z(i)=0.D0
	   fwx(i)=0.D0
	   fwy(i)=0.D0
	   fwz(i)=0.D0
	   fq0x(i)=0.D0
	   fq0y(i)=0.D0
	   fq0z(i)=0.D0
	   do k=1,4
	      fqRx(k,i)=0.D0
	      fqRy(k,i)=0.D0
	      fqRz(k,i)=0.D0
	      fqIx(k,i)=0.D0
	      fqIy(k,i)=0.D0
	      fqIz(k,i)=0.D0
	   enddo
	enddo

	Q0=0.D0
	DO I=1,4
	   QI(I)=0.D0
	   QR(I)=0.D0
	ENDDO

	WDB=0.D0
	DO 30 J=1,N
	   DO 40 K=J+1,N
	      DX=xij(j,k)
	      DY=yij(j,k)
	      DZ=zij(j,k)
	      DIST=dij(j,k)

	      if (dist.le.rcut-delta_r) then
		 weight=1.D0
		 dweight=0.D0
	      elseif (dist.ge.rcut+delta_r) then
		 weight=0.D0
		dweight=0.D0
	      else
		 weight=0.5D0+apol*(dist-rcut)**5 + bpol*(dist-rcut)**3 + cpol*(dist-rcut)
		dweight=5.D0*apol*(dist-rcut)**4 + 3.D0*bpol*(dist-rcut)**2 + cpol
	      endif
	      WDB=WDB+WEIGHT

	      w2=-dweight/dist
	      fwx(j)=fwx(j)+dx*w2
	      fwy(j)=fwy(j)+dy*w2
	      fwz(j)=fwz(j)+dz*w2
	      fwx(k)=fwx(k)-dx*w2
	      fwy(k)=fwy(k)-dy*w2
	      fwz(k)=fwz(k)-dz*w2

	      COSTH  = DZ/DIST
	      COSTH2 = COSTH  * COSTH
	      COSTH4 = COSTH2 * COSTH2
	      SINTH2 = 1. - COSTH2
	      SINTH  = SQRT(SINTH2)
	      SCTH   = SINTH  * COSTH
	      IF (SINTH .EQ. 0.) THEN
		 COSPHI(1) = 1.
		 SINPHI(1) = 1.
	      ELSE
		 COSPHI(1) = DX/DIST/SINTH
		 SINPHI(1) = DY/DIST/SINTH
	      ENDIF
	      TWOCPH    = 2.*COSPHI(1)
	      COSPHI(2) = TWOCPH*COSPHI(1)-1.
	      SINPHI(2) = TWOCPH*SINPHI(1)
	      COSPHI(3) = TWOCPH*COSPHI(2)-COSPHI(1)
	      SINPHI(3) = TWOCPH*SINPHI(2)-SINPHI(1)
	      COSPHI(4) = TWOCPH*COSPHI(3)-COSPHI(2)
	      SINPHI(4) = TWOCPH*SINPHI(3)-SINPHI(2)

!       derivatives of cos theta and cos/sin phi

	      dctdx=dx*dz/dist**3
	      dctdy=dy*dz/dist**3
	      dctdz=-(1.D0-(dz/dist)**2)/dist

	      dn=dx**2+dy**2
	      dcpdx=-SINPHI(1)*dy/dn
	      dcpdy= COSPHI(1)*dy/dn
	      dcpdz= 0.D0
	      dspdx= dcpdy
	      dspdy=-COSPHI(1)*dx/dn
	      dspdz= 0.D0

!       accumulates Q0, QR, QI and their gradients

!       Q0

	      WQ0=Q4COEF(1,0)*COSTH4+Q4COEF(2,0)*COSTH2+Q4COEF(3,0)
	      Q0=Q0+WQ0*WEIGHT

	      WQ0_1=w2*WQ0
	      dWQ0=(4.D0*Q4COEF(1,0)*COSTH2+2.D0*Q4COEF(2,0))*COSTH

	      fq0x(j)=fq0x(j)+dx*WQ0_1+WEIGHT*dWQ0*dctdx
	      fq0y(j)=fq0y(j)+dy*WQ0_1+WEIGHT*dWQ0*dctdy
	      fq0z(j)=fq0z(j)+dz*WQ0_1+WEIGHT*dWQ0*dctdz
	      fq0x(k)=fq0x(k)-dx*WQ0_1-WEIGHT*dWQ0*dctdx
	      fq0y(k)=fq0y(k)-dy*WQ0_1-WEIGHT*dWQ0*dctdy
	      fq0z(k)=fq0z(k)-dz*WQ0_1-WEIGHT*dWQ0*dctdz

!       QR(1) and QI(1)

	      WQ1=(Q4COEF(1,1)*COSTH2+Q4COEF(2,1))*SCTH
	      QR(1)=QR(1)+COSPHI(1)*WQ1*WEIGHT
	      QI(1)=QI(1)+SINPHI(1)*WQ1*WEIGHT

	      WQ1_1=w2*WQ1
	      dWQ1=SINTH*(3.D0*Q4COEF(1,1)*COSTH2+Q4COEF(2,1))
              IF (SINTH.EQ.0) THEN
                 dWQ1=0.D0
              ELSE
                 dWQ1=dWQ1-(Q4COEF(1,1)*COSTH2+Q4COEF(2,1))*COSTH2/SINTH
              ENDIF
	      dfrx=dWQ1*dctdx*COSPHI(1)+WQ1*dcpdx
	      dfry=dWQ1*dctdy*COSPHI(1)+WQ1*dcpdy
	      dfrz=dWQ1*dctdz*COSPHI(1)+WQ1*dcpdz

	      fqrx(1,j)=fqrx(1,j)+dx*COSPHI(1)*WQ1_1+WEIGHT*dfrx
	      fqry(1,j)=fqry(1,j)+dy*COSPHI(1)*WQ1_1+WEIGHT*dfry
	      fqrz(1,j)=fqrz(1,j)+dz*COSPHI(1)*WQ1_1+WEIGHT*dfrz
	      fqrx(1,k)=fqrx(1,k)-dx*COSPHI(1)*WQ1_1-WEIGHT*dfrx
	      fqry(1,k)=fqry(1,k)-dy*COSPHI(1)*WQ1_1-WEIGHT*dfry
	      fqrz(1,k)=fqrz(1,k)-dz*COSPHI(1)*WQ1_1-WEIGHT*dfrz

	      dfix=dWQ1*dctdx*SINPHI(1)+WQ1*dspdx
	      dfiy=dWQ1*dctdy*SINPHI(1)+WQ1*dspdy
	      dfiz=dWQ1*dctdz*SINPHI(1)+WQ1*dspdz

	      fqix(1,j)=fqix(1,j)+dx*SINPHI(1)*WQ1_1+WEIGHT*dfix
	      fqiy(1,j)=fqiy(1,j)+dy*SINPHI(1)*WQ1_1+WEIGHT*dfiy
	      fqiz(1,j)=fqiz(1,j)+dz*SINPHI(1)*WQ1_1+WEIGHT*dfiz
	      fqix(1,k)=fqix(1,k)-dx*SINPHI(1)*WQ1_1-WEIGHT*dfix
	      fqiy(1,k)=fqiy(1,k)-dy*SINPHI(1)*WQ1_1-WEIGHT*dfiy
	      fqiz(1,k)=fqiz(1,k)-dz*SINPHI(1)*WQ1_1-WEIGHT*dfiz

!       QR(2) and QI(2)

	      WQ2=Q4COEF(1,2)*COSTH4+Q4COEF(2,2)*COSTH2+Q4COEF(3,2)
	      QR(2)=QR(2)+COSPHI(2)*WQ2*WEIGHT
	      QI(2)=QI(2)+SINPHI(2)*WQ2*WEIGHT

	      WQ2_1=w2*WQ2
	      dWQ2=(4.D0*Q4COEF(1,2)*COSTH2+2.D0*Q4COEF(2,2))*COSTH

	      dc2pdx=4.D0*COSPHI(1)*dcpdx
	      dc2pdy=4.D0*COSPHI(1)*dcpdy
	      dc2pdz=4.D0*COSPHI(1)*dcpdz

	      dfrx=dWQ2*dctdx*COSPHI(2)+WQ2*dc2pdx
	      dfry=dWQ2*dctdy*COSPHI(2)+WQ2*dc2pdy
	      dfrz=dWQ2*dctdz*COSPHI(2)+WQ2*dc2pdz

	      fqrx(2,j)=fqrx(2,j)+dx*COSPHI(2)*WQ2_1+WEIGHT*dfrx
	      fqry(2,j)=fqry(2,j)+dy*COSPHI(2)*WQ2_1+WEIGHT*dfry
	      fqrz(2,j)=fqrz(2,j)+dz*COSPHI(2)*WQ2_1+WEIGHT*dfrz
	      fqrx(2,k)=fqrx(2,k)-dx*COSPHI(2)*WQ2_1-WEIGHT*dfrx
	      fqry(2,k)=fqry(2,k)-dy*COSPHI(2)*WQ2_1-WEIGHT*dfry
	      fqrz(2,k)=fqrz(2,k)-dz*COSPHI(2)*WQ2_1-WEIGHT*dfrz

	      ds2pdx=2.D0*(COSPHI(1)*dspdx+SINPHI(1)*dcpdx)
	      ds2pdy=2.D0*(COSPHI(1)*dspdy+SINPHI(1)*dcpdy)
	      ds2pdz=2.D0*(COSPHI(1)*dspdz+SINPHI(1)*dcpdz)

	      dfix=dWQ2*dctdx*SINPHI(2)+WQ2*ds2pdx
	      dfiy=dWQ2*dctdy*SINPHI(2)+WQ2*ds2pdy
	      dfiz=dWQ2*dctdz*SINPHI(2)+WQ2*ds2pdz

	      fqix(2,j)=fqix(2,j)+dx*SINPHI(2)*WQ2_1+WEIGHT*dfix
	      fqiy(2,j)=fqiy(2,j)+dy*SINPHI(2)*WQ2_1+WEIGHT*dfiy
	      fqiz(2,j)=fqiz(2,j)+dz*SINPHI(2)*WQ2_1+WEIGHT*dfiz
	      fqix(2,k)=fqix(2,k)-dx*SINPHI(2)*WQ2_1-WEIGHT*dfix
	      fqiy(2,k)=fqiy(2,k)-dy*SINPHI(2)*WQ2_1-WEIGHT*dfiy
	      fqiz(2,k)=fqiz(2,k)-dz*SINPHI(2)*WQ2_1-WEIGHT*dfiz

!       QR(3) and QI(3)

	      WQ3=Q4COEF(1,3)*SCTH*SINTH2
	      QR(3)=QR(3)+COSPHI(3)*WQ3*WEIGHT
	      QI(3)=QI(3)+SINPHI(3)*WQ3*WEIGHT

	      WQ3_1=w2*WQ3
	      dWQ3=Q4COEF(1,3)*SINTH*(1.D0-4.D0*COSTH2)

	      dcp=4.D0*COSPHI(1)**2-1.D0
	      dc3pdx=3.D0*dcpdx*dcp
	      dc3pdy=3.D0*dcpdy*dcp
	      dc3pdz=3.D0*dcpdz*dcp

	      dfrx=dWQ3*dctdx*COSPHI(3)+WQ3*dc3pdx
	      dfry=dWQ3*dctdy*COSPHI(3)+WQ3*dc3pdy
	      dfrz=dWQ3*dctdz*COSPHI(3)+WQ3*dc3pdz

	      fqrx(3,j)=fqrx(3,j)+dx*COSPHI(3)*WQ3_1+WEIGHT*dfrx
	      fqry(3,j)=fqry(3,j)+dy*COSPHI(3)*WQ3_1+WEIGHT*dfry
	      fqrz(3,j)=fqrz(3,j)+dz*COSPHI(3)*WQ3_1+WEIGHT*dfrz
	      fqrx(3,k)=fqrx(3,k)-dx*COSPHI(3)*WQ3_1-WEIGHT*dfrx
	      fqry(3,k)=fqry(3,k)-dy*COSPHI(3)*WQ3_1-WEIGHT*dfry
	      fqrz(3,k)=fqrz(3,k)-dz*COSPHI(3)*WQ3_1-WEIGHT*dfrz

	      ds3pdx=dcp*dspdx+8.D0*SINPHI(1)*COSPHI(1)*dcpdx
	      ds3pdy=dcp*dspdy+8.D0*SINPHI(1)*COSPHI(1)*dcpdy
	      ds3pdz=dcp*dspdz+8.D0*SINPHI(1)*COSPHI(1)*dcpdz

	      dfix=dWQ3*dctdx*SINPHI(3)+WQ3*ds3pdx
	      dfiy=dWQ3*dctdy*SINPHI(3)+WQ3*ds3pdy
	      dfiz=dWQ3*dctdz*SINPHI(3)+WQ3*ds3pdz

	      fqix(3,j)=fqix(3,j)+dx*SINPHI(3)*WQ3_1+WEIGHT*dfix
	      fqiy(3,j)=fqiy(3,j)+dy*SINPHI(3)*WQ3_1+WEIGHT*dfiy
	      fqiz(3,j)=fqiz(3,j)+dz*SINPHI(3)*WQ3_1+WEIGHT*dfiz
	      fqix(3,k)=fqix(3,k)-dx*SINPHI(3)*WQ3_1-WEIGHT*dfix
	      fqiy(3,k)=fqiy(3,k)-dy*SINPHI(3)*WQ3_1-WEIGHT*dfiy
	      fqiz(3,k)=fqiz(3,k)-dz*SINPHI(3)*WQ3_1-WEIGHT*dfiz

!       QR(4) and QI(4)

	      WQ4=Q4COEF(1,4)*SINTH2*SINTH2
	      QR(4)=QR(4)+COSPHI(4)*WQ4*WEIGHT
	      QI(4)=QI(4)+SINPHI(4)*WQ4*WEIGHT

	      WQ4_1=w2*WQ4
	      dWQ4=-Q4COEF(1,4)*4.D0*SINTH2*COSTH

	      dcp=4.D0*COSPHI(1)*(2.D0*COSPHI(1)**2-1.D0)
	      dc4pdx=4.D0*dcpdx*dcp
	      dc4pdy=4.D0*dcpdy*dcp
	      dc4pdz=4.D0*dcpdz*dcp

	      dfrx=dWQ4*dctdx*COSPHI(4)+WQ4*dc4pdx
	      dfry=dWQ4*dctdy*COSPHI(4)+WQ4*dc4pdy
	      dfrz=dWQ4*dctdz*COSPHI(4)+WQ4*dc4pdz

	      fqrx(4,j)=fqrx(4,j)+dx*COSPHI(4)*WQ4_1+WEIGHT*dfrx
	      fqry(4,j)=fqry(4,j)+dy*COSPHI(4)*WQ4_1+WEIGHT*dfry
	      fqrz(4,j)=fqrz(4,j)+dz*COSPHI(4)*WQ4_1+WEIGHT*dfrz
	      fqrx(4,k)=fqrx(4,k)-dx*COSPHI(4)*WQ4_1-WEIGHT*dfrx
	      fqry(4,k)=fqry(4,k)-dy*COSPHI(4)*WQ4_1-WEIGHT*dfry
	      fqrz(4,k)=fqrz(4,k)-dz*COSPHI(4)*WQ4_1-WEIGHT*dfrz

	      dcp2=4.D0*SINPHI(1)*(6.D0*COSPHI(1)**2-1.D0)
	      ds4pdx=dspdx*dcp+dcpdx*dcp2
	      ds4pdy=dspdy*dcp+dcpdy*dcp2
	      ds4pdz=dspdz*dcp+dcpdz*dcp2

	      dfix=dWQ4*dctdx*SINPHI(4)+WQ4*ds4pdx
	      dfiy=dWQ4*dctdy*SINPHI(4)+WQ4*ds4pdy
	      dfiz=dWQ4*dctdz*SINPHI(4)+WQ4*ds4pdz

	      fqix(4,j)=fqix(4,j)+dx*SINPHI(4)*WQ4_1+WEIGHT*dfix
	      fqiy(4,j)=fqiy(4,j)+dy*SINPHI(4)*WQ4_1+WEIGHT*dfiy
	      fqiz(4,j)=fqiz(4,j)+dz*SINPHI(4)*WQ4_1+WEIGHT*dfiz
	      fqix(4,k)=fqix(4,k)-dx*SINPHI(4)*WQ4_1-WEIGHT*dfix
	      fqiy(4,k)=fqiy(4,k)-dy*SINPHI(4)*WQ4_1-WEIGHT*dfiy
	      fqiz(4,k)=fqiz(4,k)-dz*SINPHI(4)*WQ4_1-WEIGHT*dfiz

 40	   CONTINUE
 30	CONTINUE


!	Q0=0.D0
!	DO I=4,4
!	   QI(I)=0.D0
!	   QR(I)=0.D0
!	ENDDO

	WNB = DSQRT( Q0 * Q0 + ( QR(1)*QR(1) + QI(1)*QI(1) ) + ( QR(2)*QR(2) + QI(2)*QI(2) ) + ( QR(3)*QR(3) + QI(3)*QI(3) ) + ( QR(4)*QR(4) + QI(4)*QI(4) ) )

	Q4=WNB/WDB

	do i=1,n
	   fq4x(i)=Q0*fq0x(i)
	   fq4y(i)=Q0*fq0y(i)
	   fq4z(i)=Q0*fq0z(i)
	   do k=1,4
	      fq4x(i)=fq4x(i)+QR(k)*fqrx(k,i)+QI(k)*fqix(k,i)
	      fq4y(i)=fq4y(i)+QR(k)*fqry(k,i)+QI(k)*fqiy(k,i)
	      fq4z(i)=fq4z(i)+QR(k)*fqrz(k,i)+QI(k)*fqiz(k,i)
	   enddo
	   fq4x(i)=fq4x(i)/WDB/WNB-WNB*fwx(i)/WDB**2
	   fq4y(i)=fq4y(i)/WDB/WNB-WNB*fwy(i)/WDB**2
	   fq4z(i)=fq4z(i)/WDB/WNB-WNB*fwz(i)/WDB**2
	enddo

	return
end Subroutine


Function q4(N,xij,yij,zij,dij)
    use T_kind_param_m, ONLY:  double
 	implicit none
        real(double)  OUT,gmaxmax,gminmin,rcut,delta_r,apol,bpol,cpol
        integer NEMAX,nwmax,N,nsize,i,k,j
	parameter(OUT=1D40,NEMAX=2500,nwmax=50)
	parameter(gmaxmax=1e10,gminmin=1.D-50)
        parameter(rcut=1.392,delta_r=0.1,apol=-3.D0/16.D0/delta_r**5)
        parameter(bpol=5.D0/8.D0/delta_r**3,cpol=-15.D0/16.D0/delta_r)

	real(double) x(N),y(N),z(N)
	real(double) q4coef(3,0:4)
	real(double) qr(4),qi(4)
        real(double) WDB,WNB,DX,DY,DZ,DIST,weight
	common/q4coef/q4coef
	real(double) xij(N,N)
	real(double) yij(N,N)
	real(double) zij(N,N)
	real(double) dij(N,N)
!	common/pairdist/xij,yij,zij,dij
        real(double) WQ4,Q4,Q0

	Q0=0.D0
	DO I=1,4
	   QI(I)=0.D0
	   QR(I)=0.D0
	ENDDO

      ! write (*,*) n
   
	WDB=0.D0
	DO  J=1,N
	   DO  K=J+1,N
	      DX=xij(j,k)
             ! write (*,*) 'dx', dx
	      DY=yij(j,k)
	      DZ=zij(j,k)
	      DIST=dij(j,k)
              !write (*,*) 'dx', dist
	      if (dist.le.rcut-delta_r) then
		 weight=1.D0
	      elseif (dist.ge.rcut+delta_r) then
		 weight=0.D0
	      else
		 weight=0.5D0+apol*(dist-rcut)**5 + bpol*(dist-rcut)**3 + cpol*(dist-rcut)
	      endif
  	      WDB=WDB+WEIGHT
              
	      CALL EVASH4(DX/DIST,DY/DIST,DZ/DIST, Q0, QR, QI, WEIGHT)


            end do
          end do


!	Q0=0.D0
!	DO I=4,4
!	   QI(I)=0.D0
!	   QR(I)=0.D0
!	ENDDO

!write (*,*) 'wdb', wdb, wnb, q0, qr(:), qi(:)
	Q4 = WQ4(WDB,WNB,Q0,QR,QI)


	return

end Function


Subroutine SHINIT
        use T_kind_param_m, ONLY:  double
!       *** Calculate coefficients for Q4, Q6, W4 and W6
!       Notebook JvD A92
 
	COMMON/Q4COEF/ Q4COEF
 	COMMON/Q6COEF/ Q6COEF
	INTEGER    M, I
	double precision FACT, FACS
 
	double precision  Q4COEF(3,0:4),Q6COEF(4,0:6)
 

!       *** It would be nice indeed to know the exact formula!
 
	DATA      Q4COEF /4.375,  -3.75,  0.375,     -17.5,     7.5,   0.0,    -52.5,    60.0,  -7.5,    -105.0,     0.0,   0.0,    105.0,     0.0,   0.0/
	!       *** It would be nice indeed to know the exact formula!
 
	DATA      Q6COEF / 14.4375, -19.6875, 6.5625, -0.3125,-86.625,78.75,-13.125,0.0,433.125,-236.25,13.125,0.0,-1732.5,472.5,0.0,0.0,5197.5,-472.5,0.0,0.0,-10395.0,0.0,0.0,0.0,10395.0,0.0, 0.0,0.0/

!       *** Now multiply the coefficients with sqrt(2.0* (l-m)!/(l+m)! )
	
	DO 99 M=1, 4
	   DO 98 I=1, 3
	      FACS=DSQRT(2.0*FACT(4-M)/FACT(4+M))
	      Q4COEF(I,M)=Q4COEF(I,M)*FACS
 98	   CONTINUE
 99	CONTINUE	

	DO 89 M=1, 6
	   DO 88 I=1, 4
	      FACS=DSQRT(2.0*FACT(6-M)/FACT(6+M))
	      Q6COEF(I,M)=Q6COEF(I,M)*FACS
 88	   CONTINUE
 89	CONTINUE

	RETURN
end subroutine


Subroutine EVASH4(X, Y, Z, Q0, QR, QI, WEIGHT)
        USE T_kind_param_m, ONLY:  double
!       *** Evaluate the spherical harmonics of degree 4 ********************
!       
!       REAL(double)   X,Y,Z    A vector on the unit-sphere to be processed
!       DOUBLE PRECISION
!       Q0       Accumulator of P40(cos(THETA))
!       QR(4)    Accumulator of P4m(cos(THETA))*cos(PHI)
!       QI(4)    Accumulator of P4m(cos(THETA))*sin(PHI)
!       
!       Note that the actual spherical harmonics Y4m differ from Q0
!       and Q (= QR + i*QI) by a factor
!       sqrt((9/8*pi)
!       This will be taken care of in FUNCTION Q4
!       
!       Angles THETA and PHI are defined in the regular way:
!       
!       X = cos(PHI)sin(THETA)
!       Y = sin(PHI)sin(THETA)
!       Z = cos(THETA)
!       
!       Notebook JvD A38 and A85
!       
!
	
	COMMON/Q4COEF/ Q4COEF
	double precision   Q4COEF(3,0:4)
	
	double precision  X, Y, Z,WEIGHT
	DOUBLE PRECISION  Q0, QR(4), QI(4)
	
	double precision COSTH, COSTH2, COSTH4, SINTH, SINTH2, SCTH
	double precision TWOCPH, COSPHI(4), SINPHI(4), TEMP
	
	COSTH  = Z
	COSTH2 = COSTH  * COSTH
	COSTH4 = COSTH2 * COSTH2
	SINTH2 = 1. - COSTH2
	SINTH  = SQRT(SINTH2)
	SCTH   = SINTH  * COSTH
	
!       *** Is THETA = 0 ? Then PHI is irrelevant *****************
	
	IF (SINTH .EQ. 0.) THEN
	   COSPHI(1) = 1.
	   SINPHI(1) = 1.
	ELSE
	   COSPHI(1) = X/SINTH
	   SINPHI(1) = Y/SINTH
	ENDIF
	TWOCPH    = 2.*COSPHI(1)
	COSPHI(2) = TWOCPH*COSPHI(1)-1.
	SINPHI(2) = TWOCPH*SINPHI(1)
	COSPHI(3) = TWOCPH*COSPHI(2)-COSPHI(1)
	SINPHI(3) = TWOCPH*SINPHI(2)-SINPHI(1)
	COSPHI(4) = TWOCPH*COSPHI(3)-COSPHI(2)
	SINPHI(4) = TWOCPH*SINPHI(3)-SINPHI(2)
	
!       *** Now the spherical harmonics are calculated and summed with Q
!       
!       This part of the subroutine would have been easier to understand
!       when complex numbers would have been used. However, COMPLEX*16
!C       variables are not included in the FORTRAN 77 standard.
	
	Q0    = Q0 + (Q4COEF(1,0)*COSTH4+Q4COEF(2,0)*COSTH2 + Q4COEF(3,0))*WEIGHT

!	TEMP  = (Q4COEF(1,1) * COSTH2 + Q4COEF(2,1))
	TEMP  = (Q4COEF(1,1) * COSTH2 + Q4COEF(2,1)) * SCTH
!	TEMP  = 1E6*SCTH
 

	QR(1) = QR(1) + COSPHI(1)*TEMP*WEIGHT

!	QR(1) = QR(1) + TEMP*WEIGHT
	QI(1) = QI(1) + SINPHI(1)*TEMP*WEIGHT
!	QI(1) = QI(1) + TEMP*WEIGHT

	TEMP  = Q4COEF(1,2)*COSTH4+Q4COEF(2,2)*COSTH2+Q4COEF(3,2)
	QR(2) = QR(2) + COSPHI(2)*TEMP*WEIGHT
	QI(2) = QI(2) + SINPHI(2)*TEMP*WEIGHT

	TEMP  = Q4COEF(1,3) * SCTH * SINTH2
	QR(3) = QR(3) + COSPHI(3)*TEMP*WEIGHT
	QI(3) = QI(3) + SINPHI(3)*TEMP*WEIGHT

	TEMP  = Q4COEF(1,4) * SINTH2 * SINTH2
	QR(4) = QR(4) + COSPHI(4)*TEMP*WEIGHT
	QI(4) = QI(4) + SINPHI(4)*TEMP*WEIGHT

	RETURN
end Subroutine


!	double precision 
Function WQ4(WDB,WNB,Q0,QR,QI)
          USE T_kind_param_m, ONLY:  double
!       *** The order parameter Q4 is calculated from Q0 and Q *************
!       
!       INTEGER  NB      Number of bonds accumulated
!       DOUBLE PRECISION
!       Q0      Accumulator of P40(cos(THETA))
!       QR(4)   Accumulator of P4m(cos(THETA))*cos(PHI)
!       QI(4)   Accumulator of P4m(cos(THETA))*sin(PHI)
!       
!       In the expression for Q4, the factor pi*9/4 disappears
!       The summation has to be performed over the squares of Q4m, with
!       m running from -4 to +4. However, Q4m and Q4-m are equal when
!       squared. So the summation is done for positive m and a factor
!       2 is introduced.
!       
!       ********************************************************************
	
	INTEGER   NB
	real(double)  Q0, QR(4), QI(4),WNB,WDB,WQ4

	WNB = DSQRT( Q0 * Q0    + ( QR(1)*QR(1) + QI(1)*QI(1) )   + ( QR(2)*QR(2) + QI(2)*QI(2) )    + ( QR(3)*QR(3) + QI(3)*QI(3) )	    + ( QR(4)*QR(4) + QI(4)*QI(4) ) )

	WQ4=WNB/WDB

	RETURN
	End function
 

!double precision 
Function FACT(N)
	INTEGER I, N
        real*8 FACT
	FACT = 1.d0
	IF(N.LT.2)RETURN
	DO 100 I = 2, N
	   FACT = FACT*I
 100	CONTINUE
	RETURN
End Function


Subroutine fq6(N,xij,yij,zij,dij,fq6x,fq6y,fq6z)
          USE T_kind_param_m, ONLY:  double
          implicit none
!	implicit double precision(a-h,o-z)

        real(double)  OUT,gmaxmax,gminmin,rcut,delta_r,apol,bpol,cpol
        integer NEMAX,nwmax,N,nsize,i,k,j

	parameter(OUT=1D40,NEMAX=2500,nwmax=50)
	parameter(gmaxmax=1D10,gminmin=1.D-50)
        parameter(rcut=1.24, delta_r=0.1,apol=-3.D0/16.D0/delta_r**5)
       !parameter(rcut=1.392,delta_r=0.1,apol=-3.D0/16.D0/delta_r**5)
       parameter(bpol=5.D0/8.D0/delta_r**3,cpol=-15.D0/16.D0/delta_r)

	real(double) x(N),y(N),z(N)
	real(double) fwx(N),fwy(N),fwz(N)
	real(double) fQ0x(N),fQ0y(N),fQ0z(N)
	real(double) fQRx(6,N),fQRy(6,N),fQRz(6,N)
	real(double) fQIx(6,N),fQIy(6,N),fQIz(6,N)
	real(double) fq6x(N),fq6y(N),fq6z(N)
	real(double) COSPHI(6), SINPHI(6)
	double precision q6coef(4,0:6)
	double precision qr(6),qi(6)
	common/q6coef/q6coef
	real(double) xij(N,N)
	real(double) yij(N,N)
	real(double) zij(N,N)
	real(double) dij(N,N)
        real(double) Q0,DX,DY,DZ,WDB,dist,weight,dweight,w2
        real(double) costh,costh2,costh4,sinth,sinth2,scth,twocph
        real(double) dctdx,dctdy,dctdz,dn
        real(double) dcpdx,dcpdy,dcpdz
        real(double) dspdx,dspdy,dspdz
        real(double) wq0,dwq0,wq0_1,wq1,dwq1,wq1_1,wq2,dwq2,wq2_1,wq3,dwq3,wq3_1
        real(double) dfrx,dfry,dfrz
        real(double) dc2pdx,dc2pdy,dc2pdz
        real(double) dfix,dfiy,dfiz
        real(double) ds2pdx,ds2pdy,ds2pdz
        real(double) dcp,dcp2
        real(double) dc3pdx,dc3pdy,dc3pdz
        real(double) dc4pdx,dc4pdy,dc4pdz
        real(double) ds4pdx,ds4pdy,ds4pdz
        real(double) ds3pdx,ds3pdy,ds3pdz
        real(double) WQ6,WQ6_1,dWQ6,WNB,Q6,dc6pdx,dc6pdy,dc6pdz,ds6pdx,ds6pdy,ds6pdz
        real(double) SINTH4,DSTDCT,WQ4,WQ4_1,dWQ4
        real(double) WQ5,WQ5_1,dWQ5,dc5pdx,dc5pdy,dc5pdz,ds5pdx,ds5pdy,ds5pdz


	do i=1,N
	   fq6x(i)=0.D0
	   fq6y(i)=0.D0
	   fq6z(i)=0.D0
	   fwx(i)=0.D0
	   fwy(i)=0.D0
	   fwz(i)=0.D0
	   fq0x(i)=0.D0
	   fq0y(i)=0.D0
	   fq0z(i)=0.D0
	   do k=1,6
	      fqRx(k,i)=0.D0
	      fqRy(k,i)=0.D0
	      fqRz(k,i)=0.D0
	      fqIx(k,i)=0.D0
	      fqIy(k,i)=0.D0
	      fqIz(k,i)=0.D0
	   enddo
	enddo

	Q0=0.D0
	DO I=1,6
	   QI(I)=0.D0
	   QR(I)=0.D0
	ENDDO

	WDB=0.D0
	DO 30 J=1,N
	   DO 40 K=J+1,N
	      DX=xij(j,k)
	      DY=yij(j,k)
	      DZ=zij(j,k)
	      DIST=dij(j,k)

	      if (dist.le.rcut-delta_r) then
		 weight=1.D0
		 dweight=0.D0
	      elseif (dist.ge.rcut+delta_r) then
		 weight=0.D0
		 dweight=0.D0
	      else
		 weight=0.5D0+apol*(dist-rcut)**5 + bpol*(dist-rcut)**3 + cpol*(dist-rcut)
		 dweight=5.D0*apol*(dist-rcut)**4 + 3.D0*bpol*(dist-rcut)**2 + cpol
	      endif
	      WDB=WDB+WEIGHT

	      w2=-dweight/dist
	      fwx(j)=fwx(j)+dx*w2
	      fwy(j)=fwy(j)+dy*w2
	      fwz(j)=fwz(j)+dz*w2
	      fwx(k)=fwx(k)-dx*w2
	      fwy(k)=fwy(k)-dy*w2
	      fwz(k)=fwz(k)-dz*w2

	      COSTH  = DZ/DIST
	      COSTH2 = COSTH  * COSTH
	      COSTH4 = COSTH2 * COSTH2
	      SINTH2 = 1. - COSTH2
	      SINTH4 = SINTH2 * SINTH2
	      SINTH  = SQRT(SINTH2)
	      SCTH   = SINTH  * COSTH
	      IF (SINTH .EQ. 0.) THEN
		 COSPHI(1) = 1.
		 SINPHI(1) = 1.
	      ELSE
		 COSPHI(1) = DX/DIST/SINTH
		 SINPHI(1) = DY/DIST/SINTH
	      ENDIF
	      TWOCPH    = 2.*COSPHI(1)
	      COSPHI(2) = TWOCPH*COSPHI(1)-1.
	      SINPHI(2) = TWOCPH*SINPHI(1)
	      COSPHI(3) = TWOCPH*COSPHI(2)-COSPHI(1)
	      SINPHI(3) = TWOCPH*SINPHI(2)-SINPHI(1)
	      COSPHI(4) = TWOCPH*COSPHI(3)-COSPHI(2)
	      SINPHI(4) = TWOCPH*SINPHI(3)-SINPHI(2)
	      COSPHI(5) = TWOCPH*COSPHI(4)-COSPHI(3)
	      SINPHI(5) = TWOCPH*SINPHI(4)-SINPHI(3)
	      COSPHI(6) = TWOCPH*COSPHI(5)-COSPHI(4)
	      SINPHI(6) = TWOCPH*SINPHI(5)-SINPHI(4)

!       derivatives of cos theta and cos/sin phi

	      dctdx=dx*dz/dist**3
	      dctdy=dy*dz/dist**3
	      dctdz=-(1.D0-(dz/dist)**2)/dist

	      dn=dx**2+dy**2
	      dcpdx=-SINPHI(1)*dy/dn
	      dcpdy= COSPHI(1)*dy/dn
	      dcpdz= 0.D0
	      dspdx= dcpdy
	      dspdy=-COSPHI(1)*dx/dn
	      dspdz= 0.D0

	      DSTDCT=0.D0
	      IF (SINTH.NE.0.) DSTDCT=-COSTH/SINTH

!       accumulates Q0, QR, QI and their gradients

!       Q0

	WQ0=(Q6COEF(1,0)*COSTH2+Q6COEF(2,0))*COSTH4 + Q6COEF(3,0)*COSTH2+Q6COEF(4,0)
	      Q0=Q0+WQ0*WEIGHT

	      WQ0_1=w2*WQ0
	      dWQ0=2.D0*COSTH*(3.D0*Q6COEF(1,0)*COSTH4 + 2.D0*Q6COEF(2,0)*COSTH2+Q6COEF(3,0))

	      fq0x(j)=fq0x(j)+dx*WQ0_1+WEIGHT*dWQ0*dctdx
	      fq0y(j)=fq0y(j)+dy*WQ0_1+WEIGHT*dWQ0*dctdy
	      fq0z(j)=fq0z(j)+dz*WQ0_1+WEIGHT*dWQ0*dctdz
	      fq0x(k)=fq0x(k)-dx*WQ0_1-WEIGHT*dWQ0*dctdx
	      fq0y(k)=fq0y(k)-dy*WQ0_1-WEIGHT*dWQ0*dctdy
	      fq0z(k)=fq0z(k)-dz*WQ0_1-WEIGHT*dWQ0*dctdz

!       QR(1) and QI(1)

	      WQ1=(Q6COEF(1,1)*COSTH4+Q6COEF(2,1)*COSTH2 + Q6COEF(3,1))*SCTH
	      QR(1)=QR(1)+COSPHI(1)*WQ1*WEIGHT
	      QI(1)=QI(1)+SINPHI(1)*WQ1*WEIGHT

	      WQ1_1=w2*WQ1
	      dWQ1=(5.D0*Q6COEF(1,1)*COSTH4+3.D0*Q6COEF(2,1)*COSTH2  + Q6COEF(3,1))*SINTH  + (Q6COEF(1,1)*COSTH4+Q6COEF(2,1)*COSTH2  + Q6COEF(3,1))*COSTH*DSTDCT

	      dfrx=dWQ1*dctdx*COSPHI(1)+WQ1*dcpdx
	      dfry=dWQ1*dctdy*COSPHI(1)+WQ1*dcpdy
	      dfrz=dWQ1*dctdz*COSPHI(1)+WQ1*dcpdz

	      fqrx(1,j)=fqrx(1,j)+dx*COSPHI(1)*WQ1_1+WEIGHT*dfrx
	      fqry(1,j)=fqry(1,j)+dy*COSPHI(1)*WQ1_1+WEIGHT*dfry
	      fqrz(1,j)=fqrz(1,j)+dz*COSPHI(1)*WQ1_1+WEIGHT*dfrz
	      fqrx(1,k)=fqrx(1,k)-dx*COSPHI(1)*WQ1_1-WEIGHT*dfrx
	      fqry(1,k)=fqry(1,k)-dy*COSPHI(1)*WQ1_1-WEIGHT*dfry
	      fqrz(1,k)=fqrz(1,k)-dz*COSPHI(1)*WQ1_1-WEIGHT*dfrz

	      dfix=dWQ1*dctdx*SINPHI(1)+WQ1*dspdx
	      dfiy=dWQ1*dctdy*SINPHI(1)+WQ1*dspdy
	      dfiz=dWQ1*dctdz*SINPHI(1)+WQ1*dspdz

	      fqix(1,j)=fqix(1,j)+dx*SINPHI(1)*WQ1_1+WEIGHT*dfix
	      fqiy(1,j)=fqiy(1,j)+dy*SINPHI(1)*WQ1_1+WEIGHT*dfiy
	      fqiz(1,j)=fqiz(1,j)+dz*SINPHI(1)*WQ1_1+WEIGHT*dfiz
	      fqix(1,k)=fqix(1,k)-dx*SINPHI(1)*WQ1_1-WEIGHT*dfix
	      fqiy(1,k)=fqiy(1,k)-dy*SINPHI(1)*WQ1_1-WEIGHT*dfiy
	      fqiz(1,k)=fqiz(1,k)-dz*SINPHI(1)*WQ1_1-WEIGHT*dfiz

!       QR(2) and QI(2)

	      WQ2=(Q6COEF(1,2)*COSTH4 + Q6COEF(2,2)*COSTH2+Q6COEF(3,2))*SINTH2
	      QR(2)=QR(2)+COSPHI(2)*WQ2*WEIGHT
	      QI(2)=QI(2)+SINPHI(2)*WQ2*WEIGHT

	      WQ2_1=w2*WQ2
	      dWQ2=(-6.D0*Q6COEF(1,2)*COSTH4+4.D0*(Q6COEF(1,2) - Q6COEF(2,2))*COSTH2+2.D0*(Q6COEF(2,2)-Q6COEF(3,2)))*COSTH

	      dc2pdx=4.D0*COSPHI(1)*dcpdx
	      dc2pdy=4.D0*COSPHI(1)*dcpdy
	      dc2pdz=4.D0*COSPHI(1)*dcpdz

	      dfrx=dWQ2*dctdx*COSPHI(2)+WQ2*dc2pdx
	      dfry=dWQ2*dctdy*COSPHI(2)+WQ2*dc2pdy
	      dfrz=dWQ2*dctdz*COSPHI(2)+WQ2*dc2pdz

	      fqrx(2,j)=fqrx(2,j)+dx*COSPHI(2)*WQ2_1+WEIGHT*dfrx
	      fqry(2,j)=fqry(2,j)+dy*COSPHI(2)*WQ2_1+WEIGHT*dfry
	      fqrz(2,j)=fqrz(2,j)+dz*COSPHI(2)*WQ2_1+WEIGHT*dfrz
	      fqrx(2,k)=fqrx(2,k)-dx*COSPHI(2)*WQ2_1-WEIGHT*dfrx
	      fqry(2,k)=fqry(2,k)-dy*COSPHI(2)*WQ2_1-WEIGHT*dfry
	      fqrz(2,k)=fqrz(2,k)-dz*COSPHI(2)*WQ2_1-WEIGHT*dfrz

	      ds2pdx=2.D0*(COSPHI(1)*dspdx+SINPHI(1)*dcpdx)
	      ds2pdy=2.D0*(COSPHI(1)*dspdy+SINPHI(1)*dcpdy)
	      ds2pdz=2.D0*(COSPHI(1)*dspdz+SINPHI(1)*dcpdz)

	      dfix=dWQ2*dctdx*SINPHI(2)+WQ2*ds2pdx
	      dfiy=dWQ2*dctdy*SINPHI(2)+WQ2*ds2pdy
	      dfiz=dWQ2*dctdz*SINPHI(2)+WQ2*ds2pdz

	      fqix(2,j)=fqix(2,j)+dx*SINPHI(2)*WQ2_1+WEIGHT*dfix
	      fqiy(2,j)=fqiy(2,j)+dy*SINPHI(2)*WQ2_1+WEIGHT*dfiy
	      fqiz(2,j)=fqiz(2,j)+dz*SINPHI(2)*WQ2_1+WEIGHT*dfiz
	      fqix(2,k)=fqix(2,k)-dx*SINPHI(2)*WQ2_1-WEIGHT*dfix
	      fqiy(2,k)=fqiy(2,k)-dy*SINPHI(2)*WQ2_1-WEIGHT*dfiy
	      fqiz(2,k)=fqiz(2,k)-dz*SINPHI(2)*WQ2_1-WEIGHT*dfiz

!       QR(3) and QI(3)

	      WQ3=(Q6COEF(1,3)*COSTH2+Q6COEF(2,3))*SCTH*SINTH2
	      QR(3)=QR(3)+COSPHI(3)*WQ3*WEIGHT
	      QI(3)=QI(3)+SINPHI(3)*WQ3*WEIGHT

	      WQ3_1=w2*WQ3
	      dWQ3=(-5.D0*Q6COEF(1,3)*COSTH4+3.D0*(Q6COEF(1,3) - Q6COEF(2,3))*COSTH2+Q6COEF(2,3))*SINTH  + (Q6COEF(1,3)*COSTH2+Q6COEF(2,3))*COSTH*SINTH2*DSTDCT

	      dcp=4.D0*COSPHI(1)**2-1.D0
	      dc3pdx=3.D0*dcpdx*dcp
	      dc3pdy=3.D0*dcpdy*dcp
	      dc3pdz=3.D0*dcpdz*dcp

	      dfrx=dWQ3*dctdx*COSPHI(3)+WQ3*dc3pdx
	      dfry=dWQ3*dctdy*COSPHI(3)+WQ3*dc3pdy
	      dfrz=dWQ3*dctdz*COSPHI(3)+WQ3*dc3pdz

	      fqrx(3,j)=fqrx(3,j)+dx*COSPHI(3)*WQ3_1+WEIGHT*dfrx
	      fqry(3,j)=fqry(3,j)+dy*COSPHI(3)*WQ3_1+WEIGHT*dfry
	      fqrz(3,j)=fqrz(3,j)+dz*COSPHI(3)*WQ3_1+WEIGHT*dfrz
	      fqrx(3,k)=fqrx(3,k)-dx*COSPHI(3)*WQ3_1-WEIGHT*dfrx
	      fqry(3,k)=fqry(3,k)-dy*COSPHI(3)*WQ3_1-WEIGHT*dfry
	      fqrz(3,k)=fqrz(3,k)-dz*COSPHI(3)*WQ3_1-WEIGHT*dfrz

	      ds3pdx=dcp*dspdx+8.D0*SINPHI(1)*COSPHI(1)*dcpdx
	      ds3pdy=dcp*dspdy+8.D0*SINPHI(1)*COSPHI(1)*dcpdy
	      ds3pdz=dcp*dspdz+8.D0*SINPHI(1)*COSPHI(1)*dcpdz

	      dfix=dWQ3*dctdx*SINPHI(3)+WQ3*ds3pdx
	      dfiy=dWQ3*dctdy*SINPHI(3)+WQ3*ds3pdy
	      dfiz=dWQ3*dctdz*SINPHI(3)+WQ3*ds3pdz

	      fqix(3,j)=fqix(3,j)+dx*SINPHI(3)*WQ3_1+WEIGHT*dfix
	      fqiy(3,j)=fqiy(3,j)+dy*SINPHI(3)*WQ3_1+WEIGHT*dfiy
	      fqiz(3,j)=fqiz(3,j)+dz*SINPHI(3)*WQ3_1+WEIGHT*dfiz
	      fqix(3,k)=fqix(3,k)-dx*SINPHI(3)*WQ3_1-WEIGHT*dfix
	      fqiy(3,k)=fqiy(3,k)-dy*SINPHI(3)*WQ3_1-WEIGHT*dfiy
	      fqiz(3,k)=fqiz(3,k)-dz*SINPHI(3)*WQ3_1-WEIGHT*dfiz

!       QR(4) and QI(4)

	      WQ4=(Q6COEF(1,4)*COSTH2+Q6COEF(2,4))*SINTH4
	      QR(4)=QR(4)+COSPHI(4)*WQ4*WEIGHT
	      QI(4)=QI(4)+SINPHI(4)*WQ4*WEIGHT

	      WQ4_1=w2*WQ4
	      dWQ4=(6.D0*Q6COEF(1,4)*COSTH4+4.D0*(Q6COEF(2,4) - 2.D0*Q6COEF(1,4))*COSTH2 + 2.D0*(Q6COEF(1,4)-2.D0*Q6COEF(2,4)))*COSTH

	      dcp=4.D0*COSPHI(1)*(2.D0*COSPHI(1)**2-1.D0)
	      dc4pdx=4.D0*dcpdx*dcp
	      dc4pdy=4.D0*dcpdy*dcp
	      dc4pdz=4.D0*dcpdz*dcp

	      dfrx=dWQ4*dctdx*COSPHI(4)+WQ4*dc4pdx
	      dfry=dWQ4*dctdy*COSPHI(4)+WQ4*dc4pdy
	      dfrz=dWQ4*dctdz*COSPHI(4)+WQ4*dc4pdz

	      fqrx(4,j)=fqrx(4,j)+dx*COSPHI(4)*WQ4_1+WEIGHT*dfrx
	      fqry(4,j)=fqry(4,j)+dy*COSPHI(4)*WQ4_1+WEIGHT*dfry
	      fqrz(4,j)=fqrz(4,j)+dz*COSPHI(4)*WQ4_1+WEIGHT*dfrz
	      fqrx(4,k)=fqrx(4,k)-dx*COSPHI(4)*WQ4_1-WEIGHT*dfrx
	      fqry(4,k)=fqry(4,k)-dy*COSPHI(4)*WQ4_1-WEIGHT*dfry
	      fqrz(4,k)=fqrz(4,k)-dz*COSPHI(4)*WQ4_1-WEIGHT*dfrz

	      dcp2=4.D0*SINPHI(1)*(6.D0*COSPHI(1)**2-1.D0)
	      ds4pdx=dspdx*dcp+dcpdx*dcp2
	      ds4pdy=dspdy*dcp+dcpdy*dcp2
	      ds4pdz=dspdz*dcp+dcpdz*dcp2

	      dfix=dWQ4*dctdx*SINPHI(4)+WQ4*ds4pdx
	      dfiy=dWQ4*dctdy*SINPHI(4)+WQ4*ds4pdy
	      dfiz=dWQ4*dctdz*SINPHI(4)+WQ4*ds4pdz

	      fqix(4,j)=fqix(4,j)+dx*SINPHI(4)*WQ4_1+WEIGHT*dfix
	      fqiy(4,j)=fqiy(4,j)+dy*SINPHI(4)*WQ4_1+WEIGHT*dfiy
	      fqiz(4,j)=fqiz(4,j)+dz*SINPHI(4)*WQ4_1+WEIGHT*dfiz
	      fqix(4,k)=fqix(4,k)-dx*SINPHI(4)*WQ4_1-WEIGHT*dfix
	      fqiy(4,k)=fqiy(4,k)-dy*SINPHI(4)*WQ4_1-WEIGHT*dfiy
	      fqiz(4,k)=fqiz(4,k)-dz*SINPHI(4)*WQ4_1-WEIGHT*dfiz

!       QR(5) and QI(5)

	      WQ5=Q6COEF(1,5)*SINTH4*SCTH
	      QR(5)=QR(5)+COSPHI(5)*WQ5*WEIGHT
	      QI(5)=QI(5)+SINPHI(5)*WQ5*WEIGHT

	      WQ5_1=w2*WQ5
	      dWQ5=Q6COEF(1,5)*(SINTH*(1.D0-6.D0*COSTH2+5.D0*COSTH4) + SINTH4*COSTH*DSTDCT)
	      dcp=16.D0*COSPHI(1)**4-12.D0*COSPHI(1)**2+1.D0
	      dc5pdx=5.D0*dcpdx*dcp
	      dc5pdy=5.D0*dcpdy*dcp
	      dc5pdz=5.D0*dcpdz*dcp

	      dfrx=dWQ5*dctdx*COSPHI(5)+WQ5*dc5pdx
	      dfry=dWQ5*dctdy*COSPHI(5)+WQ5*dc5pdy
	      dfrz=dWQ5*dctdz*COSPHI(5)+WQ5*dc5pdz

	      fqrx(5,j)=fqrx(5,j)+dx*COSPHI(5)*WQ5_1+WEIGHT*dfrx
	      fqry(5,j)=fqry(5,j)+dy*COSPHI(5)*WQ5_1+WEIGHT*dfry
	      fqrz(5,j)=fqrz(5,j)+dz*COSPHI(5)*WQ5_1+WEIGHT*dfrz
	      fqrx(5,k)=fqrx(5,k)-dx*COSPHI(5)*WQ5_1-WEIGHT*dfrx
	      fqry(5,k)=fqry(5,k)-dy*COSPHI(5)*WQ5_1-WEIGHT*dfry
	      fqrz(5,k)=fqrz(5,k)-dz*COSPHI(5)*WQ5_1-WEIGHT*dfrz

	      dcp2=4.D0*SINPHI(2)*(8.D0*COSPHI(1)**2-3.D0)
	      ds5pdx=dspdx*dcp+dcpdx*dcp2
	      ds5pdy=dspdy*dcp+dcpdy*dcp2
	      ds5pdz=dspdz*dcp+dcpdz*dcp2

	      dfix=dWQ5*dctdx*SINPHI(5)+WQ5*ds5pdx
	      dfiy=dWQ5*dctdy*SINPHI(5)+WQ5*ds5pdy
	      dfiz=dWQ5*dctdz*SINPHI(5)+WQ5*ds5pdz

	      fqix(5,j)=fqix(5,j)+dx*SINPHI(5)*WQ5_1+WEIGHT*dfix
	      fqiy(5,j)=fqiy(5,j)+dy*SINPHI(5)*WQ5_1+WEIGHT*dfiy
	      fqiz(5,j)=fqiz(5,j)+dz*SINPHI(5)*WQ5_1+WEIGHT*dfiz
	      fqix(5,k)=fqix(5,k)-dx*SINPHI(5)*WQ5_1-WEIGHT*dfix
	      fqiy(5,k)=fqiy(5,k)-dy*SINPHI(5)*WQ5_1-WEIGHT*dfiy
	      fqiz(5,k)=fqiz(5,k)-dz*SINPHI(5)*WQ5_1-WEIGHT*dfiz

!       QR(6) and QI(6)

	      WQ6=Q6COEF(1,6)*SINTH4*SINTH2
	      QR(6)=QR(6)+COSPHI(6)*WQ6*WEIGHT
	      QI(6)=QI(6)+SINPHI(6)*WQ6*WEIGHT

	      WQ6_1=w2*WQ6
	      dWQ6=Q6COEF(1,6)*(-6.D0*COSTH4+12.D0*COSTH2-6.D0)*COSTH

	      dcp=COSPHI(1)*(32.D0*COSPHI(1)**4-32.D0*COSPHI(1)**2+6.D0)
	      dc6pdx=6.D0*dcpdx*dcp
	      dc6pdy=6.D0*dcpdy*dcp
	      dc6pdz=6.D0*dcpdz*dcp

	      dfrx=dWQ6*dctdx*COSPHI(6)+WQ6*dc6pdx
	      dfry=dWQ6*dctdy*COSPHI(6)+WQ6*dc6pdy
	      dfrz=dWQ6*dctdz*COSPHI(6)+WQ6*dc6pdz

	      fqrx(6,j)=fqrx(6,j)+dx*COSPHI(6)*WQ6_1+WEIGHT*dfrx
	      fqry(6,j)=fqry(6,j)+dy*COSPHI(6)*WQ6_1+WEIGHT*dfry
	      fqrz(6,j)=fqrz(6,j)+dz*COSPHI(6)*WQ6_1+WEIGHT*dfrz
	      fqrx(6,k)=fqrx(6,k)-dx*COSPHI(6)*WQ6_1-WEIGHT*dfrx
	      fqry(6,k)=fqry(6,k)-dy*COSPHI(6)*WQ6_1-WEIGHT*dfry
	      fqrz(6,k)=fqrz(6,k)-dz*COSPHI(6)*WQ6_1-WEIGHT*dfrz

	      dcp2=SINPHI(1)*(160.D0*COSPHI(1)**4-96.D0*COSPHI(1)**2+6.D0)
	      ds6pdx=dspdx*dcp+dcpdx*dcp2
	      ds6pdy=dspdy*dcp+dcpdy*dcp2
	      ds6pdz=dspdz*dcp+dcpdz*dcp2

	      dfix=dWQ6*dctdx*SINPHI(6)+WQ6*ds6pdx
	      dfiy=dWQ6*dctdy*SINPHI(6)+WQ6*ds6pdy
	      dfiz=dWQ6*dctdz*SINPHI(6)+WQ6*ds6pdz

	      fqix(6,j)=fqix(6,j)+dx*SINPHI(6)*WQ6_1+WEIGHT*dfix
	      fqiy(6,j)=fqiy(6,j)+dy*SINPHI(6)*WQ6_1+WEIGHT*dfiy
	      fqiz(6,j)=fqiz(6,j)+dz*SINPHI(6)*WQ6_1+WEIGHT*dfiz
	      fqix(6,k)=fqix(6,k)-dx*SINPHI(6)*WQ6_1-WEIGHT*dfix
	      fqiy(6,k)=fqiy(6,k)-dy*SINPHI(6)*WQ6_1-WEIGHT*dfiy
	      fqiz(6,k)=fqiz(6,k)-dz*SINPHI(6)*WQ6_1-WEIGHT*dfiz

 40	   CONTINUE
 30	CONTINUE

	WNB = DSQRT(    Q0 * Q0  + ( QR(1)*QR(1) + QI(1)*QI(1) ) + ( QR(2)*QR(2) + QI(2)*QI(2)) + (QR(3)*QR(3) + QI(3)*QI(3) ) + ( QR(4)*QR(4) + QI(4)*QI(4) ) + ( QR(5)*QR(5) + QI(5)*QI(5) ) + ( QR(6)*QR(6) + QI(6)*QI(6) ) )

	Q6=WNB/WDB

	!print *,'Q6=',q6

	do i=1,n
	   fq6x(i)=Q0*fq0x(i)
	   fq6y(i)=Q0*fq0y(i)
	   fq6z(i)=Q0*fq0z(i)
	   do k=1,6
	      fq6x(i)=fq6x(i)+QR(k)*fqrx(k,i)+QI(k)*fqix(k,i)
	      fq6y(i)=fq6y(i)+QR(k)*fqry(k,i)+QI(k)*fqiy(k,i)
	      fq6z(i)=fq6z(i)+QR(k)*fqrz(k,i)+QI(k)*fqiz(k,i)
	   enddo
	   fq6x(i)=fq6x(i)/WDB/WNB-WNB*fwx(i)/WDB**2
	   fq6y(i)=fq6y(i)/WDB/WNB-WNB*fwy(i)/WDB**2
	   fq6z(i)=fq6z(i)/WDB/WNB-WNB*fwz(i)/WDB**2
	enddo

	!print *,'fq6x ',fq6x

	return
End Subroutine


Function q6(N,xij,yij,zij,dij)
    USE T_kind_param_m, ONLY:  double
    implicit none
        real(double)  OUT,gmaxmax,gminmin,rcut,delta_r,apol,bpol,cpol
        integer NEMAX,nwmax,N,nsize,i,k,j
	parameter(OUT=1d40,NEMAX=2500,nwmax=50)
	parameter(gmaxmax=1e10,gminmin=1.D-50)
       !parameter(rcut=1.392,delta_r=0.1,apol=-3.D0/16.D0/delta_r**5)
        parameter(rcut=1.24, delta_r=0.1,apol=-3.D0/16.D0/delta_r**5)
        parameter(bpol=5.D0/8.D0/delta_r**3,cpol=-15.D0/16.D0/delta_r)

	real(double) x(N),y(N),z(N)
	real(double) q6coef(4,0:6)
	real(double) qr(6),qi(6)
	common/q6coef/q6coef
        real(double) WDB,WNB,DX,DY,DZ,DIST,weight
	real(double) xij(N,N)
	real(double) yij(N,N)
	real(double) zij(N,N)
	real(double) dij(N,N)
!	common/pairdist/xij,yij,zij,dij
        real(double) WQ6,Q6,Q0
	Q0=0.D0
	DO I=1,6
	   QI(I)=0.D0
	   QR(I)=0.D0
	ENDDO

	WDB=0.D0
	DO 30 J=1,N
	   DO 40 K=J+1,N
	      DX=xij(j,k)
	      DY=yij(j,k)
	      DZ=zij(j,k)
	      DIST=dij(j,k)
	      if (dist.le.rcut-delta_r) then
		 weight=1.D0
	      elseif (dist.ge.rcut+delta_r) then
		 weight=0.D0
	      else
		 weight=0.5D0+apol*(dist-rcut)**5 +bpol*(dist-rcut)**3 + cpol*(dist-rcut)
	      endif
  	      WDB=WDB+WEIGHT
	      CALL EVASH6(DX/DIST,DY/DIST,DZ/DIST, Q0, QR, QI, WEIGHT)
 40	   CONTINUE
 30	CONTINUE

!cccccccccccccccccccccccccc
!c	Q0=0.D0
!c	DO I=5,5
!c	   QI(I)=0.D0
!c	   QR(I)=0.D0
!c	ENDDO
!cccccccccccccccccccccccccc

	Q6 = WQ6(WDB,WNB,Q0,QR,QI)

	return
End Function



Subroutine EVASH6(X, Y, Z, Q0, QR, QI, WEIGHT)
 
!       *** Evaluate the spherical harmonics of degree 6 ********************
!C       
!C       Notebook JvD A38 and A85
!C       
!C       *********************************************************************
	
	COMMON/Q6COEF/ Q6COEF
	double precision   Q6COEF(4,0:6)
	
	double precision  X, Y, Z,WEIGHT
	DOUBLE PRECISION  Q0, QR(6), QI(6)
	
	double precision COSTH, COSTH2, COSTH4, SINTH, SINTH2, SCTH
	double precision TWOCPH, COSPHI(6), SINPHI(6), TEMP, SINTH4
	
	COSTH  = Z
	COSTH2 = COSTH  * COSTH
	COSTH4 = COSTH2 * COSTH2
	SINTH2 = 1. - COSTH2
	SINTH4 = SINTH2 * SINTH2
	SINTH  = SQRT(SINTH2)
	SCTH   = SINTH  * COSTH
	
!       *** Is THETA = 0 ? Then PHI is irrelevant *****************
	
	IF (SINTH .EQ. 0.) THEN
	   COSPHI(1) = 1.
	   SINPHI(1) = 1.
	ELSE
	   COSPHI(1) = X/SINTH
	   SINPHI(1) = Y/SINTH
	ENDIF
	TWOCPH    = 2.*COSPHI(1)
	COSPHI(2) = TWOCPH*COSPHI(1)-1.
	SINPHI(2) = TWOCPH*SINPHI(1)
	COSPHI(3) = TWOCPH*COSPHI(2)-COSPHI(1)
	SINPHI(3) = TWOCPH*SINPHI(2)-SINPHI(1)
	COSPHI(4) = TWOCPH*COSPHI(3)-COSPHI(2)
	SINPHI(4) = TWOCPH*SINPHI(3)-SINPHI(2)
	COSPHI(5) = TWOCPH*COSPHI(4)-COSPHI(3)
	SINPHI(5) = TWOCPH*SINPHI(4)-SINPHI(3)
	COSPHI(6) = TWOCPH*COSPHI(5)-COSPHI(4)
	SINPHI(6) = TWOCPH*SINPHI(5)-SINPHI(4)

!       *** Now the spherical harmonics are calculated and summed with Q
!       
!       This part of the subroutine would have been easier to understand
!       when complex numbers would have been used. However, COMPLEX*16
!       variables are not included in the FORTRAN 77 standard.
	
	TEMP=(Q6COEF(1,0)*COSTH2+Q6COEF(2,0))*COSTH4 + Q6COEF(3,0)*COSTH2+Q6COEF(4,0)
	Q0  = Q0 + WEIGHT*TEMP

	TEMP=(Q6COEF(1,1)*COSTH4+Q6COEF(2,1)*COSTH2 + Q6COEF(3,1))*SCTH
      QR(1) = QR(1) + COSPHI(1)*TEMP*WEIGHT
      QI(1) = QI(1) + SINPHI(1)*TEMP*WEIGHT

      TEMP  = (Q6COEF(1,2)*COSTH4 + Q6COEF(2,2)*COSTH2+Q6COEF(3,2))*SINTH2
      QR(2) = QR(2) + COSPHI(2)*TEMP*WEIGHT
      QI(2) = QI(2) + SINPHI(2)*TEMP*WEIGHT

      TEMP  = (Q6COEF(1,3)*COSTH2+Q6COEF(2,3))*SCTH*SINTH2
      QR(3) = QR(3) + COSPHI(3)*TEMP*WEIGHT
      QI(3) = QI(3) + SINPHI(3)*TEMP*WEIGHT

      TEMP  = (Q6COEF(1,4)*COSTH2+Q6COEF(2,4))*SINTH4
      QR(4) = QR(4) + COSPHI(4)*TEMP*WEIGHT
      QI(4) = QI(4) + SINPHI(4)*TEMP*WEIGHT

      TEMP  = Q6COEF(1,5)*SINTH4*SCTH
      QR(5) = QR(5) + COSPHI(5)*TEMP*WEIGHT
      QI(5) = QI(5) + SINPHI(5)*TEMP*WEIGHT

      TEMP  = Q6COEF(1,6)*SINTH4*SINTH2
      QR(6) = QR(6) + COSPHI(6)*TEMP*WEIGHT
      QI(6) = QI(6) + SINPHI(6)*TEMP*WEIGHT


	RETURN
end Subroutine

!	double precision 

Function WQ6(WDB,WNB,Q0,QR,QI)



       USE T_kind_param_m, ONLY:  double
!C       *** The order parameter Q6 is calculated from Q0 and Q *************
!C       
!C       INTEGER  NB      Number of bonds accumulated
!C       DOUBLE PRECISION
!C       Q0      Accumulator of P60(cos(THETA))
!C       QR(6)   Accumulator of P6m(cos(THETA))*cos(PHI)
!C       QI(6)   Accumulator of P6m(cos(THETA))*sin(PHI)!
!C       
!C       In the expression for Q4, the factor pi*13/4 disappears
!C       The summation has to be performed over the squares of Q4m, with
!C       m running from -6 to +6. However, Q6m and Q6-m are equal when
!C       squared. So the summation is done for positive m and a factor
!C       2 is introduced.
!C       
!C       ********************************************************************
	
	INTEGER   NB
	real(double)  Q0, QR(6), QI(6),WNB,WDB,WQ6

	WNB = SQRT(    Q0 * Q0 + ( QR(1)*QR(1) + QI(1)*QI(1) ) + ( QR(2)*QR(2) + QI(2)*QI(2) ) + ( QR(3)*QR(3) + QI(3)*QI(3) ) + ( QR(4)*QR(4) + QI(4)*QI(4) ) + ( QR(5)*QR(5) + QI(5)*QI(5) )  + ( QR(6)*QR(6) + QI(6)*QI(6) ) )

	WQ6=WNB/WDB

	RETURN

End function

