module sundae_module
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use var_pot
  use jqmod
  use random_art
  use lanczos_defs
  !-----------------------------------------------
  ! test
  !   _c -> courant
  !   _s -> selection
  !   _d -> depart

      integer,parameter :: dpkind=selected_real_kind(13)
      real(double) :: t0,t1,elaps_1,tbuffer1,tbuffer2,t_lanczos,t_propag
      integer,     dimension(:),allocatable,save  :: ipovois

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

      integer                                       :: it_tele_vac
      integer                                       :: lupout 
      integer                                       :: luwout
      integer                                       :: lufout 
      integer                                       :: lusout
      integer                                       :: luvout
      integer                                       :: luhout
      integer                                       :: ncomptinter
      integer                                       :: idistance
      integer :: iteration
      integer :: iloop
      integer :: icheck
!      integer :: nchemin
      integer :: nfrequence,im1!,it_art

      real(double),dimension(:),allocatable,save    :: tmass_tele_vac
      real(double) ::  tstep_tele_vac,usdh_tele_vac, convert_tele_vac
      real(double) :: sigma
      real(double) :: xdist0
      real(double) :: xdist1
      real(double) :: kappas2 , potistadd
      real(double) :: barbetaq
      real(double) :: barbeta,barbetaE
      real(double) :: betaweff
      real(double) :: Ecinetique
      real(double) :: Ecinetique0
      real(double) :: dlambda
      real(double) :: lambdax
      real(double) :: lambday
      real(double) :: lambdaz
      real(double) :: deltawork
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
      logical      :: ldeter
      logical      :: ldistance
      logical      :: lta
!     parametre du lennard-jones
      character*80 :: fnamfin
      character*80 :: fnamhin

      integer N
      PARAMETER(N=127)   !1023  ! 127
      integer NHZ
      PARAMETER(NHZ=60)
      integer nfenetre,lufin,luhin
      parameter(nfenetre=100)
      real(double) :: rsig_lj      
      real(double) :: sig_lj    
      real(double) :: eps_4lj
      real(double) :: xp4
      real(double) :: pot_auxi(-10:nfenetre+10)
      real(double) :: pot_auxiliary
      real(double) :: contour_auxiliary
      real(double) :: kappaE,kappaEs2
      real(double) :: potistaddE
      real(double) :: alphadd(2)
      real(double) :: pot_contour(-10:nfenetre+10,0:nfenetre)
      real(double) :: energie_min,energie_max
      real(double) :: ergtoev
      integer :: nl_iter
      parameter(ergtoev=6.2415d11)

! cosmin added:
      integer    :: it_langevin_deter=0, it_langevin=0,it_trajectory     
contains

subroutine allocate_tele_vac ! ****
  implicit none
!    allocate ( ielat_c(imm),iwmax_c(imm), ityp_c(imm),xp_c(3,imm),  xpp_c(3,imm), &
!               vp_c(3,imm), ax_c(3,imm), fp_c(3,imm),                             &
!               ielat_s(imm),iwmax_s(imm), ityp_s(imm),xp_s(3,imm),  xpp_s(3,imm), &
!               vp_s(3,imm), ax_s(3,imm), fp_s(3,imm),                             &
!  	       ielat_d(imm),iwmax_d(imm), ityp_d(imm),xp_d(3,imm),  xpp_d(3,imm), &
!               vp_d(3,imm), ax_d(3,imm), fp_d(3,imm),                             &
    allocate ( m_i(3,imm),gau(6,imm),sig_i(3,imm),rga_i(3,imm),fadd(3,imm),       &
               xbarini(3))

  allocate ( tmass_tele_vac(imm) )
  allocate (ipovois(15), d2vois(15),xpvois(3,15),xpvoisini(3,15),xtransla(3))

  end subroutine allocate_tele_vac

  subroutine init_tele_vac (xp, vp, fp, ielat, iwmax, ityp)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    implicit none
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer       :: ielat(imm)
    integer       :: iwmax(imm)
    integer       :: ityp(imm)
    real(double)  :: xp(3,imm)
!    real(double)  :: xpp(3,imm)
    real(double)  :: vp(3,imm)
!   real(double)  :: ax(3,imm)
    real(double)  :: fp(3,imm)
    real(double)  :: gamma, text_teledyn

   !-----------------------------------------------
     integer   :: ic_local,iatom
   !-----------------------------------------------

   text_teledyn=0.19
   it_tele_vac=0
   tstep_tele_vac=tstep / utemps
   convert_tele_vac=9.6485d0/10000.d0
   im1=im+1
   do ic_local=1,im
      iatom = ic_local
      if (iatom==0) iatom=im     
      tmass_tele_vac(ic_local) = tstep_tele_vac**2*convert_tele_vac/(cm(ityp(iatom))/umass)
   enddo
   
   usdh_tele_vac= 1.d0/(2.d0*tstep_tele_vac)

   sigma=1.d-2 
   write(6,*) ' kappa ',kappa 
   kappas2 = kappa/two

   ! initialisation des parametres du langevin
   gamma=one/(tstep*1d2)
   m_i(1:3,1:im) = cm(1)

   rga_i(:,:) = exp(-gamma*tstep/two)
   sig_i(:,:) = sqrt(m_i(:,:)*text_teledyn*bk*(one-exp(-gamma*tstep)))

   write(6,*) ' rga_i  ', rga_i(:,1:1)
   write(6,*) ' sig_i ', sig_i(:,1:1)
   write(6,*) ' bk text_teledyn',bk,text_teledyn,one-exp(-gamma*tstep/two)
!   write(6,*) ' cm ',m_i(1,1:im)
!   stop

!   parametres additionnels

   m_tot = sum(m_i(1,1:im))

   fnampout = fnam(1:lenfnam)//'.pout'
   fnamwout = fnam(1:lenfnam)//'.wout'
   fnamfout = fnam(1:lenfnam)//'.fout'
   fnamsout = fnam(1:lenfnam)//'.sout'
   fnamvout = fnam(1:lenfnam)//'.vout'

   write(6,*) ' fnampout',fnampout
   lupout  = 17
   luwout  = 18
   lufout  = 19
   lusout  = 20
   luvout  = 21

   write(6,*) ' xbarini =',xbarini, im
   
!   x111        = sum(xp(1:3,1))
!  ldeter    = .true.
   return

 end subroutine init_tele_vac


 subroutine langevin(dt,temperature,rga)!,ielat,iwmax, ityp)!)!(xp, vp, fp,dt)!, ielat,iwmax, ityp)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double
    use gen_com_m
    use tab_imm_m

    implicit none
    integer ic
    real(double) :: pp(3,im)
    real(double), dimension (3,N)::sig
    real(double) :: xbar(3)
    real(double) :: Ecin0, Ecin1, Ecin3, Ecin4
    real(double) :: dt
    real(double) :: vbar(3)
    real(double), dimension(6,N+1)::gau
    real(double) ::temperature,rga
    sig(:,:) = sqrt(m_i(:,:)*temperature*(one-rga**2))
    call genere_bruit2(sig,gau)
    sig_i(:,:)=sig(:,:)
    rga_i(:,:)=rga

    Ecin0 = zero
    Ecin1 = zero
    Ecin3 = zero
    Ecin4 = zero

    im=N
    
    pp(1:3,1:im)=vp(1:3,1:im)*m_i(1:3,1:im)

    do ic=1,3
       vbar(ic)=sum(vp(ic,1:im))/dble(im) ! barycentre sur les particules
    enddo

	  if (itab/=0) then
	    if (mod(it_langevin,itab)==0) then
	       call caltabt
	    endif
	 endif
	 if (ltabvois.and.mod(it_langevin,itetabvois)==0) call caltabi
	 call calfo

     
     do ic =1,3
       Ecin0 = Ecin0   + DOT_PRODUCT(pp(ic,1:im),vp(ic,1:im))/two
     enddo
    pp(1:3,1:im)=pp(1:3,1:im)*rga_i(1:3,1:im) + gau(1:3,1:im)
    do ic=1,3
       vbar(ic)=sum(pp(ic,1:im))/dble(im) ! vit barycentre sur les particules
       pp(ic,1:im)=pp(ic,1:im)-vbar(ic)
    enddo
    !prot(1:3,1:im)=pp(1:3,1:im)
    call control_angular_momenta(pp,xp)
    !pp(1:3,1:im) = prot(1:3,1:im)
    vp(1:3,1:im)=pp(1:3,1:im)/m_i(1:3,1:im)


    do ic =1,3
       Ecin1 = Ecin1   + DOT_PRODUCT(pp(ic,1:im),vp(ic,1:im))/two
    enddo
    
    pp(1:3,1:im) =  pp(1:3,1:im) + (fp(1:3,1:im))*dt/two
    
    do ic=1,3
       xbar(ic)    = sum(xp(ic,1:im))/dble(im) ! barycentre sur les particules
    enddo
    xp(1:3,1:im)=  xp(1:3,1:im) + pp(1:3,1:im)*dt/m_i(1:3,1:im)
    

    do ic=1,3
       xbar(ic)    = sum(xp(ic,1:im))/dble(im)-xbar(ic) ! déplacement du barycentre 
       xp(ic,1:im) = xp(ic,1:im) - xbar(ic)  ! on recentre tout le systeme
    enddo

	  if (itab/=0) then
	    if (mod(it_langevin,itab)==0) then
	       call caltabt
	    endif
	 endif
	 if (ltabvois.and.mod(it_langevin,itetabvois)==0) call caltabi
	 call calfo
  
!debugC    write(*,*) 'Langevin.........:', itab, ltabvois, itetabvois, it_langevin
    
    pp(1:3,1:im) = pp(1:3,1:im) + (fp(1:3,1:im))*dt/two
    vp(1:3,1:im) = pp(1:3,1:im) / m_i(1:3,1:im)
    
    
    do ic=1,3
       Ecin3 = Ecin3 + DOT_PRODUCT(pp(ic,1:im),vp(ic,1:im))/two
    enddo

    pp(1:3,1:im) = pp(1:3,1:im)*rga_i(1:3,1:im) + gau(4:6,1:im)
    do ic=1,3
       vbar(ic)=sum(pp(ic,1:im))/dble(im) ! vit barycentre sur les particules
       pp(ic,1:im)=pp(ic,1:im)-vbar(ic)
    enddo
    !prot(1:3,1:im)=pp(1:3,1:im)
    call control_angular_momenta(pp,xp)!,prot)
    !pp(1:3,1:im) = prot(1:3,1:im)
!    call control_angular_momenta(pp,xp)
    vp(1:3,1:im) = pp(1:3,1:im)/m_i(1:3,1:im)

    do ic=1,3
      Ecin4 = Ecin4 + DOT_PRODUCT(pp(ic,1:im),vp(ic,1:im))/two
    enddo

    Ecinetique = Ecin4
   
    return

 end subroutine langevin


Subroutine LyapLanczos_vac! (xp)

  use lanczos_defs!, only: eigenvalue
 !USE T_kind_param_m, ONLY:  double
    !use jqmod
    use tab_imm_m
implicit none

integer, dimension (:), allocatable :: Nb
integer, dimension (:), allocatable :: Number
integer, parameter ::N=127 ! 127  !1023  numero atomes du cluster
integer :: jl 
real (double) ::eta    !paramètre d'énergie du potentiel
real (double) ::sigma ! dist. éq. pour le potentiel LJ 

real (double),dimension (:,:),allocatable::absdmax,xpq6,ener
integer:: kl 

real (double), dimension (3,N):: xppro
real (double), dimension(1:N):: xref
real (double), dimension(1:N):: yref
real (double), dimension(1:N):: zref
real (double), dimension(1:3,1:N):: qref,qrot,prot

real (double), dimension(3,1:im):: q,p
real (double), dimension(1:N,6):: qtemp 


real (double), dimension(1:3,1:N):: q1s2
real (double), dimension(1:3):: xbar

real (double) ::   xalea,xcumul

type Trajectoire
    real (double), dimension(1:3,1:N):: q  ! vector of position
    real (double), dimension(1:3,1:N):: p  ! vector of impulsion
    real(double), dimension(3*N) :: project
    real (double)  Lyap
    real (double)  eigenvalue
 endtype

type (Trajectoire), dimension (:,:), allocatable :: Path
type (Trajectoire), dimension (:,:), allocatable :: Pshoot
type (Trajectoire), dimension (:,:), allocatable :: Pshift

type (Trajectoire) :: Pcourant

real (double) ::mcconf,ranf, l0, alpha,tempo,teq,z, a_sto, p3
real (double) :: kinetotdt
real(double)::fftot, fftotdt,laputot,laputotdt,uproject
real (double) ::dt,TotalTime,gamma, Temperature ! Lyapunov=0,LyapunovAv=0,mu=0,
integer :: nbClones, j,l,k,a,i, nq4, ltot,mcmoves,totalmcmoves, newtraj, it_art
logical ::  new_projection
logical ::  waste_recycling
character(len=128) :: sortie
character(len=128) :: posfinal
character(len=128) :: data_mbar
character(len=128) :: dada_mbar
character(len=128) :: data_mbar_std
character(len=128) :: moyennes_mbar
character(len=128) :: moyennes_mbar_denom

character(len=128) :: kappaF
character(len=128) :: kappaFd
character(len=128) :: fnamtin

real (double)::rga,rien,kine,kinetot, pi
real (double) ::ss,delta_x, kapa, p1,p2
integer, parameter::nfenetre=100
integer ::  iter,kappa,ix, scrivi, continue,maxvec

real (double):: stat(0:nfenetre, 0:nfenetre)
real (double):: stat6(0:nfenetre, 0:nfenetre)
real(double)::  enprmoy

real(double) :: tab_contour(-10:nfenetre+10,0:nfenetre)
real(double) :: tab_contour_q6(-10:nfenetre+10,0:nfenetre)
real(double) :: tab_cont_q4q6(-10:nfenetre+10,-10:nfenetre+10)
real(double) :: cumul_contour(-10:nfenetre+10,0:nfenetre)
real(double) :: cumul_contour_q6(-10:nfenetre+10,0:nfenetre)
real(double) :: cumul_contour_q4q6(-10:nfenetre+10,-10:nfenetre+10)
real (double), dimension (:), allocatable::norm0

real (double):: tequilib, e,timefsh
real (double) :: ekin,epot
real (double) absdmax_current
integer  :: pix,lanczos_iter
integer:: iterfw, iterbw,nmax,totiter,   tprim, Nbclones_mbar, depart_boucle_nbclones, atom_bouge_abs, tprimo,lutin

real (double) ,dimension(:),   allocatable::poids,alpha_bias,dist,absdist
real (double) ,dimension(:),   allocatable:: ener0 
real (double) ,dimension(:),   allocatable:: enerK
real (double) ,dimension(:),   allocatable::triallyap
real (double) ,dimension(:),   allocatable::oldlyap 
real (double) ,dimension(:),   allocatable::rapport
real (double) ,dimension(:),   allocatable::h_A
real (double) ,dimension(:),   allocatable::tau
real (double) ,dimension(:),   allocatable:: hamilt
real (double) ,dimension(:,:), allocatable:: dh
real (double) ,dimension(:,:), allocatable:: dhx
real (double) ,dimension(:,:), allocatable:: Psel
real (double) ,dimension(:,:), allocatable:: S
real (double) ,dimension(:,:,:), allocatable:: u_kln
real (double) ,dimension(:,:,:), allocatable:: ustd_kln
real (double) ,dimension(:,:,:), allocatable:: umoy_kln
real (double) ,dimension(:,:,:), allocatable:: u2moy_kln

integer ,dimension(:), allocatable::acc

real (double) ,dimension(:), allocatable::h_F
real (double) ,dimension(:), allocatable::react_F
real (double) ,dimension(:), allocatable::h_Fd
real (double) ,dimension(:), allocatable::react_Fd
real (double) ,dimension(:), allocatable::h_dI
real (double) ,dimension(:), allocatable::react_dI
real (double) ,dimension(:), allocatable::h_Fg
real (double) ,dimension(:), allocatable::react_Fg
real (double) ,dimension(:), allocatable::h_FI
real (double) ,dimension(:), allocatable::react_FI
!#real (double) ,dimension(:), allocatable::react_FCC_col
real (double) :: tempiter,tempvar,eigenvalue_old
real(double), parameter :: KtoERG=1.3791946308724831d-16
real(double) ::  h_A_max,h_ba_max,h_ba_min,h_ba,h_ba_I,h_temp
real(double) :: pav(3)
real(double) :: genrand

namelist /input_sundae/Totalmcmoves,TotalTime,dt,NbClones,temperature,gamma,sortie,alpha, teq, delta_x, a_sto, kapa, continue, Nbclones_mbar, depart_boucle_nbclones,tprimo,maxvec,h_A_max,h_ba_min,h_ba_max,h_ba,h_ba_I

 h_A_max  = 7.d-1
 h_ba_max = 1.8d0
 h_ba_min = 0.8d0
 h_ba     = 1.3d0

 h_A_max  = 4.5d-1
 h_ba_max = 1.6d0
 h_ba_min = 1.0d0
 h_ba     = 1.3d0
 h_ba_I   = 1.8d0

fnamtin = fnam(1:lenfnam)//'.tin'
  ! variables de dynamique
  write(*,*) 'file name', fnamtin
  lutin = 777
  open(unit=lutin, file=fnamtin, status='unknown')
  read (lutin, nml=input_sundae)



write(6,*)'usage %s:\n'

!    Initialisation of the variables
write(6,*)'TotalMcmoves = ', Totalmcmoves ! Total duration of the simulation 

write(6,*)'No of steps = ',  TotalTime      ! Total duration of the simulation 

write(6,*)'dt = ', dt ! Time step
totiter = int(TotalTime)
TotalTime=real(TotalTime)*dt

write(6,*)'NbClones = le numero du canal : ' , NbClones    ! Number of the clone (channel)
write(6,*)'temperature = '  , temperature  ! temperature (kT)
!   
temperature=temperature*KtoERG
!
write(6,*)'friction*dt = '  , gamma  ! friction*dt
! interval between 2 data records
write(6,*)'sortie data_mbar = ',sortie  ! Name of the file where data are stored
write(6,*)'alpha max = ', alpha 
write(6,*)'t_equilib', teq
write(6,*)'delta_X = ', delta_x
a_sto = 1.d0-2.d0*((1.d1**(-2.d0-2.d0*dble(nbclones)/dble(nbclones_mbar))))
write(6,*)'a_sto = ', a_sto
write(6,*)'k ressort = ', kapa
write(6,*)'continue = ', continue
write(6,*)'nbclones_mbar = ', Nbclones_mbar
write(6,*)'cluster? 0 no, nbclones yes'
write (*,*) ' depart_boucle_nbclones',depart_boucle_nbclones
write(6,*)'tprimo = ', tprimo
write(*,*) 'maxvec for Lanczos = ', maxvec

ss=sqrt(temperature*1.66*1e-24*55.845)
write(*,*) 'ss  = ', ss

lenfnam   = index(sortie,' ')-1
!open (unit=11, FILE=sortie, action='write', status='replace')
!open (unit=31, action='write', status='replace')
moyennes_mbar  =sortie(1:lenfnam)//'.data_moy'
moyennes_mbar_denom  =sortie(1:lenfnam)//'.data_moy2'
data_mbar      = sortie(1:lenfnam)//'.data'
dada_mbar      = sortie(1:lenfnam)//'.dada'
data_mbar_std  = sortie(1:lenfnam)//'.data_std'

!histotot  =sortie(1:lenfnam)//'.thout'
!histopart =sortie(1:lenfnam)// '.phout'
posfinal = sortie(1:lenfnam)//'.cin'

kappaF = sortie(1:lenfnam)//'.corfunc_WR'
kappaFd = sortie(1:lenfnam)//'.corfunc_ST'
!kappadI = sortie(1:lenfnam)//'.kappadI'
!kappaFI= sortie(1:lenfnam)//'.kappaFI'
!kappaFg= sortie(1:lenfnam)//'.kappaFg'
!corbias= sortie(1:lenfnam)//'.corbias'

!!!!!!!!!! output per calcolo correl function reactivity
open(unit=30, file=kappaF,  action='write', status='replace')
open(unit=31, file=kappaFd,  action='write', status='replace')
!open(unit=32, file=kappadI,  action='write', status='replace')
!open(unit=33, file=kappaFI,  action='write', status='replace')
!open(unit=34, file=kappaFg,  action='write', status='replace')
!open(unit=333,file=corbias,  action='write', status='replace')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!open(unit=27, file=data_mbar, action='write', status='replace')
open(unit=24, file=data_mbar, action='write', status='replace')
open(unit=244, file=dada_mbar, action='write', status='replace')
!open(unit=44, file=data_mbar_std, action='write', status='replace')

write(244,'(i,i,i)') Nbclones_mbar+1,Nbclones_mbar+1,totalmcmoves

!open(unit=25, file=moyennes_mbar, action='write', status='replace')
!open(unit=23, file=moyennes_mbar_denom, action='write', status='replace')
!nmax=NbClones*3.0+1    ! Max number of clones
!TimeStore = TotalTime-AverageTime ! Next time at which data will be stored 

nmax=NbClones_mbar

write(*,*) ' nmax ', nmax
 
! Array where the weight of every clone is stored
allocate (acc(0:nmax))
allocate (tau(0:nmax))
!allocate (betaq(0:nmax))
allocate (ener0(0:nmax))
allocate (enerK(0:nmax))
!allocate(enerpro(0:nmax-1))
allocate(hamilt(0:nmax))
! Array which contains the numbers from 1 to nmax
allocate(Nb(0:nmax)) 
allocate(norm0(0:nmax)) 

kappa=1
scrivi=0
ltot=int(totaltime/real(kappa*dt))
!totiter=nint(totaltime/dt)

write(*,*)'totiter', totiter
gamma=gamma/(dt)
write(*,*)'gamma' , gamma, gamma*dt,ltot
waste_recycling=.true.
!waste_recycling=.false.
write(*,*) 'waste_recycling?',  waste_recycling
allocate (h_A(0:2*totiter+10))
allocate (absdmax(0:nmax,0:2*totiter+10))
allocate (xpq6(0:nmax,0:2*totiter+10))
allocate(ener(0:nmax,0:totiter+10))
allocate(S(0:nmax,0:2*totiter+10))
allocate(u_kln(0:nmax,0:nmax,0:Totalmcmoves))
allocate(ustd_kln(0:nmax,0:nmax,0:Totalmcmoves))
allocate(umoy_kln(0:nmax,0:nmax,0:Totalmcmoves))
allocate(u2moy_kln(0:nmax,0:nmax,0:Totalmcmoves))

allocate(old_projection(3*N))
!allocate(projection(3*N))
allocate(first_projection(3*N))
allocate(old_before_sc_projection(3*N))
!allocate(project(0:totiter+10,3*N))

allocate (h_F(0:2*totiter+10))
allocate (react_F(0:nmax))
allocate (h_Fd(0:2*totiter+10))
allocate (react_Fd(0:nmax))
allocate (h_Fg(0:2*totiter+10))
allocate (react_Fg(0:nmax))
allocate (h_dI(0:2*totiter+10))
allocate (react_dI(0:nmax))
allocate (h_FI(0:2*totiter+10))
allocate (react_FI(0:nmax))

!#allocate (react_FCC_col(0:nbclones_MBAR))

ener(:,:)=0
stat(:,:)=0
stat6(:,:)=0
sigma=1
eta=1
!enerpro(:)=0
absdmax(:,:)=0
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
it_art=0
tau(:)=0
enprmoy=0
nq4=0
e=0.000001
eigenvalue=0
!ss=0.05
l0=0.0
acc(:)=0
p1=0
p2=0
p3=0

  if (maxvec.lt.4) maxvec=4
  write(*,*) 'The maxvec in lanczos is set to .........:', maxvec

     tab_contour(-10:nfenetre+10,0:nfenetre)=zero
     tab_contour_q6(-10:nfenetre+10,0:nfenetre)=zero
     tab_cont_q4q6(-10:nfenetre+10,-10:nfenetre+10)=zero
     cumul_contour_q4q6(-10:nfenetre+10,-10:nfenetre+10)=0
     cumul_contour(-10:nfenetre+10,0:nfenetre)=0
     cumul_contour_q6(-10:nfenetre+10,0:nfenetre)=0

allocate(Number(0:NbClones-1)) 

!!!!!!!!!!!!T EQUILIBRAGE
tequilib=(dt)*teq
write(*,*) 'No of teq steps', teq
write(*,*) 'tquilib',tequilib
!!!!!!!!!!!! T EQUILIBRAGE
  
write(*,*) 'sortie affichee', 0.6*ev2erg

allocate(Path(0:nmax,0:totiter+10)) 
allocate(Pshoot(0:nmax,0:totiter+10))
allocate(Pshift(0:nmax,0:2*totiter+10))

allocate (dh(0:nmax,0:2*totiter+10))
allocate(oldLyap(0:nmax))
allocate(triallyap(0:nmax))
allocate(rapport(0:nmax))
allocate(dist(0:N))
allocate(absdist(0:N))
allocate (dhx(0:nmax,0:2*totiter+10))
allocate(poids(1:totiter))
allocate(Psel(0:nbclones_mbar,0:totiter))
allocate (alpha_bias(0:Nbclones_mbar))!!EQUILIBRAGE
xppro(:,:)=0
rapport(:)=0
iter=0
rien=0
Q=0
uproject=0
kappa=1
pi=4*atan(1._dpkind)
if ((gamma*dt).le.100000) then
    rga = exp(-gamma*dt/two)
else 
  rga=0
endif
sig(:,:) = sqrt(temperature*(one-rga**2))
write(*,*) 'rga',gamma, gamma*dt, rga, sig(1,1)
write(*,*) 'temp', temperature*erg2eV, erg2eV

dh(:,:)=0
acc(:)=0
!! initialisation paramètres de bias alfa pour reconstruction
do j= 0, Nbclones_mbar
 alpha_bias(j)=1.d12*(real(j*(alpha/real(Nbclones_mbar))))
 write(*,*)'bias', alpha_bias(j)  ! ,real(j*(alpha/real(Nbclones_mbar)))
enddo

!write(*,*)'bias', alpha_bias(:),real(j*(alpha/real(nbclones_mbar)))


xref(1:n)=xp(1,1:n)
yref(1:n)=xp(2,1:n)
zref(1:n)=xp(3,1:n)

qref(1:3,1:im)=xp(1:3,1:im)

 rang = 1 
tempo=0.0

 if (continue.ne.2) then 
   j = depart_boucle_nbclones
   absdmax(j,:)=-9999.0
   atom_bouge_abs=0
   iter=0
   do i=1,N
     absdist(i)=sqrt((xp(1,i)-qref(1,i))**2+(xp(2,i)-qref(2,i))**2+(xp(3,i)-qref(3,i))**2)
      if (absdist(i).gt.absdmax(j,iter)) then
       absdmax(j,iter) = absdist(i)
       atom_bouge_abs  = i
      endif
!   write(*,*) ' absdist(i) ' ,   absdist(i)*1.d8
   enddo
   absdmax_current = absdmax(j,iter)*1.d8
   write(*,*) ' absdmax(1,0)   ',absdmax(j,iter)*1.d8
   write(*,*) ' atom_bouge_abs ',atom_bouge_abs

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! INITIALISATION !!!!!!!!!!!!!!!!!!!!!!
  write(*,*) 'conditions initiales pour traj de reference avc distrib stoch à temperature T=',temperature
  write(*,*) 'temps equilibrage', Tequilib, tempo
  !!!!!!!!!!!!!!!!!!!!!!!! EQUILIBRAGE initial (STOCH DYN)!!!!!!!!!!!!!!!!!!
  write (*,'("Equilibrage on ..................:",i5)') nint(Tequilib/dt) 
  
   do it_langevin = 0,int(Tequilib/dt)
    call langevin(dt,temperature, rga)
    !debugC	write(*,*) 'Lyaplanc.................:', it_langevin

    absdmax(j,:)=-9999.0
    atom_bouge_abs=0
    iter=0
    do i=1,N
      absdist(i)=sqrt((xp(1,i)-qref(1,i))**2+(xp(2,i)-qref(2,i))**2+(xp(3,i)-qref(3,i))**2)
      if (absdist(i).gt.absdmax(j,iter)) then
       absdmax(j,iter) = absdist(i)
       atom_bouge_abs  = i
      endif
    !write(*,*) ' absdist(i) ' ,   absdist(i)*1.d8
    enddo
    absdmax_current = absdmax(j,iter)*1.d8
    write(*,*) ' absdmax(1,0)   ',absdmax(j,iter)*1.d8
    write(*,*) ' atom_bouge_abs ',atom_bouge_abs
   end do

   absdmax(j,:)=-9999.0
   atom_bouge_abs=0
   iter=0
   do i=1,N
     absdist(i)=sqrt((xp(1,i)-qref(1,i))**2+(xp(2,i)-qref(2,i))**2+(xp(3,i)-qref(3,i))**2)
      if (absdist(i).gt.absdmax(j,iter)) then
       absdmax(j,iter) = absdist(i)
       atom_bouge_abs  = i
      endif
!   write(*,*) ' absdist(i) ' ,   absdist(i)*1.d8
   enddo
   absdmax_current = absdmax(j,iter)*1.d8
   write(*,*) ' absdmax(1,0)   ',absdmax(j,iter)*1.d8
   write(*,*) ' atom_bouge_abs ',atom_bouge_abs


  write(*,*) 'point de depart traj de reference deterministe' 
  ! NbClones the index number of bias
  ! The startingpoint in the bias series

   do  jl = 0,NbClones 
     Path(jl,0)%q=xp(1:3,1:N)
     Path(jl,0)%p=vp(1:3,1:N)*m_i(1:3,1:N)
   enddo
  !!!!!!!!!!!!!! mettre a zero le reste avant ??
  Path(:,:)%Lyap = 0.d0

  iter=0
  it_trajectory=0
  new_projection=.true.
  it_art=1
  lanczos_iter=0
  eigenvalue = 0.0

  do icheck=1,10

  Pcourant=Path(j,iter)

  call mapping_P_Verlet(Pcourant%q,Pcourant%p,dt,N,q1s2)

  call lanczos(N,maxvec,q1s2,new_projection,Path(j,iter)%project)!!!positions avant propagation

  Path(j,iter)%eigenvalue = eigenvalue
  write(*,*) 'eigenvalue ',eigenvalue
  new_projection=.false.  
   if (eigenvalue.lt.0.0) then 
                   Path(j,iter)%Lyap =  asin(dt*sqrt(-eigenvalue)/2.d0)*2.d0 
   else 
                   Path(j,iter)%Lyap = 0.d0  ! 
   endif
  enddo

 endif  ! continue.ne.2 

  write(*,*) 'initialisation faite: go with Lanczos'
  write(*,*) 'itab', itab
  write(*,*) 'itetabvois', itetabvois
  write(*,*) 'ltabvois', ltabvois

if (continue.eq.2) then
   write(*,*) 'posfinal ',posfinal
  open(unit=27, file=posfinal, status='old')
  do j=depart_boucle_nbclones,NbClones
   write(*,*) 'j depart_boucle_nbclones NbClones  ', j ,depart_boucle_nbclones,NbClones
   do i=0, totiter-1                      ! attention i=totiter doit etre pris en compte
    read (27,*) qtemp(1:N,1)
    read (27,*) qtemp(1:N,2)
    read (27,*) qtemp(1:N,3)
    read (27,*) 
    read (27,*) qtemp(1:N,4)
    read (27,*) qtemp(1:N,5)
    read (27,*) qtemp(1:N,6)
    read (27,*)
    read (27,*) Path(j,i)%Lyap
    read (27,*)
    read (27,*) Path(j,i)%project

    Path(j,i)%q(1,1:N)=qtemp(1:N,1)
    Path(j,i)%q(2,1:N)=qtemp(1:N,2)
    Path(j,i)%q(3,1:N)=qtemp(1:N,3)
    Path(j,i)%p(1,1:N)=qtemp(1:N,4)
    Path(j,i)%p(2,1:N)=qtemp(1:N,5)
    Path(j,i)%p(3,1:N)=qtemp(1:N,6)
    !Path(j,i)%p
   enddo
   close(27)
   Path(j,totiter) =  Path(j,totiter-1) 
   oldLyap(j)=sum(Path(j,0:totiter-1)%Lyap)/real(totiter)
!   oldLyap(j)=sum(Path(j,:)%Lyap)/real(totiter)
   write(*,*) 'oldLyap(',j,') = ', oldLyap(j)


   q=Path(j,0)%q
   do i=1,3
   write(*,*) 'i       ',i
   write(*,*) 'bary    ',sum(q(i,1:im))/dble(im)
   write(*,*) 'bary ref',sum(qref(i,1:im))/dble(im)
   enddo


   rang=0
   do i=0,totiter
   !write(*,*) 'i  ',i
   
   q    = Path(j,i)%q
   prot = Path(j,i)%p

   pav(1)=sum(prot(1,1:im))/dble(im)
   pav(2)=sum(prot(2,1:im))/dble(im)
   pav(3)=sum(prot(3,1:im))/dble(im)

!   write(*,*) 'prot(:,1):', prot

!   write(*,*) 'pav',pav

   prot(1,1:im)= prot(1,1:im)-pav(1)
   prot(2,1:im)= prot(2,1:im)-pav(2)
   prot(3,1:im)= prot(3,1:im)-pav(3)

   pav(1)=sum(prot(1,1:im))/dble(im)
   pav(2)=sum(prot(2,1:im))/dble(im)
   pav(3)=sum(prot(3,1:im))/dble(im)
!   write(*,*) 'pav',pav

   call control_angular_momenta(prot,q)

!   write(*,*) 'q en 1 ',q(1,1)
!   write(*,*) 'somme de q(1,:) ',sum(q(1,1:im))
  
   qrot=q-qref
   call control_angular_momenta(qrot,qref)
   call control_angular_momenta(qrot,qref)
   q = qrot + qref
   rang=1
   call control_angular_momenta(qrot,qref)
   rang=1

   Path(j,i)%q=q

   call control_angular_momenta(prot,q)

   Path(j,i)%p=prot

   enddo

   absdmax(j,:)=-9999.0
   atom_bouge_abs=0
   iter=0
   do i=1,N
     pav(1:3)=(Path(j,iter)%q(1:3,i)-qref(1:3,i))**2
     !write(*,*) 'pav ',pav
     absdist(i)=sqrt(sum(pav(1:3)))
      if (absdist(i).gt.absdmax(j,iter)) then
       absdmax(j,iter) = absdist(i)
       atom_bouge_abs  = i
      endif
!   write(*,*) ' absdist(i) ' ,   absdist(i)*1.d8
   enddo
   absdmax_current = absdmax(j,iter)*1.d8
   write(*,*) ' absdmax(1,0)   ',absdmax(j,iter)*1.d8
   write(*,*) ' atom_bouge_abs ',atom_bouge_abs

  enddo

 endif ! block with continue == 2
! test du moment angulaire de la force
 call caltabt
 call caltabi
 call calfo

 p(1:3,1:im) = fp(1:3,1:im)
 q (1:3,1:im) = xp(1:3,1:im)
 rang=0
 call control_angular_momenta(p,q)
 call control_angular_momenta(p,q)
 call control_angular_momenta(p,q)
 ! stop

! test de la subroutine de controle angulaire
 j=depart_boucle_nbclones

 q= Path(j,0)%q
 p= Path(j,0)%p

 temp=0.d0
 do i=1,20
 call OU_control(p,q,a_sto,ss)
 call mapping_P_Verlet(q,p,dt,N,q1s2)

 temp = temp + (sum(p(1,1:im)**2)+sum(p(2,1:im)**2)+sum(p(3,1:im)**2))/dble(3*im-6)/m_i(1,1)
 write(*,*) 'température ',temp/dble(i)/KtoERG,i,alpha_bias(j)
 enddo
 !stop
 write (*,'("1st traj before shooting ........:",i5)') nint(Totaltime/dt)
 rang=1
 j=depart_boucle_nbclones

 call cpu_time(t0)
 lanczos_iter=0
 eigenvalue_old = eigenvalue
 do icheck=1,1
   Pcourant = Path(j,0)
   temp=0.d0
   eigenvalue = Pcourant%eigenvalue
   do iter=0,totiter-1
     call mapping_P_Verlet(Pcourant%q,Pcourant%p,dt,N,q1s2)
     call cpu_time(tbuffer1)
     lanczos_iter=lanczos_iter+nl_iter
     call lanczos(N,maxvec,q1s2,new_projection,Pcourant%project)!!!
     call cpu_time(tbuffer2)
     t_lanczos=t_lanczos + (tbuffer2-tbuffer1)
     Pcourant%Lyap = 0.d0  
     if (eigenvalue.lt.0.d0) then
           Pcourant%Lyap = asin(dt*sqrt(-eigenvalue)/2.d0)*2.d0
     endif
     Path(j,iter+1)    = Pcourant
     Path(j,iter)%Lyap = Pcourant%Lyap
     temp = temp +  Pcourant%Lyap
     new_projection = .false.
   enddo  !!!!!
   oldLyap(j)=sum(Path(j,0:totiter-1)%Lyap)/real(totiter)
   write(*,*) " oldLyap(j)= ", oldLyap(j),temp/real(totiter)
 enddo

  pav(1)=sum(Path(j,totiter)%p(1,1:im))/dble(im)
  pav(2)=sum(Path(j,totiter)%p(2,1:im))/dble(im)
  pav(3)=sum(Path(j,totiter)%p(3,1:im))/dble(im)

  write(*,*) 'pav',pav

 call cpu_time(t1)
 elaps_1=t1-t0

 write(*,*) 'The first trajectory..:', elaps_1
 write(*,*) 'The Lanczos time .....:', t_lanczos
 write(*,*) 'The Propag time.......:', elaps_1-t_lanczos
 write(*,*) 'Lanczos interations...:', lanczos_iter
 write(*,*) '            forces....:', lanczos_iter*maxvec*2

!   write(*,*) 'pos fin du shooting', Path(j,totiter)%q(0), Path(j,totiter)%q(0), Path(j,totiter)%q(0)
!   write(*,*) 'mom fin du shooting', Path(j,totiter)%p(0), Path(j,totiter)%p(0), Path(j,totiter)%p(0)
!   write(*,*) 'pos fin du shooting', Path(j,0)%q(10), Path(j,0)%q(10), Path(j,0)%q(10)
!   write(*,*) 'mom fin du shooting', Path(j,0)%p(10), Path(j,0)%p(10), Path(j,0)%p(10)
 absdmax(j,:)=-9999.d0
  do iter=0,totiter-1
    atom_bouge_abs=0
    do i=1,N
      pav(1:3) = (Path(j,iter)%q(1:3,i)-qref(1:3,i))**2
      absdist(i)=sqrt(sum(pav(1:3)))
       if (absdist(i).gt.absdmax(j,iter)) then
        absdmax(j,iter) = absdist(i)
        atom_bouge_abs  = i
       endif
    enddo
 enddo

  temp=0.0000
  tempvar=temp
  j= depart_boucle_nbclones

  do iterbw=0,totiter
  tempiter = SUM((Path(j,iterbw)%p(1,1:N)**2+Path(j,iterbw)%p(2,1:N)**2+Path(j,iterbw)%p(3,1:N)**2))/m_i(1,1)/3.000/dble(im-2)/KtoERG ! *2/3 /2
   !debugC write(*,*) 'tempiter ',iterbw,tempiter
   temp=temp+tempiter
   tempvar=tempvar+tempiter**2
  enddo
  temp=temp/real(totiter+1)
  tempvar=sqrt(tempvar/real(totiter+1)-temp**2)
  write(*,*) 'temp cin traject, std et cible ',temp,tempvar,temperature/KtoERG
! calcul de la position initiale 

  call cal_hamilton(Path(j,0)%q,Path(j,0)%p,N,ekin,epot)
  hamilt(j)=(ekin+epot)/temperature
  ener0(j) = epot/temperature

! vérification de la dérive 
 xbar(1) = SUM(Path(j,0)%q(1,1:N))/dble(N)
 xbar(2) = SUM(Path(j,0)%q(2,1:N))/dble(N)
 xbar(3) = SUM(Path(j,0)%q(3,1:N))/dble(N)

 write(*,*) ' centre de masse référence ', xbar(1:3)

 xbar(1) = SUM(qref(1,1:N))/dble(N)
 xbar(2) = SUM(qref(2,1:N))/dble(N)
 xbar(3) = SUM(qref(3,1:N))/dble(N)

 write(*,*) ' centre de masse pos 0     ', xbar(1:3)

 xbar(1) = SUM(Path(j,0)%p(1,1:N))/dble(N)
 xbar(2) = SUM(Path(j,0)%p(2,1:N))/dble(N)
 xbar(3) = SUM(Path(j,0)%p(3,1:N))/dble(N)

 write(*,*) ' moments translationnels ',   xbar(1:3)
 do i=1,3
  write(*,*) ' temperature translation ',  xbar(i)**2/m_i(1,1)/dble(im)/KtoERG*dble(N)
 enddo
 iterbw = 0
 tempiter = SUM(Path(j,iterbw)%p(1,1:N)**2+Path(j,iterbw)%p(2,1:N)**2+Path(j,iterbw)%p(3,1:N)**2)/m_i(1,1)/3.000/dble(im)/KtoERG ! *2/3 /2
 write(*,*) 'temp avant ',iterbw,tempiter
 enerK(j) =tempiter
 do iterfw=0,totiter
 Path(j,iterfw)%p(1,1:N)=Path(j,iterfw)%p(1,1:N)-xbar(1)
 Path(j,iterfw)%p(2,1:N)=Path(j,iterfw)%p(2,1:N)-xbar(2)
 Path(j,iterfw)%p(3,1:N)=Path(j,iterfw)%p(3,1:N)-xbar(3)
 enddo

 tempiter = SUM(Path(j,iterbw)%p(1,1:N)**2+Path(j,iterbw)%p(2,1:N)**2+Path(j,iterbw)%p(3,1:N)**2)/m_i(1,1)/3.000/dble(im)/KtoERG ! *2/3 /2
 write(*,*) 'tempiter après ',iterbw,tempiter
 write(*,*) 'diff ',tempiter-enerK(j) 
! stop



 do mcmoves = 1 , Totalmcmoves                       !mcmoves = 1,M in the paper
  it_art=mcmoves

 ix=int((totiter+1)*genrand())
 write(*,*) 
 write(*,*) 
 write(*,'(" Itération..:",i5," sur un total de ",i5)') ,mcmoves,Totalmcmoves
 write(*,'(" Shooting......................:",3i5)') mcmoves,totiter,ix

 !!!! shooting deterministique
 ! Check: if depart_boucle_nbclones,NbClones of if there are different WHY?

!  write(*,*) 'depart_boucle_nbclones NbClones  ', depart_boucle_nbclones,NbClones
!  write(*,*)'curr bias', alpha_bias(j)

 do j=depart_boucle_nbclones,NbClones

!!!!! initialization backward

! verif manuel : Massimiliano avait oublie de mettre a jour la variable oldlyap
!  oldLyapunov=SUM(Path(j,0:totiter-1)%Lyap)/real(totiter)
!  write(*,*) ' verif oldLyap',oldLyapunov,oldLyap(j)

! modif manuel

    Pshoot(j,ix) = Path(j,ix)

    eigenvalue = Path(j,ix)%eigenvalue

    temp = SUM((Path(j,ix)%p(1,1:N)**2+Path(j,ix)%p(2,1:N)**2+Path(j,ix)%p(3,1:N)**2))/m_i(1,1)/3.000/dble(im)/KtoERG ! *2/3 /2
    write(*,*) 'temp cin et temperature cible',temp,temperature/KtoERG

   !!!! shooting  time "a la Stoltz"

    p(1:3,1:im) = Pshoot(j,ix)%p(1:3,1:im)
    q(1:3,1:im) = Pshoot(j,ix)%q(1:3,1:im)

    call OU_control(p,q,a_sto,ss)  ! d'amplitude a_sto 

    Pshoot(j,ix)%p(1:3,1:N)= p(1:3,1:N)

    temp = SUM((Pshoot(j,ix)%p(1,1:N)**2+Pshoot(j,ix)%p(2,1:N)**2+Pshoot(j,ix)%p(3,1:N)**2))/m_i(1,1)/dble(3*im-6)/KtoERG ! *2/3 /2
    write(*,*) 'temp cin apres perturb et temp cible ',temp,temperature/KtoERG
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! propagation BACKWARD a partir du point de shooting (sans lyapunov)

    write(*,*) 'shooting ..ix->0...........bw ',ix,' -> ', 0
    new_projection=.false.
    it_trajectory=0
    Pcourant =  Pshoot(j,ix)
!    write(*,*) 'ix   ',ix
    if (ix.gt.0) then 
     Pcourant%project     = Path(j,ix)%project
     Pcourant%eigenvalue  = Path(j,ix)%eigenvalue
    else
     Pcourant%project     = Path(j,0)%project
     Pcourant%eigenvalue  = Path(j,0)%eigenvalue
    endif
!   question de Manuel pour Cosmin : on stocke en iterbw-1 la valeur propre calculée à iterbw-1/2
    lanczos_iter=0
    do iterbw=ix,1,-1 ! on stocke  en iterbw-1 la valeur propre calculée à iterbw-1/2
    !
     call mapping_P_Verlet(Pcourant%q,Pcourant%p,-dt,N,q1s2) !! on met -dt pour le backward !!!
!    modif Manuel
     call lanczos(N,maxvec,q1s2,new_projection,Pcourant%project)
     lanczos_iter=lanczos_iter+nl_iter
     Pcourant%eigenvalue = eigenvalue 
      if (eigenvalue.lt.0.0) then 
!    modif Manuel
           Pcourant%Lyap=asin(dt*sqrt(-eigenvalue)/2.d0)*2.d0 !
       else 
           Pcourant%Lyap=0.d0
      endif
!    modif Manuel
     Pshoot(j,iterbw-1)=Pcourant
     it_trajectory = it_trajectory + 1
     new_projection = .false.
     enddo  
!!!!!!!!!!!!!!!!!!!!!  la trajectoire bw est arrivee au temps t=0
!!!!!!!!!!!!           trajectoire forward        

    write(*,*) 'shooting ..ix->totiter-1...fw ',ix,' -> ', totiter-1
    Pcourant       = Pshoot(j,ix)
    new_projection =.false. ! Question Manuel: on doit recalculer?
    it_trajectory=0
    eigenvalue = Pshoot(j,ix)%eigenvalue

! Manuel : avec mapping_P_Verlet on calcule les forces en iterfw+1/2 et on stocke positions et moments en iterfw+1 
!                 mais on stocke le lyapunov avant en iterfw 
    do iterfw = ix,totiter-1  ! on recalcule en ix car la position en ix+1/2 est modifiée
!    modif Manuel (texte+commentaires) 
     call mapping_P_Verlet(Pcourant%q,Pcourant%p,dt,N,q1s2)
     call lanczos(N,maxvec,q1s2,new_projection,Pcourant%project)
     lanczos_iter=lanczos_iter+nl_iter
     Pcourant%eigenvalue = eigenvalue
     if (eigenvalue.le.0.0) then 
      Pcourant%Lyap = asin(dt*sqrt(-eigenvalue)/2.d0)*2.d0  ! on stocke avant !!!
     else 
      Pcourant%Lyap = 0.d0      ! on stocke toujours le lyapunov avant !!!
     endif
     Pshoot(j,iterfw)%Project=Pcourant%Project
     Pshoot(j,iterfw)%eigenvalue=Pcourant%eigenvalue
     Pshoot(j,iterfw)%Lyap=Pcourant%Lyap
     Pshoot(j,iterfw+1)   = Pcourant
     it_trajectory  = it_trajectory + 1
     new_projection = .false.
    enddo
      write(*,*) 'Lanczos iterations...:',lanczos_iter
      write(*,*) '           forces....:',lanczos_iter*maxvec*2
! 
    triallyap(j)=SUM(Pshoot(j,0:totiter-1)%Lyap)/real(totiter)
    write(*,*) 'oldlyap(j)   = ', oldLyap(j)
    write(*,*) 'triallyap(j) = ', triallyap(j)
    !!!!!!! fonction heaviside pour maintenir le point initial de la trajectoire dans 
    !!!!!!! le bassin de depart: calcul dist max pour le point x(0) de la trajectoire
    absdmax(j,iterbw)=-9999.0
    atom_bouge_abs=0
     do i=1,N

      absdist(i)=sqrt((Pshoot(j,iterbw)%q(1,i)-qref(1,i))**2 + &
                      (Pshoot(j,iterbw)%q(2,i)-qref(2,i))**2 + &
                      (Pshoot(j,iterbw)%q(3,i)-qref(3,i))**2)
      if (absdist(i).gt.absdmax(j,iterbw)) then
        absdmax(j,iterbw) = absdist(i)
        atom_bouge_abs    = i
      endif

     enddo
     absdmax(j,iterbw)      = absdmax(j,iterbw)*1.0d8
     write(*,*) 'iterbw absdmax ',iterbw,absdmax(j,iterbw)
!!!!!!!!!!!! ATTENTION AUX RELATIONS DU BILAN DETAILLE

  if (absdmax(j,iterbw).lt.h_A_max) then
    rapport(j)=exp(alpha_bias(j)*(triallyap(j)-oldLyap(j)))
  else
   rapport(j)=0.d0
  endif

  ranf=genrand()
  mcconf=min(1.d0,rapport(j))
  if (ranf.lt.mcconf) then
    acc(j)=acc(j)+1
    Path(j,0:totiter) = Pshoot(j,0:totiter)
    oldLyap(j)=triallyap(j)
    
    write(*,*) 'traj. acceptée',mcmoves,'rapport =',rapport(j)
    absdmax_current = absdmax(j,iterbw)  
    write(*,*) 'absdmax_current', absdmax_current
  else
        write(*,*) 'traj. refusée ',mcmoves,'rapport =',rapport(j)
  endif
  write (*,*) 'Shooting ....taux d''acceptation',real(acc(j))/real(mcmoves) 
  write (*,*) 

! write(*,*) 'pos fin du shooting', Path(j,totiter)%qx(0), Path(j,totiter)%qy(0), Path(j,totiter)%qz(0)
! write(*,*) 'mom fin du shooting', Path(j,totiter)%px(0), Path(j,totiter)%py(0), Path(j,totiter)%pz(0)
! write(*,*) 'pos fin du shooting', Path(j,0)%qx(10), Path(j,0)%qy(10), Path(j,0)%qz(10)
! write(*,*) 'mom fin du shooting', Path(j,0)%px(10), Path(j,0)%py(10), Path(j,0)%pz(10)

enddo!! fin du shooting
!enddo  !! fin de icheck
!write(*,*) 'stop here'
!stop

 do j=depart_boucle_nbclones,Nbclones
!  pix=nint(genrand()*totiter) probleme aux limites
   pix = int(genrand()*(totiter+1)) ! le cas pix=totiter+1 est impossible par construction
                                    ! car genrand < 1
                                    ! pix est compris entre 0 et L soit L+1 valeur possibles 
  write(*,'(" Shifting............................:",3i5)') mcmoves,totiter,pix

  write(*,'(" Shift.rec...pix+1->pix+totinter....:",i5," to ",i5)') pix+1,pix+totiter
  ! we recycle the pix,pix+1, pix+2, ...., pix+totiter because there are totiter+1 points 
  do k=0,totiter
   Pshift(j,pix+k) = Path(j,k)
  enddo

  write(*,'(" shift.bw ......pix->0..............:",i5," to 0")') pix
    ! we calculate the pix-1,....,0
  new_projection=.false.
  it_trajectory=0

  Pcourant   = Pshift(j,pix)
  eigenvalue = Pcourant%eigenvalue

  do iterbw=pix,1,-1
    call mapping_P_Verlet(Pcourant%q,Pcourant%p,-dt,N,q1s2)  ! -dt car backward en position verlet
    call lanczos(N,maxvec,q1s2,new_projection,Pcourant%project)

    Pcourant%eigenvalue = eigenvalue

    if (eigenvalue.lt.0.0) then
      Pcourant%Lyap= asin(dt*sqrt(-eigenvalue)/2.d0)*2.d0 ! on stocke Lyap en 1/2 avant
     else 
      Pcourant%Lyap=0.d0  !! autoval max del Lanczos
    endif
    Pshift(j,iterbw-1)=Pcourant
    it_trajectory=it_trajectory+1 
    new_projection=.false.
   enddo

  write(*,'(" shift.fw..pix+totiter->2*totinter..:" ,i5, " to ", i5)') pix+totiter,2*totiter
  ! we calculate the pix+toiter+1, pix+2, ...., 2*toiter
  ! Check: there are two points 2*toiter-1 and 2*toiter which
  !         are computed for nothing. Never used after !

  timefsh=dt*real(totiter+pix)
  new_projection=.true.
  iterfw=nint(timefsh/dt)
  it_trajectory=0

  Pcourant = Pshift(j,totiter+pix)
  eigenvalue = Pcourant%eigenvalue

  do iterfw=totiter+pix,2*totiter-1
    call mapping_P_Verlet(Pcourant%q,Pcourant%p,dt,N,q1s2)
    call lanczos(N,maxvec,q1s2,new_projection,Pcourant%project)!!!!!
    Pcourant%eigenvalue = eigenvalue 

    if (eigenvalue.le.0.d0) then 
        Pcourant%Lyap=asin(dt*sqrt(-eigenvalue)/2.d0)*2.d0  ! Lyapunov en 1/2 stocké avant 
       else 
        Pcourant%Lyap=0.d0 ! Pshift(j,iterfw-1)%Lyap ! Lyapunov en 1/2 stocké avant 
    endif

    Pshift(j,iterfw)%Lyap       = Pcourant%Lyap
    Pshift(j,iterfw)%Project    = Pcourant%Project
    Pshift(j,iterfw)%eigenvalue = eigenvalue

    Pshift(j,iterfw+1)=Pcourant

    it_trajectory=it_trajectory+1 
    new_projection=.false.

  enddo

  hamilt(j)=ener0(j)

  S(j,:)=0
   do l=0,totiter  ! sommation sur les 1+L "path proposals" possibles
     do a=l,l+totiter-1  ! dans une action S il y a bien totiter=L Lyap stockes de 0 à totiter-1 
       S(j,l) = S(j,l)+Pshift(j,a)%Lyap
     enddo
   enddo
  S(j,:)=S(j,:)/(real(totiter))

  !!!!!!!!!!!! calcul de la fonction echelon pour le point de départ 
  !!!!!!!!!!!! de la trajectoire N sui 2N passi della traj shiftata

  absdmax(j,:)=-9999.0
  atom_bouge_abs=0

   do l=0,2*totiter ! il y a 1+L (1+totiter) proposals + totiter positions décalée en totiter
    do i=1,N
     absdist(i)=sqrt((Pshift(j,l)%q(1,i)-xref(i))**2 + & 
                     (Pshift(j,l)%q(2,i)-yref(i))**2 + & 
                     (Pshift(j,l)%q(3,i)-zref(i))**2) 
     if (absdist(i).gt.absdmax(j,l)) then 
        absdmax(j,l)=absdist(i) 
     endif 

! Manuel : les coordonnées dans Pshift(j,l) sont celles avant l'application du mapping P_Verlet dans le sens forward
! on a donc maintenant correspondance entre Pshift(j,l) et absdmax(j,l) contrairement à auparavant

    enddo
  enddo
  absdmax(j,:)=absdmax(j,:)*1.0e8
  write(*,*) 'absdmax : '
  write(*,*) absdmax(j,0) 
  !!!!!!!!!!!! ESSAY POUR MBAR: calcul des divers poids u_kln avec waste recycling 
  z = 0.0

  !  if (waste_recycling) then
  ! Check: Why many jl (bias) if we compute 
  !          only one in this case  ?
  ! Check: Can be a real problem here !!!! 
  ! Manuel : c'est normal, il faut calculer toutes ces grandeurs pour MBAR. 
  !   do jl=0,NbClones

    do jl=0,Nbclones_mbar
    z=0.d0
    do l=0,totiter ! boucles sur les chemins proposés possibles
  ! Here is P_sel from the paper
      if (absdmax(j,l).lt.h_A_max) then  ! on a correspondance maintenant entre absdmax et Pshift
        Psel(jl,l) = exp(alpha_bias(jl)*S(j,l))
             else
        Psel(jl,l)= 0.d0
      endif
      z = z +  Psel(jl,l)
    enddo
! This is S_alpha(zeta^m) Eq.48
    u_kln(j,jl,mcmoves) = -log(z) + log(real(totiter+1))
    Psel(jl,0:totiter)  =  Psel(jl,0:totiter)/z
!   write(*,*) 'Psel(',jl,')=',Psel(jl,totiter/2),Psel(jl,totiter-1)
   enddo

   umoy_kln(j,0:nmax,mcmoves)=0.d0
   u2moy_kln(j,0:nmax,mcmoves)=0.d0

   do jl=0,Nbclones_mbar
     do i=0,totiter !
       umoy_kln(j,jl,mcmoves)  = umoy_kln(j,jl,mcmoves)  - S(j,i)*Psel(jl,i)
       u2moy_kln(j,jl,mcmoves) = u2moy_kln(j,jl,mcmoves) + (S(j,i)**2)*Psel(jl,i)
     enddo
   enddo

! selection de la trajectoire indiciée newtraj
 
 xalea    = genrand()
 xcumul=0.d0

 do k=0,totiter !  boucle sur les "proposals"
   xcumul=xcumul+Psel(j,k)
   if (xalea.lt.xcumul) goto 234
 enddo
 234 continue
 newtraj=k
 write(*,*) 'xalea, Psel, poids cumulé et newtraj = ',xalea,Psel(j,k), xcumul ,newtraj
 write(*,*) ' absdmax(j,newtraj)', absdmax(j,newtraj)
 write(*,*) 'oldlyap(j)  = ', oldLyap(j)

!!!!!!!! on copie la trajectoire selectionnée avec le shifting
 if (absdmax(j,newtraj).le.h_A_max) then 
 Path(j,0:totiter-1) = Pshift(j,newtraj:totiter+newtraj-1) ! on a Path(j,0) = Pshift(j,0) pour newtraj=0
 oldLyap(j)       = SUM(Path(j,0:totiter-1)%Lyap)/real(totiter)
 else
 write(*,*) " Problème avec le shifting "
 endif

 write(*,*) 'newlyap(j)  = ', oldLyap(j)

!  calcul des fonctions indicatrices pour l'etat B qui correspond a la barriere dans le cas lacune 

  h_F(:)      = 0
  h_Fg(:)     = 0
  h_Fd(:)     = 0
  h_FI(:)     = 0
  h_dI(:)     = 0

  ! write(*,*) 'h_ba_min...',h_ba_min,h_ba,h_ba_max
  do k=0,2*totiter 
     !passagge F-d et F-I
     if ((absdmax(j,k).ge.h_ba_min).and.(absdmax(j,k).lt.h_ba)) then
        h_Fg(k) = 1.0
     endif
     if ((absdmax(j,k).ge.h_ba).and.(absdmax(j,k).lt.h_ba_max) ) then   !!!FCC -> DEFAULT FCC
        h_Fd(k) = 1.0
     endif
     if ((absdmax(j,k).ge.h_ba_max).and.(absdmax(j,k).lt.h_ba_I)) then !!FCC -> FCC 
        h_F(k)  = 1.0
     endif
     if (absdmax(j,k).ge.h_ba_I ) then 
        h_FI(k) = 1.0
     endif
!  write(*,*) ' h_ba... ', h_ba_min,  h_ba , h_ba_max, h_ba_I
     !
!  write(*,*) 'absdmax(j,k)',k, absdmax(j,k),h_F(k),h_Fg(k),h_Fd(k),h_FI(k)

  enddo !k

! calcul des constantes de reaction avec waste recycling

  do tprim=1,100!125!,700,100
!                 modif de Manuel : j'ai déplacé l'initialization à zero avant la première utilisation. C'est plus sûr.

   react_F(0)  = 0.0
   react_Fd(0) = 0.0
   react_dI(0) = 0.0
   react_FI(0) = 0.0
   react_Fg(0) = 0.0


  do k=0,totiter  ! manuel  attention au décalage par rapport aux h_F...
     kl=k + tprim*tprimo  ! -1 tient compte du décalage entre les h_* et Psel
     react_Fg(0) = react_Fg(0) +  h_Fg(kl) * Psel(0,k)
     react_Fd(0) = react_Fd(0) +  h_Fd(kl) * Psel(0,k)
     react_F(0)  = react_F(0)  +  h_F(kl)  * Psel(0,k)
     react_FI(0) = react_FI(0) +  h_FI(kl) * Psel(0,k)
     !write(*,*) ' totiter k kl ', totiter, k , kl
  enddo  
  write (30,'(i,i,i,4(e12.4))') tprim, j, mcmoves, react_F(0),   react_Fg(0),   react_Fd(0),   react_FI(0)
  kl = newtraj + tprim*tprimo 
  write (31,'(i,i,i,4(e12.4))') tprim, j, mcmoves, h_F(kl), h_Fg(kl), h_Fd(kl), h_FI(kl)
  enddo ! fin boucle tprimo

!  write (*,*) ' ',tprim, react_Fd(0),Psel(0,totiter-1)
!  on calcule les correlations avec t=300=tprimo*100 en fonction du biais

  react_Fd(:)= 0.0
  do jl=0,Nbclones_mbar
    do k=0,totiter
      kl = k + totiter
      react_Fd(jl) = react_Fd(jl) + (h_Fd(kl)+ h_F(kl)+h_FI(kl))*Psel(jl,k) ! attention c'est 3 contributions !!
      !write(*,*) ' verif ',  h_Fd(k) , Psel(jl,k)
      ! write(*,*) ' totiter k kl ', k , kl

    enddo
  enddo
  h_temp = h_Fd(newtraj+totiter) + h_F(newtraj+totiter) + h_FI(newtraj+totiter)

 !write(*,*) ' h_temp react_Fd ', h_temp , react_Fd(0)
 do l=0,NbClones_mbar
   ustd_kln(j,l,mcmoves) = -alpha_bias(l)*S(j,newtraj)
 enddo

! écriture des poids dans MBAR pour le fichier dada

 k=depart_boucle_nbclones

  call cal_hamilton(Path(j,0)%q,Path(j,0)%p,N,ekin,epot)
  hamilt(j)=(ekin+epot)/temperature
  ener0(j) = epot/temperature
  enerK(j) = ekin/temperature
  write(24,'(i,i,4(e18.7e3))') k,mcmoves,hamilt(j),ener0(j),enerK(j),-S(j,newtraj)
  
  
  
 do l=0,Nbclones_mbar
      write(244,'(i,i,i,6(e15.4e3))') k,l,mcmoves,u_kln(k,l,mcmoves),ustd_kln(k,l,mcmoves),umoy_kln(k,l,mcmoves),u2moy_kln(k,l,mcmoves),react_Fd(l),h_temp
 enddo
 
 if(.not.waste_recycling) then ! calcul des poids pour MBAR  post-shifting sans le waste-recycling 
!!    do l=0,NbClones ! Commentaire de Manuel : il y avait une erreur monumentale mais sans grosses consequences sur la statistique, 
!!    il faut tenir compte des canaux avec des numéros > nbclones Massimiliano passait des valeurs nulles 
!!!!!!!!!!!!!!!!!!!!!  calcul des distances 

 absdmax(j,:)=-9999.0
 atom_bouge_abs=0

  do iter=0,totiter-1
   !if (mod(iter,10).eq.0) then
   atom_bouge_abs=0
    q=Path(j,iter)%q
    do i=1,N
          absdist(i)=sqrt((q(1,i)-qref(1,i))**2+(q(2,i)-qref(2,i))**2+(q(3,i)-qref(3,i))**2)
       if (absdist(i).gt.absdmax(j,iter)) then
         absdmax(j,iter) = absdist(i)
         atom_bouge_abs  = i
       endif
    enddo
!debugC   write (*,'(a,i,i,i,e12.4)') 'bougee_abs', mcmoves, iter, atom_bouge_abs, absdmax(j,iter)*1.0e8
  enddo
  absdmax(j,:)=absdmax(j,:)*1.0e8
  write(*,*) ' coucou WR ' , waste_recycling
  z=0.0
 endif

!  mesure de la température
  temp=0.0000
  tempvar=temp
  do iterbw=totiter,0,-1
  tempiter = SUM((Path(j,iterbw)%p(1,1:N)**2+Path(j,iterbw)%p(2,1:N)**2+Path(j,iterbw)%p(3,1:N)**2))/m_i(1,1)/dble(3*im-6)/KtoERG ! *2/3 /2
   temp=temp+tempiter
   tempvar=tempvar+tempiter**2
  enddo
  temp=temp/dble(totiter+1)
  tempvar=sqrt(tempvar/dble(totiter+1)-temp**2)
  write(*,'("température cinétique trajectorielle et std :",2(e15.4e3)," et cible",e15.4e3)'),temp,tempvar,temperature/KtoERG
  write(*,'("température cinétique initiale              :",e15.4e3)'),tempiter

enddo ! boucle sur le clone
enddo !! mcmoves avec incrément + 1 clones

!write(*,*)'tau',  real(acc)/real(mcmoves)

do j=depart_boucle_nbclones, NbClones
  write(*,*) 'tau',j,alpha_bias(j),real(acc(j))/real(totalmcmoves)
  !write(*,*) 'acceptance ratio',j,alpha_bias(j),real(acc(j))/real(totalmcmoves)
enddo

posfinal = sortie(1:lenfnam)//'.cout'

open(unit=27, file=posfinal, status='replace')

do j=depart_boucle_nbclones,NbClones
 do i=0,totiter
!   write (27,*) Path(j,i)%q(1,1:im)
!   write (27,*) Path(j,i)%q(2,1:im)
!   write (27,*) Path(j,i)%q(3,1:im)
!   write (27,*)
!   write (27,*) Path(j,i)%p(1,1:im)
!   write (27,*) Path(j,i)%p(2,1:im)
!   write (27,*) Path(j,i)%p(3,1:im)
   write (27,*) Path(j,i)%q
   write (27,*)
   write (27,*) Path(j,i)%p
   write (27,*)
   write (27,*) Path(j,i)%Lyap
   write (27,*)
   write (27,*) Path(j,i)%project
   write (27,*)
   write (27,*) Path(j,i)%eigenvalue
enddo 
enddo

close(27)

stop

end subroutine LyapLanczos_vac

subroutine lanczos(N,maxvec,q1s2,new_projection,projection) !!!!!!!!!!!!! SERVE POTENZIALE PER LA LACUNA!!!!!!!!!!!!!!!!!!!!
  !use defs
  use gen_com_m, ONLY : erg2ev,ev2erg
  use lanczos_defs
  !use random_art
  !use art_in_ndm_module
  implicit none
  integer :: maxvec
  integer :: N
  logical ::  new_projection
  logical :: lanczos_failed=.false.
  real(8), dimension( 2 * maxvec -1 ) :: scratcha
  real(8), dimension(maxvec) :: diag
  real(8), dimension(maxvec-1) :: offdiag
  real(8), dimension(maxvec, maxvec) :: vector
  real(8), dimension(3*N, maxvec), target :: lanc
  real(8):: sum_forcenew, sum_force
  ! Vectors used to build the matrix for Lanzcos algorithm 
  real(8), dimension(:), pointer :: z0, z1, z2
  ! Projection direction based on lanczos computations of lowest eigenvalues

  integer :: i,k, i_err, ivec, nl_failed, it_art,evalf_number
  real(8) :: a1,a0,b2,b1,increment!, eigenvalue
  real(8) :: excited_energy,c1,norm
  real(8) :: xsum, ysum, zsum, sum2, invsum
  real(8), dimension(3*N) :: pos, newpos,newforce,ref_force
  real(8), dimension(3*N) :: newforcep1,newforcem1,projection
  real (double), dimension (1:3,1:N) :: q1s2  ! vector of position
  real(double), parameter:: cmTOang=1.0d8
  real(double)  :: ran3

!  do i =1,N
!  pos(i)=qx(i-1)*cmTOang
!  pos(i+N)=qy(i-1)*cmTOang
!  pos(i+2*N)=qz(i-1)*cmTOang
!  enddo

  pos(1:N)       = q1s2(1,1:N)*cmTOang
  pos(N+1:2*N)   = q1s2(2,1:N)*cmTOang
  pos(2*N+1:3*N) = q1s2(3,1:N)*cmTOang



  !lanczos_step =0.001 ! in Angstroems  mettre en parametre d'entree.
  !boxl(:) = box(:) * scala
  increment = lanczos_step  ! Increment, convert in box units
  overlap=0.0
  evalf_number=0
  nl_iter=0
  nl_failed=0

  if(.not. new_projection ) then
    old_before_sc_projection = projection*cmTOang*cmTOang
    projection=old_before_sc_projection  ! Vectorial operation
    eigenvalue=eigenvalue/(ev2erg*cmTOang*cmTOang)
  end if

  ! We now take the current position as the reference point and will make 
  ! a displacement in a random direction or using the previous direction as
  ! the starting point.
  it_art=it_trajectory
  !debugC write(*,*) 'lanczos',it_art,it_trajectory,mcmoves
  34 continue
   !d1 call calcforce(N,pos,ref_force,total_energy,it_art)
   !d1 evalf_number = evalf_number + 1
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
    call calcforce(N,newpos,newforcep1,excited_energy,it_art)
    evalf_number = evalf_number + 1

    !d4 newpos = pos +2.d0* z0 * increment
    !d4 call calcforce(N,newpos,newforcep2,excited_energy,it_art)
    !d4 evalf_number = evalf_number + 1

    newpos = pos - z0 * increment
    call calcforce(N,newpos,newforcem1,excited_energy,it_art)
    evalf_number = evalf_number + 1

    !d4 newpos = pos - 2.d0*z0 * increment
    !d4 call calcforce(N,newpos,newforcem2,excited_energy,it_art)
    !d4 evalf_number = evalf_number + 1
     
     !d4 newforce(:)=(newforcem2(:)-8.d0*newforcem1(:)+8.d0*newforcep1(:)-newforcep2(:))/12.d0
     newforce(:)=(newforcep1(:) - newforcem1(:))/2.0d0 
     !d2 newforce(:)= newforcep1(:)-ref_force(:)
    ! We extract lanczos(1)

!write(*,*) '44'!,maxval(newforce(:)), maxval(newforce1(:)), minval(newforce2(:))

  ! We get a0
  a0 = 0.0d0
  do i=1, 3*N
    a0 = a0 + z0(i) * newforce(i)
  end do
  diag(1) = a0
  z1 => lanc(:,2)
  z1 = newforce - a0 * z0    ! Vectorial operation
!write (*,*) z1

!stop
  b1=0.0d0
  do i=1, 3*N
    b1 = b1 + z1(i) * z1(i)
  end do
  offdiag(1) = sqrt(b1)

  invsum = 1.0d0 / sqrt ( b1 )
  z1 = z1 * invsum           ! Vectorial operation

  ! We can now repeat this game for the next vectors
  do ivec = 2, maxvec-1
    z1 => lanc(:,ivec)
    newpos = pos + z1 * increment
    call calcforce(N,newpos,newforcep1,excited_energy,it_art)
    evalf_number = evalf_number + 1

    !d4 newpos = pos + 2.d0*z1 * increment
    !d4 call calcforce(N,newpos,newforcep2,excited_energy,it_art)
    !d4 evalf_number = evalf_number + 1

     newpos = pos - z1 * increment
     call calcforce(N,newpos,newforcem1,excited_energy,it_art)
     evalf_number = evalf_number + 1

    !newpos = pos - 2.d0*z1 * increment
    !call calcforce(N,newpos,newforcem2,excited_energy,it_art)
    !evalf_number = evalf_number + 1


     !d4 newforce(:)=(newforcem2(:)-8.d0*newforcem1(:)+8.d0*newforcep1(:)-newforcep2(:))/12.d0
     newforce(:)=(newforcep1(:) - newforcem1(:))/2.0d0  
     !d1 newforce(:)=newforcep1(:)-ref_force(:)
 
    a1 = 0.0d0
    do i=1,3*N
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
!write(*,*) 'cinque', newpos(1)
  ! We now consider the last line of our matrix
  ivec = maxvec
  z1 => lanc(:,maxvec)
    newpos = pos + z1 * increment
    call calcforce(N,newpos,newforcep1,excited_energy,it_art)
    evalf_number = evalf_number + 1

    !d4 newpos = pos + 2.d0*z1 * increment
    !d4 call calcforce(N,newpos,newforcep2,excited_energy,it_art)
    !d4 evalf_number = evalf_number + 1

     newpos = pos - z1 * increment
     call calcforce(N,newpos,newforcem1,excited_energy,it_art)
     evalf_number = evalf_number + 1

    !d4 newpos = pos - 2.d0*z1 * increment
    !d4 call calcforce(N,newpos,newforcem2,excited_energy,it_art)
    !d4 evalf_number = evalf_number + 1


      !d4 newforce(:)=(newforcem2(:)-8.d0*newforcem1(:)+8.d0*newforcep1(:)-newforcep2(:))/12.d0
      newforce(:)=(newforcep1(:) - newforcem1(:))/2.0d0  
      !d1 newforce(:)=newforcep1(:)-ref_force(:)

    sum_force = 0.0 
    sum_forcenew = 0.0

    do i = 1 , 3*N
       sum_force = sum_force + ref_force(i)
       sum_forcenew =  sum_forcenew + newforce(i) 
    end do
!   write(*,*) 'the sum of the forces before the move', sum_force,sum_forcenew

!   write(*,*) 'the sum of the forces After the move', sum_forcenew
! newforce = newforce - ref_force
  sum_forcenew = 0.0
  a1 = 0.0d0
  do i=1, 3*N
    a1 = a1 + z1(i) * newforce(i)
    sum_forcenew =  sum_forcenew + newforce(i)
  end do
  !write(*,*) 'the difference between the forces ' , sum_forcenew
!  write(*,*)
  diag(maxvec) = a1

  ! We now have everything we need in order to diagonalise and find the
  ! eigenvectors.

  diag = -1.0d0 * diag
  offdiag = -1.0d0 * offdiag

  ! We now need the routines from Lapack. We define a few values
  i_err = 0

  ! We call the routine for diagonalizing a tridiagonal  matrix
  call dstev('V',maxvec,diag,offdiag,vector,maxvec,scratcha,i_err)

  ! We now reconstruct the eigenvectors in the real space
  ! Of course, we need only the first 5*N elements of vec

  projection = 0.0d0    ! Vectorial operation
  do k=1, maxvec
    z1 => lanc(:,k)
    a1 = vector(k,1)
    projection = projection + a1 * z1   ! Vectorial operation
  end do 
  c1=0.0d0
  do i=1, 3*N
     c1 = c1 + projection(i) * projection(i)
  end do

     norm = 1.0/sqrt(c1)
     projection = projection * norm 

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
overlap=a1 
!d if(a1<0.0d0) then
!d    projection = -1.0d0 * projection
!d    overlap=-overlap
!d end if    
!ooo
!ooo  
!d  call center(projection,3*N)
  
! write(*,*) 'eigenv', eigenvalue-old_eigenvalue, a1

!write(*,'("lanczos nliter",i5,f5.2,2f12.3)') nl_iter,overlap, eigenvalue,old_eigenvalue-eigenvalue
! 10-3 eV/A² is the same thing as 6.15 erg/cm². For the reasons we  will take 5.  
if ( dabs((old_eigenvalue-eigenvalue)) .gt. 1.d-3) then  ! mettre lanczos_threshold en parametre d'entree
   self_consistent=.true.
   lanczos_failed=.false.
   nl_iter=nl_iter+1
   nl_failed=nl_failed+1
   if (nl_failed.gt.30) then
      nl_failed=0
      !lanczos_failed=.true.
      write(*,*) 'WARNING: LANCZOS FAILED ... convergence not reached'
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

  eigenvalue=eigenvalue*ev2erg*cmTOang*cmTOang
  projection=projection/(cmTOang*cmTOang)
  overlap=a1 
  lanczos_iter=nl_iter*maxvec

!   The condition on the scalar product , we reject the point where we loose the eigen value and try to reduce the step size
!debug    if(abs(a1)<=0.2) reject = .true. 

 overlap=a1

!debug if(a1<0.0d0) then
!debug    projection = -1.0d0 * projection
    !overlap=-overlap
!debug end if    

 call center(projection,3*N)

!write(*,'("lanczos nliter",i5,f5.2,2f12.3)') nl_iter,overlap, eigenvalue,old_eigenvalue
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

  xtotal   = xtotal / natoms
  ytotal   = ytotal / natoms
  ztotal   = ztotal / natoms

  do i = 1, natoms
    x(i)   = x(i) - xtotal
    y(i)   = y(i) - ytotal
    z(i)   = z(i) - ztotal
  end do
end subroutine

Subroutine mapping_P_Verlet(qx,px,dt,N,q1s2) !!!!!!!!!!! MODIFIER AVEC BON POTENTIEL

 USE T_kind_param_m, ONLY:  double
    use gen_com_m
    !use jqmod
    use tab_imm_m
 implicit none
 integer :: ic
 integer   ::N
 real (double), dimension (3,1:N) :: qx   ! vector of position
 real (double), dimension (3,1:N) :: px   ! vector of impulsion
 real (double), dimension (3,1:N) :: q1s2 ! vector of intermediate position 
 real (double) ::dt
 real (double) :: pp(3,N),xbar(3),dxbar(3)

 xp(1:3,1:N)=qx(1:3,1:N)
 pp(1:3,1:N)=px(1:3,1:N)
 
 do ic=1,3
 pp(ic,1:im) = pp(ic,1:im) - sum(pp(ic,1:im))/dble(im)
 enddo
 call control_angular_momenta(pp,qx) 

 do ic=1,3
   xbar(ic)    = sum(xp(ic,1:im))/dble(im) ! barycentre sur les particules de même masse
 enddo

 vp(1:3,1:im)  = pp(1:3,1:im)/m_i(1:3,1:im)

 xp(1:3,1:im)  = xp(1:3,1:im) + vp(1:3,1:im)*dt/two

 do ic=1,3
   dxbar(ic)    = sum(xp(ic,1:im))/dble(im) - xbar(ic) ! déplacement du barycentre 
   xp(ic,1:im) = xp(ic,1:im) - dxbar(ic)  ! on recentre tout le systeme
 enddo

 q1s2(1:3,1:N) = xp(1:3,1:N)

 call calfo_teledyn(it_trajectory)
 
!  ecin=  SUM((pp(1:3,1:im) + fp(1:3,1:im)*dt)**2/m_i(1:3,1:im))/2.00
!  write(*,*) 'potist ' ,potist,ecin,potist+ecin
    
 pp(1:3,1:im) =  pp(1:3,1:im) + fp(1:3,1:im)*dt

 do ic=1,3
 pp(ic,1:im) = pp(ic,1:im) - sum(pp(ic,1:im))/dble(im)
 enddo

 call control_angular_momenta(pp,qx) 


 vp(1:3,1:im) =  pp(1:3,1:im)/m_i(1:3,1:im)

 xp(1:3,1:im) =  xp(1:3,1:im) + vp(1:3,1:im)*dt/two

!debug write (*,*) 'lang1', xp (1:3,1),fp(1,1:3)

 do ic=1,3
   dxbar(ic)    = sum(xp(ic,1:im))/dble(im)-xbar(ic) ! déplacement du barycentre 
   xp(ic,1:im)  = xp(ic,1:im) - dxbar(ic)  ! on recentre tout le systeme
 enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 qx(1:3,1:N)=xp(1:3,1:N)
 px(1:3,1:N)=pp(1:3,1:N)

! write(*,*) 'two dt ',two,dt  test : dt change bien de signe pour le backward

end subroutine mapping_P_Verlet


Subroutine cal_hamilton(qx,px,N,ekin,epot) !!!!!!!!!!! MODIFIER AVEC BON POTENTIEL

 USE T_kind_param_m, ONLY:  double
    use gen_com_m
    !use jqmod
    use tab_imm_m
 implicit none
 integer   ::N
 real (double), dimension (1:3,1:N) :: qx   ! vector of position
 real (double), dimension (1:3,1:N) :: px   ! vector of impulsion
 real (double) :: pp(1:3,1:N),ekin,epot

 xp(1:3,1:N)=qx(1:3,1:N)
 pp(1:3,1:N)=px(1:3,1:N) 

 call calfo_teledyn(it_trajectory)

 ekin = SUM((pp(1:3,1:im))**2/m_i(1:3,1:im))/2.00
 epot = potist

! write(*,*) 'potist ' ,potist,ecin,potist+ecin

end subroutine cal_hamilton

Subroutine gauss(ss, l0,l)
real (double) :: ss
real (double) :: l0
real (double) :: l
real (double) :: r,v1,v2
real (double) :: x1,x2

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

Subroutine OU_control(p,q,a_sto,v0)
real(double), dimension(3,im) :: p,q,a
real(double):: z1,z2,z3,z4,v0,a_sto,r
integer:: i,ic

do i=1,im 
 do ic=1,3
 r=2.d0
 do while (r.ge.1.d0)
   call random_number(z1)
   call random_number(z2)
   z3=2.0*z1-1.0
   z4=2.0*z2-1.0
   r=z3**2+z4**2
 enddo
 a(ic,i)=z3*sqrt(-2.0*log(r)/(r))
 enddo
enddo


a(1,1:im) = a(1,1:im) - sum(a(1,1:im))/dble(im)
a(2,1:im) = a(2,1:im) - sum(a(2,1:im))/dble(im)
a(3,1:im) = a(3,1:im) - sum(a(3,1:im))/dble(im)

call control_angular_momenta(a,q)


p(1,1:im) = p(1,1:im) - sum(p(1,1:im))/dble(im)
p(2,1:im) = p(2,1:im) - sum(p(2,1:im))/dble(im)
p(3,1:im) = p(3,1:im) - sum(p(3,1:im))/dble(im)
call control_angular_momenta(p,q)

!p = p*a_sto + a*v0*sqrt(2.d0-2.d0*a_sto**2)  pour van brutzel 
p = p*a_sto + a*v0*sqrt(1.d0-1.d0*a_sto**2)


end subroutine OU_control

Subroutine control_angular_momenta(p,q)

implicit none
real(double) :: rx, ry, rz, r2x, r2y, r2z, r2
real(double) :: prx, pry, prz, px, py, pz, vrx, vry, vrz
real(double) :: omegax, omegay, omegaz
real(double), dimension(3) ::   scom
real(double), dimension(3,3) :: ainer, aineri
real(double), dimension(3,im) :: p,q
integer i,ic,ib

ainer = 0.d0
prx = 0.d0
pry = 0.d0
prz = 0.d0

!rang =0

! provisoire: doit tenir compte des masses

do ic=1,3
 scom(ic) = sum(q(ic,1:im))/dble(im)
enddo

!write(*,*) 'p : ', p(1:3,1:im) 

do i = 1, im
                           
 rx = q(1,i)-scom(1)
 ry = q(2,i)-scom(2)
 rz = q(3,i)-scom(3)
             
 r2x = rx*rx
 r2y = ry*ry
 r2z = rz*rz
 r2  = r2x+r2y+r2z
 
 ainer(1,1) = ainer(1,1)+m_i(1,i)*(r2-r2x)
 ainer(2,2) = ainer(2,2)+m_i(1,i)*(r2-r2y)
 ainer(3,3) = ainer(3,3)+m_i(1,i)*(r2-r2z)
 ainer(2,3) = ainer(2,3)-m_i(1,i)*ry*rz
 ainer(3,1) = ainer(3,1)-m_i(1,i)*rz*rx
 ainer(1,2) = ainer(1,2)-m_i(1,i)*rx*ry
 px  = p(1,i)
 py  = p(2,i)
 pz  = p(3,i)
 prx = prx+ry*pz-rz*py
 pry = pry+rz*px-rx*pz
 prz = prz+rx*py-ry*px
enddo

ainer(3,2) = ainer(2,3)
ainer(1,3) = ainer(3,1)
ainer(2,1) = ainer(1,2)

if (rang==0) then
       write(6,*) 
       write(6,997) (ainer(1,ib),ib=1,3),prx
       write(6,997) (ainer(2,ib),ib=1,3),pry
       write(6,997) (ainer(3,ib),ib=1,3),prz
endif
997    format('Inertia/anglm = ',3e15.6,4x,e14.6e3)

     !     calculate  angular velocity

call matinv(ainer,aineri)
omegax = aineri(1,1)*prx+aineri(1,2)*pry+aineri(1,3)*prz
omegay = aineri(2,1)*prx+aineri(2,2)*pry+aineri(2,3)*prz
omegaz = aineri(3,1)*prx+aineri(3,2)*pry+aineri(3,3)*prz

           !         shift velocities to make the angular momentum zero
do i = 1, im
 rx = q(1,i)-scom(1)
 ry = q(2,i)-scom(2)
 rz = q(3,i)-scom(3)
 vrx = omegay*rz-omegaz*ry
 vry = omegaz*rx-omegax*rz
 vrz = omegax*ry-omegay*rx

 p(1,i) = p(1,i)-vrx*m_i(1,i)
 p(2,i) = p(2,i)-vry*m_i(1,i)
 p(3,i) = p(3,i)-vrz*m_i(1,i)
enddo

end subroutine control_angular_momenta

subroutine calfo_teledyn(it_counter)
 USE T_kind_param_m, ONLY:  double 
 use gen_com_m
!use art_in_ndm_module
 use tab_imm_m
 integer , intent(in) :: it_counter
 
        !write(*,*) ltabvois,it_counter ,itab,itetabvois
        if (itab/=0) then
            if (mod(it_counter,itab)==0) then
               call caltabt
            endif
         endif
         if (ltabvois.and.mod(it_counter,itetabvois)==0) call caltabi 
         call calfo 

return
end subroutine

subroutine calcforce(N,pos,nforce,potist, it_art) !calcul des forces POUR LANCZOS 
! cette subroutine prends N, les positions, et doit rendre position, forces et l'energie de la configuration
 !use  
      USE T_kind_param_m, ONLY:  double
      use gen_com_m
      !use art_in_ndm_module
      use tab_imm_m

  implicit none
  integer, intent(in):: N
  real(double), dimension(3*N), target, intent(in) :: pos
  real(double), dimension(:), pointer :: x , y , z
  real(double), intent(out)  :: nforce(3*N)
  real(double)  :: potist
  integer  :: i,ic, it_art
  real(double), parameter:: cmTOang=1.d8
  real(double), parameter :: ev2erg=1.602d-12, erg2eV=1.d0/eV2erg


    x => pos(1:N)
    y => pos(N+1:2*N)
    z => pos(2*N+1:3*N)
 
do i=1,N
  xp(1,i)=x(i)/cmTOang
  xp(2,i)=y(i)/cmTOang
  xp(3,i)=z(i)/cmTOang
enddo

!write(*,*) 'force:itart', it_art

        if (itab/=0) then
            if (mod(it_art,itab)==0) then
               call caltabt
            endif
         endif
         if (ltabvois.and.mod(it_art,itetabvois)==0) call caltabi 
         call calfo 

do i=1,N
  nforce(i)=fp(1,i)*ergTOev/cmTOang
  nforce(i+N)=fp(2,i)*ergTOev/cmTOang
  nforce(i+2*N)=fp(3,i)*ergTOev/cmTOang
enddo

end subroutine calcforce


end module sundae_module



