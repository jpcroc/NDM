module sundae_module
	!-----------------------------------------------
	USE T_kind_param_m, ONLY: double
	use gen_com_m
	use var_pot, ONLY: cm
	use jqmod
	use random_art
	use lanczos_defs
#if(PARASUN) 
	use mpi
#endif(PARASUN)
	!-----------------------------------------------
	! test
	!   _c -> courant
	!   _s -> selection
	!   _d -> depart

		integer,parameter :: dpkind=selected_real_kind(13)
		real(double) :: t0,t1,elaps_1,tbuffer1,tbuffer2,t_lanczos
		integer,     dimension(:),allocatable,save  :: ipovois
		
		
		integer :: continue_sundae, maxvec, tprimo, totiter
		integer :: Totalmcmoves, Nbclones, Nbclones_mbar, depart_boucle_nbclones
		real(double) :: h_A_max, h_ba_max, h_ba_min, h_ba, h_ba_I, h_temp
		real(double) :: dt, TotalTime, gamma_sundae, Temperature
		real(double) :: alpha_max, teq, delta_x, a_sto, kapa, ss, tequilib
		character(len=128) :: sortie
		character(len=128) :: fnamtin
		character(len=128) :: posfinal
		character(len=128) :: data_abf
		character(len=128) :: data_mbar
		character(len=128) :: dada_mbar
		character(len=128) :: data_mbar_std
		character(len=128) :: moyennes_mbar
		character(len=128) :: moyennes_mbar_denom
		character(len=128) :: kappaF
		character(len=128) :: kappaFd
	
		real(double), parameter :: KtoERG=1.3791946308724831d-16

		real(double),dimension(:,:),allocatable,save  :: m_i
		real(double),dimension(:,:),allocatable,save  :: gau
		real(double),dimension(:,:),allocatable,save  :: rga_i


		real(double),dimension(:),allocatable,save    :: d2vois
		real(double),dimension(:,:),allocatable,save  :: xpvois
		real(double),dimension(:,:),allocatable,save  :: xpvoisini
		real(double),dimension(:),allocatable,save    :: xtransla


		integer :: iteration
		integer :: icheck

		integer N
		PARAMETER(N=127)   !1023 
		integer nfenetre
		parameter(nfenetre=100)
		integer :: nl_iter

		! cosmin added:
		integer    :: it_langevin_deter=0, it_langevin=0,it_trajectory    
		
		!
		!real (double), dimension(3,1:im):: q,p 
		real (double), dimension(:,:), allocatable, save :: q,p 
		type Trajectoire
			real (double), dimension(1:3,1:N):: q  ! vector of position
			real (double), dimension(1:3,1:N):: p  ! vector of impulsion
			real (double), dimension(3*N) :: project
			real (double)  Lyap
			real (double)  eigenvalue
		endtype

		type (Trajectoire), dimension (:,:), allocatable :: Path
		type (Trajectoire), dimension (:,:), allocatable :: Pshoot
		type (Trajectoire), dimension (:,:), allocatable :: Pshift
		type (Trajectoire) :: Pcourant
		
		real (double), dimension(1:3,1:N):: qref,qrot,prot
		
		integer :: i,j,l,k,iter, atom_bouge_abs
		real (double), dimension(1:3,1:N):: q1s2
		real (double), dimension(:), allocatable:: alpha_bias
		logical :: new_projection
		real(double) :: pav(3)
		real (double), dimension(:),   allocatable :: oldLyap 
		real (double), dimension(:,:), allocatable :: absdmax
		
		
		integer :: jl, kl


		real (double), dimension(1:N):: xref
		real (double), dimension(1:N):: yref
		real (double), dimension(1:N):: zref

		real (double), dimension(1:N,6):: qtemp 

		real (double), dimension(1:3):: xbar

		real (double) :: xalea,xcumul

	

		real (double) ::mcconf,ranf, z
		integer :: a, ltot,mcmoves, newtraj, it_art

		logical :: waste_recycling

		real(double) :: rga
		!real(double) :: pi
		integer :: ix


		real(double) :: e, timefsh
		real(double) :: ekin,epot
		real(double) :: absdmax_current
		integer :: pix
		integer :: iterfw, iterbw,nmax, tprim

		real (double) ,dimension(:),   allocatable:: absdist
		real (double) ,dimension(:),   allocatable:: ener0 
		real (double) ,dimension(:),   allocatable:: enerK
		real (double) ,dimension(:),   allocatable::triallyap
		real (double) ,dimension(:),   allocatable::rapport
		real (double) ,dimension(:),   allocatable:: hamilt
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
		real (double) :: tempiter,tempvar

		! Declarations ABF
		
		real(double) :: delta, delta_bin_theta_s2
		real(double) :: AR_MH
		real(double) :: theta_n, theta_temp
		integer      :: itheta_n, itheta_temp, theta_tilde, N_extra
		real(double) :: somme
		real(double) :: Lyap, gmax, gmin, P_n, pi_n, pi_n1, sum_A
		real(double), dimension(:), allocatable, save :: theta
		real(double), dimension(:), allocatable, save :: u_A
		real(double), dimension(:), allocatable, save :: A_n
		real(double), dimension(:), allocatable, save :: P_A
		real(double), dimension(:), allocatable, save :: sum_P_A
		real(double), dimension(:), allocatable, save :: A_prime
		real(double), dimension(:), allocatable, save :: A_prime_num
		real(double), dimension(:), allocatable, save :: O_moy_num
		real(double), dimension(:), allocatable, save :: O_moy_estim
		real(double), dimension(:), allocatable, save :: O_estim
		real(double) ,dimension(:), allocatable, save :: histo_theta	


		! Complements MPI
		
#if(PARASUN)
		integer rank, numproc, ierror
		real(double) ,dimension(:), allocatable, save :: MPI_histo_theta
		real(double) ,dimension(:), allocatable, save :: MPI_P_A
		real(double) ,dimension(:), allocatable, save :: MPI_sum_P_A 
		real(double) ,dimension(:), allocatable, save :: MPI_A_prime 
		real(double) ,dimension(:), allocatable, save :: MPI_A_prime_num
		real(double), dimension(:), allocatable, save :: MPI_O_moy_num
		real(double), dimension(:), allocatable, save :: MPI_O_moy_estim
		real(double), dimension(:), allocatable, save :: MPI_O_estim
#endif(PARASUN)
		
CONTAINS


!!!!!!!!!!!!!!!!!!! 21.05.14
!
! allocate_tele_vac()
!
! init_tele_vac()
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!  





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine LyapLanczos_shooting

	
	implicit none
	
	real(double) :: genrand

	j = itheta_n
	write(*,*) 'indice theta', itheta_n
	!!!!! initialization backward

	!  verif : Massimiliano avait oublie de mettre a jour la variable oldlyap
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

	if (ix.gt.0) then 
		Pcourant%project     = Path(j,ix)%project
		Pcourant%eigenvalue  = Path(j,ix)%eigenvalue
	else
		Pcourant%project     = Path(j,0)%project
		Pcourant%eigenvalue  = Path(j,0)%eigenvalue
	endif
	lanczos_iter=0
	do iterbw=ix,1,-1 ! on stocke  en iterbw-1 la valeur propre calculée à iterbw-1/2
		call mapping_P_Verlet(Pcourant%q,Pcourant%p,-dt,N,q1s2) !! on met -dt pour le backward !!!
		call lanczos(N,maxvec,q1s2,new_projection,Pcourant%project)
		lanczos_iter=lanczos_iter+nl_iter
		Pcourant%eigenvalue = eigenvalue 
		if (eigenvalue.lt.0.0) then 
			Pcourant%Lyap=asin(dt*sqrt(-eigenvalue)/2.d0)*2.d0 !
		else 
			Pcourant%Lyap=0.d0
		endif
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

	triallyap(j)=SUM(Pshoot(j,0:totiter-1)%Lyap)/real(totiter)
	write(*,*) 'oldlyap(j)   = ', oldLyap(j)*300.0/sqrt(9.270914743200000e-023)
	write(*,*) 'triallyap(j) = ', triallyap(j)*300.0/sqrt(9.270914743200000e-023)
	
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

	!! fin du shooting

end subroutine LyapLanczos_shooting


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine LyapLanczos_shifting

	implicit none
	
	real(double) :: genrand
	
	
	j = itheta_n

	pix = int(genrand()*(totiter+1))  ! le cas pix=totiter+1 est impossible par construction
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
	
	
	z=0.d0
	do l=0,totiter ! boucles sur les chemins proposés possibles
		! Here is P_sel from the paper
		if (absdmax(j,l).lt.h_A_max) then  ! on a correspondance maintenant entre absdmax et Pshift
			Psel(j,l) = exp(alpha_bias(j)*S(j,l))
		else
			Psel(j,l)= 0.d0
		endif
		z = z +  Psel(j,l)
	enddo

	Psel(j,0:totiter)  =  Psel(j,0:totiter)/z

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
	write(*,*) 'absdmax(j,newtraj)', absdmax(j,newtraj)
	write(*,*) 'oldlyap(j)  = ', oldLyap(j)*300.0/sqrt(9.270914743200000e-023)

	!!!!!!!! on copie la trajectoire selectionnée avec le shifting
	if (absdmax(j,newtraj).le.h_A_max) then 
		Path(j,0:totiter-1) = Pshift(j,newtraj:totiter+newtraj-1) ! on a Path(j,0) = Pshift(j,0) pour newtraj=0
		oldLyap(j)       = SUM(Path(j,0:totiter-1)%Lyap)/real(totiter)
	else
		write(*,*) " Problème avec le shifting "
	endif

	write(*,*) 'newlyap(j)  = ', oldLyap(j)*300.0/sqrt(9.270914743200000e-023)
	
	
	!!!!!!!!!!!! ESSAY POUR MBAR: calcul des divers poids u_kln avec waste recycling 
	
	!
	!  call essai_mbar
	!
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	

end subroutine LyapLanczos_shifting


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine LyapLanczos_ABF

	use lanczos_defs
	use tab_imm_m

	implicit none

	real(double) :: genrand
	

	! Proposition du nouveau theta 
	
	xalea = genrand()
	theta_temp = theta_n + 2.*delta*(xalea-0.5d0)
	
	
	itheta_n = nint((theta_n)/(alpha_max)*real(Nbclones_mbar))
	itheta_temp = nint((theta_temp)/(alpha_max)*real(Nbclones_mbar))
	

	! Acceptation/rejet 
	if ( (itheta_temp.ge.0) .and. (itheta_temp.le.Nbclones_mbar) ) then
		AR_MH = max( 1.0d0 , exp( -(theta_temp-theta_n)*Lyap - (A_n(itheta_temp) - A_n(itheta_n)) ) )
	endif
	ranf = genrand()
	if ( (ranf.lt.AR_MH) .and. (theta_temp.gt.0.0d0) .and. (theta_temp.lt.alpha_max) ) then
		theta_n = theta_temp
		itheta_n = nint((theta_n)/(alpha_max)*real(Nbclones_mbar))
	endif

	! Calcul du nouveau biais
	
	Lyap          = oldLyap(itheta_n)
	u_A(:)        = Lyap*theta(:) - A_n(:)     
	gmax          = maxval(u_A)
	u_A(:)        = u_A(:) - gmax
	P_A(:)        = exp(u_A(:))
	somme         = sum(P_A(-N_extra:Nbclones_mbar+N_extra))
	P_A           = P_A/somme


#if(PARASUN)
	! Mise en commun parallele ---- 1/2 !	
	call MPI_REDUCE(P_A,MPI_P_A,Nbclones_mbar+2*N_extra+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, MPI_COMM_WORLD,ierror)

	MPI_P_A = MPI_P_A/dble(numproc)
	if (rank.eq.0) then
		MPI_histo_theta(0:Nbclones_mbar) = MPI_histo_theta(0:Nbclones_mbar) + MPI_P_A(0:Nbclones_mbar)
	end if

	call MPI_BCAST(MPI_P_A,Nbclones_mbar+2*N_extra+1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
	! Mise en commun parallele ---- 1/2 !
#else
	histo_theta(0:Nbclones_mbar) = histo_theta(0:Nbclones_mbar) + P_A(0:Nbclones_mbar)
#endif(PARASUN)
	
	
	do theta_tilde=-N_extra,Nbclones_mbar+N_extra
		P_n        = P_A(theta_tilde)
		pi_n       = sum_P_A(theta_tilde)
		pi_n1      = pi_n + P_n

		if (pi_n1.eq.0.d0) then 
			pi_n1 = 1.0d-10
		endif
		
		sum_P_A(theta_tilde)       =   pi_n1

		A_prime_num(theta_tilde)   =   P_n * Lyap + A_prime_num(theta_tilde)
		!A_prime(theta_tilde)       =   A_prime_num(theta_tilde) / pi_n1
	enddo
	
	do k = 0,totiter 
		O_moy_num(k)     =   P_A(0) * O_estim(k) + O_moy_num(k)
	enddo
	

#if(PARASUN)
	! Mise en commun parallele ---- 2/2 !
	call MPI_REDUCE(A_prime_num,MPI_A_prime_num,Nbclones_mbar+2*N_extra+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, MPI_COMM_WORLD,ierror)
	call MPI_REDUCE(O_moy_num,MPI_O_moy_num,totiter,MPI_DOUBLE_PRECISION,MPI_SUM,0, MPI_COMM_WORLD,ierror)
	call MPI_REDUCE(sum_P_A,MPI_sum_P_A,Nbclones_mbar+2*N_extra+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, MPI_COMM_WORLD,ierror)

	MPI_A_prime = MPI_A_prime_num/ (MPI_sum_P_A + 1.d-6)
	
	do k = 0,totiter 
		O_moy_estim(k)  =   MPI_O_moy_num(k) / (MPI_sum_P_A(0) + 1.d-6)
	enddo
	
	if (N_extra.gt.0) then 
		MPI_A_prime(-N_extra:0) = 0.d0
		MPI_A_prime(Nbclones_mbar:Nbclones_mbar+N_extra) = 0.d0
	endif

	call MPI_BCAST(MPI_A_prime,Nbclones_mbar+2*N_extra+1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
	! Mise en commun parallele ---- 2/2 !

	A_n = 0.d0
	A_n(-N_extra) = 0.d0
	do theta_tilde  = -N_extra+1,Nbclones_mbar + N_extra
		A_n(theta_tilde) = A_n(theta_tilde-1)  + (MPI_A_prime(theta_tilde-1) + MPI_A_prime(theta_tilde))*delta_bin_theta_s2
	enddo
#else
	do theta_tilde=-N_extra,Nbclones_mbar+N_extra
		A_prime(theta_tilde)  =  A_prime_num(theta_tilde) / (pi_n1 + 1.d-6)
	enddo
	do k = 0,totiter 
		O_moy_estim(k)  =  O_moy_num(k) / (sum_P_A(0) + 1.d-6)
	enddo
	A_n = 0.d0
	A_n(-N_extra) = 0.d0
	do theta_tilde  = -N_extra+1,Nbclones_mbar + N_extra
		A_n(theta_tilde) = A_n(theta_tilde-1)  + (A_prime(theta_tilde-1) + A_prime(theta_tilde))*delta_bin_theta_s2
	enddo
#endif(PARASUN)


	gmin=minval(A_n)
	A_n = A_n - gmin

	! normalisation du generateur biaisant

	sum_A = log(sum(exp(-A_n(-N_extra:Nbclones_mbar+N_extra))))
	A_n = A_n + sum_A
	
	write(*,*) 'theta_n  = ', theta_n*sqrt(9.270914743200000e-023)/300.0


end subroutine LyapLanczos_ABF



subroutine LyapLanczos_Obs

	use lanczos_defs
	use tab_imm_m

	implicit none

	j = itheta_n

	!  calcul des fonctions indicatrices pour l'etat B qui correspond a la barriere dans le cas lacune 
	h_F(:)      = 0
	h_Fg(:)     = 0
	h_Fd(:)     = 0
	h_FI(:)     = 0
	h_dI(:)     = 0

	do k = 0,totiter 
		!passagge F-d et F-I
		if ((absdmax(j,k).ge.h_ba_min).and.(absdmax(j,k).lt.h_ba)) then
			h_Fg(k) = 1.0
		endif
		if ((absdmax(j,k).ge.h_ba).and.(absdmax(j,k).lt.h_ba_max) ) then   !!!FCC -> DEFAULT FCC
			h_Fd(k) = 1.0
			write(*,*) "*************************************************"
			write(*,*) "************                         ************"
			write(*,*) "************  passage dans l'etat B  ************"
			write(*,*) "************                         ************"
			write(*,*) "*************************************************"
		endif
		if ((absdmax(j,k).ge.h_ba_max).and.(absdmax(j,k).lt.h_ba_I)) then !!FCC -> FCC 
			h_F(k)  = 1.0
		endif
		if (absdmax(j,k).ge.h_ba_I ) then 
			h_FI(k) = 1.0
		endif
	enddo !k
	
	do k = 0,totiter
		O_estim(k) = h_Fd(k)
	enddo
	

end subroutine LyapLanczos_Obs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




subroutine LyapLanczos_vac! (xp)

	use lanczos_defs
	use tab_imm_m

	implicit none

	real(double) :: genrand
	
	
	!!!!!!!!!!!!!!!!!!!!!!      Phase d'initialisation       !!!!!!!!!!!!!!!!!!!!!!!!

	call read_sundae                         !!! Modif 21.05.14

    call LyapLanczos_allocate                !!! Modif 02.06.14

	call LyapLanczos_equilibrage             !!! Modif 02.06.14
	
	call LyapLanczos_init_tests              !!! Modif 23.05.14
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if(PARASUN) 
	
	call MPI_INIT(ierror)     !! On initialise MPI

	call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror) !! Attribue a rank le numero de processeur

	call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)  !! Attribue a numproc le nombre de processeur
	
	call init_random_seed(100000*(rank+1))
	
#else

	call init_random_seed(100000)

#endif(PARASUN)


	do mcmoves = 1 , Totalmcmoves                       !mcmoves = 1,M in the paper
	
		
		!write(*,'(" Proc number ", i5, "says hello")'), rank
	
	
		it_art=mcmoves

		ix=int((totiter+1)*genrand())
		write(*,*) 
		write(*,*) 
		write(*,'(" Itération..:",i5," sur un total de ",i5)') ,mcmoves,Totalmcmoves
		write(*,'(" Shooting......................:",3i5)') mcmoves,totiter,ix

		!!!! shooting deterministique
		! Check: if depart_boucle_nbclones,NbClones of if there are different WHY?


		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!
		
		call LyapLanczos_shooting              !!! Modif 02.06.14
	
	
		call LyapLanczos_shifting              !!! Modif 03.06.14
	
	
		call LyapLanczos_Obs				   !!! Modif 13.06.14
	
	
		call LyapLanczos_ABF				   !!! Modif 10.06.14
		
		!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	enddo !! mcmoves avec incrément + 1 clones


#if(PARASUN) 

	if (rank.eq.0) then
		do j=0,Nbclones_mbar
			write(111,*) MPI_histo_theta(j), MPI_A_prime(j), O_moy_estim(j)
		enddo
	endif
	call MPI_FINALIZE(ierror)	 !! Pour finir l'appel MPI
	
#else

	do j=0,Nbclones_mbar
		write(111,*) histo_theta(j), A_prime(j), O_moy_estim(j)
	enddo

#endif(PARASUN) 


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	do j=depart_boucle_nbclones, NbClones
		write(*,*) 'tau',j,alpha_bias(j),real(acc(j))/real(Totalmcmoves)
	enddo

	posfinal = sortie(1:lenfnam)//'.cout'

	open(unit=27, file=posfinal, status='replace')

	do j=depart_boucle_nbclones,NbClones
		do i=0,totiter
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine LyapLanczos_init_tests

	!use lanczos_defs
	use tab_imm_m
	implicit none

	! test du moment angulaire de la force
	call caltabt
	call caltabi
	call calfo

	p(1:3,1:im) = fp(1:3,1:im)
	q(1:3,1:im) = xp(1:3,1:im)
	rang=0
	call control_angular_momenta(p,q)
	call control_angular_momenta(p,q)
	call control_angular_momenta(p,q)
	! stop


	! test de la subroutine de controle angulaire
	!do j=Nbclones_mbar,Nbclones_mbar
	do j = itheta_n, itheta_n
	
		q= Path(j,0)%q
		p= Path(j,0)%p

		temp=0.d0
		do i=1,20
		call OU_control(p,q,a_sto,ss)
		call mapping_P_Verlet(q,p,dt,N,q1s2)

		temp = temp + (sum(p(1,1:im)**2)+sum(p(2,1:im)**2)+sum(p(3,1:im)**2))/dble(3*im-6)/m_i(1,1)
		write(*,*) 'température ',temp/dble(i)/KtoERG,i,alpha_bias(j)
		enddo
	
	
		write (*,'("1st traj before shooting ........:",i5)') nint(TotalTime/dt)
		rang=1

		call cpu_time(t0)
		lanczos_iter=0

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

		do iterbw=0,totiter
			tempiter = SUM((Path(j,iterbw)%p(1,1:N)**2+Path(j,iterbw)%p(2,1:N)**2+Path(j,iterbw)%p(3,1:N)**2))/m_i(1,1)/3.000/dble(im-2)/KtoERG ! *2/3 /2
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
	
	enddo

end subroutine LyapLanczos_init_tests


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine LyapLanczos_equilibrage

	use tab_imm_m

	implicit none
	
	if (continue_sundae.ne.2) then 
	
		do j = 0, Nbclones_mbar
			absdmax(j,:)=-9999.0
			atom_bouge_abs=0
			iter=0
			do i=1,N
				absdist(i)=sqrt((xp(1,i)-qref(1,i))**2+(xp(2,i)-qref(2,i))**2+(xp(3,i)-qref(3,i))**2)
				if (absdist(i).gt.absdmax(j,iter)) then
					absdmax(j,iter) = absdist(i)
					atom_bouge_abs  = i
				endif
			enddo
			absdmax_current = absdmax(j,iter)*1.d8
			write(*,*) ' absdmax(1,0)   ',absdmax(j,iter)*1.d8
			write(*,*) ' atom_bouge_abs ',atom_bouge_abs

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! INITIALISATION !!!!!!!!!!!!!!!!!!!!!!
			write(*,*) 'conditions initiales pour traj de reference avc distrib stoch à temperature T=',temperature
			write(*,*) 'temps equilibrage', Tequilib
			write (*,'("Equilibrage on ..................:",i5)') nint(Tequilib/dt) 
			!!!!!!!!!!!!!!!!!!!!!!!! EQUILIBRAGE initial (STOCH DYN)!!!!!!!!!!!!!!!!!!

			do it_langevin = 0,int(Tequilib/dt)
				call langevin(dt,temperature, rga)
				absdmax(j,:)=-9999.0
				atom_bouge_abs=0
				iter=0
				do i=1,N
					absdist(i)=sqrt((xp(1,i)-qref(1,i))**2+(xp(2,i)-qref(2,i))**2+(xp(3,i)-qref(3,i))**2)
					if (absdist(i).gt.absdmax(j,iter)) then
						absdmax(j,iter) = absdist(i)
						atom_bouge_abs  = i
					endif
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
			enddo
			absdmax_current = absdmax(j,iter)*1.d8
			write(*,*) ' absdmax(1,0)   ',absdmax(j,iter)*1.d8
			write(*,*) ' atom_bouge_abs ',atom_bouge_abs


			write(*,*) 'point de depart traj de reference deterministe' 
			! NbClones the index number of bias
			! The starting-point in the bias series

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
				call lanczos(N,maxvec,q1s2,new_projection,Path(j,iter)%project) !!!positions avant propagation
				Path(j,iter)%eigenvalue = eigenvalue
				write(*,*) 'eigenvalue ',eigenvalue
				new_projection=.false.  
				if (eigenvalue.lt.0.0) then 
					Path(j,iter)%Lyap =  asin(dt*sqrt(-eigenvalue)/2.d0)*2.d0 
				else 
					Path(j,iter)%Lyap = 0.d0  ! 
				endif
			enddo
		
		enddo
	endif  ! continue_sundae.ne.2 

	write(*,*) 'initialisation faite: go with Lanczos'
	write(*,*) 'itab', itab
	write(*,*) 'itetabvois', itetabvois
	write(*,*) 'ltabvois', ltabvois
  
  
  
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	if (continue_sundae.eq.2) then
		
		do j = 0, Nbclones_mbar
		
			write(*,*) 'posfinal ',posfinal
			open(unit=27, file=posfinal, status='old')
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
			enddo
			close(27)
			Path(j,totiter) =  Path(j,totiter-1) 
			oldLyap(j)=sum(Path(j,0:totiter-1)%Lyap)/real(totiter)
			write(*,*) 'oldLyap(',j,') = ', oldLyap(j)


			q=Path(j,0)%q
			do i=1,3
				write(*,*) 'i       ',i
				write(*,*) 'bary    ',sum(q(i,1:im))/dble(im)
				write(*,*) 'bary ref',sum(qref(i,1:im))/dble(im)
			enddo

			rang=0
			
			do i=0,totiter

				q    = Path(j,i)%q
				prot = Path(j,i)%p

				pav(1)=sum(prot(1,1:im))/dble(im)
				pav(2)=sum(prot(2,1:im))/dble(im)
				pav(3)=sum(prot(3,1:im))/dble(im)

				prot(1,1:im)= prot(1,1:im)-pav(1)
				prot(2,1:im)= prot(2,1:im)-pav(2)
				prot(3,1:im)= prot(3,1:im)-pav(3)

				pav(1)=sum(prot(1,1:im))/dble(im)
				pav(2)=sum(prot(2,1:im))/dble(im)
				pav(3)=sum(prot(3,1:im))/dble(im)

				call control_angular_momenta(prot,q)

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
				absdist(i)=sqrt(sum(pav(1:3)))
				if (absdist(i).gt.absdmax(j,iter)) then
					absdmax(j,iter) = absdist(i)
					atom_bouge_abs  = i
				endif
			enddo
			absdmax_current = absdmax(j,iter)*1.d8
			write(*,*) ' absdmax(1,0)   ',absdmax(j,iter)*1.d8
			write(*,*) ' atom_bouge_abs ',atom_bouge_abs

		enddo

	endif ! block with continue_sundae == 2
	

end subroutine LyapLanczos_equilibrage


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine LyapLanczos_allocate

	use tab_imm_m

	implicit none

	open(unit=30, file=kappaF,  action='write', status='replace')
	open(unit=31, file=kappaFd,  action='write', status='replace')
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	open(unit=111, file=data_abf, action='write', status='replace')
	open(unit=24, file=data_mbar, action='write', status='replace')
	open(unit=244, file=dada_mbar, action='write', status='replace')

	write(244,'(i,i,i)') Nbclones_mbar+1,Nbclones_mbar+1,Totalmcmoves


	nmax=NbClones_mbar

	write(*,*) ' nmax ', nmax
	 
	! Array where the weight of every clone is stored
	allocate (acc(0:nmax))
	allocate (ener0(0:nmax))
	allocate (enerK(0:nmax))
	allocate(hamilt(0:nmax))
	! Array which contains the numbers from 1 to nmax

	ltot=int(TotalTime/real(dt))

	write(*,*)'totiter', totiter
	gamma_sundae=gamma_sundae/(dt)
	write(*,*)'gamma_sundae' , gamma_sundae, gamma_sundae*dt,ltot
	waste_recycling=.true.
	!waste_recycling=.false.
	write(*,*) 'waste_recycling ?',  waste_recycling

	allocate(absdmax(0:nmax,0:2*totiter+10))
	allocate(S(0:nmax,0:2*totiter+10))
	allocate(u_kln(0:nmax,0:nmax,0:Totalmcmoves))
	allocate(ustd_kln(0:nmax,0:nmax,0:Totalmcmoves))
	allocate(umoy_kln(0:nmax,0:nmax,0:Totalmcmoves))
	allocate(u2moy_kln(0:nmax,0:nmax,0:Totalmcmoves))

	allocate(old_projection(3*N))
	allocate(first_projection(3*N))
	allocate(old_before_sc_projection(3*N))

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


	absdmax(:,:)=0
	hamilt(:)=0
	mcmoves=0
	it_art=0
	e=0.000001
	eigenvalue=0
	acc(:)=0

	if (maxvec.lt.4) maxvec=4
	write(*,*) 'The maxvec in lanczos is set to .........:', maxvec



	allocate(Path(0:nmax,0:totiter+10)) 
	allocate(Pshoot(0:nmax,0:totiter+10))
	allocate(Pshift(0:nmax,0:2*totiter+10))

	allocate(oldLyap(0:nmax))
	allocate(triallyap(0:nmax))
	allocate(rapport(0:nmax))
	allocate(absdist(0:N))
	allocate(Psel(0:nbclones_mbar,0:totiter))
	allocate(alpha_bias(0:Nbclones_mbar))!!EQUILIBRAGE
	allocate(q(1:3,1:im))
	allocate(p(1:3,1:im))
	rapport(:)=0
	iter=0
	Q=0
	!pi=4*atan(1._dpkind)
	
	if ((gamma_sundae*dt).le.100000) then
		rga = exp(-gamma_sundae*dt/two)
	else 
		rga=0
	endif
	
	sig(:,:) = sqrt(temperature*(one-rga**2))
	write(*,*) 'rga',gamma_sundae, gamma_sundae*dt, rga, sig(1,1)
	write(*,*) 'temp', temperature*erg2eV, erg2eV

	acc(:)=0
	!! initialisation paramètres de bias alpha pour reconstruction
	do j= 0, Nbclones_mbar
		alpha_bias(j)=(real(j*(alpha_max/real(Nbclones_mbar))))!*sqrt(9.270914743200000e-023)/300.0
		write(*,*)'bias', alpha_bias(j)  
	enddo



	xref(1:n)=xp(1,1:n)
	yref(1:n)=xp(2,1:n)
	zref(1:n)=xp(3,1:n)

	qref(1:3,1:im)=xp(1:3,1:im)

	rang = 1 
	
	! Allocation ABF
	
	delta = 5.d12
	delta_bin_theta_s2 = 0.5d0*real(alpha_max/real(Nbclones_mbar))
	theta_n = 50.0d12
	itheta_n = nint((theta_n)/(alpha_max)*real(Nbclones_mbar))
	N_extra = 10

	allocate(theta(-N_extra:nmax+N_extra))
	allocate(u_A(-N_extra:nmax+N_extra))
	allocate(A_n(-N_extra:nmax+N_extra))
	allocate(P_A(-N_extra:nmax+N_extra))
	allocate(sum_P_A(-N_extra:nmax+N_extra))
	allocate(A_prime(-N_extra:nmax+N_extra))
	allocate(A_prime_num(-N_extra:nmax+N_extra))
	allocate(O_moy_num(0:totiter))
	allocate(O_moy_estim(0:totiter))
	allocate(O_estim(0:totiter))
	allocate(histo_theta(0:nmax))
	
#if(PARASUN)
	allocate(MPI_histo_theta(0:nmax))
	allocate(MPI_P_A(-N_extra:nmax+N_extra))
	allocate(MPI_sum_P_A(-N_extra:nmax+N_extra))
	allocate(MPI_A_prime(-N_extra:nmax+N_extra))
	allocate(MPI_A_prime_num(-N_extra:nmax+N_extra))
	allocate(MPI_O_moy_num(0:totiter))
	allocate(MPI_O_moy_estim(0:totiter))
	allocate(MPI_O_estim(0:totiter))
#endif(PARASUN)
	
	theta = 0.d0
	do i=0,Nbclones_mbar
		theta(i) = alpha_bias(i)
	enddo



end subroutine LyapLanczos_allocate 



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!              OLD FUNCTIONS (MBAR)             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine essai_mbar

	use tab_imm_m

	implicit none
	
	real(double) :: genrand

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
		write(*,*) 'oldlyap(j)  = ', oldLyap(j)*300.0/sqrt(9.270914743200000e-023)

		!!!!!!!! on copie la trajectoire selectionnée avec le shifting
		if (absdmax(j,newtraj).le.h_A_max) then 
			Path(j,0:totiter-1) = Pshift(j,newtraj:totiter+newtraj-1) ! on a Path(j,0) = Pshift(j,0) pour newtraj=0
			oldLyap(j)       = SUM(Path(j,0:totiter-1)%Lyap)/real(totiter)
		else
			write(*,*) " Problème avec le shifting "
		endif

		write(*,*) 'newlyap(j)  = ', oldLyap(j)*300.0/sqrt(9.270914743200000e-023)

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
			!            modif de Manuel : j'ai déplacé l'initialization à zero avant la première utilisation. C'est plus sûr.

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


end subroutine essai_mbar





end module sundae_module



