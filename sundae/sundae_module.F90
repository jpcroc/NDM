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
#endif
	!-----------------------------------------------
	! test
	!   _c -> courant
	!   _s -> selection
	!   _d -> depart

	integer, parameter :: dpkind=selected_real_kind(13)
	real(double) :: t0,t1,elaps_1,tbuffer1,tbuffer2,t_lanczos
	integer, dimension(:), allocatable, save  :: ipovois
	
	
	integer :: continue_sundae, reprise_A, maxvec, totiter
	integer :: Totalmcmoves, depart_boucle_nbclones
	real(double) :: h_A_max, h_ba_max, h_ba_min, h_ba, h_ba_I, h_temp
	real(double) :: dt, TotalTime, gamma_sundae, Temperature
	real(double) :: alpha_max, teq, delta_x,ss, tequilib

	real(double) :: drift,dts2racinem,omegadts2

	character(len=128) :: recup
	character(len=128) :: srank
	character(len=128) :: sortie
	character(len=128) :: fnamtin
	character(len=128) :: posfinal
	character(len=128) :: data_abf
!m	character(len=128) :: data_mbar
!m	character(len=128) :: dada_mbar
!m	character(len=128) :: data_mbar_std
!m	character(len=128) :: moyennes_mbar
!m	character(len=128) :: moyennes_mbar_denom
!m	character(len=128) :: kappaF
!m	character(len=128) :: kappaFd

	real(double), parameter :: KtoERG=1.3791946308724831d-16

	real(double),dimension(:,:),allocatable,save  :: m_i
	real(double),dimension(:,:),allocatable,save  :: gau
	real(double),dimension(:,:),allocatable,save  :: rga_i


	real(double),dimension(:),allocatable,save    :: d2vois
	real(double),dimension(:,:),allocatable,save  :: xpvois
	real(double),dimension(:,:),allocatable,save  :: xpvoisini
	real(double),dimension(:),allocatable,save    :: xtransla


	integer :: iteration,isauvegarde
	integer :: icheck

	integer N
	PARAMETER(N=127)   !1023 
	integer nfenetre
	parameter(nfenetre=100)
	integer :: nl_iter

	! cosmin added:
	integer    :: it_langevin_deter=0, it_langevin=0,it_trajectory    
	
	!
	real (double), dimension(:,:), allocatable, save :: q,p 
	type Trajectoire
		real (double), dimension(1:3,1:N) :: q  ! vector of position
		real (double), dimension(1:3,1:N) :: p  ! vector of impulsion
		real (double), dimension(3*N) :: project
		real (double)  Lyap
		real (double)  eigenvalue
	endtype

	type (Trajectoire), dimension (:), allocatable :: Path
	type (Trajectoire), dimension (:), allocatable :: Pshoot
	type (Trajectoire), dimension (:), allocatable :: Pshift
	type (Trajectoire) :: Pcourant
	
	real (double), dimension(1:3,1:N):: qref,qrot,prot
	
	integer :: i,j,l,k,iter, atom_bouge_abs
	real (double), dimension(1:3,1:N):: q1s2
	logical :: new_projection
	real(double) :: pav(3)
	real (double) :: oldLyap 
	real (double) ,dimension(:),   allocatable:: absdmax
	
	
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

	real (double) :: ener0 
	real (double) :: enerK
	real (double) :: triallyap
	real (double) :: rapport
	real (double) :: hamilt
	real (double) ,dimension(:), allocatable:: absdist
	real (double) ,dimension(:), allocatable:: Psel
	real (double) ,dimension(:), allocatable:: S

	integer :: acc
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
	real(double) :: Lyap_max = 20.d0
	integer :: iLyap
	real(double) :: somme
	real(double) :: Lyap, gmax, gmin, P_n, pi_n, pi_n1, sum_A

      	real(double),dimension(:),  allocatable,save  :: a_sto
	real(double), dimension(:), allocatable, save :: theta
	real(double), dimension(:), allocatable, save :: u_A
	real(double), dimension(:), allocatable, save :: A_n
	real(double), dimension(:), allocatable, save :: P_A
	real(double), dimension(:), allocatable, save :: sum_P_A
	real(double), dimension(:), allocatable, save :: A_prime
	real(double), dimension(:), allocatable, save :: A_prime_num
	real(double), dimension(:), allocatable, save :: L2_moy_num
	real(double), dimension(:,:,:), allocatable, save :: O_moy_num
	real(double), dimension(:,:,:), allocatable, save :: O_moy_estim
	real(double), dimension(:,:,:), allocatable, save :: O_estim
	real(double) ,dimension(:), allocatable, save :: histo_theta	
	real(double), dimension(:), allocatable, save :: histo_Lyap
	real(double), dimension(:), allocatable, save :: theta_traj
	real(double), dimension(:), allocatable, save :: Lyap_traj
	real(double), dimension(:), allocatable, save :: ha_hb_traj

	! Complements MPI
	integer :: rank = 0, numproc = 1
#if(PARASUN)
	integer :: ierror
	real(double) ,dimension(:), allocatable, save :: MPI_histo_theta
	real(double) ,dimension(:), allocatable, save :: MPI_histo_Lyap
	real(double) ,dimension(:), allocatable, save :: MPI_P_A
	real(double) ,dimension(:), allocatable, save :: MPI_sum_P_A 
	real(double) ,dimension(:), allocatable, save :: MPI_L2
	real(double) ,dimension(:), allocatable, save :: MPI_L2_moy_num
	real(double) ,dimension(:), allocatable, save :: MPI_A_prime 
	real(double) ,dimension(:), allocatable, save :: MPI_A_prime_num
	real(double), dimension(:,:,:), allocatable, save :: MPI_O_moy_num
	real(double), dimension(:,:,:), allocatable, save :: MPI_O_moy_estim
#endif

		
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
	
	real(double) :: genrand,asto

	j = itheta_n
	!!write(*,*) 'indice theta', itheta_n
	!!!!! initialization backward

	! modif manuel

	Pshoot(ix) = Path(ix)

	eigenvalue = Path(ix)%eigenvalue

	temp = SUM((Path(ix)%p(1,1:N)**2+Path(ix)%p(2,1:N)**2+Path(ix)%p(3,1:N)**2))/m_i(1,1)/3.000/dble(im)/KtoERG ! *2/3 /2
	!!write(*,*) 'temp cin et temperature cible',temp,temperature/KtoERG

	!!!! shooting  time "a la Stoltz"

	p(1:3,1:im) = Pshoot(ix)%p(1:3,1:im)
	q(1:3,1:im) = Pshoot(ix)%q(1:3,1:im)

        asto = a_sto(itheta_n)

	call OU_control(p,q,asto,ss)  ! d'amplitude a_sto 

	Pshoot(ix)%p(1:3,1:N)= p(1:3,1:N)

	temp = SUM((Pshoot(ix)%p(1,1:N)**2+Pshoot(ix)%p(2,1:N)**2+Pshoot(ix)%p(3,1:N)**2))/m_i(1,1)/dble(3*im-6)/KtoERG ! *2/3 /2
	!!write(*,*) 'temp cin apres perturb et temp cible ',temp,temperature/KtoERG
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! propagation BACKWARD a partir du point de shooting (sans lyapunov)

	!!write(*,*) 'shooting ..ix->0...........bw ',ix,' -> ', 0
	new_projection=.false.
	it_trajectory=0
	Pcourant =  Pshoot(ix)

	if (ix.ge.0) then 
		Pcourant%project     = Path(ix)%project
		Pcourant%eigenvalue  = Path(ix)%eigenvalue
	endif
	lanczos_iter=0
	do iterbw=ix,1,-1 ! on stocke  en iterbw-1 la valeur propre calculée à iterbw-1/2
		call mapping_P_Verlet(Pcourant%q,Pcourant%p,-dt,N,q1s2) !! on met -dt pour le backward !!!
		call lanczos(N,maxvec,q1s2,new_projection,Pcourant%project)
		lanczos_iter=lanczos_iter+nl_iter
		Pcourant%eigenvalue = eigenvalue 
		if (eigenvalue.lt.0.0) then
                        omegadts2 = dts2racinem*sqrt(-eigenvalue)
			Pcourant%Lyap=asinh(omegadts2)*2.d0 !
		else 
			Pcourant%Lyap=0.d0
		endif
		Pshoot(iterbw-1)=Pcourant
		it_trajectory = it_trajectory + 1
		new_projection = .false.
	enddo  
	!!!!!!!!!!!!!!!!!!!!!  la trajectoire bw est arrivee au temps t=0
	
	!!!!!!!!!!!!           trajectoire forward        

	!!write(*,*) 'shooting ..ix->totiter-1...fw ',ix,' -> ', totiter-1
	Pcourant       = Pshoot(ix)
	new_projection =.false. ! Question Manuel: on doit recalculer?
	it_trajectory=0
	eigenvalue = Pshoot(ix)%eigenvalue

	! Manuel : avec mapping_P_Verlet on calcule les forces en iterfw+1/2 et on stocke positions et moments en iterfw+1 
	!                 mais on stocke le lyapunov avant en iterfw 
	do iterfw = ix,totiter-1  ! on recalcule en ix car la position en ix+1/2 est modifiée
		!    modif Manuel (texte+commentaires) 
		call mapping_P_Verlet(Pcourant%q,Pcourant%p,dt,N,q1s2)
		call lanczos(N,maxvec,q1s2,new_projection,Pcourant%project)
		lanczos_iter=lanczos_iter+nl_iter
		Pcourant%eigenvalue = eigenvalue
		if (eigenvalue.le.0.0) then 
!m	 	Pcourant%Lyap = asin(dt*sqrt(-eigenvalue)/2.d0)*2.d0  ! on stocke avant !!
                omegadts2 = dts2racinem*sqrt(-eigenvalue)
	        Pcourant%Lyap=asinh(omegadts2)*2.d0 !
!
		else 
			Pcourant%Lyap = 0.d0      ! on stocke toujours le lyapunov avant !!!
		endif
		Pshoot(iterfw)%Project=Pcourant%Project
		Pshoot(iterfw)%eigenvalue=Pcourant%eigenvalue
		Pshoot(iterfw)%Lyap=Pcourant%Lyap
		Pshoot(iterfw+1)   = Pcourant
		it_trajectory  = it_trajectory + 1
		new_projection = .false.
	enddo
	!!write(*,*) 'Lanczos iterations...:',lanczos_iter
	!!write(*,*) '           forces....:',lanczos_iter*maxvec*2

	triallyap=SUM(Pshoot(0:totiter-1)%Lyap)       !m    /real(totiter)
	!!write(*,*) 'oldLyap   = ', oldLyap*300.0/sqrt(9.270914743200000e-023)
	!!write(*,*) 'triallyap = ', triallyap*300.0/sqrt(9.270914743200000e-023)
	
	!!!!!!! fonction heaviside pour maintenir le point initial de la trajectoire dans 
	!!!!!!! le bassin de depart: calcul dist max pour le point x(0) de la trajectoire
	absdmax(iterbw)=-9999.0
	atom_bouge_abs=0
	do i=1,N
		absdist(i)=sqrt((Pshoot(iterbw)%q(1,i)-qref(1,i))**2 + &
		(Pshoot(iterbw)%q(2,i)-qref(2,i))**2 + &
		(Pshoot(iterbw)%q(3,i)-qref(3,i))**2)
		if (absdist(i).gt.absdmax(iterbw)) then
			absdmax(iterbw) = absdist(i)
			atom_bouge_abs    = i
		endif
	enddo
	absdmax(iterbw)      = absdmax(iterbw)*1.0d8
	!!write(*,*) 'iterbw absdmax ',iterbw,absdmax(iterbw)
	!!!!!!!!!!!! ATTENTION AUX RELATIONS DU BILAN DETAILLE

	if (absdmax(iterbw).lt.h_A_max) then
		rapport=exp(theta(j)*(triallyap-oldLyap))
	else
		rapport=0.d0
	endif

	write (*,*) ' itrialLyap oldLyap = ', triallyap , oldLyap 

	ranf=genrand()
	mcconf=min(1.d0,rapport)

	if (ranf.lt.mcconf) then
		acc=acc+1
		Path(0:totiter) = Pshoot(0:totiter)
		oldLyap=triallyap

		write(*,*) 'traj. acceptée',mcmoves,'rapport =',rapport
		absdmax_current = absdmax(iterbw)  
		!!write(*,*) 'absdmax_current', absdmax_current
	else
		write(*,*) 'traj. refusée ',mcmoves,'rapport =',rapport
	endif
	write (*,*) 'Shooting ....taux d''acceptation',real(acc)/real(mcmoves) 
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
	write(*,'(" Shifting............................:",3i6)') mcmoves,totiter,pix

	write(*,'(" Shift.rec...pix+1->pix+totinter....:",i6," to ",i6)') pix+1,pix+totiter
	! we recycle the pix,pix+1, pix+2, ...., pix+totiter because there are totiter+1 points 
	do k=0,totiter
		Pshift(pix+k) = Path(k)
	enddo

	write(*,'(" shift.bw ......pix->0..............:",i6," to 0")') pix
	! we calculate the pix-1,....,0
	new_projection=.false.
	it_trajectory=0

	Pcourant   = Pshift(pix)
	eigenvalue = Pcourant%eigenvalue

	do iterbw=pix,1,-1
		call mapping_P_Verlet(Pcourant%q,Pcourant%p,-dt,N,q1s2)  ! -dt car backward en position verlet
		call lanczos(N,maxvec,q1s2,new_projection,Pcourant%project)

		Pcourant%eigenvalue = eigenvalue

		if (eigenvalue.lt.0.0) then
!m		 Pcourant%Lyap= asin(dt*sqrt(-eigenvalue)/2.d0)*2.d0 ! on stocke Lyap en 1/2 avant
                 omegadts2 = dts2racinem*sqrt(-eigenvalue)
	         Pcourant%Lyap=asinh(omegadts2)*2.d0 !
		else 
		 Pcourant%Lyap=0.d0  !! autoval max del Lanczos
		endif
		Pshift(iterbw-1)=Pcourant
		it_trajectory=it_trajectory+1 
		new_projection=.false.
	enddo

	write(*,'(" shift.fw..pix+totiter->2*totinter..:" ,i6, " to ", i6)') pix+totiter,2*totiter
	! we calculate the pix+toiter+1, pix+2, ...., 2*toiter
	! Check: there are two points 2*toiter-1 and 2*toiter which
	!         are computed for nothing. Never used after !

	timefsh=dt*real(totiter+pix)
	new_projection=.true.
	iterfw=nint(timefsh/dt)
	it_trajectory=0

	Pcourant = Pshift(totiter+pix)
	eigenvalue = Pcourant%eigenvalue

	do iterfw=totiter+pix,2*totiter-1
		call mapping_P_Verlet(Pcourant%q,Pcourant%p,dt,N,q1s2)
		call lanczos(N,maxvec,q1s2,new_projection,Pcourant%project)!!!!!
		Pcourant%eigenvalue = eigenvalue 

		if (eigenvalue.le.0.d0) then 
!m			Pcourant%Lyap=asin(dt*sqrt(-eigenvalue)/2.d0)*2.d0  ! Lyapunov en 1/2 stocké avant
                 omegadts2 = dts2racinem*sqrt(-eigenvalue)
	         Pcourant%Lyap=asinh(omegadts2)*2.d0 ! 
		else 
		 Pcourant%Lyap=0.d0 ! Pshift(iterfw-1)%Lyap ! Lyapunov en 1/2 stocké avant 
		endif

		Pshift(iterfw)%Lyap       = Pcourant%Lyap
		Pshift(iterfw)%Project    = Pcourant%Project
		Pshift(iterfw)%eigenvalue = eigenvalue

		Pshift(iterfw+1)=Pcourant

		it_trajectory=it_trajectory+1 
		new_projection=.false.
	enddo

	hamilt=ener0

	S(:)=0
	do l=0,totiter  ! sommation sur les 1+L "path proposals" possibles
		do a=l,l+totiter-1  ! dans une action S il y a bien totiter=L Lyap stockes de 0 à totiter-1 
			S(l) = S(l)+Pshift(a)%Lyap
		enddo
	enddo
!m  	S(:)=S(:)/(real(totiter))

	!!!!!!!!!!!! calcul de la fonction echelon pour le point de départ 
	!!!!!!!!!!!! de la trajectoire N sui 2N passi della traj shiftata

	absdmax(:)=-9999.0
	atom_bouge_abs=0

	do l=0,2*totiter ! il y a 1+L (1+totiter) proposals + totiter positions décalée en totiter
		do i=1,N
			absdist(i)=sqrt((Pshift(l)%q(1,i)-xref(i))**2 + & 
				 (Pshift(l)%q(2,i)-yref(i))**2 + & 
				 (Pshift(l)%q(3,i)-zref(i))**2) 
			if (absdist(i).gt.absdmax(l)) then 
				absdmax(l)=absdist(i) 
			endif 

			! Manuel : les coordonnées dans Pshift(l) sont celles avant l'application du mapping P_Verlet dans le sens forward
			! on a donc maintenant correspondance entre Pshift(l) et absdmax(l) contrairement à auparavant

		enddo
	enddo
	absdmax(:)=absdmax(:)*1.0e8
	!!write(*,*) 'absdmax : '
	!!write(*,*) absdmax(0) 
	
	
	z=0.d0
	do l=0,totiter ! boucles sur les chemins proposés possibles
		! Here is P_sel from the paper
		if (absdmax(l).lt.h_A_max) then  ! on a correspondance maintenant entre absdmax et Pshift
			Psel(l) = exp(theta(j)*S(l))
		else
			Psel(l)= 0.d0
		endif
		z = z +  Psel(l)
	enddo

	Psel(0:totiter)  =  Psel(0:totiter)/z

	! selection de la trajectoire indiciée newtraj

	xalea    = genrand()
	xcumul=0.d0

	do k=0,totiter !  boucle sur les "proposals"
		xcumul=xcumul+Psel(k)
		if (xalea.lt.xcumul) goto 234
	enddo
	234 continue
	newtraj=k
	!!write(*,*) 'xalea, Psel, poids cumulé et newtraj = ',xalea,Psel(k), xcumul ,newtraj
	!!write(*,*) 'absdmax(newtraj)', absdmax(newtraj)
	!!write(*,*) 'oldLyap  = ', oldLyap*300.0/sqrt(9.270914743200000e-023)

	!!!!!!!! on copie la trajectoire selectionnée avec le shifting
	if (absdmax(newtraj).le.h_A_max) then 
		Path(0:totiter-1) = Pshift(newtraj:totiter+newtraj-1) ! on a Path(0) = Pshift(0) pour newtraj=0
		oldLyap       = SUM(Path(0:totiter-1)%Lyap)       !m /real(totiter)
	else
		!!write(*,*) " Problème avec le shifting "
	endif

	!!write(*,*) 'newlyap = ', oldLyap*300.0/sqrt(9.270914743200000e-023)	
	Lyap_traj(mcmoves) = oldLyap  !  *300.0/sqrt(9.270914743200000e-023)	

end subroutine LyapLanczos_shifting


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine LyapLanczos_ABF

	use lanczos_defs
	use tab_imm_m

	implicit none

	real(double) :: genrand,theta_o
	integer :: grid

	! Proposition du nouveau theta 
	xalea = genrand()
	grid = nint(10*(xalea-0.5d0))
	theta_temp = theta_n + delta*real(grid)
	theta_o    = theta_n 

        write(*,*) ' theta_n ' , theta_n
	itheta_n = nint(theta_n/alpha_max*real(nmax))
	itheta_temp = nint((theta_temp)/(alpha_max)*real(nmax))
	
	! Acceptation/rejet 
	if ( (itheta_temp.ge.0) .and. (itheta_temp.le.nmax) ) then
		AR_MH = min( 1.0d0 , exp((theta_temp-theta_n)*oldLyap - A_n(itheta_temp) + A_n(itheta_n)))
	endif
	ranf = genrand()
	if ( (ranf.lt.AR_MH) .and. (theta_temp.gt.0.0d0) .and. (theta_temp.lt.alpha_max) ) then
		theta_n = theta_temp
		itheta_n = itheta_temp
	endif

	write(*,*) 'theta n / o ' ,theta_n , theta_o
	theta_traj(mcmoves) =  theta_n

	! Calcul du nouveau biais
	Lyap          = oldLyap
	u_A(:)        = Lyap*theta(:) - A_n(:)     
	gmax          = maxval(u_A)
	u_A(:)        = u_A(:) - gmax
	P_A(:)        = exp(u_A(:))
	somme         = sum(P_A(-N_extra:nmax+N_extra))
	P_A           = P_A/somme

#if(PARASUN)
	! Mise en commun parallele ---- 1/2 !	
	call MPI_REDUCE(P_A,MPI_P_A,nmax+2*N_extra+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, MPI_COMM_WORLD,ierror)

	if (rank.eq.0) then
		MPI_P_A = MPI_P_A/dble(numproc)
		MPI_histo_theta(0:nmax) = MPI_histo_theta(0:nmax) + MPI_P_A(0:nmax)
	end if

	call MPI_BCAST(MPI_P_A,nmax+2*N_extra+1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
	! Mise en commun parallele ---- 1/2 !
#else
	histo_theta(0:nmax) = histo_theta(0:nmax) + P_A(0:nmax)
#endif
	
	! Pre-calcul des moyennes ergodiques
	do theta_tilde=-N_extra,nmax+N_extra
		P_n        = P_A(theta_tilde)
		pi_n       = sum_P_A(theta_tilde)
		pi_n1      = pi_n + P_n

		if (pi_n1.eq.0.d0) then 
			pi_n1 = 1.0d-10
		endif
		
		sum_P_A(theta_tilde)       =   pi_n1

		A_prime_num(theta_tilde)   =   P_n * Lyap + A_prime_num(theta_tilde)
		L2_moy_num(theta_tilde)    =   P_n * Lyap*Lyap + L2_moy_num(theta_tilde)
		!A_prime(theta_tilde)       =   A_prime_num(theta_tilde) / pi_n1
	enddo
	
	do k = 0,totiter 
		O_moy_num(k,1:3,1:4)    =   P_A(0) * O_estim(k,1:3,1:4) + O_moy_num(k,1:3,1:4)
	enddo


#if(PARASUN)
	! Mise en commun parallele ---- 2/2 !
	call MPI_REDUCE(A_prime_num,MPI_A_prime_num,nmax+2*N_extra+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, MPI_COMM_WORLD,ierror)
	call MPI_REDUCE(O_moy_num,MPI_O_moy_num,4*3*(totiter+1),MPI_DOUBLE_PRECISION,MPI_SUM,0, MPI_COMM_WORLD,ierror)
	call MPI_REDUCE(L2_moy_num,MPI_L2_moy_num,nmax+2*N_extra+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, MPI_COMM_WORLD,ierror)
	call MPI_REDUCE(sum_P_A,MPI_sum_P_A,nmax+2*N_extra+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, MPI_COMM_WORLD,ierror)

	MPI_A_prime = MPI_A_prime_num / (MPI_sum_P_A + 1.d-6)
	MPI_L2  = MPI_L2_moy_num / (MPI_sum_P_A + 1.d-6)
	
	do k = 0,totiter 
		O_moy_estim(k,1:3,1:4)  =   MPI_O_moy_num(k,1:3,1:4) / (MPI_sum_P_A(0) + 1.d-6)
	enddo
	
	if (N_extra.gt.0) then 
		MPI_A_prime(-N_extra:-1) = 0.d0
		MPI_A_prime(nmax+1:nmax+N_extra) = 0.d0
		MPI_L2(-N_extra:-1) = 0.d0
		MPI_L2(nmax+1:nmax+N_extra) = 0.d0
	endif

	call MPI_BCAST(MPI_A_prime,nmax+2*N_extra+1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
	call MPI_BCAST(MPI_L2,nmax+2*N_extra+1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
	! Mise en commun parallele ---- 2/2 !

	! Calcul du biais A
	if ( (continue_sundae.ne.1) .or. (reprise_A.ne.1) ) then
		A_n = 0.d0
		A_n(-N_extra) = 0.d0
		do theta_tilde  = -N_extra+1,nmax + N_extra
			A_n(theta_tilde) = A_n(theta_tilde-1)  + (MPI_A_prime(theta_tilde-1) + MPI_A_prime(theta_tilde))*delta_bin_theta_s2
		enddo
	endif
#else
	do theta_tilde=-N_extra,nmax+N_extra
		A_prime(theta_tilde)  =  A_prime_num(theta_tilde) / (pi_n1 + 1.d-6)
	enddo
	do k = 0,totiter 
		O_moy_estim(k,1:3,1:4)  =  O_moy_num(k,1:3,1:4) / (sum_P_A(0) + 1.d-6)
	enddo
	A_n = 0.d0
	A_n(-N_extra) = 0.d0
	do theta_tilde  = -N_extra+1,nmax + N_extra
		A_n(theta_tilde) = A_n(theta_tilde-1)  + (A_prime(theta_tilde-1) + A_prime(theta_tilde))*delta_bin_theta_s2
	enddo
#endif

	! normalisation du generateur biaisant
	if ( (continue_sundae.ne.1) .or. (reprise_A.ne.1) ) then
		gmin=minval(A_n)
		A_n = A_n - gmin
		sum_A = log(sum(exp(-A_n(-N_extra:nmax+N_extra))))
		A_n = A_n + sum_A
	endif
		
	!!write(*,*) 'theta_n  = ', theta_n*sqrt(9.270914743200000e-023)/300.0


end subroutine LyapLanczos_ABF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine LyapLanczos_output

	use lanczos_defs
	use tab_imm_m

	implicit none
	real(double) theta_f, A_prime_f, L2_f, somme
	character(len=128) :: fic1, fic2, fic3, fic4, fic5 , fic6, smcmoves
	
	! ----------- Sortie des observables ----------- !
#if(PARASUN)
	! Mise en commun parallele !	
	call MPI_REDUCE(histo_Lyap,MPI_histo_Lyap,nmax+2*N_extra+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, MPI_COMM_WORLD,ierror)

	MPI_histo_Lyap = MPI_histo_Lyap/dble(numproc)

	call MPI_BCAST(MPI_histo_Lyap,nmax+2*N_extra+1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
	! Mise en commun parallele !
#endif
	
	if (rank.eq.0) then
		write( smcmoves, '(i6)' )  mcmoves
		fic1 = trim(data_abf)//'_histo'
		fic2 = trim(data_abf)//'_obs_1_'//trim(adjustl(smcmoves))
		fic3 = trim(data_abf)//'_obs_2_'//trim(adjustl(smcmoves))
		open(unit=111, file=data_abf, action='write', status='replace')
		open(unit=1110, file=fic1, action='write', status='replace')
		open(unit=1111, file=fic2, action='write', status='replace')
		open(unit=1112, file=fic3, action='write', status='replace')

		fic4 = trim(sortie)//'.pathinit'
                fic5 = trim(sortie)//'.eigen'
		fic6 = trim(sortie)//'.pmf'
                write(*,*) 'fic4 fic5 fic6 ',fic4 , fic5 , fic6
		open(unit=772, file=fic4, action='write', status='replace')
		open(unit=773, file=fic5, action='write', status='replace')
		open(unit=774, file=fic6, action='write', status='replace')

	endif
#if(PARASUN) 
	if (rank.eq.0) then
		somme = sum(MPI_histo_theta)
		MPI_histo_theta = MPI_histo_theta/somme*real(nmax)
		somme = sum(MPI_histo_Lyap)
		MPI_histo_Lyap = MPI_histo_Lyap/somme*real(nmax)
		do j=-N_extra,nmax+N_extra
			! Renormalisation des variables pour les sorties fichier
			theta_f = theta(j)           !  *sqrt(9.270914743200000e-023)/300.0
			A_prime_f = MPI_A_prime(j)   !  *300/sqrt(9.270914743200000e-023)
			L2_f = MPI_L2(j)             !   *300/sqrt(9.270914743200000e-023)*300/sqrt(9.270914743200000e-023)
			! Sortie fichier : theta, histo_theta, Lyap_moyen (A_prime)
			write(111,'(3(E15.6E3))') theta_f, A_prime_f, L2_f
			write(1110,'(2(E15.6E3))') MPI_histo_theta(j), MPI_histo_Lyap(j)
		enddo
		do j=0,totiter
			write(1111,'(5(E15.6E3))') O_moy_estim(j,1,1), O_moy_estim(j,1,2), O_moy_estim(j,1,3), O_moy_estim(j,1,4), MPI_sum_P_A(0)
			write(1112,'(5(E15.6E3))') O_moy_estim(j,2,1), O_moy_estim(j,2,2), O_moy_estim(j,2,3), O_moy_estim(j,2,4), MPI_sum_P_A(0)
		enddo
	endif	
#else
	do j=-N_extra,nmax+N_extra
		! Renormalisation des variables pour les sorties fichier
		theta_f = theta(j)          !        *sqrt(9.270914743200000e-023)/300.0
		A_prime_f = A_prime(j)      !        *300/sqrt(9.270914743200000e-023)
		! Sortie fichier : theta, histo_theta, Lyap_moyen (A_prime)
		write(111,*) theta_f, histo_theta(j), A_prime_f, histo_Lyap(j)
	enddo
	do j=0,totiter
		write(1110,'(3(E15.6E3))') O_moy_estim(j,1,1), O_moy_estim(j,1,1), O_moy_estim(j,1,1)
	enddo
#endif
	if (rank.eq.0) then
		close(111)
		close(1110)
		close(1111)
		close(1112)
	endif
	! ----------- Sortie des observables ----------- !
	
	
	! ----------- Sortie fichier de recuperation ----------- !
	write( srank, '(i2)' )  rank
	recup = sortie(1:lenfnam)//trim(adjustl(srank))//'.in'
	if (continue_sundae.eq.1) then
		open(unit=112,file=recup,status='replace')		
	else if (continue_sundae.eq.2) then
		open(unit=112,file=recup,status='replace')
	endif
	

!!!!!!!!!!!!!!!!! nouveau format  

	write (772,*) Path(0)%q
	write (772,*)
	write (772,*) Path(0)%p
	write (772,*)
	write (772,*) theta_n 

	do i=0, totiter-1                      ! attention i=totiter doit etre pris en compte
		write (773,*) Path(i)%Lyap,Path(i)%eigenvalue
        enddo
	do i = -N_extra, nmax+N_extra 
		write (774,*) theta(i), A_n(i)
	enddo

	do i=0, totiter-1   
		write (775,*) Path(i)%project
		write (775,*)
	enddo


	close(772)
	close(773)
	close(774)
	close(775)

!!!!!!!!!!!!!!!!! fin nouveau format


	!!write(*,*) 'recup', recup
        	
!	do i=0, totiter-1                      ! attention i=totiter doit etre pris en compte
!		write (112,*) Path(i)%q
!		write (112,*)
!		write (112,*) Path(i)%p
!		write (112,*)
!		write (112,*) Path(i)%Lyap
!		write (112,*)
!		write (112,*) Path(i)%project
!		write (112,*)
!		write (112,*) Path(i)%eigenvalue
!		write (112,*)
!	enddo
!	do i = -N_extra, nmax+N_extra 
!		write (112,*) A_n(i)
!	enddo
!		write(112,*)
!		write(112,*) theta_n 
!	close(112)




	! ----------- Sortie fichier de recuperation ----------- !
	
	
	! ----------- Sortie fichier Lyapunov-trajectoire ----------- !
	write( srank, '(i2)' )  rank
	recup = sortie(1:lenfnam)//trim(adjustl(srank))//'.Lyapunov'
	open(unit=rank,file=recup,status='replace')		
		
	do k = 1, mcmoves
		write(rank,'(3(E15.6E3))') theta_traj(k),Lyap_traj(k), ha_hb_traj(k)
	enddo

	close(rank)
	! ----------- Sortie fichier Lyapunov-trajectoire ----------- !
	
	

end subroutine LyapLanczos_output


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine LyapLanczos_vac! (xp)

	use lanczos_defs
	use tab_imm_m

	implicit none

	real(double) :: genrand
	

#if(PARASUN) 
	
	call MPI_INIT(ierror)     !! On initialise MPI

	call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror) !! Attribue a rank le numero de processeur

	call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)  !! Attribue a numproc le nombre de processeur
	
	call init_random_seed(100000*(rank+1))
	
#else

	call init_random_seed(100000)

#endif


	!!!!!!!!!!!!!!!!!!!!!!      Phase d'initialisation       !!!!!!!!!!!!!!!!!!!!!!!!

	call read_sundae                         !!! Modif 21.05.14

    call LyapLanczos_allocate                !!! Modif 02.06.14

	call LyapLanczos_equilibrage             !!! Modif 02.06.14
	
	call LyapLanczos_init_tests              !!! Modif 23.05.14
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




	do mcmoves = 1 , Totalmcmoves                       !mcmoves = 1,M in the paper
	
		
		!write(*,'(" Proc number ", i6, "says hello")'), rank
	
	
		it_art=mcmoves

		ix=int((totiter)*genrand())   ! le cas ix = totiter ne doit pas etre possible car on stocke jusqu'à totiter. 
		!!write(*,*) 
		!!write(*,*) 
		write(*,'(" Itération..:",i6," sur un total de ",i6)') ,mcmoves,Totalmcmoves
		write(*,'(" Shooting......................:",3i6)') mcmoves,totiter,ix

		!!!! shooting deterministique
		! Check: if depart_boucle_nbclones,NbClones of if there are different WHY?


		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!
		
		call LyapLanczos_shooting              !!! Modif 02.06.14
	
	
		call LyapLanczos_shifting              !!! Modif 03.06.14
	
	
		call LyapLanczos_Obs				   !!! Modif 13.06.14
	
		call LyapLanczos_ABF				   !!! Modif 10.06.14
		
		write (*,*) isauvegarde
		if (0.eq.mod(mcmoves,isauvegarde)) then 
			call LyapLanczos_output
		endif
		
		!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	enddo !! mcmoves avec incrément + 1 clones


	call LyapLanczos_output

#if(PARASUN)
	call MPI_FINALIZE(ierror)	 !! Pour finir l'appel MPI
#endif

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!	do j=depart_boucle_nbclones, NbClones
!		!!write(*,*) 'tau',j,theta(j),real(acc)/real(Totalmcmoves)
!	enddo

	stop

end subroutine LyapLanczos_vac


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine LyapLanczos_Obs

	use lanczos_defs
	use tab_imm_m

	implicit none

	real(double) :: num

	!  calcul des fonctions indicatrices pour l'etat B qui correspond a la barriere dans le cas lacune 
	h_Fg(:)     = 0
	h_Fd(:)     = 0
	h_F(:)      = 0
	h_FI(:)     = 0
	h_dI(:)     = 0
	num = 0
	
	
	! Estimation simple sur la partie [0,L] de la trajectoire
	do k = 0,totiter 
		!passagge F-d et F-I
		if ((absdmax(newtraj+k).ge.h_ba_min) .and. (absdmax(newtraj+k).lt.h_ba) .and. (absdmax(newtraj).lt.h_A_max)) then
			h_Fg(k) = 1.0
		endif
		if ((absdmax(newtraj+k).ge.h_ba) .and. (absdmax(newtraj+k).lt.h_ba_max) .and. (absdmax(newtraj).lt.h_A_max)) then   
		!!!FCC -> DEFAULT FCC
			h_Fd(k) = 1.0
		endif
		if ((absdmax(newtraj+k).ge.h_ba_max) .and. (absdmax(newtraj+k).lt.h_ba_I) .and. (absdmax(newtraj).lt.h_A_max)) then 
		!!FCC -> FCC 
			h_F(k)  = 1.0
		endif
		if ((absdmax(newtraj+k).ge.h_ba_I) .and. (absdmax(newtraj).lt.h_A_max)) then 
			h_FI(k) = 1.0
		endif
	enddo !k
	
	do k = 0,totiter
		O_estim(k,1,1) = h_Fg(k)
		O_estim(k,1,2) = h_Fd(k)
		O_estim(k,1,3) = h_F(k)
		O_estim(k,1,4) = h_FI(k)
	enddo
	ha_hb_traj(mcmoves) = h_Fd(totiter)+h_F(totiter)+h_FI(totiter)
	
	! Estimation en moyenne glissante sur [0,2L] trajectoire shiftee
	h_Fg(:) = 0
	h_Fd(:) = 0
	h_F(:)  = 0
	h_FI(:) = 0
	do k = 0,totiter 
		do j = 0,totiter 
			if ((absdmax(j+k).ge.h_ba_min) .and. (absdmax(j+k).lt.h_ba) .and. (absdmax(j).lt.h_A_max)) then 
				h_Fg(k) = h_Fg(k) + 1.0
			endif
			if ((absdmax(j+k).ge.h_ba) .and. (absdmax(j+k).lt.h_ba_max) .and. (absdmax(j).lt.h_A_max)) then 
				h_Fd(k) = h_Fd(k) + 1.0
			endif
			if ((absdmax(j+k).ge.h_ba_max) .and. (absdmax(j+k).lt.h_ba_I) .and. (absdmax(j).lt.h_A_max)) then 
				h_F(k) = h_F(k) + 1.0
			endif
			if ((absdmax(j+k).ge.h_ba_I) .and. (absdmax(j).lt.h_A_max)) then 
				h_FI(k) = h_FI(k) + 1.0
			endif
			if (absdmax(j).lt.h_A_max) then 
				num = num + 1.0
			endif
		enddo !j
		if (num.eq.0) then
			O_estim(k,2,1) = 0.
			O_estim(k,2,2) = 0.
			O_estim(k,2,3) = 0.
			O_estim(k,2,4) = 0.
		else
			O_estim(k,2,1) = h_Fg(k)/num
			O_estim(k,2,2) = h_Fd(k)/num
			O_estim(k,2,3) = h_F(k)/num
			O_estim(k,2,4) = h_FI(k)/num
		endif
		num = 0
	enddo !k
	
	
	! Estimation en moyenne glissante avec Waste Recycling sur [0,2L] trajectoire shiftee
	h_Fd(:) = 0
	do k = 0,totiter 
		do j = 0,totiter 
			if ((absdmax(j+k).ge.h_ba) .and. (absdmax(j+k).lt.h_ba_max) .and. (absdmax(j).lt.h_A_max)) then 
				h_Fd(k) = h_Fd(k) + Psel(j) 
				!!!write(*,*) "*************************************************"
				!!!write(*,*) "************                         ************"
				!!!write(*,*) "************  passage dans l'etat B  ************"
				!!!write(*,*) "************                         ************"
				!!!write(*,*) "*************************************************"
			endif
		enddo !j
	enddo !k
	
	do k = 0,totiter
		O_estim(k,3,1) = h_Fd(k)
		O_estim(k,3,2) = h_Fd(k)
		O_estim(k,3,3) = h_Fd(k)
		O_estim(k,3,4) = h_Fd(k)
	enddo
	
	! -------- Histogramme du Lyap --------- !
	iLyap = nint(Lyap/Lyap_max*real(nmax))
	if ((iLyap.ge.0).and.(iLyap.le.nmax)) then
		histo_Lyap(iLyap) = histo_Lyap(iLyap) + 1.0d0
	endif

end subroutine LyapLanczos_Obs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine LyapLanczos_init_tests

	!use lanczos_defs
	use tab_imm_m
	implicit none
        real(double) :: asto
	

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
	!do j=nmax,nmax
	do j = itheta_n, itheta_n
	
		q= Path(0)%q
		p= Path(0)%p
	        asto = 0.d0
		temp=0.d0
		do i=1,20
		call OU_control(p,q,asto,ss)
		call mapping_P_Verlet(q,p,dt,N,q1s2)

		temp = temp + (sum(p(1,1:im)**2)+sum(p(2,1:im)**2)+sum(p(3,1:im)**2))/dble(3*im-6)/m_i(1,1)
		!!write(*,*) 'température ',temp/dble(i)/KtoERG,i,theta(j)
		enddo
	
	
		write (*,'("1st traj before shooting ........:",i6)') nint(TotalTime/dt)
		rang=1

		call cpu_time(t0)
		lanczos_iter=0

		do icheck=1,1
			Pcourant = Path(0)
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
                                 omegadts2 = dts2racinem*sqrt(-eigenvalue)
	                         Pcourant%Lyap=asinh(omegadts2)*2.d0 
				temp = temp +  asin(dt*sqrt(-eigenvalue)/2.d0)*2.d0/sqrt(9.270914743200000e-023)
				endif
				Path(iter+1)    = Pcourant
				Path(iter)%Lyap = Pcourant%Lyap
				new_projection = .false.
			enddo  !!!!!
!			oldLyap=sum(Path(0:totiter-1)%Lyap) !m /real(totiter)
			write(*,*) " oldLyap= ", oldLyap,temp !m /real(totiter)
		enddo
		pav(1)=sum(Path(totiter)%p(1,1:im))/dble(im)
		pav(2)=sum(Path(totiter)%p(2,1:im))/dble(im)
		pav(3)=sum(Path(totiter)%p(3,1:im))/dble(im)

		!!write(*,*) 'pav',pav

		call cpu_time(t1)
		elaps_1=t1-t0

		!!write(*,*) 'The first trajectory..:', elaps_1
		!!write(*,*) 'The Lanczos time .....:', t_lanczos
		!!write(*,*) 'The Propag time.......:', elaps_1-t_lanczos
		!!write(*,*) 'Lanczos interations...:', lanczos_iter
		!!write(*,*) '            forces....:', lanczos_iter*maxvec*2


		absdmax(:)=-9999.d0
		do iter=0,totiter-1
			atom_bouge_abs=0
			do i=1,N
				pav(1:3) = (Path(iter)%q(1:3,i)-qref(1:3,i))**2
				absdist(i)=sqrt(sum(pav(1:3)))
				if (absdist(i).gt.absdmax(iter)) then
					absdmax(iter) = absdist(i)
					atom_bouge_abs  = i
				endif
			enddo
		enddo

		temp=0.0000
		tempvar=temp

		do iterbw=0,totiter
			tempiter = SUM((Path(iterbw)%p(1,1:N)**2+Path(iterbw)%p(2,1:N)**2+Path(iterbw)%p(3,1:N)**2))/m_i(1,1)/3.000/dble(im-2)/KtoERG ! *2/3 /2
			temp=temp+tempiter
			tempvar=tempvar+tempiter**2
		enddo
		temp=temp/real(totiter+1)
		tempvar=sqrt(tempvar/real(totiter+1)-temp**2)
		!!write(*,*) 'temp cin traject, std et cible ',temp,tempvar,temperature/KtoERG
		! calcul de la position initiale 

		call cal_hamilton(Path(0)%q,Path(0)%p,N,ekin,epot)
		hamilt=(ekin+epot)/temperature
		ener0 = epot/temperature

		! vérification de la dérive 
		xbar(1) = SUM(Path(0)%q(1,1:N))/dble(N)
		xbar(2) = SUM(Path(0)%q(2,1:N))/dble(N)
		xbar(3) = SUM(Path(0)%q(3,1:N))/dble(N)

		!!write(*,*) ' centre de masse référence ', xbar(1:3)

		xbar(1) = SUM(qref(1,1:N))/dble(N)
		xbar(2) = SUM(qref(2,1:N))/dble(N)
		xbar(3) = SUM(qref(3,1:N))/dble(N)

		!!write(*,*) ' centre de masse pos 0     ', xbar(1:3)

		xbar(1) = SUM(Path(0)%p(1,1:N))/dble(N)
		xbar(2) = SUM(Path(0)%p(2,1:N))/dble(N)
		xbar(3) = SUM(Path(0)%p(3,1:N))/dble(N)

		!!write(*,*) ' moments translationnels ',   xbar(1:3)
		do i=1,3
			!!write(*,*) ' temperature translation ',  xbar(i)**2/m_i(1,1)/dble(im)/KtoERG*dble(N)
		enddo
		iterbw = 0
		tempiter = SUM(Path(iterbw)%p(1,1:N)**2+Path(iterbw)%p(2,1:N)**2+Path(iterbw)%p(3,1:N)**2)/m_i(1,1)/3.000/dble(im)/KtoERG ! *2/3 /2
		!!write(*,*) 'temp avant ',iterbw,tempiter
		enerK =tempiter
		do iterfw=0,totiter
			Path(iterfw)%p(1,1:N)=Path(iterfw)%p(1,1:N)-xbar(1)
			Path(iterfw)%p(2,1:N)=Path(iterfw)%p(2,1:N)-xbar(2)
			Path(iterfw)%p(3,1:N)=Path(iterfw)%p(3,1:N)-xbar(3)
		enddo

		tempiter = SUM(Path(iterbw)%p(1,1:N)**2+Path(iterbw)%p(2,1:N)**2+Path(iterbw)%p(3,1:N)**2)/m_i(1,1)/3.000/dble(im)/KtoERG ! *2/3 /2
		!!write(*,*) 'tempiter après ',iterbw,tempiter
		!!write(*,*) 'diff ',tempiter-enerK 
	
	enddo

end subroutine LyapLanczos_init_tests


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine LyapLanczos_equilibrage

	use tab_imm_m

	implicit none
	
	if ( (continue_sundae.ne.2) .and. (continue_sundae.ne.1)) then 
          write(*,*) 'pb car continue_sundae = ',continue_sundae
          stop	
		absdmax(:)=-9999.0
		atom_bouge_abs=0
		iter=0
		do i=1,N
			absdist(i)=sqrt((xp(1,i)-qref(1,i))**2+(xp(2,i)-qref(2,i))**2+(xp(3,i)-qref(3,i))**2)
			if (absdist(i).gt.absdmax(iter)) then
				absdmax(iter) = absdist(i)
				atom_bouge_abs  = i
			endif
		enddo
		absdmax_current = absdmax(iter)*1.d8
		!!write(*,*) ' absdmax(1,0)   ',absdmax(iter)*1.d8
		!!write(*,*) ' atom_bouge_abs ',atom_bouge_abs

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! INITIALISATION !!!!!!!!!!!!!!!!!!!!!!
		!!write(*,*) 'conditions initiales pour traj de reference avc distrib stoch à temperature T=',temperature
		!!write(*,*) 'temps equilibrage', Tequilib
		write (*,'("Equilibrage on ..................:",i6)') nint(Tequilib/dt) 
		!!!!!!!!!!!!!!!!!!!!!!!! EQUILIBRAGE initial (STOCH DYN)!!!!!!!!!!!!!!!!!!

		do it_langevin = 0,int(Tequilib/dt)
			call langevin(dt,temperature, rga)
			absdmax(:)=-9999.0
			atom_bouge_abs=0
			iter=0
			do i=1,N
				absdist(i)=sqrt((xp(1,i)-qref(1,i))**2+(xp(2,i)-qref(2,i))**2+(xp(3,i)-qref(3,i))**2)
				if (absdist(i).gt.absdmax(iter)) then
					absdmax(iter) = absdist(i)
					atom_bouge_abs  = i
				endif
			enddo
			absdmax_current = absdmax(iter)*1.d8
			!!write(*,*) ' absdmax(1,0)   ',absdmax(iter)*1.d8
			!!write(*,*) ' atom_bouge_abs ',atom_bouge_abs
		end do

		absdmax(:)=-9999.0
		atom_bouge_abs=0
		iter=0
		do i=1,N
			absdist(i)=sqrt((xp(1,i)-qref(1,i))**2+(xp(2,i)-qref(2,i))**2+(xp(3,i)-qref(3,i))**2)
			if (absdist(i).gt.absdmax(iter)) then
				absdmax(iter) = absdist(i)
				atom_bouge_abs  = i
			endif
		enddo
		absdmax_current = absdmax(iter)*1.d8
		!!write(*,*) ' absdmax(1,0)   ',absdmax(iter)*1.d8
		!!write(*,*) ' atom_bouge_abs ',atom_bouge_abs


		!!write(*,*) 'point de depart traj de reference deterministe' 
		! NbClones the index number of bias
		! The starting-point in the bias series

		Path(0)%q=xp(1:3,1:N)
		Path(0)%p=vp(1:3,1:N)*m_i(1:3,1:N)
		!!!!!!!!!!!!!! mettre a zero le reste avant ??
		Path(:)%Lyap = 0.d0

		iter=0
		it_trajectory=0
		new_projection=.true.
		it_art=1
		lanczos_iter=0
		eigenvalue = 0.0

		do icheck=1,10
			Pcourant=Path(iter)
			call mapping_P_Verlet(Pcourant%q,Pcourant%p,dt,N,q1s2)
			call lanczos(N,maxvec,q1s2,new_projection,Path(iter)%project) !!!positions avant propagation
			Path(iter)%eigenvalue = eigenvalue
			!!write(*,*) 'eigenvalue ',eigenvalue
			new_projection=.false.  
			if (eigenvalue.lt.0.0) then 
		!m	Path(iter)%Lyap =  asin(dt*sqrt(-eigenvalue)/2.d0)*2.d0
                         omegadts2 = dts2racinem*sqrt(-eigenvalue)
	                 Pcourant%Lyap=asinh(omegadts2)*2.d0 !
			else 
	  		 Path(iter)%Lyap = 0.d0  ! 
			endif
		enddo
		
	endif  ! continue_sundae.ne.2 

	!!write(*,*) 'initialisation faite: go with Lanczos'
	!!write(*,*) 'itab', itab
	!!write(*,*) 'itetabvois', itetabvois
	!!write(*,*) 'ltabvois', ltabvois
  
  
  
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	if (continue_sundae.eq.2) then
		
		
		!!write(*,*) 'posfinal ',posfinal
		open(unit=27, file=posfinal, status='old')
		!!write(*,*) 'j depart_boucle_nbclones NbClones  ', j ,depart_boucle_nbclones,NbClones
		
		do i=0, totiter-1                      ! attention i=totiter doit etre pris en compte
			read (27,*) qtemp(1:N,1)
			read (27,*) qtemp(1:N,2)
			read (27,*) qtemp(1:N,3)
			read (27,*) 
			read (27,*) qtemp(1:N,4)
			read (27,*) qtemp(1:N,5)
			read (27,*) qtemp(1:N,6)
			read (27,*)
			read (27,*) Path(i)%Lyap
			read (27,*)
			read (27,*) Path(i)%project

			Path(i)%q(1,1:N)=qtemp(1:N,1)
			Path(i)%q(2,1:N)=qtemp(1:N,2)
			Path(i)%q(3,1:N)=qtemp(1:N,3)
			Path(i)%p(1,1:N)=qtemp(1:N,4)
			Path(i)%p(2,1:N)=qtemp(1:N,5)
			Path(i)%p(3,1:N)=qtemp(1:N,6)
		enddo
		close(27)
		Path(totiter) =  Path(totiter-1) 
		oldLyap=sum(Path(0:totiter-1)%Lyap)      !m /real(totiter)
		write(*,*) 'oldLyap = ', oldLyap


		q=Path(0)%q
		do i=1,3
			!!write(*,*) 'i       ',i
			!!write(*,*) 'bary    ',sum(q(i,1:im))/dble(im)
			!!write(*,*) 'bary ref',sum(qref(i,1:im))/dble(im)
		enddo

		rang=0
		
		do i=0,totiter

			q    = Path(i)%q
			prot = Path(i)%p

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

			Path(i)%q=q

			call control_angular_momenta(prot,q)

			Path(i)%p=prot

		enddo

		absdmax(:)=-9999.0
		atom_bouge_abs=0
		iter=0
		do i=1,N
			pav(1:3)=(Path(iter)%q(1:3,i)-qref(1:3,i))**2
			absdist(i)=sqrt(sum(pav(1:3)))
			if (absdist(i).gt.absdmax(iter)) then
				absdmax(iter) = absdist(i)
				atom_bouge_abs  = i
			endif
		enddo
		absdmax_current = absdmax(iter)*1.d8
		!!write(*,*) ' absdmax(1,0)   ',absdmax(iter)*1.d8
		!!write(*,*) ' atom_bouge_abs ',atom_bouge_abs

	endif ! block with continue_sundae == 2
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	if (continue_sundae.eq.1) then
		
		
		write( srank, '(i2)' )  rank
		recup = sortie(1:lenfnam)//trim(adjustl(srank))//'.in'
		open(unit=112,file=recup,status='old')		
		
		!!write(*,*) 'recup', recup
		
		do i=0, totiter-1                      ! attention i=totiter doit etre pris en compte
			read (112,*) Path(i)%q
			read (112,*)
			read (112,*) Path(i)%p
			read (112,*)
			read (112,*) Path(i)%Lyap
			read (112,*)
			read (112,*) Path(i)%project
			read (112,*)
			read (112,*) Path(i)%eigenvalue
			read (112,*)
		enddo
		do i = -N_extra, nmax+N_extra 
			read (112,*) A_n(i)
		enddo
			read (112,*) 
			read (112,*) theta_n
		close(112)
              !  write(*,*) ' alpha_max ', alpha_max,theta_n
                if ((theta_n.le.alpha_max).and.(theta_n.ge.0)) then
                 write(*,*) ' alpha_max ', alpha_max,theta_n
 	  	  itheta_n = nint(theta_n/alpha_max*real(nmax))
                else
                  itheta_n = nmax/2
                  theta_n  = alpha_max/2

                endif
		Path(totiter) =  Path(totiter-1) 
		oldLyap=sum(Path(0:totiter-1)%Lyap) !m  /real(totiter)
		!!write(*,*) 'oldLyap = ', oldLyap
		
		

		q=Path(0)%q
		do i=1,3
			!!write(*,*) 'i       ',i
			!!write(*,*) 'bary    ',sum(q(i,1:im))/dble(im)
			!!write(*,*) 'bary ref',sum(qref(i,1:im))/dble(im)
		enddo

		rang=0
		
		do i=0,totiter

			q    = Path(i)%q
			prot = Path(i)%p

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

			Path(i)%q=q

			call control_angular_momenta(prot,q)

			Path(i)%p=prot

		enddo

		absdmax(:)=-9999.0
		atom_bouge_abs=0
		iter=0
		do i=1,N
			pav(1:3)=(Path(iter)%q(1:3,i)-qref(1:3,i))**2
			absdist(i)=sqrt(sum(pav(1:3)))
			if (absdist(i).gt.absdmax(iter)) then
				absdmax(iter) = absdist(i)
				atom_bouge_abs  = i
			endif
		enddo
		absdmax_current = absdmax(iter)*1.d8
		!!write(*,*) ' absdmax(1,0)   ',absdmax(iter)*1.d8
		!!write(*,*) ' atom_bouge_abs ',atom_bouge_abs

	endif ! block with continue_sundae == 1
	

end subroutine LyapLanczos_equilibrage


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine LyapLanczos_allocate

	use tab_imm_m

	implicit none

!m	open(unit=30, file=kappaF,  action='write', status='replace')
!m	open(unit=31, file=kappaFd,  action='write', status='replace')
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!m	open(unit=24, file=data_mbar, action='write', status='replace')
!m	open(unit=244, file=dada_mbar, action='write', status='replace')

!m	write(244,'(i,i,i)') nmax+1,nmax+1,Totalmcmoves


        dts2racinem = dt/2.0/sqrt(9.270914743200000e-023)

	!!write(*,*) ' nmax ', nmax
	 
	! Array where the weight of every clone is stored
	acc = 0.d0
	ener0 = 0.d0
	enerK = 0.d0
	hamilt = 0.d0

	! Array which contains the numbers from 1 to nmax

	ltot=int(TotalTime/real(dt))

	!!write(*,*)'totiter', totiter
	gamma_sundae=gamma_sundae/(dt)
	!!write(*,*)'gamma_sundae' , gamma_sundae, gamma_sundae*dt,ltot
	waste_recycling=.true.
	!waste_recycling=.false.
	!!write(*,*) 'waste_recycling ?',  waste_recycling

	allocate(absdmax(0:2*totiter+10))
	allocate(S(0:2*totiter+10))

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


	absdmax(:)=0
	mcmoves=0
	it_art=0
	e=0.000001
	eigenvalue=0

	if (maxvec.lt.4) maxvec=4
	!!write(*,*) 'The maxvec in lanczos is set to .........:', maxvec



	allocate(Path(0:totiter+10)) 
	allocate(Pshoot(0:totiter+10))
	allocate(Pshift(0:2*totiter+10))

	oldLyap = 0.d0
	triallyap = 0.d0
	rapport = 0.d0

	allocate(absdist(0:N))
	allocate(Psel(0:totiter))
	allocate(q(1:3,1:im))
	allocate(p(1:3,1:im))
	iter=0
	Q=0
	!pi=4*atan(1._dpkind)
	
	if ((gamma_sundae*dt).le.100000) then
		rga = exp(-gamma_sundae*dt/two)
	else 
		rga=0
	endif
	
	sig(:,:) = sqrt(temperature*(one-rga**2))
	!!write(*,*) 'rga',gamma_sundae, gamma_sundae*dt, rga, sig(1,1)
	!!write(*,*) 'temp', temperature*erg2eV, erg2eV

	xref(1:n)=xp(1,1:n)
	yref(1:n)=xp(2,1:n)
	zref(1:n)=xp(3,1:n)

	qref(1:3,1:im)=xp(1:3,1:im)

	rang = 1 
	
	! Allocation ABF
	
	delta = real(alpha_max/real(nmax))
	delta_bin_theta_s2 = 0.5d0*real(alpha_max/real(nmax))
	theta_n = real(rank+1)/real(numproc+1)*alpha_max
	itheta_n = nint((theta_n)/(alpha_max)*real(nmax))
	N_extra = 10

	allocate(a_sto(0:nmax))
	allocate(theta(-N_extra:nmax+N_extra))
	allocate(u_A(-N_extra:nmax+N_extra))
	allocate(A_n(-N_extra:nmax+N_extra))
	allocate(P_A(-N_extra:nmax+N_extra))
	allocate(sum_P_A(-N_extra:nmax+N_extra))
	allocate(A_prime(-N_extra:nmax+N_extra))
	allocate(A_prime_num(-N_extra:nmax+N_extra))
	allocate(L2_moy_num(-N_extra:nmax+N_extra))
	allocate(O_moy_num(0:totiter,1:3,1:4))
	allocate(O_moy_estim(0:totiter,1:3,1:4))
	allocate(O_estim(0:totiter,1:3,1:4))
	allocate(histo_theta(-N_extra:nmax+N_extra))
	allocate(histo_Lyap(-N_extra:nmax+N_extra))
	allocate(theta_traj(1:Totalmcmoves))
	allocate(Lyap_traj(1:Totalmcmoves))
	allocate(ha_hb_traj(1:Totalmcmoves))
	
	A_prime_num = 0
	L2_moy_num  = 0
	O_moy_num   = 0
	histo_theta = 0.d0
	histo_Lyap  = 0.d0
	Lyap_traj   = 0.d0
	ha_hb_traj  = 0.d0
	
#if(PARASUN)
	allocate(MPI_histo_theta(-N_extra:nmax+N_extra))
	allocate(MPI_histo_Lyap(-N_extra:nmax+N_extra))
	allocate(MPI_P_A(-N_extra:nmax+N_extra))
	allocate(MPI_sum_P_A(-N_extra:nmax+N_extra))
	allocate(MPI_A_prime(-N_extra:nmax+N_extra))
	allocate(MPI_A_prime_num(-N_extra:nmax+N_extra))
	allocate(MPI_L2(-N_extra:nmax+N_extra))
	allocate(MPI_L2_moy_num(-N_extra:nmax+N_extra))
	allocate(MPI_O_moy_num(0:totiter,1:3,1:4))
	allocate(MPI_O_moy_estim(0:totiter,1:3,1:4))
	
	MPI_histo_theta = 0.d0
	MPI_histo_Lyap = 0.d0
#endif

	!! initialisation paramètres de bias alpha pour reconstruction
	do j= -N_extra, nmax+N_extra
		theta(j)=(real(j*(alpha_max/real(nmax))))!*sqrt(9.270914743200000e-023)/300.0
		!!write(*,*)'bias', theta(j)  
	enddo

        do i = 0,nmax
	a_sto(i) = 1.d0-2.d0*((1.d1**(-2.d0-2.d0*dble(i)/dble(nmax))))
        enddo
!	write(6,*)'a_sto = ', a_sto


	
!	theta = 0.d0
!	do i=0,nmax
!		theta(i) = theta(i)
!	enddo



end subroutine LyapLanczos_allocate 



end module sundae_module



