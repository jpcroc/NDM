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

!INTERFACE 
!   FUNCTION asinhsqrt (x)
!     use T_kind_param_m, ONLY:  double
!     REAL(double) :: asinhsqrt
!     REAL(double), INTENT(IN) :: x
!   END FUNCTION asinhsqrt
!END INTERFACE

	integer, parameter :: dpkind=selected_real_kind(13)
	real(double) :: t0,t1,elaps_1,tbuffer1,tbuffer2,t_lanczos
	integer, dimension(:), allocatable, save  :: ipovois
	
	
	integer :: continue_sundae, reprise_A, maxvec, totiter
	integer :: Totalmcmoves, depart_boucle_nbclones
	real(double) :: h_A_max, h_ba_max, h_ba_min, h_ba, h_ba_I, h_temp
	real(double) :: dt, TotalTime, gamma_sundae, Temperature
	real(double) :: alpha_max, teq, delta_x,ss, tequilib
	real(double) :: drift,dts2racinem ,ra1
!        real(double) :: asinhsqrt

	character(len=128) :: recup
	character(len=128) :: srank
        character(len=128) :: entree
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
		real (double)  Lyapu(1:4)
		real (double)  eigenvals(1:4)
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
	integer :: iterfw, iterbw,nmax, tprim,nthetamax

	real (double) :: ener0 
	real (double) :: enerK
	real (double) :: triallyap
	real (double) :: rapport
	real (double) :: hamilt
	real (double) ,dimension(:), allocatable:: absdist
	real (double) ,dimension(:), allocatable:: Psel
	real (double) ,dimension(:), allocatable:: PaPsel
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
	real(double) :: Lyap_max = 9.d0  ! anciennement 20.0
        real(double) :: Lyap_min = -9.d0
	integer :: iLyap
	real(double) :: somme
	real(double) :: Lyap, gmin, P_n, pi_n, pi_n1, sum_A
	real(double) biais,oldbiais
	real(double) a_stol
      	real(double),dimension(:),  allocatable,save  :: accept_theta
      	real(double),dimension(:),  allocatable,save  :: MPI_accept_theta

      	real(double),dimension(:),  allocatable,save  :: refus_theta
      	real(double),dimension(:),  allocatable,save  :: MPI_refus_theta

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
	real(double), dimension(:,:), allocatable, save :: Lyap_traj
	real(double), dimension(:), allocatable, save :: ha_hb_traj
	real(double), dimension(:), allocatable, save :: eigmax_traj
	real(double), dimension(:), allocatable, save :: Lyapmax_traj

	! Complements MPI
	integer :: rank = 0, numproc = 1
	integer :: ificout
        character*120 :: ficout	
	real(double) :: estim_deno,cumul_deno=0.d0,MPI_cumul_deno

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
        real(double) :: miasinhsqrt

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
		Pcourant%Lyap       = miasinhsqrt(eigenvalue)
	        do i=1,4 
		  Pcourant%Lyapu(i)= miasinhsqrt(eigenvals(i))
                  Pcourant%eigenvals(i)=eigenvals(i)
                enddo
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
	        Pcourant%Lyap=miasinhsqrt(eigenvalue)
	        do i=1,4 
		 Pcourant%Lyapu(i)=miasinhsqrt(eigenvals(i))
		enddo
		Pshoot(iterfw)%Project=Pcourant%Project
		Pshoot(iterfw)%eigenvalue=Pcourant%eigenvalue
		Pshoot(iterfw)%Lyap=Pcourant%Lyap
               	Pshoot(iterfw)%Lyapu=Pcourant%Lyapu
               	Pshoot(iterfw)%eigenvals(1:4)=Pcourant%eigenvals(1:4)


		Pshoot(iterfw+1)   = Pcourant
		it_trajectory  = it_trajectory + 1
		new_projection = .false.
	enddo
	!!write(*,*) 'Lanczos iterations...:',lanczos_iter
	!!write(*,*) '           forces....:',lanczos_iter*maxvec*2

	triallyap=max(Lyap_min,SUM(Pshoot(0:totiter-1)%Lyap))       !m    /real(totiter)
	
	!!!!!!! fonction heaviside pour maintenir le point initial de la trajectoire dans 
	!!!!!!! le bassin de depart: calcul dist max pour le point x(0) de la trajectoire
	absdmax(0)=-9999.0
	atom_bouge_abs=0
	do i=1,N
		absdist(i)=sqrt((Pshoot(0)%q(1,i)-qref(1,i))**2 + &
		(Pshoot(0)%q(2,i)-qref(2,i))**2 + &
		(Pshoot(0)%q(3,i)-qref(3,i))**2)
		if (absdist(i).gt.absdmax(0)) then
			absdmax(iterbw) = absdist(i)
			atom_bouge_abs    = i
		endif
	enddo
!        write(*,*) ' iterbw ', iterbw 
!        stop
	absdmax(0)      = absdmax(0)*1.0d8
	!!write(*,*) 'iterbw absdmax ',iterbw,absdmax(iterbw)
	!!!!!!!!!!!! ATTENTION AUX RELATIONS DU BILAN DETAILLE

	if (absdmax(0).lt.h_A_max) then
		ra1=exp(theta(j)*(-triallyap+oldLyap))
                call marginal_calcul(triallyap)
		rapport=exp(biais-oldbiais)
                write(*,*) rank,'rapports ',ra1,rapport
	else
		rapport=0.d0
	endif

	write (ificout,*) ' itrialLyap oldLyap = ', triallyap , oldLyap 

	ranf=genrand()
	mcconf=min(1.d0,rapport)

	if (ranf.lt.mcconf) then
                accept_theta(j/10) = accept_theta(j/10) +1
		acc=acc+1
		Path(0:totiter) = Pshoot(0:totiter)
		oldLyap=triallyap
		oldbiais= biais
		write(ificout,*) 'traj. acceptée',mcmoves,'rapport =',rapport
		absdmax_current = absdmax(0)  
		!!write(ificout,*) 'absdmax_current', absdmax_current
	else
                refus_theta(j/10) = refus_theta(j/10) + 1
		write(ificout,*) 'traj. refusée ',mcmoves,'rapport =',rapport
	endif
	write (ificout,*) 'Shooting ....taux d''acceptation',real(acc)/real(mcmoves) 
        call marginal_calcul(oldLyap)
       	write (ificout,*) 'theta moyen ', sum(P_A(0:totiter)*theta(0:totiter))

!! fin du shooting

end subroutine LyapLanczos_shooting

subroutine marginal_calcul(L)
	implicit none
	real(double) :: L

	u_A(:)        = L*theta(:) - A_n(:)     
	gmin          = minval(u_A)
	u_A(:)        = u_A(:) - gmin
	P_A(:)        = exp(-u_A(:))
	somme         = sum(P_A(0:nmax))
	P_A           = P_A/somme
        biais         = log (somme)-gmin
	
end subroutine marginal_calcul

subroutine LyapLanczos_theta
	implicit none
	real(double) :: genrand,theta_o, cumul 
	integer :: i

        theta_o = theta_n

	ranf = genrand()
!        write(*,*) 'ranf ',ranf
	itheta_n = 0
        cumul = 0.d0
        do i=0,nmax
        cumul = cumul + P_A(i) 
         if (ranf.gt.cumul) then 
	  itheta_n = i 
 !         write(*,*) 'cumul ',cumul 
         endif
        enddo
!        write(*,*) 'cumul ',cumul
!	theta_n = dble(itheta_n)*alpha_max/dble(nmax)
       
!        write(*,*) ' P_A =',P_A
 	theta_n = theta(itheta_n) 
	write(ificout,*) 'theta n / o ' ,theta(itheta_n) , theta_o
	theta_traj(mcmoves) =  theta(itheta_n)

end subroutine LyapLanczos_theta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine LyapLanczos_shifting

	implicit none
	
	real(double) :: genrand,miasinhsqrt,Ly
	
	
	j = itheta_n

	pix = int(genrand()*(totiter+1))  ! le cas pix=totiter+1 est impossible par construction
			    			          ! car genrand < 1
			          				  ! pix est compris entre 0 et L soit L+1 valeur possibles 
	write(ificout,'(" Shifting............................:",3i6)') mcmoves,totiter,pix

	write(ificout,'(" Shift.rec...pix+1->pix+totinter....:",i6," to ",i6)') pix+1,pix+totiter
	! we recycle the pix,pix+1, pix+2, ...., pix+totiter because there are totiter+1 points 
	do k=0,totiter
		Pshift(pix+k) = Path(k)
	enddo

	write(ificout,'(" shift.bw ......pix->0..............:",i6," to 0")') pix
	! we calculate the pix-1,....,0
	new_projection=.false.
	it_trajectory=0

	Pcourant   = Pshift(pix)
	eigenvalue = Pcourant%eigenvalue

	do iterbw=pix,1,-1
		call mapping_P_Verlet(Pcourant%q,Pcourant%p,-dt,N,q1s2)  ! -dt car backward en position verlet
		call lanczos(N,maxvec,q1s2,new_projection,Pcourant%project)
		Pcourant%eigenvalue     = eigenvalue
		Pcourant%eigenvals(1:4) = eigenvals(1:4)
	        Pcourant%Lyap=miasinhsqrt(eigenvalue)
	        do i=1,4 
		  Pcourant%Lyapu(i)=miasinhsqrt(eigenvals(i))
		enddo
		Pshift(iterbw-1)=Pcourant
		it_trajectory=it_trajectory+1 
		new_projection=.false.
	enddo

	write(ificout,'(" shift.fw..pix+totiter->2*totinter..:" ,i6, " to ", i6)') pix+totiter,2*totiter
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
	        Pcourant%Lyap=miasinhsqrt(eigenvalue)
	        do i=1,4 
		 Pcourant%Lyapu(i)=miasinhsqrt(eigenvals(i))
		enddo

!		if (eigenvalue.le.0.d0) then 
!m			Pcourant%Lyap=asin(dt*sqrt(-eigenvalue)/2.d0)*2.d0  ! Lyapunov en 1/2 stocké avant
!                 omegadts2 = dts2racinem*sqrt(-eigenvalue)
!	         Pcourant%Lyap=asinh(omegadts2)*2.d0 ! 
!		else 
!		 Pcourant%Lyap=0.d0 ! Pshift(iterfw-1)%Lyap ! Lyapunov en 1/2 stocké avant 
!		endif

		Pshift(iterfw)%Lyap            = Pcourant%Lyap
		Pshift(iterfw)%Lyapu           = Pcourant%Lyapu
		Pshift(iterfw)%Project         = Pcourant%Project
		Pshift(iterfw)%eigenvalue      = eigenvalue
		Pshift(iterfw)%eigenvals(1:4)  = eigenvals(1:4)

		Pshift(iterfw+1)=Pcourant

		it_trajectory=it_trajectory+1 
		new_projection=.false.
	enddo

	hamilt=ener0

	S(:)=0
	do l=0,totiter  ! sommation sur les 1+L "path proposals" possibles
            Ly = max(Lyap_min,sum(Pshift(l:l+totiter-1)%Lyap))
            call marginal_calcul(Ly)
            S(l) = Biais
	enddo

	!!!!!!!!!!!! calcul de la fonction echelon pour le point de départ 
	!!!!!!!!!!!! de la trajectoire N sui 2N passi della traj shiftata

	absdmax(:)=-9999.0
	atom_bouge_abs=0

	do l=0,2*totiter ! il y a 1+L (1+totiter) proposals + totiter positions décalée en totiter
		do i=1,N
			absdist(i)=sqrt((Pshift(l)%q(1,i)-qref(1,i))**2 + & 
				 (Pshift(l)%q(2,i)-qref(2,i))**2 + & 
				 (Pshift(l)%q(3,i)-qref(3,i))**2) 
			if (absdist(i).gt.absdmax(l)) then 
				absdmax(l)=absdist(i) 
			endif 
		enddo
	enddo
	absdmax(:)=absdmax(:)*1.0e8
	!!write(*,*) 'absdmax : '
	!!write(*,*) absdmax(0) 
	
	z=0.d0
	do l=0,totiter ! boucles sur les chemins proposés possibles
		! Here is P_sel from the paper
		if (absdmax(l).lt.h_A_max) then  ! on a correspondance maintenant entre absdmax et Pshift
!			Psel(l) = exp(theta(j)*S(l))
                    	Psel(l) = exp(S(l))
                        PaPsel(l) =exp(A_n(0))
		else
			PaPsel(l) = 0.d0
			Psel(l)   = 0.d0
		endif
		z = z +  Psel(l)
	enddo
        PaPsel(0:totiter) = PaPsel(0:totiter)/z
	Psel(0:totiter)   = Psel(0:totiter)/z

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
		oldLyap       = max(Lyap_min,SUM(Path(0:totiter-1)%Lyap))
	        call marginal_calcul(oldLyap)
                oldbiais = biais
	else
		write(*,*) " Problème avec le shifting "
                stop
	endif

        write(ificout,*) 'shift  lyap = ', newtraj-pix, oldLyap,oldbiais
        write(ificout,*) 
!	Lyap_traj(mcmoves) = oldLyap  !  *300.0/sqrt(9.270914743200000e-023)
	
        forall (i=1:4) Lyap_traj(mcmoves,i) =  SUM(Path(0:totiter-1)%Lyapu(i))

        eigmax_traj(mcmoves)  = minval(Path(0:totiter-1)%eigenvalue)
        Lyapmax_traj(mcmoves) = maxval(Path(0:totiter-1)%Lyap)

!	write(*,*) 'min val',rank, eigmax_traj(mcmoves),Lyapmax_traj(mcmoves) 

!        do i = 0,totiter-1
! 	 write(*,*) 'eigen ', i, Path(i)%eigenvalue,  Path(i)%Lyap
!        enddo
	

end subroutine LyapLanczos_shifting


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine LyapLanczos_ABF

	use lanczos_defs
	use tab_imm_m

	implicit none

	! Calcul du nouveau biais
	Lyap          = oldLyap
	u_A(:)        = Lyap*theta(:) - A_n(:)     
	gmin          = minval(u_A)
	u_A(:)        = u_A(:) - gmin
	P_A(:)        = exp(-u_A(:))
	somme         = sum(P_A(0:nmax))
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

		A_prime_num(theta_tilde)   =   P_n * Lyap      + A_prime_num(theta_tilde)
		L2_moy_num(theta_tilde)    =   P_n * Lyap*Lyap + L2_moy_num(theta_tilde)
	enddo
	
#if(PARASUN)
	! Mise en commun parallele ---- 2/2 !

	call MPI_REDUCE(A_prime_num,MPI_A_prime_num,nmax+2*N_extra+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, MPI_COMM_WORLD,ierror)
	call MPI_REDUCE(L2_moy_num,MPI_L2_moy_num,nmax+2*N_extra+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, MPI_COMM_WORLD,ierror)
	call MPI_REDUCE(sum_P_A,MPI_sum_P_A,nmax+2*N_extra+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, MPI_COMM_WORLD,ierror)

	MPI_A_prime = MPI_A_prime_num / (MPI_sum_P_A + 1.d-6)
	MPI_L2  = MPI_L2_moy_num / (MPI_sum_P_A + 1.d-6)
	
!	do k = 0,totiter 
!		O_moy_estim(k,1:3,1:4)  =   MPI_O_moy_num(k,1:3,1:4) / (MPI_sum_P_A(0) + 1.d-6)
!	enddo
	
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
	write(*,*) 'pas fini '
	stop
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

subroutine LyapLanczos_input

	use lanczos_defs
	use tab_imm_m
	implicit none
	real(double) theta_f, A_prime_f, L2_f
	character(len=128) :: fic_abf, fic4, fic5 , fic6, fic7!,srank
        integer :: ific4,ific5,ific7	
        write(*,*) 'coucou 0 '
	if (rank.eq.0) then
                fic_abf =  trim(entree)//'.data_abf'
		open(unit=111,file=fic_abf,action='read',status='old')
		fic6 = trim(entree)//'.pmf'
		open(unit=774, file=fic6, action='read', status='old')
        endif
		write( srank, '(i2)' )  rank
                write(*,*) ' proc ', srank
		fic4 = trim(entree)//'.'//trim(adjustl(srank))//'.pathinit'
                fic5 = trim(entree)//'.'//trim(adjustl(srank))//'.eigen'
	        fic7 = trim(entree)//'.'//trim(adjustl(srank))//'.project'
!                read(*,*) 'fic4 fic5 fic6 ',fic4 , fic5 , fic6 , fic7
                ific4 = 800+rank
                ific5 = 820+rank
                ific7 = 840+rank

		open(unit=ific4, file=fic4, action='read', status='old')
		open(unit=ific5, file=fic5, action='read', status='old')
		open(unit=ific7, file=fic7, action='read', status='old')

#if(PARASUN)
	if (rank.eq.0) then
               write(*,*) ' fic_abf = ' , fic_abf
	       do j=-N_extra,nmax+N_extra
			read(111,'(3(E15.6E3))') theta_f, A_prime_f, L2_f
			write(*,*) 'theta ', theta(j), theta_f
			MPI_A_prime(j) = A_prime_f
	       enddo
         
            !   do i = - N_extra , nmax + N_extra
            !    read (774,*) theta(i), A_n(i)
            !   enddo
	endif	

#endif
	if (rank.eq.0) then
	 close(111)
         close(774)
	endif

	read (ific4,*) Path(0)%q
	read (ific4,*)
	read (ific4,*) Path(0)%p
	read (ific4,*)
	read (ific4,*) theta_n 

!      verification des rotations angulaires / qref
!       construction du chemin 
	q    = Path(0)%q
	qrot=q-qref
        rang=-1
	call control_angular_momenta(qrot,qref)
	call control_angular_momenta(qrot,qref)
	q = qrot + qref
	rang=-1
	call control_angular_momenta(qrot,qref)
	rang=-1
	Path(0)%q=q

	prot = Path(0)%p

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

	call control_angular_momenta(prot,q)
	Path(0)%p=prot

        do i = 1,totiter-1
	 Path(i) = Path(0)
        enddo

	do i=0, totiter-1                      ! attention i=totiter doit etre pris en compte
	 	read (ific5,*) Path(i)%Lyap,Path(i)%eigenvalue
        enddo
	do i=0, totiter-1   
		read(ific7,*) Path(i)%project
		read(ific7,*)
	enddo

	close(ific4)
	close(ific5)
	close(ific7)

end subroutine LyapLanczos_input

subroutine LyapLanczos_output

	use lanczos_defs
	use tab_imm_m

	implicit none
	real(double) theta_f, A_prime_f, L2_f, somme
	character(len=128) :: fictheta, fic1, fic2, fic3, fic4, fic5 , fic6, fic7, smcmoves
        integer :: ific4,ific5,ific7	
	! ----------- Sortie des observables ----------- !
#if(PARASUN)
	! Mise en commun parallele !	
	call MPI_REDUCE(histo_Lyap,MPI_histo_Lyap,nmax+2*N_extra+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, MPI_COMM_WORLD,ierror)

	MPI_histo_Lyap = MPI_histo_Lyap/dble(numproc)

	call MPI_BCAST(MPI_histo_Lyap,nmax+2*N_extra+1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
	! Mise en commun parallele !


	! Mise en commun parallele !	
	call MPI_REDUCE(accept_theta,MPI_accept_theta,nthetamax+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, MPI_COMM_WORLD,ierror)

	MPI_accept_theta = MPI_accept_theta/dble(numproc)

	call MPI_BCAST(MPI_accept_theta,nthetamax+1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
	! Mise en commun parallele !

	call MPI_REDUCE(refus_theta,MPI_refus_theta,nthetamax+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, MPI_COMM_WORLD,ierror)

	MPI_refus_theta = MPI_refus_theta/dble(numproc)

	call MPI_BCAST(MPI_refus_theta,nthetamax+1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)

#endif

	if (rank.eq.0) then
		write( smcmoves, '(i6)' )  mcmoves
		fictheta = trim(sortie)//'.accept_theta'
		fic1 = trim(sortie)//'.histo_abf'
		fic2 = trim(sortie)//'.obs_1_'//trim(adjustl(smcmoves))
		fic3 = trim(sortie)//'.obs_2_'//trim(adjustl(smcmoves))
		open(unit=111, file=data_abf, action='write', status='replace')
		open(unit=1110, file=fic1, action='write', status='replace')
		open(unit=1111, file=fic2, action='write', status='replace')
		open(unit=1112, file=fic3, action='write', status='replace')
		open(unit=1113, file=fictheta, action='write', status='replace')
		fic6 = trim(sortie)//'.pmf'
		open(unit=774, file=fic6, action='write', status='replace')
        endif
		write( srank, '(i2)' )  rank
		fic4 = trim(sortie)//'.'//trim(adjustl(srank))//'.pathinit'
                fic5 = trim(sortie)//'.'//trim(adjustl(srank))//'.eigen'
	        fic7 = trim(sortie)//'.'//trim(adjustl(srank))//'.project'
!                write(*,*) 'fic4 fic5 fic6 ',fic4 , fic5 , fic6 , fic7i
                ific4 = 772+10*rank
                ific5 = 773+10*rank
                ific7 = 775+10*rank

		open(unit=ific4, file=fic4, action='write', status='replace')
		open(unit=ific5, file=fic5, action='write', status='replace')
		open(unit=ific7, file=fic7, action='write', status='replace')

#if(PARASUN)

	call MPI_REDUCE(O_moy_num,MPI_O_moy_num,4*3*(totiter+1),MPI_DOUBLE_PRECISION,MPI_SUM,0, MPI_COMM_WORLD,ierror)
	call MPI_REDUCE(cumul_deno,MPI_cumul_deno,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierror)

	if (rank.eq.0) then
		somme = sum(MPI_histo_theta)
		MPI_histo_theta = MPI_histo_theta/somme*real(nmax)
		somme = sum(MPI_histo_Lyap)
		MPI_histo_Lyap = MPI_histo_Lyap/somme
		do j=-N_extra,nmax+N_extra
			! Renormalisation des variables pour les sorties fichier
			theta_f = theta(j)           !  *sqrt(9.270914743200000e-023)/300.0
			A_prime_f = MPI_A_prime(j)   !  *300/sqrt(9.270914743200000e-023)
			L2_f = MPI_L2(j)             !   *300/sqrt(9.270914743200000e-023)*300/sqrt(9.270914743200000e-023)
			! Sortie fichier : theta, histo_theta, Lyap_moyen (A_prime)
			write(111,'(4(E15.6E3))') theta_f, A_prime_f, L2_f,MPI_histo_theta(j)
			write(1110,'(2(E15.6E3))') - dble(j)/dble(nmax)*Lyap_max, MPI_histo_Lyap(j)
		enddo
                do j=0,nthetamax
                        somme = MPI_accept_theta(j)+MPI_refus_theta(j)
			if (somme.eq.0.d0) somme = 1.d-5
		 	write(1113,*) theta(j*10),MPI_accept_theta(j)/somme,somme
                enddo

		O_moy_estim(0:totiter,2,1:4) = MPI_O_moy_num(0:totiter,2,1:4)/MPI_cumul_deno
		O_moy_estim(0:totiter,1,1:4) = MPI_O_moy_num(0:totiter,1,1:4)/MPI_sum_P_A(0)


		do j=0,totiter
			write(1111,'(5(E15.6E3))') O_moy_estim(j,1,1), O_moy_estim(j,1,2) , O_moy_estim(j,1,3), sum(O_moy_estim(j,1,1:4)), MPI_sum_P_A(0)
			write(1112,'(5(E15.6E3))') O_moy_estim(j,2,1), O_moy_estim(j,2,2) , O_moy_estim(j,2,3), sum(O_moy_estim(j,2,1:4)), MPI_cumul_deno
		enddo
                do i = -N_extra, nmax+N_extra
                 write (774,*) theta(i), A_n(i)
                enddo

	endif	
#else
        write(*,*) ' Programme séquentiel non terminé '
	stop
#endif
	if (rank.eq.0) then
		close(111)
		close(1110)
		close(1111)
		close(1112)
		close(1113)
                close(774)
	endif
	! ----------- Sortie des observables ----------- !
		

!!!!!!!!!!!!!!!!! nouveau format  
  
	write (ific4,*) Path(0)%q
	write (ific4,*)
	write (ific4,*) Path(0)%p
	write (ific4,*)
	write (ific4,*) theta_n 

	do i=0, totiter-1                      ! attention i=totiter doit etre pris en compte
		write (ific5,*) Path(i)%Lyap,Path(i)%eigenvalue
        enddo
	do i=0, totiter-1   
		write (ific7,*) Path(i)%project
		write (ific7,*)
	enddo

	close(ific4)
	close(ific5)
	close(ific7)

	recup = sortie(1:lenfnam)//'.'//trim(adjustl(srank))//'.Lyapunov'
        
	open(unit=300+rank,file=recup,status='replace')		
		
	do k = 1, mcmoves
		write(300+rank,'(5(E15.6E3))') ha_hb_traj(k),theta_traj(k),Lyap_traj(k,1),Lyap_traj(k,2),eigmax_traj(k)
	enddo

	close(300+rank)
	! ----------- Sortie fichier Lyapunov-trajectoire ----------- !



!!!!!!!!!!!!!!!!! fin nouveau format


	!!write(*,*) 'recup', recup
! ----------- Sortie fichier de recuperation ----------- !
!	write( srank, '(i2)' )  rank
!	recup = sortie(1:lenfnam)//trim(adjustl(srank))//'.in'
!	if (continue_sundae.eq.1) then
!		open(unit=112,file=recup,status='replace')		
!	else if (continue_sundae.eq.2) then
!		open(unit=112,file=recup,status='replace')
!	endif
!      	
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

end subroutine LyapLanczos_output


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine LyapLanczos_vac! (xp)

	use lanczos_defs
	use tab_imm_m

	implicit none

	real(double) :: genrand
#if(PARASUN)
	
	call MPI_INIT(ierror)     !! On initialise MPI

	call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)     !! Attribue a rank le numero de processeur

	call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)  !! Attribue a numproc le nombre de processeur
	
	call init_random_seed(100000*(rank+1))

#else
	call init_random_seed(100000)
#endif


	!!!!!!!!!!!!!!!!!!!!!!      Phase d'initialisation       !!!!!!!!!!!!!!!!!!!!!!!!

	call read_sundae                         !!! Modif 21.05.14

        call LyapLanczos_allocate                !!! Modif 02.06.14

	call LyapLanczos_equilibrage             !!! Modif 22.01.15
	
	call LyapLanczos_init_tests              !!! Modif 22.01.15
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	ficout  = trim(sortie)//'.'//trim(adjustl(srank))//'.out'
        ificout = 860+rank
	open(unit=ificout, file=ficout, action='write', status='replace')

	do mcmoves = 1 , Totalmcmoves                       !mcmoves = 1,M in the paper

		!write(*,'(" Proc number ", i6, "says hello")'), rank
	
	
		it_art=mcmoves

		ix=int((totiter)*genrand())   ! le cas ix = totiter ne doit pas etre possible car on stocke jusqu'à totiter. 
		!!write(*,*) 
		!!write(*,*) 
		write(ificout,'(" Itération..:",i6," sur un total de ",i6)') ,mcmoves,Totalmcmoves
		write(ificout,'(" Shooting......................:",3i6)') mcmoves,totiter,ix

		!!!! shooting deterministique
		! Check: if depart_boucle_nbclones,NbClones of if there are different WHY?


		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!
		
		call LyapLanczos_shooting              !!! Modif 02.06.14

	        call LyapLanczos_theta	               !!! Modif 09.01.15
	
		call LyapLanczos_shifting              !!! Modif 29.01.15

		call LyapLanczos_Obs		  !!! Modif 13.06.14
	
		call LyapLanczos_ABF		  !!! Modif 10.06.14
		
!		write (*,*) isauvegarde
		if (0.eq.mod(mcmoves,isauvegarde)) then 
			call LyapLanczos_output
		endif
		
		!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	enddo !! mcmoves avec incrément + 1 clones


!	call LyapLanczos_output

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
	ha_hb_traj(mcmoves) =  h_Fg(totiter)+h_Fd(totiter)+h_F(totiter)+h_FI(totiter)
	
	! Estimation en moyenne glissante sur [0,2L] trajectoire shiftee
	h_Fg(:) = 0.d0
	h_Fd(:) = 0.d0
	h_F(:)  = 0.d0
	h_FI(:) = 0.d0
	O_estim(0:totiter,2,1:4) = 0.d0


	do k = 0,totiter 
		do j = 0,totiter 
			if ((absdmax(j+k).ge.h_ba_min) .and. (absdmax(j+k).lt.h_ba) .and. (absdmax(j).lt.h_A_max)) then 
				h_Fg(k) = h_Fg(k) + PaPsel(j) ! + 1.0
			endif
			if ((absdmax(j+k).ge.h_ba) .and. (absdmax(j+k).lt.h_ba_max) .and. (absdmax(j).lt.h_A_max)) then 
				h_Fd(k) = h_Fd(k) + PaPsel(j) ! + 1.0
			endif
			if ((absdmax(j+k).ge.h_ba_max) .and. (absdmax(j+k).lt.h_ba_I) .and. (absdmax(j).lt.h_A_max)) then 
				h_F(k) = h_F(k)   + PaPsel(j) ! + 1.0
			endif
			if ((absdmax(j+k).ge.h_ba_I) .and. (absdmax(j).lt.h_A_max)) then 
				h_FI(k) = h_FI(k) + PaPsel(j) ! + 1.0
			endif
		enddo !j
		O_estim(k,2,1) = h_Fg(k) 
		O_estim(k,2,2) = h_Fd(k)
		O_estim(k,2,3) = h_F(k)
		O_estim(k,2,4) = h_FI(k)
	enddo !k

	cumul_deno                   = cumul_deno + sum(PaPsel(0:totiter))
	O_moy_num(0:totiter,1:3,1:4) = O_estim(0:totiter,1:3,1:4) + O_moy_num(0:totiter,1:3,1:4)
	
	! -------- Histogramme du Lyap --------- !
	iLyap = nint(-Lyap/Lyap_max*real(nmax))
	if ((iLyap.ge.0).and.(iLyap.le.nmax)) then
		histo_Lyap(iLyap) = histo_Lyap(iLyap) + 1.0d0
	endif

end subroutine LyapLanczos_Obs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine LyapLanczos_init_tests

	!use lanczos_defs
	use tab_imm_m
	implicit none
        real(double) :: miasinhsqrt
	
	write (*,'("1st traj before shooting ........:",i6)') nint(TotalTime/dt)
	rang=1
	call cpu_time(t0)
	lanczos_iter=0
	Pcourant = Path(0)
	temp=0.d0
	eigenvalue = Pcourant%eigenvalue
        new_projection = .true.

	do iter=0,totiter-1
		call mapping_P_Verlet(Pcourant%q,Pcourant%p,dt,N,q1s2)
		call cpu_time(tbuffer1)
		lanczos_iter=lanczos_iter+nl_iter
		call lanczos(N,maxvec,q1s2,new_projection,Pcourant%project)!!!
		call cpu_time(tbuffer2)
		t_lanczos=t_lanczos + (tbuffer2-tbuffer1)
		Pcourant%Lyap = 0.d0 
 	        Pcourant%Lyap=miasinhsqrt(eigenvalue)
	        do i=1,4 
		  Pcourant%Lyapu(i)=miasinhsqrt(eigenvals(i))
		enddo
		Path(iter+1)    = Pcourant
		Path(iter)%Lyap = Pcourant%Lyap
		new_projection = .false.
                write(*,*) ' eigenvalue ',eigenvalue
                write(*,*) ' eigenvals  ',eigenvals


	enddo  !!!!!
!			oldLyap=sum(Path(0:totiter-1)%Lyap) !m /real(totiter)
	write(*,*) " oldLyap= ", oldLyap,temp !m /real(totiter)
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
		write(*,*) 'temp cin traject, std et cible ',temp,tempvar,temperature/KtoERG
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
		!!do i=1,3
			!!write(*,*) ' temperature translation ',  xbar(i)**2/m_i(1,1)/dble(im)/KtoERG*dble(N)
		!!enddo
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

	write(*,*) 'tempiter après ',iterbw,tempiter
	write(*,*) 'diff ',tempiter-enerK

end subroutine LyapLanczos_init_tests


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine LyapLanczos_equilibrage

	use tab_imm_m
	implicit none

	real(double) miasinhsqrt
		
        SELECT CASE (continue_sundae)
        CASE (1)
         
        		write( srank, '(i2)' )  rank
		recup = trim(entree)//'.'//trim(adjustl(srank))//'.in'
		open(unit=112,file=recup,status='old')		
		
		!!write(*,*) 'recup', recup
		
		do i=0, 250-1  ! totiter -1                  ! attention i=totiter doit etre pris en compte
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
	        do i=250-1,totiter -1 
                      Path(i) = Path(250-1)
	        enddo
		do i = -N_extra, nmax+N_extra 
			read (112,*) A_n(i)
		enddo
			read (112,*) 
			read (112,*) theta_n
		close(112)
              !  write(*,*) ' alpha_max ', alpha_max,theta_n
!                if ((theta_n.le.alpha_max).and.(theta_n.ge.0)) then
!                 write(*,*) ' alpha_max ', alpha_max,theta_n
! 	  	  itheta_n = nint(theta_n/alpha_max*real(nmax))
!                else
                  itheta_n = nmax
                  theta_n  = alpha_max
!               endif
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

        CASE (2)

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
		oldLyap=max(Lyap_min,sum(Path(0:totiter-1)%Lyap))      !m /real(totiter)
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

        CASE (3)

        call LyapLanczos_input

        CASE DEFAULT  !   cas pas encore traité
                 write(*,*) 'pb car continue_sundae = ',continue_sundae
          	
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
			new_projection      = .false.  
		        Pcourant%Lyap       = miasinhsqrt(eigenvalue)
		        do i=1,4 
			  Pcourant%Lyapu(i)=miasinhsqrt(eigenvals(i))
			enddo
!			if (eigenvalue.lt.0.0) then 
!		!m	Path(iter)%Lyap =  asin(dt*sqrt(-eigenvalue)/2.d0)*2.d0
!                         omegadts2 = dts2racinem*sqrt(-eigenvalue)
!	                 Pcourant%Lyap=asinh(omegadts2)*2.d0 !
!			else 
!	  		 Path(iter)%Lyap = 0.d0  ! 
!			endif

		enddo
	!!write(*,*) 'initialisation faite: go with Lanczos'
	!!write(*,*) 'itab', itab
	!!write(*,*) 'itetabvois', itetabvois
	!!write(*,*) 'ltabvois', ltabvois

        END SELECT
        WRITE(*,*)  'case Done'

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
	allocate(PaPsel(0:totiter))
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

	qref(1:3,1:im)=xp(1:3,1:im)

	rang = 1 
	
	! Allocation ABF
	
	delta = real(alpha_max/real(nmax))
	delta_bin_theta_s2 = 0.5d0*real(alpha_max/real(nmax))
	theta_n = real(rank+1)/real(numproc+1)*alpha_max
	itheta_n = nint((theta_n)/(alpha_max)*real(nmax))
	N_extra = 10
	nthetamax = nmax/10

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
	allocate(Lyap_traj(1:Totalmcmoves,1:4))
	allocate(ha_hb_traj(1:Totalmcmoves))
	allocate(eigmax_traj(1:Totalmcmoves))
	allocate(Lyapmax_traj(1:Totalmcmoves))


	allocate(accept_theta(0:nthetamax))
	allocate(refus_theta(0:nthetamax))
	
	A_prime_num = 0
	L2_moy_num  = 0
	O_moy_num   = 0
	histo_theta = 0.d0
	histo_Lyap  = 0.d0
	Lyap_traj   = 0.d0
	ha_hb_traj  = 0.d0
        eigmax_traj = 0.d0
        Lyapmax_traj = 0.d0

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
!	allocate(MPI_O_moy_estim(0:totiter,1:3,1:4))
	allocate(MPI_accept_theta(0:nthetamax))
	allocate(MPI_refus_theta(0:nthetamax))	


	MPI_histo_theta = 0.d0
	MPI_histo_Lyap = 0.d0
	MPI_accept_theta = 0.d0
	MPI_refus_theta = 0.d0

#endif

	!! initialisation paramètres de bias alpha pour reconstruction
	do j= -N_extra, nmax+N_extra
		theta(j)=(real(j*(alpha_max/real(nmax))))!*sqrt(9.270914743200000e-023)/300.0
		!!write(*,*)'bias', theta(j)  
	enddo

        do i = 0,nmax
!	a_sto(i) = 1.d0-2.d0*((1.d1**(-2.d0-2.d0*dble(i)/dble(nmax))))
!	a_sto(i) = 1.d0-1.d0*((1.d1**(-2.d0-2.5d0*dble(i)/dble(nmax))))
!	a_sto(i) = 1.d0-4.d0*((1.d1**(-2.d0-4.0d0*dble(i)/dble(nmax))))
!	a_sto(i) = 1.d0-1.d1**(-1.d0-5.0d0*dble(i)/dble(nmax))  ! bon parametrage pour alpha_max  = 2.d0
!	a_sto(i) = 1.d0-1.d1**(-1.d0-2.5d0*alpha_max*dble(i)/dble(nmax))  ! parametrage test pour alpha_max = 1.5d0
	a_sto = 0.99d0
	a_sto = a_stol


        enddo
!	write(6,*)'a_sto = ', a_sto


	
!	theta = 0.d0
!	do i=0,nmax
!		theta(i) = theta(i)
!	enddo



end subroutine LyapLanczos_allocate 



end module sundae_module



