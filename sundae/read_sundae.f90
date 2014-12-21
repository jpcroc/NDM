subroutine read_sundae()
!
! This subroutine reads initials paramaters
!

	!-----------------------------------------------
	!   M o d u l e s
	!-----------------------------------------------
	use gen_com_m,      ONLY : fnam, lenfnam
	use sundae_module,  ONLY : h_A_max, h_ba_min, h_ba_max, h_ba, h_ba_I, fnamtin, totiter, ss,          &
							 dt, Totalmcmoves, TotalTime, gamma_sundae, Temperature, alpha_max, theta,   &
							 teq, delta_x, a_sto, isauvegarde, KtoERG, sortie, continue_sundae,       &
							 maxvec, nmax, depart_boucle_nbclones, posfinal, tequilib,  &
							 data_abf, reprise_A

	implicit none


	namelist /input_sundae/ Totalmcmoves,TotalTime,dt,isauvegarde,Temperature,gamma_sundae,sortie,alpha_max, teq, delta_x, continue_sundae, reprise_A, nmax,maxvec,h_A_max,h_ba_min,h_ba_max,h_ba,h_ba_I



	h_A_max  = 4.5d-1
	h_ba_max = 1.6d0
	h_ba_min = 1.0d0
	h_ba     = 1.3d0
	h_ba_I   = 1.8d0

        isauvegarde = 100

	fnamtin = fnam(1:lenfnam)//'.tin'
	! variables de dynamique
	write(*,*) 'file name', fnamtin
	open(unit=777, file=fnamtin, status='unknown')
	read (777, nml=input_sundae)



	write(6,*)'usage %s:\n'

	!    Initialisation of the variables
	write(6,*)'TotalMcmoves = ', Totalmcmoves ! Total duration of the simulation 
	write(6,*)'No of steps = ',  TotalTime      ! Total duration of the simulation 

	write(6,*)'dt = ', dt ! Time step
	totiter = int(TotalTime)
	TotalTime=real(TotalTime)*dt

	write(6,*)'temperature = '  , temperature  ! temperature (kT)
	!   
	temperature=temperature*KtoERG
	!
	write(6,*)'friction*dt = '  , gamma_sundae  ! friction*dt
	! interval between 2 data records
	write(6,*)'sortie data_mbar = ',sortie  ! Name of the file where data are stored
!m	alpha_max = alpha_max*1.d12
	write(6,*)'alpha max = ', alpha_max 
	write(6,*)'t_equilib', teq
	write(6,*)'delta_X = ', delta_x

 	write(6,*)'continue_sundae = ', continue_sundae
	write(6,*)' nmax = ', nmax
!	write(6,*)'cluster? 0 no, nbclones yes'
!	write(*,*)'depart_boucle_nbclones',depart_boucle_nbclones
	write(*,*)'maxvec for Lanczos = ', maxvec


	ss=sqrt(temperature*1.66*1e-24*55.845)
	write(*,*) 'ss  = ', ss
	
	
	lenfnam   = index(sortie,' ')-1

!m	moyennes_mbar  =sortie(1:lenfnam)//'.data_moy'
!m	moyennes_mbar_denom  =sortie(1:lenfnam)//'.data_moy2'
	data_abf       = sortie(1:lenfnam)//'.data_abf'
!m	data_mbar      = sortie(1:lenfnam)//'.data'
!m 	 dada_mbar      = sortie(1:lenfnam)//'.dada'
!m	data_mbar_std  = sortie(1:lenfnam)//'.data_std'

	posfinal = sortie(1:lenfnam)//'.cin'

!m	kappaF = sortie(1:lenfnam)//'.corfunc_WR'
!m	kappaFd = sortie(1:lenfnam)//'.corfunc_ST'
	
	tequilib=(dt)*teq
	write(*,*) 'No of teq steps', teq
	write(*,*) 'tquilib',tequilib
	


end subroutine read_sundae
