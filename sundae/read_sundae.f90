subroutine read_sundae()
  use T_kind_param_m, ONLY : double
  use gen_com_m,      ONLY : fnam, lenfnam
  use sundae_module,  ONLY : h_A_max, h_ba_min, h_ba_max, h_ba, h_ba_I, fnamtin, totiter,     &
  							 dt, Totalmcmoves, TotalTime, gamma_sundae, Temperature,          &
  							 alpha_max, teq, delta_x, a_sto, kapa, Nbclones, KtoERG,          &
  							 sortie, continue_sundae, maxvec, tprimo, Nbclones_mbar, depart_boucle_nbclones
implicit none

namelist /input_sundae/Totalmcmoves,TotalTime,dt,Nbclones,Temperature,gamma_sundae,sortie,alpha_max, teq, delta_x, a_sto, kapa, continue_sundae, Nbclones_mbar, depart_boucle_nbclones,tprimo,maxvec,h_A_max,h_ba_min,h_ba_max,h_ba,h_ba_I



 h_A_max  = 4.5d-1
 h_ba_max = 1.6d0
 h_ba_min = 1.0d0
 h_ba     = 1.3d0
 h_ba_I   = 1.8d0

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

write(6,*)'NbClones = le numero du canal : ' , NbClones    ! Number of the clone (channel)
write(6,*)'temperature = '  , temperature  ! temperature (kT)
!   
temperature=temperature*KtoERG
!
write(6,*)'friction*dt = '  , gamma_sundae  ! friction*dt
! interval between 2 data records
write(6,*)'sortie data_mbar = ',sortie  ! Name of the file where data are stored
write(6,*)'alpha max = ', alpha_max 
write(6,*)'t_equilib', teq
write(6,*)'delta_X = ', delta_x

a_sto = 1.d0-2.d0*((1.d1**(-2.d0-2.d0*dble(nbclones)/dble(nbclones_mbar))))

write(6,*)'a_sto = ', a_sto
write(6,*)'k ressort = ', kapa
write(6,*)'continue_sundae = ', continue_sundae
write(6,*)'nbclones_mbar = ', Nbclones_mbar
write(6,*)'cluster? 0 no, nbclones yes'
write (*,*) ' depart_boucle_nbclones',depart_boucle_nbclones
write(6,*)'tprimo = ', tprimo
write(*,*) 'maxvec for Lanczos = ', maxvec







end subroutine read_sundae
