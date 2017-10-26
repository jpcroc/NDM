!----------------------------------------------------------------------------------------------------
!-------SUBROUTINEs contained: Free_energy_ABF(), create_files() Deconvolution
!-------This program calculate Free energy (Free energy mollified)files, and create mean force files
!-------This program has to be launched at the end of the simulation--------------------------------
!----------------------------------------------------------------------------------------------------
! A_ee       is A unnormalized
! A_bar_ee   is A_ee normalized corrected
! Free_energy is the free energy 
subroutine Free_energy_ABF()
 USE T_kind_param_m, ONLY:  double
 USE gen_com_m, ONLY: erg2eV,A2cm
 USE tab_imm_m
 USE mab_in_ndm_module, ONLY: delta_z,nhisto,nhisto1,nhisto2, & 
                              histo,histo1,histo2,Free_energy,&
                              mean_force1,abf_type,abf_mode,x_mol,temperature,&
                              A_ee,exp_A_bar,A_bar_ee,eta_ABFee,xi_min, &
                              KtoERG,equit,unit_histo2,  rangmab

implicit none
integer :: ii
integer::i_loop,i_iter
real(double)::Free_temp(-nhisto2:nhisto+nhisto2)
real(double)::renorm_alch(-nhisto2:nhisto+nhisto2)
real(double)::renorm_f,sum_histo
Free_temp(:)=0.d0
renorm_alch(:)=0.d0
if (abf_type .NE. 5) then ! pour ABFee, on va calculer autrement l'energie libre 
  Free_energy(-nhisto1)=0.d0
   do i_loop=-nhisto1+1,nhisto+nhisto1
    !
    Free_energy(i_loop)=Free_energy(i_loop-1)+0.5d0*delta_z*(mean_force1(i_loop-1)+mean_force1(i_loop))
   enddo

   forall(i_loop=-nhisto1:nhisto+nhisto1) Free_temp(i_loop)=exp(-Free_energy(i_loop)/temperature)
   !page 5 in my notes ...
    renorm_f=temperature*log(sum(Free_temp)*delta_z)
   !renormalise par rapport a l'aire
   forall(i_loop=-nhisto1:nhisto+nhisto1) Free_energy(i_loop)=Free_energy(i_loop)+renorm_f
endif
 

!  call set_the_unit_histo()
 select case (abf_type)
 case(1)

  write(*,*)'Free energy computation....Langevin Dynamics'
   open(unit=992,file='Free_energy_Langevin',status='unknown')
    do i_loop=-nhisto1+1,nhisto+nhisto1
     write(992,*) unit_histo2(i_loop),Free_energy(i_loop)*erg2eV
    enddo
   close(992)

 case(2)

    write(*,*)'Free energy computation....ABF BIN'
    open(unit=993,file='Free_energy_ABFBIN',status='unknown')
   
    if (.not.(abf_mode==22)) then
      do i_loop=-nhisto1+1,nhisto+nhisto1
        write(993,*)  unit_histo2(i_loop), Free_energy(i_loop)*erg2eV
      enddo
    else 

      do i_loop=-nhisto1+1,nhisto+nhisto1
        write(993,*)  1.d0/(unit_histo2(i_loop)/(temperature/KtoERG)), &
             (Free_energy(i_loop)/unit_histo2(i_loop)+renorm_f/unit_histo2(i_loop))*erg2eV
      enddo
    end if
    close(993)


 case(3)


    write(*,*)'Free energy computation....ABF BIN OMEGA'
    open(unit=994,file='Free_energy_ABFBIN_OMEGA',status='unknown')
    do i_loop=-nhisto1+1,nhisto+nhisto1
      write(994,*)  unit_histo2(i_loop), Free_energy(i_loop)*erg2eV
    enddo
    close(994)

 case(4)
  
    write(*,*)'Free energy computation....ABF Gaussian'
    open(unit=995,file='Free_energy_ABFGaussian',status='unknown')
    do i_loop=-nhisto1+1,nhisto+nhisto1
      write(995,*)  unit_histo2(i_loop), Free_energy(i_loop)*erg2eV
    enddo
    close(995)
 case(5,8)
  if (rangmab==0) write(*,*)'Free energy computation....ABF ee'
  !--------Pour ABFee, on a besoin des A_ee pour calculer l'energie libre.
  !-----------------------A_tilde-------------------------------
  !-----------------------------------------------------------------

  ! take the remor_f from A_ee
  forall(i_loop=-nhisto1:nhisto+nhisto1) Free_temp(i_loop)=exp(-A_ee(i_loop)/temperature)
  renorm_f=temperature*log(sum(Free_temp(-nhisto1:nhisto+nhisto1))*delta_z)


  if (abf_mode==2) then
    forall(i_loop=-nhisto1:nhisto+nhisto1) Free_energy(i_loop)=A_ee(i_loop) + renorm_f
    !do ii=-nhisto1,nhisto+nhisto1
 !
  !  write(*,*) ii, Free_energy(ii)
   ! end do
    !stop 
  end if 

  if (abf_mode==22) then
  
   do i_loop=-nhisto1,nhisto+nhisto1
   !renormalization by changing the temperature see the page 17 of my notes ...
   ! this was the old version from abf_type=22
   ! do i_loop= -nhisto1,nhisto+nhisto1
   !  do i_iter= -nhisto1,nhisto+nhisto1
   !   renorm_alch(i_loop)=renorm_alch(i_loop) + exp(-A_ee(i_iter)*unit_histo2(i_iter)/temperature)*delta_z
   !  end do
   !   renorm_alch(i_loop)=temperature*log(renorm_alch(i_loop))
   ! end do
   !
   ! 
   !       Free_energy(i_loop) = A_ee(i_loop)/unit_histo2(i_loop)           &
   !                     -renorm_alch(i_loop)/unit_histo2(i_loop)   &
   !                    +renorm_f/unit_histo2(i_loop)
 
   ! the temperature is alreaby puted in the renormalization constant. 
   Free_energy(i_loop)=A_ee(i_loop)/unit_histo2(i_loop)+renorm_f/unit_histo2(i_loop) 
   end do    
  end if  !abf_mode=22


! why this is here !!!!!
  forall(i_loop=-nhisto1:nhisto+nhisto1) A_ee(i_loop)=A_ee(i_loop)+renorm_f
  !minfreeval=minval(A_ee(-nhisto1:nhisto+nhisto1))
  !forall(i_loop=-nhisto1:nhisto+nhisto1) A_ee(i_loop)=A_ee(i_loop)-minfreeval
  open(unit=996,file='Free_energy_mollifiee_ABFee',status='unknown')!A_tilde
    ! 
   do i_loop=-nhisto1, nhisto+nhisto1
       write(996,*) unit_histo2(i_loop),A_ee(i_loop)*erg2eV
   enddo
    !
  close(996)
  if (abf_mode==1) then 
   ! THOSE LINES ARE OBSOLTE: mcmCHECK
     forall(i_loop=-nhisto1:nhisto+nhisto1) Free_temp(i_loop)=exp(-Free_energy(i_loop)/temperature)
      renorm_f=temperature*log(sum(Free_temp)*delta_z)
     forall(i_loop=-nhisto1:nhisto+nhisto1) Free_energy(i_loop)=Free_energy(i_loop)+renorm_f
   !END OF OBSOLETE LINES
!-----------------------------------------------------------------------------------
!--------------------------A_bar----------------------------------------------------
!------------------------------------------------------------------------------------ 

      !
      do i_loop=-nhisto1,nhisto+nhisto1
        exp_A_bar(i_loop)=0.d0
       do i_iter=-nhisto1,nhisto+nhisto1
         exp_A_bar(i_loop)=exp_A_bar(i_loop)  +  &
         exp(-(x_mol(i_loop)-x_mol(i_iter))**2/   &
        (2.d0*eta_ABFee*temperature)+A_ee(i_iter)/temperature)*delta_z
       enddo
      enddo
      renorm_f=sum(exp_A_bar(-nhisto1:nhisto+nhisto1))*delta_z
      forall(i_loop=-nhisto1:nhisto+nhisto1) A_bar_ee(i_loop)=log(exp_A_bar(i_loop)/renorm_f)*temperature
      forall(i_loop=-nhisto1:nhisto+nhisto1) Free_temp(i_loop)=exp(-A_bar_ee(i_loop)/temperature)
      renorm_f=temperature*log(sum(Free_temp)*delta_z)
      A_bar_ee(:)=A_bar_ee(:)+renorm_f
!forall(i_loop=-nhisto1:nhisto+nhisto1) A_bar_ee(i_loop)=A_bar_ee(i_loop)-minfreeval


      open(unit=998,file='Free_energy_biais_ABFee',status='unknown')
      do i_loop=-nhisto1,nhisto+nhisto1
        write(998,*) x_mol(i_loop)/A2cm,A_bar_ee(i_loop)*erg2eV
      enddo
      close(998)
!-----------------------------------------------------------------------
!-----------------------Free_energy!!!!---------------------------------
!-----------------------------------------------------------------------
 
     sum_histo=sum(histo1)*delta_z
     forall(i_loop=-nhisto1:nhisto+nhisto1)  Free_energy(i_loop)=A_bar_ee(i_loop)-temperature*log(histo1(i_loop)/sum_histo)
     forall(i_loop=-nhisto1:nhisto+nhisto1) Free_temp(i_loop)=exp(-Free_energy(i_loop)/temperature)
     renorm_f=temperature*log(sum(Free_temp)*delta_z)
     forall(i_loop=-nhisto1:nhisto+nhisto1) Free_energy(i_loop)=Free_energy(i_loop)+renorm_f
 end if  !abf_mode==1



     open(unit=79,file='Free_energy_ABFee',status='unknown')

     if (.NOT.(abf_mode==22)) then
       do i_loop=-nhisto1,nhisto+nhisto1
          write(79,*)unit_histo2(i_loop),Free_energy(i_loop)*erg2eV
       enddo
     end if 

     if (abf_mode==22) then
      do i_loop=nhisto+nhisto1,-nhisto1,-1
        !write(999,*) temperature/KtoERG/unit_histo2(i_loop),Free_energy(i_loop)*erg2eV+equit*erg2ev
        write(79,*) 1.d0/(unit_histo2(i_loop)/(temperature/KtoERG)),Free_energy(i_loop)*erg2eV+equit*erg2ev
      enddo
     end if 
  
     close(79, status='keep')


end select


end subroutine Free_energy_ABF


subroutine create_files()!!---------Create files for the programm


  use mpi
  use mod_mpi_mab

 USE T_kind_param_m, ONLY:  double
 USE gen_com_m, ONLY: erg2eV,A2cm
 USE tab_imm_m
 USE mab_in_ndm_module, ONLY: nhisto,nhisto1,nhisto2,       & 
                              histo,histo1,histo2,delta_z,  &
                              histo_equi,abf_type, &
                              A_ee,A_bar_ee,histo_zeta,x_mol,&
                              mean_force1,mean_force_ABFee_dyn,abf_type,&
                              Free_energy,temperature,rangmab, unit_histo2



implicit none
integer::i_iter,i_loop
real(double)::sum_histo1,sum_histo_zeta,fnorm1,fnorm2,temp_T,temp_A_bar,temp_A
real(double),dimension(:),allocatable::Inter_meanforce,Inter_meanforce1
real(double),dimension(:),allocatable::A_bar_corrige, A_corrige
real(double)::Free_temp(-nhisto1:nhisto+nhisto1)
allocate(Inter_meanforce(-nhisto1:nhisto+nhisto1), Inter_meanforce1(-nhisto1:nhisto+nhisto1))
allocate(A_bar_corrige(-nhisto1:nhisto+nhisto1),A_corrige(-nhisto1:nhisto+nhisto1))

temp_T=temperature*erg2eV
Free_temp(:)=0
Inter_meanforce(:)=0
Inter_meanforce1(:)=0


!---------data_A_bar: \xi(q), A_bar_ee, A_bar_ee corrected by substracting log(P(\xi))
!---------data_A_tilde: zeta, A_ee, histo_zeta, A_ee corrected by substracting log(P(\zeta))

sum_histo1=sum(histo1)
histo1(:)=histo1(:)/sum_histo1

if ((abf_type==5) .or. (abf_type == 8)) then
  sum_histo_zeta=sum(histo_zeta)
  histo_zeta(:)=histo_zeta(:)/sum_histo_zeta

  A_bar_corrige(:)=A_bar_ee(-nhisto1:nhisto+nhisto1)*erg2eV-log(histo1(:))*temperature*erg2eV
  !
  forall(i_loop=-nhisto1:nhisto+nhisto1) Free_temp(i_loop)=exp(-A_bar_corrige(i_loop)/temperature)

  temp_A_bar=temperature*log(sum(Free_temp)*delta_z)
 
  if (rangmab==0) write(6,'("MAB: create files temp_A_bar......:",f15.6)') temp_A_bar
  A_bar_corrige(:)=A_bar_corrige(:)+temp_A_bar
  
  A_corrige(:)=A_ee(-nhisto1:nhisto+nhisto1)*erg2eV-log(histo_zeta(-nhisto1:nhisto+nhisto1))*temperature*erg2eV


  forall(i_loop=-nhisto1:nhisto+nhisto1) Free_temp(i_loop)=exp(-A_corrige(i_loop)/temperature)

  temp_A=temperature*log(sum(Free_temp)*delta_z)

  if (rangmab==0) write(6,'("MAB: create_files temp_A .........:",f15.6)') temp_A

  A_corrige(:)=A_corrige(:)+temp_A
  call mpi_barrier (MPI_COMM_WORLD,codeph)
  !r if (rangmab==0) then
  !----1. bucket of \xi. 2. bar A, 3. A_bar_ee corrected by substracting log(P(\xi))
  open(unit=968,file='data_A_bar',status='unknown')

  !----1. bucket of \zeta. 2.tilde A, 3. A_ee corrected by substracting log(P(\zeta))
  open(unit=967,file='data_A_tilde',status='unknown')
   do i_iter=-nhisto1,nhisto+nhisto1
     write(968,*) unit_histo2(i_iter), A_bar_ee(i_iter)*erg2eV, A_bar_corrige(i_iter)
     write(967,*) unit_histo2(i_iter), A_ee(i_iter)*erg2ev, A_corrige(i_iter)
  enddo
  !
  close(968)
  close(967)
  !r end if
else
  A_corrige(:)=Free_energy(-nhisto1:nhisto+nhisto1)*erg2eV-log(histo1(:))*temperature*erg2eV


  forall(i_loop=-nhisto1:nhisto+nhisto1) Free_temp(i_loop)=exp(-A_corrige(i_loop)/temperature)
 
  temp_A=temperature*log(sum(Free_temp)*delta_z)

  if (rangmab==0) write(6,'("MAB: create files temp_A",f15.6)') temp_A
   A_corrige(:)=A_corrige(:)+temp_A

   !r if (rangmab==0) then 
   !----1. bucket of \xi. 2. Free energy, 3. Free energy corrected by substracting log(P(\xi))
   open(unit=966,file='data_A',status='unknown')

   do i_iter=-nhisto1,nhisto+nhisto1
     write(966,*) unit_histo2(i_iter), Free_energy(i_iter)*erg2eV, A_corrige(i_iter)
   enddo
   close(966)
   !r end if
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--------mean force and integration of meanforce--------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if ((abf_type==5).or.(abf_type==8)) then
   ! meanforce1: local force in ABFee 
   ! mean_force_ABFee_dyn (old mean_force2): mean force used in the dynamics, derivative of bar_A in ABFee
   !r if (rangmab==0) then
   !meanforce file - mean force used in dynamics, derivative of bar_A
   open(unit=991,file='meanforce',status='unknown')
   do i_iter=-nhisto2,nhisto+nhisto2
       write(991,*) i_iter,mean_force_ABFee_dyn(i_iter)
   end do
   close(991)
   !r end if

   do i_iter=-nhisto1+1,nhisto+nhisto1
     Inter_meanforce(i_iter)=Inter_meanforce(i_iter-1)+0.5d0*delta_z*(mean_force_ABFee_dyn(i_iter-1)+mean_force_ABFee_dyn(i_iter))
   enddo
   forall(i_loop=-nhisto1:nhisto+nhisto1) Free_temp(i_loop)=exp(-Inter_meanforce(i_loop)/temperature)
   fnorm1=temperature*log(sum(Free_temp)*delta_z)
   Inter_meanforce(:)=Inter_meanforce(:)+fnorm1

   !r if (rangmab==0) then
   ! Integration_meanforce file - integration of mean force used in dynamics mean_force_ABFee_dyn (old mean_force2 for ABFee)
   open(unit=75,file='Integration_meanforce',status='unknown')
   do i_iter=-nhisto1,nhisto+nhisto1
       write(75,'(f18.4,d18.10,f18.4)')unit_histo2(i_iter), Inter_meanforce(i_iter)*erg2eV, histo_zeta(i_iter)
   enddo
   close(75, status='keep')

   ! data_force_ABFee: a collection of mean forces with 
   ! meanforce1: local force in ABFee 
   ! mean_force_ABFee_dyn (old meanforce2): mean force used in the dynamics in ABFee

   open(unit=964,file='data_force_ABFee',status='unknown')
   do i_iter=-nhisto1, nhisto+nhisto1 
          ! data_force_ABFee file
          write(964,*)i_iter,mean_force1(i_iter),mean_force_ABFee_dyn(i_iter)   
   end do
   close(964)

   !r endif !rangmab

   !data_Integration_force_ABFee file
   ! a collection of \zeta 
   ! integration of local force in ABFee (meanforce1) which should be equal to bar_A
   ! A_bar_ee (\zeta)  
   ! the column 2 and 3 should be very close is the simulation is converged. 

   do i_iter=-nhisto1+1,nhisto+nhisto1
     Inter_meanforce1(i_iter)=Inter_meanforce1(i_iter-1)+0.5d0*delta_z*(mean_force1(i_iter-1)+mean_force1(i_iter))
   enddo
   forall(i_loop=-nhisto1:nhisto+nhisto1) Free_temp(i_loop)=exp(-Inter_meanforce1(i_loop)/temperature)
   fnorm2=temperature*log(sum(Free_temp)*delta_z)
   Inter_meanforce1(:)=Inter_meanforce1(:)+fnorm2


   !r if (rangmab==0) then
   open(unit=965,file='data_Integration_force_ABFee',status='unknown')
   do i_iter=-nhisto1,nhisto+nhisto1
       write(965,*) unit_histo2(i_iter), Inter_meanforce1(i_iter)*erg2eV,A_bar_ee(i_iter)*erg2eV
   enddo
   close(965)
   !r end if

else 
   !   meanforce1 - mean force used in dynamics (is different that in the case of abf_type=5) 

   !  meanforce file 
   ! contains the meanforce used in dynamics  (meanforce1) for ABFbin 


   !r if (rangmab==0) then
   open(unit=991,file='meanforce',status='unknown')
   do i_iter=-nhisto1,nhisto+nhisto1
    write(991,*) i_iter,mean_force1(i_iter)
   enddo
   close(991)
   !r end if

   do i_iter=-nhisto1+1,nhisto+nhisto1
      Inter_meanforce1(i_iter)=Inter_meanforce1(i_iter-1)+0.5d0*delta_z*(mean_force1(i_iter-1)+mean_force1(i_iter))
   enddo

   forall(i_loop=-nhisto1:nhisto+nhisto1) Free_temp(i_loop)=exp(-Inter_meanforce1(i_loop)/temperature)
   fnorm1=temperature*log(sum(Free_temp)*delta_z)
   Inter_meanforce1(:)=Inter_meanforce1(:)+fnorm1



   !r if (rangmab==0) then
   ! Intr_meanforce1 - integration of the meanforce used in dynamics for ABFbin
   ! Integration_meanforce file 
   open(unit=966,file='Integration_meanforce',status='unknown')


   do i_iter=-nhisto1,nhisto+nhisto1
       write(966,*)unit_histo2(i_iter), Inter_meanforce1(i_iter)*erg2eV, histo1(i_iter)
   enddo
   close(966)
   !r end  if
end if 


return
end subroutine create_files


subroutine deconvolution_ABFee()! Deconvolution of A_ee, it reads the file 'data_A_tilde', and writes out Deconvo_free_energy

 USE T_kind_param_m, ONLY:  double
 USE gen_com_m, ONLY: erg2eV,A2cm
 USE tab_imm_m
 USE mab_in_ndm_module, ONLY: nhisto,nhisto1,nhisto2,       & 
                              histo,histo1,histo2,delta_z,  &
                              histo_zeta,x_mol,&
                              mean_force1,abf_type,&
                              Free_energy,temperature,nom_deconvo,&
                              eta_ABFee



implicit none

real(double)::lest1,lest2
real(double),dimension(:),allocatable::A_tilde,A_tilde_temp,Deconvo_exp_pre
real(double),dimension(:),allocatable::Deconvo_denom,Deconvo_inter,Deconvo_exp,Deconvo_exp_temp
integer::i_iter,iter,optime
real(double)::temp_calcul,C_norm
logical :: dir_at
allocate(A_tilde(-nhisto1:nhisto+nhisto1),A_tilde_temp(-nhisto1:nhisto+nhisto1))
allocate(Deconvo_exp_pre(-nhisto1:nhisto+nhisto1),Deconvo_denom(-nhisto1:nhisto+nhisto1))
allocate(Deconvo_inter(-nhisto1:nhisto+nhisto1),Deconvo_exp(-nhisto1:nhisto+nhisto1))
allocate(Deconvo_exp_temp(-nhisto1:nhisto+nhisto1))

Deconvo_inter(:)=0

A_tilde(:)=0
A_tilde_temp(:)=0


inquire(file='data_A_tilde',exist=dir_at)
if (.not.dir_at) then
  write(6,*) 'MAB: The file data_A_tilde is not present.The program will stop'
  write(6,*) 'MAB: stop in sub deconvolution_ABFee'
  stop
end if 
open(unit=967,file='data_A_tilde',action='read')

do i_iter=-nhisto1,nhisto+nhisto1
read(967,*)lest1,A_tilde_temp(i_iter),lest2
enddo

A_tilde_temp(:)=A_tilde_temp(:)/erg2eV ! change unit for calculating

A_tilde(:)=dexp(-A_tilde_temp(:)/temperature)!

!------Renormalise A_tilde (:)

C_norm=sum(A_tilde(:)*delta_z)

A_tilde(:)=A_tilde(:)/C_norm

Deconvo_exp_pre(:)=A_tilde(:)
!--------begin-----------

do optime=1,nom_deconvo
 
Deconvo_exp(:)=0
Deconvo_denom(:)=0
Deconvo_inter(:)=0
Deconvo_exp_temp(:)=0

  do i_iter=-nhisto1,nhisto+nhisto1
    
    do iter=-nhisto1,nhisto+nhisto1
       temp_calcul=dexp(-(x_mol(iter)-x_mol(i_iter))**2/(2.d0*eta_ABFee*temperature))*Deconvo_exp_pre(iter)*delta_z
       Deconvo_denom(i_iter)=Deconvo_denom(i_iter)+temp_calcul
    enddo
    
    if (Deconvo_denom(i_iter)==0.d0) then
    write(*,*)'error!!! Deconvo_denom is zero!'
    stop
    endif 
    
  enddo
  !------normaliser Deconvo_denom 
   

  Deconvo_inter(:)=A_tilde(:)/Deconvo_denom(:)

  do iter=-nhisto1,nhisto+nhisto1
     
     do i_iter=-nhisto1,nhisto+nhisto1
        temp_calcul=dexp(-(x_mol(iter)-x_mol(i_iter))**2/(2.d0*eta_ABFee*temperature))*Deconvo_inter(i_iter)*delta_z
        Deconvo_exp_temp(iter)=Deconvo_exp_temp(iter)+temp_calcul  
     enddo
     if (Deconvo_exp_temp(iter) .eq. 0.d0) then
     write(*,*)'error!!!Deconvo_exp_temp is zero!'
     stop
     endif

  enddo

 

  Deconvo_exp(:)=Deconvo_exp_temp(:)*Deconvo_exp_pre(:)

C_norm=sum(Deconvo_exp(:)*delta_z)

Deconvo_exp=Deconvo_exp/C_norm


Deconvo_exp_pre(:)=Deconvo_exp(:)

enddo

Deconvo_exp(:)=-log(Deconvo_exp_pre(:))*temperature*erg2eV


open(unit=1111,file='Deconvo_Free_energy',status='unknown')

do iter=-nhisto1,nhisto+nhisto1

write(1111,*) x_mol(iter)/A2cm,Deconvo_exp(iter)
enddo





end subroutine deconvolution_ABFee



subroutine on_run_writting_histo_freetemp(itapp)
 USE T_kind_param_m, ONLY:  double
 USE gen_com_m, ONLY: erg2eV,A2cm,lenfnam,fnam
 USE tab_imm_m
 USE mab_in_ndm_module, ONLY : A_ee, &
                              P_ee_num, P_ee_denom, &
                              delta_z, & 
                              unit_histo2,nhisto,nhisto1,nhisto2, &
                              abf_type, abf_mode,        &
                              histo, histo1, histo_zeta, &
                              temperature, KtoErg, &
                              abf_mode_reaction,abf_mode_temperature, abf_mode_alchemical
 implicit none
 integer, intent(in) :: itapp
! local variables .... 
 integer      :: i_iter
 real(double) :: renorm_f
 real(double)::Free_temp(-nhisto2:nhisto+nhisto2)
 character :: extension*9, fnamhisto0*80, fnamhisto1*80
 character :: fnamhisto_zeta*80,fnamfreetemp*80,fnamaeerestart*80
 character(len=1) :: diez

 diez=char(35)
 write(extension,'(i9.9)') itapp
 fnamhisto0=fnam(1:lenfnam)//'.histogram0.'//extension
 fnamhisto1=fnam(1:lenfnam)//'.histogram1.'//extension
if (.NOT.(abf_type==5)) then
    open(unit=989,file=fnamhisto0,status='unknown')
    open(unit=990,file=fnamhisto1,status='unknown')

     do i_iter=-nhisto1,nhisto+nhisto1
      write(989,*)i_iter, histo1(i_iter)
     enddo
     do i_iter=1,nhisto
       write(990,*)i_iter, histo(i_iter)
     enddo

    close (989)
    close (990)
   !

end if 

if (abf_type==5) then
   !
    fnamhisto_zeta=fnam(1:lenfnam)//'.temp_run_histogram_zeta.'//extension
    open(unit=990,file=fnamhisto_zeta,status='unknown')
    if ((abf_mode==abf_mode_reaction).or.(abf_mode==abf_mode_alchemical)) then
      do i_iter=-nhisto1,nhisto+nhisto1
       write(990,'(i9,E25.12)') i_iter, histo_zeta(i_iter)
      enddo
    end if 

    if (abf_mode==abf_mode_temperature) then
     do i_iter=-nhisto1,nhisto+nhisto1
       write(990,'(i9,E25.12,2f15.7)') i_iter, histo_zeta(i_iter), (temperature/KtoErg)/unit_histo2(i_iter),unit_histo2(i_iter)
     enddo
    end if 
    close(990)
 
    fnamaeerestart=fnam(1:lenfnam)//'.temp_run_a_ee_restartin.'//extension
    open(unit=991,file=fnamaeerestart,status='unknown')
    if ((abf_mode==abf_mode_reaction).or.(abf_mode==abf_mode_alchemical)) then
      !writting the en-tete ... in order to have a vague trace of the simulation
      write(991,'(a,"  ",2i3,i7,3i5,E25.12E3)') diez, abf_mode, abf_type, itapp, nhisto, nhisto1, nhisto2, delta_z   
      do i_iter=-nhisto2,nhisto+nhisto2
        write(991,'(E25.10,3E25.15)') dble(i_iter)*delta_z, A_ee(i_iter)*erg2ev, P_ee_num(i_iter),P_ee_denom(i_iter)
      end do
      close(991)
     end if 

    fnamfreetemp  =fnam(1:lenfnam)//'.temp_run_free_energyyyy.'//extension
    open(unit=990,file=fnamfreetemp,status='unknown')
    if ((abf_mode==abf_mode_reaction).or.(abf_mode==abf_mode_alchemical)) then
      forall(i_iter=-nhisto1:nhisto+nhisto1) Free_temp(i_iter)=exp(-A_ee(i_iter)/temperature)
      renorm_f=temperature*log(sum(Free_temp(-nhisto1:nhisto+nhisto1))*delta_z)
      do i_iter=-nhisto1+1,nhisto+nhisto1
        write(990,'(E25.10,E25.15)') i_iter*delta_z, (A_ee(i_iter)+renorm_f)*erg2ev
      end do
     end if 

     if (abf_mode==abf_mode_temperature) then
      do i_iter=-nhisto1+1,nhisto+nhisto1
         write(990,'(f14.5,3f17.5)')  (temperature/KtoErg)/unit_histo2(i_iter),   &
          (A_ee(i_iter)/unit_histo2(i_iter)+renorm_f/unit_histo2(i_iter))*erg2ev, &
          A_ee(i_iter)/unit_histo2(i_iter)*erg2ev, &
          renorm_f/unit_histo2(i_iter)*erg2ev
      end do
      close(990)
     end if
    
end if 

return

end subroutine on_run_writting_histo_freetemp





