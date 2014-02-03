subroutine calfo_mab()
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use tab_imm_m
  USE mab_in_ndm_module, only: it_mab,abf_type,block,histo_equi,n_equilibre
  implicit none
  integer :: it_langevin
 
  it_langevin=it_mab
! one force calculation ....

! NDM part ...

         if (itab/=0) then
          if (mod(it_langevin,itab)==0) then
           call caltabt
          endif
         endif
         if (ltabvois.and.mod(it_langevin,itetabvois)==0) call caltabi
        call calfo
        if (block) call calfoblock()

!ABF part ...

        call fill_histo()

!---if we only want to fill histogram after n_equilibre steps-------------
        if ((histo_equi == .true.) .And. (it_mab >= n_equilibre)) then
        call fill_histo_equilibre()
        endif
!--------------------------------------------------------------------------

      
    select case (abf_type)
           case (1)
                  continue
           case (2)  
                  call calfo_ABF_BIN
           case (3) 
                  call calfo_ABF_BIN_OMEGA
          case (4)  
                  call calfo_ABF_Gaussien ! ABF Gaussian
          case (5) 
                  call calfo_ABFee
          case (6)
                  call calfo_ABFee_const_biais
          case (7)
                  call calfo_ABF_BIN_const_biais
          case (8)
                  call calfo_ABFee_iter
        end select 

return
end subroutine calfo_mab



subroutine calfo_ABF_BIN()
 USE T_kind_param_m, ONLY:  double
 USE gen_com_m, ONLY: zero,im,imm,low_limit,angst,ev2erg,erg2ev
 USE tab_imm_m
 USE mab_in_ndm_module, ONLY:dcsi,icsi,rfilac,histo,     &
                             mean_force,cumul_force1,nhisto,nhisto1, &
                             mean_force1,histo1
 implicit none

 real(double), dimension(3,imm) :: fpabf
 real(double) :: force

 fpabf(:,:) = zero
! Computing the forces from the ABF bins ...
force = - DOT_PRODUCT(fp(:,7),rfilac(:))

if ((icsi >= -nhisto1).and.(icsi <= nhisto+nhisto1)) then 
 cumul_force1(icsi) =  cumul_force1(icsi) + force
mean_force1 (icsi) = cumul_force1(icsi)/histo1(icsi)
fpabf(1:3,7)=rfilac(1:3)*mean_force1(icsi)

end if

! Updating the forces ...
!write (41,*) fp(1,7),fpabf(1,7)
!write (42,*) fp(2,7),fpabf(2,7)
!write (43,*) fp(3,7),fpabf(3,7)

  fp(1:3,:)=fp(1:3,:)+fpabf(1:3,:)
  

return
end subroutine calfo_ABF_BIN


subroutine calfo_ABF_BIN_const_biais()! This subroutine fixes the biais of ABF bin, then launch the dynamics and save histogram files
 USE T_kind_param_m, ONLY:  double
 USE gen_com_m, ONLY: zero,im,imm,low_limit,angst,ev2erg,erg2ev
 USE tab_imm_m
 USE mab_in_ndm_module, ONLY:dcsi,icsi,rfilac,histo,     &
                             mean_force,cumul_force1,nhisto,nhisto1, &
                             mean_force1,histo1
 implicit none

 real(double), dimension(3,imm) :: fpabf
 real(double) :: force

 fpabf(:,:) = zero

! Computing the forces from the ABF bins ...
force = - DOT_PRODUCT(fp(:,7),rfilac(:))

if ((icsi >= -nhisto1).and.(icsi <= nhisto+nhisto1)) then 
fpabf(1:3,7)=rfilac(1:3)*mean_force1(icsi)
end if

! Updating the forces ...
!write (41,*) fp(1,7),fpabf(1,7)
!write (42,*) fp(2,7),fpabf(2,7)
!write (43,*) fp(3,7),fpabf(3,7)

  fp(1:3,:)=fp(1:3,:)+fpabf(1:3,:)
  

return
end subroutine calfo_ABF_BIN_const_biais

subroutine calfo_ABF_BIN_OMEGA()
 USE T_kind_param_m, ONLY:  double
 USE gen_com_m, ONLY: zero,im,imm,low_limit,angst,ev2erg,erg2ev
 USE tab_imm_m
 USE mab_in_ndm_module, ONLY:dcsi,icsi,rfilac,histo,     &
                             mean_force,cumul_force1,nhisto,nhisto1, &
                             mean_force1,histo1,omega_abf,&
                             a_Fermi,xi_min,xi_max
  
 real(double), dimension(3,imm) :: fpabf
 real(double) :: force,omega
 real(double) :: Fermi_manuel

 
 fpabf(:,:) = zero
 ! Computing the forces from the ABF bins ...
 force = - DOT_PRODUCT(fp(:,7),rfilac(:))

 if ((icsi >= -nhisto1).and.(icsi <= nhisto+nhisto1)) then 
  cumul_force1(icsi) =  cumul_force1(icsi) + force
  mean_force1(icsi) = cumul_force1(icsi)/(1.d0/omega_abf+histo1(icsi))
  fpabf(1:3,7)=rfilac(1:3)*mean_force1(icsi)
 end if

! Updating the forces ...
!write (41,*) fp(1,7),fpabf(1,7)
!write (42,*) fp(2,7),fpabf(2,7)
!write (43,*) fp(3,7),fpabf(3,7)

  fp(1:3,:)=fp(1:3,:)+fpabf(1:3,:)
  
return
end subroutine calfo_ABF_BIN_OMEGA


subroutine calfo_ABF_Gaussien()
 USE T_kind_param_m, ONLY:  double
 USE gen_com_m, ONLY: zero,im,imm,low_limit,angst,ev2erg,erg2ev
 USE tab_imm_m
 USE mab_in_ndm_module, ONLY:dcsi,icsi,rfilac,histo,     &
                             mean_force,cumul_force1,nhisto,nhisto1, &
                             mean_force1,histo1,omega_abf,&
                             a_Fermi,xi_min,xi_max,ecart_eta,eta_mab,&
                             cumul_force_denom1,x_mol,sigma_carre

 real(double), dimension(3,imm) :: fpabf
 real(double) :: force,omega
 real(double) :: Fermi_manuel
 integer:: indice_gaussian

 
 fpabf(:,:) = zero
 ! Computing the forces from the ABF bins ...
 force = - DOT_PRODUCT(fp(:,7),rfilac(:))
do indice_gaussian=icsi-ecart_eta,icsi+ecart_eta
 if ((indice_gaussian >=-nhisto1).and.(indice_gaussien <= nhisto+nhisto1)) then 
  cumul_force1(indice_gaussian)=cumul_force1(indice_gaussian) +force*dexp(-(dcsi-x_mol(indice_gaussian))**2/2.d0*sigma_carre)

  cumul_force_denom1(indice_gaussian)=cumul_force_denom1(indice_gaussian)+dexp(-(dcsi-x_mol(indice_gaussian))**2/2.d0*sigma_carre)
  mean_force1 (indice_gaussian) = cumul_force1(indice_gaussian)/cumul_force_denom1(indice_gaussian)

  if(indice_gaussian == icsi) then
  fpabf(1:3,7)=rfilac(1:3)*mean_force1(icsi)
  endif
 end if
enddo
! Updating the forces ...
!write (41,*) fp(1,7),fpabf(1,7)
!write (42,*) fp(2,7),fpabf(2,7)
!write (43,*) fp(3,7),fpabf(3,7)

  fp(1:3,:)=fp(1:3,:)+fpabf(1:3,:)
  
return
end subroutine calfo_ABF_Gaussien



subroutine calfo_ABFee()
 USE T_kind_param_m, ONLY:  double
 USE gen_com_m, ONLY: zero,im,imm,low_limit,angst,ev2erg,erg2ev
 USE tab_imm_m
 USE mab_in_ndm_module, ONLY:dcsi,icsi,rfilac,histo,     &
                             mean_force,cumul_force1,nhisto,nhisto1, &
                             nhisto2,mean_force1,histo1,cumul_force1,   &
                             A_dev_ee,A_ee,P_ee,P_ee_num,P_ee_denom,delta_z,&
                             x_mol,eta_ABFee,temperature,omega_abf,exp_A_bar,&
                             mean_force2,it_mab,histo_zeta,n_equilibre,histo_equi,&
                             it_mab

 implicit none
integer::iter
real(double), dimension(3,imm) :: fpabf
real(double),dimension(:),allocatable::temp_log
real(double),dimension(:),allocatable::temp_exp
real(double),dimension(:),allocatable::temp_num_f
real(double)::temp_log_max,force
real(double)::denom,num_f,mean_force_ABF

allocate(temp_log(-nhisto2:nhisto+nhisto2),temp_exp(-nhisto2:nhisto+nhisto2),temp_num_f(-nhisto2:nhisto+nhisto2))
 
!write(*,*),nhisto2,nhisto,nhisto1
temp_log(:)=0
temp_num_f(:)=0!
temp_exp(:)=0


fpabf(:,:) = zero

! Computing the forces as an observable ...
force = - DOT_PRODUCT(fp(:,7),rfilac(:))

if ((icsi >= -nhisto1).and.(icsi <= nhisto+nhisto1)) then 
 cumul_force1(icsi) =  cumul_force1(icsi) + force
mean_force1 (icsi) = cumul_force1(icsi)/histo1(icsi)
end if



!---------BEGIN PRINCIPAL PROGRAM-----------------------------------------------------

!---------1. Compute A_ee(\zeta) for every \zeta, Méthode des trapèzes
A_ee(-nhisto2)=0.d0
do iter=-nhisto2+1,nhisto+nhisto2
A_ee(iter)=A_ee(iter-1)+delta_z*0.5d0*(A_dev_ee(iter-1)+A_dev_ee(iter))
end do

!--------eta= eta_ABFee for extended dynamics-----------
forall(iter=-nhisto2:nhisto+nhisto2) temp_log(iter)=-(x_mol(iter)-dcsi)**2/(2.d0*eta_ABFee*temperature)+A_ee(iter)/temperature


temp_log_max=maxval(temp_log)

forall(iter=-nhisto2:nhisto+nhisto2) temp_log(iter)=temp_log(iter)-temp_log_max

forall(iter=-nhisto2:nhisto+nhisto2) temp_exp(iter)=exp(temp_log(iter))! 

forall(iter=-nhisto2:nhisto+nhisto2) temp_num_f(iter)=temp_exp(iter)*(x_mol(iter)-dcsi)

denom=sum(temp_exp)

num_f=sum(temp_num_f)

mean_force_ABF=1.d0/eta_ABFee*num_f/(denom)

mean_force2(icsi)=mean_force_ABF

fpabf(1:3,7)=rfilac(1:3)*mean_force_ABF

fp(1:3,:)=fp(1:3,:)+fpabf(1:3,:)


!write (41,*) fp(1,7),fpabf(1,7)
!write (42,*) fp(2,7),fpabf(2,7)
!write (43,*) fp(3,7),fpabf(3,7)


do iter=-nhisto2, nhisto+nhisto2
P_ee(iter)=temp_exp(iter)/denom
P_ee_denom(iter)=P_ee_denom(iter)+P_ee(iter)
P_ee_num(iter)=P_ee_num(iter)+(x_mol(iter)-dcsi)*P_ee(iter)
enddo

!----------If histo_equi is true when it_mab > n_equilibre, or histo_equi is false, we fill the histogram of zeta
if(((histo_equi == .true.) .AND. (it_mab > n_equilibre)) .OR. (histo_equi == .False.)) then

do iter= -nhisto1, nhisto+nhisto1 
histo_zeta(iter)=histo_zeta(iter)+P_ee(iter)
enddo

endif

!---------Calcule A_dev if necessaire
do iter=-nhisto1,nhisto+nhisto1
A_dev_ee(iter)=1.d0/eta_ABFee*(P_ee_num(iter))/(P_ee_denom(iter)+1.d0/(omega_abf*dble(nhisto)))
enddo

if(nhisto1 < nhisto2) then
!---------------en dehors de nhisto1 
do iter=-nhisto2,-nhisto1
A_dev_ee(iter)=0.d0
enddo

do iter=nhisto1+nhisto,nhisto+nhisto2
A_dev_ee(iter)=0.d0
enddo
endif

return

end subroutine calfo_ABFee

subroutine calfo_ABFee_const_biais()! This subroutine fixes biais \bar A of ABF ee, then launch the dynamics and save histogram files

 USE T_kind_param_m, ONLY:  double
 USE gen_com_m, ONLY: zero,im,imm,low_limit,angst,ev2erg,erg2ev
 USE tab_imm_m
 USE mab_in_ndm_module, ONLY:dcsi,icsi,rfilac,histo,     &
                             mean_force,cumul_force1,nhisto,nhisto1, &
                             nhisto2,mean_force1,histo1,cumul_force1,   &
                             A_dev_ee,A_ee,P_ee,P_ee_num,P_ee_denom,delta_z,&
                             x_mol,eta_ABFee,temperature,omega_abf,exp_A_bar,&
                             mean_force2,it_mab,histo_zeta,n_equilibre,histo_equi,&
                             it_mab,mean_force_ABFee

 implicit none
integer::iter
real(double), dimension(3,imm) :: fpabf
real(double),dimension(:),allocatable::temp_log
real(double),dimension(:),allocatable::temp_exp
real(double),dimension(:),allocatable::temp_num_f
real(double)::temp_log_max,force
real(double)::denom,num_f,mean_force_ABF

allocate(temp_log(-nhisto2:nhisto+nhisto2),temp_exp(-nhisto2:nhisto+nhisto2),temp_num_f(-nhisto2:nhisto+nhisto2))
 
!write(*,*),nhisto2,nhisto,nhisto1
temp_log(:)=0
temp_num_f(:)=0!
temp_exp(:)=0


fpabf(:,:) = zero

! Computing the forces as an observable ...
force = - DOT_PRODUCT(fp(:,7),rfilac(:))

if ((icsi >= -nhisto1).and.(icsi <= nhisto+nhisto1)) then 
 cumul_force1(icsi) =  cumul_force1(icsi) + force
mean_force1 (icsi) = cumul_force1(icsi)/histo1(icsi)
end if



!---------BEGIN PRINCIPAL PROGRAM-----------------------------------------------------

!A_ee is given in the main program

!--------eta= eta_ABFee for extended dynamics-----------
forall(iter=-nhisto2:nhisto+nhisto2) temp_log(iter)=-(x_mol(iter)-dcsi)**2/(2.d0*eta_ABFee*temperature)+A_ee(iter)/temperature

temp_log_max=maxval(temp_log)

forall(iter=-nhisto2:nhisto+nhisto2) temp_log(iter)=temp_log(iter)-temp_log_max

forall(iter=-nhisto2:nhisto+nhisto2) temp_exp(iter)=exp(temp_log(iter))! 

forall(iter=-nhisto2:nhisto+nhisto2) temp_num_f(iter)=temp_exp(iter)*(x_mol(iter)-dcsi)

denom=sum(temp_exp)

num_f=sum(temp_num_f)

mean_force2(icsi)=1.d0/eta_ABFee*num_f/(denom)

mean_force_ABF=mean_force_ABFee(icsi)

fpabf(1:3,7)=rfilac(1:3)*mean_force_ABF

fp(1:3,:)=fp(1:3,:)+fpabf(1:3,:)



do iter=-nhisto2, nhisto+nhisto2
P_ee(iter)=temp_exp(iter)/denom
enddo

!----------If histo_equi is true when it_mab > n_equilibre, or histo_equi is false, we fill the histogram of zeta
if(((histo_equi == .true.) .AND. (it_mab > n_equilibre)) .OR. (histo_equi == .False.)) then

do iter= -nhisto1, nhisto+nhisto1 
histo_zeta(iter)=histo_zeta(iter)+P_ee(iter)
enddo

endif


return
end subroutine calfo_ABFee_const_biais



subroutine calfo_ABFee_iter()! this subroutine calculate iterally A bar of ABFee given A_bar_0. ( given mean_force_ABFee)

 USE T_kind_param_m, ONLY:  double
 USE gen_com_m, ONLY: zero,im,imm,low_limit,angst,ev2erg,erg2ev
 USE tab_imm_m
 USE mab_in_ndm_module, ONLY:dcsi,icsi,rfilac,histo,     &
                             mean_force,cumul_force1,nhisto,nhisto1, &
                             nhisto2,mean_force1,histo1,cumul_force1,   &
                             A_dev_ee,A_ee,P_ee,P_ee_num,P_ee_denom,delta_z,&
                             x_mol,eta_ABFee,temperature,omega_abf,exp_A_bar,&
                             mean_force2,it_mab,histo_zeta,n_equilibre,histo_equi,&
                             it_mab,mean_force_ABFee

 implicit none
integer::iter
real(double), dimension(3,imm) :: fpabf,fpabf1
real(double),dimension(:),allocatable::temp_log
real(double),dimension(:),allocatable::temp_exp
real(double),dimension(:),allocatable::temp_num_f
real(double)::temp_log_max,force
real(double)::denom,num_f,mean_force_ABF

allocate(temp_log(-nhisto2:nhisto+nhisto2),temp_exp(-nhisto2:nhisto+nhisto2),temp_num_f(-nhisto2:nhisto+nhisto2))
 
!write(*,*),nhisto2,nhisto,nhisto1
temp_log(:)=0
temp_num_f(:)=0!
temp_exp(:)=0


fpabf(:,:) = zero
fpabf1(:,:)= zero


! Computing the forces as an observable ...
force = - DOT_PRODUCT(fp(:,7),rfilac(:))

if ((icsi >= -nhisto1).and.(icsi <= nhisto+nhisto1)) then 
 cumul_force1(icsi) =  cumul_force1(icsi) + force
mean_force1 (icsi) = cumul_force1(icsi)/histo1(icsi)
end if



!---------BEGIN PRINCIPAL PROGRAM-----------------------------------------------------

!---------1. Compute A_ee(\zeta) for every \zeta, Méthode des trapèzes
A_ee(-nhisto2)=0.d0
do iter=-nhisto2+1,nhisto+nhisto2
A_ee(iter)=A_ee(iter-1)+delta_z*0.5d0*(A_dev_ee(iter-1)+A_dev_ee(iter))
end do

!--------eta= eta_ABFee for extended dynamics-----------
forall(iter=-nhisto2:nhisto+nhisto2) temp_log(iter)=-(x_mol(iter)-dcsi)**2/(2.d0*eta_ABFee*temperature)+A_ee(iter)/temperature

temp_log_max=maxval(temp_log)

forall(iter=-nhisto2:nhisto+nhisto2) temp_log(iter)=temp_log(iter)-temp_log_max

forall(iter=-nhisto2:nhisto+nhisto2) temp_exp(iter)=exp(temp_log(iter))! 

forall(iter=-nhisto2:nhisto+nhisto2) temp_num_f(iter)=temp_exp(iter)*(x_mol(iter)-dcsi)

denom=sum(temp_exp)

num_f=sum(temp_num_f)

mean_force_ABF=1.d0/eta_ABFee*num_f/(denom)

mean_force2(icsi)=mean_force_ABF

fpabf(1:3,7)=rfilac(1:3)*mean_force_ABF

fpabf1(1:3,7)=rfilac(1:3)*mean_force_ABFee(icsi)


fp(1:3,:)=fp(1:3,:)+fpabf(1:3,:)+fpabf1(1:3,:)

!write (41,*) fp(1,7),fpabf(1,7),fpabf1(1,7)
!write (42,*) fp(2,7),fpabf(2,7),fpabf1(2,7)
!write (43,*) fp(3,7),fpabf(3,7),fpabf1(3,7)


do iter=-nhisto2, nhisto+nhisto2
P_ee(iter)=temp_exp(iter)/denom
P_ee_denom(iter)=P_ee_denom(iter)+P_ee(iter)
P_ee_num(iter)=P_ee_num(iter)+(x_mol(iter)-dcsi)*P_ee(iter)
enddo

!----------If histo_equi is true when it_mab > n_equilibre, or histo_equi is false, we fill the histogram of zeta
if(((histo_equi == .true.) .AND. (it_mab > n_equilibre)) .OR. (histo_equi == .False.)) then

do iter= -nhisto1, nhisto+nhisto1 
histo_zeta(iter)=histo_zeta(iter)+P_ee(iter)
enddo

endif

!---------Calcule A_dev if necessaire
do iter=-nhisto1,nhisto+nhisto1
A_dev_ee(iter)=1.d0/eta_ABFee*(P_ee_num(iter))/(P_ee_denom(iter)+1.d0/(omega_abf*dble(nhisto)))
enddo

if(nhisto1 < nhisto2) then
!---------------en dehors de nhisto1 
do iter=-nhisto2,-nhisto1
A_dev_ee(iter)=0.d0
enddo

do iter=nhisto1+nhisto,nhisto+nhisto2
A_dev_ee(iter)=0.d0
enddo
endif

return

end subroutine calfo_ABFee_iter
