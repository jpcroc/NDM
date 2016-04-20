! This file contains:
!  - calfo_ABFee
!  

subroutine calfo_ABFee()

 USE T_kind_param_m, ONLY:  double
 USE gen_com_m, ONLY: zero,im,imm,low_limit,angst,ev2erg,erg2ev
 USE tab_imm_m
 USE mab_in_ndm_module, ONLY:dcsi,icsi,rfilac,histo,     &
                             mean_force,cumul_force1,nhisto,nhisto1, &
                             nhisto2,mean_force1,histo1,cumul_force1,   &
                             A_dev_ee,A_ee,P_ee,P_ee_num,P_ee_denom,delta_z,&
                             x_mol,eta_ABFee,temperature,omega_abf,exp_A_bar,&
                             mean_force2,it_mab,histo_zeta,& 
                             neq_lang, n_equilibre,histo_equi,&
                             it_mab,abf_mode,potist,ene_einstein,ene0,fpeinstein, &
                             ha_mix,equit,atom_to_jump,itype_reaction, &
                             abf_mode_reaction, abf_mode_alchemical, abf_mode_temperature

 implicit none
integer::iter,ia,jx
integer, save :: itcount=0
real(double), dimension(3,imm) :: fpabf
real(double),dimension(:),allocatable::temp_log
real(double),dimension(:),allocatable::temp_exp,U_Aee
real(double),dimension(:),allocatable::temp_num_f
real(double)::temp_log_max,force,tmp_num
real(double)::denom,num_f,mean_force_ABF
real(double) :: omega_abf_i

allocate(temp_log(-nhisto2:nhisto+nhisto2),  &
         temp_exp(-nhisto2:nhisto+nhisto2),  &
         temp_num_f(-nhisto2:nhisto+nhisto2),&
         U_Aee(-nhisto2:nhisto+nhisto2)   )
 
!write(*,*),nhisto2,nhisto,nhisto1
temp_log(:)=0
temp_num_f(:)=0!
temp_exp(:)=0


fpabf(:,:) = zero

if (abf_mode==abf_mode_reaction) then
! Computing the forces as an observable ...
if (itype_reaction==0) force = - DOT_PRODUCT(fp(:,atom_to_jump),rfilac(:))

 if ((icsi >= -nhisto1).and.(icsi <= nhisto+nhisto1)) then 
  cumul_force1(icsi) =  cumul_force1(icsi) + force
  mean_force1 (icsi) =  cumul_force1(icsi)/histo1(icsi)
 end if

end if 

!---------BEGIN PRINCIPAL PROGRAM-----------------------------------------------------

!---------1. Compute A_ee(\zeta) for every \zeta, from A_dev_ee(\zeta).  Methode des trapezes
! If delta_z is OK can be used  for any ABFee type 1 or 2
A_ee(-nhisto2)=0.d0
do iter=-nhisto2+1,nhisto+nhisto2
A_ee(iter)=A_ee(iter-1)+delta_z*0.5d0*(A_dev_ee(iter-1)+A_dev_ee(iter))
end do
if (it_mab<4) A_ee=0.d0

select case (abf_mode)

case  (abf_mode_reaction) 
 !--2. Compute pi_A_ee (csi,q)
 ! xmol is zeta, eta = eta_ABFee for extended dynamics 
 ! pi_Aee = num / denom = \exp{ -\frac{beta}{2*eta} (x_mol(iter) - csi)**2  /denom  }
 ! I think that the sign is in the opposite direction ??? mcmCHECK
 forall(iter=-nhisto2:nhisto+nhisto2) temp_log(iter)=-(x_mol(iter)-dcsi)**2/(2.d0*eta_ABFee*temperature)+A_ee(iter)/temperature
 ! Trick to avoid NaN errors ...
 temp_log_max=MAXVAL(temp_log(:) )
 temp_log(:)=temp_log(:)-temp_log_max

  do iter=-nhisto2,nhisto+nhisto2
     temp_exp(iter)=exp(temp_log(iter))!
      if (temp_exp(iter) /= temp_exp(iter)) then
        write(*,*) 'WARNING:  NaN detected look in fort.333 file'
        write(333,'(2i5,6D21.8)') it_mab, iter, x_mol(iter)-dcsi,A_ee(iter), temp_log(iter),temp_exp(iter),potist,ene0,ene_einstein
      end if  
  end do


denom=SUM(temp_exp(:))


 !--3. Compute the E(grad U(.,q)|q) 
 forall(iter=-nhisto2:nhisto+nhisto2)   temp_num_f(iter)=temp_exp(iter)*(x_mol(iter)-dcsi)/eta_ABFee
 ! strange way to have the integral but it is corect because we have only the fraction
 num_f=SUM(temp_num_f(:))
 mean_force_ABF=num_f/denom
 mean_force2(icsi)=mean_force_ABF
! write(*,*) atom_to_jump
 fpabf(1:3,atom_to_jump)=rfilac(1:3)*mean_force_ABF
 fp(1:3,:)=fp(1:3,:)+fpabf(1:3,:)

!----------------end case abf_mode==abf_mode_reaction=1 


case (abf_mode_alchemical)
 ! xmol is zeta
 ! U_Aee is U(zeta, q) = (1-zeta)*(ene_einstein + ene0) + zeta * potist
 ! Aee   is A_n(zeta)
 ! temp_log is - beta * U_A_n or -beta (U_Aee - Aee)
 ! 
 !2. compute the  conditional probability 
 !    p_A_n (zeta | q_n )  = \exp{U_A_n(zeta,q_n) / int_\zeta_min^\zeta_max{\exp{U_a_n(zeta,q_n) d\zeta}
 !--2.a num is \exp{U_A_n (zeta,q_n)}
  do iter=-nhisto1,nhisto+nhisto1
    !good_old  U_Aee(iter) = (1.d0-x_mol(iter))*(ene_einstein) + x_mol(iter)*(potist-ene0)
     U_Aee(iter) = (1.d0-x_mol(iter))*(ene_einstein+ene0) + x_mol(iter)*(potist)
     temp_log(iter)=-(U_Aee(iter)-A_ee(iter))/temperature
  end do

 ! add an extra trick to avoid the NaN errors. 
 ! No influence in the p_A_n (zeta | q_n) pi_Aee(\zeta | q) temp_log_max appears 
 ! in the same way in the numerator and denumerator ... 
  temp_log_max=MAXVAL(temp_log(:) )
  temp_log(:)=temp_log(:)-temp_log_max

  do iter=-nhisto1,nhisto+nhisto1
     temp_exp(iter)=exp(temp_log(iter))
      if ( temp_exp(iter) /= temp_exp(iter)) then
        write(*,*) 'WARNING:  NaN detected look in fort.333 file'
        write(333,'(2i5,7D21.8)') it_mab, iter, U_Aee(iter),A_ee(iter), temp_log(iter),temp_exp(iter),potist,ene0,ene_einstein
        stop
      end if  
  end do

 !--2.b denom: int_\zeta_min^\zeta_max{\exp{U_A_n(zeta,q_n) d\zeta}
  denom=0.d0
  do iter=-nhisto1+1,nhisto+nhisto1
   !C_n in my notes
   denom=denom + 0.5d0*(temp_exp(iter-1)+temp_exp(iter) )*delta_z
  end do


 !--3. Compute the F_A_n(q_n)= E(\grad_q U(.,q)|q)= 
 !                             int_\zeta_min^\zeta_max  \grad_q U(zeta,q_n) p_A_n(zeta,q_n) d zeta
 !    In this case \grad_q U(zeta,q) = (1-zeta) \grad_q U^HA(q) + zeta \grad_q ( U(q) - U0) 
 !    just remind that \grad_q U^HA(q) = -fpeinstein (q) ; \grad_q U(q) = - fp(q)
 !    tmp_num
 do ia=1,im
  do jx=1,3
   do iter=-nhisto1,nhisto+nhisto1
    temp_num_f(iter)=-((1.d0-x_mol(iter))*fpeinstein(jx,ia)+x_mol(iter)*fp(jx,ia))*temp_exp(iter)
   end do
   tmp_num=0.d0
   do iter=-nhisto1+1,nhisto+nhisto1
    tmp_num=tmp_num  + 0.5d0*(temp_num_f(iter-1)+temp_num_f(iter) )*delta_z
   end do
   
  if ((histo_equi .eqv. .true.) .And. (it_mab - neq_lang >= n_equilibre)) then
    fp(jx,ia)=-tmp_num/denom
   else 
    fp(jx,ia)=fp(jx,ia)+fpeinstein(jx,ia)
  end if 
 end do
end do
!----------------------------end case abf_mode==abf_mode_alchemical=2 

case(abf_mode_temperature)
 ! xmol is zeta
 ! U(zeta, q) = zeta*(potist-ene0+ha_mix*U_HA-equit)
 !--2. compute the  pi_A_ee(\zeta | q ) = \exp{U(zeta,q) / int_\zeta_min^\zeta_max{\exp{U(zeta,q) d\zeta}
 !--2.a num \exp{U(zeta,q)


  do iter=-nhisto1,nhisto+nhisto1
     U_Aee(iter) = x_mol(iter)*(potist-ene0+ha_mix*ene_einstein-equit)
     temp_log(iter)=-(U_Aee(iter)-A_ee(iter))/temperature
  end do
! add an extra trick to avoid the NaN errors. No influence in the pi_Aee(\zeta | q) appears in the same way
! in the numerator and denumerator ... 
  temp_log_max=MAXVAL(temp_log(:))
  temp_log(:)=temp_log(:)-temp_log_max

  do iter=-nhisto1,nhisto+nhisto1
     temp_exp(iter)=exp(temp_log(iter))!
      if ((temp_exp(iter)+1.0).eq.temp_exp(iter)) then
       itcount=itcount+1
        write(*,*) 'WARNING:  NaN detected look in fort.333 file'
        write(333,'(2i5,9D21.8)') it_mab, iter, U_Aee(iter),A_ee(iter), temp_log(iter),temp_exp(iter),potist,&
               ene0,ene_einstein,equit, (potist-ene0+ha_mix*ene_einstein-equit)*erg2ev
!debug        if (itcount==5) stop
      end if  
  end do

 !--2.b denom: int_\zeta_min^\zeta_max{\exp{U(zeta,q) d\zeta}
  denom=0.d0
  do iter=-nhisto1+1,nhisto+nhisto1
   ! Cn in my notes 
   denom=denom + 0.5d0*(temp_exp(iter-1)+temp_exp(iter) )*delta_z
  end do
 
 !--4. Compute the E(\grad_q U(.,q)|q)
 !In this case E[\grad_q U(.,q) | q ]= E[ \zeta \grad_q ( U(q) + ha_mix * U_HA(q) - E0 + equit) ]
 !We should just remind that  \grad_q U(q) = - fp(q) and  \grad_q U_HA(q)=- fpeinstein (q)
 do ia=1,im
  do jx=1,3
   do iter=-nhisto1,nhisto+nhisto1
    temp_num_f(iter)=-x_mol(iter)*(fp(jx,ia)+ha_mix*fpeinstein(jx,ia))*temp_exp(iter)
   end do
   tmp_num=0.d0
   do iter=-nhisto1+1,nhisto+nhisto1
    tmp_num=tmp_num  + 0.5d0*(temp_num_f(iter-1)+temp_num_f(iter) )*delta_z
   end do
!debug    if ((jx==1).and.(ia==7)) then
!debug    write(*,*) -tmp_num/denom, fp(jx,ia)
!debug    end if
  if ((histo_equi .eqv. .true.) .And. (it_mab- neq_lang >= n_equilibre)) then
    fp(jx,ia)=-tmp_num/denom
  else 
    fp(jx,ia)=fp(jx,ia)+fpeinstein(jx,ia)
  end if 
 end do
end do

!------------------end case abf_mode==abf_mode_temperature=22 
end select 

! 1. The computation of the mean force A_dev_ee (\zeta)
  
!1a We prepa denom and num  in order to compute A_dev_ee (\zeta}
! num= \sum_{all_md_steps) \grad_\zeta U(\zeta,q)*pi_A_ee(\zeta | q )
! denom = \tau + \sum_{all_md_steps} pi_A_ee(\zeta | q ) (below the pi_A_ee si denoted by is P_ee) 
 do iter=-nhisto1,nhisto+nhisto1
 ! How to handle zeros or/and NaN
  !zero 
  if (denom==0.d0) then
     P_ee(iter)= 0.d0
   else 
     P_ee(iter)=temp_exp(iter)/denom
  end if 
  !NaN
  if (P_ee(iter)/=P_ee(iter)) P_ee(iter)=0.d0
 end do
 !
 !
 select case (abf_mode)
  !
  case (abf_mode_reaction)
    !
    do iter=-nhisto1, nhisto+nhisto1
     P_ee_denom(iter)=P_ee_denom(iter)+P_ee(iter)
     P_ee_num(iter)=P_ee_num(iter)+(x_mol(iter)-dcsi)*P_ee(iter)/eta_ABFee
    enddo

  case (abf_mode_alchemical)
    !
    do iter=-nhisto1, nhisto+nhisto1
      P_ee_denom(iter)=P_ee_denom(iter)+P_ee(iter)
      P_ee_num(iter)=P_ee_num(iter)+(potist-ene_einstein-ene0)*P_ee(iter)
    enddo

  case (abf_mode_temperature)
    !
    do iter=-nhisto1, nhisto+nhisto1
     P_ee_denom(iter)=P_ee_denom(iter)+P_ee(iter)
     P_ee_num(iter)=P_ee_num(iter)+(potist-ene0+ha_mix*ene_einstein-equit)*P_ee(iter)
   enddo
   !
 end select 


!----------If histo_equi is true when it_mab > n_equilibre, or histo_equi is false, we fill the histogram of zeta
if(((histo_equi .eqv. .true.) .AND. (it_mab > n_equilibre)) .OR. (histo_equi .eqv. .False.)) then

    do iter= -nhisto1, nhisto+nhisto1 
      histo_zeta(iter)=histo_zeta(iter)+P_ee(iter)
    enddo
endif

 if(nhisto1 < nhisto2) then
    A_dev_ee(-nhisto2:-nhisto1)=0.d0
    A_dev_ee(nhisto1+nhisto:nhisto+nhisto2)=0.d0
 endif


!1b Final step to have  A_dev_ee (\zeta} = num/denom)
! num= \sum_{all_md_steps) \grad_\zeta U(\zeta,q)*pi_A_ee(\zeta | q )
! denom = \tau + \sum_{all_md_steps} pi_A_ee(\zeta | q ) (below the pi_A_ee si denoted by is P_ee) 

 ! 
 ! from where comes this nhisto !!!! mcmCHECK
 omega_abf_i = 1.d0/omega_abf
 select case(abf_mode)
   case(abf_mode_reaction) 
      forall(iter=-nhisto1:nhisto+nhisto1) A_dev_ee(iter)=P_ee_num(iter)/(P_ee_denom(iter)+omega_abf_i/dble(nhisto))
   case(abf_mode_alchemical) 
      forall(iter=-nhisto1:nhisto+nhisto1) A_dev_ee(iter)=P_ee_num(iter)/(P_ee_denom(iter)+omega_abf_i)
   case(abf_mode_temperature) 
      forall(iter=-nhisto1:nhisto+nhisto1) A_dev_ee(iter)=P_ee_num(iter)/(P_ee_denom(iter)+omega_abf_i)
  end select

end subroutine calfo_ABFee


