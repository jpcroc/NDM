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
                             mean_force2,it_mab,histo_zeta,n_equilibre,histo_equi,&
                             it_mab,abf_mode,potist,ene_einstein,ene0,fpeinstein, &
                             ha_mix,equit,atom_to_jump,itype_reaction

 implicit none
integer::iter,ia,jx
integer, save :: itcount=0
real(double), dimension(3,imm) :: fpabf
real(double),dimension(:),allocatable::temp_log
real(double),dimension(:),allocatable::temp_exp,U_Aee
real(double),dimension(:),allocatable::temp_num_f
real(double)::temp_log_max,force,tmp_num
real(double)::denom,num_f,mean_force_ABF

allocate(temp_log(-nhisto2:nhisto+nhisto2),  &
         temp_exp(-nhisto2:nhisto+nhisto2),  &
         temp_num_f(-nhisto2:nhisto+nhisto2),&
         U_Aee(-nhisto2:nhisto+nhisto2)   )
 
!write(*,*),nhisto2,nhisto,nhisto1
temp_log(:)=0
temp_num_f(:)=0!
temp_exp(:)=0


fpabf(:,:) = zero

if (abf_mode==1) then
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
! mcmCHECK if that is OK ? 
A_ee(-nhisto2)=0.d0
do iter=-nhisto2+1,nhisto+nhisto2
A_ee(iter)=A_ee(iter-1)+delta_z*0.5d0*(A_dev_ee(iter-1)+A_dev_ee(iter))
end do
if (it_mab<10) A_ee=0.d0

if (abf_mode==1) then
 !--2. Compute pi_A_ee (csi,q)
 ! xmol is zeta, eta = eta_ABFee for exteded dynamics 
 ! pi_Aee = num / denom = \exp{ -\frac{beta}{2*eta} (xmol(iter) - csi)**2  /denom  }
 ! I think that the sign is in the opposite direction ??? mcmCHECK
 forall(iter=-nhisto2:nhisto+nhisto2) temp_log(iter)=-(x_mol(iter)-dcsi)**2/(2.d0*eta_ABFee*temperature)+A_ee(iter)/temperature
 ! What the hack is that !!!
 temp_log_max=MAXVAL(temp_log(:) )
 temp_log(:)=temp_log(:)-temp_log_max

 temp_exp(-nhisto2:nhisto+nhisto2)=exp(temp_log(nhisto2:nhisto+nhisto2))! 
 ! strange way to have a integral :
 denom=SUM(temp_exp(:))


 !--3. Compute the E(grad U(.,q)|q) 
 temp_num_f(-nhisto2:nhisto+nhisto2)=temp_exp(nhisto2:nhisto+nhisto2)*(x_mol(-nhisto2:nhisto+nhisto2)-dcsi)
 ! strange way to have the integral but it is corect because we have only the fraction
 num_f=SUM(temp_num_f(:))
 ! 
 mean_force_ABF=1.d0/eta_ABFee*num_f/denom
 !
 mean_force2(icsi)=mean_force_ABF

 fpabf(1:3,7)=rfilac(1:3)*mean_force_ABF

 fp(1:3,:)=fp(1:3,:)+fpabf(1:3,:)
end if !abf_mode==1

if (abf_mode==2) then
 ! xmol is zeta
 ! U(zeta, q) = zeta*potist + (1-zeta)*(ene_einstein + ene0) 
 !--2. compute the  pi_A_ee(\zeta) = \exp{U(zeta,q) / int_\zeta_min^\zeta_max{\exp{U(zeta,q) d\zeta}
 !--2.a num \exp{U(zeta,q)
  do iter=-nhisto1,nhisto+nhisto1
     U_Aee(iter) = (1.d0-x_mol(iter))*(ene_einstein) + x_mol(iter)*(potist-ene0)
     temp_log(iter)=-(U_Aee(iter)-A_ee(iter))/temperature
  end do
  temp_log_max=MAXVAL(temp_log(:) )
  temp_log(:)=temp_log(:)-temp_log_max

     
  do iter=-nhisto1,nhisto+nhisto1
     temp_exp(iter)=exp(temp_log(iter))!
      if ((temp_exp(iter)+1.0).eq.temp_exp(iter)) then
        write(*,*) 'WARNING:  NaN detected look in fort.333 file'
        write(333,'(2i5,7D21.8)') it_mab, iter, U_Aee(iter),A_ee(iter), temp_log(iter),temp_exp(iter),potist,ene0,ene_einstein
      end if  
  end do

 !--2.b denom: int_\zeta_min^\zeta_max{\exp{U(zeta,q) d\zeta}
  denom=0.d0
  do iter=-nhisto1+1,nhisto+nhisto1
   !C_n in my notes
   denom=denom + 0.5d0*(temp_exp(iter-1)+temp_exp(iter) )*delta_z
  end do

 !--3. Compute the E(\grad_q U(.,q)|q)
 !In this case E[\grad_q U(.,q)|q]=E[ (1-\zeta) \grad_q U^HA(q) + \zeta \grad_q U(q) ]
 !We should just remind that \grad_q U^HA(q) = -fpeinstein (q) ; \grad_q U(q) = - fp(q)
 do ia=1,im
  do jx=1,3
   do iter=-nhisto1,nhisto+nhisto1
    temp_num_f(iter)=-((1.d0-x_mol(iter))*fpeinstein(jx,ia)+x_mol(iter)*fp(jx,ia))*temp_exp(iter)
   end do
   tmp_num=0.d0
   do iter=-nhisto1+1,nhisto+nhisto1
    tmp_num=tmp_num  + 0.5d0*(temp_num_f(iter-1)+temp_num_f(iter) )*delta_z
   end do
   if ((jx==1).and.(ia==7)) then
   
   !write(*,*) -tmp_num/denom, fp(jx,ia)
   end if 
  if ((histo_equi == .true.) .And. (it_mab >= n_equilibre)) then
    fp(jx,ia)=-tmp_num/denom
   else 
    fp(jx,ia)=fp(jx,ia)+fpeinstein(jx,ia)
  end if 
 end do
end do

end if !abf_mode==2 

if (abf_mode==22) then
 ! xmol is zeta
 ! U(zeta, q) = zeta*(potist-ene0 - ha_mix*U_HA)
 !--3. compute the  pi_A_ee(\zeta | q ) = \exp{U(zeta,q) / int_\zeta_min^\zeta_max{\exp{U(zeta,q) d\zeta}
 !--3.a num \exp{U(zeta,q)


  do iter=-nhisto1,nhisto+nhisto1
     U_Aee(iter) = x_mol(iter)*(potist-ene0+ha_mix*ene_einstein-equit)
     temp_log(iter)=-(U_Aee(iter)-A_ee(iter))/temperature
     temp_exp(iter)=exp(temp_log(iter))!
      if ((temp_exp(iter)+1.0).eq.temp_exp(iter)) then
       itcount=itcount+1
        write(*,*) 'WARNING:  NaN detected look in fort.333 file'
        write(333,'(2i5,9D21.8)') it_mab, iter, U_Aee(iter),A_ee(iter), temp_log(iter),temp_exp(iter),potist,&
               ene0,ene_einstein,equit, (potist-ene0+ha_mix*ene_einstein-equit)*erg2ev
        if (itcount==5) stop
      end if  
  end do

 !--3.b denom: int_\zeta_min^\zeta_max{\exp{U(zeta,q) d\zeta}
  denom=0.d0
  do iter=-nhisto1+1,nhisto+nhisto1
   denom=denom + 0.5d0*(temp_exp(iter-1)+temp_exp(iter) )*delta_z
  end do
 
 ! if (it_mab > n_equilibre+100) write(*,'("calfo_ABFee",3g15.4)') ene_einstein, fp(1,7),fpeinstein(1,7) 
 !--4. Compute the E(\grad_q U(.,q)|q)
 !In this case E[\grad_q U(.,q)|q]=E[ \zeta \grad_q U(q) ]
 !We should just remind that  \grad_q U(q) = - fp(q)
 do ia=1,im
  do jx=1,3
   do iter=-nhisto1,nhisto+nhisto1
    temp_num_f(iter)=-x_mol(iter)*(fp(jx,ia)+ha_mix*fpeinstein(jx,ia))*temp_exp(iter)
   end do
   tmp_num=0.d0
   do iter=-nhisto1+1,nhisto+nhisto1
    tmp_num=tmp_num  + 0.5d0*(temp_num_f(iter-1)+temp_num_f(iter) )*delta_z
   end do
!debug<
   if ((jx==1).and.(ia==7)) then
   !write(*,*) -tmp_num/denom, fp(jx,ia)
   end if
!>debug 
  if ((histo_equi == .true.) .And. (it_mab >= n_equilibre)) then
    fp(jx,ia)=-tmp_num/denom
  else 
    fp(jx,ia)=fp(jx,ia)+fpeinstein(jx,ia)
  end if 
 end do
end do

end if !abf_mode==22 


!write (41,*) fp(1,7),fpabf(1,7)
!write (42,*) fp(2,7),fpabf(2,7)
!write (43,*) fp(3,7),fpabf(3,7)
! 1. The computationb of the mean force A_dev_ee (\zeta)
  
!1a We prepa denom and num  in order to compute A_dev_ee (\zeta}
! num= \sum_{all_md_steps) \grad_\zeta U(\zeta,q)*pi_A_ee(\zeta | q )
! denom = \tau + \sum_{all_md_steps} pi_A_ee(\zeta | q ) (below the pi_A_ee si denoted by is P_ee) 
if (abf_mode==1) then
 do iter=-nhisto2, nhisto+nhisto2
  P_ee(iter)=temp_exp(iter)/denom
  P_ee_denom(iter)=P_ee_denom(iter)+P_ee(iter)
  P_ee_num(iter)=P_ee_num(iter)+(x_mol(iter)-dcsi)*P_ee(iter)
 enddo
  P_ee_num(:)=1.d0/eta_ABFee*P_ee_num(:)
end if 


if (abf_mode==2) then
 do iter=-nhisto1, nhisto+nhisto1
  P_ee(iter)=temp_exp(iter)/denom
  P_ee_denom(iter)=P_ee_denom(iter)+P_ee(iter)
  P_ee_num(iter)=P_ee_num(iter)+(potist-ene_einstein-ene0)*P_ee(iter)
 enddo
end if 


if (abf_mode==22) then
 do iter=-nhisto1, nhisto+nhisto1
  P_ee(iter)=temp_exp(iter)/denom
  P_ee_denom(iter)=P_ee_denom(iter)+P_ee(iter)
  P_ee_num(iter)=P_ee_num(iter)+(potist-ene0+ha_mix*ene_einstein-equit)*P_ee(iter)
 enddo
end if 


!----------If histo_equi is true when it_mab > n_equilibre, or histo_equi is false, we fill the histogram of zeta
if(((histo_equi == .true.) .AND. (it_mab > n_equilibre)) .OR. (histo_equi == .False.)) then

    do iter= -nhisto1, nhisto+nhisto1 
      histo_zeta(iter)=histo_zeta(iter)+P_ee(iter)
    enddo

endif

 if(nhisto1 < nhisto2) then
  !---------------en dehors de nhisto1 
   do iter=-nhisto2,-nhisto1
    A_dev_ee(iter)=0.d0
   enddo

   do iter=nhisto1+nhisto,nhisto+nhisto2
    A_dev_ee(iter)=0.d0
   enddo
 endif


!1b Final step to have  A_dev_ee (\zeta} = num/denom)
! num= \sum_{all_md_steps) \grad_\zeta U(\zeta,q)*pi_A_ee(\zeta | q )
! denom = \tau + \sum_{all_md_steps} pi_A_ee(\zeta | q ) (below the pi_A_ee si denoted by is P_ee) 

 ! 
 ! from where comes this nhisto !!!! mcmCHECK
 if (abf_mode==1) then
 do iter=-nhisto1,nhisto+nhisto1
     A_dev_ee(iter)=P_ee_num(iter)/(P_ee_denom(iter)+1.d0/(omega_abf*dble(nhisto)))
 enddo
end if 

 if (abf_mode==2) then
 do iter=-nhisto1,nhisto+nhisto1
     A_dev_ee(iter)=P_ee_num(iter)/(P_ee_denom(iter)+1.d0/omega_abf)
 enddo
end if 

 if (abf_mode==22) then
 do iter=-nhisto1,nhisto+nhisto1
     A_dev_ee(iter)=P_ee_num(iter)/(P_ee_denom(iter)+1.d0/omega_abf)
 enddo
end if 

return

end subroutine calfo_ABFee


