! This part should be adapted for the problem investigated.
! There is no universal solution. 
! You must customize your problem in order to have your 
!             -reaction coordinate
!             -extra constraints in your system


subroutine reaction ()
  
 USE T_kind_param_m, ONLY:  double
 USE gen_com_m, ONLY: lenfnam,fnam,im,imm,ev2erg
 USE tab_imm_m
 USE mab_in_ndm_module, ONLY: m_i,m_tot,normxlac, rfilac,xlaci,xbarini,xbar,& 
                              xi_min,dcsi,icsi,delta_z,abf_mode,atom_to_jump,itype_reaction
 implicit none

!  local variables ...
   integer :: ic

 do ic=1,3
    xbar(ic)  = sum(xp(ic,1:im)*m_i(ic,1:im))/m_tot
 enddo


if (itype_reaction==0) then
 ! In this case the vacancy in the atom no 7 ...
  if (abf_mode==1) then 
    dcsi=DOT_PRODUCT(rfilac(:),xp(:,atom_to_jump)-xlaci(:)-xbar(:)+xbarini(:))
    !write (*,*) dcsi ,delta_z
    !stop
  end if
end if 

if (itype_reaction==1) then
  if (abf_mode==1) then

  end if 
end if 
    icsi=nint((dcsi-xi_min)/delta_z)
    
    !write (*,*) 'icsi.(reaction) ...', dcsi,delta_z,icsi 

 return
  end subroutine reaction


subroutine calfoblock()
 USE T_kind_param_m, ONLY:  double
 USE gen_com_m, ONLY: zero,im,imm,low_limit,angst,ev2erg,erg2ev
 USE tab_imm_m
 USE mab_in_ndm_module, ONLY:radiussph,deltasph,xbar,xbarini,xp0,maxforce,nsite_block,isite_block
 
 real(double), dimension(3,imm) :: fpblock
 integer :: ic
 real(double) :: rp1,dFerDir

 fpblock(:,:) = zero
! Computing the forces on the protectives spheres...
 
 do ii=1,nsite_block
   ic = isite_block(ii)
!debug  if (ic/=7) then
   
   rp1=dsqrt(SUM((xp(:,ic)-xp0(:,ic)-xbar(:)+xbarini(:))**2))
    if (dabs(rp1).gt.low_limit) then
     fpblock(:,ic)=(xp(:,ic)-xp0(:,ic)-xbar(:)+xbarini(:))*maxforce*ev2erg*dFerDir(rp1,radiussph,deltasph)/rp1
     else
     fpblock(:,ic)=zero
    end if
    if (rp1 > (radiussph+deltasph)) then
     write(6,*) 'WARNING the atom ic ESCAPED from the protective domains',ic,rp1*angst 
    ! write(6,*) 'The force is ...', dsqrt(SUM(xp(:,ic)**2))*angst, dsqrt(SUM(fpblock(:,ic)**2))*erg2ev/angst
    end if

!debug  end if
 end do

! Updating the forces ...
 do ii=1,nsite_block
  ic=isite_block(ii)
!debug if (ic/=7) then 
  fp(1:3,ic)=fp(1:3,ic)+fpblock(1:3,ic)
!debug   end if
 end do

end subroutine calfoblock

