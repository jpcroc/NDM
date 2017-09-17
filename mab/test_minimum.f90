! ************************************************
!           Sous-programme test_minimum
!           MCM - 08/2012
! ************************************************

subroutine test_minimum !(xp, xpp, vp, ax, fp, ielat, iwmax, iwmax2, ityp)
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use tab_imm_m
  use var_pot
  use mab_in_ndm_module
  implicit none
  !-----------------------------------------------
  !   G l o b a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  ! real(double) :: fps(3,imm)
  real(double) :: forctot, formax
  !-----------------------------------------------
  !
  

  ! MPI
  if (rangph==0) write (6, *) 'MAB: Testing if the forces are zero.'

  imd = im
  nad(:ntyp) = na(:ntyp)

  ! appel de la routine generale des forces
  !position de départ
  vp=0.0
  !      do i=1,im
  !         write(6,*)i,xp(1,i),xp(2,i),xp(3,i)
  !      end do
  fp=0.
  call calfo_atomic_forces(1)
  epot0=potist*erg2eV
  if (rangph==0)   write (6, '(" MAB: *Epot (eV) = ",D21.12)') epot0
  ene0 = epot0
  fp0=fp
!debugCOS  do i=1,im
!debugCOS     write(6,'(I7,3D21.12)')i,fp0(1,i),fp0(2,i),fp0(3,i)
!debugCOS  end do
  !      fps=fp


  IF (lFrozen) THEN
     forctot = sqrt( SUM( fp0(:,1:im)**2, .NOT.Frozen(:,1:im) ) )
     formax = MaxVal( Abs(fp0(:,1:im)), .NOT.Frozen(:,1:im) ) 
  ELSE
     forctot=sqrt( SUM(fp0(1:3,1:im)**2) )
     formax = MaxVal( Abs(fp0(:,1:im)) )
  END IF
 
if (rangph==0) then
 if (lEev) then 
   forctot = forctot*erg2eV/angst
   formax  = formax*erg2eV/angst
   if (fpstop>0) then   
       if (formax.gt.fpstop) then
   if (rangph==0)        write(6,*) 'PHONDY: WARNING - The forces are not converged compared to &
        din file'
    if (rangph==0)       write(6,'(" PHONDY: WARNING - formax is greater than fpstop...:",2D12.5)') & 
        formax,fpstop
     if (rangph==0)      write(6,*) 'PHONDY: WARNING - Check the input'        
       end if
   end if       

   if (fsumstop>0) then   
      if (forctot.gt.fsumstop) then
    if (rangph==0)       write(6,*) 'PHONDY: WARNING - The forces are not converged compared to &
        din file'
    if (rangph==0)       write(6,*) 'PHONDY: WARNING - forctot greater than fsumstop', &
       forctot,fsumstop
   if (rangph==0)       write(6,*) 'PHONDY: WARNING - Check the input'        
      end if
   end if   
 end if
end if  
  
  end subroutine test_minimum



