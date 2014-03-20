! ************************************************
!           Sous-programme test_minimum_abf
!           MCM - 02/2014
! ************************************************

subroutine test_minimum_abf
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use tab_imm_m
  use var_pot
  use mab_in_ndm_module, ONLY : ene0
  implicit none
  real(double) :: forctot, formax
  !-----------------------------------------------
  real(double), dimension(3,imm) :: fp0
  


  imd = im
  nad(:ntyp) = na(:ntyp)
  vp=0.0
  fp=0.
  call caltabt
  call caltabi
  call calfo
  ene0=potist
  write (6, '(" MAB: ene0 or  *Epot (eV) = ",D21.12)') ene0*erg2ev
  fp0=fp


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
   if (rangph==0)        write(6,*) 'MAB: WARNING - The forces are not converged compared to &
        din file'
    if (rangph==0)       write(6,'(" MAB: WARNING - formax is greater than fpstop...:",2D12.5)') & 
        formax,fpstop
     if (rangph==0)      write(6,*) 'MAB: WARNING - Check the input'        
       end if
   end if       

   if (fsumstop>0) then   
      if (forctot.gt.fsumstop) then
    if (rangph==0)       write(6,*) 'MAB: WARNING - The forces are not converged compared to &
        din file'
    if (rangph==0)       write(6,*) 'MAB: WARNING - forctot greater than fsumstop', &
       forctot,fsumstop
   if (rangph==0)       write(6,*) 'MAB: WARNING - Check the input'        
      end if
   end if   
 end if
end if  
  
  end subroutine test_minimum_abf



