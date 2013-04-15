subroutine mab
!-----------------------------------------------

      USE T_kind_param_m, ONLY:  double
      use gen_com_m
      use var_pot
      use tab_imm_m
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      !use DEFS
      use mab_in_ndm_module
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      implicit none
!      integer  :: ielat(imm)
!      integer  :: iwmax(imm)
!      integer  :: ityp(imm)
!      real(double)  :: xp(3,imm)
!      real(double)  :: xpp(3,imm)
!      real(double)  :: vp(3,imm)
!      real(double)  :: ax(3,imm)



! 
! Copyright LL Cao and all NDM band, April- 2013

  write(6,*)
  write(6,*)
  write(6,*)'************ DEBUT DE MAB ****************'
  write(6,*)
  write(6,*)
  call allocate_mab()
  !call force_constant(xp, xpp, vp, ax, fp, ielat, iwmax, ityp)


  write(6,*)
  write(6,*)
  write(6,*)'************  FIN  DE MAB ****************'
  write(6,*)
  write(6,*)
 

  return

  
  end subroutine mab
