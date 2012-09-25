subroutine sundae
!-----------------------------------------------
      USE T_kind_param_m, ONLY:  double
      use gen_com_m
      use tab_imm_m
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
!      use DEFS
!      use random_art
!      use lanczos_defs

#ifdef V1
       use teledyn_module
#else
       use sundae_module
#endif

!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      implicit none
! This is the main program for TELEDYN nouveau version 2010-11
! 
! Copyleft M. Athenes & M.-C. Marinica 
! last french  touch 03.03.2008, 15h20m 
! last italian touch: 
!
!  integer :: i, ierror
  integer ic
!  parameter (nfenetre=100)
  real(double) :: xalea
!  character(len=20) :: fname
!  character(4) :: scounter
!  character(len=150) :: format_output
!  real(8), dimension(:), allocatable, target  :: diff
  integer :: igerme
  character (len=80) :: fnamtin
real (double), dimension(3,im)::xpepe

#ifdef V1
  integer :: lutin

  namelist /input_teledyn/ lta,lljc,text_teledyn,niteration,nchemin,kappa,kappaE,alphadd,massadd,tstep,bargamma,igerme,nfrequence,xtempmin

  igerme = 0
  xtempmin = 1.0e-9
  lutin=54
  fnamtin = fnam(1:lenfnam)//'.tin'
  ! variables de dynamique
  write(*,*) 'file name', fnamtin

  open(unit=lutin, file=fnamtin, status='unknown', err=567)
  read (lutin, nml=input_teledyn)
 
  write(*,*) 'lta', lta
  write(*,*) 'text_teledyn', text_teledyn
  write(*,*) ' alphadd ',alphadd
  write(*,*) ' massadd ',massadd
  write(6,*) 'kappa  ',kappa
  write(6,*) 'kappaE ',kappaE
  write(6,*) 'igerme ',igerme
  write(6,*)
  write(6,*)'************ DEBUT DE TELEDYN ****************'
  write(6,*)
  write(6,*)
  write(6,*) 'im imm',im,imm
  write(6,*) 'lljc',lljc
  
  do ic=1,igerme
    call random_number(xalea)
  enddo
  write(6,*) ' xalea ', xalea

  if (lljc) then 
    write(*,*) 'coucouc lljc  ',lljc
    !write (*,*) 'avant',xp(1:3,1:7)
   xpepe=xp
  
    call init_ljc(xp,vp,fp)
 

   !  write (*,*) 'apres', xp(1:3,1:7)
    write(*,*) 'potistadd ',potistadd
    !call calfoljq4(xp,fp)
    write(*,*) 'potistadd ',potistadd,potistaddE
    write(*,*) 'coucouc ' 
   !call loop_ljc(xp, vp, fp)
   !call neufdoublewell(xpepe)   
   !call largedeviations(xpepe)  
   !call LWDTPS(xpepe)
    call LyapLanczos(xpepe)
   
  write (*,*) 'voilà les clones'
    stop
  endif
 
 return
 
567 print *,'Erreur lors de la lecture du fichier .tin, Pas de tin pas de teledyn!'


#else  
 
!namelist /input_teledyn/ lta,text_teledyn,alphadd

 write(6,*)'************ DEBUT DE tele_vacancy ****************'
 write(6,*)
 write(6,*)
 write(6,*)
 write(6,*)
 write(6,*) 'im imm',im,imm
  lta=.true.
  
igerme = 5368
  !xtempmin = 1.0e-9
  
 ! lutin=54
  fnamtin = fnam(1:lenfnam)//'.tin'
  ! variables de dynamique
  write(*,*) 'file name', fnamtin

itab=10
 text_teledyn=0.0
  !open(unit=lutin, file=fnamtin, status='unknown', err=567)
 ! read (lutin, nml=input_teledyn)
 
 do ic=1,igerme
    call random_number(xalea)
  enddo
  write(6,*) ' xalea ', xalea
 write(*,*) 'coucouc   ' 

 xpepe=xp

  write (*,*) xp(1,1:3)
 call allocate_tele_vac
! call ndm_into_depart  (xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
write (*,*) xp(1,1:3)
 call index_premier_voisin(xp)

write (*,*) xp(1,1:3)
 call distance_premier_voisin(xp)

 write(*,*) 'coucouc' 

write (*,*) xp(1,1:3)
call init_tele_vac  (xp, vp, fp, ielat, iwmax, ityp)

write (*,*) xp(1,1:3)
 call distance_premier_voisin (xp) 
write (*,*) xp(1,1:3)
!stop
! write(*,*) 'potist0 potistadd ',potist0/(bk*text_tele_vac),potistadd/(bk*text_tele_vac)

 call LyapLanczos_vac !(xp)
stop

#endif

end subroutine  sundae
