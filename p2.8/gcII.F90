! *************************************************************
subroutine gcII(xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use var_pot
  use work_cgII
  use tab_imm_m, ONLY : bruitmd
  ! *************************************************************
  ! xp positions des atomes
  ! xpp previous positions
  ! vp  velocities
  ! ax starting positions
  ! fp forces
  !iwmax =  indice du dernier voisin de chaque atome
  !ityp  tableau des types
  implicit none
  !-----------------------------------------------
  !   G l o b a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  integer  :: ielat(imm)
  integer  :: iwmax(imm)
  integer  :: ityp(imm)
  real(double)  :: xp(3,imm)
  real(double)  :: xpp(3,imm)
  real(double)  :: vp(3,imm)
  real(double)  :: ax(3,imm)
  real(double)  :: fp(3,imm)
  !-----------------------------------------------
  !   L o c a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: n,  i, igc
  real(double) :: efinal
  !GC settings
  integer  :: NGC,criterion,NCALLS,IER
  double precision, dimension(:), allocatable:: X,G,W
  double precision     :: ACC,F

  !-----------------------------------------------
  !
  !



     imd = im
     nad(:ntyp) = na(:ntyp)

   if (mdcg_noise /= 0 ) then
     call bruit_xp
   end if 


  !New GC settings ....:
  NGC=3*imm

  allocate (X(NGC),G(NGC),W(6*NGC))


          do i=1,im
             IF (dmtype.EQ.30) THEN
                     ! Variables = reduced coordinates
                   if (mdcg_noise==0) then
                     
                     X(3*i-2:3*i) = MatMul( ax(1:3,i), bg)
                   else
                     xp(1:3,i)= ax(1:3,i)+bruitmd(1:3,i)
                     X(3*i-2:3*i) = MatMul( xp(1:3,i), bg)
                   end if
             ELSE
                     ! Variables = cartesian coordinates (in A)
                     if (mdcg_noise==0) then
                      X(3*i-2:3*i) = ax(1:3,i)*angst
                     else
                      X(3*i-2:3*i) = (ax(1:3,i)+bruitmd(1:3,i))*angst
                     end if
             END IF
          end do
          X(3*im+1:3*imm)=0.d0

  ! internal units
  !      epsilon_force=epsilon_force/6.251d3
  criterion=0
  !!$ACC=1.d-8
  ACC=0.d0


  it = 0
  efinal = 0.0

  CALL ZXCGRII(FUNCT,NGC,ACC,itmax,X,G,F,W,IER,criterion,NCALLS, &
       xp, xpp, vp, ax, fp,  ielat, iwmax, ityp)       

  deallocate (X,G,W)
  
  write(*,*) 'DEBUG after GCII ... just SAY HALLO'


#ifdef ART
  write (*, *) 'energie ', potist,'erg',potist*erg2eV,'eV'
  return

#else
  call endrun 
  write (*, *) 'energie ', potist,'erg',potist*erg2eV,'eV'
  return
#endif

end subroutine gcII


