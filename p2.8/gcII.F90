! *************************************************************
subroutine gcII(xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use work_cgII
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



  !New GC settings ....:
  !!$NGC=3*imm
  IF (lFrozen) THEN
          NGC = 3*Count(Free(1:im))
  ELSE
          NGC = 3*imm
  END IF

  allocate (X(NGC),G(NGC),W(6*NGC))


  IF (lFrozen) THEN
          iGC=0
          do i=1 ,im
             IF (.Not.Free(i)) Cycle
             iGC=iGC+1
             X(3*iGC-2:3*iGC)=ax(1:3,i)
          end do
          IF (3*iGC.NE.nGC) THEN
                  WRITE(0,'(a,i0)') "3*iGC = ", 3*iGC
                  WRITE(0,'(a,i0)') "nGC   = ", nGC
                  STOP "< gcII >"
          END IF
  ELSE
          do i=1,imm
             X(3*i-2:3*i)=ax(1:3,i)
          end do
  END IF

  ! internal units
  !      epsilon_force=epsilon_force/6.251d3
  criterion=0
  !!$ACC=1.d-8
  ACC=0.d0


  it = 0
  efinal = 0.0

  ! X into from internal units to Ang

!debugGC     do i=1,5
!debugGC      write(*,'("gcII ax(*,1) ",3f15.8)')ax(1,i)/at(1,1),ax(2,i)/at(1,1),ax(3,i)/at(1,1)
!debugGC     end do
!debugGC            write(*,*) '                   '
!debugGC
!debugGC     do i=1,5
!debugGC      write(*,'("gcII xp(*,1) ",3f15.8)')xp(1,i)/at(1,1),xp(2,i)/at(1,1),xp(3,i)/at(1,1)
!debugGC     end do

  X=X*angst
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


