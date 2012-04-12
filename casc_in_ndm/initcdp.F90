! **************************************************************
subroutine initcdp
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m, iseed_glob=>iseed
  use tab_imm_m
  use defcdp
  implicit none
  integer, dimension(2) :: iseedt
  !-----------------------------------------------
  !   L o c a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: i,itapp
  !-----------------------------------------------
  namelist /inputcdp/itecdp,nintrodp,imin,imax,nposI,iseed,dminins,typint,ioxdef,Ed,ideftyp,rsphdef,centresphdef



  itecdp=-1      ! introduction de DP tout les itecdp pas
  ideftyp=0      ! type de def ; 0=PF; 1=int
  nintrodp=1    ! nombre de DP intrduit à chaque fois
  imin=1        ! indice minimal possible pour les atomes déplacés
  imax=im       ! indice maximal possible pour les atomes déplacés
  nposI=-1      ! nombre de positions interstitielles
  iseed=-1      ! graine pour la generation aleatoire si <0 tirage avec SECNDS
  dminins=1.3   ! distance minimum entre nouvel interstitiel et atomes deja present
  typint=0     ! type d'introduction des Intestitiels : 0 dans les sites prédéfinis, 1 aléatoirement
  ioxdef=0
  Ed(:)=0.
  rsphdef=0.
  centresphdef(:)=0.5

  open(unit=73, file='creaDPin', status='unknown')
  read (73, nml=inputcdp)
  dminins=dminins*1d-8
  rsphdef=rsphdef*1d-8


  if((typint.lt.0).or.(typint.GT.1)) then
     write(6,*)'mauvaise introduction des interstitiels stop'
     call arret_ndm
  end if

  if(ideftyp==1) then
     write(6,*)'introduction d interstitiels de type 1'
  end if

  if (typint==0) then


     allocate(xposint(3,nposI))
     do i=1,nposI
        read(73,*)xposint(1,i),xposint(2,i),xposint(3,i)
     end do


     where (xposint(:,:)<0.0)
        xposint(:,:)=xposint(:,:)+1.
     end where
     where (xposint(:,:)>1.0)
        xposint(:,:)=xposint(:,:)-1.
     end where

     !         do i=1,nposI
     !            xposint(:,i)=(xposint(:,i)-0.5)*zl(:)
     !         end do
  end if
  if (iseed.le.0) then
     call system_clock (iseed) 
     write(6,*)'rang iseed ',rang,iseed
  end if
  iseedt(1)=iseed
  call random_seed (iseedt(1))
! call    random_seed (put=iseedt)

  return
end subroutine initcdp
