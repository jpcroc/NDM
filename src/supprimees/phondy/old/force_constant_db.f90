! ************************************************
!           Sous-programme force_constant
! ************************************************

subroutine force_constant (xp, xpp, vp, ax, fp,  ielat, iwmax, ityp)
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use var_pot
  use phondy_in_ndm_module
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
  real(double)  :: xpref(3,imm)
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
  integer :: i, ia,j,jb,u,v
  real(double) :: diff
  real(double) :: deltaE 

  real(double),allocatable, dimension(:,:)  :: fpp1,fpp2,fpm1,fpm2
 !-----------------------------------------------
  !
  !
  ! calcul des phonons
  ! construit la matrice dynamique par calcul des forces pour de petits deplacements des atomes autour de positions
  ! d'equilibre. 
  ! diagonalise cette matrice. La diagonalisation sort les carre des pulsations
  ! les frequences sont les racines carrées divisées par 2pi des valeurs propres
  !
  !
  !

  
! Here we apply a coefficient in order to put the Hessian in eV/ A²
  !coeff=
  if (rang==0) write (6, *) 'PHONDY: Computing the force constants .... '

  nad(:ntyp) = na(:ntyp)
  xpref(:,:) = xp(:,:)

 !debugCOS write(*,*) 'the mass', umass, cm(ityp(20)), cm(ityp(20))/umass
  
  if (HessianOrder.eq.1) then
 !
  do j=1,im
     do jb=1,3
        u=3*(j-1)+jb
         
        xp(jb,j)=xpref(jb,j)+deltax
        call calfo
        
        do i=1,im
           do ia=1,3
              v=3*(i-1)+ia
              !( j, jb ; i , ia)
              matfor(u,v)=-1.*(fp(ia,i)-fp0(ia,i))/ &
              (deltax*dsqrt(cm(ityp(i)/umass))*dsqrt(cm(ityp(j)/umass)))
           
             end do
        end do

     end do
  end do
 !
 end if

!debugCOSwrite(*,*)erg2eV/angst

 if (HessianOrder.eq.2) then
 allocate(fpp1(3,imm),fpm1(3,imm))
 !
  do j=1,im
     do jb=1,3
        u=3*(j-1)+jb
        !other_good_way_to_pack: u=j+im*(jb-1) 
        xp(jb,j)=xpref(jb,j)+deltax
        call calfo
        !debugCOS write(*,*)  potist*erg2eV-epot0, epot0,  potist*erg2eV 
        if (dabs(potist*erg2eV-epot0).gt.0.5) then 
         write(*,*) 'WARNING the energy is greater than 0.5. Decrease deltax', (potist*erg2eV-epot0)
        end if
        fpp1(:,:)=fp(:,:)*erg2eV/angst
        xp(jb,j)=xpref(jb,j)-deltax
        call calfo
        if (dabs(potist*erg2eV-epot0).gt.0.5) then 
         write(*,*) 'WARNING the energy is greater than 0.5. Decrease deltax',(potist*erg2eV-epot0)
        end if
        fpm1(:,:)=fp(:,:)*erg2eV/angst
        
        do i=1,im
           do ia=1,3
              v=3*(i-1)+ia
              !other_good_way_to_pack: v=i+im*(ia-1)
              !( j, jb ; i , ia)
              matfor(u,v)=-1.*(fpp1(ia,i)-fpm1(ia,i))/ ( (2.d0*deltax/ang2cm) &  
              * dsqrt(cm(ityp(i))/umass)*dsqrt(cm(ityp(j))/umass) )
              !if (dabs(matfor(u,v)).lt.1.d-10) then
              ! matfor(u,v)=0.d0
              !end if
           
             end do
        end do

     end do
  end do
 !
 deallocate(fpp1,fpm1)
 end if


 if (HessianOrder.eq.4) then
 allocate(fpp1(3,imm),fpp2(3,imm),fpm1(3,imm), fpm2(3,imm))
 !
  do j=1,im
     do jb=1,3
        u=3*(j-1)+jb
        !other_good_way_to_pack u=j+im*(jb-1) 
         
        xp(jb,j)=xpref(jb,j)+deltax
        call calfo
        fpp1(:,:)=fp(:,:)*erg2eV/angst

        xp(jb,j)=xpref(jb,j)+2.d0*deltax
        call calfo
        fpp2(:,:)=fp(:,:) *erg2eV/angst
 

        xp(jb,j)=xpref(jb,j)-deltax
        call calfo
        fpm1(:,:)=fp(:,:)*erg2eV/angst

        xp(jb,j)=xpref(jb,j)-2.d0*deltax
        call calfo
        fpm2(:,:)=fp(:,:)*erg2eV/angst

         
        do i=1,im
           do ia=1,3
              v=3*(i-1)+ia
              !other_good_way_to_pack v=i+im*(ia-1)
              !( j, jb ; i , ia)
              matfor(u,v)=-1.*(fpm2(ia,i)-8.d0*fpm1(ia,i)+8.d0*fpp1(ia,i)-fpp2(ia,i) &
                     )/ ( (12.d0*deltax/ang2cm) &  
              * dsqrt(cm(ityp(i))/umass)*dsqrt(cm(ityp(j))/umass) )

           
             end do
        end do

     end do
  end do
 !
 deallocate(fpp1,fpp2,fpm1,fpm2)
 end if


  write(*,*) 'MATFOR', matfor (1,2) ,matfor(2,1)
  write(*,*) 'MATFOR', matfor (1,5) ,matfor(5,1)
  write(*,*) 'PHONDY: ...the force constants were  filled'

end subroutine force_constant


