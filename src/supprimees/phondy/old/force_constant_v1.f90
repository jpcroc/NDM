! ************************************************
!           Sous-programme force_constant
! ************************************************

subroutine force_constant (xp, xpp, vp, ax, fp,  ielat, iwmax, iwmax2,ityp)
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
  integer  :: iwmax2(imm)
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
  integer :: iw, iw1,iw2,iti
  real(double) :: diff
  real(double) :: deltaE 

  
  real(double),allocatable, dimension(:,:)  :: fpp1,fpp2,fpm1,fpm2
  real(double), dimension(:,:), allocatable :: xpnp
  real(double)::rue,rue2,rtest2,time
  REAL(double), dimension(1:3) :: dxp
  integer :: count1,count2,count_rate,count_max 
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

    rue=rue_pot(ipotentiel)
    rue=2.d0*5.5d0*A2cm
    rue2=rue**2
    write(*,*) 'RUE.....(SHOULD BE  AT LEAST 2 x Rcut ...: ', rue 
    write(*,*) 'RUE.POTENTIEL............................: ', rue_pot(ipotentiel) 


  ALLOCATE(xpnp(3,imm))
  if (lperiod) then
   xpnp(:,:)=xp(:,:)
  else
   call notperiod(xp,xpnp)
  end if
  call cryst_to_cart (imm, xpnp, bg, -1)
  
! Here we apply a coefficient in order to put the Hessian in eV/ A²
  !coeff=
  if (rang==0) write (6, *) 'PHONDY: Computing the force constants .... '

  nad(:ntyp) = na(:ntyp)
  xpref(:,:) = xp(:,:)
 !debug write(*,*) 'x,y,z....2 ', xp(1:3,2)
 !debugCOS write(*,*) 'the mass', umass, cm(ityp(20)), cm(ityp(20))/umass
  call system_clock (count1,count_rate,count_max)
  if (HessianOrder.eq.1) then
 !
  do j=1,im
   iti= ityp(j)
   iw1=iw2+1
   iw2=iwmax(j)
     do jb=1,3
        u=3*(j-1)+jb
         
        xp(jb,j)=xpref(jb,j)+deltax
        call calfo
        
        loopvois: do iw=iw1,iw2
           i=indi(iw)
           do ia=1,3
              v=3*(i-1)+ia
              !( j, jb ; i , ia)
              matfor(u,v)=-1.*(fp(ia,i)-fp0(ia,i))/ &
              (deltax*dsqrt(cm(ityp(i)/umass))*dsqrt(cm(ityp(j)/umass)))
           
             end do
        end do loopvois

     end do
  end do
 !
 end if

!debugCOSwrite(*,*)erg2eV/angst
 matfor(:,:)=0.d0

 if (HessianOrder.eq.2) then
 allocate(fpp1(3,imm),fpm1(3,imm))
 !
  iw2=0
  do j=1,im
  iti=ityp(j)
  iw1=iw2+1
  iw2=iwmax2(j)
   if ((mod (j,400)==0)) write(*,*) '.......j ',j,iw2-iw1,deltax

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
        
          do iw=iw1,iw2
             i=indi2(iw)
  
             dxp(1:3) = xpnp(1:3,i) - xpnp(1:3,j)

             WHERE ( (dxp.GT.0.5d0).OR.(dxp.LT.-0.5d0) )
               dxp(1:3) = dxp(1:3) - Dble(Nint(dxp(1:3)))
             END WHERE
             dxp = MatMul(at,dxp)
             rtest2 = Sum( dxp(1:3)**2 )
            if (rtest2<rue2)  then

             do ia=1,3
              v=3*(i-1)+ia
              !other_good_way_to_pack: v=i+im*(ia-1)
              !( j, jb ; i , ia)
              matfor(u,v)=-1.*(fpp1(ia,i)-fpm1(ia,i))/ ( (2.d0*deltax/ang2cm) &  
              * dsqrt(cm(ityp(i))/umass)*dsqrt(cm(ityp(j))/umass) )
              
              !matfor(v,u)=matfor(u,v)
              !if (dabs(matfor(u,v)).lt.1.d-10) then
              ! matfor(u,v)=0.d0
              !end if
           
             end do
            else
              do ia=1,3 
               v=3*(i-1)+ia
               matfor(u,v)=0.d0
              end do 
            end if
        end do

     end do
  end do
  



 !
 deallocate(fpp1,fpm1)
 end if


 if (HessianOrder.eq.4) then
 allocate(fpp1(3,imm),fpp2(3,imm),fpm1(3,imm), fpm2(3,imm))
 !
  iw2=0
  do j=1,im
  iti=ityp(j)
  iw1=iw2+1
  iw2=iwmax2(j)
   if ((mod (j,100)==0)) write(*,*) '.......j ',j,iw2-iw1,deltax
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

         
        do iw=iw1,iw2
           i=indi2(iw)
            !
             dxp(1:3) = xpnp(1:3,i) - xpnp(1:3,j)
             WHERE ( (dxp.GT.0.5d0).OR.(dxp.LT.-0.5d0) )
               dxp(1:3) = dxp(1:3) - Dble(Nint(dxp(1:3)))
             END WHERE
             dxp = MatMul(at,dxp)
             rtest2 = Sum( dxp(1:3)**2 )
             !
             if (rtest2<rue2)  then
              !
              do ia=1,3
               v=3*(i-1)+ia
               matfor(u,v)=-1.*(fpm2(ia,i)-8.d0*fpm1(ia,i)+8.d0*fpp1(ia,i)-fpp2(ia,i) &
                     )/ ( (12.d0*deltax/ang2cm) &  
               * dsqrt(cm(ityp(i))/umass)*dsqrt(cm(ityp(j))/umass) )
              end do
              !
             else
              ! 
              do ia=1,3
               v=3*(i-1)+ia
               matfor(u,v)=0.d0
              end do
              !
             end if
             !
        end do

     end do
  end do
 !
 deallocate(fpp1,fpp2,fpm1,fpm2)
 end if

  call system_clock (count2,count_rate,count_max)
  time=real((count2-count1))/real(count_rate)
  write(*,"(' PHONDY: MATFOR  was filled in.......:  ',f16.8,' s')") time
  write(*,*) 'PHONDY: MATFOR', matfor (1,2) ,matfor(2,1)
  write(*,*) 'PHONDY: MATFOR', matfor (1,5) ,matfor(5,1)


  do i=1,3*im
    do j=i,3*im
    matfor(i,j)=0.5d0*(matfor(i,j)+matfor(j,i))
    matfor(j,i)=matfor(i,j)
    end do
  end do
  call system_clock (count1,count_rate,count_max)
  time=real((count1-count2))/real(count_rate)
  write(*,"(' PHONDY: MATFOR  was symmetrized in...:  ',f16.8,' s')") time

  write(*,*) 'PHONDY: ...the force constants were  filled'

end subroutine force_constant


