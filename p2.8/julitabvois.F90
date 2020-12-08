module calfojuli_mod
        USE notperiod_mod,only: notperiod
        USE cryst_to_cart_mod,only: cryst_to_cart
        USE calfocommon
        USE gen_com_m, ONLY:nvat,fnemd,lcalcjq,lnemd,lperiod,&
             zero
        implicit none
      contains


!----------------------------------------------------------------------
SUBROUTINE calfojuli(im,imm,xp, fp, iwmax, ityp,indi,at,bg,volu)
  !tentaive de calfoeam avec une seule grande boucle sur i
  USE T_kind_param_m

  USE var_pot, ONLY:ipotentiel,potisglue,potisrep,rhomax,rhomin,rue_pot,ngrid,eamrho,ipo,npair,eamglue,typ_pot_pair,eamrep
  USE SMjuli
  USE jqmod
  implicit none

  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  ! eam variables
  integer,intent(in)::im,imm
  integer , intent(in),allocatable :: iwmax(:),ityp(:),indi(:)
  real(double),intent(inout),allocatable  :: xp(:,:)
  real(double) , intent(inout),allocatable :: fp(:,:)
    real(double),intent(in),dimension(3,3)::at,bg
    real(double),intent(in)::volu
 
  !local variables
  integer :: i,j,l !atomes
  integer ::iti,itj,itl,ic !types
  integer :: ll !paires
  integer :: k ! position dans les splines
  real(double) :: rk, drk,ktor, ktorho !pour splines
  real(double) :: cv(1,3)


  real(double) :: rue2 !coupure**2
  real(double) :: Erep,dErep ! potentiel et gradient de la repulsion de paire ij
  real(double) :: dDensityi,dDensityj,gradDensityi(3) ! gradient de la densité sur i et j
  real(double) :: dEembi,dEembj,Eembi ! potentiel et gradient de l'immersion
  real(double) :: rhoi,rhoj,drhoj,rho ! densite de i sur j et j sur i
  real(double) :: rholsi,drholsi,rholsj,drholsj! densite de i sur j et j sur i
  real(double) :: densityi !densite totale sur i


  integer ::nvi                     !nb de voisins de i
  integer, dimension(nvat) :: jvi   !indices j des voisins de i
  real(double), dimension (nvat) ::Vrij, & !distances 
       Vc1ij,Vc2ij,Vc3ij,& ! xpi-xpj
       rhojsi,drhojsi,& ! rho de j sur i et dérivée
       sij,&            ! fonction intermédiaire dans d'écrantage
       ecrsij,&         !écrantage 
       rhotildjsi        ! rho écrantée


  real(double) :: a1i,a2i,a3i,a1j,a2j,a3j !delta x y z
  real(double) :: r2ij,r2il,r2jl !distance carree  i-j i-l j-k
  real(double) :: rij,c1ij,c2ij,c3ij ! distance et delta X ij
  real(double) :: ril,c1il,c2il,c3il ! distance et delta X il
  real(double) :: rjl,c1jl,c2jl,c3jl ! distance et delta X jl
  real(double) :: costetlij, costetijl
  real(double) :: sijl
  real(double) :: aux1,aux2(3),aux3(3),aux4(3) ! aux. pour la force
  real(double) :: gradij(3), gradil(3),gradjl(3) ! deltaX/r pour ij,il et lj
  integer::  iw1,iw2, iw,iw1j,iw2j, iwj, iwl1,iwl2,iwl ! indices des voisins j et l de i



  !    real(double) ::  alphaPbeta,beta
  real(double) :: rcut2 (npair)
  real(double) :: rhoitot,fpi
  real(double) :: tdepcos

  real(double)::rue
  real(double), dimension(:,:), allocatable :: xpnp


  real(double):: fpnemd(3,im),fpnemdmoy(3), XijdotF,XildotF,XjldotF

  rue=rue_pot(ipotentiel)
  fpnemd=0

  do l=1,npair
     rcut2(l)=(reppairjl(l)%rc*1.0d-8)**2
     !       write(6,*)l,sqrt(rcut2(l))
  end do

  ktor=rue/ngrid
  ktorho=(rhomax-rhomin)/ngrid
!  fp(:,:) = 0.0

  jq(:)=0.
  !  if (allocated(eat)) eat(:)=0.


  potist = zero
  potisrep=0.;potisglue=0.
  rue2=rue**2
  !    iw2=0

  ALLOCATE(xpnp(3,imm))
  if (lperiod) then
   xpnp(:,:)=xp(:,:)
  else
   call notperiod(imm,xp,xpnp,at,bg)
  end if
   
  call cryst_to_cart (imm, xpnp, bg, -1)    !cart vers cryst


  loop1at1: do i=1,im

     ! --- Calcul de la densite sur i ---    
     nvi=0.
     rhoitot=0.
     iti = ityp(i)
!            write(6,*)'i iti = ',i,iti
     densityi=0.0 ;Eembi=0.0; dEembi=0.0
     if (i==1) then
        iw1 = 1
     else
        iw1=iwmax(i-1)+1
     end if
     iw2 = iwmax(i)
     loopvois :do iw = iw1, iw2
        j = indi(iw)
        itj=ityp(j)


        c1ij = xpnp(1,i)-xpnp(1,j)
        c2ij = xpnp(2,i)-xpnp(2,j)
        c3ij = xpnp(3,i)-xpnp(3,j)
        if (c1ij>0.5) c1ij = c1ij-1.
        if (c1ij<(-0.5)) c1ij = c1ij+1.
        if (c2ij>0.5) c2ij = c2ij-1.
        if (c2ij<(-0.5)) c2ij = c2ij+1.
        if (c3ij>0.5) c3ij = c3ij-1.
        if (c3ij<(-0.5)) c3ij = c3ij+1.
        cv(1,1) = c1ij
        cv(1,2) = c2ij
        cv(1,3) = c3ij
        call cryst_to_cart (1, cv, at, 1) !cryst vers cart sur cv
        r2ij = cv(1,1)*cv(1,1)+cv(1,2)*cv(1,2)+cv(1,3)*cv(1,3)             
        c1ij= cv(1,1) 
        c2ij =cv(1,2) 
        c3ij =cv(1,3) 

        itj=ityp(j)
        if (r2ij>rcut2(ipo(iti,itj))) cycle

        rij=sqrt(r2ij)
        !          write(6,*)'i iti j  itj r ',i,iti,j,itj
        k=Int(rij/ktor)
        !       write(6,*)'k ',k
        drk=rij-k*ktor
        ll = ipo(iti,itj)
        rhoj=eamrho(1,ll,k)+eamrho(2,ll,k)*drk+eamrho(3,ll,k)*drk**2+eamrho(4,ll,k)*drk**3  !rho de j sur i
        drhoj=eamrho(2,ll,k)+2.0*eamrho(3,ll,k)*drk+3.0*eamrho(4,ll,k)*drk**2  !rho' de j sur i

        !          write(6,*)'i j rij rh',i,j,rij,rhoj
        !          if (iti==1.and.itj==1) then
        !             trh=exp(czz*(dzz-rij/1.0d-8)+kzz/(rij/1.0d-8-rczz))
        !          else
        !             trh=exp(2.0*(czc*(dzc-rij/1.0d-8)+kzc/(rij/1.0d-8-rczc)))
        !          end if
        !          write(6,*)'i jrij trh',i,j,rij,trh



        nvi=nvi+1
        Vc1ij(nvi)=c1ij ; Vc2ij(nvi)=c2ij ; Vc3ij(nvi)=c3ij
        jvi(nvi)=j
        Vrij(nvi)=rij
        rhojsi(nvi)=rhoj
        drhojsi(nvi)=drhoj
        rhoitot=rhoitot+rhojsi(nvi)
     end do loopvois

     !boucle sur les voisins de i
     !       write(6,*)'i nvi', i,nvi
     sij(:)=0.
     loopvj1 : do iw=1,nvi
        j=jvi(iw)
        itj=ityp(j)
        !          if((itj==2).and.(iti==2))cycle
        if(ityp(j).ne.iti) then
           itj=ityp(j)
           !        write(6,*)'i j dis ',i,j,Vrij(iw)
           c1ij=Vc1ij(iw) ;c2ij= Vc2ij(iw); c3ij=Vc3ij(iw)
           rij=Vrij(iw)


           !calcul terme l voisins de i
           !voisins l de i de type itj
           loopvli1 : do iwl=1,nvi
              l=jvi(iwl)
              itl=ityp(l)
              !
              IF (L==J)cycle
              !
              if(itl.eq.iti) cycle 
              !                write(6,*)'L i iti j itj l itl ',i,iti,j,itj,l,itl
              c1il=Vc1ij(iwl) ;c2il= Vc2ij(iwl); c3il=Vc3ij(iwl)
              ril=Vrij(iwl)
              rholsi=rhojsi(iwl)

              costetlij=(c1ij*c1il+c2ij*c2il+c3ij*c3il)/(rij*ril)             
              !                write(6,*)'costetlij',costetlij
              tdepcos=1.0+costetlij
              if (tdepcos.le.0.0) tdepcos=1.0d-10
              !                write(6,*)'tdepcos',tdepcos
              sijl=rholsi*(tdepcos**beta)/alphaPbeta          
              sij(iw)=sij(iw)+sijl
              !                write(6,*)'Lsijl',i,iw,j,l,sij(iw),sijl
           end do loopvli1

           ! fin terme l debut terme "k"

           !initialisations pour voisins de j

           if (j==1) then
              iw1j=1
           else
              iw1j=iwmax(j-1)+1
           endif
           iw2j=iwmax(j)

           loopvjk1 :do iwl = iw1j, iw2j          
              l = indi(iwl)
              !                write(6,*)'i j l',i,j,l
              itl=ityp(l)
              if(itl==itj) cycle 
              !
              IF (L==I)cycle
              !

              c1jl = xpnp(1,j)-xpnp(1,l)
              c2jl = xpnp(2,j)-xpnp(2,l)
              c3jl = xpnp(3,j)-xpnp(3,l)
              if (c1jl>0.5) c1jl = c1jl-1.
              if (c1jl<(-0.5)) c1jl = c1jl+1.
              if (c2jl>0.5) c2jl = c2jl-1.
              if (c2jl<(-0.5)) c2jl = c2jl+1.
              if (c3jl>0.5) c3jl = c3jl-1.
              if (c3jl<(-0.5)) c3jl = c3jl+1.
              cv(1,1) = c1jl
              cv(1,2) = c2jl
              cv(1,3) = c3jl
              call cryst_to_cart (1, cv, at, 1) !cryst vers cart sur cv
              r2jl = cv(1,1)*cv(1,1)+cv(1,2)*cv(1,2)+cv(1,3)*cv(1,3)             
              c1jl=cv(1,1)
              c2jl=cv(1,2) 
              c3jl=cv(1,3)



              if (r2jl>rcut2(ipo(itj,itl))) cycle               
              !                write(6,*)'K i iti j itj l itl ',i,iti,j,itj,l,itl                
              !                if (r2jl>rue2) cycle
              rjl=sqrt(r2jl)
              !       write(6,*)'k ',k
              k=Int(rjl/ktor)
              drk=rjl-k*ktor
              ll = ipo(itj,itl)
              rholsj=eamrho(1,ll,k)+eamrho(2,ll,k)*drk+eamrho(3,ll,k)*drk**2+eamrho(4,ll,k)*drk**3  !rho de l sur j
              costetijl=-(c1ij*c1jl+c2ij*c2jl+c3ij*c3jl)/(rij*rjl)
              !                write(6,*)'costetijk',costetijl
              tdepcos=1+costetijl

              if (tdepcos.le.0.0) tdepcos=1.0d-10
              !                write(6,*)'tdepcos',tdepcos
              sijl=rholsj*(tdepcos**beta)/alphaPbeta          

              !                sijl=rholsj*((1+costetijl)**beta)/alphaPbeta          
              sij(iw)=sij(iw)+sijl 
              !                write(6,*)'Ksijl',i,iw,j,l,sij(iw),sijl
           end do loopvjk1
           if (abs(sij(iw)/rhojsi(iw)).gt.50.0) then
              ecrsij(iw)=0.0
           else
              ecrsij(iw)=exp(-2.0*sqrt(sij(iw)/rhojsi(iw)))
           end if
           rhotildjsi(iw)=rhojsi(iw)*ecrsij(iw)
           !if (i==344)  write(6,*)'ecr ',i,j,sij(iw),rhojsi(iw),sij(iw)/rhojsi(iw),ecrsij(iw)
        else
           sij(iw)=0.0
           ecrsij(iw)=1.0

           rhotildjsi(iw)=rhojsi(iw)*ecrsij(iw)
        end if

        densityi=densityi+rhotildjsi(iw)

     end do loopvj1


     ! calcul et stockage de Eembi et dEembi
     k=Int((densityi-rhomin)/ktorho)
     k=max(k,3) ; k=min(k,ngrid-3)
     !       write(6,*)'i rhoitot densityi ',i,rhoitot,densityi
     drk=densityi-(rhomin+k*ktorho)
     !       write(6,*)i,iti,k
     Eembi=eamglue(1,iti,k)+eamglue(2,iti,k)*drk+eamglue(3,iti,k)*drk**2+eamglue(4,iti,k)*drk**3
     dEembi=eamglue(2,iti,k)+2.0*eamglue(3,iti,k)*drk+3.0*eamglue(4,iti,k)*drk**2
     !       if(i==344)write(6,*)'dembi ', dEembi,k,drk,eamglue(2,iti,k),eamglue(3,iti,k),eamglue(4,iti,k)
     !             write(6,*)'i Eembi', i ,Eembi
     !       call extrapolateEam(embtyp(iti),density(i),Embf=Eembi, dembF=dEembi)
     !densityi, Eembi et dEembi sont des scalaires associés au i courant

     !******************************************************************************************************************

     !boucle des forces
     !    write(6,*)
     !    write(6,*)'REP'
     !    write(6,*)
     loopvj2 :do iw = 1, nvi 
        j = jvi(iw)
        c1ij=Vc1ij(iw) ;c2ij= Vc2ij(iw); c3ij=Vc3ij(iw)
        rij=Vrij(iw)
        gradij(1)=c1ij/rij ; gradij(2)=c2ij/rij ; gradij(3)=c3ij/rij
        itj=ityp(j)

        k=Int(rij/ktor)
        drk=rij-k*ktor
        if (typ_pot_pair(ll).ne.12) cycle
        if (i.gt.j) then  !terme de repulsion deja calculé
           ll = ipo(iti,itj)

           !             Erep=0. ; dErep=0.
           Erep=2.0*eamrep(1,ll,k)+eamrep(2,ll,k)*drk+eamrep(3,ll,k)*drk**2+eamrep(4,ll,k)*drk**3
           dErep=2.0*(eamrep(2,ll,k)+2.0*eamrep(3,ll,k)*drk+3.0*eamrep(4,ll,k)*drk**2)
           !             write(6,*)'i iti j itj rij Erep',i,iti,j,itj,rij,Erep
           !             write(6,*)
           !          write(6,*)'i iti j  itj r ',i,iti,j,itj
           !          write(6,*)'i j rij erep',i,j,rij,erep
           !          if (iti==1.and.itj==1) then
           !             trh=exp(azz*(bzz-rij/1.0d-8)+kzz/(rij/1.0d-8-rczz))*ev2erg
           !          else
           !             trh=exp((azc*(bzc-rij/1.0d-8)+kzc/(rij/1.0d-8-rczc)))*ev2erg
           !          end if
           !          write(6,*)'i jrij terp',i,j,rij,trh
           

           potist=potist+Erep
           potisrep=potisrep+Erep
           fp(1:3,i)=fp(1:3,i)-dErep*gradij(1:3)
           fp(1:3,j)=fp(1:3,j)+dErep*gradij(1:3)
           if (lnemd) then
              XijdotF=c1ij*Fnemd
              do ic=1,3
                 fpnemd(ic,i)=fpnemd(ic,i)+0.5*dErep*gradij(ic)*XijdotF
                 fpnemd(ic,j)=fpnemd(ic,j)+0.5*dErep*gradij(ic)*XijdotF
              end do
           end if

           sig(1:3,1) = sig(1:3,1) -dErep*gradij(1:3)*c1ij/volu
           sig(1:3,2) = sig(1:3,2) -dErep*gradij(1:3)*c2ij/volu
           sig(1:3,3) = sig(1:3,3) -dErep*gradij(1:3)*c3ij/volu
     


           ! A commenter qd lcalcjq=false pour ne pas perdre de temps dans le test
!!$           if (lcalcjq) then
!!$              jqf=0.0
!!$              eat(i)=eat(i)+0.5*Erep
!!$              eat(j)=eat(j)+0.5*Erep
!!$              do ic=1,3
!!$                 jqf=jqf-0.5*(dErep*gradij(ic)*(vp(ic,i)+vp(ic,j)))
!!$              end do
!!$              do ic=1,3
!!$                 jq(ic)=jq(ic)+jqf*gradij(ic)*rij
!!$              end do
!!$           end if


           !if ((i==1).or.(j==1))  write(6,*)'f1 ',fp(1,1),fp(2,1),fp(3,1)
        end if  !i>j

        !          if((itj==2).and.(iti==2))cycle

        rhoj=rhojsi(iw)
        drhoj=drhojsi(iw)
        if(itj==iti) then
           ! terme standard
           fp(1:3,i)=fp(1:3,i)-dEembi*drhoj*gradij(1:3)
           fp(1:3,j)=fp(1:3,j)+dEembi*drhoj*gradij(1:3)

           sig(1:3,1) = sig(1:3,1) -dEembi*drhoj*gradij(1:3)*c1ij/volu
           sig(1:3,2) = sig(1:3,2) -dEembi*drhoj*gradij(1:3)*c2ij/volu
           sig(1:3,3) = sig(1:3,3) -dEembi*drhoj*gradij(1:3)*c3ij/volu

           if (lnemd) then
              XijdotF=c1ij*Fnemd
              do ic=1,3
                 fpnemd(ic,i)=fpnemd(ic,i) +0.5*dEembi*drhoj*gradij(ic)*XijdotF
                 fpnemd(ic,j)=fpnemd(ic,j) +0.5*dEembi*drhoj*gradij(ic)*XijdotF
              end do
           end if

!!$
!!$           if(lcalcjq) then
!!$              jqf=0.0
!!$              do ic=1,3
!!$                 jqf=jqf+dEembi*drhoj*gradij(ic)*vp(ic,j)
!!$              end do
!!$              do ic=1,3
!!$                 jq(ic)=jq(ic)-jqf*gradij(ic)*rij
!!$              end do
!!$           end if

           !if ((i==1).or.(j==1))  write(6,*)'f 2 ',fp(1,1),fp(2,1),fp(3,1)
           !             auxt=dEembi*drhoj*gradij(1)
           !             write(6,*)'auxt 1 ', auxt
        else
           !terme quasi-standard sans d(sij)
           if (ecrsij(iw).ne.0) then
              aux1=ecrsij(iw)*(1.0+sqrt(sij(iw)/rhoj))
              fp(1:3,i)=fp(1:3,i)-dEembi*aux1*drhoj*gradij(1:3)
              fp(1:3,j)=fp(1:3,j)+dEembi*aux1*drhoj*gradij(1:3)

              sig(1:3,1) = sig(1:3,1) -dEembi*aux1*drhoj*gradij(1:3)*c1ij/volu
              sig(1:3,2) = sig(1:3,2) -dEembi*aux1*drhoj*gradij(1:3)*c2ij/volu
              sig(1:3,3) = sig(1:3,3) -dEembi*aux1*drhoj*gradij(1:3)*c3ij/volu

              if (lnemd) then
                 XijdotF=c1ij*Fnemd
                 do ic=1,3
                    fpnemd(ic,i)=fpnemd(ic,i) +0.5*dEembi*aux1*drhoj*gradij(ic)*XijdotF
                    fpnemd(ic,j)=fpnemd(ic,j) +0.5*dEembi*aux1*drhoj*gradij(ic)*XijdotF
                 end do
              end if


!!$              if(lcalcjq) then
!!$                 jqf=0.0
!!$                 do ic=1,3
!!$                    jqf=jqf+dEembi*aux1*drhoj*gradij(ic)*vp(ic,j)
!!$                 end do
!!$                 do ic=1,3
!!$                    jq(ic)=jq(ic)-jqf*gradij(ic)*rij
!!$                 end do
!!$              end if


           end if
        end if

        !terme d(sij)
        if(itj==iti) cycle
        if (ecrsij(iw)==0) cycle
        !       
        !          if(i==344) write(6,*)'i iw j sij(iw) ',i,iw,j,sij(iw)
        aux1=-1.0*dEembi*ecrsij(iw)*sqrt(rhoj/sij(iw))/alphaPbeta  !chgt de signe par rapport à d(E)-> force
        !          if(i==344) write(6,*)'aux1',aux1,dEembi,ecrsij(iw),rhoj,sij(iw)
        !calcul terme l voisins de i
        !voisins l de i de type itj
        loopvli2 : do iwl=1,nvi
           l=jvi(iwl)
           itl=ityp(l)

           if(itl.eq.iti) cycle 
           IF (L==J)cycle
           c1il=Vc1ij(iwl) ;c2il= Vc2ij(iwl); c3il=Vc3ij(iwl)
           ril=Vrij(iwl)
           rholsi=rhojsi(iwl)             
           drholsi=drhojsi(iwl)
           costetlij=(c1ij*c1il+c2ij*c2il+c3ij*c3il)/(rij*ril)             
           gradil(1)=c1il/ril ; gradil(2)=c2il/ril ; gradil(3)=c3il/ril             
           tdepcos=1.0+costetlij
           if (tdepcos.le.0.0) tdepcos=1.0d-10

           aux2(1:3)=((tdepcos)**beta)*ril*drholsi*gradil(1:3)
           aux3(1:3)=beta*((tdepcos)**(beta-1.))*rholsi*(gradij(1:3)-costetlij*gradil(1:3))
           aux4(1:3)=beta*((tdepcos)**(beta-1.))*rholsi*(gradil(1:3)-costetlij*gradij(1:3))

           fp(1:3,i)=fp(1:3,i)-aux1*(aux2(1:3)+aux3(1:3))/ril
           fp(1:3,l)=fp(1:3,l)+aux1*(aux2(1:3)+aux3(1:3))/ril
           fp(1:3,i)=fp(1:3,i)-aux1*aux4(1:3)/rij
           fp(1:3,j)=fp(1:3,j)+aux1*aux4(1:3)/rij



           if (lnemd) then
              XijdotF=c1ij*Fnemd
              XildotF=c1il*Fnemd
              do ic=1,3
                 fpnemd(ic,i)=fpnemd(ic,i)+0.5*(aux1*aux4(ic)/rij)*XijdotF
                 fpnemd(ic,j)=fpnemd(ic,j)+0.5*(aux1*aux4(ic)/rij)*XijdotF
                 fpnemd(ic,i)=fpnemd(ic,i)+0.5*(aux1*(aux2(ic)+aux3(ic))/ril)*XildotF
                 fpnemd(ic,l)=fpnemd(ic,l)+0.5*(aux1*(aux2(ic)+aux3(ic))/ril)*XildotF

              end do
           end if



           
           sig(1:3,1) = sig(1:3,1) -(aux1*(aux2(1:3)+aux3(1:3))/ril)*c1il/volu
           sig(1:3,2) = sig(1:3,2) -(aux1*(aux2(1:3)+aux3(1:3))/ril)*c2il/volu
           sig(1:3,3) = sig(1:3,3) -(aux1*(aux2(1:3)+aux3(1:3))/ril)*c3il/volu
           
           sig(1:3,1) = sig(1:3,1) -(aux1*aux4(1:3)/rij )*c1ij/volu
           sig(1:3,2) = sig(1:3,2) -(aux1*aux4(1:3)/rij )*c2ij/volu
           sig(1:3,3) = sig(1:3,3) -(aux1*aux4(1:3)/rij )*c3ij/volu


!!$           if(lcalcjq) then
!!$              jqf=0.0
!!$              do ic=1,3
!!$                 jqf=jqf+aux1*vp(ic,l)*(aux2(ic)+aux3(ic))/ril
!!$              end do
!!$              do ic=1,3
!!$                 jq(ic)=jq(ic)-jqf*gradil(ic)*ril
!!$              end do
!!$              jqf=0.0
!!$              do ic=1,3
!!$                 jqf=jqf+aux1*vp(ic,j)*aux4(ic)/rij
!!$              end do
!!$              do ic=1,3
!!$                 jq(ic)=jq(ic)-jqf*gradij(ic)*rij
!!$              end do
!!$           end if

        end do loopvli2


        ! fin terme l debut terme "k"
        !initialisations pour voisins de j
        if (j==1) then
           iw1j=1
        else
           iw1j=iwmax(j-1)+1
        endif
        iw2j=iwmax(j)             
        loopvjk2 :do iwl = iw1j, iw2j          
           l = indi(iwl)
           itl=ityp(l)
           IF (L==I)cycle
           if(itl==itj) cycle                 
           c1jl = xpnp(1,j)-xpnp(1,l)
           c2jl = xpnp(2,j)-xpnp(2,l)
           c3jl = xpnp(3,j)-xpnp(3,l)
           if (c1jl>0.5) c1jl = c1jl-1.
           if (c1jl<(-0.5)) c1jl = c1jl+1.
           if (c2jl>0.5) c2jl = c2jl-1.
           if (c2jl<(-0.5)) c2jl = c2jl+1.
           if (c3jl>0.5) c3jl = c3jl-1.
           if (c3jl<(-0.5)) c3jl = c3jl+1.
           cv(1,1) = c1jl
           cv(1,2) = c2jl
           cv(1,3) = c3jl
           call cryst_to_cart (1, cv, at, 1) !cryst vers cart sur cv
           r2jl = cv(1,1)*cv(1,1)+cv(1,2)*cv(1,2)+cv(1,3)*cv(1,3)             
           c1jl=cv(1,1)
           c2jl=cv(1,2) 
           c3jl=cv(1,3)                   

           if (r2jl>rcut2(ipo(itj,itl))) cycle               
           !                if (r2jl>rue2) cycle
           rjl=sqrt(r2jl)
           gradjl(1)=c1jl/rjl ; gradjl(2)=c2jl/rjl ; gradjl(3)=c3jl/rjl
           !       write(6,*)'k ',k
           k=Int(rjl/ktor)
           drk=rjl-k*ktor
           ll = ipo(itj,itl)
           rholsj=eamrho(1,ll,k)+eamrho(2,ll,k)*drk+eamrho(3,ll,k)*drk**2+eamrho(4,ll,k)*drk**3  !rho de l sur 
           drholsj=eamrho(2,ll,k)+2.0*eamrho(3,ll,k)*drk+3.0*eamrho(4,ll,k)*drk**2  !rho de l sur j



           costetijl=-(c1ij*c1jl+c2ij*c2jl+c3ij*c3jl)/(rij*rjl)
           tdepcos=1.0+costetijl
           if (tdepcos.le.0.0) tdepcos=1.0d-10

           aux2(1:3)=((tdepcos)**beta)*rjl*drholsj*gradjl(1:3)
           aux3(1:3)=beta*((tdepcos)**(beta-1.))*rholsj*(-1.*gradij(1:3)-costetijl*gradjl(1:3))
           aux4(1:3)=-1.0*beta*((tdepcos)**(beta-1.))*rholsj*(gradjl(1:3)+costetijl*gradij(1:3))

           fp(1:3,j)=fp(1:3,j)-aux1*(aux2(1:3)+aux3(1:3))/rjl
           fp(1:3,l)=fp(1:3,l)+aux1*(aux2(1:3)+aux3(1:3))/rjl
           fp(1:3,i)=fp(1:3,i)-aux1*aux4(1:3)/rij
           fp(1:3,j)=fp(1:3,j)+aux1*aux4(1:3)/rij

           if (lnemd) then
              XijdotF=c1ij*Fnemd
              XjldotF=c1jl*Fnemd
              do ic=1,3
                 fpnemd(ic,j)=fpnemd(ic,j)+0.5*(aux1*(aux2(ic)+aux3(ic))/rjl)*XjldotF
                 fpnemd(ic,l)=fpnemd(ic,l)+0.5*(aux1*(aux2(ic)+aux3(ic))/rjl)*XjldotF
                 fpnemd(ic,i)=fpnemd(ic,i)+0.5*(aux1*aux4(ic)/rij)*XijdotF
                 fpnemd(ic,j)=fpnemd(ic,j)+0.5*(aux1*aux4(ic)/rij)*XijdotF

              end do
           end if




           sig(1:3,1) = sig(1:3,1) -(aux1*(aux2(1:3)+aux3(1:3))/rjl)*c1jl/volu
           sig(1:3,2) = sig(1:3,2) -(aux1*(aux2(1:3)+aux3(1:3))/rjl)*c2jl/volu
           sig(1:3,3) = sig(1:3,3) -(aux1*(aux2(1:3)+aux3(1:3))/rjl)*c3jl/volu


           sig(1:3,1) = sig(1:3,1) -(aux1*aux4(1:3)/rij )*c1ij/volu
           sig(1:3,2) = sig(1:3,2) -(aux1*aux4(1:3)/rij )*c2ij/volu
           sig(1:3,3) = sig(1:3,3) -(aux1*aux4(1:3)/rij )*c3ij/volu

!!$
!!$           if(lcalcjq) then
!!$              jqf=0.0
!!$              do ic=1,3
!!$                 jqf=jqf+(aux1*aux4(ic)/rij-aux1*(aux2(ic)+aux3(ic)/rjl))*vp(ic,j)
!!$              end do
!!$              do ic=1,3
!!$                 jq(ic)=jq(ic)-jqf*gradij(ic)*rij
!!$              end do
!!$              jqf=0.0
!!$              do ic=1,3
!!$                 jqf=jqf+(aux1*(aux2(ic)+aux3(ic))/rjl)*vp(ic,l)
!!$              end do
!!$              do ic=1,3
!!$                 jq(ic)=jq(ic)-jqf*(gradij(ic)*rij+gradjl(ic)*rjl)
!!$              end do
!!$           end if


        end do loopvjk2
        !if (i==1)  write(6,*)'fb1 ',fp(1,1),fp(2,1),fp(3,1)
     end do loopvj2
     !if (i==1)  write(6,*)'fb2 ',fp(1,1),fp(2,1),fp(3,1)

     potist=potist+Eembi
     if (lcalcjq) eat(i)=eat(i)+Eembi
     potisglue=potisglue+Eembi

     !    write(6,*)


  end do loop1at1
  if (lnemd) then
     fpnemdmoy=0
     do i=1,im   
        do l=1,3
           fpnemdmoy(l)=fpnemdmoy(l)+fpnemd(l,i)/float(im)
        enddo
     end do

     do i=1,im
!        write(6,*)'A',i,fp(:,i)
        do l=1,3
           fp(l,i)=fp(l,i)-fpnemdmoy(l)
           fp(l,i)=fp(l,i)+fpnemd(l,i)
        enddo
!        write(6,*)'B',i,fp(:,i)
     end do
  end if






  deALLOCATE(xpnp)

  !  write(6,*)'f8 ',fp(1,1),fp(2,1),fp(3,1)
  return
end SUBROUTINE calfojuli


end module
