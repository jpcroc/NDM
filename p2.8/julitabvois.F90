module calfojuli_mod
  USE calfocommon
  USE gen_com_m, ONLY:nvat,fnemd,lcalcjq,lnemd,lperiod, zero
  use vect_dist_mod,only:vect_dist
  USE atomconfig,only : atom_config,atom_config_d,atom_config_e
  USE cellconfig, only : cell_config
  use boxconfig,only: box_config

        implicit none
      contains


!----------------------------------------------------------------------
SUBROUTINE calfojuli(atcf,celcf,boxcf)
  !tentaive de calfoeam avec une seule grande boucle sur i
  USE T_kind_param_m

  USE var_pot, ONLY:ipotentiel,potisglue,potisrep,rhomax,rhomin,rue_pot,ngrid,eamrho,ipo,npair,eamglue,typ_pot_pair,eamrep,ntyp
  USE SMjuli
  USE jqmod
  implicit none
  class(atom_config),intent(inout)::atcf
  type(cell_config),intent(in)::celcf
  type(box_config),intent(in)::boxcf
 
  !local variables
  integer :: i,j,l !atomes
  integer ::iti,itj,itl,ic !types
  integer :: ll !paires
  integer :: k ! position dans les splines
  real(double) :: rk, drk,ktor
  real(double),allocatable::ktorho(:) !pour splines
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
  real(double) :: rcut2 (npair),rcut(npair)
  real(double) :: rhoitot,fpi
  real(double) :: tdepcos
  real(double)::rue
  real(double):: fpnemd(3,atcf%im),fpnemdmoy(3), XijdotF,XildotF,XjldotF
  real(double)::dxp(3),dxpjl(3)
  logical::linter,linterjl

  rue=rue_pot(ipotentiel)
  fpnemd=0

  do l=1,npair
     rcut2(l)=(reppairjl(l)%rc*1.0d-8)**2
     rcut(l)=(reppairjl(l)%rc*1.0d-8)
  end do
  allocate(ktorho(ntyp))
  ktor=rue/ngrid
  ktorho(:)=(rhomax(:)-rhomin(:))/ngrid

!  fp(:,:) = 0.0

  jq(:)=0.
  !  if (allocated(eat)) eat(:)=0.
  potist = zero
  potisrep=0.;potisglue=0.
  rue2=rue**2
  !    iw2=0
  loop1at1: do i=1,atcf%im

     ! --- Calcul de la densite sur i ---    
     nvi=0.
     rhoitot=0.
     iti = atcf%ityp(i)
!            write(6,*)'i iti = ',i,iti
     densityi=0.0 ;Eembi=0.0; dEembi=0.0
     if (i==1) then
        iw1 = 1
     else
        iw1=atcf%iwmax(i-1)+1
     end if
     iw2 = atcf%iwmax(i)
     loopvois :do iw = iw1, iw2
        j = atcf%indi(iw)
        itj=atcf%ityp(j)
        ll=ipo(iti,itj)
        call vect_dist(atcf,celcf,boxcf,i,j,VJI=dxp, lperiod=boxcf%lperiod,rum=rcut(ll)&
             &,linter=linter,dist=rij)
        if (.not.linter) cycle
        
        !          write(6,*)'i iti j  itj r ',i,iti,j,itj
        k=Int(rij/ktor)
        !       write(6,*)'k ',k
        drk=rij-k*ktor
        ll = ipo(iti,itj)
        rhoj=eamrho(1,ll,k)+eamrho(2,ll,k)*drk+eamrho(3,ll,k)*drk**2+eamrho(4,ll,k)*drk**3  !rho de j sur i
        drhoj=eamrho(2,ll,k)+2.0*eamrho(3,ll,k)*drk+3.0*eamrho(4,ll,k)*drk**2  !rho' de j sur i

        nvi=nvi+1
        Vc1ij(nvi)=dxp(1) ; Vc2ij(nvi)=dxp(2) ; Vc3ij(nvi)=dxp(3)
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
        itj=atcf%ityp(j)
        !          if((itj==2).and.(iti==2))cycle
        if(atcf%ityp(j).ne.iti) then
           itj=atcf%ityp(j)
           !        write(6,*)'i j dis ',i,j,Vrij(iw)
           c1ij=Vc1ij(iw) ;c2ij= Vc2ij(iw); c3ij=Vc3ij(iw)
           rij=Vrij(iw)


           !calcul terme l voisins de i
           !voisins l de i de type itj
           loopvli1 : do iwl=1,nvi
              l=jvi(iwl)
              itl=atcf%ityp(l)
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
              iw1j=atcf%iwmax(j-1)+1
           endif
           iw2j=atcf%iwmax(j)

           loopvjk1 :do iwl = iw1j, iw2j          
              l = atcf%indi(iwl)
              !                write(6,*)'i j l',i,j,l
              itl=atcf%ityp(l)
              if(itl==itj) cycle 
              IF (L==I)cycle
              call vect_dist(atcf,celcf,boxcf,j,l,VJI=dxpjl, lperiod=boxcf%lperiod&
                   &,rum=rcut(ipo(itj,itl)),linter=linterjl,dist=rjl)
              if (.not.linterjl)cycle
              c1jl=dxpjl(1);              c2jl=dxpjl(2);              c3jl=dxpjl(3);

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
     k=Int((densityi-rhomin(iti))/ktorho(iti))
     k=max(k,3) ; k=min(k,ngrid-3)
     !       write(6,*)'i rhoitot densityi ',i,rhoitot,densityi
     drk=densityi-(rhomin(iti)+k*ktorho(iti))
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
        itj=atcf%ityp(j)

        k=Int(rij/ktor)
        drk=rij-k*ktor
        ll = ipo(iti,itj)
        if (typ_pot_pair(ll).ne.12) cycle
        if (i.gt.j) then  !terme de repulsion deja calculé


           !             Erep=0. ; dErep=0.
           Erep=2.0*eamrep(1,ll,k)+eamrep(2,ll,k)*drk+eamrep(3,ll,k)*drk**2+eamrep(4,ll,k)*drk**3
           dErep=2.0*(eamrep(2,ll,k)+2.0*eamrep(3,ll,k)*drk+3.0*eamrep(4,ll,k)*drk**2)
           potist=potist+Erep
           potisrep=potisrep+Erep
           atcf%fp(1:3,i)=atcf%fp(1:3,i)-dErep*gradij(1:3)
           atcf%fp(1:3,j)=atcf%fp(1:3,j)+dErep*gradij(1:3)
           if (lnemd) then
              XijdotF=c1ij*Fnemd
              do ic=1,3
                 fpnemd(ic,i)=fpnemd(ic,i)+0.5*dErep*gradij(ic)*XijdotF
                 fpnemd(ic,j)=fpnemd(ic,j)+0.5*dErep*gradij(ic)*XijdotF
              end do
           end if

           sig(1:3,1) = sig(1:3,1) -dErep*gradij(1:3)*c1ij/boxcf%volu
           sig(1:3,2) = sig(1:3,2) -dErep*gradij(1:3)*c2ij/boxcf%volu
           sig(1:3,3) = sig(1:3,3) -dErep*gradij(1:3)*c3ij/boxcf%volu
        end if  !i>j

        !          if((itj==2).and.(iti==2))cycle

        rhoj=rhojsi(iw)
        drhoj=drhojsi(iw)
        if(itj==iti) then
           ! terme standard
           atcf%fp(1:3,i)=atcf%fp(1:3,i)-dEembi*drhoj*gradij(1:3)
           atcf%fp(1:3,j)=atcf%fp(1:3,j)+dEembi*drhoj*gradij(1:3)

           sig(1:3,1) = sig(1:3,1) -dEembi*drhoj*gradij(1:3)*c1ij/boxcf%volu
           sig(1:3,2) = sig(1:3,2) -dEembi*drhoj*gradij(1:3)*c2ij/boxcf%volu
           sig(1:3,3) = sig(1:3,3) -dEembi*drhoj*gradij(1:3)*c3ij/boxcf%volu

           if (lnemd) then
              XijdotF=c1ij*Fnemd
              do ic=1,3
                 fpnemd(ic,i)=fpnemd(ic,i) +0.5*dEembi*drhoj*gradij(ic)*XijdotF
                 fpnemd(ic,j)=fpnemd(ic,j) +0.5*dEembi*drhoj*gradij(ic)*XijdotF
              end do
           end if
        else
           !terme quasi-standard sans d(sij)
           if (ecrsij(iw).ne.0) then
              aux1=ecrsij(iw)*(1.0+sqrt(sij(iw)/rhoj))
              atcf%fp(1:3,i)=atcf%fp(1:3,i)-dEembi*aux1*drhoj*gradij(1:3)
              atcf%fp(1:3,j)=atcf%fp(1:3,j)+dEembi*aux1*drhoj*gradij(1:3)

              sig(1:3,1) = sig(1:3,1) -dEembi*aux1*drhoj*gradij(1:3)*c1ij/boxcf%volu
              sig(1:3,2) = sig(1:3,2) -dEembi*aux1*drhoj*gradij(1:3)*c2ij/boxcf%volu
              sig(1:3,3) = sig(1:3,3) -dEembi*aux1*drhoj*gradij(1:3)*c3ij/boxcf%volu

              if (lnemd) then
                 XijdotF=c1ij*Fnemd
                 do ic=1,3
                    fpnemd(ic,i)=fpnemd(ic,i) +0.5*dEembi*aux1*drhoj*gradij(ic)*XijdotF
                    fpnemd(ic,j)=fpnemd(ic,j) +0.5*dEembi*aux1*drhoj*gradij(ic)*XijdotF
                 end do
              end if

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
           itl=atcf%ityp(l)

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

           atcf%fp(1:3,i)=atcf%fp(1:3,i)-aux1*(aux2(1:3)+aux3(1:3))/ril
           atcf%fp(1:3,l)=atcf%fp(1:3,l)+aux1*(aux2(1:3)+aux3(1:3))/ril
           atcf%fp(1:3,i)=atcf%fp(1:3,i)-aux1*aux4(1:3)/rij
           atcf%fp(1:3,j)=atcf%fp(1:3,j)+aux1*aux4(1:3)/rij



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
           
           sig(1:3,1) = sig(1:3,1) -(aux1*(aux2(1:3)+aux3(1:3))/ril)*c1il/boxcf%volu
           sig(1:3,2) = sig(1:3,2) -(aux1*(aux2(1:3)+aux3(1:3))/ril)*c2il/boxcf%volu
           sig(1:3,3) = sig(1:3,3) -(aux1*(aux2(1:3)+aux3(1:3))/ril)*c3il/boxcf%volu
           
           sig(1:3,1) = sig(1:3,1) -(aux1*aux4(1:3)/rij )*c1ij/boxcf%volu
           sig(1:3,2) = sig(1:3,2) -(aux1*aux4(1:3)/rij )*c2ij/boxcf%volu
           sig(1:3,3) = sig(1:3,3) -(aux1*aux4(1:3)/rij )*c3ij/boxcf%volu


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
           iw1j=atcf%iwmax(j-1)+1
        endif
        iw2j=atcf%iwmax(j)             
        loopvjk2 :do iwl = iw1j, iw2j          
           l = atcf%indi(iwl)
           itl=atcf%ityp(l)
           IF (L==I)cycle
           if(itl==itj) cycle                 

           call vect_dist(atcf,celcf,boxcf,j,l,VJI=dxpjl, lperiod=boxcf%lperiod,rum=rcut(ipo(itj,itl))&
                &,linter=linterjl,dist=rjl)
              if (.not.linterjl)cycle
              c1jl=dxpjl(1);              c2jl=dxpjl(2);              c3jl=dxpjl(3);
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

           atcf%fp(1:3,j)=atcf%fp(1:3,j)-aux1*(aux2(1:3)+aux3(1:3))/rjl
           atcf%fp(1:3,l)=atcf%fp(1:3,l)+aux1*(aux2(1:3)+aux3(1:3))/rjl
           atcf%fp(1:3,i)=atcf%fp(1:3,i)-aux1*aux4(1:3)/rij
           atcf%fp(1:3,j)=atcf%fp(1:3,j)+aux1*aux4(1:3)/rij

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




           sig(1:3,1) = sig(1:3,1) -(aux1*(aux2(1:3)+aux3(1:3))/rjl)*c1jl/boxcf%volu
           sig(1:3,2) = sig(1:3,2) -(aux1*(aux2(1:3)+aux3(1:3))/rjl)*c2jl/boxcf%volu
           sig(1:3,3) = sig(1:3,3) -(aux1*(aux2(1:3)+aux3(1:3))/rjl)*c3jl/boxcf%volu


           sig(1:3,1) = sig(1:3,1) -(aux1*aux4(1:3)/rij )*c1ij/boxcf%volu
           sig(1:3,2) = sig(1:3,2) -(aux1*aux4(1:3)/rij )*c2ij/boxcf%volu
           sig(1:3,3) = sig(1:3,3) -(aux1*aux4(1:3)/rij )*c3ij/boxcf%volu

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
!     if (lcalcjq) eat(i)=eat(i)+Eembi
     potisglue=potisglue+Eembi

     !    write(6,*)


  end do loop1at1
  if (lnemd) then
     fpnemdmoy=0
     do i=1,atcf%im   
        do l=1,3
           fpnemdmoy(l)=fpnemdmoy(l)+fpnemd(l,i)/float(atcf%im)
        enddo
     end do

     do i=1,atcf%im
!        write(6,*)'A',i,fp(:,i)
        do l=1,3
           atcf%fp(l,i)=atcf%fp(l,i)-fpnemdmoy(l)
           atcf%fp(l,i)=atcf%fp(l,i)+fpnemd(l,i)
        enddo
!        write(6,*)'B',i,fp(:,i)
     end do
  end if








  !  write(6,*)'f8 ',fp(1,1),fp(2,1),fp(3,1)
  return
end SUBROUTINE calfojuli


end module
