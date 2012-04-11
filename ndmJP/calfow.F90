! ***************************************************************
      SUBROUTINE CALFOW(xp,xpp,vp,ax,fp,ielat,iwmax,ityp)
!     calcule des forces a 2 corps dans le pot de Watanabe
!     version Avril 2001  
! ***************************************************************

       USE T_kind_param_m 
       use gen_com_m 
!      include 'gen.com'
!      include 'pos.com'

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
!--------------------------------------------
!  L o c a l   v a r i a b l e s
!-------------------------------------------
      integer iw2,iti,l,iw1,i,koo,i1,ko1,j, &
     & i2,itj,k,ic,kz,&
     & ncelvois,	&	!nonbre de cellules =0 ou 26
     & itO,itSi         ! types du O et du Si

      parameter(itSi=1,itO=2)

      real(double) :: xpnp(3,imm)

      real(double) SCALAR(imM),z(imM),dzdx(imM,3),spotr(imM),&
     & virdzdx(imM,3,3),cv(1,3)
      real(double) wclau,aux,alp,f1,f2,f3,&
      x1,x2,x3,auy,auz,c1,c2,c3,ddsq,r2,sk,dk,r,&
      phu,ra,sz,&
      hbn2,&		!||g||**2
      dfp,fdp,&
      potr,foncgz,dpotr,dfoncgz,dzdxpart   ! intermediaires de calcul 

! spline
      real(double) dr,potpart,dz     

! sig=0
!       write(6,*)'entree calfw'

      if (any(free).NEQV..true.)then
         write(6,*)'free +SW =pas code'
         stop
      end if

      if (lperiod) then
         xpnp(:,:)=xp(:,:)
      else 
         call notperiod(xp,xpnp)
      end if


       if (itesigma.gt.0) then
         if (mod(it,itesigma).eq.0) then
         do i1=1,3
          do i2=1,3
            do i=1,im
               virdzdx(i,i1,i2)=0.0
!               sigat(i,i1,i2)=0.0
            enddo
            sig(i1,i2)=0.0
                 if (lTPcel.EQV..true.) then
                    do koo=1,noxyz
                       sigc(i1,i2,koo)=0.0
                       sigc(i1,i2,koo)=0.0
                    enddo
                 end if
          enddo
         enddo
         endif
       endif
   
! Epot=0      
      potis1=zero
      POTCP=0.0 
      POTIST=ZERO
! F=0/
      DO 3333 I=1,Im
      z(i)=ZERO
      spotr(i)=ZERO
      dzdx(i,1)=ZERO
      dzdx(i,2)=ZERO
      dzdx(i,3)=ZERO
      FP(1,i)=ZERO
      FP(2,i)=ZERO
 3333 FP(3,i)=ZERO

      AUX=23.06134575D-20 
      ALP=ALPHA/SQRT(PI)*AUX
      IW2=0 

!      if (iterdf.gt.0)  nrdf=nrdf+1



! calcul de la cordination des atomes d'O (z)
! plus some des fonctions de SW pour les forces qui suivent

      do i=1,imd
!         write(6,*)i
      iti=ityp(i)
      if(iti.eq.itO) then

      koo=ielat(i)   ! Numero de la cellule
      ncelvois=min(noxyz,27)-1
      do 11 i1=0,ncelvois
         ko1=ncel(koo,i1)
!         write(6,*)'ko1',ko1
         c1p = xpnp(1,i)+sum(at(1,:)*deltadist(:,i1,koo))
         c2p = xpnp(2,i)+sum(at(2,:)*deltadist(:,i1,koo))
         c3p = xpnp(3,i)+sum(at(3,:)*deltadist(:,i1,koo))

         do 12 i2=1,nato(ko1)
            j=last(i2,ko1)
!         write(6,*)'j ', j
            itj=ityp(j)
            if(itj.eq.itSi)then
!         write(6,*)'i j 1er ',i,j
!               c1=xpnp(1,i)-xpnp(1,j)
!               c2=xpnp(2,i)-xpnp(2,j)
!               c3=xpnp(3,i)-xpnp(3,j)

!               if (ltriclin) then
!                  do ic=1,3
!                     c1=c1+at(1,ic)*deltadist(ic,i1,koo)
!                     c2=c2+at(2,ic)*deltadist(ic,i1,koo)
!                     c3=c3+at(3,ic)*deltadist(ic,i1,koo)
!                  enddo
!               else
!                  c1=c1+zl(1)*deltadist(1,i1,koo)
!                  c2=c2+zl(2)*deltadist(2,i1,koo)
!                  c3=c3+zl(3)*deltadist(3,i1,koo)
!               endif

           if (noxyz.ne.1) then
              c1 = c1p-xpnp(1,j)
              c1abs=abs(c1)
              if(c1abs>rawat(l)) cycle
              
              c2 = c2p-xpnp(2,j)
              c2abs=abs(c2)
              if(c2abs>rawat(l)) cycle
              c3 = c3p-xpnp(3,j)
              c3abs=abs(c3)
              if(c3abs>rawat(l)) cycle
              cv(1,1) = c1
              cv(1,2) = c2
              cv(1,3) = c3

           else
              c1 = c1p-xpnp(1,j)
              c2 = c2p-xpnp(2,j)
              c3 = c3p-xpnp(3,j)
              cv(1,1) = c1
              cv(1,2) = c2
              cv(1,3) = c3
              call cryst_to_cart (1, cv, bg, -1) !cryst vers cart sur cv
              WHERE ( (cv.GT.0.5d0).OR.(cv.LT.-0.5d0) )
                 cv(:,1:3) = cv(:,1:3) - Dble(Nint(cv(:,1:3)))
              END WHERE
              call cryst_to_cart (1, cv, at, 1) !cryst vers cart sur cv
              c1=cv(1,1)
              c2=cv(1,2)
              c3=cv(1,3)

           end if


               r2=c1*c1+c2*c2+c3*c3
               l=IPO(iti,itj)
               if(r2.gt.rawat(l)**2) goto 12 ! r2<r_coupure2
               r=sqrt(r2)

               sk=r/csive
               k=sk
               dr=r-float(k)*csive
!         write(6,*)'i j 2eme ',i,j
!         write(6,*)'r sk k dr ',r,sk,k,dr
               z(i)= z(i)+ fcr(k)+    &                ! calul de la coordination
                bspf(k)*dr+cspf(k)*dr**2+dspf(k)*dr**3
               spotr(i)= spotr(i)+     &               ! calcul de la somme des termes de SW
               potw(l,k)+bspw(l,k)*dr+cspw(l,k)*dr**2+dspw(l,k)*dr**3
!         write(6,*)'tata'
               dzdxpart=(bspf(k)+2.0*cspf(k)*dr+3.0*dspf(k)*dr**2)/r
               dzdx(i,1)= dzdx(i,1) +dzdxpart *c1
               dzdx(i,2)= dzdx(i,2) +dzdxpart *c2
               dzdx(i,3)= dzdx(i,3) +dzdxpart *c3
!         write(6,*)'toto'
!               virdzdx(i,1,1)= virdzdx(i,1,1) +dzdxpart *c1 *c1
!               virdzdx(i,1,2)= virdzdx(i,1,2) +dzdxpart *c1 *c2
!               virdzdx(i,1,3)= virdzdx(i,1,3) +dzdxpart *c1 *c3
!               virdzdx(i,2,1)= virdzdx(i,2,1) +dzdxpart *c2 *c1
!               virdzdx(i,2,2)= virdzdx(i,2,2) +dzdxpart *c2 *c2
!               virdzdx(i,2,3)= virdzdx(i,2,3) +dzdxpart *c2 *c3
!               virdzdx(i,3,1)= virdzdx(i,3,1) +dzdxpart *c3 *c1
!               virdzdx(i,3,2)= virdzdx(i,3,2) +dzdxpart *c3 *c2
!               virdzdx(i,3,3)= virdzdx(i,3,3) +dzdxpart *c3 *c3
!         write(6,*)'titi'
            endif
   12 continue
   11 continue
      endif
      enddo


! boucle de calcul des forces

      DO 699 I=1,IMD
!         write(6,*)i
      KOO=IELAT(I)   ! Numero de la cellule
      ITI=ITYP(I)
      ncelvois=min(noxyz,27)-1

      DO 61 I1=0,ncelvois
      KO1=NCEL(KOO,I1)
      do ic=1,3
         c1=c1+at(1,ic)*deltadist(ic,i1,koo)
         c2=c2+at(2,ic)*deltadist(ic,i1,koo)
         c3=c3+at(3,ic)*deltadist(ic,i1,koo)
      enddo

      
      DO 62 I2=1,NATO(KO1)
      j=last(i2,ko1)
      if(i.eq.j) goto 62
  !    C1=xpnp(1,i)-xpnp(1,J)
  !    C2=xpnp(2,i)-xpnp(2,J)
  !    C3=xpnp(3,i)-xpnp(3,J)

      if (noxyz.ne.1) then
         c1 = c1p-xpnp(1,j)
         c1abs=abs(c1)
         if(c1abs>rawat(l)) cycle
         
         c2 = c2p-xpnp(2,j)
         c2abs=abs(c2)
         if(c2abs>rawat(l)) cycle
         c3 = c3p-xpnp(3,j)
         c3abs=abs(c3)
         if(c3abs>rawat(l)) cycle
         cv(1,1) = c1
         cv(1,2) = c2
         cv(1,3) = c3
         
      else
         c1 = c1p-xpnp(1,j)
         c2 = c2p-xpnp(2,j)
         c3 = c3p-xpnp(3,j)
         cv(1,1) = c1
         cv(1,2) = c2
         cv(1,3) = c3
         call cryst_to_cart (1, cv, bg, -1) !cryst vers cart sur cv
         WHERE ( (cv.GT.0.5d0).OR.(cv.LT.-0.5d0) )
            cv(:,1:3) = cv(:,1:3) - Dble(Nint(cv(:,1:3)))
         END WHERE
         call cryst_to_cart (1, cv, at, 1) !cryst vers cart sur cv
         c1=cv(1,1)
         c2=cv(1,2)
         c3=cv(1,3)
         
      end if



      R=sqrt(C1*C1+C2*C2+C3*C3)
      ITJ=ITYP(J)
      L=IPO(ITI,ITJ)
!      if(r.gt.rue) goto 62
      if(r.gt.rawat(l)) goto 62
      SK=R/csive
      K=SK

! spline

      dr=r-float(k)*csive
!      write(6,*)l,k
      potr=potw(l,k)+&                          ! partie Stillinger Weber (r)
     &     bspw(l,k)*dr+cspw(l,k)*dr**2+dspw(l,k)*dr**3
      dpotr=bspw(l,k)+2.0*cspw(l,k)*dr+3.0*dspw(l,k)*dr**2 ! derivation de la partie SW

      if(iti.ne.itj)then                      ! liaison SiO
         if(iti.eq.itO) then                  ! O est l'atome considere
           sz=z(i)/csive_g                    ! on prend le z de O
           kz=sz                              ! voila un entier
           dz=z(i)-float(kz)*csive_g 
           foncgz=  gz(kz)+     &              ! g(z)
           bspg(kz)*dz+cspg(kz)*dz**2+dspg(kz)*dz**3
           dfoncgz= bspg(kz)+2.0*cspg(kz)*dz+3.0*dspg(kz)*dz**2 ! derivee de g
           potis1=potis1+0.5*foncgz*potr ! energie potentielle
           dfp=dfoncgz*potr
           fdp=foncgz*dpotr/r
           F1=-1.0*(dzdx(i,1)*dfp+fdp*c1) ! force x
           F2=-1.0*(dzdx(i,2)*dfp+fdp*c2) ! force y
           F3=-1.0*(dzdx(i,3)*dfp+fdp*c3) ! force z
         else
           dzdxpart=bspf(k)+2.0*cspf(k)*dr+3.0*dspf(k)*dr**2
           sz=z(j)/csive_g                      ! on prend le z de O 
           kz=sz                                ! voila un entier
           dz=z(j)-float(kz)*csive_g
           foncgz=  gz(kz)+  &                   ! g(z)
             bspg(kz)*dz+cspg(kz)*dz**2+dspg(kz)*dz**3
           dfoncgz= bspg(kz)+2.0*cspg(kz)*dz+3.0*dspg(kz)*dz**2 ! derivee de g
           potis1=potis1+0.5*foncgz*potr        ! energie potentielle  
           F1=-1.0*((dzdxpart*dfoncgz*spotr(j)+foncgz*dpotr)*c1/r) ! force x
           F2=-1.0*((dzdxpart*dfoncgz*spotr(j)+foncgz*dpotr)*c2/r) ! force y
           F3=-1.0*((dzdxpart*dfoncgz*spotr(j)+foncgz*dpotr)*c3/r) ! force z        
        endif
      else                      ! autres liaison Si-Si ou O-O
         potis1=potis1+0.5*potr
         phu=-1.0*dpotr/r
         F1=PHU*C1
         F2=PHU*C2
         F3=PHU*C3           
      endif
      FP(1,i)=FP(1,i)+F1
      FP(2,i)=FP(2,i)+F2
      FP(3,i)=FP(3,i)+F3
 !     if(i.eq.4) write(6,*)'it,fp1itesigma=',it,fp(1,i)

! calcul des contraintes
       if (itesigma.gt.0) then
         if (mod(it,itesigma).eq.0) then
            if(iti.ne.itj.and.iti.eq.itO) then ! O est l'atome considere
               sig(1,1)=sig(1,1)-0.5*(virdzdx(i,1,1)*dfp+fdp*c1*c1)/volu
               sig(1,2)=sig(1,2)-0.5*(virdzdx(i,1,2)*dfp+fdp*c1*c2)/volu
               sig(1,3)=sig(1,3)-0.5*(virdzdx(i,1,3)*dfp+fdp*c1*c3)/volu
               sig(2,1)=sig(2,1)-0.5*(virdzdx(i,2,1)*dfp+fdp*c2*c1)/volu
               sig(2,2)=sig(2,2)-0.5*(virdzdx(i,2,2)*dfp+fdp*c2*c2)/volu
               sig(2,3)=sig(2,3)-0.5*(virdzdx(i,2,3)*dfp+fdp*c2*c3)/volu
               sig(3,1)=sig(3,1)-0.5*(virdzdx(i,3,1)*dfp+fdp*c3*c1)/volu
               sig(3,2)=sig(3,2)-0.5*(virdzdx(i,3,2)*dfp+fdp*c3*c2)/volu
               sig(3,3)=sig(3,3)-0.5*(virdzdx(i,3,3)*dfp+fdp*c3*c3)/volu
               !contraintes atomiques statique 
!               sigat(i,1,1)=sigat(i,1,1)-&
!                   0.5*(virdzdx(i,1,1)*dfp+fdp*c1*c1)
!               sigat(i,1,2)=sigat(i,1,2)-&
!                   0.5*(virdzdx(i,1,2)*dfp+fdp*c1*c2)
!               sigat(i,1,3)=sigat(i,1,3)-&
!                   0.5*(virdzdx(i,1,3)*dfp+fdp*c1*c3)
!               sigat(i,2,1)=sigat(i,2,1)-&
!                   0.5*(virdzdx(i,2,1)*dfp+fdp*c2*c1)
!               sigat(i,2,2)=sigat(i,2,2)-&
!                   0.5*(virdzdx(i,2,2)*dfp+fdp*c2*c2)
!               sigat(i,2,3)=sigat(i,2,3)-&
!                   0.5*(virdzdx(i,2,3)*dfp+fdp*c2*c3)
!               sigat(i,3,1)=sigat(i,3,1)-&
!                   0.5*(virdzdx(i,3,1)*dfp+fdp*c3*c1)
!               sigat(i,3,2)=sigat(i,3,2)-&
!                   0.5*(virdzdx(i,3,2)*dfp+fdp*c3*c2)
!               sigat(i,3,3)=sigat(i,3,3)-&
!                   0.5*(virdzdx(i,3,3)*dfp+fdp*c3*c3)
               !contraintes locales
                 if (lTPcel.EQV..true.) then
                    sigc(1,1,koo)=sigc(1,1,koo)-&
                         0.5*(virdzdx(i,1,1)*dfp+fdp*c1*c1)*noxyz/volu
                    sigc(1,2,koo)=sigc(1,2,koo)-&
                         0.5*(virdzdx(i,1,2)*dfp+fdp*c1*c2)*noxyz/volu
                    sigc(1,3,koo)=sigc(1,3,koo)-&
                         0.5*(virdzdx(i,1,3)*dfp+fdp*c1*c3)*noxyz/volu
                    sigc(2,1,koo)=sigc(2,1,koo)-&
                         0.5*(virdzdx(i,2,1)*dfp+fdp*c2*c1)*noxyz/volu
                    sigc(2,2,koo)=sigc(2,2,koo)-&
                         0.5*(virdzdx(i,2,2)*dfp+fdp*c2*c2)*noxyz/volu
                    sigc(2,3,koo)=sigc(2,3,koo)-&
                         0.5*(virdzdx(i,2,3)*dfp+fdp*c2*c3)*noxyz/volu
                    sigc(3,1,koo)=sigc(3,1,koo)-&
                         0.5*(virdzdx(i,3,1)*dfp+fdp*c3*c1)*noxyz/volu
                    sigc(3,2,koo)=sigc(3,2,koo)-&
                         0.5*(virdzdx(i,3,2)*dfp+fdp*c3*c2)*noxyz/volu
                    sigc(3,3,koo)=sigc(3,3,koo)-&
                         0.5*(virdzdx(i,3,3)*dfp+fdp*c3*c3)*noxyz/volu
                 end if
            else
               sig(1,1)=sig(1,1)+0.5*F1*c1/volu
               sig(1,2)=sig(1,2)+0.5*F1*c2/volu
               sig(1,3)=sig(1,3)+0.5*F1*c3/volu
               sig(2,1)=sig(2,1)+0.5*F2*c1/volu
               sig(2,2)=sig(2,2)+0.5*F2*c2/volu
               sig(2,3)=sig(2,3)+0.5*F2*c3/volu
               sig(3,1)=sig(3,1)+0.5*F3*c1/volu
               sig(3,2)=sig(3,2)+0.5*F3*c2/volu
               sig(3,3)=sig(3,3)+0.5*F3*c3/volu
               !contraintes atomiques
!               sigat(i,1,1)=sigat(i,1,1)+0.5*F1*c1
!               sigat(i,1,2)=sigat(i,1,2)+0.5*F1*c2
!               sigat(i,1,3)=sigat(i,1,3)+0.5*F1*c3
!               sigat(i,2,1)=sigat(i,2,1)+0.5*F2*c1
!               sigat(i,2,2)=sigat(i,2,2)+0.5*F2*c2
!               sigat(i,2,3)=sigat(i,2,3)+0.5*F2*c3
!               sigat(i,3,1)=sigat(i,3,1)+0.5*F3*c1
!               sigat(i,3,2)=sigat(i,3,2)+0.5*F3*c2
!               sigat(i,3,3)=sigat(i,3,3)+0.5*F3*c3
               !contraintes locales   
               if (lTPcel.EQV..true.) then
                  sigc(1,1,koo)=sigc(1,1,koo)+0.5*F1*c1*noxyz/volu
                  sigc(1,2,koo)=sigc(1,2,koo)+0.5*F1*c2*noxyz/volu
                  sigc(1,3,koo)=sigc(1,3,koo)+0.5*F1*c3*noxyz/volu
                  sigc(2,1,koo)=sigc(2,1,koo)+0.5*F2*c1*noxyz/volu
                  sigc(2,2,koo)=sigc(2,2,koo)+0.5*F2*c2*noxyz/volu
                  sigc(2,3,koo)=sigc(2,3,koo)+0.5*F2*c3*noxyz/volu
                  sigc(3,1,koo)=sigc(3,1,koo)+0.5*F3*c1*noxyz/volu
                  sigc(3,2,koo)=sigc(3,2,koo)+0.5*F3*c2*noxyz/volu
                  sigc(3,3,koo)=sigc(3,3,koo)+0.5*F3*c3*noxyz/volu
               end if
            endif
         endif
       endif
   62 CONTINUE
   61 CONTINUE
  699 CONTINUE


      POTIST=potis1




      return
      end


