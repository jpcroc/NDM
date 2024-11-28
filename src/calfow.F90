module calfow_mod
  USE arret_ndm_mod,only:arret_ndm
  USE gen_com_m, ONLY: lperiod,pi,potcp,zero
  USE atomconfig,only : atom_config,atom_config_d,atom_config_e
  USE cellconfig, only : cell_config
  use boxconfig,only: box_config
  USE calfocommon
  use vect_dist_mod,only:vect_dist

  implicit none
contains
  ! ***************************************************************
  SUBROUTINE CALFOw(atcf,celcf,boxcf)
    ! ***************************************************************

    USE T_kind_param_m 
    USE var_pot, ONLY:alpha,csive,csive_g,ipo,dspf,ipo,rawat,dspw,dspf,dspg,cspg,cspf,bspf,bspg,fcr,cspw,gz,bspw,potw,potis1
    class(atom_config),intent(inout)::atcf
    type(cell_config),intent(in)::celcf
    type(box_config),intent(in)::boxcf

    !--------------------------------------------
    !  L o c a l   v a r i a b l e s
    !-------------------------------------------
    integer iw2,iti,l,koo,i1,ko1,j, &
         & i2,itj,k,kz,i,&
         & itO,itSi         ! types du O et du Si

    parameter(itSi=1,itO=2)


    real(double) z(atcf%imM),dzdx(atcf%imM,3),spotr(atcf%imM),&
         & virdzdx(atcf%imM,3,3)
    real(double) aux,alp,f1,f2,f3,&
         auy,sk,r,&
         phu,sz,&
         dfp,fdp,&
         potr,foncgz,dpotr,dfoncgz,dzdxpart   ! intermediaires de calcul 

    ! spline
    real(double) dr,dz     

    logical::linter
    real(double)::dxp(3)

    ! sig=0
    !       write(6,*)'entree calfw'


    if (test_sigma)  then
       do i1=1,3
          do i2=1,3
             do i=1,atcf%im
                virdzdx(i,i1,i2)=0.0
                !               sigat(i,i1,i2)=0.0
             enddo
             sigcalfo(i1,i2)=0.0
             if (lcalcsigc.EQV..true.) then
                do koo=1,celcf%noxyz
                   sigc(i1,i2,koo)=0.0
                   sigc(i1,i2,koo)=0.0
                enddo
             end if
          enddo
       enddo
    end if


    ! Epot=0      
    potis1=zero
    POTCP=0.0 
    POTISTcalfo=ZERO
    ! F=0/
    DO  I=1,Atcf%Im
       z(i)=ZERO
       spotr(i)=ZERO
       dzdx(i,1)=ZERO
       dzdx(i,2)=ZERO
       dzdx(i,3)=ZERO
       atcf%FP(1,i)=ZERO
       ATCF%FP(2,i)=ZERO
       ATCF%FP(3,i)=ZERO
    end DO

       AUX=23.06134575D-20 
       ALP=ALPHA/SQRT(PI)*AUX
       IW2=0 

       !      if (iterdf.gt.0)  nrdf=nrdf+1



       ! calcul de la cordination des atomes d'O (z)
       ! plus some des fonctions de SW pour les forces qui suivent

       do i=1,atcf%im
          !         write(6,*)i
          iti=atcf%ityp(i)
          if(iti.eq.itO) then

             koo=atcf%ielat(i)   ! Numero de la cellule
             do  i1=0,celcf%ncelvois(koo)
                ko1=celcf%ncel(koo,i1)
                !         write(6,*)'ko1',ko1

                do  i2=1,celcf%nato(ko1)
                   j=celcf%atincel(i2,ko1)
                   !         write(6,*)'j ', j
                   itj=atcf%ityp(j)
                   if(itj.eq.itSi)then
                      l=ipo(iti,itj)
                      call vect_dist(atcf,celcf,boxcf,i,j,VJI=dxp,indcv=i1, lperiod=boxcf%lperiod,rum=rawat(l),linter=linter,dist=r)
                      if (.not.linter) cycle
                      
                      sk=r/csive
                      k=int(sk)
                      dr=r-float(k)*csive
                      !         write(6,*)'i j 2eme ',i,j
                      !         write(6,*)'r sk k dr ',r,sk,k,dr
                      z(i)= z(i)+ fcr(k)+    &                ! calul de la coordination
                           bspf(k)*dr+cspf(k)*dr**2+dspf(k)*dr**3
                      spotr(i)= spotr(i)+     &               ! calcul de la somme des termes de SW
                           potw(l,k)+bspw(l,k)*dr+cspw(l,k)*dr**2+dspw(l,k)*dr**3
                      !         write(6,*)'tata'
                      dzdxpart=(bspf(k)+2.0*cspf(k)*dr+3.0*dspf(k)*dr**2)/r
                      dzdx(i,1)= dzdx(i,1) +dzdxpart *dxp(1)
                      dzdx(i,2)= dzdx(i,2) +dzdxpart *dxp(2)
                      dzdx(i,3)= dzdx(i,3) +dzdxpart *dxp(3)
                   endif
                end do
             end do
          endif
       enddo


       ! boucle de calcul des forces

       DO 699 I=1,ATCF%IM
          !         write(6,*)i
          KOO=atcf%IELAT(I)   ! Numero de la cellule
          ITI=atcf%ITYP(I)

          DO 61 I1=0,celcf%ncelvois(koo)
             KO1=celcf%NCEL(KOO,I1)


             DO 62 I2=1,celcf%NATO(KO1)
                j=celcf%atincel(i2,ko1)
                if(i.eq.j) goto 62
                ITJ=atcf%ITYP(J)
                L=IPO(ITI,ITJ)
                call vect_dist(atcf,celcf,boxcf,i,j,VJI=dxp,indcv=i1, lperiod=boxcf%lperiod,rum=rawat(l),linter=linter,dist=r)
 
                if (.not.linter) cycle

                SK=R/csive
                K=int(SK)

                ! spline

                dr=r-float(k)*csive
                !      write(6,*)l,k
                potr=potw(l,k)+&                          ! partie Stillinger Weber (r)
                     &     bspw(l,k)*dr+cspw(l,k)*dr**2+dspw(l,k)*dr**3
                dpotr=bspw(l,k)+2.0*cspw(l,k)*dr+3.0*dspw(l,k)*dr**2 ! derivation de la partie SW

                if(iti.ne.itj)then                      ! liaison SiO
                   if(iti.eq.itO) then                  ! O est l'atome considere
                      sz=z(i)/csive_g                    ! on prend le z de O
                      kz=int(sz)                              ! voila un entier)
                      dz=z(i)-float(kz)*csive_g 
                      foncgz=  gz(kz)+     &              ! g(z)
                           bspg(kz)*dz+cspg(kz)*dz**2+dspg(kz)*dz**3
                      dfoncgz= bspg(kz)+2.0*cspg(kz)*dz+3.0*dspg(kz)*dz**2 ! derivee de g
                      potis1=potis1+0.5*foncgz*potr ! energie potentielle
                      dfp=dfoncgz*potr
                      fdp=foncgz*dpotr/r
                      F1=-1.0*(dzdx(i,1)*dfp+fdp*dxp(1)) ! force x
                      F2=-1.0*(dzdx(i,2)*dfp+fdp*dxp(2)) ! force y
                      F3=-1.0*(dzdx(i,3)*dfp+fdp*dxp(3)) ! force z
                   else
                      dzdxpart=bspf(k)+2.0*cspf(k)*dr+3.0*dspf(k)*dr**2
                      sz=z(j)/csive_g                      ! on prend le z de O 
                      kz=sz                                ! voila un entier
                      dz=z(j)-float(kz)*csive_g
                      foncgz=  gz(kz)+  &                   ! g(z)
                           bspg(kz)*dz+cspg(kz)*dz**2+dspg(kz)*dz**3
                      dfoncgz= bspg(kz)+2.0*cspg(kz)*dz+3.0*dspg(kz)*dz**2 ! derivee de g
                      potis1=potis1+0.5*foncgz*potr        ! energie potentielle  
                      F1=-1.0*((dzdxpart*dfoncgz*spotr(j)+foncgz*dpotr)*dxp(1)/r) ! force x
                      F2=-1.0*((dzdxpart*dfoncgz*spotr(j)+foncgz*dpotr)*dxp(2)/r) ! force y
                      F3=-1.0*((dzdxpart*dfoncgz*spotr(j)+foncgz*dpotr)*dxp(3)/r) ! force z        
                   endif
                else                      ! autres liaison Si-Si ou O-O
                   potis1=potis1+0.5*potr
                   phu=-1.0*dpotr/r
                   F1=PHU*DXP(1)
                   F2=PHU*DXP(2)
                   F3=PHU*DXP(3)           
                endif
                ATCF%FP(1,i)=ATCF%FP(1,i)+F1
                ATCF%FP(2,i)=ATCF%FP(2,i)+F2
                ATCF%FP(3,i)=ATCF%FP(3,i)+F3
                !     if(i.eq.4) write(6,*)'it,fp1itesigma=',it,fp(1,i)

                ! calcul des contraintes
                if (test_sigma) then
                   if(iti.ne.itj.and.iti.eq.itO) then ! O est l'atome considere
                      sigcalfo(1,1)=sigcalfo(1,1)-0.5*(virdzdx(i,1,1)*dfp+fdp*dxp(1)*dxp(1))/boxcf%volu
                      sigcalfo(1,2)=sigcalfo(1,2)-0.5*(virdzdx(i,1,2)*dfp+fdp*dxp(1)*dxp(2))/boxcf%volu
                      sigcalfo(1,3)=sigcalfo(1,3)-0.5*(virdzdx(i,1,3)*dfp+fdp*dxp(1)*dxp(3))/boxcf%volu
                      sigcalfo(2,1)=sigcalfo(2,1)-0.5*(virdzdx(i,2,1)*dfp+fdp*dxp(2)*dxp(1))/boxcf%volu
                      sigcalfo(2,2)=sigcalfo(2,2)-0.5*(virdzdx(i,2,2)*dfp+fdp*dxp(2)*dxp(2))/boxcf%volu
                      sigcalfo(2,3)=sigcalfo(2,3)-0.5*(virdzdx(i,2,3)*dfp+fdp*dxp(2)*dxp(3))/boxcf%volu
                      sigcalfo(3,1)=sigcalfo(3,1)-0.5*(virdzdx(i,3,1)*dfp+fdp*dxp(3)*dxp(1))/boxcf%volu
                      sigcalfo(3,2)=sigcalfo(3,2)-0.5*(virdzdx(i,3,2)*dfp+fdp*dxp(3)*dxp(2))/boxcf%volu
                      sigcalfo(3,3)=sigcalfo(3,3)-0.5*(virdzdx(i,3,3)*dfp+fdp*dxp(3)*dxp(3))/boxcf%volu
                      if (lcalcsigc.EQV..true.) then
                         sigc(1,1,koo)=sigc(1,1,koo)-&
                              0.5*(virdzdx(i,1,1)*dfp+fdp*dxp(1)*dxp(1))*celcf%noxyz/boxcf%volu
                         sigc(1,2,koo)=sigc(1,2,koo)-&
                              0.5*(virdzdx(i,1,2)*dfp+fdp*dxp(1)*dxp(2))*celcf%noxyz/boxcf%volu
                         sigc(1,3,koo)=sigc(1,3,koo)-&
                              0.5*(virdzdx(i,1,3)*dfp+fdp*dxp(1)*dxp(3))*celcf%noxyz/boxcf%volu
                         sigc(2,1,koo)=sigc(2,1,koo)-&
                              0.5*(virdzdx(i,2,1)*dfp+fdp*dxp(2)*dxp(1))*celcf%noxyz/boxcf%volu
                         sigc(2,2,koo)=sigc(2,2,koo)-&
                              0.5*(virdzdx(i,2,2)*dfp+fdp*dxp(2)*dxp(2))*celcf%noxyz/boxcf%volu
                         sigc(2,3,koo)=sigc(2,3,koo)-&
                              0.5*(virdzdx(i,2,3)*dfp+fdp*dxp(2)*dxp(3))*celcf%noxyz/boxcf%volu
                         sigc(3,1,koo)=sigc(3,1,koo)-&
                              0.5*(virdzdx(i,3,1)*dfp+fdp*dxp(3)*dxp(1))*celcf%noxyz/boxcf%volu
                         sigc(3,2,koo)=sigc(3,2,koo)-&
                              0.5*(virdzdx(i,3,2)*dfp+fdp*dxp(3)*dxp(2))*celcf%noxyz/boxcf%volu
                         sigc(3,3,koo)=sigc(3,3,koo)-&
                              0.5*(virdzdx(i,3,3)*dfp+fdp*dxp(3)*dxp(3))*celcf%noxyz/boxcf%volu
                      end if
                   else
                      sigcalfo(1,1)=sigcalfo(1,1)+0.5*F1*dxp(1)/boxcf%volu
                      sigcalfo(1,2)=sigcalfo(1,2)+0.5*F1*dxp(2)/boxcf%volu
                      sigcalfo(1,3)=sigcalfo(1,3)+0.5*F1*dxp(3)/boxcf%volu
                      sigcalfo(2,1)=sigcalfo(2,1)+0.5*F2*dxp(1)/boxcf%volu
                      sigcalfo(2,2)=sigcalfo(2,2)+0.5*F2*dxp(2)/boxcf%volu
                      sigcalfo(2,3)=sigcalfo(2,3)+0.5*F2*dxp(3)/boxcf%volu
                      sigcalfo(3,1)=sigcalfo(3,1)+0.5*F3*dxp(1)/boxcf%volu
                      sigcalfo(3,2)=sigcalfo(3,2)+0.5*F3*dxp(2)/boxcf%volu
                      sigcalfo(3,3)=sigcalfo(3,3)+0.5*F3*dxp(3)/boxcf%volu
                      if (lcalcsigc.EQV..true.) then
                         sigc(1,1,koo)=sigc(1,1,koo)+0.5*F1*dxp(1)*celcf%noxyz/boxcf%volu
                         sigc(1,2,koo)=sigc(1,2,koo)+0.5*F1*dxp(2)*celcf%noxyz/boxcf%volu
                         sigc(1,3,koo)=sigc(1,3,koo)+0.5*F1*dxp(3)*celcf%noxyz/boxcf%volu
                         sigc(2,1,koo)=sigc(2,1,koo)+0.5*F2*dxp(1)*celcf%noxyz/boxcf%volu
                         sigc(2,2,koo)=sigc(2,2,koo)+0.5*F2*dxp(2)*celcf%noxyz/boxcf%volu
                         sigc(2,3,koo)=sigc(2,3,koo)+0.5*F2*dxp(3)*celcf%noxyz/boxcf%volu
                         sigc(3,1,koo)=sigc(3,1,koo)+0.5*F3*dxp(1)*celcf%noxyz/boxcf%volu
                         sigc(3,2,koo)=sigc(3,2,koo)+0.5*F3*dxp(2)*celcf%noxyz/boxcf%volu
                         sigc(3,3,koo)=sigc(3,3,koo)+0.5*F3*dxp(3)*celcf%noxyz/boxcf%volu
                      end if
                   endif
                endif
62              CONTINUE
61              CONTINUE
699             CONTINUE


                POTISTcalfo=potis1




                return
              end subroutine

            end module calfow_mod
