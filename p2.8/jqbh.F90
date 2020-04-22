module jqbh_mod
        USE tempinst_mod
        USE cryst_to_cart_mod
        USE gen_com_m, ONLY:at,bg,epcoud,epsil,erg2ev,erg2joule,im,it,ittherm,kthg,njqbh,ntr,&
             &rang,rulayer,tstep,nzl,zl,zls2

        implicit none
        contains
subroutine jqbh (xp,xpp,vp,ityp)

  USE T_kind_param_m, ONLY:  double
  USE var_pot, ONLY:cm
#ifdef PARA
  USE mod_para
#endif
  implicit none
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  integer  :: ityp(imm)
  real(double)  :: xp(3,imm)
  real(double)  :: xpp(3,imm)
  real(double)  :: vp(3,imm)
  !-----------------------------------------------
  !   L o c a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: nacou1,nacou2,i,indtr,lenfn2,ii,itr
  integer,save :: imesureT
  integer, save :: iunite
  real(double):: tcou1,tcou2,dec1,dec2,delt,ecou1,ecou2,alph1,alph2,xdec
  real(double), save ::tfcou1,tfcou2,tmcou1,tmcou2,csup,cinf,crul
  real(double),allocatable,save :: temptr(:)
  real(double),allocatable,save :: temptra(:)
  integer,allocatable,save::nattr(:)
#ifdef PARA 
  real(double) :: ecou1_tot,ecou2_tot
  integer::nacou1_tot,nacou2_tot
  real(double),allocatable,save :: temptra_tot(:)
  integer,allocatable,save::nattr_tot(:)

#endif
  character*15:: fnamtr
  character(len=2) :: extension
  logical :: loc(imm)
  !real(double):: dTtot,tempact,dTloc,tempinst,crulinv
  real(double):: dTtot,tempact,dTloc,crulinv
  ! tempinst does not need declaration becaUSE it is declared in tempinst_mod
  if(it.eq.1) then
     if(rang==0) write(6,*)'condutivité thermique méthode directe'

     if(rang==0) write(6,*)'METHODE n° ',njqbh
     if(rang==0) write(6,*)'epaisseur', epcoud
     if(rang==0) write(6,*)'flux au bord ', epsil, ' erg',epsil*erg2joule,' joules', epsil*erg2eV,' eV' 
     if ((njqbh==5).and.(rang==0)) write(6,*)'kthg ', kthg
     if (njqbh.lt.4) then
        if(rang==0) write(6,*) 'You d rather not...'
        stop
     end if
     open(file='DeltaE',unit=58)
     epcoud=epcoud*1.0d-8
     allocate (temptr(ntr))
     allocate (temptra(ntr))
     allocate (nattr(ntr))
#ifdef PARA
     allocate (temptra_tot(ntr))
     allocate (nattr_tot(ntr))
#endif
     temptr(:)=0. ; nattr(:)=0 ; temptra(:)=0
     if (rang==0) then
        do i=1,ntr
          write(extension,'(i2.2)') i

200        format(i2)
101        format(a1)
201        format(a2)
           fnamtr = 'tranche.'//extension
           ii=i+1000
           open(ii, file=fnamtr, form='formatted')
        end do
     end if

     if (njqbh==5) then
        if(rang==0) write(6,*)'tempinst',tempinst(vp,ityp)
        dTtot=epsil*(nzl(1)-2*rulayer)/(nzl(2)*nzl(3)*kthg*tstep)
        if(rang==0) write(6,*)'dTtot', dTtot
        tempact=tempinst(vp,ityp)
        temptra(:)=0.
        nattr(:)=0

        do itr=1,ntr
           dTloc=-2*(float(itr)/ntr-0.5)*dTtot+tempact


           loc(:)=.false.

           call cryst_to_cart (imm, xp, bg, -1) !cart vers cryst
           crul=rulayer/nzl(1)           
           crulinv=1-crul
!           write(6,*)'rulayer, nzl,crul',rulayer, nzl(1),crul
           do i=1,imd
              if((xp(1,i).lt.crul).or.(xp(1,i).gt.crulinv))cycle
              indtr=1+Int(ntr*(xp(1,i)-crul)/(1-2*crul))
              if (indtr==itr) then
                 loc(i)=.true.
                 nattr(indtr)=nattr(indtr)+1
                 temptra(indtr)=temptra(indtr)+ (vp(1,i)**2+vp(2,i)**2+vp(3,i)**2)*cm(ityp(i))/(3.*bk)
              end if
           end do

           temptra(itr) = temptra(itr)/float(nattr(itr))
           if(rang==0) write(6,*)'itr dTloc ',itr,dTloc,temptra(itr),nattr(itr)
           alph2=sqrt(dTloc/temptra(itr))
!           if(rang==0) write(6,*)'alph2',alph2
           call cryst_to_cart (imm, xp, at, 1)  !cryst vers cart

           do i = 1, imd
              if (loc(i))then
                 xpp(:,i) = xp(:,i)-(xp(:,i)-xpp(:,i))*alph2
                 vp(:,i)=vp(:,i)*alph2
              end if
           end do

        end do
     end if
!     tcou2=tempinst(vp,ityp)
!     if(rang==0) write(6,*)'tempinst',tcou2
  end if
  
  select case (njqbh)
  case(4,5)             ! imposition d'un flux en +/-Z 
     loc(:)=.false.
     !     write(6,*)'tempinst BCL',tempinst(vp,ityp)
     !write(6,*)'RRRRRRRRRRRRRR',xp(1,38603)
     nacou1 = 0 ; nacou2=0
     if (it==1)then 
        if(rang==0) write(6,*)'boite non périodique, flux de - à + (ZL/2 -rumax)'
        epcoud=epcoud+rulayer           
     end if

     ecou1 = 0.d0;      ecou2 = 0.d0
     nacou1=0; nacou2=0

     call cryst_to_cart (imm, xp, bg, -1) !cart vers cryst

     csup = 1-epcoud/nzl(1)
     cinf = epcoud/nzl(1)
     crul=rulayer/nzl(1)
!           write(6,*)'rulayer, nzl,crul',rulayer, nzl,crul
     if ((it==1).and.(rang==0)) write(6,*)'cinf,csup crul',cinf,csup,crul


#ifdef PARA 
     do i = 1, im
        if(free(i))then
           if (xp(1,i)<cinf) then
              nacou1 = nacou1+1
              ecou1 = ecou1+0.5*(vp(1,i)**2+vp(2,i)**2+vp(3,i)**2)*cm(ityp(&
                   i))
              loc(i)=.true.
           end if
        end if
     end do
!     write(6,*)'ecou1', rang,nacou1,ecou1
     call MPI_ALLREDUCE(ecou1,ecou1_tot,1,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
     ecou1=ecou1_tot
     call MPI_ALLREDUCE(nacou1,nacou1_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
     nacou1=nacou1_tot
!     write(6,*)'ecou1B', rang,nacou1,ecou1

#else

     do i = 1, imd
        if (xp(1,i)<cinf) then
           nacou1 = nacou1+1
           ecou1 = ecou1+0.5*(vp(1,i)**2+vp(2,i)**2+vp(3,i)**2)*cm(ityp(&
                i))
           loc(i)=.true.
        end if
     end do


#endif



     call cryst_to_cart (imm, xp, at, 1)  !cryst vers cart

     tcou1=ecou1*2./(3.*bk*nacou1)
     alph1=sqrt(1+epsil/ecou1)
!     if(rang==0) write(6,*)'nacou1 alph1 tcou1', nacou1,alph1,tcou1
     do i = 1, imd
        if (loc(i)) then
           xpp(:,i) = xp(:,i)-(xp(:,i)-xpp(:,i))*alph1
           vp(:,i)=vp(:,i)*alph1
        end if
     end do




     loc(:)=.false.
     call cryst_to_cart (imm, xp, bg, -1) !cart vers cryst

#ifdef PARA
     do i = 1, im
        if (free(i))then
           if (xp(1,i).gt.csup)then
              nacou2 = nacou2+1
              ecou2 = ecou2+(vp(1,i)**2+vp(2,i)**2+vp(3,i)**2)*cm(ityp(&
                   i))*0.5
              loc(i)=.true.
           end if
        end if
     end do
!     write(6,*)'ecou2', rang,nacou2,ecou2
     call MPI_ALLREDUCE(ecou2,ecou2_tot,1,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
     ecou2=ecou2_tot
     call MPI_ALLREDUCE(nacou2,nacou2_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
     nacou2=nacou2_tot
!     write(6,*)'ecou2B', rang,nacou2,ecou2
#else

     do i = 1, imd
        if (xp(1,i).gt.csup)then
           nacou2 = nacou2+1
           ecou2 = ecou2+(vp(1,i)**2+vp(2,i)**2+vp(3,i)**2)*cm(ityp(&
                i))*0.5
           loc(i)=.true.
        end if
     end do
#endif



     call cryst_to_cart (imm, xp, at, 1)  !cryst vers cart
     tcou2=ecou2*2./(3.*bk*nacou2)
     alph2=sqrt(1-epsil/ecou2)
!                      write(6,*)'nacou2 alph2 tcou2', nacou2,alph2,tcou2

     do i = 1, imd
        if (loc(i))then
           xpp(:,i) = xp(:,i)-(xp(:,i)-xpp(:,i))*alph2
           vp(:,i)=vp(:,i)*alph2
        end if
     end do



     call cryst_to_cart (imm, xp, bg, -1) !cart vers cryst

#ifdef PARA
     temptra(:)=0.
     nattr(:)=0
     do i=1,imd
	        if(free(i))then	
        indtr=1+Int(ntr*(xp(1,i)-crul)/(1-2*crul))
        !          write(6,*)i,indtr,xp(1,i), (xp(1,i)-crul)/(1-2*crul)
        nattr(indtr)=nattr(indtr)+1
        temptra(indtr)=temptra(indtr)+ (vp(1,i)**2+vp(2,i)**2+vp(3,i)**2)*cm(ityp(i))/(3.*bk*ittherm)
		endif
     end do
!     write(6,*)'temptra', rang,temptra,nattr
     call MPI_ALLREDUCE(temptra,temptra_tot,ntr,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
     temptra=temptra_tot
     call MPI_ALLREDUCE(nattr,nattr_tot,ntr,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
     nattr=nattr_tot
!     write(6,*)'temptra2', rang,temptra,nattr

#else
     temptra(:)=0.
     nattr(:)=0
     crulinv=1-crul
     do i=1,imd
        if((xp(1,i).lt.crul).or.(xp(1,i).gt.crulinv))cycle
        indtr=1+Int(ntr*(xp(1,i)-crul)/(1-2*crul))
!              write(6,*)i,indtr,xp(1,i), (xp(1,i)-crul)/(1-2*crul)
        nattr(indtr)=nattr(indtr)+1
        temptra(indtr)=temptra(indtr)+ (vp(1,i)**2+vp(2,i)**2+vp(3,i)**2)*cm(ityp(i))/(3.*bk*ittherm)
     end do
#endif

     temptr(:)=temptr(:)+temptra(:)/nattr(:)

     if (rang==0)then
     if (mod(it,ittherm).eq.0)then     
        do i=1,ntr
           ii=i+1000
           write(ii,*)it,temptr(i),nattr(i)
        end do
        temptr(:)=0
        write(58,'(I6,4D15.7,I11)')it,alph1,alph2,tcou1,tcou2,sum(nattr)
!         write(6,'(A,I6,4D15.7,I6)')'!!KTH!! ',it,alph1,alph2,tcou1,tcou2,sum(nattr)
     end if
  endif
     call cryst_to_cart (imm, xp, at, 1)  !cryst vers cart

     ! AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
     !FAUX FAUX FAUX
     !A REECRIRE

  case(3)    ! chauffage aux bord refroidissement au centre
     nacou1 = 0 ; nacou2=0
     ecou1 = 0.d0;      ecou2 = 0.d0
     do i = 1, im
        if (zls2(1)-abs(xp(1,i))<epcoud) then
           nacou1 = nacou1+1
           ecou1 = ecou1+0.5*(vp(1,i)**2+vp(2,i)**2+vp(3,i)**2)*cm(ityp(&
                i))
        end if
     end do
     tcou1=ecou1*2./(3.*bk*nacou1)
     alph1=sqrt(1+epsil/ecou1)
     !            write(6,*)'nacou1 alph1 tcou1', nacou1,alph1,tcou1
     do i = 1, imd
        if (zls2(1)-abs(xp(1,i))<epcoud) then
           xpp(:,i) = xp(:,i)-(xp(:,i)-xpp(:,i))*alph1
           vp(:,i)=vp(:,i)*alph1
        end if
     end do


     do i = 1, im
        if (abs(xp(1,i)).lt.epcoud)then
           nacou2 = nacou2+1
           ecou2 = ecou2+(vp(1,i)**2+vp(2,i)**2+vp(3,i)**2)*cm(ityp(&
                i))*0.5
        end if
     end do
     tcou2=ecou2*2./(3.*bk*nacou2)
     alph2=sqrt(1-epsil/ecou1)
     !            write(6,*)'nacou2 alph2 tcou2', nacou2,alph2,tcou2
     do i = 1, im
        if (abs(xp(1,i)).lt.epcoud) then
           xpp(:,i) = xp(:,i)-(xp(:,i)-xpp(:,i))*alph2
           vp(:,i)=vp(:,i)*alph2
        end if
     end do


     write(58,'(I6,4D15.7)')it,alph1,alph2,tcou1,tcou2

     temptra(:)=0.
     nattr(:)=0
     do i=1,im
        indtr=1+Int(ntr*(xp(1,i)+zls2(1))/zl(1))
        !            write(6,*)indtr
        nattr(indtr)=nattr(indtr)+1
        temptra(indtr)=temptra(indtr)+ (vp(1,i)**2+vp(2,i)**2+vp(3,i)**2)*cm(ityp(i))/(3.*bk*ittherm)
     end do
     temptr(:)=temptr(:)+temptra(:)/nattr(:)

     if (mod(it,ittherm).eq.0)then

        do i=1,ntr
           ii=i+1000
           write(ii,*)it,temptr(i)
        end do
        temptr(:)=0
     end if

     ! autres cas : imposition de l etmperature te mesure du flux ! methode abandonnée


  case(1)


     if(.not.(mod(it,ittherm).eq.0))return
     nacou1 = 0 ; nacou2=0
     tcou1 = 0.d0;      tcou2 = 0.d0
     do i = 1, im
        !         delt=zls2(1)-abs(xp(1,i))
        !         write(6,*)delt,epcoud
        if (zls2(1)-abs(xp(1,i))<epcoud) then
           nacou1 = nacou1+1
           tcou1 = tcou1+(vp(1,i)**2+vp(2,i)**2+vp(3,i)**2)*cm(ityp(&
                i))/(3.0*bk)
        end if
     end do
     tcou1 = tcou1/float(nacou1)
     dec1=1.5*bk*(tcou1-tfcou1)
     if(rang==0) write(6,*)'nacou1 tcou1 ', nacou1,tcou1
     do i = 1, im
        if (zls2(1)-abs(xp(1,i))<epcoud) &
             &         xpp(:,i) = xp(:,i)-(xp(:,i)-xpp(:,i))*sqrt(tfcou1/tcou1)
     end do


     do i = 1, im
        if (abs(xp(1,i)).lt.epcoud)then
           nacou2 = nacou2+1
           tcou2 = tcou2+(vp(1,i)**2+vp(2,i)**2+vp(3,i)**2)*cm(ityp(&
                i))/(3.0*bk)
        end if
     end do
     tcou2 = tcou2/float(nacou2)
     dec2=1.5*bk*(tcou2-tfcou2)
!     if(rang==0) write(6,*)'nacou2 tcou2 ', nacou2,tcou2
     do i = 1, im
        if (abs(xp(1,i)).lt.epcoud) &
             &          xpp(:,i) = xp(:,i)-(xp(:,i)-xpp(:,i))*sqrt(tfcou2/tcou2)
     end do


     write(58,'(I6,4D18.10)')it,dec1,dec2,tcou1,tcou2





  case(2)

     nacou1 = 0 ; nacou2=0
     tcou1 = 0.d0; tcou2 = 0.d0
     do i = 1, im
        if (zls2(1)-abs(xp(1,i))<epcoud) then
           nacou1 = nacou1+1
           tcou1 = tcou1+(vp(1,i)**2+vp(2,i)**2+vp(3,i)**2)*cm(ityp(&
                i))/(3.0*bk)
        end if
     end do
     tcou1 = tcou1/float(nacou1)
     tmcou1=tmcou1+tcou1/ittherm
     !               write(6,*)'tcou1 ',tcou1            

     do i = 1, im
        if (abs(xp(1,i)).lt.epcoud)then
           nacou2 = nacou2+1
           tcou2 = tcou2+(vp(1,i)**2+vp(2,i)**2+vp(3,i)**2)*cm(ityp(&
                i))/(3.0*bk)
        end if
     end do
     tcou2 = tcou2/float(nacou2)
     tmcou2=tmcou2+tcou2/ittherm
     !               write(6,*)'tcou2 ',tcou2            

     if (mod(it,ittherm).eq.0)then

        dec1=1.5*bk*(tmcou1-tfcou1)
        if(rang==0) write(6,*)'tcou1 tmcou1 nacou1 ',it, tcou1,tmcou1,nacou1
        do i = 1, im
           if (zls2(1)-abs(xp(1,i))<epcoud) &
                &         xpp(:,i) = xp(:,i)-(xp(:,i)-xpp(:,i))*sqrt(tfcou1/tmcou1)
        end do

        dec2=1.5*bk*(tcou2-tfcou2)
        if(rang==0) write(6,*)'tcou2 mtcou2 nacou2', it,tcou2,tmcou2,nacou2
        do i = 1, im
           if (abs(xp(1,i)).lt.epcoud) &
                &          xpp(:,i) = xp(:,i)-(xp(:,i)-xpp(:,i))*sqrt(tfcou2/tmcou2)
        end do



        write(58,'(I6,4D15.7)')it,dec1,dec2,tmcou1,tmcou2
        tmcou1=0. ; tmcou2=0.

     end if




  end select

end subroutine jqbh

end module
