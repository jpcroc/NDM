module controle_mod
        use endrun_mod 
        use dynalloccell
        use tempinst_mod
        use jqbh_mod
        use caltabt_mod
        use desinteg_insert_mod
        use period_mod
        use caltabi_mod
        use heat_mod
        use creadp_mod
        implicit none
        contains
! ***********************************************************
!           sous-programme controle.f
! ***********************************************************

subroutine controle
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use var_pot
  use tab_imm_m
  use suivinonpbc
#ifdef PARA
  use mod_para
#endif
  use defcdp, ONLY :itecdp
  implicit none
  !-----------------------------------------------
  !   G l o b a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: nacou, i, ic, iti
  real(double) :: vv, a1, a2, a3, c1, c2, c3
  real(double), dimension(1,3) :: g1,aux
  real(double) :: ltc, ctime, tdev, tcool, epc1, epc2, epc3,masstot,massa,tclt
  real(double), dimension(1,3) :: xtr, cv
  real(double) :: fpmax,fpn,forctot,formax
  real(double) :: xpnp(3,imm)
  real(double) :: potistmean,potistdif
  real(double),save :: potist1000
  real, pointer,save :: potiststock(:)
#ifdef PARA
  real(double) :: tcou_glob
  integer      :: nacou_glob
  real(double) :: fpmax_glob
  real(double) :: forctot_glob
  real(double) :: formax_glob
#endif
  save ltc
  !-----------------------------------------------
  !
  !
  !

!    write(6,*) 'entree controle',it,im


  ! last iteration ?


  !dewrite(*,*) 'IN CONTROLE: dmtype,  tfcou, idorcou', dmtype, tfcou, ibordcou
  !stop


  if (it>=itmax) then
     if (rang==0) write (6, *) '*******Derniere iteration **** '
     call endrun
     call DeallocateAll

     call arret_ndm

  endif
  
  if (timel>=timemax) then
     if (rang==0) write (6, *) '*******max time reached **** ',timel,timemax
     call endrun
     call DeallocateAll

     call arret_ndm

  endif
  if(lEev) then
     unitE=erg2eV
     cunitE='  eV'
  else
     unitE=1.0
     cunitE=' erg'
  end if
  if(lPkbar) then
     unitP=1.0d-9
     cunitP='kbar'
  else
     unitP=1.0
     cunitP='d/cm2'
  endif

  ! temperature is down enough ?
  if (itetemp>0) then
     if (mod(it,itetemp)==0) then
        if (temp<=tempstop) then
	if (rang==0)  write (6, *) 'temperature < tempstop '
           call endrun
           call DeallocateAll

	   call arret_ndm

        endif
        if (tempstopcel.gt.0) then
           if (maxtcel<=tempstopcel) then
              write (6, *) 'temperature dans toutes les cels < tempstopcel '
              call endrun
              call DeallocateAll
              
              call arret_ndm
           end if
        endif
     endif
  endif



  if (lcdp) then 
     !     write(6,*)'it',it
     if(itecdp.gt.0) then
        if (mod(it,itecdp).ne.0) goto 666


        call creadp (xp, xpp, ityp,vp)
        call caltabt
        if (ltabvois) call caltabi
        if (lperiod) then 
           call period
        else 
           write(*,*) 'WARNING .... Not implemented for lperiod  FALSE nad lcdp TRUE'
           write(*,*) 'FIX THAT! Until there the program will stop'
           stop      
        end if
     end if
  end if
666 continue


  ! cell dispatching
  if (dmtype.ne.4) then
     if (itab/=0) then
        if (mod(it,itab)==0) then
           call caltabt
        endif
     endif
  endif !dmtype

  if (ltabvois.and.mod(it,itetabvois)==0) call caltabi 

  !write(6,*)it


  if (dmtype==1) then

     if (lperiod) then
        xpnp(:,:)=xp(:,:)
     else
        call notperiod(xp,xpnp)
     end if

     ! Scaling temperature if intolerable ?
     if (itetemp/=0) then
        if (mod(it,itetemp)==0) then
           if (ttol>0.) then
              tdev = dabs(temp-tfroi)/tfroi*1.d2
              if (tdev>=ttol) then
	         if (rang==0) then
                    write (6, *) 'pas :',it,'refroidir de T=', temp,&
                         ' a T=', tfroi
		 endif
                 vv = sqrt(tfroi/temp)
                 xpp(:,:im) = xp(:,:im)-(xp(:,:im)-xpp(:,:im))*vv
                 if (lperiod) call period              !Conditions periodiques
              endif  !end tdev>=toll
           endif !end toll>0
        endif !end mod(it,itetemp)==0
     endif ! end itetemp/=0

     ! change in time step ? itetimestep >0
     if (itetimestep>0) then
        if (mod(it,itetimestep)==0) call deftimestep 
     endif


     ! controle temperature en bord de boite ?  tcou >0

     if (tfcou > 0.0) then
        if (itetemp/=0) then
           if (mod(it,itetemp)==0) then

              !                                                
              call cryst_to_cart (imm, xp, bg, -1)          !cart vers cryst
              call cryst_to_cart (imm, xpnp, bg, -1)        !cart vers cryst
              epc1 = epcou/zl(1)
              epc2 = epcou/zl(2)
              epc3 = epcou/zl(3)
              nacou = 0
              tcou = 0.d0



              select case(ibordcou)
              case(0)
                 do i = 1, im
                    if (.not.(xpnp(1,i)<epc1.or.xpnp(1,i)>1.0-epc1.or.xpnp(2,i)<&
                         epc2.or.xpnp(2,i)>1.0-epc2.or.xpnp(3,i)<epc3.or.xpnp(3,i)>&
                         1.0-epc3)) cycle
                    nacou = nacou+1
                    tcou = tcou+(vp(1,i)**2+vp(2,i)**2+vp(3,i)**2)*cm(ityp(&
                         i))/(3.0*bk)
                 end do
                 tcou = tcou/float(nacou)

                 do i = 1, im
                    if (.not.(xpnp(1,i)<epc1.or.xpnp(1,i)>1.0-epc1.or.xpnp(2,i)<&
                         epc2.or.xpnp(2,i)>1.0-epc2.or.xpnp(3,i)<epc3.or.xpnp(3,i)>&
                         1.0-epc3)) cycle
                    xtr(1,:) = xp(:,i)

                    call cryst_to_cart (1, xtr, at, 1)
                    xpp(:,i) = xtr(1,:)-(xtr(1,:)-xpp(:,i))*sqrt(tfcou/&
                         tcou)
                 end do
              case(3)
                 do i = 1, im
                    if (.not.(xpnp(3,i)<epc3.or.xpnp(3,i)>&
                         1.0-epc3)) cycle
                    nacou = nacou+1
                    tcou = tcou+(vp(1,i)**2+vp(2,i)**2+vp(3,i)**2)*cm(ityp(&
                         i))/(3.0*bk)
                 end do

                 tcou = tcou/float(nacou)

                 do i = 1, im
                    if (.not.(xpnp(3,i)<epc3.or.xpnp(3,i)>&
                         1.0-epc3)) cycle
                    xtr(1,:) = xp(:,i)

                    call cryst_to_cart (1, xtr, at, 1)
                    xpp(:,i) = xtr(1,:)-(xtr(1,:)-xpp(:,i))*sqrt(tfcou/&
                         tcou)
                 end do
              end select
              !     !cryst vers cart
              call cryst_to_cart (imm, xp, at, 1)

           endif                             ! itetemp
        endif                                ! mod(it,itetemp)
     endif                                   !controle couche externe




     ! test du flux de chaleur a decommenter avec precaution






  endif                                      !dmtype=1 end here




  if (ljqbh) then
     if (lperiod) then
        call jqbh(xp,xpp,vp,ityp)
     else 
        write(*,*) 'No implementation for ljqbh .true. and lperiod .false.'
        write(*,*) 'Stop in controle'
        stop        
     end if
  end if

  ! change in time step ? itetimestep >0
  if (itetimestep>0) then
     if (mod(it,itetimestep)==0) call deftimestep 
  endif


  if (dmtype==4) then

     !debug    write(*,*) 'IN CONTROLE: dmtype,  tfcou, idorcou, landerscou', dmtype, tfcou, ibordcou, landerscou

     if (lperiod) then
        xpnp(:,:)=xp(:,:)
     else       
        call notperiod(xp,xpnp)
     end if


     ! Temperature constante a  la Andersen
     if (LTandersen) then
        do i=1,im
           massa=cm(ityp(i))
        end do
     end if

     if (mod(it,itetemp)==0) then
        if (ttol>0.) then
           tdev = dabs(temp-tfroi)/tfroi*1.d2
           if (tdev>=ttol) then
              write (6, *) 'refroidir de T=', temp, ' a T=', tfroi
              vv = sqrt(tfroi/temp)
              vp(:,:im) = vp(:,:im)*vv
           endif
        endif
     endif



     ! controle temperature en bord de boite ?  tcou >0 AVEC L'ALGORITHME DE RESCALING
     if ((tfcou>0.0).and.(.not.landerscou)) then
        if (itetemp/=0) then
           if (mod(it,itetemp)==0) then

              !  
              call cryst_to_cart (imm, xp, bg, -1) !cart vers cryst
              call cryst_to_cart (imm, xpnp, bg, -1) !cart vers cryst
              epc1 = epcou/zl(1)
              epc2 = epcou/zl(2)
              epc3 = epcou/zl(3)
              nacou = 0
              tcou = 0.d0
              !      write (6,*) 'it, epc1a3',it,epc1,epc2,epc3
              do i = 1, im
                 if (.not.(xpnp(1,i)<epc1.or.xpnp(1,i)>1.0-epc1.or.xpnp(2,i)<&
                      epc2.or.xpnp(2,i)>1.0-epc2.or.xpnp(3,i)<epc3.or.xpnp(3,i)>&
                      1.0-epc3)) cycle
                 nacou = nacou+1
                 !      write(6,*) 'xp(1a3,i)',i,xp(1,i),xp(2,i),xp(3,i)
                 tcou = tcou+(vp(1,i)**2+vp(2,i)**2+vp(3,i)**2)*cm(ityp(&
                      i))/(3.0*bk)
              end do
#ifdef PARA
              call MPI_ALLREDUCE(tcou,tcou_glob,1,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
              tcou = tcou_glob
              call MPI_ALLREDUCE(nacou,nacou_glob,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
              nacou = nacou_glob
#endif

              !     write(6,*) 'it, nacou, tcou=',it,nacou,tcou
              tcou = tcou/float(nacou)
              !     write(6,*) 'tcou apres division=',tcou

              do i = 1, im
                 if (.not.(xpnp(1,i)<epc1.or.xpnp(1,i)>1.0-epc1.or.xpnp(2,i)<&
                      epc2.or.xpnp(2,i)>1.0-epc2.or.xpnp(3,i)<epc3.or.xpnp(3,i)>&
                      1.0-epc3)) cycle

                 vp(:,i) = vp(:,i)*sqrt(tfcou/tcou)
                 !       write(6,*) 'xpp(ic,i)apres',i,ic,xpp(ic,i)
              end do

              !                                                
              call cryst_to_cart (imm, xp, at, 1)  !cryst vers cart



           endif                             ! itetemp
        endif                                ! mod(it,itetemp)
     endif                                   !controle couche externe 


     ! controle temperature en bord de boite ?  tcou >0 AVEC L'ALGORITHME DE ANDERSEN

     if ((tfcou>0.0).and.(landerscou)) then
        call cryst_to_cart (imm, xp, bg, -1)  !cart vers cryst
        call cryst_to_cart (imm, xpnp, bg, -1)  !cart vers cryst
	epc1 = epcou/zl(1)
        epc2 = epcou/zl(2)
        epc3 = epcou/zl(3)
        nacou = 0
        tcou = 0.d0
        do i = 1, im
           if (.not.(xpnp(1,i)<epc1.or.xpnp(1,i)>1.0-epc1.or.xpnp(2,i)<&
                epc2.or.xpnp(2,i)>1.0-epc2.or.xpnp(3,i)<epc3.or.xpnp(3,i)>&
                1.0-epc3)) cycle
           massa=cm(ityp(i))
           call andersenth(vp(:,i),massa,Text,nuandersen,tstep,bk)
        end do
        call cryst_to_cart (imm, xp, at, 1)    !cryst vers cart

     end if
  endif      ! end of dmtype=4                                    

  ! cooling of the atoms ? tcooling >0
  if (tcooling>0) then
     if (itetemp/=0) then
        if (mod(it,itetemp)==0) then
           ctime = timel-lastcool
           if (lastcool==0) then
              tcool = temp-ctime*tcooling
           else
              tcool = ltc-ctime*tcooling
           endif
           tclt=tcool-tcooling*tstep*itetemp
           if (tclt.le.0) &
                & call endrun
           if (tcool<temp) then
              ltc = tcool
              lastcool = timel
              write (6, *) 'refroidir de T=', temp, ' a T=', tcool
              vv = sqrt(tcool/temp)
              if (dmtype==4) then
                 vp(:,:im) = vp(:,:im)*vv
              else
                 xpp(:,:im) = xp(:,:im)-(xp(:,:im)-xpp(:,:im))*vv
              endif
           endif
        endif
     endif    ! end of itetemp /= 0
  endif       !end of tcooling > 0






  ! *** correction de la derive ***
  if (itederive>0) then
     if (.not.lperiod) then
        write(*,*) '----------------WARNING-----------------------------'
        write(*,*) 'there is no implemantation for itederive > 0 and lperiod=.false.'
        write(*,*) 'However, you are free to implement that'
        write(*,*) 'After that please call 2 9185 and/or ask for Jean-Paul in Bat. 520'
     end if
     if (mod(it,itederive)==0) then

        aux(1,:) = 0.0                        ! initialisation du sommateur de delta X


        call cryst_to_cart (imm, xp, bg, -1)    !cart vers cristal
        call cryst_to_cart (imm, ax, bg, -1) !cart vers cristal

        g1(1,:) = 0.0
        masstot=0.
        do i=1,im
           if (i==iko) cycle
           do ic =1,3
              g1(1,ic)=g1(1,ic)+xp(ic,i)*cm(ityp(i))
           end do
           masstot=masstot+cm(ityp(i))
        end do
        if (iko>0) then
           g1(1,:) = g1(1,:)/(im*masstot)
        else
           g1(1,:) = g1(1,:)/((im-1)*masstot)
        end if

        !centre de masse des ax
        aux(1,:) = 0.0
        do i = 1, im
           if (i==iko) cycle
           do ic=1,3
              aux(1,ic)=aux(1,ic)+ax(ic,i)*cm(ityp(i))
           end do
        end do

        if (iko>0) then
           aux(1,:) = aux(1,:)/(im*masstot)
        else
           aux(1,:) = aux(1,:)/((im-1)*masstot)
        end if

        write(6,*)'g1 avant ', g1
        write(6,*)'auxavant ', aux
        call cryst_to_cart (imm, ax, at, 1)
        call cryst_to_cart (imm, xp, at, 1)
        call cryst_to_cart (1, g1, at, 1)
        call cryst_to_cart (1, aux, at, 1)
        write(6,*)'g1 apres ', g1
        write(6,*)'auxapres ', aux
        write(6,*)
        !correction de la derive
        do i=1,im
           xp(:,i) = xp(:,i)+aux(1,:)-g1(1,:)
           xpp(:,i) = xpp(:,i)+aux(1,:)-g1(1,:)
        end do


        !                                                !Conditions periodiques
        if (lperiod) call period

     endif   ! end of mod(it,itederive)==0
  endif      ! end of itederive > 0




  if (deltaestop.ne.0) then

     if (it==1) then

        allocate(potiststock(0:nbmoye-1))
        potiststock(:)=0.
        potistmean=potist
        potiststock(it)=potist
     elseif (it<nbmoye) then

        potiststock(it)=potist
        potistmean=SUM(potiststock)/it
     else
        potiststock(mod(it,nbmoye))=potist
        potistmean=SUM(potiststock)/nbmoye
     end if

     potistdif=0


     if (it==1000)  potist1000=potistmean
     if(it>1000) then
        potistdif=(potistmean-potist1000)*erg2ev
        if(potistdif.le.deltaestop)then
           write(6,*)'delta e***, stop',it, timel, potistdif
           stop
        end if
     end if
     write(823,*)potist*erg2ev,potistmean*erg2ev,potistdif
  end if   !end of deltaestop.ne.0
!  write(6,* )'dmtype',dmtype


  select case (dmtype)

  case(2,10,8)
     if ((lprtrp.EQV..false.).and.((dmtype==10).or.(dmtype==8))) goto 123
     if ((fpstop>0.0).AND.(it.GE.1)) then 
        IF (lFrozen) THEN
           !fpmax=sqrt( MAXVAL( Sum(fp(1:3,1:im)**2,1), Free(1:im) ) )
           fpmax = MaxVal( Abs(fp(:,1:im)), .NOT.Frozen(:,1:im) )
        ELSE
           !fpmax=sqrt( MAXVAL( Sum(fp(1:3,1:im)**2,1) ) )
           fpmax = MaxVal( Abs(fp(:,1:im)) )
        END IF
#ifdef PARA
        call MPI_ALLREDUCE(fpmax,fpmax_glob,1,NDM_MPI_REAL_DOUBLE,MPI_MAX,MPI_COMM_WORLD,ierr)
        fpmax=fpmax_glob
#endif

        fpn=fpmax*erg2eV/angst
        !if (rang==0)     write(6,*)
        if (rang==0)      write(6,'("TR: force max, energy",i6,3E20.10)') it,fpn, potist*erg2eV
        if ( rang==0)     write (6, *) 'energie ',potist*erg2eV
        !if (rang==0)     write(6,'(a,2g20.12)')'force max cgs  ev/Ang ',fpmax, fpn
        if((rang==0).and.(sigstop.ge.0))write(6,*)'sigma max kbar', 1d-9*maxval(abs(sigtot))
        if (fpn.le.fpstop)then
           if (sigstop.ge.0) then
              if(maxval(abs(sigtot)).le.sigstop/1d-9) call endrun
           else
              call endrun
           end if

           if (rang==0)      write(6,*)
        end if
     end if

     if ((fsumstop>0.0).AND.(it.GE.1)) then 
        IF (lFrozen) THEN
           !fpmax=sqrt( Sum( SUM(fp(1:3,1:im)**2,1), Free(1:im) ) )
           fpmax=sqrt( SUM( fp(:,1:im)**2, .NOT.Frozen(:,1:im) ) )
        ELSE
           fpmax=sqrt( SUM(fp(:,1:im)**2) )
        END IF
#ifdef PARA
        fpmax=fpmax**2
        call MPI_ALLREDUCE(fpmax,fpmax_glob,1,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
        fpmax=sqrt(fpmax_glob)
#endif


        fpn=fpmax*erg2eV/angst
        if (rang==0)      write(6,*)
        if (rang==0)      write(6,*)'  sqrt ( sum_f F_i^2 ):  cgs  ev/Ang ',fpmax, fpn
        if((rang==0).and.(sigstop.ge.0))write(6,*)'sigma max kbar',1d-9* maxval(abs(sigtot))
        if (fpn.le.fsumstop) then
           if (sigstop.ge.0) then
              if(maxval(abs(sigtot)).le.sigstop/1d-9) call endrun
           else
              call endrun
           end if
           if (rang==0)      write(6,*)
        end if
     end if
!     if (rang==0) write (6, '(I10,A,D21.12,A)') it,  '*Epot = ', potist*unitE, cunitE
123  continue

  case(3,30) ! Gradient conjugue sur coordonnee cartesiennes (3) ou reduites (30)


     if (it==1) then
        if (lEev.EQV..true.) then 
           if (rang==0) write(6,*)'Resultats en eV, Ang'
        else
           if (rang==0) write(6,*)'Resultats en cgs'
        end if
        if (rang==0)      write(*,'(70("="))')
        if (rang==0)      write(*,'("CG:     ","iter",10(" "),"epsi",14(" "),"Fmax",14(" "), "Energy")')
        if (rang==0)      write(*,'(70("="))')
     end if


     !debug     write(*,*) 'DEBUG ALL IT IN CONTROLE',it

     IF (it.GE.1) THEN
        IF (lFrozen) THEN
           !forctot=sqrt( Sum( SUM(fp(1:3,1:im)**2,1), Free(1:im) ) )
           !formax=sqrt( MAXVAL( Sum(fp(1:3,1:im)**2,1), Free(1:im) ) )
           forctot = sqrt( SUM( fp(:,1:im)**2, .NOT.Frozen(:,1:im) ) )
           formax = MaxVal( Abs(fp(:,1:im)), .NOT.Frozen(:,1:im) ) 
        ELSE
           forctot=sqrt( SUM(fp(1:3,1:im)**2) )
           !formax=sqrt( MAXVAL( Sum(fp(1:3,1:im)**2,1) ) )
           formax = MaxVal( Abs(fp(:,1:im)) )
        END IF
#ifdef PARA
        call MPI_ALLREDUCE(formax,fpmax_glob,1,NDM_MPI_REAL_DOUBLE,MPI_MAX,MPI_COMM_WORLD,ierr)
        formax=fpmax_glob
        forctot=forctot**2
        call MPI_ALLREDUCE(forctot,fpmax_glob,1,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
        forctot=sqrt(fpmax_glob)
#endif

        if (lEev.EQV..true.) then 
           forctot = forctot*erg2eV/angst
           formax  = formax*erg2eV/angst
           if (rang==0) write(*,'("GC: ",i6,3E20.10)') it,forctot, formax, potist*erg2eV
           if (fpstop>0) then   
              if (formax.le.fpstop) then
                 if (rang==0) write(6,*)'force par atome  max  ev/Ang ', formax
                 if (rang==0) write (6, *) 'energie ', potist*erg2eV
                 if (it.le.1) xp(:,:)=ax(:,:)
                 call endrun
              end if
           end if
           if (fsumstop>0) then   
              if (forctot.le.fsumstop) then
                 if (rang==0) write(6,*)'  sqrt ( sum_f F_i^2 ):   ev/Ang ', forctot
                 if (rang==0) write(6, *) 'energie ', potist*erg2eV
                 if (it.le.1) xp(:,:)=ax(:,:)
                 call endrun
              end if
           end if

        else
           if (rang==0) write(*,'("GC: ",i6,3E20.10)') it,forctot, formax, potist
           if (fpstop>0) then   
              if (formax.le.fpstop) then
                 if (rang==0) write(6,*)'force par atome  max cgs ',formax
                 if (rang==0) write (6, *) 'energie ', potist
                 if (it.le.1) xp(:,:)=ax(:,:)
                 call endrun

              end if
           end if

           if (fsumstop>0) then   
              if (forctot.le.fsumstop) then
                 if (rang==0) write(6,*)'  sqrt ( sum_f F_i^2 ):  cgs  ', forctot
                 if (rang==0) write (6, *) 'energie ', potist
                 if (it.le.1) xp(:,:)=ax(:,:)
                 call endrun
              end if
           end if

        end if
     end if ! it .ge.1

  case default
  end select



  if ((ldesinteg.EQV..true.).and.(itdes==nstepdes))call desinteg_insert
  if ((iteheat.gt.0).and.(mod(it,iteheat)==0))call heat




!
  return
end subroutine controle
end module
