module analyse_mod
  USE temp_com,only:h0, volu,zls2,nzl, normat, nvois, celsize,eatomtotm,patcel,patcelmax,tempc,tcp,pmc,celpp,&
       &pm1,sigat,pmc,tempc,sigc,celpm1,tm1,ltabvois ! A EFFACER
  USE Mat_utils_mod
  USE spebc_fin_mod,only: spebc_fin
  USE adf_mod,only: adf
  USE calctemp_mod,only: calctemp
  USE calcdepla_mod,only: calcdepla
  USE calcdepla2_mod,only: calcdepla2
  USE calccoordo_mod,only: calccoordo
  USE calcdigr_mod,only: calcdigr
  USE calcangle_mod,only: calcangle
  USE bondval_mod,only: bondval
  USE rasmol_mod,only: rasmol
  USE rdf_mod,only: rdf
  USE prtplz_mod,only: prtplz


  use var_pot, only: iewald,l3c,npotmax,potisglue,potisrep,lpotentiel
  use gen_com_m, only:bk,cunite,deltaespr,deltaf,ecellpr,espr,flag_fin,fnose,iteanapos,iteangle,itebdv,&
       &itecfg,itecoordo,itedepla,itefcc,iterasmol,iterdf,itesigma,itetemp,itetemp2,kcell,kine,kinemean,knose,&
       &lambdades,leev,leparat,linstantfda,lpr,lprteattotm,lsigatcel,lthoover,ltnose,ltpcel,lucell,&
       &nfda,parallele,pist,pmean,potcp,potis1,potis2,potis3,potist,potistersoff,potiszbl,sigatcel,&
       &tcou,temp,tempep,tfcou,tmean,ucell,unite,unose,zhoover,sig,sigkine,lprtcel,lprtcel,natchk,natchk,tpseuils,&
       &sigtot,unitP,tdepla2,nrdf,lprtsigat,lprteat,lpkbar,linstantrdf,&
       &ldesinteg,itmax,cunitp,erg2ev,iteplz,itespebcout

  USE cellconfig,only:cell_config, ndm2cellconfig, cellconfig2ndm,caltabtC
  USE atomconfig,only:atom_config,atom_config_d,atom_config_e, ndm2config, config2ndm
  use boxconfig,only: box_config,ndm2boxconfig,boxconfig2ndm
  implicit none
contains
  ! ************************************************
  !         Sous-programme analyse.f
  ! ************************************************

  subroutine analyse
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double

    use tab_imm_m

    USE fcc_module
    USE cfg_module
    USE posana
    use elec_cell, only : Eelec,Teavg,Tecmax,ietm,eleccellmol
    use eloss, only : ibrake, elosselectot1, elosselectot
    implicit none
    !-----------------------------------------------
    !   G l o b a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s        
    !-----------------------------------------------
    integer :: i, iti, ic, ko,kx,ky,kz
    real(double), dimension(ntyp) :: temptyp


    type(atom_config_e)::atdml,attyp
    type(cell_config):: celndm,celtyp
    type(box_config)::boxndm

    real(double) :: ppot, pkin
    real(double) :: a1, a2, a3, b1, b2, b3, c1, c2, c3
    real(double) :: fteta, tbc, tca, tab, amod, bmod, cmod,kinetyp
    real(double) ::  alat
    real(double), external :: tempinst
    real(double), save :: volumean,amodmean,bmodmean,cmodmean,tcamean,tabmean,tbcmean
    real(double), dimension(3,3) :: transformation, strain, rotation, invh0

    INTEGER, dimension(:), allocatable :: fcc_nVoisins, fcc_cluster
    INTEGER, dimension(:,:), allocatable :: aux_int
    REAL(kind(0.d0)), dimension(:,:), allocatable :: aux_real
    CHARACTER(len=20), dimension(:), allocatable :: aux_title
    CHARACTER(len=100) :: out_file
    integer,save::ncalceattotm=0, nposmoy=0
    integer::ipot, nAux_real, n

    !  real(double)::celpP,celpp2,Tcp,Tcp2  
    real(double)::ptest
    integer::luvisuc=888,lenfn2,nprt
    character :: extension*9
    real(double)::xb(3),minp,maxp,mint,maxt
    real(double)::minpP,maxpP,mintP,maxtP
    real(double),save::Cminp,Cmaxp,Cmint,Cmaxt
    real(double),save::CminpP,CmaxpP,CmintP,CmaxtP
    real(double),save::CminpP2,CmaxpP2,CmintP2,CmaxtP2
    real(double),save::timelm1=0
    integer::koo
    !-----------------------------------------------
    !
    !
    !      logical:: lEev=.false., lPkbar=.false.
    !        write(6,*)'analyse',itetemp,it

    call ndm2cellconfig(celndm,noxyz,nox,noy,noz,natperc,nato,ncel,atincel,deltadist,celsize)
    call ndm2config(atdml,im,imm,xp,fp,ityp,ielat,num_at_glob=num_at_glob,ltabvois=ltabvois,&
         &iwmax=iwmax,indi=indi,nvois=nvois,vp=vp,xpp=xpp)
    call ndm2boxconfig(at,bg,zl,zls2,nzl,volu,normat,boxndm)

    if(it==1) then
       Cminp=100000
       Cmaxp=-100000
       Cmint=100000
       Cmaxt=-100000

       CminpP=100000
       CmaxpP=-100000
       CmintP=100000
       CmaxtP=-100000

       CminpP2=100000
       CmaxpP2=-100000
       CmintP2=100000
       CmaxtP2=-100000
    end if
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

    ! calcul de la temperature
    ! calcul de la temperature
    !  if (itetemp2>0) then
    !     if (mod(it,itetemp2)==0) then
    !        temp2=tempinst(vp,ityp)
    !        write(112,'(I10,G15.6,F12.2)')it,timel,temp2
    !     end if
    !  end if
 
    if (itetemp>0) then
       if (mod(it,itetemp)==0) then
          call calctemp (temp,kine,atdml,celndm)
          
          do iti=1,ntyp
             if (na(iti)==0) cycle
             atdml%lgul=.false.
             celtyp=celndm
             where(atdml%ityp==iti)
                atdml%lgul=.true.
             end where
             call atdml%fab(attyp)
             !write(6,*)'nbat',count(atdml%ityp==iti),na(iti),attyp%im,attyp%ityp
             call caltabtC(celtyp,attyp,lperiod,boxndm)
             call calctemp(temptyp(iti),kinetyp,attyp,celtyp)
             !temptyp=0
             call celtyp%dealloc ; call attyp%dealloc
          !        call calctemp (temptyp)
          end do
          ! MPI
          !remarque 1erg = 6.24d11 eV
          if (mod(it,itetemp2)==0) then
             if (rang==0) then

                if (ibrake.GT.0)  write(6,*)'electronic losses ', elosselectot, elosselectot1

                write (6, *)
                write (6, *)
                write (6, '(A,I7,A,G15.7)') '<<<<<<<<  ITERATION =', it, &
                     '  time = ', timel
                write (6, *)
                write (6, *)
                write (6, *) '--------- Temperatures et Pression--------'
                write (6, *) '----------valeurs instantanees------------'


                write (6, '(I10,G10.3,A,G21.12,A)') it, timel, '*Epot = ', potist*unitE, cunitE
                do ipot=1,npotmax
                   if (lpotentiel(ipot).eqv..true.) then
                      select case (ipot)
                      case (1:9)
                         write(6,'(A,G21.12,A)')'    *energie paire 2 corps = ',(potis1+potis2)*unitE, cunitE
                         if(l3c) write(6,'(A,G21.12,A)')'    *energie pot 3 corps = ',potcp*unitE, cunitE
                         if(iewald.gT.0) write (6, '(A,G21.12,A)') '    *energie pot coul recip = ', potis3*unitE, cunitE
                      case(10:12)
                         write(6,'(A,G21.12,A)')'    *energie PAIRE EAM = ',potisrep*unitE, cunitE
                         write(6,'(A,G21.12,A)')'    *energie GLUE  = ',potisglue*unitE, cunitE
                      case(13)
                         write(6,'(A,G21.12,A)')'    *energie Tersoff = ',potisTersoff*unitE, cunitE
                         If (potisZBL.ne.0)write(6,'(A,G21.12,A)')'    *energie ZBL = ',potisZBL*unitE, cunitE
                      end select
                   end if
                end do
                write (6,'(I10,G10.3,A,G21.12,A,a,f0.3,a)') it,timel,'*Ec = ',kine*unitE, cunitE, &
                     '  (', temp, ' K)'

                write (6,'(I10,G10.3,A,G21.12,A)') it,timel,'*Etot = ',(kine+potist)*unitE, cunitE
                If (l2T) then

                   write (6,'(I10,G10.3,A,G21.12,A)') it,timel,'*Eelec = ',Eelec*unitE, cunitE                 
                   write (6,'(I10,G10.3,A,G21.12,A)') it,timel,'*IE_Et = ',&
                        &(kine+potist+Eelec)*unitE, cunitE                 
                   write (6,'(I10,G10.3,A,G21.12,A)') it,timel,'*Telec = ',Teavg
                   write (6,'(I10,G10.3,A,G21.12,A)') it,timel,'*TempEP = ',TempEP
                   write (6,'(I10,G10.3,A,G21.12,3I5)') it,timel,'*maxTe = ',Tecmax,ietm(:)

                end If

                write (6, *)

                !	IF(it==1) Open(unit=774, file='energie_it.dat', status='unknown', action='write')
                !	IF(modulo(it,10)==0) WRITE(774,'(i6, f)') it, (kine+potist)*unitE  ! controle .. !*!

                !*! ---------------------------------- APPEL DE LA ROUTINE D'ECRITURE -------------
                ! ------------------------------------ DES FICHIERS DE SORTIE, DANS LE ------------
                ! ------------------------------------ CAS DES CL CONTROLEES EN CONTRAINTE --------
                IF (ibound==2 .OR. ibound==3) THEN
                   IF(flag_fin.EQV..true.) THEN
                      Call spebc_fin (.true.)
                   ELSE IF(itespebcout > 0. .AND. (mod(it,itespebcout)==0 .OR. it==1)) THEN
                      Call spebc_fin (.false.)
                   END IF
                END IF
                !*!


                if (lEparat) then
                   write (6, '(I10,G10.3,A,G21.12,A)') it, timel, '*Epot/at = ', potist*unitE/im_glob, cunitE
                   write(6,*)
                end if

                write (6, '(I10,G10.3,A,f0.3)') it,timel, &
                     '*Temp instantanee = ',temp

                if (tfcou>0.0) write (6, '(A,G15.4)') '*temperature externe = ', tcou

                IF (lpr) THEN
                   write(6,*) 'NPT With Parrinello-Rahman'
                   write(6,'(I10,G10.3,A,G21.12,A,a,f0.3,a)') it,timel,'*Kcell = ',Kcell*unitE,cunitE, &
                        '  (', 2.d0*Kcell/(9.d0*bk), ' K)'
                   IF (lUcell) THEN
                      write(6,'(I10,G10.3,A,G21.12,A)') it,timel,'*Ucell = ',Ucell*unitE,cunitE
                      write(6,'(I10,G10.3,A,G21.12,A)') it,timel,'*Ecell = ',EcellPR*unitE,cunitE
                      write(6,'(I10,G10.3,A,G21.12,A)') it,timel,'*Htot_PR = ', &
                           (potist+kine+EcellPR)*unitE,cunitE
                   END IF
                   write(6,*) 'Box tensor'
                   write(6,*)'a',at(1,1),at(2,1),at(3,1)
                   write(6,*)'b',at(1,2),at(2,2),at(3,2)
                   write(6,*)'c',at(1,3),at(2,3),at(3,3)

                   Call MatInv(h0, invh0)
                   Transformation=MatMul(at,invh0)
                   ! Strain tensor (Lagrange definition)
                   strain = 0.5d0*MatMul(Transformation,Transpose(Transformation))
                   DO i=1, 3
                      strain(i,i) = strain(i,i) - 0.5d0
                   END DO
                   rotation = 0.5d0*(Transformation - Transpose(Transformation))
                   WRITE(6,'(a)') 'Strain (Lagrange def.):'
                   WRITE(6,'(a,3g14.6)') '  e(1:3,1) = ', strain(1:3,1)
                   WRITE(6,'(a,3g14.6)') '  e(1:3,2) = ', strain(1:3,2)
                   WRITE(6,'(a,3g14.6)') '  e(1:3,3) = ', strain(1:3,3)
                   WRITE(6,'(a)') 'Rotation:'
                   WRITE(6,'(a,3g14.6)') '  r(1:3,1) = ', rotation(1:3,1)
                   WRITE(6,'(a,3g14.6)') '  r(1:3,2) = ', rotation(1:3,2)
                   WRITE(6,'(a,3g14.6)') '  r(1:3,3) = ', rotation(1:3,3)



                   a1 = at(1,1)
                   a2 = at(2,1)
                   a3 = at(3,1)
                   b1 = at(1,2)
                   b2 = at(2,2)
                   b3 = at(3,2)
                   c1 = at(1,3)
                   c2 = at(2,3)
                   c3 = at(3,3)
                   amod = dsqrt(a1**2+a2**2+a3**2)
                   bmod = dsqrt(b1**2+b2**2+b3**2)
                   cmod = dsqrt(c1**2+c2**2+c3**2)
                   fteta = 360d0/2d0/pi
                   tbc = fteta*dacos((b1*c1+b2*c2+b3*c3)/(bmod*cmod))
                   tca = fteta*dacos((c1*a1+c2*a2+c3*a3)/(cmod*amod))
                   tab = fteta*dacos((a1*b1+a2*b2+a3*b3)/(amod*bmod))
                   amod=amod*1.0d8 ; bmod=bmod*1.0d8 ; cmod=cmod*1.0d8
                   if(it<=1) then
                      volumean=volu
                      amodmean=amod ; bmodmean=bmod ; cmodmean=cmod
                      tcamean=tca; tabmean=tab; tbcmean=tbc
                   else
                      volumean=(volumean*(it/itetemp-1)+volu)/(it/itetemp)
                      amodmean=(amodmean*(it/itetemp-1)+amod)/(it/itetemp)
                      bmodmean=(bmodmean*(it/itetemp-1)+bmod)/(it/itetemp)
                      cmodmean=(cmodmean*(it/itetemp-1)+cmod)/(it/itetemp)
                      tcamean=(tcamean*(it/itetemp-1)+tca)/(it/itetemp)
                      tbcmean=(tbcmean*(it/itetemp-1)+tbc)/(it/itetemp)
                      tabmean=(tabmean*(it/itetemp-1)+tab)/(it/itetemp)
                   end if
                   write(6,'(I10,G10.3,A,2G18.9)') it,timel,'*l_a,m  ',amod,amodmean
                   write(6,'(I10,G10.3,A,2G18.9)') it,timel,'*l_b,m  ',bmod,bmodmean
                   write(6,'(I10,G10.3,A,2G18.9)') it,timel,'*l_c,m  ',cmod,cmodmean
                   write(6,'(I10,G10.3,A,2F11.4)') it,timel,'*ang_bc,m  ',tbc,tbcmean
                   write(6,'(I10,G10.3,A,2F11.4)') it,timel,'*ang_ca,m  ',tca,tcamean
                   write(6,'(I10,G10.3,A,2F11.4)') it,timel,'*ang_ab,m  ',tab,tabmean
                   write(6,'(I10,G10.3,A,2G21.12)') it,timel,'*volume  ',volu*1d24,volumean*1d24

                endif    ! if (lpr)

                IF (lTNose) THEN
                   WRITE(6,'(a)') 'Thermostat de Nose'
                   WRITE(6,'(I10,G10.3,A,G21.12,A,a,f0.3,a)') it,timel, &
                        '*KNose = ',KNose*unitE,cunitE, '  (', 2.d0*KNose/(bk), ' K)'
                   WRITE(6,'(i7,G10.3,a,g22.12)') it,timel,'*fNose = ', fNose
                   WRITE(6,'(i7,G10.3,a,g22.12,a)') it,timel,'*Htot_Nose = ', (kine+potist+EcellPR+KNose+UNose)*unitE,cunitE

                ELSEIF (lTHoover) THEN
                   WRITE(6,'(a)') 'Thermostat de Nose-Hoover'
                   WRITE(6,'(I10,G10.3,A,G21.12,A,a,f0.3,a)') it,timel, &
                        '*KNose = ',KNose*unitE,cunitE, '  (', 2.d0*KNose/(bk), ' K)'
                   WRITE(6,'(i7,G10.3,a,g22.12)') it,timel,'*zHoover(1) = ', zHoover(1)
                   !WRITE(6,'(i7,G10.3,a,g22.12,a)') it,timel,'*Htot_Hoover = ', (kine+potist+Ecell+KNose+UNose)*unitE,cunitE

                END IF   ! if lTNose / lTHoover

             endif                                ! rang=0
          end if
          IF (it<=1) THEN
             tmean = temp
             kinemean = kine
          ELSE
             tmean = (tmean*(it/itetemp-1)+temp)/(it/itetemp)
             kinemean = (kinemean*(it/itetemp-1)+kine)/(it/itetemp)
          END IF


          if (rang==0) then
             write (6, *)
             do iti = 1, ntyp
                if (na(iti)==0) cycle
                if (mod(it,itetemp2)==0) then
                   write (6, '(A,I2,A,F12.2)') &
                        '*temp instantanee des atomes de type', iti, ' = ', &
                        temptyp(iti)
                end if
             end do

             if (iteSigma>0) then
                if (mod(it,iteSigma)==0) then
                   if (mod(it,itetemp2)==0) then
                      write (6, *)
                      write (6, *) '* stress en ', cunitP
                   end if
                   ppot = 0.0
                   pkin = 0.0
                   do ic = 1, 3
                      if (mod(it,itetemp2)==0) then
                         write (6, '(I1,3(A,I1),A,3G18.10)') ic,' sigma potentiel (1,', ic, ') (2,', ic, &
                              ') (3,', ic, ') =',sig(1:3,ic)*unitP
                      end if
                      ppot = ppot+1.0/3.0*sig(ic,ic)
                   end do
                   if (mod(it,itetemp2)==0) then
                      write (6, *)
                   end if
                   do ic = 1, 3
                      if (mod(it,itetemp2)==0) then
                         write (6, '(I1,3(A,I1),A,3G18.10)') ic,' sigma cinetique (1,', ic, ') (2,', ic, &
                              ') (3,', ic, ') =',sigkine(1:3,ic)*unitP
                      end if
                      pkin = pkin+1.0/3.0*sigkine(ic,ic)
                   end do
                   if (mod(it,itetemp2)==0) then
                      write (6, *)
                      do ic = 1, 3
                         write (6, '(I1,3(A,I1),A,3G18.10)') ic,' sigma total (1,', ic, ') (2,', ic, &
                              ') (3,', ic, ') =',sigtot(1:3,ic)*unitP
                      end do

                      write (6, *)
                   end if


                   pist = ppot+pkin
                   if (mod(it,itetemp2)==0) then
                      write (6, '(I10,G10.3,A,G14.5)') it, timel, '*Pression = '&
                           , pist*unitP
                   end if
                   if(it<=1) then
                      pmean = pist
                   else
                      pmean = (pmean*(it/itesigma-1)+pist)/(it/itesigma)
                   end if
                   if (mod(it,itetemp2)==0) then
                      write (6, *) 'pression totale  potentiel = cinetique ='
                      write (6, '(3g14.5)') pist*unitP, ppot*unitP, pkin*unitP
                   end if




                endif
             endif                             !  itesigma>0

             IF (itefcc>0) THEN
                IF (mod(it,itefcc)==0) THEN
                   alat=(4.d0*volu/dble(im))**(1.d0/3.d0)   ! Lattice parameter
                   Allocate(fcc_nVoisins(1:im))
                   Allocate(fcc_cluster(1:im))
                   CALL Count_Fcc_Neighbours(xp, iwmax, alat, fcc_nVoisins, fcc_cluster, 6)
                   IF (Allocated(aux_int)) DeAllocate(aux_int)
                   Allocate(aux_int(2,1:im))
                   IF (Allocated(aux_title)) DeAllocate(aux_title)
                   Allocate(aux_title(2))
                   aux_title(1) = 'fcc neighbours'
                   aux_int(1,1:im) = fcc_nVoisins(1:im)
                   aux_title(2) = 'cluster index'
                   aux_int(2,1:im) = fcc_cluster(1:im)
                   Deallocate(fcc_cluster)
                   WRITE(out_file,'(2a,i0,a)') fnam(1:lenfnam),'.', it, '.fcc.cfg'
                   OPEN(file=out_file, unit=60, action='write')
                   CALL WriteCfg(xp, ityp, im, at, 60, fcc_nVoisins(1:im).NE.12, 2, aux_int, aux_title=aux_title)
                   CLOSE(60)
                   Deallocate(fcc_nVoisins)
                   Deallocate(aux_int, aux_title)
                END IF
             END IF

             if (mod(it,itetemp2)==0) then
                !               write (6, *)
                write (6, *) '----------valeurs moyennes------------'
                write (6, '(A,F12.2)') '*temperature moyenne = ', tmean
                !               write (6, '(A,G21.12,A)') '*energie cinetique moyenne = ', &
                !                  kinemean*unitE, cunitE
                if (itesigma>0) then
                   if (mod(it,itesigma)==0) write (6, '(A,G14.5)') &
                        '*pression moyenne = ', pmean*unitP
                endif
                write (6, *) '--------------------------------------'
             end if




             if (ltpcel) then
                minp=100000
                maxp=-100000
                mint=100000
                maxt=-100000
                !           write(extension,'(i9.9)') it
                !           lenfn2 = 9
                !           if (it.ge.3) then
                !              open(luvisuc, file=fnam(1:lenfnam)//'.'//extension(1:lenfn2)//'.CEL.mol', form='formatted', &
                !                   status='unknown')
                !              write (luvisuc, '(I9,A,I7,A,D15.6)') 2511 , ' IT =', it, ' Time = ', timel
                !              at=at*1.d8
                !              write (luvisuc,'(9F12.6)')at(1,1),at(2,1),at(3,1),at(1,2),at(2,2),at(3,2),at(1,3),at(2,3),at(3,3)
                !              at=at/1.d8
                !           end if



                nprt=0
                do kx=0,nox-1
                   do ky=0,noy-1
                      do kz=0,noz-1
                         ko=1+kx+nox*(ky+noy*kz)
                         xb(1)=float(kx)/float(nox)*at(1,1)+float(ky)/float(noy)*at(1,2)+float(kz)/float(noz)*at(1,3)
                         xb(2)=float(kx)/float(nox)*at(2,1)+float(ky)/float(noy)*at(2,2)+float(kz)/float(noz)*at(2,3)
                         xb(3)=float(kx)/float(nox)*at(3,1)+float(ky)/float(noy)*at(3,2)+float(kz)/float(noz)*at(3,3)
                         xb=xb*1d8
                         pmc(ko)=0.0
                         if (itesigma.gt.0) then
                            if (mod(it,itesigma).eq.0) then
                               if (ltabvois) then
                                  continue
                               else
                                  !                                write(6,*)'dans la celulle ',ko,' sigma  '
                                  do ic =1,3
                                     pmc(ko)=pmc(ko)+sigc(ic,ic,ko)/3.0



                                     !                                write(6,'(3g14.5)')sigc(1,ic,ko),sigc(2,ic,ko),sigc(3,ic,ko)
                                     !                                write(6,'(A,I2,I2,I2,I2,G14.5,G14.5,G14.5)')'CEL-SIG ',ic,kx,ky&
                                     !                                     &,kz,sigc(1,ic,ko),sigc(2,ic,ko),sigc(3,ic,ko)
                                  enddo
                                  !                             write(6,'(A,I5,3I4,G14.5,A,A,I4)')'CEL-PRESS ', ko,kx,ky,kz,pmc*unitP, '  ',cunitP,nato(ko)
                                  celpP(ko)=1d-14*unitP*(pmc(ko)-celpm1(ko))/(timel-timelm1)
                                  !                             celpP2(ko)=1d-14*1d-14*unitP*(pmc(ko)-2*celpm1(ko)+celpm2(ko))/(tstep**2)
                                  !                             celpp2(ko)=1d-14*1d-14*unitP*((pmc(ko)-celpm1(ko))/tstep -(celpm1(ko)-celpm2(ko))/timelm1)/tstep

                                  tcp(ko)=1d-14*(tempc(ko)-tm1(ko))/(timel-timelm1)
                                  !                             tcp2(ko)=1d-14*1d-14*(tempc(ko)-2*tm1(ko)+tm2(ko))/(tstep**2)
                                  !                             tcp2(ko)=1d-14*1d-14*((tempc(ko)-tm1(ko))/tstep -(tm1(ko)-tm2(ko))/timelm1)/tstep
                                  !                             celpm2(ko)=celpm1(ko)
                                  celpm1(ko)=pmc(ko)
                                  !                             tm2(ko)=tm1(ko)
                                  tm1(ko)=tempc(ko)
                                  lprtcel(ko)=.false.

                                  !                             if ((kx==15).or.(ky==15).or.(kz==15).or.(kx==5).or.(ky==5).or.(kz==5)) lprt=.true.
                                  !                             if((kx.le.13).and.(kx.ge.7).and.(ky.le.13).and.(ky.ge.7).and.(kz.le.13).and.(kz.ge.7)) lprt=.true.
                                  !                             if (it.le.2) lprt=.false.
                                  if(tempc(ko).gt.tpseuils(1)) lprtcel(ko)=.true.
                                  if(abs(tcp(ko)).gt.tpseuils(2)) lprtcel(ko)=.true.
                                  !                             if(abs(tcp2(ko)).gt.tpseuils(3)) lprtcel(ko)=.true.
                                  if(abs(pmc(ko)*unitP).gt.tpseuils(3)) lprtcel(ko)=.true.
                                  if(abs(celpp(ko)).gt.tpseuils(4)) lprtcel(ko)=.true.
                                  !                             if(abs(celpp2(ko)).gt.tpseuils(6)) lprtcel(ko)=.true.
                                  if (it.le.2) lprtcel(ko)=.false.
                                  !                             lprtcel(ko)=.true.
                                  write(6,*)'tempcprt', tempc(ko),abs(tcp(ko))
                                  if (lprtcel(ko).EQV..true.) nprt=nprt+1

                                  if(it.ge.3) then
                                     minp=min(minp,pmc(ko)*unitP)
                                     maxp=max(maxp,pmc(ko)*unitP)
                                     mint=min(mint,tempc(ko))
                                     maxt=max(maxt,tempc(ko))

                                     minpP=min(minpP,celpp(ko))
                                     maxpP=max(maxpP,celpp(ko))
                                     mintP=min(mintP,tcp(ko))
                                     maxtP=max(maxtP,tcP(ko))

                                     !                                minpP2=min(minpP2,celpP2(ko))
                                     !                                maxpP2=max(maxpP2,celpP2(ko))
                                     !                                mintP2=min(mintP2,tcp2(ko))
                                     !                                maxtP2=max(maxtP2,tcp2(ko))
                                  end if

                               endif
                            endif
                         endif
                      enddo
                   end do
                end do
                timelm1=timel

                if ((mod(it,itetemp2)==0).and.(nprt.gt.0))then
                   write (6, *) '------valeurs par cellules-------',it,nprt
                   write(extension,'(i9.9)') it
                   lenfn2 = 9
                   open(luvisuc, file=fnam(1:lenfnam)//'.'//extension(1:lenfn2)//'.CEL.mol', form='formatted', &
                        status='unknown')
                   write (luvisuc, '(I9,A,I7,A,D15.6)') nprt , ' IT =', it, ' Time = ', timel
                   at=at*1.d8
                   write (luvisuc,'(9F12.6)')at(1,1),at(2,1),at(3,1),at(1,2),at(2,2),at(3,2),at(1,3),at(2,3),at(3,3)
                   at=at/1.d8
                   do kx=0,nox-1
                      do ky=0,noy-1
                         do kz=0,noz-1
                            ko=1+kx+nox*(ky+noy*kz)
                            xb(1)=float(kx)/float(nox)*at(1,1)+float(ky)/float(noy)*at(1,2)+float(kz)/float(noz)*at(1,3)
                            xb(2)=float(kx)/float(nox)*at(2,1)+float(ky)/float(noy)*at(2,2)+float(kz)/float(noz)*at(2,3)
                            xb(3)=float(kx)/float(nox)*at(3,1)+float(ky)/float(noy)*at(3,2)+float(kz)/float(noz)*at(3,3)
                            xb=xb*1d8
                            if(lprtcel(ko).EQV..true.)  write (luvisuc, 136) 'Au',kx,ky,kz,xb(1), xb(2), &
                                 xb(3),tempc(ko),tcp(ko),pmc(ko)*unitP,celpp(ko),nato(ko)
                         end do
                      end do
                   end do
                   close (luvisuc)
                end if

                if (it.ge.3)then
                   Cminp=min(Cminp,minp)
                   Cmaxp=max(Cmaxp,maxP)
                   Cmint=min(Cmint,mint)
                   Cmaxt=max(Cmaxt,maxT)

                   CminpP=min(CminpP,minpP)
                   CmaxpP=max(CmaxpP,maxpP)
                   CmintP=min(CmintP,mintP)
                   CmaxtP=max(CmaxtP,maxTP)

                   !              CminpP2=min(CminpP2,minpP2)
                   !              CmaxpP2=max(CmaxpP2,maxPP2)
                   !              CmintP2=min(CmintP2,mintP2)
                   !              CmaxtP2=max(CmaxtP2,maxTP2)

                   write(6,'(A,4G20.10)')'minmaxp', minp,maxp,mint,maxt
                   write(6,'(A,4G20.10)')'Cminmaxp', Cminp,Cmaxp,Cmint,Cmaxt
                   write(6,'(A,4G20.10)')'minmaxpP', minpP,maxpP,mintP,maxtP
                   write(6,'(A,4G20.10)')'CminmaxpP', CminpP,CmaxpP,CmintP,CmaxtP
                   !              write(6,'(A,4G20.10)')'minmaxpP2', minpP2,maxpP2,mintP2,maxtP2
                   !              write(6,'(A,4G20.10)')'CminmaxpP2', CminpP2,CmaxpP2,CmintP2,CmaxtP2
                end if
             endif

136          format(A,3I4,3E15.5,4E15.7,I4)

             if (ldesinteg)then


                write(6,'(A,3G21.12)')'lambda, deltaF',lambdades,deltaF,deltaF*erg2eV
                write(6,'(A,3G21.12)')'deltaEspr, Espr', lambdades,deltaEspr*erg2eV,Espr*erg2eV
             end if

          endif

       endif
    endif





    ! calcul des deplacements
    !  if (rang==0) write(6,*) 'PARA-T itedepla' ,itedepla
    if (itedepla>0) then
       if (mod(it,itedepla)==0) then
          call calcdepla(im,xp,ielat,ityp,ax)
          if (tdepla2>0.0) call calcdepla2(im,xp,ielat,ityp,ax)
       endif
    endif


    ! calcul des coordinences
    if (itecoordo>0) then
       if (mod(it,itecoordo)==0) call calccoordo  (im,imm,ityp,xp,ielat)
    endif


    ! calcul de la RDF et de la position moyenne
    if (iterdf>0) then
       if (mod(it,iterdf)==0) then
          call calcdigr (im,xp,ityp,ielat)
          nrdf = nrdf+1
          if (linstantrdf) then
             call rdf 
             nrdf = 0
          endif
       endif
    else if (iterdf==0) then
       if (itmax-it<nrdf) then
          call calcdigr (im,xp,ityp,ielat)
          nrdf = nrdf+1
          !        write(6,*)'sortie digr de analyse'
       endif
    endif

    !                 CALCUL DES DISTRIBUTIONS ANGULAIRES

    if (iteangle>0) then
       if (mod(it,iteangle)==0) then
          call calcangle (im,imm,ityp,xp,ielat)
          nfda=nfda+1 
          if (linstantfda) then
             call adf
             nfda=0
          endif
       endif
    else if(iteangle==0) then
       if(itmax-it<nfda) then
          call calcangle (im,imm,ityp,xp,ielat)
          nfda=nfda+1
       endif
    endif
    ! -----------------------------------------------

    ! ecriture de CFG

    if (itecfg>0) then
       if (mod(it,itecfg)==0) then

          WRITE(out_file,'(2a,i0,a)') fnam(1:lenfnam),'.', it, '.cfg'
          OPEN(file=out_file, unit=60, action='write')

          !WHAT_THE_HACK_IS_THAT        if (dmtype==17)  CALL redefine_ty()

          IF (lPrtEat.OR.lPrtSigat) THEN       ! Energy and/or stress per atom
             nAux_real=0
             IF (lPrtEat)   nAux_real = nAux_real + 1
             IF (lPrtSigat) nAux_real = nAux_real + 6
             IF (Allocated(aux_real)) DeAllocate(aux_real)
             Allocate(aux_real(nAux_real,1:im))
             IF (Allocated(aux_title)) DeAllocate(aux_title)
             Allocate(aux_title(nAux_real))
             n=0
             IF (lPrtEat) THEN
                aux_title(n+1)="Energy per atom (eV)"
                aux_real(n+1,1:im)=Eatom(1:im)*erg2eV
                n = n+1
             END IF
             IF (lPrtSigat) THEN
                aux_title(n+1) = 'Stress Sxx (' // cunitP // ')'
                aux_title(n+2) = 'Stress Syy (' // cunitP // ')'
                aux_title(n+3) = 'Stress Szz (' // cunitP // ')'
                aux_title(n+4) = 'Stress Syz (' // cunitP // ')'
                aux_title(n+5) = 'Stress Sxz (' // cunitP // ')'
                aux_title(n+6) = 'Stress Sxy (' // cunitP // ')'
                aux_real(n+1,1:im) = sigat(1,1,1:im)*unitP
                aux_real(n+2,1:im) = sigat(2,2,1:im)*unitP
                aux_real(n+3,1:im) = sigat(3,3,1:im)*unitP
                aux_real(n+4,1:im) = 0.5d0*( sigat(2,3,1:im) + sigat(3,2,1:im) )*unitP
                aux_real(n+5,1:im) = 0.5d0*( sigat(1,3,1:im) + sigat(3,1,1:im) )*unitP
                aux_real(n+6,1:im) = 0.5d0*( sigat(1,2,1:im) + sigat(2,1,1:im) )*unitP
             END IF
             CALL WriteCfg(xp, ityp, im, at, 60, nAux_real=nAux_real, aux_real=aux_real, aux_title=aux_title)
             DEALLOCATE(aux_real, aux_title)
          ELSE
             CALL WriteCfg(xp, ityp, im, at, 60)
          END IF
          CLOSE(60)
       endif

       !WHAT_THE_HACK_IS_THAT     if (dmtype==17) CALL refix_ty()

    endif                                      ! fin rang=0



    ! ecriture de rasmol

    if (iterasmol>0) then     
       if (mod(it,iterasmol)==0) then
     call ndm2cellconfig(celndm,noxyz,nox,noy,noz,natperc,nato,ncel,atincel,deltadist,celsize)
    call ndm2config(atdml,im,imm,xp,fp,ityp,ielat,num_at_glob=num_at_glob,ltabvois=ltabvois,&
         &iwmax=iwmax,indi=indi,nvois=nvois,vp=vp,xpp=xpp,eat=eatom)
    call ndm2boxconfig(at,bg,zl,zls2,nzl,volu,normat,boxndm)
         call rasmol(atdml,boxndm,it)
          if (l2T) call  eleccellmol

       end if
    endif




    !     write(6,*)'sortie ssprogramme analyse'
    !crcmit      call flush(6)



    if (itebdv>0) then
       ! Pas pris en compte en parallele
       if (.not.parallele.and.mod(it,itebdv)==0) call bondval(im,imm,xp,ityp,ielat,num_at_glob)
    endif

    if (iteanapos>0) then
       ! Pas pris en compte en parallele
       if (.not.parallele.and.mod(it,iteanapos)==0) call anapos(it)
    endif


    if(iteplz>0.and.(.not.parallele)) then
       call prtplz(xp,ityp)
    end if
    if (lposmoy.eqv..true.) then
       nposmoy=nposmoy+1
       do i=1,im
          posmoyx(1:3,i)=(xp(1:3,i)+(nposmoy-1)*posmoyx(1:3,i))/nposmoy
       end do
    end if

    if (lprteattotm.EQV..true.) then
       if (ncalceattotm==0)eatomtotm(:)=0.
       ncalceattotm=ncalceattotm+1
       do i=1,im
          eatomtotm(i)=(eatom(i)+(ncalceattotm-1)*eatomtotm(i))/ncalceattotm
       end do
       !     write(912,*)ncalceattotm,eatomtotm(1)*erg2eV,eatom(1)*erg2eV
    end if

    !  if (lbulle) then
    !     if (parallele) then
    !        write(6,*)' test pos He pas para'
    !        stop
    !     end if
    !     call test_position_He(xp,at,ityp,rang,imm,im,it,ldesinteg,num_at_glob,nstepdes,itmax)
    !  end if
    !  write(6,*) 'sortie canalyse',it,im

    if (mod(it,itesigma)==0) then
       if (lsigatcel)then
          natchk(:)=0
          sigatcel=0 ; patcel=0 ; patcelmax=0
          do i=1,im
             koo = ielat(i)                          ! Numero de la cellule
             iti = ityp(i)
             natchk(koo)=natchk(koo)+1
             sigatcel(:,:,koo)=sigatcel(:,:,koo)+sigat(:,:,i)
             ptest=0
             do ic=1,3
                ptest=ptest+sigat(ic,ic,i)/3
             end do
             ptest=abs(ptest)
             patcelmax(koo)=max(patcelmax(koo),ptest)
          end do
          do koo=1,noxyz
             sigatcel(:,:,koo)=sigatcel(:,:,koo)/natchk(koo)
             do ic=1,3
                patcel(koo)= patcel(koo)+sigatcel(ic,ic,koo)/3.0
             end do
             if (natchk(koo).ne.nato(koo)) then
                write(6,*)'nato ? koo natcchk nato', koo, natchk(koo), nato (koo)
                stop

             end if
          end do


          !        write (6, *) '------valeurs par cellules-------',it
          write(extension,'(i9.9)') it
          lenfn2 = 9
          open(luvisuc, file=fnam(1:lenfnam)//'.'//extension(1:lenfn2)//'.CEL.mol', form='formatted', &
               status='unknown')
          write (luvisuc, '(I9,A,I7,A,D15.6)') noxyz , ' IT =', it, ' Time = ', timel
          at=at*1.d8
          write (luvisuc,'(9F12.6)')at(1,1),at(2,1),at(3,1),at(1,2),at(2,2),at(3,2),at(1,3),at(2,3),at(3,3)
          at=at/1.d8
          do kx=0,nox-1
             do ky=0,noy-1
                do kz=0,noz-1
                   koo=1+kx+nox*(ky+noy*kz)
                   xb(1)=float(kx)/float(nox)*at(1,1)+float(ky)/float(noy)*at(1,2)+float(kz)/float(noz)*at(1,3)
                   xb(2)=float(kx)/float(nox)*at(2,1)+float(ky)/float(noy)*at(2,2)+float(kz)/float(noz)*at(2,3)
                   xb(3)=float(kx)/float(nox)*at(3,1)+float(ky)/float(noy)*at(3,2)+float(kz)/float(noz)*at(3,3)
                   xb=xb*1d8
                   write (luvisuc, 149) 'Au',xb(1), xb(2),xb(3),patcel(koo)*unitP,patcelmax(koo)*unitP,nato(koo)
                end do
             end do
          end do
          close (luvisuc)
148       format(A,3I4,3E15.5,1E15.7,I4)
149       format(A,3E15.5,2E15.7,I4)
       end if
    end if

    return
  end subroutine analyse
end module analyse_mod
