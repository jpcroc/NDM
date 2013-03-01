! ************************************************
!         Sous-programme analyse.f
! ************************************************

subroutine analyse
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use tab_imm_m
#if(PARA)
  use mod_mpi
#endif
  USE fcc_module
  USE cfg_module
  USE posana
  implicit none
  !-----------------------------------------------
  !   G l o b a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s        tersoff_zbl.o\
  !-----------------------------------------------
  integer :: i, iti, ic, ko,kx,ky,kz
  real(double), dimension(ntyp) :: temptyp



  real(double) :: ppot, pkin
  real(double) :: a1, a2, a3, b1, b2, b3, c1, c2, c3
  real(double) :: fteta, tbc, tca, tab, amod, bmod, cmod
  real(double) :: unitE,unitP
  character*5 :: cunitE, cunitP
  real(double) :: temp2,pmc, alat
  real(double), external :: tempinst
  real(double), save :: zlm(3)
  real(double), save :: volumean,amodmean,bmodmean,cmodmean,tcamean,tabmean,tbcmean
  real(double), dimension(3,3) :: transformation, strain, rotation, invh0

  INTEGER, dimension(:), allocatable :: fcc_nVoisins, fcc_cluster
  INTEGER, dimension(:,:), allocatable :: aux_int
  REAL(kind(0.d0)), dimension(:,:), allocatable :: aux_real
  CHARACTER(len=20), dimension(:), allocatable :: aux_title
  CHARACTER(len=100) :: out_file
  integer,save::ncalceattotm=0, nposmoy=0
  integer::ipot
  !-----------------------------------------------
  !
  !
  !      logical:: lEev=.false., lPkbar=.false.
!        write(6,*)'analyse',itetemp,it
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
  if (itetemp2>0) then
     if (mod(it,itetemp2)==0) then
        temp2=tempinst(vp,ityp)
        write(112,'(I10,G15.6,F12.2)')it,timel,temp2
     end if
  end if

  if (itetemp>0) then
     if (mod(it,itetemp)==0) then

        call calctemp (temptyp)
        !        write(6,*) 'sortie calctemp',temp,rang

        ! MPI
        !remarque 1erg = 6.24d11 eV
        if (rang==0) then
           write (6, *)
           write (6, *)
           write (6, '(A,I7,A,G10.3)') '<<<<<<<<  ITERATION =', it, &
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
                 end select
              end if
           end do
           write (6,'(I10,G10.3,A,G21.12,A,a,f0.3,a)') it,timel,'*Ec = ',kine*unitE, cunitE, &
                '  (', 2.d0*kine/(3.d0*float(im_glob)*bk), ' K)'
           write (6,'(I10,G10.3,A,G21.12,A)') it,timel,'*Etot = ',(kine+potist)*unitE, cunitE
           write (6, *)

!	IF(it==1) Open(unit=774, file='energie_it.dat', status='unknown', action='write')
!	IF(modulo(it,10)==0) WRITE(774,'(i6, f)') it, (kine+potist)*unitE  ! controle .. !*!

	   !*! ---------------------------------- APPEL DE LA ROUTINE D'ECRITURE -------------
	   ! ------------------------------------ DES FICHIERS DE SORTIE, DANS LE ------------
	   ! ------------------------------------ CAS DES CL CONTROLEES EN CONTRAINTE --------
	   IF (ibound==2 .OR. ibound==3) THEN
		IF(flag_fin==.true.) THEN
			Call spebc_fin (.true.)
		ELSE IF(itespebcout > 0. .AND. (mod(it,itespebcout)==0 .OR. it==1)) THEN
			Call spebc_fin (.false.)
		END IF
	   END IF
	   !*!
	

           if (lEparat) then
              write (6, '(I10,G10.3,A,G21.12,A)') it, timel, '*Epot/at = ', potist*unitE/im, cunitE
              write(6,*)
           end if

           write (6, '(I10,G10.3,A,F12.2)') it,timel, &
                '*Temp instantanee = ',temp

           if (tfcou>0.0) write (6, '(A,F12.2)') '*temperature externe = ', tcou

           IF (lpr) THEN
              write(6,*) 'NPT With Parrinello-Rahman'
              write(6,'(I10,G10.3,A,G21.12,A,a,f0.3,a)') it,timel,'*Kcell = ',Kcell*unitE,cunitE, &
                   '  (', 2.d0*Kcell/(9.d0*bk), ' K)'
              IF (lUcell) THEN
                 write(6,'(I10,G10.3,A,G21.12,A)') it,timel,'*Ucell = ',Ucell*unitE,cunitE
                 write(6,'(I10,G10.3,A,G21.12,A)') it,timel,'*Ecell = ',Ecell*unitE,cunitE
                 write(6,'(I10,G10.3,A,G21.12,A)') it,timel,'*Htot_PR = ', &
                      (potist+kine+Ecell)*unitE,cunitE
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
              WRITE(6,'(i7,G10.3,a,g22.12,a)') it,timel,'*Htot_Nose = ', (kine+potist+Ecell+KNose+UNose)*unitE,cunitE

           ELSEIF (lTHoover) THEN
              WRITE(6,'(a)') 'Thermostat de Nose-Hoover'
              WRITE(6,'(I10,G10.3,A,G21.12,A,a,f0.3,a)') it,timel, &
                   '*KNose = ',KNose*unitE,cunitE, '  (', 2.d0*KNose/(bk), ' K)'
              WRITE(6,'(i7,G10.3,a,g22.12)') it,timel,'*zHoover(1) = ', zHoover(1)
              !WRITE(6,'(i7,G10.3,a,g22.12,a)') it,timel,'*Htot_Hoover = ', (kine+potist+Ecell+KNose+UNose)*unitE,cunitE

           END IF   ! if lTNose / lTHoover

        endif                                ! rang=0

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
              write (6, '(A,I2,A,F12.2)') &
                   '*temp instantanee des atomes de type', iti, ' = ', &
                   temptyp(iti)
           end do

           if (itesigma>0) then
              if (mod(it,itesigma)==0) then

                 write (6, *)
                 write (6, *) '* stress en ', cunitP
                 ppot = 0.0
                 pkin = 0.0
                 do ic = 1, 3
                    write (6, '(I1,3(A,I1),A,3G12.4)') ic,' sigma potentiel (1,', ic, ') (2,', ic, &
                         ') (3,', ic, ') =',sig(1:3,ic)*unitP
                    ppot = ppot+1.0/3.0*sig(ic,ic)
                 end do
                 write (6, *)
                 do ic = 1, 3
                    write (6, '(I1,3(A,I1),A,3G12.4)') ic,' sigma cinetique (1,', ic, ') (2,', ic, &
                         ') (3,', ic, ') =',sigkine(1:3,ic)*unitP
                    pkin = pkin+1.0/3.0*sigkine(ic,ic)
                 end do
                 write (6, *)
                 do ic = 1, 3
                    write (6, '(I1,3(A,I1),A,3G12.4)') ic,' sigma total (1,', ic, ') (2,', ic, &
                         ') (3,', ic, ') =',sigtot(1:3,ic)*unitP
                 end do

                 write (6, *)


                 pist = ppot+pkin
                 write (6, '(I10,G10.3,A,G14.5)') it, timel, '*Pression = '&
                      , pist*unitP
                 if(it<=1) then
                    pmean = pist
                 else
                    pmean = (pmean*(it/itesigma-1)+pist)/(it/itesigma)
                 end if
                 write (6, *) 'pression totale  potentiel = cinetique ='
                 write (6, '(3g14.5)') pist*unitP, ppot*unitP, pkin*unitP





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
                 CALL WriteCfg(xp, ityp, 60, fcc_nVoisins(1:im).NE.12, 2, aux_int, aux_title=aux_title)
                 CLOSE(60)
                 Deallocate(fcc_nVoisins)
                 Deallocate(aux_int, aux_title)
              END IF
           END IF


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




        endif                                ! fin rang=0
        if (ltpcel) then

           write (6, *)
           write (6, *) '----------valeurs par cellules------------'

           do kx=0,nox-1
              do ky=0,noy-1
                 do kz=0,noz-1
                    ko=1+kx+nox*(ky+noy*kz)
                    pmc=0.0
                    if (itesigma.gt.0) then
                       if (mod(it,itesigma).eq.0) then
                          if (ltabvois) then
                             continue
                          else
!                             write(6,*)'dans la celulle ',ko,' sigma  '
                             do ic =1,3
                                pmc=pmc+sigc(ic,ic,ko)/3.0
!                                write(6,'(3g14.5)')sigc(1,ic,ko),sigc(2,ic,ko),sigc(3,ic,ko)
!                                write(6,'(A,I2,I2,I2,I2,G14.5,G14.5,G14.5)')'CEL-SIG ',ic,kx,ky&
!                                     &,kz,sigc(1,ic,ko),sigc(2,ic,ko),sigc(3,ic,ko)
                             enddo
                             write(6,'(A,I5,3I4,G14.5,A,A,I4)')'CEL-PRESS ', ko,kx,ky,kz,pmc*unitP, '  ',cunitP,nato(ko)


                          endif
                       endif
                    endif
                 enddo
              end do
           end do
        endif
  if (ldesinteg)then
     if (rang==0) then
        
        write(6,'(A,3G21.12)')'lambda, deltaF',lambdades,deltaF,deltaF*erg2eV
        write(6,'(A,3G21.12)')'deltaEspr, Espr', lambdades,deltaEspr*erg2eV,Espr*erg2eV
     end if
  end if


     endif
  endif
  if ((mod(it,itesigma)==0).and.(lsigtyp)) then

     write(6,*)

     do iti=1,ntyp
        pmc=0.
        do ic =1,3
           pmc=pmc+sigtyp(ic,ic,iti)/3.0 
        end do
        if (rang==0) write(6,'(A,G14.5)')'sigtyp ',iti,' = ',pmc
        pmc=0.
        do ic =1,3
           pmc=pmc+sigtyptyp(ic,ic,iti,iti)/3.0 
        end do
        if (rang==0) write(6,'(A,G14.5)')'sigtyptyp ',iti,iti,' = ',pmc
     end do
  end if

  

  ! calcul des deplacements
  !  if (rang==0) write(6,*) 'PARA-T itedepla' ,itedepla
  if (itedepla>0) then
     if (mod(it,itedepla)==0) then
        call calcdepla
        if (tdepla2>0.0) call calcdepla2
     endif
  endif


  ! calcul des coordinences
  if (itecoordo>0) then
     if (mod(it,itecoordo)==0) call calccoordo 
  endif


  ! calcul de la RDF et de la position moyenne
  if (iterdf>0) then
     if (mod(it,iterdf)==0) then
        call calcdigr 
        nrdf = nrdf+1
        if (linstantrdf) then
           call rdf 
           nrdf = 0
        endif
     endif
  else if (iterdf==0) then
     if (itmax-it<nrdf) then
        call calcdigr 
        nrdf = nrdf+1
        !        write(6,*)'sortie digr de analyse'
     endif
  endif

  !                 CALCUL DES DISTRIBUTIONS ANGULAIRES

  if (iteangle>0) then
     if (mod(it,iteangle)==0) then
        call calcangle 
        nfda=nfda+1 
        if (linstantfda) then
           call adf
           nfda=0
        endif
     endif
  else if(iteangle==0) then
     if(itmax-it<nfda) then
        call calcangle 
        nfda=nfda+1
     endif
  endif
  ! -----------------------------------------------

  ! ecriture de CFG

  if (itecfg>0) then
     if (mod(it,itecfg)==0) then

        WRITE(out_file,'(2a,i0,a)') fnam(1:lenfnam),'.', it, '.cfg'
        OPEN(file=out_file, unit=60, action='write')
        IF (lprteat) THEN       ! Energy per atom
           IF (Allocated(aux_real)) DeAllocate(aux_real)
           Allocate(aux_real(1,1:im))
           IF (Allocated(aux_title)) DeAllocate(aux_title)
           Allocate(aux_title(2))
           aux_title(1)="Energy per atom (eV)"
           aux_real(1,1:im)=Eatom(1:im)*erg2eV
           CALL WriteCfg(xp, ityp, 60, nAux_real=1, aux_real=aux_real, aux_title=aux_title)
           DEALLOCATE(aux_real, aux_title)
        ELSE
           CALL WriteCfg(xp, ityp, 60)
        END IF
        CLOSE(60)
     endif

  endif                                      ! fin rang=0



  ! ecriture de rasmol

  if (iterasmol>0) then
     if (mod(it,iterasmol)==0) call rasmol (it)
  endif




  !     write(6,*)'sortie ssprogramme analyse'
  !crcmit      call flush(6)



  if (itebdv>0) then
     ! Pas pris en compte en parallele
     if (.not.parallele.and.mod(it,itebdv)==0) call bondval
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

  if (lbulle) then
     if (parallele) then
        write(6,*)' test pos He pas para'
        stop
     end if
     
     call test_position_He(xp,at,ityp,rang,imm,im,it,ldesinteg,num_at_glob,nstepdes,itmax)
  end if
!  write(6,*) 'sortie canalyse',it,im
  return
end subroutine analyse
