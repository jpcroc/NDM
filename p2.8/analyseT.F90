module analyseT_mod
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
  USE rasmolT_mod,only: rasmolT
  USE rdf_mod,only: rdf
  USE prtplz_mod,only: prtplz
  USE sauvegardeT_mod,only:sauvegardeT
 USE sauveforce_mod,only: sauveforce

  use var_pot, only: iewald,l3c,npotmax,potisglue,potisrep,lpotentiel,ntyp,na
  use gen_com_m, only:bk,cunite,deltaespr,deltaf,ecellpr,espr,flag_fin,fnose,iteanapos,iteangle,itebdv,&
       &itecfg,itecoordo,itedepla,itefcc,iterasmol,iterdf,itesigma,itetemp,itetemp2,kcell,kine,kinemean,knose,&
       &lambdades,leev,leparat,linstantfda,lpr,lprteattotm,lsigatcel,lthoover,ltnose,ltpcel,lucell,&
       &nfda,parallele,pist,pmean,potcp,potis1,potis2,potis3,potist,potistersoff,potiszbl,&
       &tcou,temp,tempep,tfcou,tmean,ucell,unite,unose,zhoover,sig,sigkine,lprtcel,&
       &natchk,tpseuils,sigtot,unitP,tdepla2,nrdf,lprtsigat,lprteat,lpkbar,linstantrdf,&
       &ldesinteg,itmax,cunitp,erg2ev,iteplz,itespebcout,lperiod,pi,rang,timel,&
       & itesauvforce,itesauv,formatsauv,fnamcout,itesauvinter,itesauvposition,fnam,lenfnam,im_glob,it,l2T

  USE cellconfig,only:cell_config, caltabtC
  USE atomconfig,only:atom_config,atom_config_d,atom_config_e
  use boxconfig,only: box_config
  implicit none
contains
  ! ************************************************
  !         Sous-programme analyse.f
  ! ************************************************

  subroutine analyseT(atdml,celndm,boxndm)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double

!    use tab_imm_m

!    USE fcc_module
!    USE cfg_module
!    USE posana
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


    class(atom_config_d)::atdml
    type(cell_config):: celndm
    type(box_config)::boxndm

    type(atom_config_d)::attyp
    type(cell_config):: celtyp

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
    ! MPI
    if (rang==0) then
       !          write(6,*)'analyse -> sauvegarde'
       if (itesauv.GT.0) then
          if (mod(it,itesauv)==0) then 
             formatsauv = 3
             if(itesauvinter.gt.0) then
                if (mod(it,itesauvinter).eq.0) then
                   write(extension,'(i9.9)') it
                   fnamcout = fnam(1:lenfnam)//'.cout.'//extension
                else
                   fnamcout = fnam(1:lenfnam)//'.cout'
                endif
             else
                fnamcout = fnam(1:lenfnam)//'.cout'
             end if
             call sauvegardeT(atdml,celndm,boxndm,formatsauv,fnamcout)
          end if
       endif

              if (itesauvposition.GT.0) then
          if (mod(it,itesauvposition)==0) then
             formatsauv = 2
             write(extension,'(i9.9)') it
             fnamcout = fnam(1:lenfnam)//'.cout.'//extension
             call sauvegardeT(atdml,celndm,boxndm,formatsauv,fnamcout)
          end if
       endif
       
       if (itesauvforce.GT.0) then
          if (mod(it,itesauvforce)==0) call sauveforce ( it)
       endif
       !          write(6,*)'sauvposition -> control'
    endif                                   ! fin rang=0

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
             where(atdml%ityp(1:atdml%im)==iti)
                atdml%lgul(1:atdml%im)=.true.
             end where
             call atdml%fab(attyp)
             call caltabtC(celtyp,attyp,lperiod,boxndm)
             call calctemp(temptyp(iti),kinetyp,attyp,celtyp)
             call celtyp%dealloc ; call attyp%dealloc
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
!!$                IF (ibound==2 .OR. ibound==3) THEN
!!$                   IF(flag_fin.EQV..true.) THEN
!!$                      Call spebc_fin (.true.)
!!$                   ELSE IF(itespebcout > 0. .AND. (mod(it,itespebcout)==0 .OR. it==1)) THEN
!!$                      Call spebc_fin (.false.)
!!$                   END IF
!!$                END IF
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
                   write(6,*)'a',boxndm%at(1,1),boxndm%at(2,1),boxndm%at(3,1)
                   write(6,*)'b',boxndm%at(1,2),boxndm%at(2,2),boxndm%at(3,2)
                   write(6,*)'c',boxndm%at(1,3),boxndm%at(2,3),boxndm%at(3,3)

                   Call MatInv(boxndm%h0, invh0)
                   Transformation=MatMul(boxndm%at,invh0)
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



                   a1 = boxndm%at(1,1)
                   a2 = boxndm%at(2,1)
                   a3 = boxndm%at(3,1)
                   b1 = boxndm%at(1,2)
                   b2 = boxndm%at(2,2)
                   b3 = boxndm%at(3,2)
                   c1 = boxndm%at(1,3)
                   c2 = boxndm%at(2,3)
                   c3 = boxndm%at(3,3)
                   amod = dsqrt(a1**2+a2**2+a3**2)
                   bmod = dsqrt(b1**2+b2**2+b3**2)
                   cmod = dsqrt(c1**2+c2**2+c3**2)
                   fteta = 360d0/2d0/pi
                   tbc = fteta*dacos((b1*c1+b2*c2+b3*c3)/(bmod*cmod))
                   tca = fteta*dacos((c1*a1+c2*a2+c3*a3)/(cmod*amod))
                   tab = fteta*dacos((a1*b1+a2*b2+a3*b3)/(amod*bmod))
                   amod=amod*1.0d8 ; bmod=bmod*1.0d8 ; cmod=cmod*1.0d8
                   if(it<=1) then
                      volumean=boxndm%volu
                      amodmean=amod ; bmodmean=bmod ; cmodmean=cmod
                      tcamean=tca; tabmean=tab; tbcmean=tbc
                   else
                      volumean=(volumean*(it/itetemp-1)+boxndm%volu)/(it/itetemp)
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
                   write(6,'(I10,G10.3,A,2G21.12)') it,timel,'*volume  ',boxndm%volu*1d24,volumean*1d24

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



             
136          format(A,3I4,3E15.5,4E15.7,I4)


          endif

       endif
    endif








    ! ecriture de rasmol

    if (iterasmol>0) then     
       if (mod(it,iterasmol)==0) then
         call rasmolT(atdml,boxndm,it)
          if (l2T) call  eleccellmol

       end if
    endif



    !     write(6,*)'sortie ssprogramme analyse'
    !crcmit      call flush(6)

!!$
!!$    if (lprteattotm.EQV..true.) then
!!$       if (ncalceattotm==0)eatomtotm(:)=0.
!!$       ncalceattotm=ncalceattotm+1
!!$       do i=1,im
!!$          eatomtotm(i)=(eatom(i)+(ncalceattotm-1)*eatomtotm(i))/ncalceattotm
!!$       end do
!!$       !     write(912,*)ncalceattotm,eatomtotm(1)*erg2eV,eatom(1)*erg2eV
!!$    end if

    !  if (lbulle) then
    !     if (parallele) then
    !        write(6,*)' test pos He pas para'
    !        stop
    !     end if
    !     call test_position_He(xp,at,ityp,rang,imm,im,it,ldesinteg,num_at_glob,nstepdes,itmax)
    !  end if
    !  write(6,*) 'sortie canalyse',it,im

!!$    if (mod(it,itesigma)==0) then
!!$       if (lsigatcel)then
!!$          natchk(:)=0
!!$          sigatcel=0 ; patcel=0 ; patcelmax=0
!!$          do i=1,im
!!$             koo = ielat(i)                          ! Numero de la cellule
!!$             iti = ityp(i)
!!$             natchk(koo)=natchk(koo)+1
!!$             sigatcel(:,:,koo)=sigatcel(:,:,koo)+sigat(:,:,i)
!!$             ptest=0
!!$             do ic=1,3
!!$                ptest=ptest+sigat(ic,ic,i)/3
!!$             end do
!!$             ptest=abs(ptest)
!!$             patcelmax(koo)=max(patcelmax(koo),ptest)
!!$          end do
!!$          do koo=1,noxyz
!!$             sigatcel(:,:,koo)=sigatcel(:,:,koo)/natchk(koo)
!!$             do ic=1,3
!!$                patcel(koo)= patcel(koo)+sigatcel(ic,ic,koo)/3.0
!!$             end do
!!$             if (natchk(koo).ne.nato(koo)) then
!!$                write(6,*)'nato ? koo natcchk nato', koo, natchk(koo), nato (koo)
!!$                stop
!!$
!!$             end if
!!$          end do
!!$
!!$
!!$          !        write (6, *) '------valeurs par cellules-------',it
!!$          write(extension,'(i9.9)') it
!!$          lenfn2 = 9
!!$          open(luvisuc, file=fnam(1:lenfnam)//'.'//extension(1:lenfn2)//'.CEL.mol', form='formatted', &
!!$               status='unknown')
!!$          write (luvisuc, '(I9,A,I7,A,D15.6)') noxyz , ' IT =', it, ' Time = ', timel
!!$          at=at*1.d8
!!$          write (luvisuc,'(9F12.6)')at(1,1),at(2,1),at(3,1),at(1,2),at(2,2),at(3,2),at(1,3),at(2,3),at(3,3)
!!$          at=at/1.d8
!!$          do kx=0,nox-1
!!$             do ky=0,noy-1
!!$                do kz=0,noz-1
!!$                   koo=1+kx+nox*(ky+noy*kz)
!!$                   xb(1)=float(kx)/float(nox)*at(1,1)+float(ky)/float(noy)*at(1,2)+float(kz)/float(noz)*at(1,3)
!!$                   xb(2)=float(kx)/float(nox)*at(2,1)+float(ky)/float(noy)*at(2,2)+float(kz)/float(noz)*at(2,3)
!!$                   xb(3)=float(kx)/float(nox)*at(3,1)+float(ky)/float(noy)*at(3,2)+float(kz)/float(noz)*at(3,3)
!!$                   xb=xb*1d8
!!$                   write (luvisuc, 149) 'Au',xb(1), xb(2),xb(3),patcel(koo)*unitP,patcelmax(koo)*unitP,nato(koo)
!!$                end do
!!$             end do
!!$          end do
!!$          close (luvisuc)
!!$148       format(A,3I4,3E15.5,1E15.7,I4)
!!$149       format(A,3E15.5,2E15.7,I4)
!!$       end if
!!$    end if
    return
  end subroutine analyseT
end module analyseT_mod
