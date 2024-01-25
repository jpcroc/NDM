module analyseT_mod
  USE Mat_utils_mod
  USE calctemp_mod,only: calctemp
  USE calccoordo_mod,only: calccoordo
  USE calcdigr_mod,only: calcdigr,initrdf,rdfT,rdf0
  USE calcangle_mod,only: calcangle,adf0,initadf,adfT
  USE bondval_mod,only: bondval
  USE rasmolT_mod,only: rasmolT
  use calcextr_mod,only:calfoextr
  USE sauvegardeT_mod,only:sauvegardeT
  use notperiod_mod,only:notperiod
  use var_pot, only: iewald,l3c,npotmax,potisglue,potisrep,lpotentiel,ntyp,nkmax,contmax,zz,potis1,&
       &cm
  use gen_com_m, only:bk,cunite,fnose,iteanapos,iteangle,itebdv,ecellpr,&
       &itecoordo,iterasmol,iterdf,iteprtsigma,itetemp,itetemp2,kcell,kine,kinemean,knose,&
       &leev,leparat,linstantfda,lprahman,lprteattotm,lsigatcel,lthoover,ltnose,ltpcel,lucell,&
       &nfda,pist,pmean,potcp,potis2,potis3,potist,potistersoff,potiszbl,thetamin,thetamax,&
       &tcou,temp,tempep,tfcou,tmean,ucell,unite,unose,zhoover,sig,sigkine,lprtcel,rcangle,&
       &tpseuils,sigtot,unitP,nrdf,lprtsigat,lprteat,lpkbar,linstantrdf,linstantfda,&
       &itloopmax,cunitp,erg2ev,lperiod,pi,rang,timel,latcomp,h0,rcrdf,iteangle,itedepla,tdepla,tdepla2,&
       & itesauvforce,itesauv,fnamcout,itesauvinter,itesauvposition,fnam,lenfnam,iteration,l2T,iteprtkin

  USE cellconfig,only:cell_config, caltabtC,cell_config_arps
  USE atomconfig,only:atom_config,atom_config_d,atom_config_e
  use boxconfig,only: box_config
  use Tpara,only:nprocspace,myidsp
  use calcdepla_mod,only:calcdepla
  use newunit_mod,only:newunit
  implicit none

contains
  ! ************************************************
  !         Sous-programme analyse.f
  ! ************************************************

  subroutine analyseT(atdml,celndm,boxndm,psc)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double

    use elec_cell, only : Eelec,Teavg,Tecmax,ietm,eleccellmol
    use eloss, only : ibrake, elosselectot1, elosselectot
    use Tpara,only:para_space_config    
    implicit none


    type(para_space_config)::psc    
    !-----------------------------------------------
    !   G l o b a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s        
    !-----------------------------------------------
    integer :: i, iti, ic,formatsauv
    real(double), dimension(ntyp) :: temptyp


    class(atom_config_d)::atdml
    class(cell_config):: celndm
    class(box_config)::boxndm

    type(atom_config_d)::attyp
    type(cell_config):: celtyp

    real(double) :: ppot, pkin
    real(double) :: a1, a2, a3, b1, b2, b3, c1, c2, c3
    real(double) :: fteta, tbc, tca, tab, amod, bmod, cmod,kinetyp

    real(double), save :: volumean,amodmean,bmodmean,cmodmean,tcamean,tabmean,tbcmean
    real(double), dimension(3,3) :: transformation, strain, rotation, invh0
    integer::ipot
    character :: extension*9
    real(double)::xb(3)
    integer,save::unitkin
    real(double)::atkin
    logical,save::linitrdf=.false.,linitadf=.false.,lopenkin=.false.


    if (itloopmax==0) itetemp=0
    if (itesauv.GT.0) then
       if (mod(iteration,itesauv)==0) then 
          formatsauv = 3
          if(itesauvinter.gt.0) then
             if (mod(iteration,itesauvinter).eq.0) then
                write(extension,'(i9.9)') iteration
                fnamcout = fnam(1:lenfnam)//'.cout.'//extension
             else
                fnamcout = fnam(1:lenfnam)//'.cout'
             endif
          else
             fnamcout = fnam(1:lenfnam)//'.cout'
          end if
          call sauvegardeT(atdml,celndm,boxndm,formatsauv,fnamcout,latcomp=latcomp)
       end if
    endif
    select type (atdml)
    type is (atom_config_e)
       if (itedepla.GT.0) then
          if (mod(iteration,itedepla)==0) then 
             call calcdepla(atdml,celndm,boxndm,tdepla,iteration)
          end if
       end if
       if (itedepla.GT.0) then
          if ((mod(iteration,itedepla)==0).and.(tdepla2.gt.0)) then 
             call calcdepla(atdml,celndm,boxndm,tdepla2,iteration,'2')
          end if
       end if
    end select
    
    if (itesauvposition.GT.0) then
       if (mod(iteration,itesauvposition)==0) then
          formatsauv = 2
          write(extension,'(i9.9)') iteration
          fnamcout = fnam(1:lenfnam)//'.cout.'//extension
          call sauvegardeT(atdml,celndm,boxndm,formatsauv,fnamcout,latcomp=latcomp)
       end if
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
    if (itetemp>0) then
       if (mod(iteration,itetemp)==0) then
          !          call atdml%print
          call calctemp (temp,kine,atdml,celndm)
          block
            logical:: lgs(atdml%imm)
            do iti=1,ntyp
               lgs=atdml%lgul
               atdml%lgul=.false.
               celtyp=celndm
               where(atdml%ityp(1:atdml%im)==iti)
                  atdml%lgul(1:atdml%im)=.true.
               end where
               !             if (ALL(atdml%lgul(1:atdml%im).eqv..false.)) cycle
               call atdml%fab(attyp,lback=.false.)
               call caltabtC(celtyp,attyp,lperiod,boxndm)
               call calctemp(temptyp(iti),kinetyp,attyp,celtyp)
               call celtyp%dealloc ; call attyp%dealloc
               atdml%lgul=lgs
             ! ceci est un test du calcul des forces sur un sous-ensemble des atomes
             !ne fonctionne que pour les pots de paires
!!$             if (iti==2)then
!!$                do i=1,atdml%im
!!$                   if (atdml%ityp(i)==2) write(344,*)i,atdml%xp(:,i)
!!$                   if (atdml%ityp(i)==2) write(344,*)i,atdml%fp(:,i)
!!$                end do
!!$                call calfoextr(atdml,celndm,boxndm,psc)
!!$             end if
            end do
          end block
          ! MPI
          !remarque 1erg = 6.24d11 eV
          if (mod(iteration,itetemp2)==0) then
             if (rang==0) then

                if (ibrake.GT.0)  write(6,*)'electronic losses ', elosselectot, elosselectot1

                write (6, *)
                write (6, *)
                write (6, '(A,I7,A,G15.7)') '<<<<<<<<  ITERATION =', iteration, &
                     '  time = ', timel
                write (6, *)
                write (6, *)
                write (6, *) '--------- Temperatures et Pression--------'
                write (6, *) '----------valeurs instantanees------------'


                write (6, '(I10,G10.3,A,G25.16,A)') iteration, timel, '*Epot = ', potist*unitE, cunitE
                do ipot=1,npotmax
                   if (lpotentiel(ipot).eqv..true.) then
                      select case (ipot)
                      case (1:9)
                         write(6,'(A,G21.12,A)')'    *energie paire 2 corps = ',(potis1)*unitE, cunitE
                         if(l3c) write(6,'(A,G21.12,A)')'    *energie pot 3 corps = ',potcp*unitE, cunitE
                      case(10:12)
                         write(6,'(A,G21.12,A)')'    *energie PAIRE EAM = ',potisrep*unitE, cunitE
                         write(6,'(A,G21.12,A)')'    *energie GLUE  = ',potisglue*unitE, cunitE
                      case(16)
                         write(6,'(A,G21.12,A)')'    *energie PAIRE EAM = ',potisrep*unitE, cunitE
                         write(6,'(A,G21.12,A)')'    *energie GLUE  = ',potisglue*unitE, cunitE
                      case(13)
                         write(6,'(A,G21.12,A)')'    *energie Tersoff = ',potisTersoff*unitE, cunitE
                         If (potisZBL.ne.0)write(6,'(A,G21.12,A)')'    *energie ZBL = ',potisZBL*unitE, cunitE
                      end select
                   end if
                end do
                if( potis2.ne.0) write (6, '(A,G21.12,A)') '    *energie EWALD 2eme terme = ', (potis2)*unitE, cunitE
                if ((iewald.gt.0).and.(iewald.ne.3)) write (6, '(A,G21.12,A)') '    *energie EWALD RECIP = ', (potis3)*unitE, cunitE
                write (6,'(I10,G10.3,A,G21.12,A,a,f0.3,a)') iteration,timel,'*Ec = ',kine*unitE, cunitE, &
                     '  (', temp, ' K)'

                write (6,'(I10,G10.3,A,G25.16,A)') iteration,timel,'*Etot = ',(kine+potist)*unitE, cunitE
                If (l2T) then

                   write (6,'(I10,G10.3,A,G21.12,A)') iteration,timel,'*Eelec = ',Eelec*unitE, cunitE                 
                   write (6,'(I10,G10.3,A,G21.12,A)') iteration,timel,'*IE_Et = ',&
                        &(kine+potist+Eelec)*unitE, cunitE                 
                   write (6,'(I10,G10.3,A,G21.12,A)') iteration,timel,'*Telec = ',Teavg
                   write (6,'(I10,G10.3,A,G21.12,A)') iteration,timel,'*TempEP = ',TempEP
                   write (6,'(I10,G10.3,A,G21.12,3I5)') iteration,timel,'*maxTe = ',Tecmax,ietm(:)

                end If

                write (6, *)

                !	IF(it==1) Open(unit=774, file='energie_it.dat', status='unknown', action='write')
                !	IF(modulo(it,10)==0) WRITE(774,'(i6, f)') it, (kine+potist)*unitE  ! controle .. !*!

                !*! ---------------------------------- APPEL DE LA ROUTINE D'ECRITURE -------------
                ! ------------------------------------ DES FICHIERS DE SORTIE, DANS LE ------------
                ! ------------------------------------ CAS DES CL CONTROLEES EN CONTRAINTE --------


                if (lEparat) then
                   write (6, '(I10,G10.3,A,G21.12,A)') iteration, timel, '*Epot/at = ', potist*unitE/atdml%im_glob, cunitE
                   write(6,*)
                end if

                write (6, '(I10,G10.3,A,f0.3)') iteration,timel, &
                     '*Temp instantanee = ',temp

                if (tfcou>0.0) write (6, '(A,G15.4)') '*temperature externe = ', tcou

                if (celndm%ltpcel) then
                   call plottpcel(celndm)
                end if
                IF (lprahman) THEN
                   write(6,*) 'NPT With Parrinello-Rahman'
                   write(6,'(I10,G10.3,A,G21.12,A,a,f0.3,a)') iteration,timel,'*Kcell = ',Kcell*unitE,cunitE, &
                        '  (', 2.d0*Kcell/(9.d0*bk), ' K)'
                   IF (lUcell) THEN
                      write(6,'(I10,G10.3,A,G21.12,A)') iteration,timel,'*Ucell = ',Ucell*unitE,cunitE
                      write(6,'(I10,G10.3,A,G21.12,A)') iteration,timel,'*Ecell = ',EcellPR*unitE,cunitE
                      write(6,'(I10,G10.3,A,G21.12,A)') iteration,timel,'*Htot_PR = ', &
                           (potist+kine+EcellPR)*unitE,cunitE
                   END IF
                   write(6,*) 'Box tensor'
                   write(6,*)'a',boxndm%at(1,1),boxndm%at(2,1),boxndm%at(3,1)
                   write(6,*)'b',boxndm%at(1,2),boxndm%at(2,2),boxndm%at(3,2)
                   write(6,*)'c',boxndm%at(1,3),boxndm%at(2,3),boxndm%at(3,3)
                   write(6,*)'H0(1)',h0(1,1),h0(2,1),h0(3,1)
                   write(6,*)'H0(2)',h0(1,2),h0(2,2),h0(3,2)
                   write(6,*)'H0(3)',h0(1,3),h0(2,3),h0(3,3)
                   Call MatInv(h0, invh0)
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
                   if(iteration<=1) then
                      volumean=boxndm%volu
                      amodmean=amod ; bmodmean=bmod ; cmodmean=cmod
                      tcamean=tca; tabmean=tab; tbcmean=tbc
                   else
                      volumean=(volumean*(iteration/itetemp-1)+boxndm%volu)/(iteration/itetemp)
                      amodmean=(amodmean*(iteration/itetemp-1)+amod)/(iteration/itetemp)
                      bmodmean=(bmodmean*(iteration/itetemp-1)+bmod)/(iteration/itetemp)
                      cmodmean=(cmodmean*(iteration/itetemp-1)+cmod)/(iteration/itetemp)
                      tcamean=(tcamean*(iteration/itetemp-1)+tca)/(iteration/itetemp)
                      tbcmean=(tbcmean*(iteration/itetemp-1)+tbc)/(iteration/itetemp)
                      tabmean=(tabmean*(iteration/itetemp-1)+tab)/(iteration/itetemp)
                   end if
                   write(6,'(I10,G10.3,A,2G18.9)') iteration,timel,'*l_a,m  ',amod,amodmean
                   write(6,'(I10,G10.3,A,2G18.9)') iteration,timel,'*l_b,m  ',bmod,bmodmean
                   write(6,'(I10,G10.3,A,2G18.9)') iteration,timel,'*l_c,m  ',cmod,cmodmean
                   write(6,'(I10,G10.3,A,2F11.4)') iteration,timel,'*ang_bc,m  ',tbc,tbcmean
                   write(6,'(I10,G10.3,A,2F11.4)') iteration,timel,'*ang_ca,m  ',tca,tcamean
                   write(6,'(I10,G10.3,A,2F11.4)') iteration,timel,'*ang_ab,m  ',tab,tabmean
                   write(6,'(I10,G10.3,A,2G21.12)') iteration,timel,'*volume  ',boxndm%volu*1d24,volumean*1d24

                endif    ! if (lprahman)

                IF (lTNose) THEN
                   WRITE(6,'(a)') 'Thermostat de Nose'
                   WRITE(6,'(I10,G10.3,A,G21.12,A,a,f0.3,a)') iteration,timel, &
                        '*KNose = ',KNose*unitE,cunitE, '  (', 2.d0*KNose/(bk), ' K)'
                   WRITE(6,'(i7,G10.3,a,g22.12)') iteration,timel,'*fNose = ', fNose
                   WRITE(6,'(i7,G10.3,a,g22.12,a)') iteration,timel,'*Htot_Nose = ', (kine+potist+EcellPR+KNose+UNose)*unitE,cunitE

                ELSEIF (lTHoover) THEN
                   WRITE(6,'(a)') 'Thermostat de Nose-Hoover'
                   WRITE(6,'(I10,G10.3,A,G21.12,A,a,f0.3,a)') iteration,timel, &
                        '*KNose = ',KNose*unitE,cunitE, '  (', 2.d0*KNose/(bk), ' K)'
                   WRITE(6,'(i7,G10.3,a,g22.12)') iteration,timel,'*zHoover(1) = ', zHoover(1)
                   !WRITE(6,'(i7,G10.3,a,g22.12,a)') iteration,timel,'*Htot_Hoover = ', (kine+potist+Ecell+KNose+UNose)*unitE,cunitE

                END IF   ! if lTNose / lTHoover

             endif                                ! rang=0
          end if
          IF (iteration<=1) THEN
             tmean = temp
             kinemean = kine
          ELSE
             tmean = (tmean*(iteration/itetemp-1)+temp)/(iteration/itetemp)
             kinemean = (kinemean*(iteration/itetemp-1)+kine)/(iteration/itetemp)
          END IF


          if (rang==0) then
             write (6, *)
             do iti = 1, ntyp
                if (count(atdml%ityp==iti)==0) cycle
                if (mod(iteration,itetemp2)==0) then
                   write (6, '(A,I2,A,F12.2)') &
                        '*temp instantanee des atomes de type', iti, ' = ', &
                        temptyp(iti)
                end if
             end do

             if (iteprtSigma>0) then
                if (mod(iteration,iteprtsigma)==0) then
                   if (mod(iteration,itetemp2)==0) then
                      write (6, *)
                      write (6, *) '* stress en ', cunitP
                   end if
                   ppot = 0.0
                   pkin = 0.0
                   do ic = 1, 3
                      if (mod(iteration,itetemp2)==0) then
                         write (6, '(I1,3(A,I1),A,3G18.10)') ic,' sigma potentiel (1,', ic, ') (2,', ic, &
                              ') (3,', ic, ') =',sig(1:3,ic)*unitP
                      end if
                      ppot = ppot+1.0/3.0*sig(ic,ic)
                   end do
                   if (mod(iteration,itetemp2)==0) then
                      write (6, *)
                   end if
                   do ic = 1, 3
                      if (mod(iteration,itetemp2)==0) then
                         write (6, '(I1,3(A,I1),A,3G18.10)') ic,' sigma cinetique (1,', ic, ') (2,', ic, &
                              ') (3,', ic, ') =',sigkine(1:3,ic)*unitP
                      end if
                      pkin = pkin+1.0/3.0*sigkine(ic,ic)
                   end do
                   if (mod(iteration,itetemp2)==0) then
                      write (6, *)
                      do ic = 1, 3
                         write (6, '(I1,3(A,I1),A,3G18.10)') ic,' sigma total (1,', ic, ') (2,', ic, &
                              ') (3,', ic, ') =',sigtot(1:3,ic)*unitP
                      end do

                      write (6, *)
                   end if


                   pist = ppot+pkin
                   if (mod(iteration,itetemp2)==0) then
                      write (6, '(I10,G10.3,A,G14.5)') iteration, timel, '*Pression = '&
                           , pist*unitP
                   end if
                   if(iteration<=1) then
                      pmean = pist
                   else
                      pmean = (pmean*(iteration/iteprtsigma-1)+pist)/(iteration/iteprtsigma)
                   end if
                   if (mod(iteration,itetemp2)==0) then
                      write (6, *) 'pression totale  potentiel = cinetique ='
                      write (6, '(3g14.5)') pist*unitP, ppot*unitP, pkin*unitP
                   end if




                endif
             endif                             !  itesigma>0


             if (mod(iteration,itetemp2)==0) then
                !               write (6, *)
                write (6, *) '----------valeurs moyennes------------'
                write (6, '(A,F12.2)') '*temperature moyenne = ', tmean
                !               write (6, '(A,G21.12,A)') '*energie cinetique moyenne = ', &
                !                  kinemean*unitE, cunitE
                if (iteprtsigma>0) then
                   if (mod(iteration,iteprtsigma)==0) write (6, '(A,G14.5)') &
                        '*pression moyenne = ', pmean*unitP
                endif
                write (6, *) '--------------------------------------'
             end if




136          format(A,3I4,3E15.5,4E15.7,I4)


          endif

       endif
    endif

!!$    select type (atdml)
!!$    type is (atom_config_e)
!!$       write(6,*)'prteat'
!!$       if (lprteat) then
!!$          do i=1,atdml%im
!!$             write(6,*)atdml%eat(i)
!!$          end do
!!$       end if
!!$    end select

    ! ecriture de rasmol

    if (iterasmol>0) then     
       if (mod(iteration,iterasmol)==0) then
          call rasmolT(atdml,boxndm,iteration,latcomp=latcomp)
          if (l2T) call  eleccellmol

       end if
    endif



    if (iterdf>0) then
       if (linitrdf.eqv..false.) then
          call initrdf(rdf0,nkmax,rcrdf,linstantrdf,'00')
          linitrdf=.true.
       end if

       if (mod(iteration,iterdf)==0) then
          call calcdigr (atdml,celndm,boxndm,rdf0)
          if (rdf0%linstantrdf) then
             call rdfT(rdf0)
             rdf0%nrdf = 0
          endif
       endif
    else if (iterdf==0) then
       if (itloopmax-iteration<nrdf) then
          call calcdigr (atdml,celndm,boxndm,rdf0)
       endif
    endif

    if (iteangle>0) then
       if (linitadf.eqv..false.) then
          call initadf(adf0,contmax,rcangle,linstantfda,'00',thetamax,thetamin)
          linitadf=.true.
       end if

       if (mod(iteration,iteangle)==0) then
          call calcangle(atdml,celndm,boxndm,adf0)
          if (adf0%linstantfda) then
             call adfT(adf0)
             adf0%nfda = 0
          endif
       endif
    else if (iteangle==0) then
       if (itloopmax-iteration<nfda) then
          call calcangle(atdml,celndm,boxndm,adf0)
       endif
    endif
    if (itebdv>0) then
       if (mod(iteration,itebdv)==0) call bondval(atdml,celndm,boxndm)
    end if

    if (iteprtkin>0) then
       if (mod(iteration,iteprtkin)==0) then
          if (.not.lopenkin) then
             call newunit(unitkin)
             open(unitkin, file='prtkin')
             lopenkin=.true.
          end if
          do i=1,atdml%im
             atkin=erg2ev*(atdml%vp(1,i)**2+atdml%vp(2,i)**2+atdml%vp(3,i)**2)*0.5*cm(atdml%ityp(i))
             write(unitkin,*)atkin
          end do
       end if
    end if
    return
  end subroutine analyseT
  subroutine plottpcel(celcf,itapp)
    integer,optional::itapp
    class(cell_config)::celcf
    integer::itp,unitlt,i,koxyz(3),unitlp
    character :: extension*9
    character*80::namef
    real(double)::Pcell
    
    if (.not.celcf%ltpcel) then
       write(6,*)'coding error in plottpcel call, %ltpcel.ne.true'
       call arret_ndm
    end if
    
    if (present(itapp)) then
       itp=itapp
    else
       itp=iteration
    end if
    if (myidsp==0) then
       call newunit(unitlt)
       write(extension,'(i9.9)')itp
       namef='TEMPC.'//trim(extension)
       open(unitlt,file=namef,form='formatted')


       
       do i=1,celcf%noxyz
          koxyz=celcF%koxyz(i)
          select type (celcf)
          type is (cell_config)
             write(unitlt,'(I12,3I5,G15.5,I5)'),i,koxyz(1:3),celcf%tempc(i),celcf%nato(i)
          type is (cell_config_arps)
             write(unitlt,'(I12,3I5,G15.5,4I5)')i,koxyz(1:3),celcf%tempc(i),celcf%nato(i),celcf%nmov(0,i),celcf%nmov(1,i),celcf%nmov(2,i)
          end select
       end do
       close(unitlt)

       call newunit(unitlp)
       write(extension,'(i9.9)')itp
       namef='PRESSC.'//trim(extension)
       open(unitlp,file=namef,form='formatted')
     
       do i=1,celcf%noxyz
          koxyz=celcF%koxyz(i)
          pcell=(celcf%sigc(1,1,i)+celcf%sigc(2,2,i)+celcf%sigc(3,3,i))*unitP
          select type (celcf)
          type is (cell_config)
             write(unitlp,'(I12,3I5,G15.5,I5,9E15.5)'),i,koxyz(1:3),pcell,celcf%nato(i),celcf%sigc(:,:,i)*unitP
          type is (cell_config_arps)
             write(unitlp,'(I12,3I5,G15.5,4I5,9E15.5)')i,koxyz(1:3),pcell,celcf%nato(i),celcf%nmov(0,i),&
                  &celcf%nmov(1,i),celcf%nmov(2,i),celcf%sigc(:,:,i)*unitP
          end select
       end do
       close(unitlt)



       
    end if
  end subroutine plottpcel
end module analyseT_mod
