subroutine divid (appel)
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use var_pot
  use tab_imm_m
  implicit none
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  integer , intent(in) :: appel ! 0: appel partiel juste pour calcul nox/y/z
  ! 1: appel de la routine complete
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  real(double) :: celmin, zlmin, zlm2,rus
  integer::izonr2,izonr,ic
  real(double)::voluperat,rm2,qtot,rut

  !-----------------------------------------------
  !   E x t e r n a l   F u n c t i o n s
  !-----------------------------------------------
  real(double) , external :: distmin, calcvol
  !-----------------------------------------------

  if ((rang==0).and.(appel==0)) then
     write(6,*)
     write(6,*)' -------------------------------------------------------------------'
     write(6,*)'             definition des rayons de coupure'
  endif

  call param_det


  !  write(6,*)rumax,rue_pair,maxval(rue_pair)
  rumax = max(rumax,maxval(rue_pair))
  csive=rumax/float(ngrid)



  if((rang==0).and.(appel==0))      write(6,'(A,F12.2)') 'DIVID rumax',rumax*1d8

  !     if(l3c) then
  if (r3cm.eq.0) r3cm=5.0d-8
  itab=1
  r3cm2=r3cm**2
  !     endif
  !if((rang==0).and.(appel==0))         write (6, *) ' rvois  ', rvois

  if (ltabvois) then
     if (rumax>rvois) then
        if((rang==0).and.(appel==0)) write (6, '(A,2F12.2)') ' rvois trop petit rvois rumax ', rvois*1d8, rumax*1d8
        !cosdebug call arret_ndm
     else
        if((rang==0).and.(appel==0)) write (6,'(A,2F12.2)') ' rumax devient rvois&
             & pour le dimmensionnement en cel', rvois*1d8, rumax*1d8
        rumax=rvois
     end if


  endif
  ! MPI
  !      if ((rang==0).and.(appel==0)) write (6, *) '---- def rayon coupure ok -----'
  !calcul de zlmin
  zlmin = distmin(at(1,1),at(1,2))
  zlm2 = distmin(at(1,1),at(1,3))
  zlmin = min(zlmin,zlm2)
  zlm2 = distmin(at(1,2),at(1,3))
  zlmin = min(zlmin,zlm2)
  zlmin=zlmin*2
  rut=rumax
  if (lpotentiel(10).eqv..true.)      rut=max(rut,2*rue_pot(10))
  !     write(6,*)'BIP',rumax,rut,rue_pot(10)
  !  end if
  if (lpotentiel(11).eqv..true.) rut=max(rut,2*rue_pot(11))
  if (lpotentiel(12).eqv..true.) rut=max(rut,2*rue_pot(12))

  izonr = int(zlmin/rut)
  ! MPI
  if ((rang==0).and.(appel==0)) write (6, *) 'izonr,zlmin,rut', izonr, zlmin*1d8, rut*1d8

     if (izonr<2) then
        write (6, *) 'trop petite boite !!!'
        !cosboite  stop
#if(PHONDY || PARAPH || MAB || ML || PARAML)
        write (6, *) 'trop petite boite !!!'
#else
  if (lrctest) then
        write (6, *) 'STOP ; supprimer avec lrctest=.false. dans din'
     stop
  endif
#endif

  end if
  ! calcul du volume
  !      if ((rang==0).and.(appel==0)) write (6, *) 'avant volu'

  volu = calcvol(at(1,1),at(1,2),at(1,3))

  if ((rang==0).and.(appel==0)) then 
     write (6, '(A,D15.8,A,D15.8,A)') 'volume=', volu,' cm3 ',volu*1d24,' Ang3'
  end if


  !DETERMINATION DE NOX NOY NOZ     
#if(PHONDY || PARAPH || MAB || ML || PARAML)
#else 

  if ((rang==0).and.(appel==0)) write (6, *) 'nox,noy,noz dans .din =', nox, noy, noz
  if ((rang==0).and.(appel==1)) write (6, *) 'nox,noy,noz deuxième passage  =', nox, noy, noz
#endif 

  ! ==== MODIF CLOUET 1 ====================
  !  nzl(1:3) doivent etre calcules ici: ils etaient calcules apres l'appel a divid
  !  dans config, ce qui n'etait pas correct
  call recips (at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3))
  do ic=1,3
     normat(ic)=sqrt(sum(bg(:,ic)**2))
     nzl(ic)=1.0/normat(ic)
  enddo
  ! ==== FIN MODIF CLOUET 1 ================

  if (nox<=0.or.noy<=0.or.noz<=0) then
     ! détermination de nox noy noz qui ne sont pas donnes dans .din
     ! 
     if ((rang==0).and.(appel==0)) write (6, *) 'calcul de nox noy noz !!!'

     ! ==== MODIF CLOUET 2 ====================
     if (izonr<3) then
#ifdef PARA
        write(6,*)'trop petite boite pour para'
        call arret_ndm
#endif
        if ((rang==0).and.(appel==0)) then
           WRITE(6,'(a)') "Boite trop petite: le nombre de cellules est fixe a son minimum"
        endif
     endif
     nox = int(nzl(1)/rumax)
     noy = int(nzl(2)/rumax)
     noz = int(nzl(3)/rumax)
     if ((rang==0).and.(appel==0)) THEN
        write (6,'(a)') 'nox noy noz calcules a partir de ru'
        WRITE(6,'(2(a,g12.4),a,i0)') '  nox = Int( ', nzl(1),'/',rumax,') = ', nox
        WRITE(6,'(2(a,g12.4),a,i0)') '  noy = Int( ', nzl(2),'/',rumax,') = ', noy
        WRITE(6,'(2(a,g12.4),a,i0)') '  noz = Int( ', nzl(3),'/',rumax,') = ', noz
     END IF

     IF (nox.LT.3) nox=3
     IF (noy.LT.3) noy=3
     IF (noz.LT.3) noz=3

     IF ( (nox.LE.3).AND.(noy.LE.3).AND.(noz.LE.3) ) THEN
        nox=1 ; noy=1 ; noz=1
        !             ltabvois=.TRUE.
        !             lconstrtot=.TRUE.
        if ((rang==0).and.(appel==0)) write(6,*)'!!!!!!!!!!Envisager ltabvois = true !!!!!!!!!!!!!!'
     END IF

     if ((rang==0).and.(appel==0)) write (6,'(a,3(i0,1x))') 'nox noy noz apres correction = '&
          , nox, noy, noz
!!$     if (izonr<3) then
!!$
!!$        ltabvois = .TRUE.
!!$        lconstrtot=.TRUE.
!!$#ifdef PARA
!!$        write(6,*)'trop petite boite pour para'
!!$        call arret_ndm
!!$
!!$#endif
!!$        nox = 1 ; noy = 1 ; noz = 1
!!$        if ((rang==0).and.(appel==0)) then
!!$           write (6, *) 'plus de cellules, on considere toute la boite '
!!$        endif
!!$     else
!!$        nox = int(nzl(1)/rumax)
!!$        noy = int(nzl(2)/rumax)
!!$        noz = int(nzl(3)/rumax)
!!$        if ((rang==0).and.(appel==0)) write (6, *) 'nox noy noz calcules a partir de ru = '&
!!$             , nox, noy, noz
!!$     endif
     ! ==== FIN MODIF CLOUET 2 ================
     celsize(1) = zl(1)/float(nox)
     celsize(2) = zl(2)/float(noy)
     celsize(3) = zl(3)/float(noz)

  else

     ! *** nox noy noz sont donnes dans.din ***

     if (nox==2.or.noy==2.or.noz==2) then
        write (6, *) rang,'wrong noxyz stop'
        call arret_ndm
     endif
     !     if (nox==1.or.noy==1.or.noz==1) then
     !        if ((rang==0).and.(appel==0))  write (6, *) 'plus de cellules, on considere toute la boite'            
     !        ltabvois = .TRUE.
     !        lconstrtot=.TRUE.
     !        nox = 1
     !        noy = 1
     !        noz = 1
     !     endif

     celsize(1) = zl(1)/float(nox)
     celsize(2) = zl(2)/float(noy)
     celsize(3) = zl(3)/float(noz)
     celmin = min(celsize(1),celsize(2))
     celmin = min(celsize(3),celmin)
     ! ==== MODIF CLOUET 3 ================================
     ! Je ne comprends pas l'interet de ce test qui vient detruire ce qui a ete
     ! fait auparavant. Ne peut-on pas faire confiance aux valeurs fournies par
     ! l'utilisateur concernant nox, noy et noz. Si ce n'est pas possible, il
     ! faudrait recopier le code correspondant a MODIF CLOUET 2
!!$     if (celmin<=rumax) then
!!$        write (6, *) 'changement de nox noy noz !!!'
!!$
!!$        if (izonr<3) then
!!$           ltabvois = .TRUE.
!!$           lconstrtot=.TRUE.
!!$           nox = 1
!!$           noy = 1
!!$           noz = 1
!!$           if ((rang==0).and.(appel==0)) then
!!$              write (6, *) 'plus de cellules, on considere toute la boite'
!!$           endif
!!$        else
!!$           nox = int(nzl(1)/rumax)
!!$           noy = int(nzl(2)/rumax)
!!$           noz = int(nzl(3)/rumax)
!!$           if ((rang==0).and.(appel==0)) write (6, *) &
!!$                'nox noy noz calcules a partir de ru = ', nox, noy, noz
!!$        endif
!!$        celsize(1) = zl(1)/float(nox)
!!$        celsize(2) = zl(2)/float(noy)
!!$        celsize(3) = zl(3)/float(noz)
!!$        !        write(6,*)'celsize ok '
!!$     endif
     ! ==== FIN MODIF CLOUET 3 ============================

  endif
  ! nox noy et noz sont determines

  noxy = nox*noy
  noxyz = nox*noy*noz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (appel==0) then
     if (rang==0) then
        write(6,*)'retour à la construction de la boite'
        write(6,*)
     end  if
     return
  end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !RETURN!RETURN!RETURN!RETURN!RETURN!RETURN!RETURN!RETURN!RETURN!RETURN!RETURN apppel=0


  IF (natperc.LE.0) THEN        ! MODIF Clouet
     natperc= INT(im_glob/noxyz)
     nvat=10*natperc
     !          if(natperc.le.2)then
     !             if (rang==0) &
     !                  write(6,*) 'moins de 3 atomes par celulle -> table des voisins complete' 
     !             ltabvois=.TRUE. ; lconstrtot=.TRUE.
     !          end if
!!$natperc=max(2*natperc,10)     ! MODIF Clouet
     natperc=max(5*natperc,10)     ! MODIF Clouet
  ELSE                          ! MODIF Clouet
     nvat=10*natperc       ! MODIF Clouet
  END IF                        ! MODIF Clouet



  !  natperc= INT(im_glob/noxyz)
  !  nvat=10*natperc

  natperc=max(5*natperc,20)

#if(PHONDY || PARAPH || MAB || ML || PARAML)
#else 
  if (rang==0) &
       write(6,*) 'natperc im/noxyz', natperc, im_glob/noxyz

  if (rang==0)  write(6,*)'ltabvois,lconstrtot',ltabvois,lconstrtot
#endif 

  if (ltabvois) then
     if (rumax>rvois) then
        write (6, *) rang,' rvois trop petit rvois rumax ', rvois, rumax
        call arret_ndm
     endif
     !crc        rm2=max(rumax,2*rvois)
     rm2=max(rumax,rvois)
     izonr2 = int(zlmin/rm2)
     if (izonr2<1) then
        write (6, *) rang,'trop petite boite pour rvois !!!'
        !cosboite   call arret_ndm
#if(PHONDY || PARAPH || MAB || ML || PARAML)

        write (6, *) rang,'trop petite boite pour rvois !!!'
#else 
        call arret_ndm
#endif
     endif
     if(.not.lconstrtot)rumax=rvois
     voluperat=volu/im
     IF (nvperat.LE.0) nvperat=4*Pi*(rvois+1.0d-8)**3/(3*voluperat)

#if(PHONDY || PARAPH || MAB || ML || PARAML)

#else
     if(rang==0)         write (6, '(A,D10.3)') 'volumeperat=', voluperat
     if(rang==0)         write (6, '(A,D10.3)') 'Rvois=', RVois
     if(rang==0)         write (6, *) 'NVperat= ', nvperat
#endif 
     if (ldemitab) then
        nvois=max(Int(1.5*nvperat*im),100)
        nvat=max(Int(nvperat*1.3),10)
     else
        nvois=max(Int(1.5*nvperat*im),100)
        nvat=max(Int(nvperat*1.3),10)
     end if
#if(PHONDY || PARAPH || MAB || ML || PARAML)

#else

     if(rang==0)         write (6, *) 'Nvois= ', nvois
#endif 
     allocate(indi(nvois))
     allocate(indi2(nvois))
  else

  end if


  qtot = 0
  qtot = sum(q(:ntyp)*na(:ntyp))
  if ((qtot/=0.0).and.(rang==0)) write (6, *) ' CHARGE NON NULLE !! = ', qtot

  izonr = int(zlmin/r3cm)
#if(PHONDY || PARAPH || MAB || ML || PARAML)

#else
  if(rang==0) write(6,*)
  if(rang==0) write(6,*)' TABLEAUX DIMENSIONES POUR UNE BOITE UNIFORME !! ' 
  if(rang==0) write(6,*) '-------------------------------------------------------------------'
  if(rang==0) write(6,*)
#endif

  return
end subroutine divid

