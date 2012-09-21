
! **************************************************************
subroutine init
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use tab_imm_m
  use eam
  use SMjuli
  use jqmod
  use neb_module
  use posana
  use defcdp, ONLY :itecdp
!  use var_pot
#if(PARA)
  use mod_mpi
#endif 

  !      use eam
  ! **************************************************************

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
  integer :: i, lufilmpaf,itapp,ipotcont
  integer :: complet=1    ! flag d'appel a divid : complet : exec de la routine complete
  !-----------------------------------------------

  tmean = 0.0
  pmean = 0.0
  timel = 0.0
  kinemean = 0.0
  lufilmpaf = 79

  !     write(6,*)'entree dans init.f'
  !potentiel BKS

#if(PARA)
  temps_input_deb = MPI_Wtime()
#endif
  if (npotentiel.gt.1)then
     ipotentiel=-1
     npair=  ntyp*(ntyp+1)/2 ; ntrip= ntyp*ntyp *(ntyp+1)/2
     allocate(typ_and_pot(ntyp,npotmax))
     typ_and_pot(:,:)=.false.
     call  alloc_typ
  end if
  allocate(rue_pot(npotmax))
  rue_pot(:)=0.


  do ipotcont=0,npotmax
     if(lpotentiel(ipotcont).EQV..true.) then 
        ipotentiel=ipotcont
     else
        cycle
     end if

     if (ipotentiel.lt.10) then
        if (rang.eq.0) then
           write(6,*)
           write(6,*)'POTENTIEL DE PAIRES '
           write(6,*)
        end if

        call input_pair
     else
        select case(ipotentiel)
        case(10)
           if (rang.eq.0) then
              write(6,*)
              write(6,*)'POTENTIEL EAM'
              write(6,*)
           end if
           call inputeam(ntyp,npair,ntrip,cm,catom,ty,umass,rue_pot(ipotentiel),rumax,&
                iewald,l3c,rang,r3cm,roff1,roff2,typ_and_pot,npotmax,&
                ipotentiel,typ_pot_pair,lue_typ,lue_paire,lu_roff_pair,&
                npotentiel,ipo)
           do i=1,npair
              if (typ_pot_pair(i)==ipotentiel) rue_pair(i)=rue_pot(ipotentiel)
           end do
        case(12)
           if (rang.eq.0) then
              write(6,*)
              write(6,*)'POTENTIEL Ju Li'
              write(6,*)
           end if
           call inputeamjl(ntyp,npair,ntrip,cm,catom,ty,umass,rue_pot(ipotentiel),&
                rumax,iewald,l3c,rang,r3cm,roff1,roff2,typ_and_pot,&
                npotmax,ipotentiel,typ_pot_pair)
           rue_pair(:)=rue_pot(ipotentiel)
        case(13,14,15)
           !nguyen mettre input tersoff
           !        if (rang.eq.0) then
           if (rang==0)   write(6,*)'POTENTIEL tersoff.potin'
           !           if (ipotentiel==13) then
           !           
           !           else
           !              write(6,*)'POTENTIEL tersoff.potin COUPURE MODIFIEE !!!!!!!!!!!!!!!!!!!!!!!!'
           !           end if
           !           write(6,*)
           !        end if
           if (rang.eq.0) then
              write(6,*)
              write(6,*)'POTENTIEL Tersoff-Brenner'
              write(6,*)
           end if

           call inputtersoff

        end select
     endif
  end do
  write(6,*)'cm',cm
  usdh = 1/(two*tstep)         
  !endif



  it=0 

  !---inNEB
  if (dmtype.ne.9) then
#if(PARA)
     temps_input=MPI_Wtime()-temps_input_deb

     temps_config_deb = MPI_Wtime()
#endif
     call config 
#if(PARA)
     temps_config=MPI_Wtime()-temps_config_deb
#endif


#if(DECOUP)
     ! Pas la peine d'aller plus loin dans l'initialisation
     return
#endif
  else

     call configNEB(xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
  end if
  !...inNEB

  if (iterasmol>=0) call rasmol (it)

  if (igen==(-1)) then
     formatsauv = 2
     call sauvegarde
     if (rang==0) write (6, *) 'generation terminee'
     call arret_ndm
  endif

  if (igen==2) then
     call transf 
     formatsauv = 2
     call sauvegarde
     if (rang==0) write (6, *) 'modification terminee'
     call arret_ndm
  endif

  ! determination des tailles du nombre de cel. (nox, noy, noz)

  call divid (complet)
  call DynamicalAllocationCell
  !      allocate (tabv3(-ncoucx:ncoucx,-ncoucy:ncoucy,-ncoucz:ncoucz))

  do ipotcont=0,npotmax
     if(lpotentiel(ipotcont).EQV..true.) then 
        ipotentiel=ipotcont
     else
        cycle
     end if

     select case(ipotentiel) 
     case(0:9)
        call calpo 
     case(10:12)
        call calpoeam
     case(13,14,15)
        if (.not.parallele) then
           if (maxval(roff1).gt.0) call tersoff_zbl
        end if
     end select
  end do
  if (npotentiel.gt.1) then
     !     write(6,*)
     !     write(6,*)'BILAN DES POTENTIELS'
     !     do i=1,ntyp
     !        do j=1,ntyp
     !           l=ipo(i,j)
     !           write(6,*)
     !           write(6,*)'paire i-j l',i,j,l
     !           write(6,*)'lue paire typ_pot_pair coupure'
     !           write(6,*)lue_paire(l),typ_pot_pair(l),rue_pair(l)*1d8
     !        end do
     !     end do
     if (rang==0) then
        write(6,*)
        write(6,*)'decoupage en cellule suivant'
        write(6,*)'rumax',rumax*1d8
        write(6,*)
     end if
  end if
  call neigcel 

#if(PARA)
  call init_voisinage()

  if (rang==0)  write(6,*) 'NOMBRE DE CELLULES FRONTIERES ASSOCIEES A CHAQUE PROCESSEUR'
  write(6,*) 'Le proc ',myid,' a ',nbr_proc_voisin,' processeur voisin'
  !  do i=1,nbr_proc_voisin
  !     WRITE(6,*) 'Le proc ',myid,' envoit ',nbr_cell_frontiere(i),' vers le proc ',proc_voisin(i)
  !  enddo
  !  WRITE(6,*) 'Le proc ',myid,' recoit ',nbr_cell_ftm,' cell. fantome de ses voisins'
#endif 

  ! !!! compcr non pris en charge en parallele !!!




  imd = im
  if (ltranche) call layer


  nad(:ntyp) = na(:ntyp)
!!!  endif

  imf=im
  imana=min(imd,imd)



  !  if (itmax>0) then
  call caltabt
  ! if (rang==0)  write(6,*)'>>>>>>>>>>>apres caltabt'
  if (ltabvois) call caltabi 
  ! if (rang==0)  write(6,*)'>>>>>>>>>>>apres caltabi'
  !  end if
  if (.not.lrestart) then
     !    if (rang==0)     write(6,*)'>>>>>>>>>>>avant initspeed'
#if(PARA)
     temps_initspeed_deb = MPI_Wtime()
#endif
     call initspeed 
     if (lcorrelvp) then
        ax=vp
        write(6,*)'AX DEVIENT VP0'
        write(6,*)'AX DEVIENT VP0'
        write(6,*)'AX DEVIENT VP0'
        write(6,*)'AX DEVIENT VP0'
        write(6,*)'AX DEVIENT VP0'

        call correlvp(xp,xpp,vp,ax,fp,ityp)
     end if

     if (lHcyl) then
        call Hcyl
     end if



#if(PARA)
     temps_initspeed=MPI_Wtime()-temps_initspeed_deb
#endif

     if ((itetimestep>0).and.(.not.lcasca)) call deftimestep 

  endif

  if (lcontr) call initcontr(xp,xpp,vp,ax,ityp)


  if (lcasca) then
     call initcasca 
     if (rang==0) then
        if (lfilm) then
           write (lufilmpaf, *) '1'
           write (lufilmpaf, *) 'IT ', '0 ', 'time      0.'
           !            call type (ityp(iko))
           write (lufilmpaf, 114) 'Pb ', xp(1,iko)*1D+8, xp(2,iko)*1D+8, xp(3&
                ,iko)*1D+8, iko
        endif
114     format(a3,1x,3(f10.4,1x),i5)

     end if                                  !  fin rang=0 






     if (itmax==0) stop
     call caltabt 
     if (rang==0)     write(6,*)'>>>>>>>>>>>apres caltabt'
     if (ltabvois) call caltabi 
  end if



  if (ldislo) call at_bord(xp)



  !initialisation pour lcalcjq

  if (lcalcjq.or.lprteat) allocate(eatom(imm))
  if (lprteattotm) allocate(eatomtotm(imm))


  if(lsigat) allocate(sigat(3,3,imm))
  if(lsigtyp) then
     allocate(sigtyp(3,3,ntyp))
     allocate(sigtyptyp(3,3,ntyp,ntyp))
#if(PARA)
     allocate (sigtyp_loc(3,3,ntyp))
     allocate (sigtyptyp_loc(3,3,ntyp,ntyp))
#endif 
  end if

  if (lcdp) then
     call initcdp
     if (itecdp==0)then
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
        write(6,*)'im',im
     end if
  end if

  if ((lheat.EQV..true.).and.(iteheat==0))call heat
  
  if(iteplz>0)  call prtplz(xp,ityp)


  if (dmtype==6) then
     call anapos(it) 
     call arret_ndm
  end if
  !  call sauvegarde 
  if (itmax==0) call arret_ndm

  if(iteanapos>0)then
     itapp=0
     call sauveposition (itapp)
  end if
  if (rang==0) write(6,*)'sortie init'


  if (ldesinteg) then
     itdes=0
     allocate(xpchup(3,imm))
     allocate(xpchdeb(3,imm))
     allocate(xpchdn(3,imm))
     allocate(vpchup(3,imm))
     allocate(vpchdeb(3,imm))
     allocate(vpchdn(3,imm))
     allocate(itichup(imm))
     allocate(itichdn(imm))
     allocate(itichdeb(imm))
     xpchup=xp
     vpchup=vp
     xpchdeb=xp
     vpchdeb=vp

#if(PARA)
     allocate(itichup(imm))
     allocate(itichdn(imm))
     allocate(itichdeb(imm))
     allocate(num_at_globdesdeb(imm))
     allocate(num_at_globdesup(imm))
     allocate(num_at_globdesdn(imm))
     itichdeb=ityp
     itichup=ityp
     imdesup=im ; imdesdeb=im
     num_at_globdesup=num_at_glob
     num_at_globdesdeb=num_at_glob

#endif
  end if




  return
end subroutine init
