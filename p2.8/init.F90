module init_mod
  USE config_mod,only:config
  USE divid_mod,only:divid
  USE contrainte,only:initcontr
  USE sauveposition_mod,only:sauveposition
  USE alloc_typ_mod,only: alloc_typ
  USE input_pair_mod,only: input_pair
  USE inputtersoff_mod,only: inputtersoff
  USE calpoeam_mod,only: calpoeam
  USE calpo_mod,only: calpo
  USE transf_mod,only: transf
  USE init_spebc_mod,only: init_spebc
  USE Hcyl_mod,only: Hcyl
  USE neigcel_mod,only: neigcel
  USE tersoff_zbl_mod,only: tersoff_zbl
  USE dynalloccell
  USE initspeed_mod,only: initspeed
  USE caltabt_mod,only: caltabt
  USE sauvegarde_mod,only: sauvegarde,cin2gin
  USE heat_mod,only: heat
  USE caltabi_mod,only: caltabi
  USE creadp_mod,only: creadp
  USE correl_mod,only: correlvp
  USE layer_mod,only: layer
  USE dislo_mod,only: at_bord
  USE initcdp_mod,only: initcdp
  USE initcasca_mod,only: initcasca
  USE deftimestep_mod,only: deftimestep
  USE rasmol_mod,only: rasmol
  USE prtplz_mod,only: prtplz
  USE neb_module,only: configneb
  USE atomconfig
#ifdef PARA
  USE init_vois_mod,only: init_voisinage
#endif
#ifdef ML
  USE calfo_ml_mod,only: calfo_ml 
#endif 
  USE gen_com_m, ONLY:igen,ilangevin,imf,iteheat,lcdp,lcorrelvp,ldislo,lhcyl,lheat,itichup,itichdn,itichdeb,formatsauv,iko,&
       &imana,itdes,iteanapos,iteplz,iterasmol,itetimestep,itmax,lcalcjq,lcasca,ldesinteg,lcontr,lfilm,lprteat,&
       &lrestart,lsigtyp,ltabvois,ltranche,parallele,tmean,tstep,two,umass,usdh,vpchdeb,vpchup,xpchdeb,xpchup,sigat,sigtyp,&
       kinemean,lsigat,pmean,xpchdn,sigtyptyp,sigtyp,eatomtotm,lprteattotm,vpchdn,indi,nvois,sigtyp_loc,sigtyptyp_loc,&
       &num_at_globdesdeb,num_at_globdesup,num_at_globdesdn,imdesup,imdesdn,IMDESDEB

  
      USE var_pot, ONLY:npair,ntrip,r3cm,rumax,typ_and_pot,lpotentiel,l3c,npotmax,rue_pot,ipotentiel
  implicit none 
contains
  ! **************************************************************
  subroutine init
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double

    USE tab_imm_m
    USE eam
    USE eamerco
    USE SMjuli
    USE jqmod

    USE posana
    USE defcdp, ONLY :itecdp
    USE elec_cell,ONLY: i2t,t_cpl, readelec
    USE eloss, ONLY : ibrake,ecelec,initeloss

#ifdef PARA
    USE mod_para
#endif

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
    integer :: i, lufilmpaf,itapp,ipotcont,j,lenfn2
    integer :: complet=1    ! flag d'appel a divid : complet : exec de la routine complete
    !-----------------------------------------------
    character*2::extension
    type(atom_config_d)::atdml
    tmean = 0.0
    pmean = 0.0
    timel = 0.0
    kinemean = 0.0
    lufilmpaf = 79

    !     write(6,*)'entree dans init.f'
    !potentiel BKS
#ifdef PARAPH
    rang=rangph
#endif


#ifdef PARA
    temps_input_deb = MPI_Wtime()
#endif

    !>---------allocating the types------------------
    if (npotentiel.gt.1)then
       ipotentiel=-1
       npair=  ntyp*(ntyp+1)/2 ; ntrip= ntyp*ntyp *(ntyp+1)/2
       allocate(typ_and_pot(ntyp,npotmax))
       typ_and_pot(:,:)=.false.
       call  alloc_typ
       typ_pot_pair=0
    end if

    allocate(rue_pot(npotmax))
    rue_pot(:)=0.

    !>---------setting the potential------------------
    if (rang.eq.0) then
       write(6,*)
       write(6,*)'-*-*-*-*-*-*-*POTENTIELS*-*-*-*-*-'
       write(6,*)
    end if

#ifdef LAMMPS_VERSION
    firsttime_lammps=.true.
    if ((ipotentiel==-10).or.(ipotentiel==-11))then
       call init_potential_simple
    else
#endif  

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
             case(11)
                if (rang.eq.0) then
                   write(6,*)

                   write(6,*)'POTENTIEL EAM ERCOLESI'
                   write(6,*)
                end if
                call inputeamerco(ntyp,npair,ntrip,cm,catom,ty,umass,rue_pot(ipotentiel),rumax,&
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
                !           if (rang==0)   write(6,*)'POTENTIEL tersoff.potin'
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

#ifdef ML
                ! MiLaDy
             case(20)
                if (rang.eq.0) then
                   write(6,*)
                   write(6,*)' ML ..... set-up MiLady potential'
                   write(6,*)
                end if
                !This comes with MiLaDy package
                call md_init_potential_ml
#endif
             end select

          endif
       end do

#ifdef LAMMPS_VERSION
    endif
#endif  

    !<---------end setting the potential---------------


    !  if (rang == 0)  write(6,*)'cm',cm
    usdh = 1/(two*tstep)
    !endif

    if (ibrake.gt.0) then
       call initeloss
    end if


    it=0




    !<---------setting the configuration by reading gin / cin file --------------
    !---inNEB
    if (dmtype.ne.9) then
#ifdef PARA
       temps_input=MPI_Wtime()-temps_input_deb

       temps_config_deb = MPI_Wtime()
#endif
       call config
#ifdef PARA
       temps_config=MPI_Wtime()-temps_config_deb
#endif


#ifdef DECOUP
       ! Pas la peine d'aller plus loin dans l'initialisation
       return
#endif
    else

       call configNEB(xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
    end if
    !...inNEB
    !<---------ends etting the configuration by reading gin /  cin file ---------

#ifdef LAMMPS_VERSION

    if ((ipotentiel==-10).or.(ipotentiel==-11))then
       firsttime_lammps=.true.
       allocate (posa(3*im),  forca(3*im))
       call read_lammps()
    end if
#endif  


    if (iterasmol>=0) then
       itapp=-1
       call rasmol (itapp)
    end if

    !<---------setting the configuration by generation gin / cin file --------------
    select case (igen)
    case (-1)
       formatsauv = 2
       call sauvegarde
       if (rang==0) write (6, *) 'generation terminee'
       call arret_ndm
       !  case (0)
       !     if (rang==0) write (6, *) 'generation du crystal ; puis run'
       !  case (1)
       !     if (rang==0) write (6, *) 'run a partir du fichier .cin'
    case (2)
       call cin2gin
       call arret_ndm

    case (3)
       call transf
       formatsauv = 2
       call sauvegarde
       if (rang==0) write (6, *) 'modification terminee'
       call arret_ndm
    case default
    end select


    !<---------setting the cell division -------------------------
    ! determination des tailles du nombre de cel. (nox, noy, noz)

    call divid (complet)
    call DynamicalAllocationCell

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

#ifdef PARA
    call init_voisinage()

    if (rang==0)  write(6,*) 'NOMBRE DE CELLULES FRONTIERES ASSOCIEES A CHAQUE PROCESSEUR'
    write(6,*) 'Le proc ',myid,' a ',nbr_proc_voisin,' processeur voisin'
    !  do i=1,nbr_proc_voisin
    !     WRITE(6,*) 'Le proc ',myid,' envoit ',nbr_cell_frontiere(i),' vers le proc ',proc_voisin(i)
    !  enddo
    !  WRITE(6,*) 'Le proc ',myid,' recoit ',nbr_cell_ftm,' cell. fantome de ses voisins'
#endif

    ! !!! compcr non pris en charge en parallele !!!


    !<---------end setting the cell division ----------------------



    imd = im
    if (ltranche) call layer


    nad(:ntyp) = na(:ntyp)
!!!  endif

    imf=im
    imana=min(imd,imd)


    !computing the neighbours for the very first time ......
    !  if (itmax>0) then
    if (rang==0)write(6,*)'1ER CALL init'
    call caltabt(im,xp,ielat)

    ! if (rang==0)  write(6,*)'>>>>>>>>>>>apres caltabt'
    if (ltabvois) then
       call ndm2config(atdml,im,imm,xp,fp,ityp,ielat,num_at_glob=num_at_glob,ltabvois=ltabvois,&
            &iwmax=iwmax,indi=indi,nvois=nvois,vp=vp,xpp=xpp)
       call caltabi(atdml%atom_config)
       call config2ndm(atdml,im,imm,xp,fp,ityp,ielat,num_at_glob,ltabvois,iwmax=iwmax,indi=indi,vp=vp,xpp=xpp)
    end if
!       call caltabi
    ! if (rang==0)  write(6,*)'>>>>>>>>>>>apres caltabi'
    !  end if
    !computing the neighbours for the very first time ......



#ifdef ML
    ! MiLaDy
    if(ipotentiel==20) then
       if (rang.eq.0) then
          write(6,*)
          write(6,*)' ML  ..... configuration MiLady '
          write(6,*)
       end if
       !This comes with MiLaDy Package
       call md_init_config_ml
    end if
#endif


    if (L2T.eqv..true.) then
       call readelec
       if (rang==0) write(6,*)'!*!*!*!*! 2T MD version =', i2t,'*!*!*!*!'
       dmtype=4
       ibrake=1
       ilangevin=1
       if((ecelec==0))then
          write(6,*) 'eccelec<>0  and l2T : STOP'
          stop
       end if
       if ((i2T==0).and.(t_cpl.lt.0)) then
          write(6,*) 'i2T=0 t_cpl<0 and l2T : STOP'
          stop
       end if
       if (nox.le.0 ) then
          write(6,*) 'nox noy noz MUST be defined in .din with 2T: STOP'
          stop
       end if


    endif

    if (.not.lrestart) then
       !    if (rang==0)     write(6,*)'>>>>>>>>>>>avant initspeed'
#ifdef PARA
       temps_initspeed_deb = MPI_Wtime()
#endif

       ! input and initialization of 2T
       call initspeed
       if (iterasmol>=0) then
          itapp=0
          call rasmol (itapp)
       end if
    end if
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
    !end init the speed using Maxwell proba density-----------------



#ifdef PARA
    temps_initspeed=MPI_Wtime()-temps_initspeed_deb
#endif

    if ((itetimestep>0).and.(.not.lcasca)) call deftimestep


    if (lcontr) call initcontr(xp,xpp,vp,ax,ityp)


    if (lcasca) then
       call initcasca
#if PARA
#else

       if (lfilm) then
          write (lufilmpaf, *) '1'
          write (lufilmpaf, *) 'IT ', '0 ', 'time      0.'
          write (lufilmpaf, 114) 'Pb ', xp(1,iko)*1D+8, xp(2,iko)*1D+8, xp(3&
               ,iko)*1D+8, iko
       endif
114    format(a3,1x,3(f10.4,1x),i5)
#endif
       if(iteanapos>0)then
          itapp=0
          call sauveposition (itapp)
       end if

       if (itmax==0) stop
    call caltabt(im,xp,ielat)
       if (rang==0)     write(6,*)'>>>>>>>>>>>apres caltabt'
       if (ltabvois) then
          call ndm2config(atdml,im,imm,xp,fp,ityp,ielat,num_at_glob=num_at_glob,ltabvois=ltabvois,&
               &iwmax=iwmax,indi=indi,nvois=nvois,vp=vp,xpp=xpp)
       call caltabi(atdml%atom_config)
       call config2ndm(atdml,im,imm,xp,fp,ityp,ielat,num_at_glob,ltabvois,iwmax=iwmax,indi=indi,vp=vp,xpp=xpp)
    end if
       
    end if



    if (ldislo) call at_bord(xp)



    !initialisation pour lcalcjq

    if (lcalcjq.or.lprteat) allocate(eatom(imm))
    if (lprteattotm) allocate(eatomtotm(imm))


    if(lSigat) allocate(sigat(3,3,imm))
    if(lsigtyp) then
       allocate(sigtyp(3,3,ntyp))
       allocate(sigtyptyp(3,3,ntyp,ntyp))
#ifdef PARA
       allocate (sigtyp_loc(3,3,ntyp))
       allocate (sigtyptyp_loc(3,3,ntyp,ntyp))
#endif
    end if

    if (lcdp) then
       call initcdp
       if (itecdp==0)then
          call creadp (xp, xpp, ityp,vp)
          call caltabt(im,xp,ielat)
          
          if (ltabvois) call caltabi(atdml%atom_config)

          if (lperiod) then
             call period (imm,xp,xpp,ax)
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

#ifdef PARA
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

    if (ibound==1 .OR. ibound==2 .OR. ibound==3) call init_spebc		!*!

    if (rang==0) write(6,*)'sortie init'

    return
  end subroutine init

  subroutine init_potential_simple
    USE T_kind_param_m, ONLY:  double
    USE gen_com_m, ONLY: rang,A2cm,umass
    USE var_pot, ONLY: ntyp, npair, ntrip,cm,catom, ty,rue_pair,ipotentiel,q
#ifdef PARA
    USE mpi
#endif
    implicit none
    integer :: i,error, beggin,  endding,lupotin
    character ::  fnampotin*80
    real(double)::rue


    fnampotin = 'simple.potin'
    lupotin = 95
    open(unit=lupotin, file=fnampotin, status='old')
    read(lupotin,*)ntyp
    npair=  ntyp*(ntyp+1)/2 ; ntrip= ntyp*ntyp *(ntyp+1)/2
    call  alloc_typ
    read(lupotin,*) rue
    rue=rue*A2cm
    rue_pair(:)=rue
    select case (ipotentiel)
    case(-10)
       do i = 1, ntyp
          read (lupotin,*) cm(i),catom(i),ty(i)
          if (rang/=0) cycle
          write (6, '(I4,2F9.3,A5)') i, cm(i),catom(i),ty(i)
       end do
    case(-11)
       do i = 1, ntyp
          read (lupotin,*) cm(i),catom(i),ty(i),q(i)
          if (rang/=0) cycle
          write (6, '(I4,2F9.3,A5,F9.3)') i, cm(i),catom(i),ty(i),q(i)
       end do

    end select
    cm(:ntyp) = cm(:ntyp)*umass

    close (lupotin)
  end subroutine init_potential_simple
end module init_mod
