module init_mod
!  USE config_mod,only:config
  USE setcell,only:setcellconf
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
  USE sauvegardeT_mod,only: sauvegardeT,cin2gin
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
  USE neb_module,only: constrconfNEB,atneb,cellneb,boxneb
  USE atomconfig,only:atom_config,atom_config_d,ndm2config,config2ndm,atom_config_e
  USE cellconfig, only:cell_config,ndm2cellconfig,cellconfig2ndm,caltabtC,init_cel
  use boxconfig,only: box_config,ndm2boxconfig,boxconfig2ndm
  USE param_det_mod,only: param_det
  USE constrconf_mod, only :constrconf

#ifdef PARA
  USE init_vois_mod,only: init_voisinage
#endif
#ifdef ML
  USE calfo_ml_mod,only: calfo_ml 
#endif
#ifdef LAMMPS_VERSION
  use lammps_util_mod
  use vars_lammps
#endif

  USE gen_com_m, ONLY:igen,ilangevin,iteheat,lcdp,lcorrelvp,ldislo,lhcyl,lheat,itichup,itichdn,itichdeb,formatsauv,iko&
       &,iteanapos,iteplz,iterasmol,itetimestep,itmax,lcalcjq,lcasca,ldesinteg,lcontr,lfilm,lprteat,&
       &lrestart,ltabvois,ltranche,parallele,tmean,tstep,two,umass,usdh,vpchdeb,vpchup,xpchdeb,xpchup,sigat,&
       kinemean,lsigat,pmean,xpchdn,eatomtotm,lprteattotm,vpchdn,indi,nvois,&
       &num_at_globdesdeb,num_at_globdesup,num_at_globdesdn,imdesup,imdesdn,IMDESDEB,npath,&
       &posa,forca,firsttime_lammps,normat,nzl,volu,zls2,celsize,im_glob,fnamcout


  USE var_pot, ONLY:npair,ntrip,r3cm,rumax,typ_and_pot,lpotentiel,l3c,npotmax,rue_pot,ipotentiel,ngrid,csive
  implicit none

contains
  ! **************************************************************
  subroutine init(atdml,boxndm,celndm)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double

    USE tab_imm_m,only:num_at_glob,fp,iwmax,vp,xpp,vp,fp
    USE eam
    USE eamerco
    USE SMjuli
    USE jqmod

    USE posana
    USE defcdp, ONLY :itecdp
    USE elec_cell,ONLY: i2t,t_cpl, readelec
    USE eloss, ONLY : ibrake,ecelec,initeloss

#ifdef PARA
    use mpi
    USE mod_para,only:MPI_COMM_space,TEMPS_INPUT_DEB,TEMPS_INPUT,TEMPS_CONFIG_DEB,TEMPS_CONFIG,MYID,NBR_PROC_VOISIN,TEMPS_INITSPEED_DEB,TEMPS_INITSPEED,nprocs
#endif

    ! **************************************************************

    implicit none
    class(atom_config)::atdml
    type(cell_config),intent(out)::celndm
    type(box_config),intent(out)::boxndm

    integer :: i, lufilmpaf,itapp,ipotcont,j,lenfn2,ipath,ierr
    !-----------------------------------------------
    character*2::extension
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

    !#ifdef LAMMPS_VERSION
    firsttime_lammps=.true.
    if ((ipotentiel==-10).or.(ipotentiel==-11))then
       call init_potential_simple
    else
       !#endif  

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

       !#ifdef LAMMPS_VERSION
    endif
    !#endif  

    !<---------end setting the potential---------------


    !  if (rang == 0)  write(6,*)'cm',cm
    usdh = 1/(two*tstep)
    !endif

    if (ibrake.gt.0) then
       call initeloss
    end if


    it=0




    !<---------setting the configuration by reading gin / cin file --------------


#ifdef PARA
       temps_input=MPI_Wtime()-temps_input_deb

       temps_config_deb = MPI_Wtime()
#endif



#ifdef ML

#else
       if (rang==0)then
          write(6,*)
       write(6,*)' -------------------------------------------------------------------'
       write(6,*)'             definition des rayons de coupure'
    end if
       call param_det
       ! rumax défini en ce point
#endif 

       if (dmtype.ne.9) then
          call constrconf(atdml,boxndm,celndm)

#ifdef PARA
       temps_config=MPI_Wtime()-temps_config_deb
#endif


#ifdef DECOUP
       ! Pas la peine d'aller plus loin dans l'initialisation
       return
#endif
    else

       call constrconfNEB !(xp, xpp, vp, fp, ielat, iwmax, ityp)
       !       do i=1,npath
       !          call print(atneb(1))
       !       end do
    end if
    !...inNEB
    !<---------ends etting the configuration by reading gin /  cin file ---------
    !<---------setting the cell division -------------------------
    ! determination des tailles du nombre de cel. (nox, noy, noz)

    !   call celndm%print

    !    call DynamicalAllocationCell
    !     call celndm%print
#ifdef LAMMPS_VERSION

    if ((ipotentiel==-10).or.(ipotentiel==-11))then
       firsttime_lammps=.true.
       allocate (posa(3*im),  forca(3*im))
       call read_lammps()
    end if
#endif  

    if (dmtype.ne.9)then  !pas NEB
       if (iterasmol>=0) then
          itapp=-1
          call rasmol (atdml,boxndm,itapp)
       end if
!       call config2ndm(atdml,im,imm,xp,fp,ityp,ielat,num_at_glob,ltabvois,iwmax=iwmax,indi=indi,vp=vp,&
!            &xpp=xpp)
!       call cellconfig2ndm(celndm,noxyz,nox,noy,noz,natperc,nato,ncel,atincel,deltadist,celsize) !sans doute inutile
!       call boxconfig2ndm(at,bg,zl,zls2,nzl,volu,normat,boxndm)


       !<---------setting the configuration by generation gin / cin file --------------
       select case (igen)
       case (-1)
          formatsauv = 2 ; fnamcout= fnam(1:lenfnam)//'.cout.'
          call sauvegardeT(atdml,celndm,boxndm,formatsauv,fnamcout)
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
          formatsauv = 2 ; fnamcout= fnam(1:lenfnam)//'.cout.'
          call sauvegardeT(atdml,celndm,boxndm,formatsauv,fnamcout)
          if (rang==0) write (6, *) 'modification terminee'
          call arret_ndm
       case default
       end select

    else !NEB
!       do ipath=2,npath-1
!          call init_cel(cellneb(ipath),cellnox,noy,noz,natperc)
!       end do
    end if

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

    
    if (dmtype.ne.9) then

!   call neigcel
#ifdef PARA
        CALL MPI_BARRIER(MPI_COMM_space,ierr)

   call init_voisinage(celndm)

       

    if (rang==0)  write(6,*) 'NOMBRE DE CELLULES FRONTIERES ASSOCIEES A CHAQUE PROCESSEUR'
    write(6,*) 'Le proc ',myid,' a ',nbr_proc_voisin,' processeur voisin'
    !  do i=1,nbr_proc_voisin
    !     WRITE(6,*) 'Le proc ',myid,' envoit ',nbr_cell_frontiere(i),' vers le proc ',proc_voisin(i)
    !  enddo
    !  WRITE(6,*) 'Le proc ',myid,' recoit ',nbr_cell_ftm,' cell. fantome de ses voisins'
#endif

    ! !!! compcr non pris en charge en parallele !!!


    !<---------end setting the cell diviion ----------------------



       if (ltranche) call layer
       nad(:ntyp) = na(:ntyp)

       call caltabtC(celndm,atdml,lperiod,boxndm)

       if (ltabvois) then
          call caltabi(atdml,celndm)
       end if
       write(6,*)'post caltabi'

!       call config2ndm(atdml,im,imm,xp,fp,ityp,ielat,num_at_glob,ltabvois,iwmax=iwmax,indi=indi,vp=vp,&
!            &xpp=xpp)
!       call cellconfig2ndm(celndm,noxyz,nox,noy,noz,natperc,nato,ncel,atincel,deltadist,celsize) !sans doute inutile
!       call boxconfig2ndm(at,bg,zl,zls2,nzl,volu,normat,boxndm)

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
          select type(atdml)
             class is (atom_config_d)
             call initspeed(atdml,im_glob,boxndm)
             !      call atdml%print          
          end select
          if (iterasmol>=0) then
             itapp=0
             !         call ndm2boxconfig(at,bg,zl,zls2,nzl,volu,normat,boxndm)
             !    call ndm2config(atdml,im,imm,xp,fp,ityp,ielat,num_at_glob=num_at_glob,ltabvois=ltabvois,&
             !         &iwmax=iwmax,indi=indi,nvois=nvois,vp=vp,xpp=xpp)
             call rasmol (atdml,boxndm,itapp)


          end if
       end if
       !  
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
114       format(a3,1x,3(f10.4,1x),i5)
#endif
          if(iteanapos>0)then
             itapp=0
             call sauveposition (itapp)
          end if
       end if

       if (itmax==0) stop
       !                  call ndm2cellconfig(celndm,noxyz,nox,noy,noz,natperc,nato,ncel,atincel,deltadist,celsize)
       !           call ndm2config(atdml,im,imm,xp,fp,ityp,ielat,num_at_glob=num_at_glob,ltabvois=ltabvois,i&
       !                &wmax=iwmax,indi=indi,nvois=nvois,vp=vp,xpp=xpp)
                   call caltabtC(celndm,atdml,lperiod,boxndm)
       if (ltabvois) then
          call caltabi(atdml,celndm)
       end if
!       call config2ndm(atdml,im,imm,xp,fp,ityp,ielat,num_at_glob,ltabvois,iwmax=iwmax,indi=indi,vp=vp,&
!            &xpp=xpp)
!       call cellconfig2ndm(celndm,noxyz,nox,noy,noz,natperc,nato,ncel,atincel,deltadist,celsize) !sans doute inutile



       if (rang==0)     write(6,*)'>>>>>>>>>>>apres caltabt'



       if (dmtype==6) then
          call anapos(it)
          call arret_ndm
       end if
                 write(6,*)'RANG,im',rang,atdml%im,atdml%imm,atdml%xp(3,3)
!    call MPI_FINALIZE(ierr)
!    stop


       if (lcasca) then
          fnamcout = fnam(1:lenfnam)//'.0.cout'
          call sauvegardeT(atdml,celndm,boxndm,formatsauv,fnamcout)
       else
#ifdef PARA
          CALL MPI_BARRIER(MPI_COMM_space,ierr)
          do i=0,nprocs-1
             if (myid==i) then
                write(6,*)
                write(6,*)'PPPPPPPPPPRRRRRRRRRTTTTTT',myid

#endif
!                call atdml%print
!                call celndm%print
!                call boxndm%print
#ifdef PARA                
             end if
             CALL MPI_BARRIER(MPI_COMM_space,ierr)
          end do
#endif
          fnamcout = fnam(1:lenfnam)//'.cout'
          call sauvegardeT(atdml,celndm,boxndm,formatsauv,fnamcout)
       end if
       !  call sauvegarde
       if (itmax==0) call arret_ndm

       if(iteanapos>0)then
          itapp=0
          call sauveposition (itapp)
       end if



       if (ibound==1 .OR. ibound==2 .OR. ibound==3) call init_spebc		!*!
    end if
    if (rang==0) write(6,*)'sortie init'
!    call boxneb%print
    !    call celndm%print
    !   call atdml%print

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
