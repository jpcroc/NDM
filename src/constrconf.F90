module constrconf_mod
  USE arret_ndm_mod,only:arret_ndm
#ifdef PARA
  USE decoupage_mod,only: decoupage,decoup2im,constrandrepart
#endif  
  USE read_val,only:imm,ipbc,nox,noy,noz
  USE gen_com_m, ONLY: lenfnam, fnam,fmt_cin,igen,imm_glob,ldecoup,lperiod,lrestart,rang,&
       &lvpread,zero,low_limit,lspacendm,rang,dmtype
  USE var_pot, ONLY:ntyp,rumax,ipotentiel
  use cryst_to_cart_mod,only:cryst_to_cart
  USE arret_ndm_mod,only: arret_ndm
  USE atomconfig,only:atom_config,atom_config_d,atom_config_e
  USE cellconfig,only:cell_config
  USE boxconfig,only:box_config,periodbox
  USE setcell,only:setnox,setcellconf
  USE decoupage_mod,only: decoupage
  USE Mat_utils_mod,only: Matinv_gen,is_upper_triangular
  USE T_kind_param_m, ONLY:  double
#ifdef PARA
  USE Tpara,only:COMM_space,nprocspace
#endif
  use Tpara,only:para_space_config,nprocspace
  use coord_to_cell_mod
  use config2data_mod,only:config2data
#ifdef ML
  !use gen_com_m_ml, only: at, im
  use derived_types, only: config_real
#endif



  implicit none
  logical::lprt=.true.
  logical :: ldecalcor
  logical :: lsecondpath
  real(double),allocatable::xpd(:,:)
contains
  subroutine constrconf (atrcf,boxrcf,cellrcf,lrepart,filename,psc)
    !********************************************************************
    !             CONSTRUCTION DE LA BOITE DE SIMULATION
    !********************************************************************
    implicit none
    class(atom_config),intent(inout)::atrcf
    type(cell_config),intent(out)::cellrcf
    class(box_config),intent(out)::boxrcf
    logical::lrepart
    character(len=*),optional::filename
    character*80::filenom
    integer :: i, iti
    type(para_space_config)::psc
    !    integer,      dimension(:), allocatable   :: num_at_buff
    integer, dimension(:),allocatable     :: ibuffer
    real(double), dimension(:,:),allocatable    :: buffer
    integer::ncore,ic

#ifdef PARA
    class (atom_config),allocatable::COMPatrcf
#endif
    character :: fnamcin*80, fnamgin*80
    integer::itread
    integer::nati
    logical::lwrite
#ifdef DKIO
    character(len=6) :: format
#endif
    !-----------------------------------------------------
    ! READING FROM THE CONFIGURATION FILE
    !---------------------------------------------------
    filenom=fnam(1:lenfnam)
    if (present(filename))filenom=filename
    if ((rang==0).and.(lprt)) then
       write(6,*)
       write(6,*)' *-*-*-*-*-*CONSTRUCTION DE LA BOITE*-*-*-*-*-*-'
       write(6,*)'imm,imm_glob,nprocspace',imm,imm_glob,nprocspace
       write(6,*)
    endif

    if (igen.ge.1 .and. igen.le.3) then
       allocate (ibuffer(imm_glob))
       allocate (buffer(3,imm_glob))
       if ((rang==0).and.(lprt))  write(6,*)'********** reading configuration from file********'
       if (lrestart) then
          fnamcin = fnam(1:lenfnam)//'.cout'
       else
          fnamcin = fnam(1:lenfnam)//'.cin'
       end if

#ifdef PARA
       if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
          itread=0
          call read_cin(boxrcf,itread,fnamcin=fnamcin,fmtcin=fmt_cin) !0=at seulement; 1=complet; 2 = at, xp et num_at_glob seulement
          if ((rang==0).and.(lprt)) then
             write (6, '(A,D15.8,A,D15.8,A)') 'volume=', boxrcf%volu,' cm3 ',boxrcf%volu*1d24,' Ang3'
          end if
          call setnox(boxrcf,cellrcf,rumax,lverbose=lprt,noxr=nox,noyr=noy,nozr=noz)
          ncore=0
          call  decoupage(nprocspace,ncore,cellrcf,psc=psc,lverbose=lprt)

          if ((rang==0).and.(lprt)) then
             write (6, '(A,D15.8,A,D15.8,A)') 'volume=', boxrcf%volu,' cm3 ',boxrcf%volu*1d24,' Ang3'
          end if
          ncore=0
          if (lrepart.eqv..true.) then
             call read_cin2(boxrcf,atrcf,cellrcf,imm_glob,fnamcin,lrestart,fmt_cin,psc) !0=at seulement; 1=complet; 2 = at, xp et num_at_glob seulement
          else
             itread=1
             call read_cin(boxrcf,itread,atrcf,imm_glob,fnamcin,lrestart,fmt_cin) !0=at seulement; 1=complet; 2 = at, xp et num_at_glob seulement
          end if
       else
          itread=1
          call atrcf%init(immin=imm_glob,imin=0,ltabvois=atrcf%ltabvois,rvois=atrcf%rvois)
          call read_cin(boxrcf,itread,atrcf,imm,fnamcin,lrestart,fmt_cin) !0=at seulement; 1=complet; 2 = at, xp et num_at_glob seulement , 3 trié par num_at_buff
          atrcf%im_glob=atrcf%im
          if ((rang==0).and.(lprt)) then
             write (6, '(A,D15.8,A,D15.8,A)') 'volume=', boxrcf%volu,' cm3 ',boxrcf%volu*1d24,' Ang3'
          end if
          call setnox(boxrcf,cellrcf,rumax,lverbose=lprt,noxr=nox,noyr=noy,nozr=noz)
          ncore=0

       end if
#else


       if (ldecoup) then 
          itread=0
          call read_cin(boxrcf,itread,fnamcin=fnamcin)
          if ((rang==0).and.(lprt)) then
             write (6, '(A,D15.8,A,D15.8,A)') 'volume=', boxrcf%volu,' cm3 ',boxrcf%volu*1d24,' Ang3'
          end if
          call setnox(boxrcf,cellrcf,rumax,lverbose=lprt,noxr=nox,noyr=noy,nozr=noz)
          open(123, file='decoup.dat', status='old')
          read (123, *) nprocspace,ncore
          close(123)         
          call  decoupage(nprocspace,ncore,cellrcf,psc=psc,lverbose=lprt)
          call arret_ndm
       else
          itread=1
          call atrcf%init(immin=imm_glob,imin=0,ltabvois=atrcf%ltabvois,rvois=atrcf%rvois)
          call read_cin(boxrcf,itread,atrcf,imm,fnamcin,lrestart,fmt_cin) !0=at seulement; 1=complet; 2 = at, xp et num_at_glob seulement , 3 trié par num_at_buff

          atrcf%im_glob=atrcf%im
          if ((rang==0).and.(lprt)) then
             write (6, '(A,D15.8,A,D15.8,A)') 'volume=', boxrcf%volu,' cm3 ',boxrcf%volu*1d24,' Ang3'
          end if
          call setnox(boxrcf,cellrcf,rumax,lverbose=lprt,noxr=nox,noyr=noy,nozr=noz)
          !          CALL fin allocation CELL et FIN DIVID
       end if

#endif
       call setcellconf(cellrcf,atrcf,boxrcf,rumax)
       deallocate (ibuffer)
       deallocate (buffer)

#ifdef DKIO

    else if (igen.ge.11 .and. igen.le.29) then
       !-----------------------------------------------------
       ! BUILDING OF THE CRISTAL FROM DK-IO
       !-----------------------------------------------------
       lvpread=.false.
       select case (igen)
       case(11,21)
          fnamgin = fnam(1:lenfnam)//'.cfg'
          format='xfg'
          if (igen==21) lvpread=.true.
       case(12,22)
          fnamgin = fnam(1:lenfnam)//'.xfg'
          format='xfg'
          if (igen==22) lvpread=.true.
       case(13,23)
          fnamgin = fnam(1:lenfnam)//'.cell'
          format='castep'
          if (igen==23) lvpread=.true.
       case(14)
          fnamgin = fnam(1:lenfnam)//'.cif'
          format='cif'
       case(15,25)
          fnamgin = fnam(1:lenfnam)//'.CONFIG'
          format='dlpoly'
          if (igen==25) lvpread=.true.
       case(16)
          fnamgin = fnam(1:lenfnam)//'.gulp'
          format='gulp'
       case(17)
          fnamgin = fnam(1:lenfnam)//'.lmp'
          format='lammps'
       case(18)
          fnamgin = fnam(1:lenfnam)//'.POSCAR'
          format='vasp'
       case(19)
          fnamgin = fnam(1:lenfnam)//'.xyz'
          format='xyz'
       end select
       call dkio2ndm(atrcf,cellrcf,boxrcf,fnamgin,rumax,format,lrepart,psc)
       call periodbox (boxrcf,atrcf)

       select type(atrcf)
       class is (atom_config_e)
!          atrcf%xpp(:,1:atrcf%im)=atrcf%xp(:,1:atrcf%im)
          if (atrcf%lax) then
             atrcf%ax(:,1:atrcf%im)=atrcf%xp(:,1:atrcf%im)
          end if
       end select

#endif

       !-----------------------------------------------------
       ! BUILDING OF THE CRISTAL FROM .GIN FILE
       !-----------------------------------------------------
    else    ! igen.eq.0



       lvpread=.false.

       ! open fichier .gin
       fnamgin = fnam(1:lenfnam)//'.gin'
       call gin2ndm(atrcf,cellrcf,boxrcf,fnamgin,rumax,lrepart,psc)
       call periodbox (boxrcf,atrcf)

       select type(atrcf)
       class is (atom_config_e)
          !          atrcf%xpp(:,1:atrcf%im)=atrcf%xp(:,1:atrcf%im)
          if (atrcf%lax) then
             atrcf%ax(:,1:atrcf%im)=atrcf%xp(:,1:atrcf%im)
          end if
       end select

    end if


    if ((rang==0).and.(lprt)) then

       write(6,*)
       write (6, *) '-------- boite de simulation ------'
       write (6, *) 'nombre d atomes =', atrcf%im_glob
       do i=1,3
          write(6,'(A,I2,3F15.6)')'vecteur ',i, (boxrcf%at(ic,i)*1.0d8,ic=1,3)
       end do
    endif                                  ! fin rang=0

    do iti = 1, ntyp
       nati=count(atrcf%ityp==iti)
#ifdef PARA
       if ((nprocspace.gt.1).and.(lspaceNDM)) then
          call comm_space%sum(nati)
       end if
#endif
       if ((rang==0).and.(lprt)) then
          if (nati.ne.0) write (6, *) nati, ' atomes de type', iti
       end if

    end do



    if((ipotentiel==-10).or.(ipotentiel==-11)) then
       !       if(rang==0)then
       !          write(6,*)'write configuration to conf.lmp RANG 0 !'
       if (rang==0) then
          lwrite=.true.
       else
          lwrite=.false.
       end if
       call config2data (imm,atrcf%im,atrcf%xp,atrcf%ityp,boxrcf%at,ntyp,lwrite)
       !       end if
    end if
    !#endif     


    !       end if



    return

456 print *,'Erreur dans la lecture du fichier .cin, verifier son format&
         & et fmt_cin ATTENTION A BIG_ENDIAN !! SI COMMPILE BIG_ENDIAN NE LIT PLUS QUE CA'


  end subroutine constrconf


  subroutine gin2ndm(at2b,cel2b,box2b,fnamg,rum,lrepartition,psc,lconstrsimple,immread)
    type(para_space_config),optional::psc
    class(atom_config)::at2b
    type(cell_config)::cel2b
    class(box_config)::box2b
    character,intent(in) :: fnamg*80
    real(double),intent(in)::rum
    logical,optional,intent(in)::lrepartition,lconstrsimple
    integer,optional::immread
    integer::immr,npr
    logical::lcs
    logical::lrepart
    !    type (atom_config)::COMPatrcf
    type(atom_config)::atrgin
    type(box_config)::boxrgin
    real(double)::atg(3,3)
    integer::lat(3),ic,ncore,itread
    lrepart=.true.
    lcs=.false.
    if(present(lrepartition))lrepart=lrepartition
    if(present(lconstrsimple))lcs=lconstrsimple
    immr=imm_glob
    if (present(immread)) immr=immread
    !    write(6,*)'IMMR',immr,lcs
    if (ldecoup) then
       itread=0
    else
       itread=1
    end if

    call read_gin(boxrgin,atrgin,fnamg,lat,itread=itread,immread=immr)
    do ic=1,3
       atg(:,ic)=boxrgin%at(:,ic)*lat(ic)
    end do
    call box2b%init(atg,ipbc)

    call setnox(box2b,cel2b,rum,lverbose=lprt,noxr=nox,noyr=noy,nozr=noz)
    if ((rang==0).and.(lprt)) then
       write (6, '(2A,D15.8,A,D15.8,A)') fnamg,'volume=', box2b%volu,' cm3 ',box2b%volu*1d24,' Ang3'
    end if

    if (lcs) then ! construction simpple sans repartition en sequentiel
       call constr_2gin (at2b,atrgin,lat,immr)
       call cryst_to_cart (at2b%imm, at2b%xp, box2b%at, 1)
       at2b%im_glob=at2b%im
       call setcellconf(cel2b,at2b,box2b,rum)
       return
    end if

#ifdef PARA
    if (rang==0) then
       if (ldecoup) then
          open(123, file='decoup.dat', status='old')
          read (123, *) npr,ncore
          close(123)         
          call  decoupage(npr,ncore,cel2b,psc=psc,lverbose=lprt)
          call arret_ndm
       end if
    end if
    !    COMPatrcf%ltabvois=at2b%ltabvois; compatrcf%nvois=at2b%nvois; compatrcf%rvois=at2b%rvois
    call  decoupage(nprocspace,ncore,cel2b,psc=psc,lverbose=lprt)
    ncore=0
    at2b%imm_glob=imm_glob
    at2b%im_glob=atrgin%im*lat(1)*lat(2)*lat(3)

    if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.).and.(lrepart)) then
       call constrandrepart(atrgin,at2b,cel2b,box2b,lat,psc)
!!$       call constr_2gin (COMPatrcf,atrgin,lat,imm_glob)
!!$       call cryst_to_cart (COMPatrcf%imm, COMPatrcf%xp, box2b%at, 1)
!!$       compatrcf%imm_glob=imm_glob
!!$       call  decoup2im(nprocspace,cel2b,at2b,psc=psc,atcomp=compatrcf,boxrep=box2b)
!!$       call repartition(COMPatrcf,at2b,box2b,cel2b)
    else
       call constr_2gin (at2b,atrgin,lat,imm_glob)
       call cryst_to_cart (at2b%imm, at2b%xp, box2b%at, 1)
       at2b%imm_glob=imm_glob

       !       call compatrcf%copy_config(at2b, lrescl=.true.)
    end if
#else
    if (ldecoup) then
       open(123, file='decoup.dat', status='old')
       read (123, *) npr,ncore
       close(123)         
       call  decoupage(npr,ncore,cel2b,psc=psc,lverbose=lprt)
       call arret_ndm
    end if
    call constr_2gin (at2b,box2b,cel2b,atrgin,boxrgin,lat,imm)
    call cryst_to_cart (at2b%imm, at2b%xp, box2b%at, 1)
    at2b%im_glob=at2b%im
#endif
    call setcellconf(cel2b,at2b,box2b,rum)
    return

  end subroutine gin2ndm


  subroutine constr_2gin(atrcf,atrgin,lat,immread)

    class(atom_config),intent(inout)::atrcf

    type(atom_config)::atrgin
    integer,intent(in)::lat(3)
    integer,intent(in),optional::immread

    integer::i,ia,ib,ic,icell,imloc,immr

    real(double)::rvN
    logical::liniint
    !imtot=lat(1)*lat(2)*lat(3)*atrgin%im
    imloc=lat(1)*lat(2)*lat(3)*atrgin%im
    immr=imm_glob
    if (present(immread)) immr=immread
    if (imloc>immread) then
       write (6, *) rang,'imm trop petit',imloc,immr
       call arret_ndm
    endif

    if (atrcf%rvois.gT.0)then
       rvn=atrcf%rvois
    else
       rvn=0
    end if
    if (present(immread))then

       call atrcf%init(imloc,immin=immread,ltabvois=atrcf%ltabvois,nvois=atrcf%nvois,rvois=rvn,im_glob=imloc)

    else
       call atrcf%init(imloc,ltabvois=atrcf%ltabvois,nvois=atrcf%nvois,rvois=rvn,im_glob=imloc)
    end if

    !    imtot=imloc
    i=0
    do ia = 1,lat(1)
       do ib = 1,lat(2)
          do ic = 1,lat(3)
             do icell = 1, atrgin%im
                i  = i + 1
                atrcf%xp(1,i) = (atrgin%xp(1,icell)+float(ia-1))/float(lat(1))
                atrcf%xp(2,i) = (atrgin%xp(2,icell)+float(ib-1))/float(lat(2))
                atrcf%xp(3,i) = (atrgin%xp(3,icell)+float(ic-1))/float(lat(3))
                atrcf%num_at_glob(i)=i
                atrcf%ityp(i)=atrgin%ityp(icell)
             end do
          end do
       end do
    end do

    return
  end subroutine constr_2gin

  subroutine repartition(atcomp,atrep,boxrep,cellrep,nab)
#ifdef PARA
    use Tpara,only:myidsp
#endif
    use paraconfig,only:para_config
    class(atom_config)::atrep
    class(atom_config)::atcomp
    type(cell_config)::cellrep
    class(box_config)::boxrep
    integer,optional, dimension(:), allocatable   :: nab
    integer::i,icomp,iti,im,numcell,numproc
    real(double)::xt(3)

#ifdef PARA  


    i=0;im=0
    atrep%im_glob=atcomp%im_glob
    do icomp=1,atcomp%im


       xt(:)=atcomp%xp(:,icomp)
       iti = atcomp%ityp(icomp)
       call coord_to_cell(xt,numcell,boxrep,cellrep%nox,cellrep%noy,cellrep%noz)
       numproc=cellrep%proc_cell(numcell)
       atcomp%proc_at(icomp)=numproc
       if (numproc == myidsp) then
          i=i+1
          im=im+1
          call atcomp%copy_atom(icomp,atrep,i)
          !             atrep%num_at_glob(i)=atcomp%num_at_glob(icomp)
          if (present(nab)) nab(i)=icomp
          atrep%xp(:,i)=xt(:)
          !             atrep%ityp(i)=iti
          atrep%proc_at(i)=myidsp

       endif
    end do
    atrep%im=im
#endif
    return
  end subroutine repartition

  subroutine read_cin2(boxcin,atcinr,celcf,immr,fnamcin,lres,fmtcin,psc)
    USE gen_com_m,only: iteration,itmax,nitmax,pmean,oldtstep,timel,two,usdh,dilat,tmean,tstep
    use Tpara,only:myidsp,nprocspace
    use read_val,only:ltabvois
    type(para_space_config)::psc
    character,intent(in) :: fnamcin*80
    class(box_config)::boxcin
    class(atom_config)::atcinr
    type(cell_config)::celcf
    integer,intent(in),optional::immr,fmtcin
    logical,intent(in),optional::lres
    real(double)::at(3,3),rvois0
    real(double),allocatable::xpr(:,:)
    integer,allocatable::itypr(:),proc(:),natgr(:)
    integer::natlocm,icomp,numcell,numproc,imm_loc1,im,cellules_max,cellules_int,im0,im_glob,imm_loc,imm,ig,&
         &lucin,icintype,icintypemod,i,nvois0,it,im_gr,ii
    integer,allocatable::natloc(:) !indice de boucle
    logical::lprt=.true.
    lucin = 94
    open(unit=lucin, file=fnamcin, form='unformatted', status='old', err=431)

    read (lucin, err=432) icintype

    if ((rang==0).and.(lprt))  write (6, *) 'config type of  .cin file : ', icintype
    if (icintype>5.or.icintype<0) then
       write (6, *) rang, 'wrong icintype'
       call arret_ndm
    endif




    allocate (natloc(0:nprocspace-1))
    natloc=0

    allocate(xpr(3,immr))
    allocate(itypr(immr))
    allocate(proc(immr))
    allocate(natgr(immr))

    icintypemod = mod(icintype,2)
    read (lucin, err=433) at
    if(dilat(1).ne.0.0)then
       do i=1,3
          at(i,:)=at(i,:)*dilat(i)
       end do
    end if

    call boxcin%init(at,ipbc)


    read (lucin, err=434) im_gr                         !number of atoms in the box
    if (im_gr>immr) then
       if(rang==0)                    write (6, *) 'P2 im > imM', im_gr, immr
       call arret_ndm
    endif
    !       call atcinr%init(im_gr,immr,im_glob=im_gr)

    read (lucin, err=435) itypr   !ityp muet
    if ((rang==0).and.(lprt))  write (6, *) 'types'
    read (lucin, err=436) xpr    ! xp
    if ((rang==0).and.(lprt))  write (6, *) 'XPR'
    !       atcinr%xp(1:3,1:im_gr)=buffer(1:3,1:im_gr)

    cellules_max=0
    cellules_int=0
    do ii = 0,nprocspace-1
       cellules_max = max(cellules_max,(psc%res_cpu(ii,1)+2) * (psc%res_cpu(ii,2)+2)* (psc%res_cpu(ii,3)+2))
       cellules_int = max(cellules_int,(psc%res_cpu(ii,1)+0) * (psc%res_cpu(ii,2)+0)* (psc%res_cpu(ii,3)+0))
    enddo
    cellules_max = min (cellules_max, celcf%noxyz)

    formcin:select case (fmt_cin)
    case (0) formcin
       do i=1,im_gr
          natgr(i) = i
       enddo
    case(1) formcin
       read (lucin, err=437) natgr
       if ((rang==0).and.(lprt))  write (6, *) 'num_at_glob'
    case default  formcin
       if (rang.eq.0) write(6,*) 'precisez le format fmt_cin'
       call arret_ndm
    end select formcin
    !       atcinr%num_at_glob(1:im_gr)=ibuffer(1:im_gr)

    do i=1,im_gr
       call coord_to_cell(xpr(:,i),numcell,boxcin,celcf%nox,celcf%noy,celcf%noz)

       proc(i)=celcf%proc_cell(numcell)
       if (proc(i) == myidsp) then
          natloc(myidsp)=natloc(myidsp)+1
       endif
    end do

    call comm_space%sum(natloc)
    !             if (rang==0) write(6,*)'natloc',natloc
    natlocm=maxval(natloc)
    natlocm=int(natlocm*float(cellules_max)/cellules_int)
    imm_loc=min( imm_glob, int(1.2 * natlocm))
    imm = imm_loc
    im0=0;rvois0=0
    call atcinr%init(im0,imm,ltabvois,nvois0,rvois0,im_glob=im_gr,imm_glob=imm_glob)
    i=0
    do it=1,im_gr
       if (proc(it) == myidsp) then
          i=i+1
          atcinr%xp(:,i)=xpr(:,it)
          atcinr%ityp(i)=itypr(it)
          atcinr%num_at_glob(i)=natgr(it)
#ifdef PARA
          atcinr%proc_at(i)=myidsp
#endif
          
       endif
    end do
    atcinr%im=i
    !dans la suite xpr est un simple buffer
    select type(atcinr)
    type is (atom_config)
       if (icintypemod==1) then
          if (icintype==3)read (lucin, err=438) xpr                     !xpp
          read (lucin, err=439) xpr                     !vp

       end if

       lvpread=.false.
    type is (atom_config_d)
       if (icintypemod==1) then
          if (icintype==3) read (lucin, err=440) xpr                     !xpp
          read (lucin, err=441) xpr                     !vp
          if ((rang==0).and.(lprt))  write (6, *) 'VP'
          i=0
          do it=1,im_gr
             if (proc(it) == myidsp) then
                i=i+1
                atcinr%vp(:,i)=xpr(:,it)
             endif
          end do
       end if
    end select
    select type(atcinr)
    class is (atom_config_e)
       if (icintypemod==1) then

          if (icintype==3) then
             read (lucin, err=442) xpr                     !xpp
          end if

          read (lucin, err=443) xpr                     !vp
          if ((rang==0).and.(lprt))  write (6, *) 'VP'
          i=0
          do it=1,im_gr
             if (proc(it) == myidsp) then
                i=i+1
                atcinr%vp(:,i)=xpr(:,it)
             endif
          end do

       else
       end if
       if (atcinr%lax) then
          i=0
          do it=1,im_gr
             if (proc(it) == myidsp) then
                i=i+1
                atcinr%ax(:,i)=atcinr%xp(:,i)
             endif
          end do

       end if

    end select


    if (icintypemod==1) then
       read (lucin, err=444) oldtstep
       if (lrestart) then
          read (lucin, err=445) tmean, pmean, iteration, timel
          if (nitmax.ge.0) itmax=iteration+nitmax
          tstep = oldtstep

          if ((rang==0).and.(lprt)) then

             write (6, *) 'restart parameters'
             write (6, *) 'it =', iteration, ' time =', timel
             write (6, *) 'pmean', pmean, ' tmean =', tmean
             write (6, *) 'tstep', tstep
          endif                                ! fin rang=0
       end if
       usdh = 1.0/(two*tstep)
    endif



    close (lucin)
    return
431 print *,'Erreur 431'
432 print *,'Erreur 432'
433 print *,'Erreur 433'
434 print *,'Erreur 434'
435 print *,'Erreur 435'
436 print *,'Erreur 436'
437 print *,'Erreur 437'
438 print *,'Erreur 438'
439 print *,'Erreur 439'
440 print *,'Erreur 440'
441 print *,'Erreur 441'
442 print *,'Erreur 442'
443 print *,'Erreur 443'
444 print *,'Erreur 444'
445 print *,'Erreur 445'
    call arret_ndm




  end subroutine read_cin2

  subroutine read_cin(boxcin,itread,atcinr,immr,fnamcin,lres,fmtcin)
    !itread 0=at seulement; 1=complet; 2 = at, xp et num_at_glob seulement
    !immr : imm extrait .cin
    !lres : lrestart,
    !fmtcin=1 avec num_at_glob (optional)
    !icible tableau de taille imic qui donne les atomes à lire (utile pour para), optionel
    USE T_kind_param_m, ONLY:  double
    USE gen_com_m,only: iteration,itmax,nitmax,pmean,oldtstep,timel,two,usdh,dilat,tmean,tstep
    implicit none
    character,intent(in) :: fnamcin*80
    integer,intent(in)::itread
    class(box_config)::boxcin
    class(atom_config),optional::atcinr
    integer,intent(in),optional::immr,fmtcin
    logical,intent(in),optional::lres
    !    integer,allocatable,optional,intent(in)::icible(:)
    !    integer,optional,intent(in)::imic

    logical::lrestart=.false.
    integer :: i, ic, icintype, icintypemod , lucin,fmt_cin=1
    integer, dimension(:),allocatable     :: ibuffer
    real(double), dimension(:,:),allocatable    :: buffer
    real(double)::at(3,3)
    !      integer , dimension(imm,ntyp) :: fv
    !         integer , dimension(6000,10) :: fv    !Truc_bizarre_jmd
    integer::im_gr,i_loc
    if (present(immr))then
       allocate (ibuffer(immr))      ; allocate (buffer(3,immr))
    endif
    if (present(lres))lrestart=lres
    if (present(fmtcin))fmt_cin=fmtcin
    if ((rang==0).and.(lprt)) then
       write(6,*)
       write(6,*)' *-*-*-*-*-*READING OF CIN FILE*-*-*-*-*-*-'
       write(6,*)' *-*-*-*-*- LRESTART =',lrestart!, '*** itread',itread
    endif

    lucin = 93
    open(unit=lucin, file=fnamcin, form='unformatted', status='unknown', err=456)

    read (lucin, err=456) icintype

    if ((rang==0).and.(lprt))  write (6, *) 'config type of  .cin file : ', icintype
    if (icintype>5.or.icintype<0) then
       write (6, *) rang, 'wrong icintype'
       call arret_ndm
    endif

    icintypemod = mod(icintype,2)

!!$    if (lrestart.and.icintypemod==0) then
!!$       write (6, *) rang, 'not possible to restart from this file'
!!$       call arret_ndm
!!$    endif
    !at(vect123,xyz)
    !        if (icintype>=2) then
    read (lucin, err=456) at
    if(dilat(1).ne.0.0)then
       do i=1,3
          at(i,:)=at(i,:)*dilat(i)
       end do
    end if

    call boxcin%init(at,ipbc)

    select case(itread)
    case(0)
       close (lucin)
       return
    case(1)
       if (.not.present(atcinr))then
          write(6,*)'atcinr pas present et itread=1'
          call arret_ndm
       end if

       read (lucin, err=456) im_gr                         !number of atoms in the box
       if (im_gr>immr) then
          if(rang==0)                    write (6, *) 'P1 im > imM', im_gr, immr
          call arret_ndm
       endif
       atcinr%im=im_gr
       atcinr%im_glob=im_gr


       read (lucin, err=456) ibuffer   !ityp
       atcinr%ityp(1:im_gr)=ibuffer(1:im_gr)

       if ((rang==0).and.(lprt))  write (6, *) 'types'
       read (lucin, err=456) buffer    ! xp
       atcinr%xp(1:3,1:im_gr)=buffer(1:3,1:im_gr)
       if ((rang==0).and.(lprt))  write (6, *) 'xp'
       formcin:select case (fmt_cin)
       case (0) formcin
          do i=1,im_gr
             atcinr%num_at_glob(i) = i
          enddo
       case(1) formcin
          read (lucin, err=456) ibuffer
          atcinr%num_at_glob(1:im_gr)=ibuffer(1:im_gr)
          if ((rang==0).and.(lprt))  write (6, *) 'num_at_glob'
       case default  formcin
          if (rang.eq.0) write(6,*) 'precisez le format fmt_cin'
          call arret_ndm
       end select formcin


       select type(atcinr)
       type is (atom_config)
          if (icintypemod==1) then
             if (icintype==3)read (lucin, err=456) buffer                     !xpp
             read (lucin, err=456) buffer                     !vp
          end if
          lvpread=.false.
       type is (atom_config_d)
          if (icintypemod==1) then
             if (icintype==3) read (lucin, err=456) buffer                     !xpp
             read (lucin, err=456) buffer                     !vp
             atcinr%vp(:,1:im_gr)=buffer(:,1:im_gr)
             !             read (lucin, err=456) buffer                     !former positions
             !             atcinr%ax(:,1:im_gr)=buffer(:,1:im_gr)
             !             read (lucin, err=456) buffer                     !ax inutile
          end if
       end select
       select type(atcinr)
       class is (atom_config_e)
          if (icintypemod==1) then

             if (icintype==3) then
                read (lucin, err=456) buffer                     !xpp
             end if

!!$             read (lucin, err=456) buffer                     !xpp
!!$             atcinr%xpp(:,1:im_gr)=buffer(:,1:im_gr)
!!$             write(6,*)'xpp_e'
             read (lucin, err=456) buffer                     !vp
             atcinr%vp(:,1:im_gr)=buffer(:,1:im_gr)
             !             write(6,*)'vp_e'
             !             read (lucin, err=456) buffer                     !former positions
             !             atcinr%ax(:,1:im_gr)=buffer(:,1:im_gr)

             !             read (lucin, err=456) buffer                     !ax utile peut-être
             !             if ((icintype==3).or.(icintype==7)) then 
             !                if (lrestart) then 
             !                   if (atcinr%lax) then
             !                      atcinr%ax(:,1:im_gr)=buffer(:,1:im_gr)
             !                  end if
             !               end if
             !                  end if
          else
          end if
          if (atcinr%lax) then
             atcinr%ax(:,1:im_gr)=atcinr%xp(:,1:im_gr)

          end if

       end select


       if (icintypemod==1) then
          read (lucin, err=456) oldtstep
          if (lrestart) then
             read (lucin, err=456) tmean, pmean, iteration, timel
             if (nitmax.ge.0) itmax=iteration+nitmax
             tstep = oldtstep

             if ((rang==0).and.(lprt)) then

                write (6, *) 'restart parameters'
                write (6, *) 'it =', iteration, ' time =', timel
                write (6, *) 'pmean', pmean, ' tmean =', tmean
                write (6, *) 'tstep', tstep
             endif                                ! fin rang=0
          end if
          usdh = 1.0/(two*tstep)
       endif
    case(2) ! at xp et num_at_glob
       if (.not.present(atcinr))then
          write(6,*)'atcinr pas present et itread=2'
          call arret_ndm
       end if

       read (lucin, err=456) im_gr                         !number of atoms in the box
       if (im_gr>immr) then
          if(rang==0)                    write (6, *) 'P2 im > imM', im_gr, immr
          call arret_ndm
       endif
       call atcinr%init(im_gr,immr,im_glob=im_gr)

       read (lucin, err=456) ibuffer   !ityp muet
       read (lucin, err=456) buffer    ! xp
       atcinr%xp(1:3,1:im_gr)=buffer(1:3,1:im_gr)

       read (lucin, err=456) ibuffer
       atcinr%num_at_glob(1:im_gr)=ibuffer(1:im_gr)


    end select
    call cryst_to_cart (atcinr%im, atcinr%xp, boxcin%bg, -1) !cart vers cryst
    do ic=1,3
       if ((lperiod).or.(ipbc(ic).ne.1)) then
          do i=1,atcinr%im
             if ( (atcinr%xp(ic,i).LT.0.d0).OR.(atcinr%xp(ic,i).GE.1.d0) ) then
                atcinr%xp(ic,i)  = atcinr%xp(ic,i)  - Dble(Floor(atcinr%xp(ic,i)))
             END if
          end do
       end if
    end do
    call cryst_to_cart (atcinr%im, atcinr%xp, boxcin%at, 1) ! cryst vers cart

    close (lucin)
    return
456 print *,'Erreur dans la lecture du fichier .cin, verifier son format&
         & et fmt_cin ATTENTION A BIG_ENDIAN !! SI COMMPILE BIG_ENDIAN NE LIT PLUS QUE CA'


  end subroutine read_cin
  !******************************************************************************************************
  !******************************************************************************************************
  !******************************************************************************************************

  subroutine read_gin (boxrg,atrg,fnamgin,latr,itread,immread)
    USE T_kind_param_m, ONLY:  double




    character,intent(in) :: fnamgin*80
    type(atom_config),intent(out)::atrg
    type(box_config),intent(out)::boxrg
    integer,intent(out)::latr(3)
    integer,optional,intent(in)::itread,immread
    integer::immr

    integer::itr=1
    !    real(double)::rumax_init,alpha_init
    real(double)::at(3,3),deltx

    integer ::  lugin, imcell, ic,i

    if(present(itread))itr=itread
    if(present(immread))then
       immr=immread
    else
       immr=imm_glob
    end if
    !  si coordonnees reduites

    if ((rang==0).and.(lprt))  write (6, *) '**********construction du reseau************'
    lugin=92
    !                                                !number of cells in 3 directions
    open(unit=lugin, file=fnamgin, status='unknown')
    read (lugin, *) latr(1), latr(2), latr(3)
    if ((rang==0).and.(lprt))  write (6, *) 'repetition de mailles', latr


    !     **** coordonnes des vecteurs de maille en A dans une base orthonormee ****
    !                                                !a
    read (lugin, *) at(1,1), at(2,1), at(3,1)
    !b
    read (lugin, *) at(1,2), at(2,2), at(3,2)
    !                                                !c
    read (lugin, *) at(1,3), at(2,3), at(3,3)
    at=at*1d-8
    call boxrg%init(at,ipbc)
    read (lugin, *) imcell               !number of atoms in UC
    if (itr==0) return
    if (imcell>immr) then
       if ((rang==0).and.(lprt)) write (6, *) 'trop d_atomes dans la cel. unite',immr,imcell
       call arret_ndm
    endif
    call atrg%init(imcell)
    do i = 1, imcell
       read (lugin, *) atrg%xp(1,i), atrg%xp(2,i), atrg%xp(3,i),atrg%ityp(i)
    end do
    if (ldecalcor) then
       if (any(atrg%xp (1:3,1:imcell)==0)) then
          if ((rang==0).and.(lprt))  write(6,*)' .gin with 0 coordinates; creates FAILURES,  POSITIONS SHIFTED By +2e-7'
          atrg%xp (1:3,1:imcell)=atrg%xp (1:3,1:imcell)+2e-7
       end if
       if (any(atrg%xp (1:3,1:imcell)==1)) then
          if ((rang==0).and.(lprt))  write(6,*)' .gin with 1 coordinates; creates FAILURES,  POSITIONS SHIFTED By -1e-7'
          atrg%xp (1:3,1:imcell)=atrg%xp (1:3,1:imcell)-1e-7
       end if
    end if

    if (dmtype==9) then
       if (lsecondpath) then
          !write(6,*)xpd
          do i=1,imcell
             do ic=1,3
                deltx=atrg%xp(ic,i)-xpd(ic,i)
                if (deltx.gt.0.5) atrg%xp(ic,i)=atrg%xp(ic,i)-1.
                if (deltx.lt.-0.5) atrg%xp(ic,i)=atrg%xp(ic,i)+1.
             end do
          end do
       else
          allocate (xpd(3,imcell))
          xpd=atrg%xp
       end if
    end if


    do ic=1,3
       if ((lperiod).or.(ipbc(ic).ne.1)) then
          do i=1,imcell
             if ( (atrg%xp(ic,i).LT.0.d0).OR.(atrg%xp(ic,i).GE.1.d0) ) then
                atrg%xp(ic,i)  = atrg%xp(ic,i)  - Dble(Floor(atrg%xp(ic,i)))
             END if
          end do
       end if
    end do
    close(lugin)
  end subroutine read_gin

#ifdef DKIO
  subroutine dkio2ndm(at2b,cel2b,box2b,fnam,rum,format,lrepartition,psc,lconstrsimple,immread)
    !-----------------------------------------------------
    !  Subroutine for interfacing with the dk_io library
    !-----------------------------------------------------
    use dk_structure_io, only: read_structure, TAG_LENGTH
    USE decoupage_mod,only: constrandrepart_dkio

    type(para_space_config),optional::psc
    class(atom_config)::at2b
    type(cell_config)::cel2b
    class(box_config)::box2b
    character,intent(in) :: fnam*80,format*6
    real(double),intent(in)::rum
    logical,optional,intent(in)::lrepartition,lconstrsimple
    integer,optional::immread
    integer::immr,npr
    logical::lcs,lrepart
    real(double) :: boxrin(3,3),deltx
    real(double), dimension(:,:), allocatable :: atrin,vpin
    character(TAG_LENGTH), dimension(:), allocatable :: tags
    integer::i,ic,ncore,itread,imcell
    lrepart=.true.
    lcs=.false.
    if(present(lrepartition))lrepart=lrepartition
    if(present(lconstrsimple))lcs=lconstrsimple
    immr=imm_glob
    if (present(immread)) immr=immread

    if (ldecoup) then
       itread=0 ! Attention : itread=0 pas utilisé avec dk-io => tous les atomes seront lus malgrés itread=0
    else
       itread=1
    end if
    !call read_gin(boxrgin,atrgin,fnamg,lat,itread=itread,immread=immr)
    if (lvpread) then
       call read_structure(trim(fnam),boxrin,atrin,tags,format=format, velocities=vpin)
    else
       call read_structure(trim(fnam),boxrin,atrin,tags,format=format)
    end if
    imcell=size(atrin,2)

    if (ldecalcor) then
       if (any(atrin(1:3,1:imcell)==0)) then
          if ((rang==0).and.(lprt))  write(6,*)' atom configuration file with 0 coordinates; creates FAILURES,  POSITIONS SHIFTED By +2e-7'
          atrin(1:3,1:imcell)=atrin(1:3,1:imcell)+2e-7
       end if
       if (any(atrin(1:3,1:imcell)==1)) then
          if ((rang==0).and.(lprt))  write(6,*)' atom configuration file with 1 coordinates; creates FAILURES,  POSITIONS SHIFTED By -1e-7'
          atrin(1:3,1:imcell)=atrin(1:3,1:imcell)-1e-7
       end if
    end if

    if (dmtype==9) then
       if (lsecondpath) then
          !write(6,*)xpd
          do i=1,imcell
             do ic=1,3
                deltx=atrin(ic,i)-xpd(ic,i)
                if (deltx.gt.0.5) atrin(ic,i)=atrin(ic,i)-1.
                if (deltx.lt.-0.5) atrin(ic,i)=atrin(ic,i)+1.
             end do
          end do
       else
          allocate (xpd(3,imcell))
          xpd=atrin
       end if
    end if  
    
    do ic=1,3
       if ((lperiod).or.(ipbc(ic).ne.1)) then
          do i=1,imcell
             if ( (atrin(ic,i).LT.0.d0).OR.(atrin(ic,i).GE.1.d0) ) then
                atrin(ic,i)  = atrin(ic,i)  - Dble(Floor(atrin(ic,i)))
             END if
          end do
       end if
    end do

    boxrin=boxrin*1d-8
    call box2b%init(boxrin,ipbc)
    
    call setnox(box2b,cel2b,rum,lverbose=lprt,noxr=nox,noyr=noy,nozr=noz)
    if ((rang==0).and.(lprt)) then
       write (6, '(2A,D15.8,A,D15.8,A)') fnam,'volume=', box2b%volu,' cm3 ',box2b%volu*1d24,' Ang3'
    end if

    if (lcs) then ! construction simpple sans repartition en sequentiel
       call constr_2dkio (at2b,atrin,tags,vpin,immr)
       call cryst_to_cart (at2b%imm, at2b%xp, box2b%at, 1)
       at2b%im_glob=at2b%im
       call setcellconf(cel2b,at2b,box2b,rum)
       return
    end if
    
#ifdef PARA
    if (rang==0) then
       if (ldecoup) then
          open(123, file='decoup.dat', status='old')
          read (123, *) npr,ncore
          close(123)         
          call  decoupage(npr,ncore,cel2b,psc=psc,lverbose=lprt)
          call arret_ndm
       end if
    end if

    call  decoupage(nprocspace,ncore,cel2b,psc=psc,lverbose=lprt)
    ncore=0
    at2b%imm_glob=imm_glob
    at2b%im_glob=imcell
    
    if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.).and.(lrepart)) then
       call constrandrepart_dkio(atrin,vpin,at2b,cel2b,box2b,tags,psc)
    else
       call constr_2dkio(at2b,atrin,tags,vpin,imm_glob)
       call cryst_to_cart (at2b%imm, at2b%xp, box2b%at, 1)
       at2b%imm_glob=imm_glob
    end if
#else
    if (ldecoup) then
       open(123, file='decoup.dat', status='old')
       read (123, *) npr,ncore
       close(123)         
       call  decoupage(npr,ncore,cel2b,psc=psc,lverbose=lprt)
       call arret_ndm
    end if
    call constr_2dkio (at2b,atrin,tags,vpin,imm)
    call cryst_to_cart (at2b%imm, at2b%xp, box2b%at, 1)
    at2b%im_glob=at2b%im
#endif
    call setcellconf(cel2b,at2b,box2b,rum)
    if (allocated(atrin)) deallocate(atrin)
    if (allocated(tags)) deallocate(tags)
    return
  end subroutine dkio2ndm

  subroutine constr_2dkio(atrcf,atrin,tags,vpin,immread)
    use dk_structure_io, only: TAG_LENGTH
    use decoupage_mod,only: get_ityp

    class(atom_config),intent(inout)::atrcf 
    real(double), dimension(:,:), intent(in) :: atrin,vpin
    character(TAG_LENGTH), dimension(:), intent(in) :: tags
    integer,intent(in),optional::immread
    integer::i,imloc,immr
    real(double)::rvn

    imloc=size(atrin, 2)
    immr=imm_glob
    if (present(immread)) immr=immread
    if (imloc>immread) then
       write (6, *) rang,'imm trop petit',imloc,immr
       call arret_ndm
    endif

    if (atrcf%rvois.gT.0)then
       rvn=atrcf%rvois
    else
       rvn=0
    end if
    if (present(immread))then
       call atrcf%init(imloc,immin=immread,ltabvois=atrcf%ltabvois,nvois=atrcf%nvois,rvois=rvn,im_glob=imloc)
    else
       call atrcf%init(imloc,ltabvois=atrcf%ltabvois,nvois=atrcf%nvois,rvois=rvn,im_glob=imloc)
    end if
    

    do i = 1, imloc
       atrcf%xp(1,i) = atrin(1,i)
       atrcf%xp(2,i) = atrin(2,i)
       atrcf%xp(3,i) = atrin(3,i)
       atrcf%num_at_glob(i)=i
       atrcf%ityp(i)=get_ityp(tags(i))
    end do

    if (lvpread) then
       select type (atrcf)
       type is(atom_config)
          if (rang==0) write(6,*)'no velocity in atom-config and import asked with velocities stop'
          call arret_ndm
       class is (atom_config_d)
          do i = 1, imloc
             atrcf%vp(1,i) = vpin(1,i)*1d4 ! car *1d8*1d-12 à l'écriture
             atrcf%vp(2,i) = vpin(2,i)*1d4
             atrcf%vp(3,i) = vpin(3,i)*1d4
          end do
       end select
    end if

    return
  end subroutine constr_2dkio

#endif
  !#endif
end module constrconf_mod
