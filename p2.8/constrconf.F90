module constrconf_mod
!    USE period_mod,only:period
  !  USE divid_mod, only:divid
#ifdef PARA
  USE decoupage_mod,only: decoupage
#endif
  USE gen_com_m, ONLY: lenfnam, fnam,fmt_cin,igen,im_glob,imm,imm_glob,ldecoup,lperiod,lrestart,rang,&
       &lvpread,zero,low_limit !at,bg,zls2,tstep,oldtstep,tmean,timel,nox,noy,noz,im,imm,&
  !       &it,itmax,ldesinteg,lperiod,pmean,zl,xpspr,nzl,normat,cell_debx,cell_deby,cell_debz,&
  !       &cell_finx,cell_finy,cell_finz,low_limit,llangevin,lsuivinonpbc
  USE var_pot, ONLY:na,ntyp,rumax,ipotentiel
    use cryst_to_cart_mod,only:cryst_to_cart
    USE arret_ndm_mod,only: arret_ndm
    USE atomconfig,only:atom_config,atom_config_d,atom_config_e
  USE cellconfig,only:cell_config
  USE boxconfig,only:box_config,initbox,periodbox
!  USE read_conf,only:read_cin,read_gin
  USE setcell,only:setnox,setcellconf
  USE decoupage_mod,only: decoupage
  USE T_kind_param_m, ONLY:  double
#ifdef PARA
    use mpi
    USE mod_para,only:MPI_COMM_space,ierr,NDM_MPI_REAL_DOUBLE,status,nprocs
#endif


  implicit none
contains
  subroutine constrconf (atrcf,boxrcf,cellrcf)
    !********************************************************************
    !             CONSTRUCTION DE LA BOITE DE SIMULATION
    !********************************************************************

    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------

    !    USE tab_imm_m,only:xp,ax,glangv,vp,xpp
    !    USE suivinonpbc
    !    USE read_conf_mod

    implicit none

    class(atom_config),intent(inout)::atrcf
    type(cell_config),intent(out)::cellrcf
    type(box_config),intent(out)::boxrcf
    integer :: i, icell, iti

    type(box_config)::boxrgin
    type(atom_config)::atrgin
#ifndef PARA
    integer :: nprocs
#endif
    integer,      dimension(:), allocatable   :: num_at_buff
    integer, dimension(:),allocatable     :: ibuffer
    real(double), dimension(:,:),allocatable    :: buffer
    integer::ncore,ic

#ifdef PARA
    type (atom_config)::COMPatrcf
#endif
    character :: fnamcin*80, fnamgin*80
    integer::itread

    
    !-----------------------------------------------------
    ! READING FROM THE CONFIGURATION FILE
    !---------------------------------------------------

    if (rang==0) then
       write(6,*)
       write(6,*)' *-*-*-*-*-*CONSTRUCTION DE LA BOITE*-*-*-*-*-*-',imm,imm_glob,nprocs
       write(6,*)
    endif

    if (igen.ge.1) then
    allocate (ibuffer(imm_glob))
    allocate (buffer(3,imm_glob))
       if(rang==0)               write(6,*)'********** reading configuration from file********'
       if (lrestart) then
          fnamcin = fnam(1:lenfnam)//'.cout'
       else
          fnamcin = fnam(1:lenfnam)//'.cin'
       end if

#ifdef PARA
       ! En parallele, la lecture du fichier de position se fait en passes
       !  - la premiere pour lire toutes les positions et determiner le
       !    meilleur equilibrage/decoupage
       !  - la deuxieme pour lire uniquement les positions propres au
       !    processeur

       itread=2
       call read_cin(boxrcf,itread,COMPatrcf,imm_glob,fnamcin,lrestart,fmt_cin) !0=at seulement; 1=complet; 2 = at, xp et num_at_glob seulement
       if (rang==0) then
          write (6, '(A,D15.8,A,D15.8,A)') 'volume=', boxrcf%volu,' cm3 ',boxrcf%volu*1d24,' Ang3'
       end if
       call setnox(boxrcf,cellrcf,rumax)
       do iti=1,ntyp
          na(iti)=count(COMPatrcf%ityp(1:COMPatrcf%im)==iti)
       end do
       ncore=0
       call  decoupage(nprocs,ncore,cellrcf,atrcf)
       allocate(num_at_buff(imm_glob))
       im_glob=COMPatrcf%im
       
       call repartition(COMPatrcf,atrcf,boxrcf,cellrcf,num_at_buff)
       itread=3
       call read_cin(boxrcf,itread,atrcf,imm_glob,fnamcin,lrestart,fmt_cin,num_at_buff,imm_glob) !0=at seulement; 1=complet; 2 = at, xp et num_at_glob seulement , 3 num_at_buff masque des atomes locaux

#else


       if (ldecoup) then 
          itread=0
          call read_cin(boxrcf,itread,fnamcin=fnamcin)
          if (rang==0) then
             write (6, '(A,D15.8,A,D15.8,A)') 'volume=', boxrcf%volu,' cm3 ',boxrcf%volu*1d24,' Ang3'
          end if
          call setnox(boxrcf,cellrcf,rumax)
          open(123, file='decoup.dat', status='old')
          read (123, *) nprocs,ncore
          close(123)         
          call  decoupage(nprocs,ncore,cellrcf)
          stop
       else
          itread=1
          call read_cin(boxrcf,itread,atrcf,imm,fnamcin,lrestart,fmt_cin) !0=at seulement; 1=complet; 2 = at, xp et num_at_glob seulement , 3 trié par num_at_buff
          im_glob=atrcf%im
          if (rang==0) then
             write (6, '(A,D15.8,A,D15.8,A)') 'volume=', boxrcf%volu,' cm3 ',boxrcf%volu*1d24,' Ang3'
          end if
          call setnox(boxrcf,cellrcf,rumax)
          !          CALL fin allocation CELL et FIN DIVID
          do iti=1,ntyp
             na(iti)=count(atrcf%ityp(1:atrcf%im).eq.iti)
          end do
       end if

#endif
       call setcellconf(cellrcf,atrcf,boxrcf,im_glob,rumax)
    deallocate (ibuffer)
    deallocate (buffer)

       !-----------------------------------------------------
       ! BUILDING OF THE CRISTAL FROM .GIN FILE
       !-----------------------------------------------------
    else    ! igen.eq.0



       lvpread=.false.


       ! open fichier .gin
       fnamgin = fnam(1:lenfnam)//'.gin'

       call gin2ndm(atrcf,cellrcf,boxrcf,fnamgin,im_glob,rumax)
       do iti=1,ntyp
          na(iti)=count(atrcf%ityp(1:atrcf%im)==iti)
       end do

       if (lperiod.EQV..true.) call periodbox (boxrcf,atrcf)

       select type(atrcf)
       type is (atom_config_d)
          atrcf%xpp(:,1:atrcf%im)=atrcf%xp(:,1:atrcf%im)
       type is (atom_config_e)
          atrcf%xpp(:,1:atrcf%im)=atrcf%xp(:,1:atrcf%im)
          if (atrcf%lax) then
             atrcf%ax(:,1:atrcf%im)=atrcf%xp(:,1:atrcf%im)
          end if
       end select

    end if

    

    if (rang==0) then

       write(6,*)
       write (6, *) '-------- boite de simulation ------'
       write (6, *) 'nombre d atomes =', im_glob
       do i=1,3
          write(6,'(A,I2,3F15.6)')'vecteur ',i, (boxrcf%at(ic,i)*1.0d8,ic=1,3)
       end do
       do iti = 1, ntyp
          if (na(iti)==0) cycle
          write (6, *) na(iti), ' atomes de type', iti
       end do
    endif                                  ! fin rang=0

    ! ----------------------------------------------------------
    !  CONDITIONS PERIODIQUES : REMETTRE LES ATOMES DANS BOITE
    ! ----------------------------------------------------------




    ! SUMMARY
    !#ifdef LAMMPS_VERSION

    if((ipotentiel==-10).or.(ipotentiel==-11)) then
       if(rang==0)write(6,*)'write configuration to conf.lmp'
       call config2data (imm,atrcf%im,atrcf%xp,atrcf%ityp,boxrcf%at,ntyp)
    end if
    !#endif     


    !       end if




    return

456 print *,'Erreur dans la lecture du fichier .cin, verifier son format&
         & et fmt_cin ATTENTION A BIG_ENDIAN !! SI COMMPILE BIG_ENDIAN NE LIT PLUS QUE CA'


  end subroutine constrconf
  
  subroutine coord_to_cell(tab_coord, cell,bg,nox,noy,noz)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double
    

    implicit none

    !       version du 10 janvier 2007
    !
    ! Cette routine retourne dans cell le numero de cellule
    ! contenant les coordonnees tab_coord
    !
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    real(double)  :: tab_coord(3),bg(3,3)
    integer       :: cell,nox,noy,noz
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: kx, ky, kz
    real(double) :: aux, auy, auz
    real(double) :: coord_loc(3)   ! permet de ne pas ecraser coord
    ! lors de l'appel a cryst_to_cart                       
    !-----------------------------------------------

    coord_loc(:) = tab_coord(:)
    call cryst_to_cart (1, coord_loc, bg, -1) !cart vers cryst

    aux = coord_loc(1)*nox
    auy = coord_loc(2)*noy
    auz = coord_loc(3)*noz

    kx = int(aux)
    ky = int(auy)
    kz = int(auz)

    cell = 1+kx+nox*(ky+noy*kz)

    return
  end subroutine coord_to_cell

  subroutine constr_2gin(atrcf,boxrcf,cellrcf,atrgin,boxrgin,lat,imm)

    class(atom_config),intent(inout)::atrcf
    type(cell_config),intent(in)::cellrcf
    type(box_config),intent(in)::boxrcf    
    type(box_config)::boxrgin
    type(atom_config)::atrgin
    integer,intent(in)::lat(3)
    integer,intent(in),optional::imm

    integer::i,ia,ib,ic,icell,imloc,ncore

    !imtot=lat(1)*lat(2)*lat(3)*atrgin%im
    imloc=lat(1)*lat(2)*lat(3)*atrgin%im
    if (imloc>imm_glob) then
       write (6, *) rang,'imm trop petit',imloc,imm_glob
       call arret_ndm
    endif

    if (present(imm))then

       call atrcf%init(imloc,immin=imm,ltabvois=atrcf%ltabvois,nvois=atrcf%nvois)
       
    else
       call atrcf%init(imloc,ltabvois=atrcf%ltabvois,nvois=atrcf%nvois)
    end if
    
!    imtot=imloc
    i=0
    write(6,*)
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
                !                   write(6,*)'constr',i,atrcf%xp(:,i)
             end do
          end do
       end do
    end do

    return
  end subroutine constr_2gin

  subroutine repartition(atcomp,atrep,boxrep,cellrep,nab,div)
#ifdef PARA
    use mod_para,only:myid
#endif
    use paraconfig,only:para_config
    class(atom_config)::atrep
    class(atom_config),intent(in)::atcomp
    type(cell_config)::cellrep
    type(box_config)::boxrep
    integer,optional, dimension(:), allocatable   :: nab
    type (para_config),optional::div
    
    integer::i,icomp,k,iti,im,ic,numcell,numproc
    real(double)::xt(3),xpici,cpp
    
#ifdef PARA  


       if (atrep%im>imm_glob) then
          write (6, *) rang,'imm trop petit'
          call arret_ndm
       endif

       i=0;im=0
       do icomp=1,atcomp%im
          xt(:)=atcomp%xp(:,icomp)
          do k=1,3
             xpici=xt(k)
             if ( (xpici < 0.d0 ).OR.( xpici >= 1.d0 ) ) then
                if ( (xpici > -low_limit).and.(xpici<0.d0) ) then
                   xt(k)=zero
                else
                   cpp  = Dble(Floor(atcomp%xp(ic,i)))
                   xt(k) = xpici     - cpp
                end if
             end if
          end do
          iti = atcomp%ityp(icomp)
          call coord_to_cell(xt,numcell,boxrep%bg,cellrep%nox,cellrep%noy,cellrep%noz)
          numproc=cellrep%proc_cell(numcell)
          if (numproc == myid) then
             i=i+1
             im=im+1
             atrep%num_at_glob(i)=atcomp%num_at_glob(icomp)
             if (present(nab)) nab(i)=icomp
             atrep%xp(:,i)=xt(:)
             atrep%ityp(i)=iti
             div%proc_at(i)=myid
          endif
       end do
       atrep%im=im
#endif
       return
     end subroutine repartition
          
  
  subroutine config2data (imm,im,xp,ityp,at,ntyp)
    USE T_kind_param_m, ONLY:  double
    USE gen_com_m, ONLY : position_conversion_lammps
    USE var_pot, ONLY:q,ipotentiel
    USE Mat_utils_mod,only: Matinv,is_upper_triangular
    implicit none
    integer,intent(in)::imm,im,ntyp
    real(double),intent(in)::xp(3,imm),at(3,3)
    integer,intent(in)::ityp(imm)

    logical lrotated,upper
    real(double)::xhi,yhi,zhi,xy,xz,yz,xlo,ylo,zlo
    real(double), dimension(3,3)::at_lammps,passage,passage_inv
    integer::ic,i
    real(double), dimension(3) :: tmp_coord_i,new_tmp_coord_i

    !  if ((at(2,1).ne.0).or.(at(3,1).ne.0).or.(at(3,2).ne.0))then
    write(6,*)
    write(6,*)'LAMMPS BOX BUILT FROM NDM'
    write(6,*)'transforms NDM box to lamps shaped box if needed'
    write(6,*)
    !     stop
    ! endif

    open(63,file='conf.lmp',status='unknown')
    write(63,*)

    call is_upper_triangular(at,upper)
    write(6,*)'upper ? ',upper
    if (.not.upper) then
       call convert_cell (at,at_lammps,passage)
       call matinv(passage, passage_inv)
       do ic=1,3
          write(6,*)passage(:,ic)
       end do
       write(6,*)
       do ic=1,3
          write(6,*)passage_inv(:,ic)
       end do
    else
       at_lammps=at
       passage(:,:)=0
       do ic=1,3
          passage(ic,ic)=1
       end do
       passage_inv(:,:)=passage(:,:)
    end if


    !building  lammps header
    xlo = 0.d0
    ylo = 0.d0
    zlo = 0.d0
    xhi = at_lammps(1,1)/position_conversion_lammps
    yhi = at_lammps(2,2)/position_conversion_lammps
    zhi = at_lammps(3,3)/position_conversion_lammps
    xy = at_lammps(1,2)/position_conversion_lammps
    xz = at_lammps(1,3)/position_conversion_lammps
    yz = at_lammps(2,3)/position_conversion_lammps

    !  xhi=at(1,1)/position_conversion_lammps
    !  yhi=dsqrt(at(1,2)**2+at(2,2)**2)/position_conversion_lammps
    !  zhi=dsqrt(at(1,3)**2+at(2,3)**2+at(3,3)**2)/position_conversion_lammps
    !  xz=at(1,1)*at(1,3)/(xhi*zhi*position_conversion_lammps*position_conversion_lammps)
    !  xy=at(1,1)*at(1,2)/(xhi*yhi*position_conversion_lammps*position_conversion_lammps)
    !  yz=(at(1,2)*at(1,3)+at(2,2)*at(2,3))/(yhi*zhi*position_conversion_lammps*position_conversion_lammps)

    if(IM.lt.10)then
       write(63,"(I1,A)") IM,' atoms'
    elseif((IM.lt.100).and.(IM.gt.10))then
       write(63,"(I2,A)") IM,' atoms'
    elseif((IM.lt.1000).and.(IM.gt.100))then
       write(63,"(I3,A)") IM,' atoms' 
    elseif((IM.lt.10000).and.(IM.gt.1000))then
       write(63,"(I4,A)") IM,' atoms'       
    elseif((IM.lt.100000).and.(IM.gt.10000))then
       write(63,"(I5,A)") IM,' atoms'         
    elseif((IM.lt.1000000).and.(IM.gt.100000))then
       write(63,"(I6,A)") IM,' atoms'         
    elseif((IM.lt.10000000).and.(IM.gt.1000000))then
       write(63,"(I7,A)") IM,' atoms'        
    elseif((IM.lt.100000000).and.(IM.gt.10000000))then
       write(63,"(I8,A)") IM,' atoms'         
    else 
       stop 'ADD FORMAT'
    endif
    write(6,*)
    write(6,*) " a = (xhi-xlo,0,0); b = (xy,yhi-ylo,0); c = (xz,yz,zhi-zlo). "
    write(6,'(f22.16,a,f22.16,a)')xlo,' ',xhi,' xlo xhi'
    write(6,'(f22.16,a,f22.16,a)')ylo,' ',yhi,' ylo yhi'
    write(6,'(f22.16,a,f22.16,a)')zlo,' ',zhi,' zlo zhi'
    write(6,'(f22.16,a,f22.16,a,f22.16,a)')xy,' ',xz,' ',yz,' xy xz yz'
    write(6,*)


    write(63,"(I1,A)") ntyp, ' atom types' 
    write(63,*) 
    write(63,'(f22.16,a,f22.16,a)')xlo,' ',xhi,' xlo xhi'
    write(63,'(f22.16,a,f22.16,a)')ylo,' ',yhi,' ylo yhi'
    write(63,'(f22.16,a,f22.16,a)')zlo,' ',zhi,' zlo zhi'
    write(63,'(f22.16,a,f22.16,a,f22.16,a)')xy,' ',xz,' ',yz,' xy xz yz'
    write(63,*)
    write(63,"(A)")'Atoms'
    write(63,*)
    select case (ipotentiel)
    case(-10)
       do i=1,im
          tmp_coord_i = xp(:,i)
          new_tmp_coord_i = matmul(passage,tmp_coord_i)/position_conversion_lammps
          write(63,'(I8,a,I2,a,f22.15,a,f22.15,a,f22.15)')i,' ',ityp(i),' ',new_tmp_coord_i(1),' '&
               &,new_tmp_coord_i(2),' ',new_tmp_coord_i(3)

          !        write(63,"(I8,I6,F21.12,F20.12,F20.12)") i,ityp(i),&
          !             & xp(1,i)/position_conversion_lammps,xp(2,i)/position_conversion_lammps,xp(3,i)/position_conversion_lammps
       end do
    case(-11)
       do i=1,im
          tmp_coord_i = xp(:,i)
          new_tmp_coord_i = matmul(passage,tmp_coord_i)/position_conversion_lammps
          write(63,'(I8,a,I2,a,F20.12,a,f22.15,a,f22.15,a,f22.15)')i,' ',ityp(i),' ',q(ityp(i)),&
               &' ',new_tmp_coord_i(1),' ',new_tmp_coord_i(2),' ',new_tmp_coord_i(3)
          !        write(63,"(I8,I6,F20.12,F20.12,F20.12,F20.12)") i,ityp(i),&
          !             & q(ityp(i)), xp(1,i)/position_conversion_lammps,xp(2,i)/position_conversion_lammps,xp(3,i)/position_conversion_lammps
       end do
    end select
    close (63)
  end subroutine config2data

  !---------------------------------------------------
  subroutine convert_cell(mat_ini,new_mat,transform)
    USE Mat_utils_mod,only: Matinv,norme,cross_product,is_upper_triangular,right_hand_basis


    implicit none
    real(kind(0.d0)), dimension(3,3), intent(in) :: mat_ini
    real(kind(0.d0)), dimension(3,3), intent(out) :: new_mat, transform

    !internal
    real(kind(0.d0)), dimension(3,3) :: transit_cell,inv_mat_ini
    real(kind(0.d0)), dimension(3) :: A,B,C, Ahat,AxBhat
    real(kind(0.d0)) :: volume
    logical :: upper, right
    integer :: i

    !matrix is already transpose
    transit_cell = mat_ini
    !write(*,*) 'transpose matrix is :',transit_cell

    call is_upper_triangular(transit_cell,upper)
    if (.not.upper) then
       ! rotate bases into triangular matrix
       new_mat(:,:) = 0.d0
       A = transit_cell(:,1)
       B = transit_cell(:,2)
       C = transit_cell(:,3)
       call right_hand_basis(A,B,C,right)

       if (.not.right) then
          write(*,*)"WARNING: your reper is not right handed."
          write(*,*)"WARNING: This is a critical issue. The LAMMPS results are wrong !!!!!"
          stop
       end if

       new_mat(1,1) = norme(A)
       Ahat = A / norme(A)
       AxBhat = cross_product(A, B) / norme(cross_product(A, B))
       new_mat(1,2) = dot_product(B, Ahat)
       new_mat(2,2) = norme(cross_product(Ahat, B))
       new_mat(1,3) = dot_product(C,Ahat)
       new_mat(2,3) = dot_product(C,cross_product(AxBhat, Ahat))
       new_mat(3,3) = abs(dot_product(C, AxBhat))
       !create and save the transformation for coordinates
       !volume = matdet(mat_ini)
       !trans = np.array([np.cross(B, C), np.cross(C, A), np.cross(A, B)])
       !trans = trans / volume
       !coord_transform = np.dot(tri_mat , trans)
       call matinv(mat_ini,inv_mat_ini)
       transform = matmul(new_mat,inv_mat_ini)

    else
       new_mat = mat_ini
       transform(:,:) = 0.d0
       do i=1,3
          transform(i,i) = 1.d0
       enddo

    endif
    return
  end subroutine convert_cell

  subroutine gin2ndm(at2b,cel2b,box2b,fnamg,imtot,rum)
    class(atom_config)::at2b
    type(cell_config)::cel2b
    type(box_config)::box2b
    integer,intent(out)::imtot
    character,intent(in) :: fnamg*80
    real(double),intent(in)::rum
    type (atom_config)::COMPatrcf
    type(atom_config)::atrgin
    type(box_config)::boxrgin
    real(double)::atg(3,3)
    integer::lat(3),ic,ncore,npr,ierr
    call read_gin(boxrgin,atrgin,fnamg,lat)
    do ic=1,3
       atg(:,ic)=boxrgin%at(:,ic)*lat(ic)
    end do
    call initbox(box2b,atg)
    
    call setnox(box2b,cel2b,rum)
    if (rang==0) then
       write (6, '(2A,D15.8,A,D15.8,A)') fnamg,'volume=', box2b%volu,' cm3 ',box2b%volu*1d24,' Ang3'
    end if
#ifdef PARA
    ncore=0
    call  decoupage(nprocs,ncore,cel2b,at2b)
    call constr_2gin (COMPatrcf,box2b,cel2b,atrgin,boxrgin,lat,imm_glob)
    imtot=COMPatrcf%im
    im_glob=COMPatrcf%im
   call cryst_to_cart (COMPatrcf%imm, COMPatrcf%xp, box2b%at, 1)
    call repartition(COMPatrcf,at2b,box2b,cel2b)
    
#else
    if (ldecoup) then
       open(123, file='decoup.dat', status='old')
       read (123, *) npr,ncore
       close(123)         
       call  decoupage(npr,ncore,cel2b)
       stop
    end if
    call constr_2gin (at2b,box2b,cel2b,atrgin,boxrgin,lat,imm)
    call cryst_to_cart (at2b%imm, at2b%xp, box2b%at, 1)

    im_glob=at2b%im
#endif             


    call setcellconf(cel2b,at2b,box2b,im_glob,rum)
    return

  end subroutine gin2ndm

  subroutine read_cin(boxcin,itread,atcinr,immr,fnamcin,lres,fmtcin,icible,imic)
    !itread 0=at seulement; 1=complet; 2 = at, xp et num_at_glob seulement
    !immr : imm extrait .cin
    !lres : lrestart,
    !fmtcin=1 avec num_at_glob (optional)
    !icible tableau de taille imic qui donne les atomes à lire (utile pour para), optionel
    USE T_kind_param_m, ONLY:  double
    !    USE suivinonpbc
#ifdef PARA
    !    use mpi
    USE var_pot, ONLY:ntyp
    !    USE mod_para,only:MPI_COMM_space,ierr,NDM_MPI_REAL_DOUBLE,status,nprocs,proc_cell
#endif
    USE gen_com_m,only: it,itmax,nitmax,pmean,oldtstep,timel,two,usdh,dilat,tmean,tstep
    implicit none
    character,intent(in) :: fnamcin*80
    integer,intent(in)::itread
    type(box_config)::boxcin
    class(atom_config),optional::atcinr
    integer,intent(in),optional::immr,fmtcin
    logical,intent(in),optional::lres
    integer,allocatable,optional,intent(in)::icible(:)
    integer,optional,intent(in)::imic
    
    logical::lrestart=.false.
    integer :: i, ic, icintype, icintypemod , lucin,fmt_cin=1
    integer, dimension(:),allocatable     :: ibuffer
    real(double), dimension(:,:),allocatable    :: buffer
    real(double)::at(3,3)
#ifdef PARA
    integer,      dimension(ntyp)         :: na_loc

    integer  :: pointeur_loc
    integer  :: numcell, numproc
    integer  :: i_glob
    integer  :: cellx,celly,cellz
#endif
    !      integer , dimension(imm,ntyp) :: fv
    !         integer , dimension(6000,10) :: fv    !Truc_bizarre_jmd
    integer::im_gr,i_loc
    if (present(immr))then
       allocate (ibuffer(immr))
       allocate (buffer(3,immr))
    endif
    if (present(lres))lrestart=lres
        if (present(fmtcin))fmt_cin=fmtcin
    if (rang==0) then
       write(6,*)
       write(6,*)' *-*-*-*-*-*LECTURE DE CIN*-*-*-*-*-*-'
       write(6,*)' *-*-*-*-*- LRESTART =',lrestart
    endif

    lucin = 93
    open(unit=lucin, file=fnamcin, form='unformatted', status='unknown', err=456)

    read (lucin, err=456) icintype

    !          if (rang==0) write (6, *) 'config type de fichier .cin : ', icintype
    if (icintype>3.or.icintype<0) then
       write (6, *) rang, 'wrong icintype'
       call arret_ndm
    endif

    icintypemod = mod(icintype,2)

    if (lrestart.and.icintypemod==0) then
       write (6, *) rang, 'not possible to restart from this file'
       call arret_ndm
    endif
    !at(vect123,xyz)
    !        if (icintype>=2) then
    read (lucin, err=456) at
    if(dilat(1).ne.0.0)then
       do i=1,3
          at(i,:)=at(i,:)*dilat(i)
       end do
    end if

    call initbox(boxcin,at)

    select case(itread)
    case(0)
       return
    case(1)
       if (.not.present(atcinr))then
          write(6,*)'atcinr pas present et itread=1'
          stop
       end if

       read (lucin, err=456) im_gr                         !number of atoms in the box
       if (im_gr>immr) then
          if(rang==0)                    write (6, *) 'im > imM', im_gr, immr
          call arret_ndm
       endif
       atcinr%im=im_gr

       !                    write(6,*)im
       read (lucin, err=456) ibuffer   !ityp
       atcinr%ityp(1:im_gr)=ibuffer(1:im_gr)

       if (rang==0) write (6, *) 'types'
       read (lucin, err=456) buffer    ! xp
       atcinr%xp(1:3,1:im_gr)=buffer(1:3,1:im_gr)
       if (rang==0) write (6, *) 'xp'
       formcin:select case (fmt_cin)
       case (0) formcin
          do i=1,im_gr
             atcinr%num_at_glob(1:im_gr) = i
          enddo
       case(1) formcin
          read (lucin, err=456) ibuffer
          atcinr%num_at_glob(1:im_gr)=ibuffer(1:im_gr)
          if (rang==0) write (6, *) 'num_at_glob'
       case default  formcin
          if (rang.eq.0) write(6,*) 'precisez le format fmt_cin'
          call arret_ndm
       end select formcin


       select type(atcinr)
       type is (atom_config)
          read (lucin, err=456) buffer                     !xpp
          read (lucin, err=456) buffer                     !vp
          lvpread=.false.
       type is (atom_config_d)
          if (icintypemod==1) then
             read (lucin, err=456) buffer                     !xpp
             atcinr%xpp(:,1:im_gr)=buffer(:,1:im_gr)
             write(6,*)'xpp'
             read (lucin, err=456) buffer                     !vp
             atcinr%vp(:,1:im_gr)=buffer(:,1:im_gr)
             write(6,*)'vp'
             !             read (lucin, err=456) buffer                     !former positions
             !             atcinr%ax(:,1:im_gr)=buffer(:,1:im_gr)
             read (lucin, err=456) buffer                     !ax inutile
          end if
       type is (atom_config_e)
          if (icintypemod==1) then
             read (lucin, err=456) buffer                     !xpp
             atcinr%xpp(:,1:im_gr)=buffer(:,1:im_gr)
             write(6,*)'xpp_e'
             read (lucin, err=456) buffer                     !vp
             atcinr%vp(:,1:im_gr)=buffer(:,1:im_gr)
             write(6,*)'vp_e'
             !             read (lucin, err=456) buffer                     !former positions
             !             atcinr%ax(:,1:im_gr)=buffer(:,1:im_gr)
             read (lucin, err=456) buffer                     !ax utile peut-être
             if (lrestart) then 
                if (atcinr%lax) then
                   atcinr%ax(:,1:im_gr)=buffer(:,1:im_gr)
                end if
             end if
          else
             if (atcinr%lax) then
                atcinr%ax(:,1:im_gr)=atcinr%xp(:,1:im_gr)
                lvpread=.false.
             end if
          end if
       end select
       read (lucin, err=456) oldtstep
       !  Si l'option de redemarrage (lrestart) n'est pas activee
       !  alors les positions d'origine ax deviennent les xp du fichier .cin
       !       if (.not.lrestart) then
       !          atcinr%ax(:,:im) = atcinr%xp(:,:im)
       !          if (lsuivinonpbc) axnonpbc(:,:im)=ax(:,:im)
       !       endif


       if (lrestart) then
          read (lucin, err=456) tmean, pmean, it, timel
          if (nitmax.ge.0) itmax=it+nitmax
          tstep = oldtstep

          if (rang==0) then

             write (6, *) 'restart parameters'
             write (6, *) 'it =', it, ' time =', timel
             write (6, *) 'pmean', pmean, ' tmean =', tmean
             write (6, *) 'tstep', tstep
          endif                                ! fin rang=0
          usdh = 1.0/(two*tstep)
       endif
    case(2) ! at xp et num_at_glob
       if (.not.present(atcinr))then
          write(6,*)'atcinr pas present et itread=2'
          stop
       end if

       read (lucin, err=456) im_gr                         !number of atoms in the box
       if (im_gr>immr) then
          if(rang==0)                    write (6, *) 'im > imM', im_gr, immr
          call arret_ndm
       endif
       atcinr%im=im_gr

       read (lucin, err=456) ibuffer   !ityp muet
       read (lucin, err=456) buffer    ! xp
       atcinr%xp(1:3,1:im_gr)=buffer(1:3,1:im_gr)

       read (lucin, err=456) ibuffer
       atcinr%num_at_glob(1:im_gr)=ibuffer(1:im_gr)

    case(3)
        if (.not.present(atcinr))then
          write(6,*)'atcin pas present et itread=3'
          stop
       end if

       read (lucin, err=456) im_gr                         !number of atoms in the box
       if (im_gr>immr) then
          if(rang==0)                    write (6, *) 'im > imM', im_gr, immr
          call arret_ndm
       endif
       atcinr%im=imic

       read (lucin, err=456) ibuffer   !ityp
       do i_loc=1,imic
          atcinr%ityp(i_loc)=ibuffer(icible(i_loc))
       enddo


       if (rang==0) write (6, *) 'types'
       read (lucin, err=456) buffer    ! xp muet

       select case (fmt_cin) !num_at_glob muet
       case (0) 
       case(1) 
          read (lucin, err=456) ibuffer
       end select


       select type(atcinr)
       type is (atom_config)
          read (lucin, err=456) buffer                     !xpp
          read (lucin, err=456) buffer                     !vp
          lvpread=.false.
       type is (atom_config_d)
          if (icintypemod==1) then
             read (lucin, err=456) buffer                     !xpp
             do i_loc=1,imic
                atcinr%xpp(:,i_loc)=buffer(:,icible(i_loc))
             enddo
             read (lucin, err=456) buffer                     !vp
             do i_loc=1,imic
                atcinr%vp(:,i_loc)=buffer(:,icible(i_loc))
             enddo
             read (lucin, err=456) buffer                     !ax inutile
          end if
       type is (atom_config_e)
          if (icintypemod==1) then
             read (lucin, err=456) buffer                     !xpp
             do i_loc=1,imic
                atcinr%xpp(:,i_loc)=buffer(:,icible(i_loc))
             enddo
             read (lucin, err=456) buffer                     !vp
             do i_loc=1,imic
                atcinr%vp(:,i_loc)=buffer(:,icible(i_loc))
             enddo
             read (lucin, err=456) buffer                     !ax utile peut-être
             if (lrestart) then 
                if (atcinr%lax) then
                   do i_loc=1,imic
                      atcinr%ax(:,i_loc)=buffer(:,icible(i_loc))
                   enddo
                end if
             end if
          else
             if (atcinr%lax) then
                atcinr%ax(:,1:imic)=atcinr%xp(:,1:imic)
                lvpread=.false.
             end if
          end if
       end select
       read (lucin, err=456) oldtstep

       if (lrestart) then
          read (lucin, err=456) tmean, pmean, it, timel
          if (nitmax.ge.0) itmax=it+nitmax
          tstep = oldtstep

          if (rang==0) then

             write (6, *) 'restart parameters'
             write (6, *) 'it =', it, ' time =', timel
             write (6, *) 'pmean', pmean, ' tmean =', tmean
             write (6, *) 'tstep', tstep
          endif                                ! fin rang=0
          usdh = 1.0/(two*tstep)
       endif
      
    end select
    close (lucin)
    return
456 print *,'Erreur dans la lecture du fichier .cin, verifier son format&
         & et fmt_cin ATTENTION A BIG_ENDIAN !! SI COMMPILE BIG_ENDIAN NE LIT PLUS QUE CA'


  end subroutine read_cin
  !******************************************************************************************************
   !******************************************************************************************************
  !******************************************************************************************************
  
  subroutine read_gin (boxrg,atrg,fnamgin,latr)
    USE T_kind_param_m, ONLY:  double



    
    character,intent(in) :: fnamgin*80
    type(atom_config),intent(out)::atrg
    type(box_config),intent(out)::boxrg
    integer,intent(out)::latr(3)

!    real(double)::rumax_init,alpha_init
    integer ,     dimension(:),   allocatable :: itypc
    real(double), dimension(:,:), allocatable :: xc
    real(double),dimension(:,:),allocatable :: tmpxc

    real(double)::at(3,3)

    integer ::  lugin, imcell, la, lb, lc, icell,ic,i,ia,ib
    
    !  si coordonnees reduites

    if (rang==0) write (6, *) '**********construction du reseau************'
  lugin=92
    !                                                !number of cells in 3 directions
     open(unit=lugin, file=fnamgin, status='unknown')
    read (lugin, *) latr(1), latr(2), latr(3)
    if (rang==0) write (6, *) 'repetition de mailles', latr


    !     **** coordonnes des vecteurs de maille en A dans une base orthonormee ****
    !                                                !a
    read (lugin, *) at(1,1), at(2,1), at(3,1)
    !b
    read (lugin, *) at(1,2), at(2,2), at(3,2)
    !                                                !c
    read (lugin, *) at(1,3), at(2,3), at(3,3)
    at=at*1d-8
    call initbox(boxrg,at)
    read (lugin, *) imcell               !number of atoms in UC
    if (imcell>imm) then
       if(rang==0)               write (6, *) 'trop d_atomes dans la cel. unite'
       call arret_ndm
    endif
    call atrg%init(imcell)
    do i = 1, imcell
       read (lugin, *) atrg%xp(1,i), atrg%xp(2,i), atrg%xp(3,i),atrg%ityp(i)
    end do
    if (lperiod) then
       do i=1,imcell
          WHERE ( (atrg%xp(:,i).LT.0.d0).OR.(atrg%xp(:,i).GE.1.d0) )
             atrg%xp(:,i)  = atrg%xp(:,i)  - Dble(Floor(atrg%xp(:,i)))
          END WHERE
       end do
    end if

    close(lugin)
  end subroutine read_gin
  
  !#endif
end module constrconf_mod
