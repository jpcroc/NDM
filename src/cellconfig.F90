module cellconfig
  USE arret_ndm_mod,only:arret_ndm
  USE T_kind_param_m
  use atomconfig,only : atom_config,atom_config_d,atom_config_e
  use boxconfig,only:box_config
  use paraconfig,only:para_config
  use Tpara,only:para_space_config,mpi_communicator,myidsp,comm_space
  use gen_com_m,only:uwrt,lwrt,rang
  implicit none
  !  integer:: incr=20 ! incrément des tailles de tableau 

  type cell_config
     integer:: nox(3)=0,noxyz=0 ,noxyzact
     integer::natperc !nombre (max) d'atomes par cellule
     integer(long)::icaltabt 
     integer,allocatable::nato (:) ! nombre d'atomes dans la cellule ko
     integer,allocatable:: ncel (:,:) ! ncel(ko,i1)=numéro de la ième cellule voisine de la cellule ko
     integer,allocatable:: atincel (:,:) ! atincel(i,k)= indice du ième atome de la cellule k
     integer,allocatable:: deltadist(:,:,:) !gestion des conditions périodiques entre les cellules (voisinage de bords de boites)
     integer,allocatable:: ncelvois(:) !nombre de cellules voisines de la cellule actuelle (26 pour PBC standard; dépend de la position dans la boite pour les PBC partielles
     logical :: ltpcel
     real(double),allocatable::sigc(:,:,:),tempc(:)
     real(double):: celsize(3)
#ifdef PARA
     integer,allocatable::proc_cell(:)
     logical:: ismall(3)=.false.
     logical,allocatable::isghost(:)
     integer,allocatable::copyof(:)
     integer::ncelvmax ! maximum number of cells neighbouring a cell, for large cells=26, for small nox*noy*noz-1

#endif     

   contains
     procedure, pass::init=>init_cel
     procedure, pass::dealloc=>dealloc_cel
     procedure, pass::copy
     procedure, pass::print=>cellprint
     procedure, pass::send2proc=>cells2p
     procedure, pass::send2all=>cells2a
     procedure, pass::recv=>cellrecv
     procedure, pass::koxyz=>kox
     procedure, pass::edge=>edgec
     procedure, pass::constrcomp=>cococe
  end type cell_config

  type, extends (cell_config):: cell_config_g !
     integer,allocatable::natotot (:) ! nombre d'atomes dans la cellule ko
  end type cell_config_g

  type, extends (cell_config):: cell_config_arps !
     integer,allocatable::nmov(:,:) ! nombre d'atomes dans la cellule ko
  end type cell_config_arps


contains

  function edgec(cell,ko,box)
    class(cell_config)::cell
    class(box_config)::box
    real(double):: edgec(3)
    integer::ko
    integer::kox(3)

    integer:: ic,ic2
    real(double)::flx(3)
    kox=cell%koxyz(ko)
    do ic=1,3
       if (cell%ismall(ic)) then
          flx(ic)=kox(ic)-0.5*(cell%nox(ic)-1)
       else
          flx(ic)=float(kox(ic))/cell%nox(ic)
       end if
    end do
    edgec=0.
    do ic=1,3
       do ic2=1,3
          edgec(ic2)=edgec(ic2)+flx(ic)*box%at(ic2,ic)
       end do
    end do
  end function edgec

  function kox(cell,ko)
    class(cell_config)::cell
    integer:: kox(3),ko

    integer:: kx,ky,kz,kyz,km1,km2,nox,noy
    !    do ko=1,cell%noxyz
    nox=cell%nox(1);noy=cell%nox(2)
    km1=ko-1
    kx=mod(km1,nox)
    km2=(km1-kx)/nox
    ky=mod(km2,noy)
    kz=(km2-ky)/noy
    kox(1)=kx;kox(2)=ky;kox(3)=kz

    !   end do

  end function kox


  subroutine init_cel(cell,box,nox,noy,noz,natperc,ltpc,latomalloc)
    class(cell_config)::cell
    integer,intent(in),optional::nox,noy,noz,natperc
    logical,optional,intent(in):: ltpc
    logical::ltpcel
    logical,optional::latomalloc
    logical::lata=.true.
    class(box_config)::box
    ltpcel=.false.
    if (present(latomalloc))lata=latomalloc
    if (present (ltpc))ltpcel=ltpc
    if (present(nox)) then
       cell%nox(1)=nox;cell%noxyz=nox*noy*noz
       cell%nox(2)=noy
       cell%nox(3)=noz
    end if
    if (present(natperc)) then
       cell%natperc=natperc
    else
       cell%natperc=0
    end if
    cell%icaltabt=0
    !write(uwrt,*) 'nox', cell%nox(1)
    cell%ltpcel=ltpcel
    call dealloc_cel(cell)
    call allocatecelN(cell,lata)
    call neigcelN(cell,box)

    return

  end subroutine init_cel

  subroutine allocatecelN(cell,latomalloc)
    class(cell_config)::cell
    !    integer,intent(in)::nox,noy,noz,natperc
    integer::nsize,ic
    logical,optional::latomalloc
    logical::lata=.true.
    if (present(latomalloc))lata=latomalloc
    cell%noxyz=cell%nox(1)*cell%nox(2)*cell%nox(3)
    !    write(uwrt,*)'NOX',cell%nox,cell%noxyz
    cell%ncelvmax=1
    do ic=1,3
       if(cell%ismall(ic)) then
          cell%ncelvmax=cell%ncelvmax*cell%nox(ic)
       else
          cell%ncelvmax=cell%ncelvmax*3
       end if
    end do
    cell%ncelvmax=cell%ncelvmax-1

    nsize=cell%noxyz
    if (nsize.ne.0) then
       allocate(cell%ncel(nsize,0:cell%ncelvmax))
       cell%ncel=0
       allocate(cell%nato(nsize))
       cell%nato=0
       select type(cell)
       class is (cell_config_g)
          allocate(cell%natotot(nsize))
          cell%natotot=0
       end select
       allocate(cell%ncelvois(nsize))
       allocate(cell%deltadist(3,0:cell%ncelvmax,nsize))
       cell%deltadist=0
       allocate (cell%isghost(nsize))
       cell%isghost=.false.
       if (any(cell%ismall(:)))then
          allocate (cell%copyof(nsize))
          cell%copyof(:)=-1

       end if


       if (cell%ltpcel) then
          allocate(cell%sigc(3,3,nsize))
          cell%sigc=0
          allocate(cell%tempc(nsize))
          cell%tempc=0
       end if
       if ((cell%natperc.ne.0).and.(lata))   then
          allocate(cell%atincel(cell%natperc,nsize))
          cell%atincel=0
       end if
#ifdef PARA
       allocate(cell%proc_cell(nsize))
       cell%proc_cell=0
#endif       
    end if
    return
  end subroutine allocatecelN



  subroutine dealloc_cel(cell)
    class(cell_config)::cell

    if (allocated(cell%ncel))       deallocate(cell%ncel)
    if (allocated(cell%ncelvois))       deallocate(cell%ncelvois)
    if (allocated(cell%nato))       deallocate(cell%nato)
    if (allocated(cell%atincel))    deallocate(cell%atincel)
    if (allocated(cell%deltadist))  deallocate(cell%deltadist)
    if (allocated(cell%sigc))  deallocate(cell%sigc)
    if (allocated(cell%tempc))  deallocate(cell%tempc)
    if (allocated(cell%isghost))  deallocate(cell%isghost)
    if (allocated(cell%copyof))  deallocate(cell%copyof)
    select type(cell)
    class is (cell_config_g)
       deallocate(cell%natotot)
       !       cell%natotot=0
    end select

#ifdef PARA
    if (allocated(cell%proc_cell))  deallocate(cell%proc_cell)
#endif

    return

  end subroutine dealloc_cel

  subroutine neigcelN(cell,box)
    class(cell_config)::cell
    class(box_config),intent(in)::box
    integer :: kx, ky, kz, koo, l, lz, mz, ly, my, lx, mx, kxy!,ldx,lfx,ldy,lfy,ldz,lfz
    integer::midnox(3),mindecx(3),maxdecx(3)
    integer:: vx,vy,vz,vxyz,ic,koxyz(3),kcp,koxyzv(3)
    cell%ncelvois=-1
    do ic=1,3
       if (cell%ismall(ic)) then
          midnox(ic)=(cell%nox(ic)+1)/2
          mindecx(ic)=-(cell%nox(ic)-1)/2
          maxdecx(ic)=(cell%nox(ic)-1)/2
       else
          mindecx(ic)=-1
          maxdecx(ic)=1          
       end if
    end do
    !    write(uwrt,*)'mindecx', mindecx
    !    write(uwrt,*)'maxdecx', maxdecx
    !    write(uwrt,*)'midnox', midnox
    if (cell%noxyz==1) then
       cell%ncel(1,0)=1
       cell%deltadist=0
       cell%ncelvois(1)=0
    else
       cell%deltadist(:,:,:) = 0 ! par défaut celldeltadist=0
       do kz = 1, cell%nox(3)
          do ky = 1, cell%nox(2)
             loopext:do kx = 1, cell%nox(1)
                koo = 1+(kx-1)+cell%nox(1)*((ky-1)+cell%nox(2)*(kz-1))
                if (any(cell%ismall)) then
                   if((cell%ismall(1)).and.(kx.ne.midnox(1)))then
                      cell%isghost(koo)=.true.
                   end if
                   if((cell%ismall(2)).and.(ky.ne.midnox(2)))then
                      cell%isghost(koo)=.true.
                   end if
                   if((cell%ismall(3)).and.(kz.ne.midnox(3)))then
                      cell%isghost(koo)=.true.
                   end if
                   if (cell%isghost(koo)) then
                      if (cell%ismall(1)) then
                         vx=midnox(1)
                      else
                         vx=kx
                      end if
                      if (cell%ismall(2)) then
                         vy=midnox(2)
                      else
                         vy=ky
                      end if
                      if (cell%ismall(3)) then
                         vz=midnox(3)
                      else
                         vz=kz
                      end if
                      vxyz = 1+(vx-1)+cell%nox(1)*((vy-1)+cell%nox(2)*(vz-1))
                      cell%copyof(koo)=vxyz
                      cycle loopext 
                   end if

                end if

                cell%ncel(koo,0) = koo
                !                cell%deltadist(:,0,koo) = 0
                l = 0 

                do lz = mindecx(3),maxdecx(3)
                   do ly = mindecx(2),maxdecx(2)
                      loopin:  do lx = mindecx(1),maxdecx(1)
                         if((lx==0).and.(ly==0).and.(lz==0)) cycle loopin !cellule en cours d'analyse (lignes au dessus)
                         mz = kz+lz
                         mx = kx+lx
                         my = ky+ly
                         if(box%ipbc(1).ne.1) then
                            if ((cell%ismall(1)).and.(lx.ne.0)) then
                               cycle loopin
                            else
                               if((mx<1).or.(mx>cell%nox(1))) then !débordement
                                  cycle loopin
                               end if
                            end if
                         else
                            if (.Not.cell%ismall(1)) then
                               select case(cell%nox(1))
                               case(1)
                                  if ((lx==-1).or.(lx==1)) cycle loopin
                               case(2)
                                  if (lx==-1) cycle loopin
                               end select
                            end if
                         end if
                         if(box%ipbc(2).ne.1) then
                            if ((cell%ismall(2)).and.(ly.ne.0)) then
                               cycle loopin
                            else
                               if((my<1).or.(my>cell%nox(2))) then !débordement
                                  cycle loopin
                               end if
                            end if
                         else
                            if (.Not.cell%ismall(2)) then
                               select case(cell%nox(2))
                               case(1)
                                  if ((ly==-1).or.(ly==1)) cycle loopin
                               case(2)
                                  if (ly==-1) cycle loopin
                               end select
                            end if
                         end if

                         if(box%ipbc(3).ne.1) then
                            if ((cell%ismall(3)).and.(lz.ne.0)) then
                               cycle loopin
                            else
                               if((mz<1).or.(mz>cell%nox(3))) then !débordement
                                  cycle loopin
                               end if
                            end if
                         else
                            if (.Not.cell%ismall(3)) then
                               select case(cell%nox(3))
                               case(1)
                                  if ((lz==-1).or.(lz==1)) cycle loopin
                               case(2)
                                  if (lz==-1) cycle loopin
                               end select
                            end if
                         end if

                         l=l+1 ! on est dans une cellule voisine
!!$                         write(uwrt,*)'kx ky kz koo', kx,ky,kz,koo
!!$                         write(uwrt,*)'mx my mz ', mx,my,mz
                         if (cell%ismall(3)) then
                            cell%deltadist(3,l,koo) = -1*(mz-midnox(3))
                         else
                            if (mz<1) then ! on ne peut pas être ici si ipbc(3).ne.1
                               mz = mz+cell%nox(3)
                               cell%deltadist(3,l,koo) = 1
                            endif
                            if (mz>cell%nox(3)) then
                               mz = mz-cell%nox(3)
                               cell%deltadist(3,l,koo) = -1
                            endif
                         end if
                         if (cell%ismall(2)) then
                            cell%deltadist(2,l,koo) = -1*(my-midnox(2))
                         else
                            if (my<1) then
                               my = my+cell%nox(2)
                               cell%deltadist(2,l,koo) = 1
                            endif
                            if (my>cell%nox(2)) then
                               my = my-cell%nox(2)
                               cell%deltadist(2,l,koo) = -1
                            endif
                         end if

                         if (cell%ismall(1)) then
                            cell%deltadist(1,l,koo) = -1*(mx-midnox(1))
                         else
                            if (mx<1) then
                               mx = mx+cell%nox(1)
                               cell%deltadist(1,l,koo) = 1
                            endif
                            if (mx>cell%nox(1)) then
                               mx = mx-cell%nox(1)
                               cell%deltadist(1,l,koo) = -1
                            endif
                         end if
                         kxy = 1+(mx-1)+cell%nox(1)*((my-1)+cell%nox(2)*(mz-1))
                         !                         if (kxy==koo) cycle
                         cell%ncel(koo,l) = kxy
                         !                        write(uwrt,*)koo,lz,ly,lx,l,kxy
                         !                        if ((kz==cell%nox(3)).and.(lz==1))write(uwrt,*)koo,lz,l,kxy
                         !                        if ((kz==1).and.(lz==-1))write(uwrt,*)koo,lz,l,kxy
                      end do loopin
                   end do
                end do
                cell%ncelvois(koo)=l
             end do loopext
          end do
       end do
    end if
    if (allocated(cell%isghost)) then
       cell%noxyzact = COUNT(.NOT. cell%isghost)
    else
       cell%noxyzact = cell%noxyz
    end if
    if (rang==0)write(uwrt,*)cell%noxyzact ,' active cells among ', cell%noxyz





    return
  end subroutine neigcelN


  subroutine caltabtC (cell,atcf,lperiod,boxcf,lextr,psc,lchktrav)
    USE notperiod_mod,only: notperiod
    USE cryst_to_cart_mod,only: cryst_to_cart
    use gen_com_m,only:uwrt,lwrt,lspacendm
#ifdef PARA
    USE Tpara,only:comm_space
#endif

    class(cell_config), intent(inout):: cell
    class(atom_config),intent(inout)::atcf
    class(box_config),intent(inout)::boxcf
    type(para_space_config),optional::psc    
    logical,intent(in)::lperiod
    logical,intent(in),optional::lextr
    logical::lchktrav
    logical::lextrait=.false.

    integer :: i,  koo,kxyz(3)
    real(double) :: auxyz(3)
    real(double), dimension(:,:), allocatable :: xpnp !
    integer(long), save:: icaltabt=0
    integer::iml,midnox(3)
    integer::indor,ic



    integer :: ntravtot,itrav
    integer,allocatable:: ntrav(:)
    integer,parameter::maxtrav=10
    integer,allocatable:: indtrav(:),proccib(:)
    allocate (ntrav(0:comm_space%nproc-1)) ; ntrav(:)=0
    !
    ! --------- Initialisation --------------
    !
    midnox(:)=(cell%nox(:)+1)/2
    if (present(lextr))lextrait=lextr
    if (lextrait) then
       iml=atcf%imm
    else
       iml=atcf%im
    end if
    icaltabt=icaltabt+1
    !       write(uwrt,*)'caltabt',icaltabt

    cell%nato(1:cell%noxyz) = 0
    cell%atincel(1:cell%natperc,1:cell%noxyz) = 0

    !  -------- cas sans cellule  -----------
    if (atcf%im==0) then
       cell%nato(:)=0
       goto 101
    end if
    if (cell%noxyz==1) then
       cell%nato(1) = iml
       do i = 1, iml
          atcf%ielat(i) = 1
          cell%atincel(i,1) = i
       end do
    else

       ALLOCATE(xpnp(3,iml))
       call notperiod(iml,atcf%xp,xpnp,boxcf%at,boxcf%bg,lperiod)       

       !  -------- Initialisations  -----------

       ! -------------------------------------------
       !   1. loop: lattice-coordinates of all atoms
       ! - - - - - - - - - - - - - - - - - - - - - -

       !debug       write (*,*) 'sub caltabt 1',it,xp(1,1)
       call cryst_to_cart (iml, xpnp, boxcf%bg, -1) ! cart vers cryst
       if (any(xpnp(:,1:iml).gt.1).or.any(xpnp(:,1:iml).lt.0)) then
          write(uwrt,*)'PLANTE',rang
          !          write(300+RANG,*)'PLANTE'
!          do i=1,iml
!             if (any(xpnp(:,i).gt.1).or.any(xpnp(:,i).lt.0))  write(uwrt,*) i,xpnp(:,i)
!          end do
          write(uwrt,*)'caltabtc xpnp <0 ou >1 stop'
          call arret_ndm(.true.)
       end if
       !debug       write (*,*) 'sub caltabt 2',it,xp(1,1)

       !     if (it.gt.1000) write(uwrt,*)'CALTABT',it
       !       write(uwrt,*)'caltabt icaltabt im',icaltabt,atcf%im

       do i = 1, iml
          !     if  ((it.ge.1000).and.(i.lt.20)) write(uwrt,'(I5,3G15.7)')i, xpnp(1,i),xpnp(2,i),xpnp(3,i)
          do ic=1,3
             if  (cell%ismall(ic)) then
                kxyz(ic)=midnox(ic)-1
             else
                auxyz(ic) = xpnp(ic,i)*cell%nox(ic)
                kxyz(ic) = int(auxyz(ic))
                kxyz(ic) = Modulo(kxyz(ic),cell%nox(ic))
             end if
          end do
          koo = 1+kxyz(1)+cell%nox(1)*(kxyz(2)+cell%nox(2)*kxyz(3))

          IF ( (koo.GT.cell%noxyz).OR.(koo.LT.0) ) THEN
             WRITE(0,'(a,i0,a,3g20.12)') &
                  'Problem with atom ', i, ', x,y,z = ', atcf%xp(1:3,i)
             WRITE(0,'(2(a,i0))') ' koo = ', koo, ' - noxyz = ', cell%noxyz
             STOP
          END IF
          if ((present(psc)).and.(lspacendm).and.lchktrav) then 
             if (cell%proc_cell(koo).ne.myidsp) then
                !                write(uwrt,*)'atout',i,atcf%num_at_glob(i),rang,cell%proc_cell(koo),koo

                if(.not.(any(psc%cell_ftm(:)==koo))) then
                   write(uwrt,'(A,6I7)')'WARNING ::: attrrav:i natg ielat rangem rangf newcell',i,&
                        &atcf%num_at_glob(i),atcf%ielat(i),rang,cell%proc_cell(koo),koo
                   write(uwrt,'(A,6G20.7)')'WARNING ::: travelled from cell to cell ',cell%edge(atcf%ielat(i),boxcf)&
                        &,cell%edge(koo,boxcf)
                   ntrav(myidsp)=ntrav(myidsp)+1
                   if (ntrav(myidsp)==1) then
                      allocate(indtrav(maxtrav));allocate(proccib(maxtrav)); proccib=-1
                   end if
                   proccib(ntrav(myidsp))=cell%proc_cell(koo)
                   indtrav(ntrav(myidsp))=i
                   atcf%ielat(i) = koo
                   !                   cell%nato(koo) = cell%nato(koo)+1
                   !                   cell%atincel(cell%nato(koo),koo) = i                   
                   IF (cell%nato(koo).GT.cell%natperc) THEN
                      WRITE(0,'(a)') 'caltabtC : You need to increase the maximal number of atoms per cell'
                      WRITE(0,'(a,i0)') 'current value: natperc=', cell%natperc
                      STOP '< CaltabtC >'
                   END IF


                   !                   write(uwrt,*)'atom', i,atcf%num_at_glob(i),'in cell', koo, ' originally in proc', &
                   !                        &myidsp, 'now in ', cell%proc_cell(koo),' travelled too far. its cell is not a frontier cell',&
                   !                        &'RANG actuel = ',rang
                   !                   call arret_ndm
                else
                   !                   if (cell%proc_cell(koo).ne.myidsp)       write(uwrt,'(A,6I7)')'afrt:i natg ielat rangem rangf newcell',i,atcf%num_at_glob(i),atcf%ielat(i),rang,cell%proc_cell(koo),koo
                   atcf%ielat(i) = koo
                   cell%nato(koo) = cell%nato(koo)+1
                   cell%atincel(cell%nato(koo),koo) = i
                end if

             else

                atcf%ielat(i) = koo
                cell%nato(koo) = cell%nato(koo)+1
                IF (cell%nato(koo).GT.cell%natperc) THEN
                   WRITE(0,'(a)') 'caltabtC : You need to increase the maximal number of atoms per cell'
                   WRITE(0,'(a,i0)') 'current value: natperc=', cell%natperc
                   STOP '< CaltabtC >'
                END IF
                cell%atincel(cell%nato(koo),koo) = i
             end if
          else
             atcf%ielat(i) = koo
             cell%nato(koo) = cell%nato(koo)+1
             IF (cell%nato(koo).GT.cell%natperc) THEN
                WRITE(0,'(a)') 'caltabtC : You need to increase the maximal number of atoms per cell'
                WRITE(0,'(a,i0)') 'current value: natperc=', cell%natperc
                STOP '< CaltabtC >'
             END IF
             cell%atincel(cell%nato(koo),koo) = i
          end if



       end do

       ntravtot=ntrav(myidsp)
#ifdef PARA
       if (lchktrav) then 
          call comm_space%sum(ntravtot)
          if (ntravtot.gt.0) then
             if (.not.allocated(proccib)) then
                allocate(proccib(maxtrav))
                proccib=-1
             end if
             call transfer_atoms(atcf,cell,indtrav,ntravtot,ntrav,proccib)
          end if
       end if
#endif


       DEALLOCATE(xpnp)   ! MODIF CLOUET
    endif
    if (any(cell%ismall)) then
       do koo=1,cell%noxyz
          if (cell%isghost(koo)) then
             indor=cell%copyof(koo)
             cell%atincel(:,koo)=cell%atincel(:,indor)
             cell%nato(koo)=cell%nato(indor)
          end if
       end do
    end if
    !    do koo=1,cell%noxyz
    !       write(uwrt,*)'nato',koo,cell%nato(koo)
    !    end do

101 continue
    cell%icaltabt=icaltabt
    atcf%icaltabt=icaltabt
    boxcf%icaltabt=icaltabt

    select type(cell)
    class is (cell_config_g)
       cell%natotot=cell%nato
#ifdef PARA
       if (lspacendm) call comm_space%sum(cell%natotot)
#endif
    end select

    return
  end subroutine caltabtC


  ! copie d'une config entière vers config de base

  subroutine copy (cellsource,cellcible,box,lzeroinit,ltpccible,latomcp)
    class(cell_config)::cellsource
    class(cell_config)::cellcible
    class(box_config)::box
    logical,optional::lzeroinit,ltpccible,latomcp
    logical::lzi=.false.,ltpcel,latcp=.true.
    if (present(latomcp))latcp=latomcp
    if (present(lzeroinit))lzi=lzeroinit

    call cellcible%dealloc
    !    cellcible=cellsource
    if (present(ltpccible))then
       ltpcel=ltpccible
    else
       ltpcel=cellsource%ltpcel
    end if
    call cellcible%init(box,cellsource%nox(1),cellsource%nox(2),cellsource%nox(3),cellsource%natperc,ltpc=ltpcel,latomalloc=latcp)

    cellcible%nox(1)=cellsource%nox(1)
    cellcible%nox(2)=cellsource%nox(2)
    cellcible%nox(3)=cellsource%nox(3)
    cellcible%noxyz=cellsource%noxyz
    cellcible%natperc=cellsource%natperc
    cellcible%icaltabt=cellsource%icaltabt
    if (lzi) then
       cellcible%nato(:)=0
       if (latcp) cellcible%atincel(:,:)=0
    else
       cellcible%nato(:)=cellsource%nato(:)
       if (latcp)        cellcible%atincel(:,:)=cellsource%atincel(:,:)
       if ((cellcible%ltpcel).and.(cellsource%ltpcel))then
          cellcible%sigc=cellsource%sigc
          cellcible%tempc=cellsource%tempc
       end if
    end if
    cellcible%ncel(:,:)=cellsource%ncel(:,:)
    cellcible%deltadist(:,:,:)=cellsource%deltadist(:,:,:)
    cellcible%celsize=cellsource%celsize
  end subroutine copy

  subroutine cococe(celloc,celcomp,box,latomcp,psc)
    class(cell_config)::celloc
    class(cell_config)::celcomp
    class(box_config)::box
    logical::latomcp
    type(para_space_config),optional::psc    

    integer::iko
    call celcomp%init(box,celloc%nox(1),celloc%nox(2),celloc%nox(3),celloc%natperc,ltpc=celloc%ltpcel,latomalloc=latomcp)

#ifdef PARA
    celcomp%celsize=celloc%celsize
    celcomp%nato=0
    if(latomcp)celcomp%atincel=0
    if (allocated(celcomp%sigc))celcomp%sigc=0
    if (allocated(celcomp%tempc))celcomp%tempc=0
    do iko=1,celloc%noxyz
       if (celloc%proc_cell(iko)==myidsp) then
          celcomp%nato(iko)=celloc%nato(iko)
          if (allocated(celcomp%tempc))celcomp%tempc(iko)=celloc%tempc(iko)
          if (allocated(celcomp%sigc))celcomp%sigc(:,:,iko)=celloc%sigc(:,:,iko)
          if(latomcp)celcomp%atincel(:,iko)=celloc%atincel(:,iko)
       end if
    end do
    call comm_space%sum(celcomp%nato)
    if (allocated(celcomp%tempc))call comm_space%sum(celcomp%tempc)
    if (allocated(celcomp%sigc))call comm_space%sum(celcomp%sigc)
    if(latomcp)call comm_space%sum(celcomp%atincel)

#else
    call celloc%copy(celcomp,box,latomcp=latomcp)
#endif

  end subroutine cococe



  subroutine cellprint(cellv,unit,mess)
    class(cell_config)::cellv
    integer,intent(in),optional::unit
    integer::i,un,j,ko(3)
    character(len=*),optional::mess
    un=6
    if (present(unit))un=unit
    write(un,*)'in cellprint ',mess
    write(un,*)'nox noy noz noxyz',cellv%nox(1),cellv%nox(2),cellv%nox(3),cellv%noxyz
    write(un,*)'natperc',cellv%natperc
    write(un,*)'celsize',cellv%celsize
    write(un,*)'icaltabt',cellv%icaltabt
    write(un,*)'LTPCEL',cellv%ltpcel
    do i=1,cellv%noxyz
       write(un,*)'ncelvois',i,cellv%ncelvois(i)
    end do
    do i=1,cellv%noxyz
       write(un,*)'nato',i,cellv%nato(i)
    end do
    if (allocated(cellv%atincel))then 
       do i=1,cellv%noxyz
          write(un,*)'atincel',i,cellv%atincel(:,i)
       end do
    end if
    if (cellv%ltpcel) then
       if (allocated(cellv%tempc) )then
          do i=1,cellv%noxyz
             write(un,*)'tempc',i,cellv%tempc(i)
          end do
       end if
       if (allocated(cellv%sigc)) then
          do i=1,cellv%noxyz
             write(un,*)'sigc',i,cellv%sigc(1,1,i),cellv%sigc(1,2,i)
          end do
       end if
    end if
    write(un,*)
    write(un,*)'DELTADIST et ghost'
    do i=1,cellv%noxyz
       if (cellv%isghost(i)) then
          write(un,*)'cell i',i,'is ghost, copy of ', cellv%copyof(i)
       else
          ko(:)=cellv%koxyz(i)
          write(un,'(A,5I4)')'CELL i ko(3) ncellvois',i,ko(1),ko(2),ko(3),cellv%ncelvois(i)
          do j=0,cellv%ncelvois(i)
             write(un,'(A,3I4,A,3I4)')'VOISCELL',i,j,cellv%ncel(i,j),'deltadist',cellv%deltadist(:,j,i)
          end do
       end if
    end do
#ifdef PARA
    if (allocated(cellv%proc_cell))then 
       do i=1,cellv%noxyz
          write(un,*)'proc_cell',i,cellv%proc_cell(i)
       end do
    end if
#endif
    !    write(un,*)
    !    do i=1,cellv%noxyz
    !       if (cellv%isghost(i)) then
    !          write(un,*)'cell i',i,'is ghost, copy of ', cellv%copyof(i)
    !       end if
    !    end do


  end subroutine cellprint



  subroutine cells2p(cell,rgcib,mpic)
    class(cell_config)::cell
    type(mpi_communicator),intent(in)::mpic
    integer,intent(in)::rgcib
    integer::sizeI,sizeR,ibi,ibr,nsize,ip,ip2,ip3
    integer,allocatable:: ibuffer(:)
    real(double),allocatable::rbuffer(:)

    sizeI=6+size(cell%nato)+size(cell%ncel)+size(cell%atincel)+size(cell%deltadist)
#ifdef PARA
    sizeI=sizeI+size(cell%proc_cell)
#endif
    nsize=cell%noxyz
    sizeR=3
    if (cell%ltpcel) sizer=sizer+size(cell%sigc)+size(cell%tempc)

    allocate (ibuffer(sizeI)) ; allocate (rbuffer(sizeR))
    ibuffer(1)=cell%nox(1); ibuffer(2)=cell%nox(2) ; ibuffer(3)=cell%nox(3)
    ibuffer(4)=cell%noxyz
    ibuffer(5)=cell%natperc
    ibuffer(6)=int(cell%icaltabt)
    ibi=6
    do ip=1,nsize
       ibi=ibi+1
       ibuffer(ibi)=cell%nato(ip)
    end do
    do ip=1,nsize
       do ip2=0,cell%ncelvmax
          ibi=ibi+1
          ibuffer(ibi)=cell%ncel(ip,ip2)
       end do
    end do
    do ip=1,nsize
       do ip2=0,cell%ncelvmax
          do ip3=1,3
             ibi=ibi+1
             ibuffer(ibi)=cell%deltadist(ip3,ip2,ip)
          end do
       end do
    end do
    do ip=1,nsize
       do ip2=1,cell%natperc
          ibi=ibi+1
          ibuffer(ibi)=cell%atincel(ip2,ip)
       end do
    end do
#ifdef PARA
    do ip=1,nsize
       ibi=ibi+1
       ibuffer(ibi)=cell%proc_cell(ip)
    end do

#endif
    rbuffer(1)=cell%celsize(1);     rbuffer(2)=cell%celsize(2) ;    rbuffer(3)=cell%celsize(3)
    ibr=3
    if (cell%ltpcel) then
       do ip=1,nsize
          do ip2=1,3
             do ip3=1,3
                ibR=ibr+1
                rbuffer(ibr)=cell%sigc(ip2,ip3,ip)
             end do
          end do
       end do
       do ip=1,nsize
          ibr=ibr+1
          rbuffer(ibr)=cell%tempc(ip)
       end do
    end if
    call mpic%send(ibuffer,rgcib,201)
    call mpic%send(cell%ltpcel,rgcib,202)
    call mpic%send(rbuffer,rgcib,203)
  end subroutine cells2p


  subroutine cellrecv(cell,rgem,mpic)
    class(cell_config)::cell
    type(mpi_communicator),intent(in)::mpic
    integer,intent(in)::rgem
    integer::sizeI,sizeR,ibi,ibr,nsize,ip,ip2,ip3
    integer,allocatable:: ibuffer(:)
    real(double),allocatable::rbuffer(:)


    sizeI=6+size(cell%nato)+size(cell%ncel)+size(cell%atincel)+size(cell%deltadist)
#ifdef PARA
    sizeI=sizeI+size(cell%proc_cell)
#endif
    nsize=cell%noxyz
    sizeR=3
    if (cell%ltpcel) sizer=sizer+size(cell%sigc)+size(cell%tempc)


    allocate (ibuffer(sizeI)) ; allocate (rbuffer(sizeR))

    call mpic%recv(ibuffer,rgem,201)
    call mpic%recv(cell%ltpcel,rgem,202)
    call mpic%recv(rbuffer,rgem,203)

    cell%nox(1)=ibuffer(1); cell%nox(2)=ibuffer(2) ; cell%nox(3)=ibuffer(3)
    cell%noxyz=ibuffer(4)
    cell%natperc=ibuffer(5)
    cell%icaltabt=ibuffer(6)
    ibi=6
    do ip=1,nsize
       ibi=ibi+1
       cell%nato(ip)=ibuffer(ibi)
    end do
    do ip=1,nsize
       do ip2=0,cell%ncelvmax
          ibi=ibi+1
          cell%ncel(ip,ip2)=ibuffer(ibi)
       end do
    end do
    do ip=1,nsize
       do ip2=0,cell%ncelvmax
          do ip3=1,3
             ibi=ibi+1
             cell%deltadist(ip3,ip2,ip)=ibuffer(ibi)
          end do
       end do
    end do
    do ip=1,nsize
       do ip2=1,cell%natperc
          ibi=ibi+1
          cell%atincel(ip2,ip)=ibuffer(ibi)
       end do
    end do
#ifdef PARA
    do ip=1,nsize
       ibi=ibi+1
       cell%proc_cell(ip)=ibuffer(ibi)
    end do

#endif
    cell%celsize(1)=rbuffer(1);     cell%celsize(2)=rbuffer(2) ;    cell%celsize(3)=rbuffer(3)
    ibr=3
    if (cell%ltpcel) then
       do ip=1,nsize
          do ip2=1,3
             do ip3=1,3
                ibR=ibr+1
                cell%sigc(ip2,ip3,ip)=rbuffer(ibr)
             end do
          end do
       end do
       do ip=1,nsize
          ibr=ibr+1
          cell%tempc(ip)=rbuffer(ibr)
       end do
    end if
  end subroutine cellrecv

  subroutine cells2a(cell,rgem,mpic)
    type(mpi_communicator),intent(in)::mpic
    class(cell_config)::cell
    integer,intent(in)::rgem
    integer::sizeI,sizeR,ibi,ibr,nsize,ip,ip2,ip3
    integer,allocatable:: ibuffer(:)
    real(double),allocatable::rbuffer(:)

    sizeI=6+size(cell%nato)+size(cell%ncel)+size(cell%atincel)+size(cell%deltadist)
#ifdef PARA
    sizeI=sizeI+size(cell%proc_cell)
#endif
    nsize=cell%noxyz
    sizeR=3
    if (cell%ltpcel) sizer=sizer+size(cell%sigc)+size(cell%tempc)
    allocate (ibuffer(sizeI)) ; allocate (rbuffer(sizeR))
    ibuffer(1)=cell%nox(1); ibuffer(2)=cell%nox(2) ; ibuffer(3)=cell%nox(3)
    ibuffer(4)=cell%noxyz
    ibuffer(5)=cell%natperc
    ibuffer(6)=int(cell%icaltabt)
    ibi=6
    do ip=1,nsize
       ibi=ibi+1
       ibuffer(ibi)=cell%nato(ip)
    end do
    do ip=1,nsize
       do ip2=0,cell%ncelvmax
          ibi=ibi+1
          ibuffer(ibi)=cell%ncel(ip,ip2)
       end do
    end do
    do ip=1,nsize
       do ip2=0,cell%ncelvmax
          do ip3=1,3
             ibi=ibi+1
             ibuffer(ibi)=cell%deltadist(ip3,ip2,ip)
          end do
       end do
    end do
    do ip=1,nsize
       do ip2=1,cell%natperc
          ibi=ibi+1
          ibuffer(ibi)=cell%atincel(ip2,ip)
       end do
    end do
#ifdef PARA
    do ip=1,nsize
       ibi=ibi+1
       ibuffer(ibi)=cell%proc_cell(ip)
    end do

#endif
    rbuffer(1)=cell%celsize(1);     rbuffer(2)=cell%celsize(2) ;    rbuffer(3)=cell%celsize(3)
    ibr=3
    if (cell%ltpcel) then
       do ip=1,nsize
          do ip2=1,3
             do ip3=1,3
                ibR=ibr+1
                rbuffer(ibr)=cell%sigc(ip2,ip3,ip)
             end do
          end do
       end do
       do ip=1,nsize
          ibr=ibr+1
          rbuffer(ibr)=cell%tempc(ip)
       end do
    end if
    call mpic%bcast(rgem,ibuffer)
    call mpic%bcast(rgem,cell%ltpcel)
    call mpic%bcast(rgem,rbuffer)
    ibi=0;ibr=0
    cell%nox(1)=ibuffer(1); cell%nox(2)=ibuffer(2) ; cell%nox(3)=ibuffer(3)
    cell%noxyz=ibuffer(4)
    cell%natperc=ibuffer(5)
    cell%icaltabt=ibuffer(6)
    ibi=6
    do ip=1,nsize
       ibi=ibi+1
       cell%nato(ip)=ibuffer(ibi)
    end do
    do ip=1,nsize
       do ip2=0,cell%ncelvmax
          ibi=ibi+1
          cell%ncel(ip,ip2)=ibuffer(ibi)
       end do
    end do
    do ip=1,nsize
       do ip2=0,cell%ncelvmax
          do ip3=1,3
             ibi=ibi+1
             cell%deltadist(ip3,ip2,ip)=ibuffer(ibi)
          end do
       end do
    end do
    do ip=1,nsize
       do ip2=1,cell%natperc
          ibi=ibi+1
          cell%atincel(ip2,ip)=ibuffer(ibi)
       end do
    end do
#ifdef PARA
    do ip=1,nsize
       ibi=ibi+1
       cell%proc_cell(ip)=ibuffer(ibi)
    end do

#endif
    cell%celsize(1)=rbuffer(1);     cell%celsize(2)=rbuffer(2) ;    cell%celsize(3)=rbuffer(3)
    ibr=3
    if (cell%ltpcel) then
       do ip=1,nsize
          do ip2=1,3
             do ip3=1,3
                ibR=ibr+1
                cell%sigc(ip2,ip3,ip)=rbuffer(ibr)
             end do
          end do
       end do
       do ip=1,nsize
          ibr=ibr+1
          cell%tempc(ip)=rbuffer(ibr)
       end do
    end if

  end subroutine cells2a


  subroutine transfer_atoms(atcf,cellcf,indtrav,ntravtot,ntrav,proccib)
    class(atom_config):: atcf
    class(cell_config)::cellcf
    integer,allocatable::indtrav(:),proccib(:)
    integer::ntravtot
    integer,allocatable::ntrav(:)

    integer::nemp
    integer,allocatable::procvis(:) !,procem(:)=iproc

    !    integer, allocatable:: indtravtot(:),proccibtot(:),procemtot(:)
    integer::iproc,nprocs,natem,natrecv
    integer::iatem,j,i,k,rgcib,itrf,nato,icelj,jjj
    logical :: lwrk
    lwrk=.false.

    natem=0;natrecv=0
    nprocs=comm_space%nproc
    do iproc=0,nprocs-1
       if (ntrav(iproc).ne.0) then
          lwrk=.true.
       else
          lwrk=.false.
       end if
       call comm_space%bcast(iproc,lwrk)
       if (lwrk) then
          nemp=ntrav(iproc) ! ne sert que pour iproc, mais effacé la ligne suivant
          call comm_space%bcast(iproc,nemp)
          allocate(procvis(nemp))
          procvis=0
          procvis(1:nemp)=proccib(1:nemp) ! ne sert que pour iproc, mais effacé la ligne suivant
          call comm_space%bcast(iproc,procvis)
          do iatem=1,nemp
             if (myidsp==iproc) then  ! on est sur le proc émetteur
                natem=natem+1
                rgcib=proccib(iatem)
                itrf=indtrav(iatem)
                write(uwrt,*)'em,cib,itrf',iproc,rgcib,itrf
                icelj=atcf%ielat(atcf%im)
                call atcf%s1at2p_rm(rgcib,itrf,comm_space)
                do jjj=1,cellcf%nato(icelj)
                   if (cellcf%atincel(jjj,icelj)==atcf%im+1) cellcf%atincel(jjj,icelj)=itrf
                end do

             end if
             if (myidsp==procvis(iatem)) then
                write(uwrt,*)'recpt,em',myidsp,iproc,atcf%im
                call atcf%recv1at(iproc,comm_space,myidsp)
                atcf%proc_at(atcf%im)=myidsp
                write(uwrt,*)'ielat recv',myidsp, atcf%ielat(atcf%im),atcf%im
                cellcf%nato(atcf%ielat(atcf%im))= cellcf%nato(atcf%ielat(atcf%im))+1
                nato=cellcf%nato(atcf%ielat(atcf%im))
                cellcf%atincel(nato,atcf%ielat(atcf%im)) = atcf%im                
                natrecv=natrecv+1
             end if
          end do
          deallocate(procvis)
       end if
    end do
    call comm_space%sum(natrecv)
    call comm_space%sum(natem)
    if (natrecv.ne.natem) then
       write(uwrt,*)'natrecv<>natem',rang,myidsp,natrecv,natem
       call arret_ndm(.true.)
    end if
  end subroutine transfer_atoms


end module cellconfig


