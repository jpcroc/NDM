module cellconfig
  USE arret_ndm_mod,only:arret_ndm
  USE T_kind_param_m
  use atomconfig,only : atom_config,atom_config_d,atom_config_e
  use boxconfig,only:box_config
  use paraconfig,only:para_config
  use Tpara,only:para_space_config,mpi_communicator,myidsp,comm_space
  use gen_com_m,only:rang
  implicit none
  !  integer:: incr=20 ! incrément des tailles de tableau 

  type cell_config
     integer:: nox(3)=0,noxyz=0 
     integer::natperc !nombre (max) d'atomes par cellule
     integer(long)::icaltabt 
     integer,allocatable::nato (:) ! nombre d'atomes dans la cellule ko
     integer,allocatable:: ncel (:,:) ! ncel(ko,i1)=numéro de la ième cellule voisine de la cellule ko
     integer,allocatable:: atincel (:,:) ! atincel(i,k)= indice du ième atome de la cellule k
     integer,allocatable:: deltadist(:,:,:) !gestion des conditions périodiques entre les cellules (voisinage de bords de boites)
     integer,allocatable:: ncelvois(:) !nombre de cellules voisines de la cellule actuelle (26 pour PBC standard; dépend de la position dans la boite pour les PBC partielles
     logical :: ltpcel
     logical:: ismall(3)=.false.
     real(double),allocatable::sigc(:,:,:),tempc(:)
     real(double):: celsize(3)
!     integer:: ngx(3)=0,ngxyz=0 
#ifdef PARA
     integer,allocatable::proc_cell(:)
     
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
    flx(:)=float(kox(:))/cell%nox(:)
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
    !write(6,*) 'nox', cell%nox
    cell%ltpcel=ltpcel
    call dealloc_cel(cell)
    call allocatecelN(cell,lata)
    call neigcelN(cell,box)
    !    call cell%print
       
    return

  end subroutine init_cel

  subroutine allocatecelN(cell,latomalloc)
    class(cell_config)::cell
    integer::nsize
    logical,optional::latomalloc
    logical::lata=.true.
    if (present(latomalloc))lata=latomalloc
    cell%noxyz=cell%nox(1)*cell%nox(2)*cell%nox(3)
    nsize=cell%noxyz
    if (nsize.ne.0) then
       allocate(cell%ncel(nsize,0:26))
       cell%ncel=0
       allocate(cell%nato(nsize))
       cell%nato=0
       select type(cell)
       class is (cell_config_g)
          allocate(cell%natotot(nsize))
          cell%natotot=0
       end select
       allocate(cell%ncelvois(nsize))
       allocate(cell%deltadist(3,0:26,nsize))
       cell%deltadist=0
            
       
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
          cell%proc_cell=-1
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

    if (cell%noxyz==1) then
       cell%ncel(1,0)=1
       cell%deltadist=0
       cell%ncelvois(1)=0
    else
       cell%deltadist(:,:,:) = 0 ! par défaut celldeltadist=0
       do kz = 1, cell%nox(3)
          do ky = 1, cell%nox(2)
             do kx = 1, cell%nox(1)
                koo = 1+(kx-1)+cell%nox(1)*((ky-1)+cell%nox(2)*(kz-1))
                cell%ncel(koo,0) = koo
                !                cell%deltadist(:,0,koo) = 0
!                cell%ncelvois(koo)= min(celcf%noxyz,27)-1
                l = 0 

                do lz = -1, 1
                   do ly = -1, 1
                      loopin:  do lx = -1, 1
                         if((lx==0).and.(ly==0).and.(lz==0)) cycle loopin !cellule en cours d'analyse (lignes au dessus)
                         mz = kz+lz
                         mx = kx+lx
                         my = ky+ly
                         if(box%ipbc(1).ne.1) then
                            if((mx<1).or.(mx>cell%nox(1))) then !débordement
                               cycle loopin
                            end if
                         else
                            select case(cell%nox(1))
                            case(1)
                               if ((lx==-1).or.(lx==1)) cycle loopin
                            case(2)
                               if (lx==-1) cycle loopin
                            end select
                         end if
                         
                         if(box%ipbc(2).ne.1) then
                            if((my<1).or.(my>cell%nox(2))) then !débordement
                               cycle loopin
                            end if
                         else
                            select case(cell%nox(2))
                            case(1)
                               if ((ly==-1).or.(ly==1)) cycle loopin
                            case(2)
                               if (ly==-1) cycle loopin
                            end select
                         end if
                         
                         if(box%ipbc(3).ne.1) then
                            if((mz<1).or.(mz>cell%nox(3))) then !débordement
                               cycle loopin
                            end if
                         else
                            select case(cell%nox(3))
                            case(1)
                               if ((lz==-1).or.(lz==1)) cycle loopin
                            case(2)
                               if (lz==-1) cycle loopin
                            end select
                         end if

                         l=l+1 ! on est dans une vraie cellule voisine !
                         if (mz<1) then ! on ne peut pas être ici si ipbc(3).ne.1
                            mz = mz+cell%nox(3)
                            cell%deltadist(3,l,koo) = 1
                         endif
                         if (mz>cell%nox(3)) then
                            mz = mz-cell%nox(3)
                            cell%deltadist(3,l,koo) = -1
                         endif

                         if (my<1) then
                            my = my+cell%nox(2)
                            cell%deltadist(2,l,koo) = 1
                         endif
                         if (my>cell%nox(2)) then
                            my = my-cell%nox(2)
                            cell%deltadist(2,l,koo) = -1
                         endif

                         if (mx<1) then
                            mx = mx+cell%nox(1)
                            cell%deltadist(1,l,koo) = 1
                         endif
                         if (mx>cell%nox(1)) then
                            mx = mx-cell%nox(1)
                            cell%deltadist(1,l,koo) = -1
                         endif

                         kxy = 1+(mx-1)+cell%nox(1)*((my-1)+cell%nox(2)*(mz-1))
!                         if (kxy==koo) cycle
                         cell%ncel(koo,l) = kxy
                         !                        write(6,*)koo,lz,ly,lx,l,kxy
                         !                        if ((kz==cellnox(3)).and.(lz==1))write(6,*)koo,lz,l,kxy
                         !                        if ((kz==1).and.(lz==-1))write(6,*)koo,lz,l,kxy
                      end do loopin
                   end do
                end do
                cell%ncelvois(koo)=l
             end do
          end do
       end do
    end if
!!$    do kx=1,cell%noxyz
!!$       write(6,*)'ncelvois', cell%ncelvois(kx)
!!$       write(6,*)cell%ncel(kx,:)
!!$    end do
          

    return
  end subroutine neigcelN


  subroutine caltabtC (cell,atcf,lperiod,boxcf,lextr,psc,lcheckfrontier)
    USE notperiod_mod,only: notperiod
    USE cryst_to_cart_mod,only: cryst_to_cart
    use gen_com_m,only:lspacendm
#ifdef PARA
  USE Tpara,only:comm_space
#endif
    
    class(cell_config), intent(inout):: cell
    class(atom_config),intent(inout)::atcf
    class(box_config),intent(inout)::boxcf
    type(para_space_config),optional::psc    
    logical,intent(in)::lperiod
    logical,intent(in),optional::lextr,lcheckfrontier
    logical::lextrait=.false.,lchkftr=.true.

    integer :: i,  kx, ky, kz, koo
    real(double) :: aux, auy, auz
    real(double), dimension(:,:), allocatable :: xpnp !
    integer(long), save:: icaltabt=0
    integer::iml
    !
    ! --------- Initialisation --------------
    !
    if (present(lextr))lextrait=lextr
    if (present(lcheckfrontier))lchkftr=lcheckfrontier
    if (lextrait) then
       iml=atcf%imm
    else
       iml=atcf%im
    end if
    icaltabt=icaltabt+1
!       write(6,*)'caltabt',icaltabt

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
          write(6,*)'PLANTE',rang
          write(300+RANG,*)'PLANTE'
          do i=1,iml
             if (any(xpnp(:,i).gt.1).or.any(xpnp(:,i).lt.0))  write(6,*) i,xpnp(:,i)
          end do
          write(6,*)'caltabtc xpnp <0 ou >1 stop'
          call arret_ndm
       end if
       !debug       write (*,*) 'sub caltabt 2',it,xp(1,1)

       !     if (it.gt.1000) write(6,*)'CALTABT',it
       !       write(6,*)'caltabt icaltabt im',icaltabt,atcf%im

       do i = 1, iml
          !     if  ((it.ge.1000).and.(i.lt.20)) write(6,'(I5,3G15.7)')i, xpnp(1,i),xpnp(2,i),xpnp(3,i)
          aux = xpnp(1,i)*cell%nox(1)
          auy = xpnp(2,i)*cell%nox(2)
          auz = xpnp(3,i)*cell%nox(3)
          kx = int(aux)
          ky = int(auy)
          kz = int(auz)
          kx = Modulo(kx,cell%nox(1))
          ky = Modulo(ky,cell%nox(2))
          kz = Modulo(kz,cell%nox(3))
          !          if  ((it.ge.1000).and.(i.lt.20))  write(6,'(I5,3G15.7)')i, kx,ky,kz
          !==============================================================
          koo = 1+kx+cell%nox(1)*(ky+cell%nox(2)*kz)
          IF ( (koo.GT.cell%noxyz).OR.(koo.LT.0) ) THEN
             WRITE(0,'(a,i0,a,3g20.12)') &
                  'Problem with atom ', i, ', x,y,z = ', atcf%xp(1:3,i)
             WRITE(0,'(2(a,i0))') ' koo = ', koo, ' - noxyz = ', cell%noxyz
             STOP
          END IF
          atcf%ielat(i) = koo
          cell%nato(koo) = cell%nato(koo)+1
!                    write(6,*)'caltabt', koo,i,cell%nato(koo)        ! DEBUG
          ! ==== MODIF Clouet =====================
          IF (cell%nato(koo).GT.cell%natperc) THEN
             WRITE(0,'(a)') 'caltabtC : You need to increase the maximal number of atoms per cell'
             WRITE(0,'(a,i0)') 'current value: natperc=', cell%natperc
             STOP '< CaltabtC >'
          END IF
          ! ==== Fin MODIF Clouet =================
          cell%atincel(cell%nato(koo),koo) = i
!          write(6,*)i,koo
#ifdef PARA
          if ((present(psc)).and.(lspacendm).and.(lchkftr)) then 
             if (cell%proc_cell(koo).ne.myidsp) then
                if(.not.(any(psc%cell_ftm(:)==koo))) then
                   write(6,*)'atom', i,atcf%num_at_glob(i),'in cell', koo, ' originally in proc', &
                        &myidsp, 'now in ', cell%proc_cell(koo),' travelled too far. its cell is not a frontier cell',&
                        &'RANG actuel = ',rang
                   call arret_ndm
                end if
             end if
          end if

#endif
          
          
       end do
       !debug            call cryst_to_cart (imm, xpnp, at, 1)  !cryst vers cart


       !        open(unit=809, file='cell.csv', form='formatted', &
       !             status='unknown')
       !             do i=1,noxyz
       !                write(6,*) i, nato(i) 
       !             end do
       DEALLOCATE(xpnp)   ! MODIF CLOUET
    endif
    !    do koo=1,cell%noxyz
    !       write(6,*)'nato',koo,cell%nato(koo)
    !    end do

!!$write(6,*)'sortie caltabt'     ! DEBUG

    !      write(6,*)'maxnato', maxval(cell%nato)
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

!    call atcf%print
!    call cell%print
    return
  end subroutine caltabtC



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

    cellcible%nox=cellsource%nox
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
    integer::i,un
    character(len=*),optional::mess
    un=6
    if (present(unit))un=unit
    write(un,*)'in cellprint ',mess
    write(un,*)'nox noy noz noxyz',cellv%nox,cellv%nox(2),cellv%nox(3),cellv%noxyz
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
!    write(6,*)'deltadist',cellv%deltadist
#ifdef PARA
    if (allocated(cellv%proc_cell))then 
       do i=1,cellv%noxyz
          write(un,*)'proc_cell',i,cellv%proc_cell(i)
       end do
    end if
#endif
    
    
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
       do ip2=0,26
          ibi=ibi+1
          ibuffer(ibi)=cell%ncel(ip,ip2)
       end do
    end do
    do ip=1,nsize
       do ip2=0,26
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
!!$    ibi=ibi+1; ibuffer(ibi)=cell%cell_debx
!!$    ibi=ibi+1; ibuffer(ibi)=cell%cell_deby
!!$    ibi=ibi+1; ibuffer(ibi)=cell%cell_debz
!!$
!!$    ibi=ibi+1; ibuffer(ibi)=cell%cell_finx
!!$    ibi=ibi+1; ibuffer(ibi)=cell%cell_finy
!!$    ibi=ibi+1; ibuffer(ibi)=cell%cell_finz
!!$
!!$    ibi=ibi+1; ibuffer(ibi)=cell%nb_cell_x
!!$    ibi=ibi+1; ibuffer(ibi)=cell%nb_cell_y
!!$    ibi=ibi+1; ibuffer(ibi)=cell%nb_cell_z
   
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
    
    cell%nox=ibuffer(1); cell%nox(2)=ibuffer(2) ; cell%nox(3)=ibuffer(3)
    cell%noxyz=ibuffer(4)
    cell%natperc=ibuffer(5)
    cell%icaltabt=ibuffer(6)
    ibi=6
    do ip=1,nsize
       ibi=ibi+1
       cell%nato(ip)=ibuffer(ibi)
    end do
    do ip=1,nsize
       do ip2=0,26
          ibi=ibi+1
          cell%ncel(ip,ip2)=ibuffer(ibi)
       end do
    end do
    do ip=1,nsize
       do ip2=0,26
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
!!$    ibi=ibi+1; cell%cell_debx=ibuffer(ibi)
!!$    ibi=ibi+1; cell%cell_deby=ibuffer(ibi)
!!$    ibi=ibi+1; cell%cell_debz=ibuffer(ibi)
!!$
!!$    ibi=ibi+1; cell%cell_finx=ibuffer(ibi)
!!$    ibi=ibi+1; cell%cell_finy=ibuffer(ibi)
!!$    ibi=ibi+1; cell%cell_finz=ibuffer(ibi)
!!$
!!$    ibi=ibi+1; cell%nb_cell_x=ibuffer(ibi)
!!$    ibi=ibi+1; cell%nb_cell_y=ibuffer(ibi)
!!$    ibi=ibi+1; cell%nb_cell_z=ibuffer(ibi)
   
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
!    write(6,*)'sizes',sizer,sizei,cell%ltpcel
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
       do ip2=0,26
          ibi=ibi+1
          ibuffer(ibi)=cell%ncel(ip,ip2)
       end do
    end do
    do ip=1,nsize
       do ip2=0,26
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
!!$    ibi=ibi+1; ibuffer(ibi)=cell%cell_debx
!!$    ibi=ibi+1; ibuffer(ibi)=cell%cell_deby
!!$    ibi=ibi+1; ibuffer(ibi)=cell%cell_debz
!!$
!!$    ibi=ibi+1; ibuffer(ibi)=cell%cell_finx
!!$    ibi=ibi+1; ibuffer(ibi)=cell%cell_finy
!!$    ibi=ibi+1; ibuffer(ibi)=cell%cell_finz
!!$
!!$    ibi=ibi+1; ibuffer(ibi)=cell%nb_cell_x
!!$    ibi=ibi+1; ibuffer(ibi)=cell%nb_cell_y
!!$    ibi=ibi+1; ibuffer(ibi)=cell%nb_cell_z
   
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
       do ip2=0,26
          ibi=ibi+1
          cell%ncel(ip,ip2)=ibuffer(ibi)
       end do
    end do
    do ip=1,nsize
       do ip2=0,26
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
!!$    ibi=ibi+1; cell%cell_debx=ibuffer(ibi)
!!$    ibi=ibi+1; cell%cell_deby=ibuffer(ibi)
!!$    ibi=ibi+1; cell%cell_debz=ibuffer(ibi)
!!$
!!$    ibi=ibi+1; cell%cell_finx=ibuffer(ibi)
!!$    ibi=ibi+1; cell%cell_finy=ibuffer(ibi)
!!$    ibi=ibi+1; cell%cell_finz=ibuffer(ibi)
!!$
!!$    ibi=ibi+1; cell%nb_cell_x=ibuffer(ibi)
!!$    ibi=ibi+1; cell%nb_cell_y=ibuffer(ibi)
!!$    ibi=ibi+1; cell%nb_cell_z=ibuffer(ibi)
   
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


    

end module cellconfig


