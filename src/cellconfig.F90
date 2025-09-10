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
     integer:: nox=0,noy=0,noz=0,noxyz=0 
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
    flx(1)=float(kox(1))/cell%nox
    flx(2)=float(kox(2))/cell%noy
    flx(3)=float(kox(3))/cell%noz
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
    nox=cell%nox;noy=cell%noy
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
       cell%nox=nox; cell%noy=noy; cell%noz=noz;cell%noxyz=nox*noy*noz
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
    !    integer,intent(in)::nox,noy,noz,natperc
    integer::nsize
    logical,optional::latomalloc
    logical::lata=.true.
    if (present(latomalloc))lata=latomalloc
    cell%noxyz=cell%nox*cell%noy*cell%noz
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
       do kz = 1, cell%noz
          do ky = 1, cell%noy
             do kx = 1, cell%nox
                koo = 1+(kx-1)+cell%nox*((ky-1)+cell%noy*(kz-1))
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
                            if((mx<1).or.(mx>cell%nox)) then !débordement
                               cycle loopin
                            end if
                         else
                            select case(cell%nox)
                            case(1)
                               if ((lx==-1).or.(lx==1)) cycle loopin
                            case(2)
                               if (lx==-1) cycle loopin
                            end select
                         end if
                         
                         if(box%ipbc(2).ne.1) then
                            if((my<1).or.(my>cell%noy)) then !débordement
                               cycle loopin
                            end if
                         else
                            select case(cell%noy)
                            case(1)
                               if ((ly==-1).or.(ly==1)) cycle loopin
                            case(2)
                               if (ly==-1) cycle loopin
                            end select
                         end if
                         
                         if(box%ipbc(3).ne.1) then
                            if((mz<1).or.(mz>cell%noz)) then !débordement
                               cycle loopin
                            end if
                         else
                            select case(cell%noz)
                            case(1)
                               if ((lz==-1).or.(lz==1)) cycle loopin
                            case(2)
                               if (lz==-1) cycle loopin
                            end select
                         end if

                         l=l+1 ! on est dans une vraie cellule voisine !
                         if (mz<1) then ! on ne peut pas être ici si ipbc(3).ne.1
                            mz = mz+cell%noz
                            cell%deltadist(3,l,koo) = 1
                         endif
                         if (mz>cell%noz) then
                            mz = mz-cell%noz
                            cell%deltadist(3,l,koo) = -1
                         endif

                         if (my<1) then
                            my = my+cell%noy
                            cell%deltadist(2,l,koo) = 1
                         endif
                         if (my>cell%noy) then
                            my = my-cell%noy
                            cell%deltadist(2,l,koo) = -1
                         endif

                         if (mx<1) then
                            mx = mx+cell%nox
                            cell%deltadist(1,l,koo) = 1
                         endif
                         if (mx>cell%nox) then
                            mx = mx-cell%nox
                            cell%deltadist(1,l,koo) = -1
                         endif

                         kxy = 1+(mx-1)+cell%nox*((my-1)+cell%noy*(mz-1))
!                         if (kxy==koo) cycle
                         cell%ncel(koo,l) = kxy
                         !                        write(6,*)koo,lz,ly,lx,l,kxy
                         !                        if ((kz==cell%noz).and.(lz==1))write(6,*)koo,lz,l,kxy
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
          
    !      if (ltranche) then
    !         do kz = 1, noz
    !            do ky = 1, noy
    !               do kx = 1, nox
    !                  if (kz==noz) then
    !                     koo = 1+(kx-1)+nox*((ky-1)+noy*(kz-1))
    !                     ncel(koo,1) = 0
    !                     ncel(koo,2) = 0
    !                     ncel(koo,3) = 0
    !                     ncel(koo,4) = 0
    !                     ncel(koo,5) = 0
    !                     ncel(koo,6) = 0
    !                     ncel(koo,7) = 0
    !                     ncel(koo,8) = 0
    !                     ncel(koo,9) = 0
    !                  endif
    !                  if (kz/=1) cycle
    !                  koo = 1+(kx-1)+nox*((ky-1)+noy*(kz-1))
    !                  ncel(koo,18) = 0
    !                  ncel(koo,19) = 0
    !                  ncel(koo,20) = 0
    !                  ncel(koo,21) = 0
    !                  ncel(koo,22) = 0
    !                  ncel(koo,23) = 0
    !                  ncel(koo,24) = 0
    !                  ncel(koo,25) = 0
    !                  ncel(koo,26) = 0
    ! 
    !               end do
    !            end do
    !         end do
    !     endif

    return
  end subroutine neigcelN


  subroutine caltabtC (cell,atcf,lperiod,boxcf,lextr,psc)
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
    logical,intent(in),optional::lextr
    logical::lextrait=.false.

    integer :: i,  kx, ky, kz, koo
    real(double) :: aux, auy, auz
    real(double), dimension(:,:), allocatable :: xpnp !
    integer(long), save:: icaltabt=0
    integer::iml
    integer :: ntrav,ntravtot,itrav
    integer,parameter::maxtrav=10
    integer,allocatable:: indtrav(:),proccib(:)
    ntrav=0
    !
    ! --------- Initialisation --------------
    !
    if (present(lextr))lextrait=lextr
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
          aux = xpnp(1,i)*cell%nox
          auy = xpnp(2,i)*cell%noy
          auz = xpnp(3,i)*cell%noz
          kx = int(aux)
          ky = int(auy)
          kz = int(auz)
!!$          write(*,*) i, cell%nox,cell%noy,cell%noz, kx,ky,kz
!!$          write(*,*) i, aux,auy,auz, xpnp(1,i), xpnp(2,i), xpnp(3,i)
!!$          write(6,*)
          kx = Modulo(kx,cell%nox)
          ky = Modulo(ky,cell%noy)
          kz = Modulo(kz,cell%noz)
          !          if  ((it.ge.1000).and.(i.lt.20))  write(6,'(I5,3G15.7)')i, kx,ky,kz
          !==============================================================
          koo = 1+kx+cell%nox*(ky+cell%noy*kz)
          IF ( (koo.GT.cell%noxyz).OR.(koo.LT.0) ) THEN
             WRITE(0,'(a,i0,a,3g20.12)') &
                  'Problem with atom ', i, ', x,y,z = ', atcf%xp(1:3,i)
             WRITE(0,'(2(a,i0))') ' koo = ', koo, ' - noxyz = ', cell%noxyz
             STOP
          END IF
          if ((present(psc)).and.(lspacendm)) then 
             if (cell%proc_cell(koo).ne.myidsp) then
!                write(6,*)'atout',i,atcf%num_at_glob(i),rang,cell%proc_cell(koo),koo

                if(.not.(any(psc%cell_ftm(:)==koo))) then
                   write(6,'(A,6I7)')'WARNING ::: attrrav:i natg ielat rangem rangf newcell',i,atcf%num_at_glob(i),atcf%ielat(i),rang,cell%proc_cell(koo),koo
                   write(6,'(A,6I7)')'WARNING ::: travelled from cell to cell ',cell%edge(atcf%ielat(i),boxcf) ,cell%edge(koo,boxcf)
                   ntrav=ntrav+1
                   if (ntrav==1) then
                      allocate(indtrav(maxtrav));allocate(proccib(maxtrav)); proccib=-1
                   end if
                   proccib(ntrav)=cell%proc_cell(koo)
                   indtrav(ntrav)=i
                   atcf%ielat(i) = koo
!                   cell%nato(koo) = cell%nato(koo)+1
!                   cell%atincel(cell%nato(koo),koo) = i                   
                   IF (cell%nato(koo).GT.cell%natperc) THEN
                   WRITE(0,'(a)') 'caltabtC : You need to increase the maximal number of atoms per cell'
                   WRITE(0,'(a,i0)') 'current value: natperc=', cell%natperc
                   STOP '< CaltabtC >'
                END IF


                   !                   write(6,*)'atom', i,atcf%num_at_glob(i),'in cell', koo, ' originally in proc', &
                   !                        &myidsp, 'now in ', cell%proc_cell(koo),' travelled too far. its cell is not a frontier cell',&
                   !                        &'RANG actuel = ',rang
                   !                   call arret_ndm
                else
                   if (cell%proc_cell(koo).ne.myidsp)       write(6,'(A,6I7)')'afrt:i natg ielat rangem rangf newcell',i,atcf%num_at_glob(i),atcf%ielat(i),rang,cell%proc_cell(koo),koo
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

       ntravtot=ntrav
#ifdef PARA
       call comm_space%sum(ntravtot)
       if (ntravtot.gt.0) then
          if (.not.allocated(proccib)) then
             allocate(proccib(maxtrav))
             proccib=-1
          end if
          call transfer_atoms(atcf,cell,indtrav,ntravtot,ntrav,proccib)
          
       end if
#endif

       
       DEALLOCATE(xpnp)   ! MODIF CLOUET
    endif

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
    return
  end subroutine caltabtC


!!$  subroutine ndm2cellconfig(celndm,box,noxyz,nox,noy,noz,natperc,nato,ncel,atincel,deltadist,celsize,ltpc,sigc,tempc,proc_cell)
!!$    type(cell_config), intent(out):: celndm
!!$    type(box_config)::box
!!$    integer, intent(in):: nox,noy,noz,natperc,noxyz
!!$    integer,intent(in)::ncel(noxyz,0:26),nato(noxyz),atincel(natperc,noxyz),deltadist(3,0:26,noxyz)
!!$    real(double),intent(in)::celsize(3)
!!$    logical,optional,intent(in)::ltpc
!!$    real(double),intent(in),optional,allocatable::sigc(:,:,:),tempc(:)
!!$    logical::ltpcel=.false.
!!$    integer,optional::proc_cell(:)
!!$    if (present (ltpc))ltpcel=ltpc
!!$    if (.not.allocated(celndm%ncel))then
!!$       call init_cel(celndm,box,nox,noy,noz,natperc,ltpcel)
!!$    end if
!!$    !    celndm%nox=nox
!!$    !    celndm%noy=noy
!!$    !    celndm%noz=noz
!!$    !    celndm%natperc=natperc
!!$    !    celndm%noxyz=nox*noy*noz
!!$    !    if (ltpcel)then
!!$    !       celndm%ltpcel=.true.
!!$    !    end if
!!$    !    call allocatecelN(celndm)
!!$!    celndm%icaltabt=0
!!$    celndm%ncel(1:noxyz,0:26)=ncel(1:noxyz,0:26)
!!$    celndm%nato(1:noxyz)=nato(1:noxyz)
!!$    celndm%atincel(1:natperc,1:noxyz)=atincel(1:natperc,1:noxyz)
!!$    celndm%deltadist(1:3,0:26,1:noxyz)=deltadist(1:3,0:26,1:noxyz)
!!$    celndm%celsize(1:3)=celsize(1:3)
!!$    if (ltpcel)then
!!$       celndm%sigc=sigc
!!$       celndm%tempc=tempc
!!$    end if
!!$#ifdef PARA
!!$    celndm%proc_cell=proc_cell
!!$#endif    
!!$
!!$  endsubroutine ndm2cellconfig
!!$
!!$  subroutine cellconfig2ndm(celndm,noxyz,nox,noy,noz,natperc,nato,ncel,atincel,deltadist,celsize,ltpcel,sigc,tempc,proc_cell)
!!$    type(cell_config), intent(inout):: celndm
!!$    integer, intent(inout):: nox,noy,noz,natperc,noxyz
!!$    integer,intent(inout),allocatable::ncel(:,:),nato(:),atincel(:,:),deltadist(:,:,:)!ncel(0:noxyz,0:26),nato(0:noxyz),atincel(natperc,0:noxyz),deltadist(3,0:26,noxyz)
!!$    real(double),intent(out)::celsize(3)
!!$    logical,optional,intent(out)::ltpcel
!!$    real(double),intent(out),optional,allocatable::sigc(:,:,:),tempc(:)
!!$    integer,optional,allocatable::proc_cell(:)
!!$    integer::nsize
!!$!    if (.not.allocated(ncel))then
!!$       nox=celndm%nox ;noy=celndm%noy;noz=celndm%noz
!!$       noxyz=nox*noy*noz;nsize=noxyz
!!$       natperc=celndm%natperc
!!$       if(allocated(ncel)) deallocate(ncel)
!!$       if(allocated(atincel)) deallocate(atincel)
!!$       if(allocated(deltadist)) deallocate(deltadist)
!!$       allocate(ncel(1:noxyz,0:26));allocate(atincel(natperc,1:noxyz));allocate(deltadist(3,0:26,1:noxyz))
!!$       if (allocated(nato)) deallocate(nato)
!!$       allocate(nato(noxyz))
!!$#ifdef PARA
!!$       if (present(proc_cell))then
!!$          if(allocated(proc_cell)) deallocate(proc_cell)
!!$          allocate(proc_cell(nsize))
!!$       end if
!!$#endif       
!!$!    else
!!$!       if ((nox.ne.celndm%nox).or.(noy.ne.celndm%noy).or.(noz.ne.celndm%noz).or.(natperc.ne.celndm%natperc)) then
!!$!          write(6,*)'incohérence entre noxyz et celndm%noxyz'
!!$!          write(6,*)nox,celndm%nox,natperc,celndm%natperc
!!$!#ifdef PARA
!!$!          call MPI_finalize(ierr)
!!$!#endif         
!!$!          call arret_ndm
!!$!       end if
!!$!        end if
!!$    if (celndm%ltpcel)then
!!$       ltpcel=celndm%ltpcel
!!$       if (ltpcel) then
!!$          if (.not.allocated(sigc))allocate (sigc(3,3,celndm%noxyz))
!!$          if (.not.allocated(tempc))allocate (tempc(celndm%noxyz))
!!$       end if
!!$    end if
!!$!    celndm%icaltabt=0
!!$    ncel(1:noxyz,0:26)=celndm%ncel(1:noxyz,0:26)
!!$    nato(1:noxyz)=celndm%nato(1:noxyz)
!!$    atincel(1:natperc,1:noxyz)=celndm%atincel(1:natperc,1:noxyz)
!!$    deltadist(1:3,0:26,1:noxyz)=celndm%deltadist(1:3,0:26,1:noxyz)
!!$    celsize(1:3)=celndm%celsize(1:3)
!!$    if (celndm%ltpcel)then
!!$       sigc(:,:,:)=celndm%sigc(:,:,:)
!!$       tempc(:)=celndm%tempc(:)
!!$    end if
!!$!    call celndm%dealloc
!!$#ifdef PARA
!!$    if (present(proc_cell))proc_cell=celndm%proc_cell
!!$#endif    
!!$  end subroutine cellconfig2ndm

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
    call cellcible%init(box,cellsource%nox,cellsource%noy,cellsource%noz,cellsource%natperc,ltpc=ltpcel,latomalloc=latcp)

    cellcible%nox=cellsource%nox
    cellcible%noy=cellsource%noy
    cellcible%noz=cellsource%noz
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
    call celcomp%init(box,celloc%nox,celloc%noy,celloc%noz,celloc%natperc,ltpc=celloc%ltpcel,latomalloc=latomcp)

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
    write(un,*)'nox noy noz noxyz',cellv%nox,cellv%noy,cellv%noz,cellv%noxyz
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
    ibuffer(1)=cell%nox; ibuffer(2)=cell%noy ; ibuffer(3)=cell%noz
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
    
    cell%nox=ibuffer(1); cell%noy=ibuffer(2) ; cell%noz=ibuffer(3)
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
    ibuffer(1)=cell%nox; ibuffer(2)=cell%noy ; ibuffer(3)=cell%noz
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
        cell%nox=ibuffer(1); cell%noy=ibuffer(2) ; cell%noz=ibuffer(3)
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


  subroutine transfer_atoms(atcf,cellcf,indtrav,ntravtot,ntrav,proccib)
    class(atom_config):: atcf
    class(cell_config)::cellcf
    integer,allocatable::indtrav(:),proccib(:)
    integer::ntravtot,ntrav

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
       if (myidsp==iproc) then
          if (ntrav.ne.0) then
             lwrk=.true.
          else
             lwrk=.false.
          end if
       end if
       call comm_space%bcast(iproc,lwrk)
       if (lwrk==.true.) then
          nemp=ntrav ! ne sert que pour iproc, mais effacé la ligne suivant
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
                write(6,*)'em,cib,itrf',iproc,rgcib,itrf
                icelj=atcf%ielat(atcf%im)
                call atcf%s1at2p_rm(rgcib,itrf,comm_space)
                do jjj=1,cellcf%nato(icelj)
                   if (cellcf%atincel(jjj,icelj)==atcf%im+1) cellcf%atincel(jjj,icelj)=itrf
                end do
                
             end if
             if (myidsp==procvis(iatem)) then
                write(6,*)'recpt,em',myidsp,iproc,atcf%im
                call atcf%recv1at(iproc,comm_space,myidsp)
                atcf%proc_at(atcf%im)=myidsp
                write(6,*)'ielat recv',myidsp, atcf%ielat(atcf%im),atcf%im
                cellcf%nato(atcf%ielat(atcf%im))= cellcf%nato(atcf%ielat(atcf%im))+1
                nato=cellcf%nato(atcf%ielat(atcf%im))
                cellcf%atincel(nato,atcf%ielat(atcf%im)) = atcf%im                
                natrecv=natrecv+1
             end if
          end do
          deallocate(procvis)
       end if
    end do
!    stop
    call comm_space%sum(natrecv)
    call comm_space%sum(natem)
    if (natrecv.ne.natem) then
       write(6,*)'natrecv<>natem',rang,myidsp,natrecv,natem
       call arret_ndm
    end if
! stop   
  end subroutine transfer_atoms
    
  
end module cellconfig


