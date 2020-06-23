! nox,noy,noz, natperc doivent être connus pour l'instant

module cellconfig
  USE T_kind_param_m
  use atomconfig,only : atom_config


  implicit none
  !  integer:: incr=20 ! incrément des tailles de tableau 

  type cell_config
     integer:: nox,noy,noz,noxyz 
     integer::natperc !nombre (max) d'atomes par cellule
     integer(long)::icaltabt 
     integer,allocatable::nato (:) ! nombre d'atomes dans la cellule ko

     integer,allocatable:: ncel (:,:) ! ncel(ko,i1)=numéro de la ième cellule voisine de la cellule ko
     integer,allocatable:: atincel (:,:) ! atincel(i,k)= indice du ième atome de la cellule k
     integer,allocatable:: deltadist(:,:,:) !gestion des conditions périodiques entre les cellules (voisinage de bords de boites)

     real(double):: celsize(3)

#ifdef PARA
     integer :: cell_debx, cell_deby, cell_debz     !numero de la premiere cellule locale suivant x, y et z
     integer :: cell_finx, cell_finy, cell_finz     !numero de la derniere cellule locale  suivant x, y et z
     integer :: nb_cell_x, nb_cell_y, nb_cell_z     !nb de cel locales suivant x y z     
#endif     

   contains
     procedure, pass::init=>init_cel
     procedure, pass::dealloc=>dealloc_cel


  end type cell_config

!  type,extends(cell_config)::  cel_config_e ! des tas de tableaux annexes par toujours alloués car le plus souvent inutiles
!     integer::iewald !pilote les tableux de la sommation d'Ewald
!     real(double), allocatable:: tabv3,tabf3
!     real(double),allocatable::sigc(:,:,:)! contrainte par cellule
!     logical:: ltpcel ! écriture des quantités par cellules
!     real(double)::tempstopcel ! température d'arrêt de la cellule
!     integer::istopcel ! diverses versions de tempstopcel
!     real(double),allocatable,dimension(:):: eatcel,tempc,tempcm,celpm1,tm1,celpp,tcp,pmc
!     logical::lprtcel ! ecriture des résultats seulement sur certaines cellules.
!     
!     logical :: l2T ! cellule pour modèles à 2T
!     real(double),allocatable:: elossCel(:)
!     
!     logical :: lsigatcel ! moyenne des contraintes atomiques par cellule ?
!     integer, allocatable:: natchk(:)
!     real(double), allocatable::patcel(:),patcelmax, sigatcel(:,:,:)
!   contains
!     procedure, pass::init=>init_cel_e
!     procedure, pass::dealloc=>dealloc_cel_e
!     
!  end type cel_config_e
  
  
contains


  
  subroutine init_cel(cell)
    class(cell_config)::cell

    call allocatecelN(cell)

    call neigcelN(cell)
    return

  end subroutine init_cel
  subroutine allocatecelN(cell)
    class(cell_config)::cell
    integer::nsize
    cell%noxyz=cell%nox*cell%noy*cell%noz
    nsize=cell%noxyz
    allocate(cell%ncel(0:nsize,0:26))
    allocate(cell%nato(0:nsize))
    allocate(cell%atincel(cell%natperc,0:nsize))
    allocate(cell%deltadist(3,0:26,nsize))
    return
  end subroutine allocatecelN

  
  subroutine dealloc_cel(cell)
    class(cell_config)::cell

    deallocate(cell%ncel)
    deallocate(cell%nato)
    deallocate(cell%atincel)
    deallocate(cell%deltadist)


    return

  end subroutine dealloc_cel

  subroutine neigcelN(cell)
    class(cell_config)::cell

    integer :: kx, ky, kz, koo, l, lz, mz, ly, my, lx, mx, kxy


    if (cell%noxyz==1) then
       cell%ncel(1,0)=1
       cell%deltadist=0
    else

       do kz = 1, cell%noz
          do ky = 1, cell%noy
             do kx = 1, cell%nox
                koo = 1+(kx-1)+cell%nox*((ky-1)+cell%noy*(kz-1))

                cell%ncel(koo,0) = koo
                cell%deltadist(:,0,koo) = 0

                l = 1
                do lz = -1, 1
                   do ly = -1, 1
                      do lx = -1, 1
                         cell%deltadist(:,l,koo) = 0

                         mz = kz+lz
                         if (mz<1) then
                            mz = mz+cell%noz
                            cell%deltadist(3,l,koo) = 1
                         endif
                         if (mz>cell%noz) then
                            mz = mz-cell%noz
                            cell%deltadist(3,l,koo) = -1
                         endif

                         my = ky+ly
                         if (my<1) then
                            my = my+cell%noy
                            cell%deltadist(2,l,koo) = 1
                         endif
                         if (my>cell%noy) then
                            my = my-cell%noy
                            cell%deltadist(2,l,koo) = -1
                         endif

                         mx = kx-lx
                         if (mx<1) then
                            mx = mx+cell%nox
                            cell%deltadist(1,l,koo) = 1
                         endif
                         if (mx>cell%nox) then
                            mx = mx-cell%nox
                            cell%deltadist(1,l,koo) = -1
                         endif

                         kxy = 1+(mx-1)+cell%nox*((my-1)+cell%noy*(mz-1))
                         if (kxy==koo) cycle
                         cell%ncel(koo,l) = kxy
                         !                        write(6,*)koo,lz,ly,lx,l,kxy
                         !                        if ((kz==cell%noz).and.(lz==1))write(6,*)koo,lz,l,kxy
                         !                        if ((kz==1).and.(lz==-1))write(6,*)koo,lz,l,kxy
                         l = l+1
                      end do
                   end do
                end do
             end do
          end do
       end do
    end if
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


  subroutine caltabtC (cell,atcf,lperiod,bg)
    USE notperiod_mod,only: notperiod
    USE cryst_to_cart_mod,only: cryst_to_cart
    class(cell_config), intent(inout):: cell
    class(atom_config),intent(inout)::atcf
    logical,intent(in)::lperiod
    real(double),intent(in)::bg(3,3)

    integer :: i, ic, icell, kx, ky, kz, koo
    real(double) :: aux, auy, auz
    real(double), dimension(:,:), allocatable :: xpnp !
    integer(long), save:: icaltabt=0
    !
    ! --------- Initialisation --------------
    !
    !   write(6,*)'caltabt',it
    icaltabt=icaltabt+1
    
    cell%nato(0:cell%noxyz) = 0
    cell%atincel(1:cell%natperc,0:cell%noxyz) = 0

    !  -------- cas sans cellule  -----------
    if (cell%noxyz==1) then
       cell%nato(1) = atcf%im
       do i = 1, atcf%im
          atcf%ielat(i) = 1
          cell%atincel(i,1) = i
       end do
    else

       ALLOCATE(xpnp(3,atcf%im))        
       if (lperiod) then            
          xpnp(:,:)=atcf%xp(:,:)         
       else                         
          call notperiod(atcf%im,atcf%xp,xpnp)   
       end if
       !  -------- Initialisations  -----------

       ! -------------------------------------------
       !   1. loop: lattice-coordinates of all atoms
       ! - - - - - - - - - - - - - - - - - - - - - -

       !debug       write (*,*) 'sub caltabt 1',it,xp(1,1)
       call cryst_to_cart (atcf%im, xpnp, bg, -1) !cart vers cryst
       !debug       write (*,*) 'sub caltabt 2',it,xp(1,1)

       !     if (it.gt.1000) write(6,*)'CALTABT',it
       do i = 1, atcf%im
          !     if  ((it.ge.1000).and.(i.lt.20)) write(6,'(I5,3G15.7)')i, xpnp(1,i),xpnp(2,i),xpnp(3,i)
          aux = xpnp(1,i)*cell%nox
          auy = xpnp(2,i)*cell%noy
          auz = xpnp(3,i)*cell%noz
          kx = int(aux)
          ky = int(auy)
          kz = int(auz)
          !write(*,*) i, nox,noy,noz, kx,ky,kz
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
          atcf%ielat(i) = koo
          cell%nato(koo) = cell%nato(koo)+1
!write(*,*) MAXVAL(nato(:)),koo,i        ! DEBUG
          ! ==== MODIF Clouet =====================
          IF (cell%nato(koo).GT.cell%natperc) THEN
             WRITE(0,'(a)') 'You need to increase the maximal number of atoms per cell'
             WRITE(0,'(a,i0)') 'current value: natperc=', cell%natperc
             STOP '< Caltabt >'
          END IF
          ! ==== Fin MODIF Clouet =================
          cell%atincel(cell%nato(koo),koo) = i
       end do
       !debug            call cryst_to_cart (imm, xpnp, at, 1)  !cryst vers cart


       !        open(unit=809, file='cell.csv', form='formatted', &
       !             status='unknown')
       !    do i=1,im_glob
       !       write(809,'(A,3I6,3G22.13)')'Cel ', it,i,ielat(i),xp(:,i)
       !       write(809,'(A,3I6,3G22.13)')'Cel ', it,i,ielat(i),xpnp(:,i)
       !    end do
       !             do i=1,noxyz
       !                write(6,*) i, nato(i) 
       !             end do
       DEALLOCATE(xpnp)   ! MODIF CLOUET
    endif

!!$write(6,*)'sortie caltabt'     ! DEBUG

!      write(6,*)'maxnato', maxval(cell%nato)
    cell%icaltabt=icaltabt
    atcf%icaltabt=icaltabt
    return
  end subroutine caltabtC


  subroutine ndm2cellconfig(celndm,noxyz,nox,noy,noz,natperc,nato,ncel,atincel,deltadist,celsize)
    type(cell_config), intent(out):: celndm
    integer, intent(in):: nox,noy,noz,natperc,noxyz
    integer,intent(in)::ncel(0:noxyz,0:26),nato(0:noxyz),atincel(natperc,0:noxyz),deltadist(3,0:26,noxyz)
    real(double),intent(in)::celsize(3)
    
    celndm%nox=nox
    celndm%noy=noy
    celndm%noz=noz
    celndm%natperc=natperc
    celndm%noxyz=nox*noy*noz
    call allocatecelN(celndm)
    celndm%icaltabt=0

    celndm%ncel(0:noxyz,0:26)=ncel(0:noxyz,0:26)
    celndm%nato(0:noxyz)=nato(0:noxyz)
    celndm%atincel(1:natperc,0:noxyz)=atincel(1:natperc,0:noxyz)
    celndm%deltadist(1:3,0:26,1:noxyz)=deltadist(1:3,0:26,1:noxyz)
    celndm%celsize(1:3)=celsize(1:3)
  end subroutine ndm2cellconfig

  subroutine cellconfig2ndm(celndm,noxyz,nox,noy,noz,natperc,nato,ncel,atincel,deltadist,celsize)
    type(cell_config), intent(inout):: celndm
    integer, intent(inout):: nox,noy,noz,natperc,noxyz
    integer,intent(inout)::ncel(0:noxyz,0:26),nato(0:noxyz),atincel(natperc,0:noxyz),deltadist(3,0:26,noxyz)
    real(double),intent(out)::celsize(3)

    if ((nox.ne.celndm%nox).or.(noy.ne.celndm%noy).or.(noz.ne.celndm%noz).or.(natperc.ne.celndm%natperc)) then
       write(6,*)'incohérence entre noxyz et celndm%noxyz'
       stop
    end if
    
    celndm%icaltabt=0
    ncel(0:noxyz,0:26)=celndm%ncel(0:noxyz,0:26)
    nato(0:noxyz)=celndm%nato(0:noxyz)
    atincel(1:natperc,0:noxyz)=celndm%atincel(1:natperc,0:noxyz)
    deltadist(1:3,0:26,1:noxyz)=celndm%deltadist(1:3,0:26,1:noxyz)
    celsize(1:3)=celndm%celsize(1:3)

    call celndm%dealloc
    
  end subroutine cellconfig2ndm
end module cellconfig


