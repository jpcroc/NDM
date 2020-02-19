module neighbours_mod
        use read_poscar_mod
        use alloc_typ_ml_mod
        use caltabt_mod
        implicit none
        contains
!$---------------------------------------------
subroutine test_if_config_is_small(iconf)
!test is a configuration is small comapred to r_cut.
!Input:
!         r_cut, config_real(iconf)
!Output:
!         config_real(iconf)%small
!                           T - box is small
!                           F - box is large
!$-------------------------------------------
use T_kind_param_m, only : double
use ml_in_ndm_module, only: r_cut, rangml, debug
use derived_types, only: config_real
use recips_mod
implicit none
integer, intent(in) :: iconf
real(double) :: bval(3)
integer :: i
call recips (config_real(iconf)%cell(1,1), config_real(iconf)%cell(1,2), config_real(iconf)%cell(1,3), &
             config_real(iconf)%bg_cell(1,1), config_real(iconf)%bg_cell(1,2), config_real(iconf)%bg_cell(1,3))

do i=1,3
   bval(i)=1.d0/sqrt(sum(config_real(iconf)%bg_cell(:,i)**2))
end do
config_real(iconf)%nxCell=int(2.d0*r_cut/bval(1))+1
config_real(iconf)%nyCell=int(2.d0*r_cut/bval(2))+1
config_real(iconf)%nzCell=int(2.d0*r_cut/bval(3))+1

if (max(config_real(iconf)%nxCell, config_real(iconf)%nyCell, config_real(iconf)%nzCell) .gt.1) config_real(iconf)%small=.true.

if (debug) then
    if (rangml==0) write(*,'("ML: small or big box in test_if_config_is_small ",  i5,  3i4, l3)') iconf,  config_real(iconf)%nxCell, &
                                                                                config_real(iconf)%nyCell, &
                                                                                config_real(iconf)%nzCell, &
                                                                                config_real(iconf)%small
end if
return
end subroutine test_if_config_is_small



!$---------------------------------------------
subroutine calc_neighbours(iconf)
!$ compute neighbours using milady and  NDM style
use ml_in_ndm_module, only : debug, rangml
use derived_types, only: config_real
implicit none
integer, intent(in) :: iconf
    !<<<starting the neighbours calculations
    if (debug) then
      if (rangml==0) write(6,'("ML: the start of neighbours calc of:", a)') config_real(iconf)%filename
    end if
    call neighbours_ndm_layer(iconf)
    if (config_real(iconf)%small) then
        call neighbours_small(iconf)
      else
        call neighbours_large(iconf)
     end if
    !<<<ending the neighbours calculations

return
end subroutine calc_neighbours





!$---------------------------------------------
subroutine neighbours_ndm_layer(iconf)
! set-up the ndm using the values imported froim iconf
use ml_in_ndm_module, only: r_cut
use gen_com_m, only: dmtype, im_glob, im,imm, at,bg, A2cm, rvois, umass
use derived_types, only:config_real
use var_pot, only: ntyp, cm, rumax, na
use tab_imm_m
implicit none
integer, intent(in) :: iconf
!convert xp, at from A to cm
real(kind(0.d0)), dimension(3,imm) :: xpp_copy, ax_copy
integer, dimension(imm) :: ielat_copy
integer, dimension(config_real(iconf)%ntypes) :: na_copy
real(kind(0.d0)) :: rumax_copy

!$$d call dealloc_all_tab_imm
imm=config_real(iconf)%nat
im=imm
im_glob=im                        ! Mise à jour de im_glob pour divid
ntyp=config_real(iconf)%ntypes

if (dmtype /= 18) then
   xpp_copy(:,:)=xpp(:,:)
   ax_copy(:,:)=ax(:,:)
   ielat_copy(:)=ielat(:)
   rumax_copy = rumax
   na_copy(:ntyp) = na(:ntyp)
end if
!write(*,*) 'IN NEIGH1', im, imm, ntyp,xp(1,1), config_real(iconf)%pos_cart(1,1)
if (dmtype==18) then
   call alloc_all_tab_imm(im)
else
   call realloc_all_tab_imm(im)
end if

if (dmtype /= 18) then
  ax(:,:) = ax_copy(:,:)
  xpp(:,:)=xpp_copy(:,:)
  ielat(:)=ielat_copy(:)
end if

!debug write(*,*) 'IN NEIGH10', im, imm, ntyp,xp(1,1), config_real(iconf)%pos_cart(1,1)
ityp(:)=config_real(iconf)%itype(:)
xp(:,:)=config_real(iconf)%pos_cart(:,:)
fp(:,:)=config_real(iconf)%force(:,:)
at(:,:)=config_real(iconf)%cell(:,:)
bg(:,:)=config_real(iconf)%bg_cell(:,:)
!debug write(*,*) 'IN NEIGH2', im, imm, ntyp,xp(1,1), config_real(iconf)%pos_cart(1,1)
!C_DEBUG cutoff  for the neighbours list update.
rvois=r_cut*A2cm
call convert_A2cm(1,.true.,.true.)
call alloc_typ_ml()

if (dmtype /= 18) then
   na(:ntyp) = na_copy(:ntyp)
   rumax = rumax_copy
end if
cm(:)=config_real(iconf)%mass_per_type(:)*umass
!debug write(*,*) 'IN NEIGH3', rumax, rvois, im, imm, ntyp,xp(1,1), config_real(iconf)%pos_cart(1,1)
return
end subroutine neighbours_ndm_layer
!$---------------------------------------------

!$---------------------------------------------
subroutine neighbours_large (iconf)
use gen_com_m, only: natperc, nox, noy, noz
use ml_in_ndm_module, only: debug, rangml
use dynalloccell
use divid_mod
use neigcel_mod
use caltabi_mod
implicit none
integer, intent(in) :: iconf
if (debug) then
  if (rangml==0) write(6,'("ML: debug ...neigh for iconf:", i7)') iconf
end if
call Deallocatecel()
nox=-1;  noy=-1; noz=-1
natperc=-1                        ! Force le calcul de natperc dans divid
call divid(0)
call divid(1)
!write(*,*) 'here', rvois
call DynamicalAllocationCell()    ! Réallocation du pointeur last(natperc,:noxyz) pour caltabt
call neigcel()
!write(*,*) 'after1', rvois
call caltabt()
!write(*,*) 'after2', rvois
call caltabi()
!write(*,*) 'after3', rvois
return
end subroutine neighbours_large
!$---------------------------------------------

!$---------------------------------------------
subroutine neighbours_small(iconf)
use T_kind_param_m, ONLY:  double
use gen_com_m, only : im
use ml_in_ndm_module, only : imm_neigh, r_cut
use derived_types, only : config_real
implicit none

integer, intent (in) :: iconf

integer :: ia,ja,c, k,n1,n2,n3, ibox
double precision :: r2, r_cut2
real(double), dimension(:,:), allocatable :: xpnp
integer :: nsize1, nsize2, nsize3
real(double), dimension(1:3) :: utemp

if  (allocated(config_real(iconf)%n_neigh)) deallocate(config_real(iconf)%n_neigh) ; allocate(config_real(iconf)%n_neigh(im))
if  (allocated(config_real(iconf)%r_ij)) deallocate(config_real(iconf)%r_ij) ; allocate(config_real(iconf)%r_ij(im,imm_neigh))
if  (allocated(config_real(iconf)%u_per)) deallocate(config_real(iconf)%u_per) ; allocate(config_real(iconf)%u_per(im,imm_neigh,3))
if  (allocated(config_real(iconf)%type_neigh)) deallocate(config_real(iconf)%type_neigh) ; allocate(config_real(iconf)%type_neigh(im,imm_neigh))
if  (allocated(config_real(iconf)%kind_neigh)) deallocate(config_real(iconf)%kind_neigh) ; allocate(config_real(iconf)%kind_neigh(im,imm_neigh))

if  (allocated(config_real(iconf)%u_ij)) deallocate(config_real(iconf)%u_ij) ; allocate(config_real(iconf)%u_ij(im,imm_neigh,3))

ALLOCATE(xpnp(3,config_real(iconf)%nat))
!if (lperiod) then
 xpnp(:,:)=config_real(iconf)%pos_cart
!else
! call notperiod(xp,xpnp)
!end if
! call cryst_to_cart (imm, xpnp, config_real(iconf)%bg_cell, -1)


nsize1=int(config_real(iconf)%nxCell/2)+1
nsize2=int(config_real(iconf)%nyCell/2)+1
nsize3=int(config_real(iconf)%nzCell/2)+1


r_cut2 = r_cut**2
!debug write(*,*) 'test', r_cut, config_real(iconf)%cell(:,1)
do ia=1,config_real(iconf)%nat
   c=0
   ibox=0
   do ja=1,config_real(iconf)%nat
      !$ dxp_ji(:)=config_real(iconf)%pos_cart(:,ia) - config_real(iconf)%pos_cart(:,ja)

      !$ ds(:)=MatMul(dxp_ji,config_real(iconf)%bg_cell)
      !$ WHERE ( (ds.GT.0.5d0).OR.(ds.LT.-0.5d0) )
      !$      ds(1:3) = ds(1:3) - Dble(Nint(ds(1:3)))
      !$ END WHERE
      !$ dxp_ji(:) = MatMul(config_real(iconf)%cell,ds)
      !$ r2_ji = SUM(dxp_ji(:)**2)

      do n1=-nsize1,nsize1
         do n2=-nsize2,nsize2
            do n3=-nsize3,nsize3
              r2=0.d0
              do k=1,3
                 utemp(k)=xpnp(k,ia)-(xpnp(k,ja)+config_real(iconf)%cell(k,1)*dble(n1)+config_real(iconf)%cell(k,2)*dble(n2)+config_real(iconf)%cell(k,3)*dble(n3))
                 r2=r2+utemp(k)**2
              enddo
              if (r2.lt.1d-15) cycle
              if (r2.gt.r_cut2) cycle
              c=c+1
              !$ if ( dabs(r2 - r2_ji) .lt.1.d-15 ) then
              !$   if (n1==0).or.(n2==0)
              !$  ibox = ibox+1
              !$    ! debug write(*,'("nnn", 4i5)') ia, ja, ibox, c
              !$ end if
              config_real(iconf)%type_neigh(ia,c)=config_real(iconf)%itype(ja)
              config_real(iconf)%kind_neigh(ia,c)=ja
              config_real(iconf)%r_ij(ia,c)=dsqrt(r2)
              config_real(iconf)%u_per(ia,c,:) = utemp(:) + xpnp(:,ja) - xpnp(:,ia)
              config_real(iconf)%u_ij(ia,c,:)=-utemp(:)!/config_real(iconf)%r_ij(ia,c)
            enddo
         enddo
      enddo
   enddo
   config_real(iconf)%n_neigh(ia)=c
enddo
deallocate(xpnp)
return

end subroutine neighbours_small


subroutine  calc_volume (a1,a2,a3, volume)
 implicit none
 real(kind(0.d0)), dimension(3), intent(in) :: a1, a2, a3
 real(kind(0.d0)), intent(out) :: volume
 integer :: i, j, k, l, s, iperm
 !-----------------------------------------------
 volume = 0.0
 i = 1
 j = 2
 k = 3
 s = 1.D0
 do iperm = 1, 3
    volume = volume+s*a1(i)*a2(j)*a3(k)
    l = i
    i = j
    j = k
    k = l
 end do
 i = 2
 j = 1
 k = 3
 s = -s
 do while(s<0.D0)
    do iperm = 1, 3
       volume = volume+s*a1(i)*a2(j)*a3(k)
       l = i
       i = j
       j = k
       k = l
    end do
    i = 2
    j = 1
    k = 3
    s = -s
 end do
 volume=dabs(volume)
 return
end subroutine calc_volume
!$---------------------------------------------



subroutine deallocate_real_config(ifile)
  use ml_in_ndm_module, only: ml_type, ml_type_descriptors
  use derived_types, only: config_real
  implicit none
  integer, intent(in) :: ifile


   if (allocated(config_real(ifile)%itype)) deallocate(config_real(ifile)%itype)
   if (allocated(config_real(ifile)%pos_cart)) deallocate(config_real(ifile)%pos_cart)
   if (allocated(config_real(ifile)%pos_crst)) deallocate(config_real(ifile)%pos_crst)
   if (allocated(config_real(ifile)%force))    deallocate(config_real(ifile)%force)
   if (allocated(config_real(ifile)%atomic_spin))    deallocate(config_real(ifile)%atomic_spin)


   if (allocated(config_real(ifile)%mass_per_type))  deallocate(config_real(ifile)%mass_per_type)
   if (allocated(config_real(ifile)%weight_per_type))  deallocate(config_real(ifile)%weight_per_type)
   if (allocated(config_real(ifile)%Z_per_type))  deallocate(config_real(ifile)%Z_per_type)
   if (allocated(config_real(ifile)%fix_type_poscar_to_periodic)) deallocate(config_real(ifile)%fix_type_poscar_to_periodic)


   if (allocated(config_real(ifile)%proc_atom)) deallocate(config_real(ifile)%proc_atom)
   if (ml_type==ml_type_descriptors) then
      !if (allocated(config_real(ifile)%stress)) deallocate(config_real(ifile)%stress)
      if (allocated(config_real(ifile)%type_neigh)) deallocate(config_real(ifile)%type_neigh)
      if (allocated(config_real(ifile)%kind_neigh)) deallocate(config_real(ifile)%kind_neigh)
      if (allocated(config_real(ifile)%n_neigh)) deallocate(config_real(ifile)%n_neigh)
      if (allocated(config_real(ifile)%r_ij))  deallocate(config_real(ifile)%r_ij)
      if (allocated(config_real(ifile)%u_ij))  deallocate(config_real(ifile)%u_ij)
      if (allocated(config_real(ifile)%u_per))  deallocate(config_real(ifile)%u_per)
   end if


  return
end subroutine deallocate_real_config
end module

