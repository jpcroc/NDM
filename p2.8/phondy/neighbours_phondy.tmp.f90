!$---------------------------------------------
subroutine calc_neighbours_ph(iconf)
!$ compute neighbours using milady and  NDM style
use phondy_in_ndm_module, only : rangph
use derived_types_ph, only: config_real
implicit none
integer, intent(in) :: iconf
    !<<<starting the neighbours calculations
    if (rangph==0) write(6,'("PHONDY: the start of neighbours calc of:", a)') config_real(iconf)%filename
    call neighbours_ndm_layer(iconf)
    call neighbours_large(iconf)

return
end subroutine calc_neighbours_ph





!$---------------------------------------------
subroutine neighbours_ndm_layer(iconf)
! set-up the ndm using the values imported froim iconf
use phondy_in_ndm_module, only: rcut_ph, rcut_ph_ang
use gen_com_m, only: dmtype, im_glob, im,imm, at,bg, A2cm, rvois, umass
use derived_types_ph, only:config_real
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
!old_milady imm=config_real(iconf)%nat
!old_milady im=imm
im_glob=config_real(iconf)%im                       ! Mise à jour de im_glob pour divid
im=config_real(iconf)%im
imm=config_real(iconf)%imm
ntyp=config_real(iconf)%ntypes

!if (dmtype /= 18) then
!   xpp_copy(:,:)=xpp(:,:)
!   ax_copy(:,:)=ax(:,:)
!   ielat_copy(:)=ielat(:)
!   rumax_copy = rumax
!   na_copy(:ntyp) = na(:ntyp)
!end if
!write(*,*) 'IN NEIGH1', im, imm, ntyp,xp(1,1), config_real(iconf)%pos_cart(1,1)
!if (dmtype==18) then
   call alloc_all_tab_imm(im)
!else
!   call realloc_all_tab_imm(im)
!end if

!if (dmtype /= 18) then
!  ax(:,:) = ax_copy(:,:)
!  xpp(:,:)=xpp_copy(:,:)
!  ielat(:)=ielat_copy(:)
!end if

!write(*,*) 'IN NEIGH10', im, imm, ntyp,xp(1,1), config_real(iconf)%pos_cart(1,1)
ityp(:)=config_real(iconf)%itype(:)
xp(:,:)=config_real(iconf)%pos_cart(:,:)
fp(:,:)=config_real(iconf)%force(:,:)
at(:,:)=config_real(iconf)%cell(:,:)
bg(:,:)=config_real(iconf)%bg_cell(:,:)
!write(*,*) 'IN NEIGH2', im, imm, ntyp,xp(1,1), config_real(iconf)%pos_cart(1,1)
!C_DEBUG cutoff  for the neighbours list update.

!debug write(*,*) 'r_cut_ph', rcut_ph
!r_cut = r_cut_ph
rvois=rcut_ph_ang*A2cm
call convert_A2cm(1,.true.,.true.)
call alloc_typ_ph()

if (dmtype /= 18) then
   na(:ntyp) = na_copy(:ntyp)
   rumax = rumax_copy
end if
cm(:)=config_real(iconf)%mass_per_type(:)*umass
!write(*,*) 'IN NEIGH3', im, imm, ntyp,xp(1,1), config_real(iconf)%pos_cart(1,1)
return
end subroutine neighbours_ndm_layer
!$---------------------------------------------

!$---------------------------------------------
subroutine neighbours_large (iconf)
use gen_com_m, only: natperc, nox, noy,noz,im,rvois
use phondy_in_ndm_module, only: rangph
use tab_imm_m, ONLY : iwmax, iwmax2
implicit none
integer, intent(in) :: iconf
integer :: i
!if (debug) then
  if (rangph==0) write(6,'("PHONDY neighbours_large: debug ...neigh for iconf:", i7)') iconf
!end if
call Deallocatecel()
natperc=-1                        ! Force le calcul de natperc dans divid
nox=-1 ; noy=-1 ; noz=-1
!write(*,*) 'here 0', rvois
call divid(0)
!write(*,*) 'here 1', rvois
call divid(1)
!write(*,*) 'here 2', rvois
call DynamicalAllocationCell()    ! Réallocation du pointeur last(natperc,:noxyz) pour caltabt
call neigcel()
!write(*,*) 'after1', rvois
call caltabt()
!write(*,*) 'after2', rvois
call caltabi()
!write(*,*) 'after3', rvois

! do i=1,im
!     write(*,*) "test  neigh ...", i,iwmax(i), iwmax2(i)
! end do

return
end subroutine neighbours_large
!$---------------------------------------------


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

subroutine alloc_typ_ph

  use var_pot
  implicit none
  integer::i,j,k

     npair=ntyp*(ntyp+1)/2
     ntrip=ntyp*ntyp*(ntyp+1)/2
     rumax=0.
     if (associated(na)) deallocate(na) ; allocate(na(ntyp))

     if (associated(ipo)) deallocate(ipo); allocate(ipo(ntyp,ntyp))
     ! initialisation de ipo
     k = 0
     do i = 1, ntyp
        do j = i, ntyp
           k = k+1
           ipo(i,j) = k
           ipo(j,i) = k
        end do
     end do
      if (associated(cm))    deallocate(cm)    ; allocate(cm(ntyp))
      if (associated(catom)) deallocate(catom) ; allocate(catom(ntyp))

      if (associated(ty)) deallocate(ty); allocate(ty(ntyp))
      if (associated(pot)) deallocate(pot); allocate(pot(4,npair,0:ngrid+1))
!      if (associated(pot_d)) deallocate(pot_d); allocate(pot_d(4,npair,0:ngrid+1))

      if (associated(q)) deallocate(q); allocate(q(ntyp))
       q(:)=0
      if (associated(rc)) deallocate(rc); allocate(rc(ntyp))
     rc(1:ntyp)=rclu(1:ntyp)*1.d-8
      if (associated(zz)) deallocate(zz); allocate(zz(npair))
      if (associated(lue_paire)) deallocate(lue_paire); allocate(lue_paire(npair))
      if (associated(lu_roff_pair)) deallocate(lu_roff_pair); allocate(lu_roff_pair(npair))
      lu_roff_pair(:)=.false.

      if (associated(rue_pair)) deallocate(rue_pair); allocate(rue_pair(npair))
     rue_pair(:)=0.


      if (associated(lue_typ)) deallocate(lue_typ)
      allocate(lue_typ(npair))
      if (associated(typ_pot_pair)) deallocate(typ_pot_pair)
      allocate(typ_pot_pair(npair))
      if (associated(lue_trip)) deallocate(lue_trip)
      allocate(lue_trip(ntrip))
     lue_trip(:)=.false.;lue_paire(:)=.false.; lue_typ(:)=.false.
      if (associated(ro)) deallocate(ro) ; allocate(ro(npair))
      if (associated(dip)) deallocate(dip) ; allocate(dip(npair))
      if (associated(pm)) deallocate(pm) ; allocate(pm(npair))
      if (associated(roff1)) deallocate(roff1) ; allocate(roff1(npair))
      if (associated(roff2)) deallocate(roff2) ; allocate(roff2(npair))
      if (associated(a_factor)) deallocate(a_factor) ; allocate(a_factor(npair))
      if (associated(r8p)) deallocate(r8p) ; allocate(r8p(npair))

!     if(iterdf.ge.0) then
      if (associated(coord)) deallocate(coord)
         allocate(coord(ntyp,ntyp,nkmax))
      if (associated(digr)) deallocate(digr)
         allocate(digr(ntyp,ntyp,nkmax))
        digr=0.d0
      if (associated(fda)) deallocate(fda)
         allocate(fda(ntyp,ntyp,ntyp,contmax))

     if (associated(nad)) deallocate(nad)
     allocate(nad(ntyp))
     if (associated(nas)) deallocate(nas)
     allocate(nas(ntyp))
     if (associated(nai)) deallocate(nai)
     allocate(nai(ntyp))

!  end if



  return
end subroutine alloc_typ_ph

subroutine convert_A2cm(i,xa,xpos)
! i= 1 convert positions from A  -> cm
! i=-1 convert positions from cm -> A
use gen_com_m, only : A2cm,at,bg
use tab_imm_m, only : xp
use phondy_in_ndm_module, ONLY: rangph

implicit none

integer,intent(in) :: i
logical,intent(in) :: xa,xpos

if (i==1) then ! A2cm
   if (xa) then
         at(:,:)=at(:,:)*A2cm
         bg(:,:)=bg(:,:)/A2cm
   end if
   if (xpos) xp(:,:)=xp(:,:)*A2cm
elseif (i==-1) then ! cm2A
   if (xa) then
         at(:,:)=at(:,:)/A2cm
         bg(:,:)=bg(:,:)*A2cm
   end if
   if (xpos) xp(:,:)=xp(:,:)/A2cm
else
   if (rangph==0) write(*,*) "Wrong input : input = 1 or -1"
   stop "fatal in convert_A2cm"
endif


return
end subroutine convert_A2cm
