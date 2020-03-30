


!$-------------------------------------------------------------
subroutine md_init_config_ph
!$-------------------------------------------------------------
! Here the ml potential is initialized. This subroutine is used
! by MD program and is called in NDM's init.F90 sub after that configuration was read
use phondy_in_ndm_module, only: iconf_ini, iconf_big, nat, matfor, debug_ph,rangph
use gen_com_m, only: im, imm, volu, A2cm
use derived_types_ph, only:config_real
use var_pot, only: ntyp
implicit none
integer :: md_iconf

if (debug_ph) then
    if (rangph==0 ) write(6,*) '>...enter in md_init_config_ph'
end if

! config_real and condif_desc objects
! In phondy mode the config_real has the dimesion 2: 1 for ini, 2 for big boxes
if (allocated(config_real)) deallocate(config_real) ; allocate(config_real(2))

if (imm .lt. nat) then
!  imm=im
  write(6,'("PHONDY error: imm in din file should be bigger than number of atoms.")')
  stop 'md_init_config_ph imm lower than nat'
end if

md_iconf=1
config_real(md_iconf)%volume = volu/A2cm**3

do md_iconf=1,2

  if (allocated(config_real(md_iconf)%itype))    deallocate(config_real(md_iconf)%itype)    ;  allocate(config_real(md_iconf)%itype(imm))
  if (allocated(config_real(md_iconf)%pos_cart)) deallocate(config_real(md_iconf)%pos_cart) ;  allocate(config_real(md_iconf)%pos_cart(3,imm))
  if (allocated(config_real(md_iconf)%pos_crst)) deallocate(config_real(md_iconf)%pos_crst) ;  allocate(config_real(md_iconf)%pos_crst(3,imm))
  if (allocated(config_real(md_iconf)%force))    deallocate(config_real(md_iconf)%force)    ;  allocate(config_real(md_iconf)%force(3,imm))
  if (allocated(config_real(md_iconf)%mass_per_type))    &
                                               deallocate(config_real(md_iconf)%mass_per_type)  &
                                                                                          ;  allocate(config_real(md_iconf)%mass_per_type(ntyp))
end do

! This is initial set-up of im and imm imm_setup
config_real(iconf_ini)%im=im
config_real(iconf_ini)%imm=imm
config_real(iconf_ini)%ntypes = ntyp


if (debug_ph) then
    if (rangph==0 ) write(6,*) '>...exit in md_init_config_ph'
end if

return
end subroutine md_init_config_ph
!<-------------------------------------------------------------


!$-------------------------------------------------------------
subroutine put_ph_config_into_ndm (iconf)
!$-------------------------------------------------------------
use gen_com_m, only: at, bg , im, imm, A2cm, erg2ev, umass, volu
use tab_imm_m, only : xp,  fp,  ityp
use derived_types_ph, only:config_real
use var_pot, only: ntyp, cm
use phondy_in_ndm_module, only: debug_ph, rangph
implicit none
integer, intent(in) :: iconf


if (debug_ph) then
    if (rangph==0 ) write(6,*) '>...enter in  put_ph_config_into_ndn '
end if

imm = config_real(iconf)%imm
im= config_real(iconf)%im
ityp(1:imm)=config_real(iconf)%itype(1:imm)
xp(1:3,1:imm)=A2cm*config_real(iconf)%pos_cart(1:3,1:imm)
!xc(3,1:imm)=A2cm*config_real(iconf)%pos_crst(3,1:imm)
fp(1:3,1:imm)= config_real(iconf)%force(1:3,1:imm) / (A2cm*erg2ev)
at(:,:)=config_real(iconf)%cell(:,:) * A2cm
bg(:,:)=config_real(iconf)%bg_cell(:,:) / A2cm
cm(:)=config_real(iconf)%mass_per_type(:)*umass
volu=config_real(iconf)%volume*A2cm**3
ntyp=config_real(iconf)%ntypes
!future xpp(1:3,1:imm)=config_real(iconf)%prev_pos_cart(1:3,1:imm)*A2cm

if (debug_ph) then
    if (rangph==0 ) write(6,*) '>...exit in  put_ph_config_into_ndn '
end if


return
end subroutine put_ph_config_into_ndm
!<-------------------------------------------------------------


!$-------------------------------------------------------------
subroutine put_ndm_into_ph_config(iconf)
!$-------------------------------------------------------------
use gen_com_m, only: im, at, bg, imm, A2cm, erg2ev, umass, volu
use tab_imm_m, only : xp, fp, ityp
use var_pot, only: ntyp, cm
use derived_types_ph, only:config_real
use phondy_in_ndm_module, only : rcut_ph, rcut_ph_ang, nat,debug_ph, rangph
implicit none
integer, intent(in) :: iconf

if (debug_ph) then
    if (rangph==0 ) write(6,*) '>...enter in  put_ndm_into_ph_config '
end if


config_real(iconf)%im = im
config_real(iconf)%imm = imm
config_real(iconf)%volume = volu/A2cm**3
config_real(iconf)%ntypes=ntyp
config_real(iconf)%nat=nat

config_real(iconf)%itype(1:imm) = ityp(1:imm)
config_real(iconf)%pos_cart(1:3,1:imm) = xp(1:3,1:imm)/A2cm
!config_real(iconf)%pos_crst(3,1:imm) = xc(3,1:imm)/A2cm
config_real(iconf)%force(1:3,1:imm) = fp(1:3,1:imm)*(A2cm*erg2ev)
config_real(iconf)%cell = at(:,:)/A2cm
config_real(iconf)%bg_cell = bg(:,:)*A2cm
config_real(iconf)%mass_per_type(:) = cm(:)/umass
!future config_real(iconf)%prev_pos_cart(1:3,1:imm)=xpp(1:3,1:imm)/A2cm


rcut_ph_ang = rcut_ph/A2cm


if (debug_ph) then
    if (rangph==0 ) write(6,*) '>...exit in  put_ndm_into_ph_config '
end if



return
end subroutine put_ndm_into_ph_config
!<-------------------------------------------------------------



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
use phondy_in_ndm_module, only: rcut_ph_ang, rangph, debug_ph
use derived_types_ph, only: config_real
implicit none
integer, intent(in) :: iconf
real(double) :: bval(3)
integer :: i


if (debug_ph) then
    if (rangph==0 ) write(6,*) '>...enter in  test_if_config_is_small '
end if



call recips (config_real(iconf)%cell(1,1), config_real(iconf)%cell(1,2), config_real(iconf)%cell(1,3), &
             config_real(iconf)%bg_cell(1,1), config_real(iconf)%bg_cell(1,2), config_real(iconf)%bg_cell(1,3))

do i=1,3
   bval(i)=1.d0/sqrt(sum(config_real(iconf)%bg_cell(:,i)**2))
end do
config_real(iconf)%nxCell=int(2.d0*rcut_ph_ang/bval(1))+1
config_real(iconf)%nyCell=int(2.d0*rcut_ph_ang/bval(2))+1
config_real(iconf)%nzCell=int(2.d0*rcut_ph_ang/bval(3))+1

if (max(config_real(iconf)%nxCell, config_real(iconf)%nyCell, config_real(iconf)%nzCell) .gt.1) config_real(iconf)%small=.true.

if (rangph==0) write(6,'("PHONDY: cut off in Ang is defined to: ", f20.10)') rcut_ph_ang

!if (debug) then
if (rangph==0) write(*,'("PHONDY: small or big box in test_if_config_is_small ",  i5,  3i4, l3)') iconf,  config_real(iconf)%nxCell, &
                                                                                config_real(iconf)%nyCell, &
                                                                                config_real(iconf)%nzCell, &
                                                                                config_real(iconf)%small
!end if

if (debug_ph) then
    if (rangph==0 ) write(6,*) '>...exit in test_if_config_is_small '
end if


return
end subroutine test_if_config_is_small


subroutine ph_build_box
use gen_com_m, only : im, imm
use phondy_in_ndm_module, only : rangph,  iconf_big, iconf_ini, debug_ph
use derived_types_ph, only : config_real !, atom_ph, inv_atom_ph
implicit none
real(kind(0.d0)) :: Rtemp(3)
real(kind(0.d0)), dimension(3,3) :: at_inv
real(kind(0.d0)), dimension(:,:), allocatable :: xpnp
integer :: icount, ia, i1,i2,i3, nn1, nn2, nn3, iia
real(kind(0.d0)), dimension(:,:), allocatable :: xc_temp


if (debug_ph) then
    if (rangph==0 ) write(6,*) '>...enter in  ph_build_box '
end if



!if (debug) call print_message('Enter in subroutine...','ph_build_box_ml')
nn1=int(config_real(iconf_ini)%nxCell/2)!+1
nn2=int(config_real(iconf_ini)%nyCell/2)!+1
nn3=int(config_real(iconf_ini)%nzCell/2)!+1

config_real(iconf_big)%im=(2*nn1+1)*(2*nn2+1)*(2*nn3+1)*config_real(iconf_ini)%nat
config_real(iconf_big)%imm= config_real(iconf_big)%im


!config_real_copy(iconf) = config_real(iconf)

if (allocated(config_real(iconf_big)%pos_cart)) deallocate(config_real(iconf_big)%pos_cart) ; allocate (config_real(iconf_big)%pos_cart(3,config_real(iconf_big)%imm))
if (allocated(config_real(iconf_big)%force)) deallocate(config_real(iconf_big)%force) ; allocate (config_real(iconf_big)%force(3,config_real(iconf_big)%imm))
if (allocated(config_real(iconf_big)%itype)) deallocate(config_real(iconf_big)%itype) ; allocate (config_real(iconf_big)%itype(config_real(iconf_big)%imm))
if (allocated(config_real(iconf_big)%ia_ini)) deallocate(config_real(iconf_big)%ia_ini) ; allocate (config_real(iconf_big)%ia_ini(config_real(iconf_big)%imm))
if (allocated(config_real(iconf_big)%Rperiodic)) deallocate(config_real(iconf_big)%Rperiodic) ; allocate (config_real(iconf_big)%Rperiodic(3,config_real(iconf_big)%imm))

config_real(iconf_big)%nat=config_real(iconf_big)%im
config_real(iconf_big)%volume = config_real(iconf_ini)%volume*config_real(iconf_big)%nat/config_real(iconf_ini)%nat
config_real(iconf_big)%ntypes=config_real(iconf_ini)%ntypes


!if (allocated(xpnp)) deallocate(xpnp); allocate (xpnp(3,imm))

icount=0
do i1=0, 2*nn1
do i2=0, 2*nn2
do i3=0, 2*nn3
      do  ia=1,config_real(iconf_ini)%nat

        icount = icount + 1
        Rtemp(:) =   dble(i1)*config_real(iconf_ini)%cell(:,1) + &
                     dble(i2)*config_real(iconf_ini)%cell(:,2) + &
                     dble(i3)*config_real(iconf_ini)%cell(:,3)
        config_real(iconf_big)%pos_cart(:,icount)=Rtemp(:) + config_real(iconf_ini)%pos_cart(:,ia)
        config_real(iconf_big)%itype(icount) = config_real(iconf_ini)%itype(ia)
        config_real(iconf_big)%ia_ini(icount) = ia
        config_real(iconf_big)%Rperiodic(:,icount) = Rtemp(:)
      end do
end do
end do
end do




config_real(iconf_big)%cell(:,1) = config_real(iconf_ini)%cell(:,1)*dble(2*nn1+1)
config_real(iconf_big)%cell(:,2) = config_real(iconf_ini)%cell(:,2)*dble(2*nn2+1)
config_real(iconf_big)%cell(:,3) = config_real(iconf_ini)%cell(:,3)*dble(2*nn3+1)

if (rangph==0) then

    !write(6,'("positions of atoms 1000 in big ", 3f20.8)') config_real(iconf_big)%pos_cart(:,1000)

    write(6,*) '... The big box ... with ... ', (2*nn1+1)*(2*nn2+1)*(2*nn3+1)*config_real(iconf_ini)%nat, ' ...atoms'
    write(6,'("a1 ", 3f20.8)') config_real(iconf_big)%cell(1,1), config_real(iconf_big)%cell(2,1), config_real(iconf_big)%cell(3,1)
    write(6,'("a2 ", 3f20.8)') config_real(iconf_big)%cell(1,2), config_real(iconf_big)%cell(2,2), config_real(iconf_big)%cell(3,2)
    write(6,'("a3 ", 3f20.8)') config_real(iconf_big)%cell(1,3), config_real(iconf_big)%cell(2,3), config_real(iconf_big)%cell(3,3)
end if
config_real(iconf_big)%nat = icount
if (icount /= (2*nn1+1)*(2*nn2+1)*(2*nn3+1)*config_real(iconf_ini)%nat ) then
  write(6,*) "big problem in the box replication"
  stop 'wrong dimension in <ph_build_box>'
end if

!C_DEBUG the order in cell and bg_cell
call recips (config_real(iconf_big)%cell(1,1), config_real(iconf_big)%cell(1,2), config_real(iconf_big)%cell(1,3), &
             config_real(iconf_big)%bg_cell(1,1), config_real(iconf_big)%bg_cell(1,2), config_real(iconf_big)%bg_cell(1,3))

config_real(iconf_big)%mass_per_type(:) = config_real(iconf_ini)%mass_per_type(:)

!C_debug WRITE GIN FOR DEBUG PURPOSES ...
!C_debug if (allocated(xc_temp)) deallocate(xc_temp) ; allocate(xc_temp(3,config_real(iconf_big)%imm))
!C_debug call MatInv(config_real(iconf_big)%cell(:,:), at_inv)
!C_debug xc_temp(1:3,1:config_real(iconf_big)%imm) = MatMul( at_inv(1:3,1:3),  config_real(iconf_big)%pos_cart(1:3,1:config_real(iconf_big)%imm) )
!C_debug if (rangph==0) then
!C_debug write(55,*)'1  1  1'
!C_debug write(55,'(3f20.10)') config_real(iconf_big)%cell(1,1), config_real(iconf_big)%cell(2,1), config_real(iconf_big)%cell(3,1)
!C_debug write(55,'(3f20.10)') config_real(iconf_big)%cell(1,2), config_real(iconf_big)%cell(2,2), config_real(iconf_big)%cell(3,2)
!C_debug write(55,'(3f20.10)') config_real(iconf_big)%cell(1,3), config_real(iconf_big)%cell(2,3), config_real(iconf_big)%cell(3,3)
!C_debug write(55,'(i9)') config_real(iconf_big)%nat
!C_debug do iia=1,config_real(iconf_big)%nat
!C_debug     write(55,'(3f20.15, " 1 ")') xc_temp(:,iia)
!C_debug end do
!C_debug end if

if (debug_ph) then
    if (rangph==0 ) write(6,*) '>...exit in  ph_build_box '
end if

return
end subroutine ph_build_box



subroutine  allocate_object_ph
use phondy_in_ndm_module, only : rcut_ph_ang, NN_MAX, iconf_ini,debug_ph, rangph
use derived_types_ph, only: config_real
implicit none

!local
integer :: ia, ja, nn1, nn2, nn3, i1, i2, i3, icount
real(kind(0.d0)), parameter  :: zero=1.d-10
real(kind(0.d0))  :: Rtemp(3), norm, Rmax2
integer, dimension(:,:,:,:), allocatable :: itab


if (debug_ph) then
    if (rangph==0 ) write(6,*) '>...enter in  allocate_object_ph '
end if

Rmax2=rcut_ph_ang**2



!ca nn1=int(config_real(iconf_ini)%nxCell/2)+1
!ca nn2=int(config_real(iconf_ini)%nyCell/2)+1
!ca nn3=int(config_real(iconf_ini)%nzCell/2)+1

!ca if (allocated(itab)) deallocate(itab) ; allocate(itab(config_real(iconf_ini)%nat, -nn1:nn1, -nn2:nn2, -nn3:nn3))

!ca itab(:,:,:,:)=0

if (allocated(config_real(iconf_ini)%type_neigh)) deallocate(config_real(iconf_ini)%type_neigh) ; allocate(config_real(iconf_ini)%type_neigh(config_real(iconf_ini)%nat,NN_MAX))
if (allocated(config_real(iconf_ini)%kind_neigh)) deallocate(config_real(iconf_ini)%kind_neigh) ; allocate(config_real(iconf_ini)%kind_neigh(config_real(iconf_ini)%nat,NN_MAX))
if (allocated(config_real(iconf_ini)%r_ij))       deallocate(config_real(iconf_ini)%r_ij) ;       allocate(config_real(iconf_ini)%r_ij(config_real(iconf_ini)%nat,NN_MAX))
if (allocated(config_real(iconf_ini)%u_ij))       deallocate(config_real(iconf_ini)%u_ij) ;       allocate(config_real(iconf_ini)%u_ij(config_real(iconf_ini)%nat,NN_MAX,3))
if (allocated(config_real(iconf_ini)%uperiod_ij))       deallocate(config_real(iconf_ini)%uperiod_ij) ;       allocate(config_real(iconf_ini)%uperiod_ij(config_real(iconf_ini)%nat,NN_MAX,3))
if (allocated(config_real(iconf_ini)%n_neigh))    deallocate(config_real(iconf_ini)%n_neigh) ;    allocate(config_real(iconf_ini)%n_neigh(config_real(iconf_ini)%nat))
if (allocated(config_real(iconf_ini)%kind_neigh_big)) deallocate(config_real(iconf_ini)%kind_neigh_big) ; allocate(config_real(iconf_ini)%kind_neigh_big(config_real(iconf_ini)%nat,NN_MAX))

config_real(iconf_ini)%n_neigh(:) =0
config_real(iconf_ini)%type_neigh(:,:)=0
config_real(iconf_ini)%kind_neigh(:,:)=0
config_real(iconf_ini)%kind_neigh_big(:,:)=0
config_real(iconf_ini)%r_ij(:,:)=0.d0
config_real(iconf_ini)%u_ij(:,:,:)=0.d0

!ca do  ia=1,config_real(iconf_ini)%nat
!ca   icount=0
!ca   do  ja=1,config_real(iconf_ini)%nat
!ca       do i1=-nn1,nn1
!ca       do i2=-nn2,nn2
!ca       do i3=-nn3,nn3
!ca         Rtemp(:) = dble(i1)*config_real(iconf_ini)%cell(:,1) + &
!ca                     dble(i2)*config_real(iconf_ini)%cell(:,2) + &
!ca                     dble(i3)*config_real(iconf_ini)%cell(:,3)
!ca
!ca         Rtemp(:)=Rtemp(:) + config_real(iconf_ini)%pos_cart(:,ja) - config_real(iconf_ini)%pos_cart(:,ia)
!ca         norm=DOT_PRODUCT(Rtemp(:),Rtemp(:))
!ca         ! should be added the RCut ....
!ca         !if (norm.lt.1d-15) cycle
!ca         if (norm.gt.Rmax2) cycle
!ca            icount=icount+1
!ca            !!!o Rn(ia,icount,:)=Rtemp(:)
!ca          !itab(ja,i1,i2,i3)=icount
!ca            !!!o neigh_type(ia,icount)=ja
!ca            !!!o Rn_period(ia,icount,:)=Rn(ia,icount,:)+tau(:,ia)-tau(:,ja)
!ca            config_real(iconf_ini)%type_neigh(ia,icount)=config_real(iconf_ini)%itype(ja)
!ca            config_real(iconf_ini)%kind_neigh(ia,icount)=ja
!ca            config_real(iconf_ini)%r_ij(ia,icount)=dsqrt(norm)
!ca            config_real(iconf_ini)%u_ij(ia,icount,:)=-Rtemp(:)!/config_real(iconf_ini)%r_ij(ia,c)
!ca          itab(ja,i1,i2,i3)=icount
!ca
!ca       end do
!ca       end do
!ca       end do
!ca   end do
!ca   !!!!o ndir2(ia)=icount
!ca   config_real(iconf_ini)%n_neigh(ia)=icount
!ca end do
if (debug_ph) then
    if (rangph==0 ) write(6,*) '>...exit in  allocate_object_ph '
end if


return
end subroutine allocate_object_ph
!*********************************************************************!
!*********************************************************************!

subroutine dyn_mat_q(q_point_in_b)
use derived_types_ph, only: config_real
use phondy_in_ndm_module, only: nmat, dynmat, matfor, NN_MAX,iconf_ini,iconf_big, debug_ph,rangph, unit_nu
implicit none
real(kind(0.d0)), dimension(3), intent(in) :: q_point_in_b
complex(kind(1.d0)), dimension(:,:), allocatable :: phase_factor
complex(kind(1.d0)) :: ci,ctemp, ctemp2, ctemp1
real(kind(0.d0)) :: c1, c2, c3, pi, rtemp
integer :: ia, ja, ialpha, jbeta,imat, jmat,inn, jj
character (len=60) :: CHFMT

nmat=3*config_real(iconf_ini)%nat
dynmat(:,:)=cmplx(0.d0,0.d0)

if (debug_ph) then
    if (rangph==0 ) write(6,*) '>...enter in  dyn_mat_q '
end if


if (allocated(dynmat))       deallocate(dynmat)       ; allocate(dynmat(nmat,nmat))
if (allocated(phase_factor)) deallocate(phase_factor) ; allocate(phase_factor(config_real(iconf_ini)%nat,NN_MAX))

pi=4.d0*atan(1.d0)
ci=cmplx(0.d0,1.d0)

!C        write(6,'("g1 ", 3e20.7)') config_real(iconf_ini)%bg_cell(:,1)
!C        write(6,'("g2 ", 3e20.7)') config_real(iconf_ini)%bg_cell(:,2)
!C        write(6,'("g3 ", 3e20.7)') config_real(iconf_ini)%bg_cell(:,3)
!C        write(6,'("XK ", 3e20.7)') q_point_in_b(:)

do ia=1,config_real(iconf_ini)%nat
  do inn=1,config_real(iconf_ini)%n_neigh(ia)
   jj=config_real(iconf_ini)%kind_neigh_big(ia,inn)
   !write(*,*) rangph, config_real(iconf_ini)%n_neigh(ia), size(config_real(iconf_ini)%kind_neigh(:,:),DIM=2)
   !c1=DOT_PRODUCT(config_real(iconf_big)%Rperiodic(:,jj), config_real(iconf_ini)%bg_cell(:,1))
   !c2=DOT_PRODUCT(config_real(iconf_big)%Rperiodic(:,jj), config_real(iconf_ini)%bg_cell(:,2))
   !c3=DOT_PRODUCT(config_real(iconf_big)%Rperiodic(:,jj), config_real(iconf_ini)%bg_cell(:,3))

   !write(*,*) rangph, config_real(iconf_ini)%n_neigh(ia), size(config_real(iconf_ini)%kind_neigh(:,:),DIM=2)
   c1=DOT_PRODUCT(config_real(iconf_ini)%uperiod_ij(ia,inn,:), config_real(iconf_ini)%bg_cell(:,1))
   c2=DOT_PRODUCT(config_real(iconf_ini)%uperiod_ij(ia,inn,:), config_real(iconf_ini)%bg_cell(:,2))
   c3=DOT_PRODUCT(config_real(iconf_ini)%uperiod_ij(ia,inn,:), config_real(iconf_ini)%bg_cell(:,3))

   phase_factor(ia,inn) = exp(2.d0*pi*ci*(c1*q_point_in_b(1) + c2*q_point_in_b(2) + c3*q_point_in_b(3)))
  end do
end do

do ia=1,config_real(iconf_ini)%nat
  do ja=1,config_real(iconf_ini)%nat

    do  ialpha=1,3
      do  jbeta=1,3
        imat=ialpha+(ia-1)*3
        jmat=jbeta+(ja-1)*3
        ctemp=cmplx(0.d0,0.d0)
        do inn=1,config_real(iconf_ini)%n_neigh(ia)
                 if(config_real(iconf_ini)%kind_neigh(ia,inn).eq.ja) then
                    ctemp=ctemp+matfor(imat,3*(inn-1)+jbeta)*phase_factor(ia,inn)
                 endif
        end do
        dynmat(imat,jmat)=ctemp
      end do
    end do

  end do
 end do


if (debug_ph) then

  do imat=1,size(dynmat(:,:), DIM=1)
     write(CHFMT,*)'(i6, 1x, ',  2*size(dynmat(:,:), DIM=2)   ,'e20.10)'
     write(38, CHFMT), imat,  dynmat(imat, :)*unit_nu*1.d+3
  end do


 !testing and applying the sum rule
 do imat=1,nmat
   ctemp = SUM(dynmat(imat,:))-dynmat(imat,imat)
   rtemp = abs(ctemp + dynmat(imat,imat))
   !write(33,*) imat, rtemp
   !Cif ( rtemp.gt.1d-10) write(6,'("PHONDY warning: The line sum rule is not respected imat, deviance dynmat",i9,d20.10)') imat, rtemp
   !dynmat(imat,imat) = -ctemp
   ctemp1 = -ctemp

   ctemp2 = SUM(dynmat(:,imat))-dynmat(imat,imat)
   rtemp = abs(ctemp2 + dynmat(imat,imat))
   if ( rtemp.gt.1d-10) write(6,'("PHONDY warning: The column sum rule is not respected imat, deviance dynmat",i9,d20.10)') imat, rtemp
   !write(34,*) imat, rtemp
   !write(35,'(4e20.7)') SUM(dynmat(imat,:)) ,  SUM(dynmat(:,imat))
 end do

do imat = 1,nmat
  do jmat = imat, nmat
    !write (36,'(6e20.7)') dynmat(imat,jmat) / dynmat(jmat,imat), dynmat(imat,jmat), dynmat(jmat,imat)
  end do
end do
end if



if (debug_ph) then
    if (rangph==0 ) write(6,*) '>...exit in  dyn_mat_q '
end if

return
end subroutine dyn_mat_q


subroutine diago_q
  USE T_kind_param_m, ONLY:  double
use phondy_in_ndm_module, only: nmat, rangph, q_point, debug_ph, no_of_qpoints, q_eigenvalues, W, unit_nu, iconf_ini, dynmat
use derived_types_ph, only : config_real
implicit none
real(kind(1.d0)) :: eigenvalues(nmat)
integer :: ik,n
character (len=60) :: CHFMT
integer :: count1,count2,count3, count_rate,count_max
real(kind(1.d0)) :: time_d, time_m, time
integer :: i

  call system_clock (count1,count_rate,count_max)
  time_m=0.d0
  time_d=0.d0

if (debug_ph) then
    if (rangph==0 ) write(6,*) '>...enter in  diago_q '
end if


nmat=3*config_real(iconf_ini)%nat
 if (allocated(q_eigenvalues)) deallocate(q_eigenvalues) ; allocate(q_eigenvalues(nmat, no_of_qpoints))


 do ik=1,no_of_qpoints
  !
  if (debug_ph) write(*,*) 'ik point', ik, nmat
  !
  if (debug_ph) write(*,*) "dyn_mat_q ... in diago_q "
  call dyn_mat_q(q_point(:,ik))

  call system_clock (count2,count_rate,count_max)
  time=real((count2-count1))/real(count_rate)
  time_m = time_m + time
  if (debug_ph) write(*,*) "threading complex  ... "
  call diago_threading_complex

  call system_clock (count3,count_rate,count_max)

  time=real((count3-count2))/real(count_rate)
  time_d = time_d + time

  eigenvalues(:)=sign(1.d0,W(:))*dsqrt(dabs(W(:)))*unit_nu*1.D3
  q_eigenvalues(:,ik)=eigenvalues(:)

 end do

write(6,'("PHONDY: time for build Dq matrix and diagonalization, respectivelly: ", 2f30.4)') time_m, time_d



if (debug_ph) then
    if (rangph==0 ) write(6,*) '>...exit in  diago_q '
end if


!     call write_ldos(iam, N,W,Z,DESCZ)  ! all others variables are sended by the module phony_in_ndm_module
!     call write_lmodes(iam,N,W,Z,DESCZ)
end subroutine diago_q
