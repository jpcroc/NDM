subroutine partial_take_size (i_start,i_final,imax)
 USE T_kind_param_m, ONLY:  double
  use gen_com_m, ONLY: im,imm,A2cm,lperiod,bg,at,indi,indi2,rangph
  use var_pot
  use tab_imm_m, ONLY : ityp,iwmax,iwmax2,fp,xp
  use phondy_in_ndm_module, ONLY: rcut_ph,it_phondy,iconf_ini, iconf_big, iconf_ini, lq_points
  use derived_types_ph, only : config_real
 implicit none
 integer, intent (in) :: i_start,i_final
 integer, intent(out) :: imax
  real(double)::rcut_ph2,rtest2
  real(double), dimension(:,:), allocatable :: xpnp
  real(double), dimension(1:3) :: dxp, Rtemp, ds
  real(double) :: xpref(3,imm)
  integer :: i,j,u,v,ia,jb,iw,iw1,iw2,iti,jti, ilocal,iconf, icount


    rcut_ph2=rcut_ph**2
    !debug if (rangph==0) write(*,*) 'iRUE.....(SHOULD BE  AT LEAST 2 x Rcut ...: ', rue
    !debug if (rangph==0) write(*,*) 'iRUE.POTENTIEL............................: ', rue_pot(ipotentiel)



  ALLOCATE(xpnp(3,imm))
  if (lperiod) then
   xpnp(:,:)=xp(:,:)
  else
   call notperiod(xp,xpnp)
  end if
  !C_DEBUG call cryst_to_cart (imm, xpnp, bg, -1)

 !Here we apply a coefficient in order to put the Hessian in eV/ A²
  !coeff=
  if (rangph==0) write (6, *) 'PHONDY: buiding the index of the hessian on each proc .... '

  xpref(:,:) = xp(:,:)
  it_phondy=0
! here is aproblem with LAMMPS
  call calfo_phondy(im,im)
  ilocal=0
  !
  if (i_start==1)  iw2=0
  if (i_start > 1) iw2=iwmax2(i_start-1)
  !
  do j=i_start,i_final

!ca    if (lq_points) then
!ca      iw1=1
!ca      iw2=config_real(iconf)%n_neigh(j)
!ca      iti=config_real(iconf)%itype(j)
!ca    else
      jti= ityp(j)
      iw1=iw2+1
      iw2=iwmax2(j)
!ca    end if

     do jb=1,3
        u=3*(j-1)+jb
        !other_good_way_to_pack: u=j+im*(jb-1)
        icount = 0
        do iw=iw1,iw2

!ca            if (lq_points) then
!ca              i=config_real(iconf)%kind_neigh(j,iw)
!ca            else
          i=indi2(iw)
          iti= ityp(i)

!ca            end if

!ca            if (lq_points) then
!ca               rtest2 = (config_real(iconf)%r_ij(j,iw)*A2cm)**2
!ca               dxp(1:3) = config_real(iconf)%u_ij(j,iw,:)*A2cm
!ca            else
          !call recips (at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3))

          dxp(1:3) = xpnp(1:3,i) - xpnp(1:3,j)
          ds=MatMul(dxp,bg)
          WHERE ( (ds.GT.0.5d0).OR.(ds.LT.-0.5d0) )
            ds(1:3) = ds(1:3) - Dble(Nint(ds(1:3)))
          END WHERE
          dxp = MatMul(at,ds)
          rtest2 = Sum( dxp(1:3)**2 )
!ca            end if
          if (rtest2<rcut_ph2)  then
          !debug write(*,'(3i6,d22.10,F12.5,F12.5)') j, i, jb, rtest2, dsqrt(rtest2)/A2cm, dsqrt(rcut_ph2)/A2cm
          !write(*,'(5i6,d20.14)') j, i, jb, iw1,iw2, rtest2

            if (lq_points) then
              icount = icount + 1
              config_real(iconf_ini)%type_neigh(j,icount)=iti
              config_real(iconf_ini)%kind_neigh(j,icount)=config_real(iconf_big)%ia_ini(i)
              config_real(iconf_ini)%kind_neigh_big(j,icount)=i
              config_real(iconf_ini)%r_ij(j,icount)=dsqrt(rtest2)/A2cm
              !config_real(iconf_ini)%u_ij(j,icount,:)=-dxp(:)/A2cm
              config_real(iconf_ini)%uperiod_ij(j,icount,:)=-dxp(:)/A2cm - config_real(iconf_ini)%pos_cart(:,config_real(iconf_big)%ia_ini(j)) + config_real(iconf_ini)%pos_cart(:,config_real(iconf_big)%ia_ini(i)) !/config_real(iconf_ini)%r_ij(ia,c)
              config_real(iconf_ini)%u_ij(j,icount,:)=-dxp(:)/A2cm 
                !config_real(iconf_ini)%Rphase(j,icount) = config_real(iconf_big)%Rperiodic(:,i)
            end if

            do ia=1,3
              v=3*(i-1)+ia
              ilocal=ilocal+1
            end do
          end if

        end do !iw

     end do  !jb  x,y,z
     if (lq_points) config_real(iconf_ini)%n_neigh(j) = icount
  end do !j

!  write(6,'("coordinates for  1000 in  big box NDM ", 3f20.8)')        xp(:,1000)/A2cm

  imax=ilocal

 !write(*,*) 'debug imax', rangph, imax

end subroutine partial_take_size
