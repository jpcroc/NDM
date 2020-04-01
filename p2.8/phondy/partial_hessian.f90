subroutine partial_hessian(i_start,i_final,nlocal,matforlocal,u_local,v_local,HessianOrder)

  !use mpi
  !use mod_mpi_phondy

  USE T_kind_param_m, ONLY:  double
  use gen_com_m, ONLY: im,imm,A2cm,erg2ev,lperiod,bg,at,potist,indi,indi2,rangph
  use var_pot
  use tab_imm_m, ONLY : ityp,iwmax,iwmax2,fp,xp
  use phondy_in_ndm_module, ONLY: units_ndm,deltax,fp0,matfor,epot0,rcut_ph,it_phondy,cm_phondy_at, lq_points, iconf_ini, debug_ph
  use derived_types_ph, only : config_real

 implicit none
  integer, intent(in):: i_start,i_final,nlocal,HessianOrder
  real(double),intent(out) :: matforlocal(nlocal)
  integer, intent(out)   :: u_local(nlocal),v_local(nlocal)

  real(double),allocatable, dimension(:,:)  :: fpp1,fpp2,fpm1,fpm2
  real(double)::rcut_ph2,rtest2
  real(double), dimension(:,:), allocatable :: xpnp
  real(double), dimension(1:3) :: dxp, ds
  real(double) :: xpref(3,imm)
  integer :: i,j,u,v,ia,jb,iw,iw1,iw2,iti,icount
  integer :: ilocal,iconf

  !iconf = iconf_data_ph
  rcut_ph2=rcut_ph**2

  !debug  if (rangph==0) write(*,*) 'hRUE.....(SHOULD BE  AT LEAST 2 x Rcut ...: ', rue
  !debug  if (rangph==0) write(*,*) 'hRUE.POTENTIEL............................: ', rue_pot(ipotentiel)



  if (rangph==0) write(6,'("PHONDY: u, v,  mlocal has the size inside partial hessian: ",i12, i12, i12)') size(u_local,DIM=1), size(v_local,DIM=1),size(matforlocal,DIM=1)
  allocate(xpnp(3,imm))
  if (lperiod) then
   xpnp(:,:)=xp(:,:)
  else
   call notperiod(xp,xpnp)
  end if
!C_debug   call cryst_to_cart (imm, xpnp, bg, -1)

 !Here we apply a coefficient in order to put the Hessian in eV/ A²
  !coeff=
  if (rangph==0) write (6, *) 'PHONDY: Computing the force constants .... '

  xpref(:,:) = xp(:,:)

 it_phondy=1
 matforlocal(:)=0.d0

  ilocal=0

  if (HessianOrder.eq.1) then
  !
  if (i_start==1)  iw2=0
  if (i_start > 1) iw2=iwmax2(i_start-1)
  !
  do j=i_start,i_final
!ca   if (lq_points) then
!ca    iw1=1
!ca    iw2=config_real(iconf)%n_neigh(j)
!ca    iti=config_real(iconf)%itype(j)
!ca   else
    iti= ityp(j)
    iw1=iw2+1
    iw2=iwmax2(j)
!ca   end if
     do jb=1,3
        xp(:,:) = xpref(:,:)
        u=3*(j-1)+jb
        xp(jb,j)=xpref(jb,j)+deltax
        call calfo_phondy(j,j)

        if (dabs(potist*erg2eV-epot0).gt.0.5) then
         write(*,*) 'PHONDY: Hessian O(1) WARNING the energy is greater than 0.5. Decrease deltax', (potist*erg2eV-epot0)
        end if
        icount=0
        do iw=iw1,iw2
!ca             if (lq_points) then
!ca              i=config_real(iconf)%kind_neigh(j,iw)
!ca             else
          i=indi2(iw)
!ca             end if

!ca             if (lq_points) then
!ca               rtest2 = (config_real(iconf)%r_ij(j,iw)*A2cm)**2
!ca               dxp(1:3) = config_real(iconf)%u_ij(j,iw,:)*A2cm
!ca             else
            dxp(1:3) = xpnp(1:3,i) - xpnp(1:3,j)
            ds=MatMul(dxp,bg)
            WHERE ( (ds.GT.0.5d0).OR.(ds.LT.-0.5d0) )
               ds(1:3) = ds(1:3) - Dble(Nint(ds(1:3)))
            END WHERE
            dxp = MatMul(at,ds)
            rtest2 = Sum( dxp(1:3)**2 )
!ca             end if
            if (rtest2<rcut_ph2)  then
              icount = icount + 1

              do ia=1,3
                if (lq_points) then
                  v = 3*(icount-1) + ia
                else
                  v=3*(i-1)+ia
                end if
                ilocal=ilocal+1
                matforlocal(ilocal)=-1.*(fp(ia,i)-fp0(ia,i))/ &
                (deltax*dsqrt(cm_phondy_at(i))*dsqrt(cm_phondy_at(j)))*units_ndm
                u_local(ilocal)=u
                v_local(ilocal)=v
              end do

            end if
        end do !iw

     end do !jb
  end do    !j
 !
end if


 if (HessianOrder.eq.2) then

 if (allocated(fpp1)) deallocate(fpp1) ; allocate(fpp1(3,imm))
 if (allocated(fpm1)) deallocate(fpm1) ; allocate(fpm1(3,imm))

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
      iti= ityp(j)
      iw1=iw2+1
      iw2=iwmax2(j)
!ca    end if

    if (debug_ph) then
      if (mod(j-1,5)==0) then
        if (rangph==0) write(*,*) 'partial hessian calc on', j, iw2-iw1
      end if
    end if

    do jb=1,3
        !debug Cosmin ... should be added that ? (next line)
        xp(:,:) = xpref(:,:)
        u=3*(j-1)+jb
        !other_good_way_to_pack: u=j+im*(jb-1)

        xp(jb,j)=xpref(jb,j)+deltax


        !if (rangph==0) write(*,*) 'partial_hessian.f90 test1.........:', rangph, j

        call calfo_phondy(j,j)


        !debug if (rangph==0) write(*,*) 'partial_hessian.f90 test2.........:', j, jb, potist*erg2eV

        if (dabs(potist*erg2eV-epot0).gt.0.5) then
         write(*,*) 'WARNING the energy is greater than 0.5. Decrease deltax', (potist*erg2eV-epot0)
        end if

        fpp1(:,:)=fp(:,:)


        xp(:,:) = xpref(:,:)
        xp(jb,j)=xpref(jb,j)-deltax

        call calfo_phondy(j,j)

        !if (rangph==0) write(*,*) 'partial_hessian.f90 test3.........:', rangph, j

        if (dabs(potist*erg2eV-epot0).gt.0.5) then
         write(*,*) 'WARNING the energy is greater than 0.5. Decrease deltax',(potist*erg2eV-epot0)
        end if
        fpm1(:,:)=fp(:,:)
        icount=0
        do iw=iw1,iw2

!ca            if (lq_points) then
!ca              i=config_real(iconf)%kind_neigh(j,iw)
!ca            else
              i=indi2(iw)
!ca            end if

!ca            if (lq_points) then
!ca               rtest2 = (config_real(iconf)%r_ij(j,iw)*A2cm)**2
!ca               dxp(1:3) = config_real(iconf)%u_ij(j,iw,:)*A2cm
!ca            else
              dxp(1:3) = xpnp(1:3,i) - xpnp(1:3,j)
              !if (rangph==0) write(32,'(3i8,e15.5)') j, jb , i,  sqrt(SUM(dxp(:)**2))/A2cm
              ds=MatMul(dxp,bg)
              WHERE ( (ds.GT.0.5d0).OR.(ds.LT.-0.5d0) )
               ds(1:3) = ds(1:3) - Dble(Nint(ds(1:3)))
              END WHERE
              dxp = MatMul(at,ds)
              rtest2 = Sum( dxp(1:3)**2 )
!ca            end if

            if (rtest2<rcut_ph2)  then
              icount=icount + 1
              do ia=1,3
                if (lq_points) then
                  v=3*(icount-1)+ia
                else
                  v=3*(i-1)+ia
                end if
                ilocal=ilocal+1
                !other_good_way_to_pack: v=i+im*(ia-1)
                !if (rangph==0) write(33,'(4i8,2e15.5)') j, jb , i , ia, sqrt(rtest2)/A2cm, sqrt(rcut_ph2)/A2cm
                 matforlocal(ilocal)=-1.d0*(fpp1(ia,i)-fpm1(ia,i))/ ( (2.d0*deltax) &
                * dsqrt(cm_phondy_at(i))*dsqrt(cm_phondy_at(j)))* units_ndm
                u_local(ilocal)=u
                v_local(ilocal)=v
              end do
            end if

        end do !iw
        if (lq_points) then
          if (abs(icount - config_real(iconf_ini)%n_neigh(j)) .gt. 0) then
             write(6,*) 'WARNNNNNINNNG ', icount,config_real(iconf_ini)%n_neigh(j)
          end if
        end if
     end do  ! jb that comes from xyz
      !if (rangph==0)  write(*,*) 'partial_hessian.f90 testF.........:', rangph, j
  end do ! j the index of atom

   !if (rangph==0) write(*,*) 'partial_hessian.f90 test3.........:', rangph, i_start, i_final
 !
 deallocate(fpp1,fpm1)
 end if


 if (HessianOrder.eq.4) then


 if (allocated(fpp1)) deallocate(fpp1)
 if (allocated(fpp2)) deallocate(fpp2)
 if (allocated(fpm1)) deallocate(fpm1)
 if (allocated(fpm2)) deallocate(fpm2)
 allocate(fpp1(3,imm),fpp2(3,imm),fpm1(3,imm), fpm2(3,imm))
 !
  if (i_start==1)  iw2=0
  if (i_start > 1) iw2=iwmax2(i_start-1)
  !
  do j=i_start,i_final  ! j in atom index

!ca    if (lq_points) then
!ca      iw1=1
!ca      iw2=config_real(iconf)%n_neigh(j)
!ca      iti=config_real(iconf)%itype(j)
!ca    else
      iti= ityp(j)
      iw1=iw2+1
      iw2=iwmax2(j)
!ca    end if
!debug    if ((mod (j,100)==0)) write(*,*) '.......j ',j,iw2-iw1,deltax
    do jb=1,3
        u=3*(j-1)+jb
        !other_good_way_to_pack u=j+im*(jb-1)

        xp(:,:) = xpref(:,:)
        xp(jb,j)=xpref(jb,j)+deltax
        call calfo_phondy(1,im)
        fpp1(:,:)=fp(:,:)
        if (dabs(potist*erg2eV-epot0).gt.0.5) then
         write(*,*) 'PHONDY: Hessian O(4) 1 WARNING the energy is greater than 0.5. Decrease deltax', (potist*erg2eV-epot0)
        end if


        xp(:,:) = xpref(:,:)
        xp(jb,j)=xpref(jb,j)+2.d0*deltax
        call calfo_phondy(1,im)
        fpp2(:,:)=fp(:,:)
        if (dabs(potist*erg2eV-epot0).gt.0.5) then
         write(*,*) 'PHONDY: Hessian O(4) 2 WARNING the energy is greater than 0.5. Decrease deltax', (potist*erg2eV-epot0)
        end if



        xp(:,:) = xpref(:,:)
        xp(jb,j)=xpref(jb,j)-deltax
        call calfo_phondy(1,im)
        fpm1(:,:)=fp(:,:)
        if (dabs(potist*erg2eV-epot0).gt.0.5) then
         write(*,*) 'PHONDY: Hessian O(4) 3 WARNING the energy is greater than 0.5. Decrease deltax', (potist*erg2eV-epot0)
        end if


        xp(:,:) = xpref(:,:)
        xp(jb,j)=xpref(jb,j)-2.d0*deltax
        call calfo_phondy(1,im)
        fpm2(:,:)=fp(:,:)
        if (dabs(potist*erg2eV-epot0).gt.0.5) then
         write(*,*) 'PHONDY: Hessian O(4) 4 WARNING the energy is greater than 0.5. Decrease deltax', (potist*erg2eV-epot0)
        end if
        icount = 0
        do iw=iw1,iw2

!ca            if (lq_points) then
!ca              i=config_real(iconf)%kind_neigh(j,iw)
!ca            else
              i=indi2(iw)
!ca            end if
            !

!ca            if (lq_points) then
!ca               rtest2 = (config_real(iconf)%r_ij(j,iw)*A2cm)**2
!ca               dxp(1:3) = config_real(iconf)%u_ij(j,iw,:)*A2cm
!ca            else
              dxp(1:3) = xpnp(1:3,i) - xpnp(1:3,j)

              ds=MatMul(dxp,bg)
              WHERE ( (ds.GT.0.5d0).OR.(ds.LT.-0.5d0) )
               ds(1:3) = ds(1:3) - Dble(Nint(ds(1:3)))
              END WHERE
              dxp = MatMul(at,ds)
              rtest2 = Sum( dxp(1:3)**2 )
!ca            end if

             !
            if (rtest2<rcut_ph2)  then
              !
              icount = icount + 1
              do ia=1,3
                if (lq_points) then
                  v = 3*(icount-1) + ia
                else
                  v=3*(i-1)+ia
                end if
                ilocal=ilocal+1
                matforlocal(ilocal)=-1.*(fpm2(ia,i)-8.d0*fpm1(ia,i)+8.d0*fpp1(ia,i)-fpp2(ia,i) &
                     )/ ( (12.d0*deltax) &
                * dsqrt(cm_phondy_at(i))*dsqrt(cm_phondy_at(j)) )*units_ndm

                u_local(ilocal)=u
                v_local(ilocal)=v
              end do
              !
              !
            end if
             !
        end do  !iw

     end do !jb is office for xyz
  end do !j is the index of atom
 !
 deallocate(fpp1,fpp2,fpm1,fpm2)
 end if

 if (nlocal /= ilocal) then
  if (rangph==0) write(*,*) 'Problems in partial index. The shape of tables are wrong', ilocal,nlocal
  if (rangph==0) write(*,*) 'Stop in partial_index.f90'
  stop
 end if
  !if (rangph==0) write(*,*) 'in partial hessian', rangph

  !call MPI_BARRIER(MPI_COMM_WORLD,codeph)
deallocate(xpnp)
end subroutine partial_hessian
