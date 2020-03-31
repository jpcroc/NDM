subroutine calcangle
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m  
  use var_pot
  use tab_imm_m
#if(PARA)
  use mod_mpi
#endif

  !******************************************************************
  implicit none
  !-----------------------------------------------
  !   G l o b a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: i1, i2, i3, i4, koo, ko1, ko2, i, &
       j, k, m,m1,ka,ma,m2
  real(double) :: thetaijk, c11, c21, c31, c12, c22, c32, &
       incre,dij2,dik2,dij,dik
  real(double) :: costheta,invincre,rspace2,rc2(ntyp,ntyp),rc22(ntyp,ntyp),cv(1,3)
#if(PARA)
  real(double) :: fda_glob(ntyp,ntyp,ntyp,contmax)
#endif
  !-----------------------------------------------

  real(double),pointer :: xpnp(:,:)
  allocate (xpnp(3,imm))
  !repartition des atomes entre les petites cel.
!  if (rang==0) write(6,*) 'PARA-T entree calcangle'
  rc2=rcangle*1.0d-8
  rc22=rc2**2
  incre = (thetamax-thetamin)/cont
  invincre = 1/incre

  if (lperiod) then
     xpnp(:,:)=xp(:,:)
  else 
     call notperiod(xp,xpnp)
  end if

  do i = 1, imana-1
     koo = ielat(i)
     do i1 = 0, 26
        ko1 = ncel(koo,i1)
        do i2 = 1, nato(ko1) 
           j = last(i2,ko1)
           c11 = xpnp(1,i)-xpnp(1,j) 
           c21 = xpnp(2,i)-xpnp(2,j) 
           c31 = xpnp(3,i)-xpnp(3,j)
           c11 = c11+sum(at(1,:)*deltadist(:,i1,koo))
           c21 = c21+sum(at(2,:)*deltadist(:,i1,koo))
           c31 = c31+sum(at(3,:)*deltadist(:,i1,koo))
           
           if (noxyz.ne.1) then
              if (abs(c11)>rc2(ityp(i),ityp(j))) cycle
              if (abs(c21)>rc2(ityp(i),ityp(j))) cycle
              if (abs(c31)>rc2(ityp(i),ityp(j))) cycle
           else
              cv(1,1) = c11
              cv(1,2) = c21
              cv(1,3) = c31
              call cryst_to_cart (1, cv, bg, -1) !cryst vers cart sur cv
              WHERE ( (cv.GT.0.5d0).OR.(cv.LT.-0.5d0) )
                 cv(:,1:3) = cv(:,1:3) - Dble(Nint(cv(:,1:3)))
              END WHERE
              call cryst_to_cart (1, cv, at, 1) !cryst vers cart sur cv
              c11=cv(1,1)
              c21=cv(1,2)
              c31=cv(1,3)
           end if


              dij2=c11**2+c21**2+c31**2
              if(dij2>rc22(ityp(i),ityp(j))) cycle
           do i3 = 0,26
              ko2 = ncel(koo,i3)
              do i4 = 1, nato(ko2)
                 k = last(i4,ko2)
                 if(j==i.or.k==i.or.j==k) cycle                 
                 c12 = xpnp(1,i)-xpnp(1,k) 
                 c22 = xpnp(2,i)-xpnp(2,k) 
                 c32 = xpnp(3,i)-xpnp(3,k)
                    c12 = c12+sum(at(1,:)*deltadist(:,i1,koo))
                    c22 = c22+sum(at(2,:)*deltadist(:,i1,koo))
                    c32 = c32+sum(at(3,:)*deltadist(:,i1,koo))
                    if (noxyz.ne.1) then
                       if (abs(c12)>rc2(ityp(i),ityp(k))) cycle
                       if (abs(c22)>rc2(ityp(i),ityp(k))) cycle
                       if (abs(c32)>rc2(ityp(i),ityp(k))) cycle
                    else
                       cv(1,1) = c12
                       cv(1,2) = c22
                       cv(1,3) = c32
                       call cryst_to_cart (1, cv, bg, -1) !cryst vers cart sur cv
                       WHERE ( (cv.GT.0.5d0).OR.(cv.LT.-0.5d0) )
                          cv(:,1:3) = cv(:,1:3) - Dble(Nint(cv(:,1:3)))
                       END WHERE
                       call cryst_to_cart (1, cv, at, 1) !cryst vers cart sur cv
                       c12=cv(1,1)
                       c22=cv(1,2)
                       c32=cv(1,3)
                    end if

!                    if (abs(c12)>rc2(ityp(i),ityp(k))) cycle
!                    if (abs(c22)>rc2(ityp(i),ityp(k))) cycle
!                    if (abs(c32)>rc2(ityp(i),ityp(k))) cycle
                    if(c12**2+c22**2+c32**2>rc22(ityp(i),ityp(k))) cycle

                    costheta = (c11*c12+c21*c22+c31*c32)/ &
                      (sqrt(c11*c11+c21*c21+c31*c31)*sqrt(c12*c12+ &
                      c22*c22+c32*c32))

                 thetaijk = acos(costheta)
                 ka = int(thetaijk*invincre)
                 m = ka+1
                 fda(ityp(j),ityp(i),ityp(k),m) = fda(ityp(j),ityp(i), &
                      ityp(k),m)+1

              end do
           end do
        end do
     end do
  end do

#if(PARA)
  call MPI_ALLREDUCE(fda,fda_glob,ntyp*ntyp*ntyp*contmax,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
  fda = fda_glob
#endif
  deallocate (xpnp)
  return
end subroutine calcangle




