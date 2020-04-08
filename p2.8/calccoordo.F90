module calccoordo_mod
        USE notperiod_mod
        USE cryst_to_cart_mod
        USE gen_com_m, ONLY: rang,last,nato,noxyz,it,timel,imm,im,imd,at,deltadist,bg,lperiod,ncel
        implicit none
        contains
subroutine calccoordo
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  USE tab_imm_m
#ifdef PARA
  USE mod_para
#endif
  USE var_pot, ONLY:ntyp,rc,nad
  implicit none
  !-----------------------------------------------
  !   G l o b a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: i, iti, itj, i1, i2, koo, ko1, j, ic,nci,ip,ll
  real(double), dimension(ntyp,ntyp) :: dnco
  real(double) :: r2, c1, c2, c3,cv(1,3),x1,x2,x3
  real(double):: a1, a2, a3

#ifdef PARA
  real(double), dimension(ntyp,ntyp) :: dnco_glob
  integer :: nci_glob
#endif

  real(double),pointer :: xpnp(:,:)

  allocate(xpnp(3,imm))
  !

  !      write(6,*)'entree calcoordo'
  write (6, *)
  write (6, *) '--------- Coordinations ----------------'
  dnco(:ntyp,:ntyp) = 0

  !      write(6,*)' calcoordo2'
  !      a1(:imd-1) = -sign(zl(1),xp(1,:imd-1))
  !      a2(:imd-1) = -sign(zl(2),xp(2,:imd-1))
  !      a3(:imd-1) = -sign(zl(3),xp(3,:imd-1))


  if (lperiod) then
     xpnp(:,:)=xp(:,:)
  else 
     call notperiod(xp,xpnp)
  end if

  do i = 1, imd
     koo = ielat(i)
     nci=0
     iti=ityp(i)
     do i1 = 0, 26
        ko1 = ncel(koo,i1)

        do i2 = 1, nato(ko1)
           j = last(i2,ko1)
           if (i==j) cycle
           c1 = xpnp(1,i)-xpnp(1,j)
           c2 = xpnp(2,i)-xpnp(2,j)
           c3 = xpnp(3,i)-xpnp(3,j)
              c1 = c1+sum(at(1,:)*deltadist(:,i1,koo))
              c2 = c2+sum(at(2,:)*deltadist(:,i1,koo))
              c3 = c3+sum(at(3,:)*deltadist(:,i1,koo))
           if (noxyz.eq.1) then
              cv(1,1) = c1
              cv(1,2) = c2
              cv(1,3) = c3
              call cryst_to_cart (1, cv, bg, -1) !cryst vers cart sur cv
              WHERE ( (cv.GT.0.5d0).OR.(cv.LT.-0.5d0) )
                 cv(:,1:3) = cv(:,1:3) - Dble(Nint(cv(:,1:3)))
              END WHERE
              call cryst_to_cart (1, cv, at, 1) !cryst vers cart sur cv
              c1=cv(1,1)
              c2=cv(1,2)
              c3=cv(1,3)
              
           end if
              
           r2 = sqrt(c1*c1+c2*c2+c3*c3)

           !               if (r2.le.3.0d-15)               write(6,*)i,j,r2

           if (r2<rc(ityp(i)))then
              dnco(ityp(i),ityp(j)) = dnco(ityp(i),ityp(j))+1
              nci=nci+1
           end if
           !               if (r2>=rc(ityp(j))) cycle
           !               dnco(ityp(j),ityp(i)) = dnco(ityp(i),ityp(j))+1

        end do
     end do
           !   write(6,*)i,nci

  end do

#ifdef PARA
  call MPI_ALLREDUCE(dnco(:,:),dnco_glob(:,:),ntyp*ntyp,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
  dnco = dnco_glob
  ! nci n'est pas utilise dans la suite, je laisse en commentaire la reduction
  !  call MPI_ALLREDUCE(nci,nci_glob,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
  !  nci = nci_glob
#endif

  do i1 = 1, ntyp
     if (nad(i1)==0) cycle
     where (nad(:ntyp)/=0) dnco(i1,:ntyp) = dnco(i1,:ntyp)/nad(i1)
  enddo
  if (rang==0) then
     write (6, '(A,I5,A,D10.3)') '*  ITERATION  = ', it, '  time = ', timel
     write (6, *) 'nombre de coor autour de :'
     do i1 = 1, ntyp
        if (nad(i1)==0) cycle
        write (6, *) i1, (dnco(i1,j),j=1,ntyp)
        !   28 FORMAT (31H Nombre de coor autour de ITYP=,I2,2X,7(F5.2,1X))
     end do
     do i1 = 1, ntyp
        if (nad(i1)==0) cycle
        write (6, *) 'rayon autour des type ', i1, ' = ', rc(i1)*1D+8
     end do
  endif

  deallocate (xpnp)
  return

end subroutine calccoordo
end module
